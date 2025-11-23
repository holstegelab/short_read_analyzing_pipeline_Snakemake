#!/usr/bin/env python3
import argparse
import sys
import subprocess
from collections import deque, Counter


def open_sam_stream(path, threads=1):
    cmd = ["samtools", "view", "-h"]
    if threads and threads > 1:
        cmd += ["-@", str(threads)]
    cmd.append(path)
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1)
    except Exception as e:
        print(f"Error starting samtools for {path}: {e}", file=sys.stderr)
        raise
    return p


def sam_alignment_iter(proc):
    for raw in proc.stdout:
        if not raw:
            break
        line = raw.rstrip('\n').rstrip('\r')
        if not line:
            continue
        if line[0] == '@':
            continue
        if '\t' not in line:
            continue
        tabs = line.count('\t')
        if tabs < 10:
            continue
        yield line


def parse_sam_line(line):
    fields = line.split('\t')
    fixed = fields[:11]
    tags = fields[11:]
    parsed = []
    for t in tags:
        parts = t.split(':', 2)
        if len(parts) >= 3:
            parsed.append((parts[0], parts[1], parts[2]))
        else:
            parsed.append((t, '', ''))
    parsed.sort(key=lambda x: (x[0], x[1], x[2]))
    canonical = '\t'.join(fixed + [f"{k}:{tp}:{val}" if tp else k for (k, tp, val) in parsed])
    return fixed, parsed, canonical


def qname_of(line):
    i = line.find('\t')
    return line if i == -1 else line[:i]


def qname_base_of(line):
    q = qname_of(line)
    if len(q) > 2 and q[-2] == '/' and q[-1] in ('1', '2'):
        return q[:-2]
    return q


def canonical_norm(line):
    fixed, parsed, _ = parse_sam_line(line)
    fixed = fixed[:]
    fixed[0] = qname_base_of(line)
    return '\t'.join(fixed + [f"{k}:{tp}:{val}" if tp else k for (k, tp, val) in parsed])


class Peekable:
    def __init__(self, it):
        self._it = it
        self.queue = deque()
        self._exhausted = False

    def ensure(self, n):
        while not self._exhausted and len(self.queue) < n:
            try:
                nxt = next(self._it)
            except StopIteration:
                self._exhausted = True
                break
            self.queue.append(nxt)

    def read(self):
        if self.queue:
            return self.queue.popleft()
        if self._exhausted:
            return None
        try:
            return next(self._it)
        except StopIteration:
            self._exhausted = True
            return None

    def skip(self, n):
        for _ in range(n):
            if self.queue:
                self.queue.popleft()
            else:
                # consume from iterator to maintain position
                try:
                    next(self._it)
                except StopIteration:
                    self._exhausted = True
                    break

    def drain_for_qname(self, qn, max_scan):
        self.ensure(max_scan)
        res = []
        newq = deque()
        scanned = 0
        while self.queue and scanned < max_scan:
            l = self.queue.popleft()
            scanned += 1
            if qname_base_of(l) == qn:
                res.append(l)
            else:
                newq.append(l)
        while self.queue:
            newq.append(self.queue.popleft())
        self.queue = newq
        return res


def find_resync(ahead1, ahead2):
    names1 = {}
    names2 = {}
    for i, l in enumerate(ahead1):
        nm = qname_base_of(l)
        if nm not in names1:
            names1[nm] = i
    for j, l in enumerate(ahead2):
        nm = qname_base_of(l)
        if nm not in names2:
            names2[nm] = j
    common = set(names1.keys()) & set(names2.keys())
    if not common:
        return None
    # choose earliest alignment by minimizing max(i,j), then min(i+j)
    best = None
    best_key = None
    for nm in common:
        i = names1[nm]
        j = names2[nm]
        key = (max(i, j), i + j)
        if best is None or key < best_key:
            best = (i, j, nm)
            best_key = key
    return best  # (i, j, name)


def find_resync3(a0, a1, a2):
    n0 = {}
    n1 = {}
    n2 = {}
    for i, l in enumerate(a0):
        nm = qname_base_of(l)
        if nm not in n0:
            n0[nm] = i
    for j, l in enumerate(a1):
        nm = qname_base_of(l)
        if nm not in n1:
            n1[nm] = j
    for k, l in enumerate(a2):
        nm = qname_base_of(l)
        if nm not in n2:
            n2[nm] = k
    common = set(n0.keys()) & set(n1.keys()) & set(n2.keys())
    if not common:
        return None
    best = None
    best_key = None
    for nm in common:
        i = n0[nm]
        j = n1[nm]
        k = n2[nm]
        key = (max(i, j, k), i + j + k)
        if best is None or key < best_key:
            best = (i, j, k, nm)
            best_key = key
    return best  # (i, j, k, name)


def print_context(prev_pairs, out):
    if not prev_pairs:
        return
    print("-- context (last matched) --", file=out)
    for idx, (l1, l2) in enumerate(prev_pairs[-len(prev_pairs):]):
        f1, _, _ = parse_sam_line(l1)
        f2, _, _ = parse_sam_line(l2)
        print(f"{idx:2d}: {f1[0]} == {f2[0]}", file=out)


def print_context3(prev_triples, out):
    if not prev_triples:
        return
    print("-- context (last matched) --", file=out)
    for idx, (l0, l1, l2) in enumerate(prev_triples[-len(prev_triples):]):
        f0, _, _ = parse_sam_line(l0)
        f1, _, _ = parse_sam_line(l1)
        f2, _, _ = parse_sam_line(l2)
        print(f"{idx:2d}: {f0[0]} == {f1[0]} == {f2[0]}", file=out)


def collect_group_lines(pb, first_line, qn):
    lines = []
    cur = first_line
    while cur is not None and qname_base_of(cur) == qn:
        lines.append(cur)
        cur = pb.read()
    return lines, cur


def seek_and_collect_group_lines(pb, current_line, qn):
    cur = current_line
    while cur is not None and qname_base_of(cur) != qn:
        cur = pb.read()
    if cur is None:
        return [], None, False
    lines, after = collect_group_lines(pb, cur, qn)
    return lines, after, True


def backlog_suffix_for_qname(backlog, qn):
    if not backlog:
        return []
    bl = list(backlog)
    i = len(bl) - 1
    while i >= 0 and qname_base_of(bl[i]) == qn:
        i -= 1
    return bl[i+1:]


def backlog_all_for_qname(backlog, qn):
    if not backlog:
        return []
    return [l for l in backlog if qname_base_of(l) == qn]


def canonicalize_lines(lines):
    return sorted(canonical_norm(l) for l in lines)


essential_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8]  # up to TLEN; keep SEQ/QUAL full in raw print


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam_original", help="Original BAM")
    ap.add_argument("bam_python", help="BAM from Python version")
    ap.add_argument("bam_c", help="BAM from C version")
    ap.add_argument("--max-diffs", type=int, default=10)
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--max-window", type=int, default=100)
    ap.add_argument("--backlog", type=int, default=20000, help="Number of prior lines to retain per stream to reconstruct non-contiguous groups")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--context", type=int, default=3)
    ap.add_argument("--pyc-only", action="store_true")
    args = ap.parse_args()

    p0 = open_sam_stream(args.bam_original, threads=args.threads)
    p1 = open_sam_stream(args.bam_python, threads=args.threads)
    p2 = open_sam_stream(args.bam_c, threads=args.threads)

    it0 = sam_alignment_iter(p0)
    it1 = sam_alignment_iter(p1)
    it2 = sam_alignment_iter(p2)
    pb0 = Peekable(iter(it0))
    pb1 = Peekable(iter(it1))
    pb2 = Peekable(iter(it2))

    prev_matched = deque(maxlen=args.context)
    back0 = deque(maxlen=256)  # matched-only context
    back1 = deque(maxlen=256)
    back2 = deque(maxlen=256)
    # full backlogs (matched or not), up to max_window
    back_all0 = deque(maxlen=args.backlog)
    back_all1 = deque(maxlen=args.backlog)
    back_all2 = deque(maxlen=args.backlog)
    diffs = 0

    line0 = pb0.read()
    line1 = pb1.read()
    line2 = pb2.read()

    while line0 is not None and line1 is not None and line2 is not None and diffs < args.max_diffs:
        f0, t0, c0 = parse_sam_line(line0)
        f1, t1, c1 = parse_sam_line(line1)
        f2, t2, c2 = parse_sam_line(line2)
        # record current lines in full backlogs
        back_all0.append(line0)
        back_all1.append(line1)
        back_all2.append(line2)
        n0 = canonical_norm(line0)
        n1 = canonical_norm(line1)
        n2 = canonical_norm(line2)
        eq = (n1 == n2) if args.pyc_only else ((n0 == n1) and (n1 == n2))
        if eq:
            prev_matched.append((line0, line1, line2))
            back0.append(line0)
            back1.append(line1)
            back2.append(line2)
            line0 = pb0.read()
            line1 = pb1.read()
            line2 = pb2.read()
            continue

        diff_id = diffs + 1
        diff_context = list(prev_matched)
        initial_line0 = line0
        initial_line1 = line1
        initial_line2 = line2

        # resync across streams
        base_win = max(args.window, 1)
        found3 = None
        found2 = None
        ahead0 = [line0]
        ahead1 = [line1]
        ahead2 = [line2]
        w = base_win
        while (found3 is None and found2 is None) and w <= max(args.max_window, base_win):
            pb0.ensure(w - 1)
            pb1.ensure(w - 1)
            pb2.ensure(w - 1)
            ahead0 = [line0] + list(pb0.queue)[: w - 1]
            ahead1 = [line1] + list(pb1.queue)[: w - 1]
            ahead2 = [line2] + list(pb2.queue)[: w - 1]
            if args.pyc_only:
                found3 = find_resync3(ahead0, ahead1, ahead2)
                if found3 is None:
                    found2 = find_resync(ahead1, ahead2)
            else:
                found3 = find_resync3(ahead0, ahead1, ahead2)
            if (found3 is None and found2 is None):
                w *= 2
        if (found3 is None and found2 is None):
            if args.pyc_only:
                print("-- resync: no common qname found within window (python/c focus) --")
            else:
                print("-- resync: no common qname found within window across original/python/c --")
            break
        if found3 is not None:
            i0, i1, i2, nm = found3
            print(f"-- resync: matched qname {nm} at lookahead original={i0}, python={i1}, c={i2}")
            triple = True
        else:
            i1, i2, nm = found2
            print(f"-- resync: matched qname {nm} at lookahead python={i1}, c={i2}")
            triple = False

        # skip to match positions
        if triple and i0 > 0:
            pb0.ensure(i0)
            for l in list(pb0.queue)[:i0]:
                back_all0.append(l)
            pb0.skip(i0)
            line0 = pb0.read()
        if i1 > 0:
            pb1.ensure(i1)
            for l in list(pb1.queue)[:i1]:
                back_all1.append(l)
            pb1.skip(i1)
            line1 = pb1.read()
        if i2 > 0:
            pb2.ensure(i2)
            for l in list(pb2.queue)[:i2]:
                back_all2.append(l)
            pb2.skip(i2)
            line2 = pb2.read()

        if (triple and (line0 is None)) or line1 is None or line2 is None:
            break

        # collect and print full readgroup for this qname
        g0 = []
        # collect any earlier non-contiguous lines for this qname from backlog
        p0_all = backlog_all_for_qname(back_all0, nm)
        if triple:
            f0, line0 = collect_group_lines(pb0, line0, nm)
            e0 = pb0.drain_for_qname(nm, args.max_window)
            # merge unique lines preserving order preference: backlog-all, contiguous, drained
            seen = set()
            for l in (p0_all + f0 + e0):
                if l not in seen:
                    g0.append(l); seen.add(l)
        else:
            f0, line0, found0 = seek_and_collect_group_lines(pb0, line0, nm)
            e0 = pb0.drain_for_qname(nm, args.max_window)
            seen = set()
            for l in (p0_all + f0 + e0):
                if l not in seen:
                    g0.append(l); seen.add(l)

        p1_all = backlog_all_for_qname(back_all1, nm)
        f1, line1 = collect_group_lines(pb1, line1, nm)
        e1 = pb1.drain_for_qname(nm, args.max_window)
        g1 = []
        seen = set()
        for l in (p1_all + f1 + e1):
            if l not in seen:
                g1.append(l); seen.add(l)

        p2_all = backlog_all_for_qname(back_all2, nm)
        f2, line2 = collect_group_lines(pb2, line2, nm)
        e2 = pb2.drain_for_qname(nm, args.max_window)
        g2 = []
        seen = set()
        for l in (p2_all + f2 + e2):
            if l not in seen:
                g2.append(l); seen.add(l)

        resolved = False
        c1_list = canonicalize_lines(g1)
        c2_list = canonicalize_lines(g2)
        s1 = Counter(c1_list)
        s2 = Counter(c2_list)
        if args.pyc_only:
            eq_pc = (s1 == s2)
            if not eq_pc:
                # Aggressively scan further ahead for additional lines for this qname
                extra1 = pb1.drain_for_qname(nm, args.backlog)
                extra2 = pb2.drain_for_qname(nm, args.backlog)
                if extra1 or extra2:
                    # Merge unique
                    seen = set(g1)
                    for l in extra1:
                        if l not in seen:
                            g1.append(l); seen.add(l)
                    seen = set(g2)
                    for l in extra2:
                        if l not in seen:
                            g2.append(l); seen.add(l)
                    c1_list = canonicalize_lines(g1)
                    c2_list = canonicalize_lines(g2)
                    s1 = Counter(c1_list)
                    s2 = Counter(c2_list)
                eq_pc = (s1 == s2)
            if eq_pc:
                resolved = True
            else:
                diffs = diff_id
                print(f"=== Difference #{diffs}")
                print_context3(diff_context, sys.stdout)
                print("-- original --")
                print(initial_line0)
                print("--  python  --")
                print(initial_line1)
                print("--     c    --")
                print(initial_line2)
                print(f"-- readgroup: {nm} --")
                if g0:
                    print(f"--- original ({len(g0)}) ---")
                    for l in g0:
                        print(l)
                else:
                    print(f"--- original (not found after scanning) ---")
                print(f"--- python   ({len(g1)}) ---")
                for l in g1:
                    print(l)
                print(f"--- c        ({len(g2)}) ---")
                for l in g2:
                    print(l)
                print(f"--- order-insensitive canonical match (python vs c): {eq_pc} ---")
                only1 = []
                only2 = []
                keys = set(s1.keys()) | set(s2.keys())
                for k in sorted(keys):
                    c1c = s1.get(k, 0)
                    c2c = s2.get(k, 0)
                    if c1c > c2c:
                        only1.extend([k] * (c1c - c2c))
                    if c2c > c1c:
                        only2.extend([k] * (c2c - c1c))
                if only1:
                    print("--- only in python   ---")
                    for l in only1:
                        print(l)
                if only2:
                    print("--- only in c        ---")
                    for l in only2:
                        print(l)
        else:
            c0_list = canonicalize_lines(g0)
            s0 = Counter(c0_list)
            eq_all = (s0 == s1 == s2)
            if eq_all:
                resolved = True
            else:
                diffs = diff_id
                print(f"=== Difference #{diffs}")
                print_context3(diff_context, sys.stdout)
                print("-- original --")
                print(initial_line0)
                print("--  python  --")
                print(initial_line1)
                print("--     c    --")
                print(initial_line2)
                print(f"-- readgroup: {nm} --")
                if g0:
                    print(f"--- original ({len(g0)}) ---")
                    for l in g0:
                        print(l)
                else:
                    print(f"--- original (not found after scanning) ---")
                print(f"--- python   ({len(g1)}) ---")
                for l in g1:
                    print(l)
                print(f"--- c        ({len(g2)}) ---")
                for l in g2:
                    print(l)
                print(f"--- order-insensitive canonical match across all three: {eq_all} ---")
                only0 = []
                only1 = []
                only2 = []
                keys = set(s0.keys()) | set(s1.keys()) | set(s2.keys())
                for k in sorted(keys):
                    c0c = s0.get(k, 0)
                    c1c = s1.get(k, 0)
                    c2c = s2.get(k, 0)
                    m = max(c1c, c2c)
                    if c0c > m:
                        only0.extend([k] * (c0c - m))
                    m = max(c0c, c2c)
                    if c1c > m:
                        only1.extend([k] * (c1c - m))
                    m = max(c0c, c1c)
                    if c2c > m:
                        only2.extend([k] * (c2c - m))
                if only0:
                    print("--- only in original ---")
                    for l in only0:
                        print(l)
                if only1:
                    print("--- only in python   ---")
                    for l in only1:
                        print(l)
                if only2:
                    print("--- only in c        ---")
                    for l in only2:
                        print(l)

        if resolved:
            # treat regrouped lines as matched context for future diffs
            limit = min(len(g1), len(g2))
            for idx in range(limit):
                o_line = g0[idx] if idx < len(g0) else g1[idx]
                prev_matched.append((o_line, g1[idx], g2[idx]))
                back0.append(o_line)
                back1.append(g1[idx])
                back2.append(g2[idx])
            continue

    try:
        p0.stdout.close()
        p1.stdout.close()
        p2.stdout.close()
    except Exception:
        pass
    try:
        p0.terminate()
        p1.terminate()
        p2.terminate()
    except Exception:
        pass


if __name__ == "__main__":
    main()
