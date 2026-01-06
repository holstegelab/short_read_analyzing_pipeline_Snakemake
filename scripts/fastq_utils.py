import sys
import gzip
import os
import subprocess
import threading
import time
      
class PairedFastQReader(threading.Thread):
    def __init__(self, filea,fileb):
        self.filea = filea
        self.fileb = fileb
        self.buffer = {}
        self.splitchar = None
        self.buffer_size = 1000

        self.file_finished = False
        self.error = None
        self.continue_reading = threading.Event()
        self.search_read = None
        self.read_found_or_buffer_full = threading.Event()

        threading.Thread.__init__(self)


    def _raise_if_error(self):
        if self.error is not None:
            raise self.error


    def _start_pigz(self, path):
        proc = subprocess.Popen(
            ['pigz', '-dc', '--', path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            errors='replace',
            bufsize=1,
        )
        stderr_buf = []

        def _drain():
            try:
                for line in proc.stderr:
                    stderr_buf.append(line)
                    if len(stderr_buf) > 2000:
                        stderr_buf.pop(0)
            except Exception:
                pass

        t = threading.Thread(target=_drain, daemon=True)
        t.start()
        return proc, stderr_buf, t


    def _read_fastq_record(self, handle, path, recno):
        h = handle.readline()
        if not h:
            return None
        s = handle.readline()
        p = handle.readline()
        q = handle.readline()
        if not s or not p or not q:
            raise RuntimeError(f"FASTQ appears truncated/corrupted while reading {path} at record {recno}")
        h = h.strip()
        s = s.strip()
        p = p.strip()
        q = q.strip()
        if not h.startswith('@'):
            raise RuntimeError(f"FASTQ parsing error in {path}: expected header starting with '@' but got: {h[:200]}")
        if not p.startswith('+'):
            raise RuntimeError(f"FASTQ parsing error in {path}: expected '+' line but got: {p[:200]}")
        return (h, s, q)


    def _fastq_name(self, header):
        token = header.split()[0]
        if token.startswith('@'):
            token = token[1:]
        if token.endswith('/1') or token.endswith('/2'):
            token = token[:-2]
        return token


    def run(self):
        proc1 = None
        proc2 = None
        t1 = None
        t2 = None
        stderr1 = None
        stderr2 = None
        try:
            proc1, stderr1, t1 = self._start_pigz(self.filea)
            proc2, stderr2, t2 = self._start_pigz(self.fileb)
            f1 = proc1.stdout
            f2 = proc2.stdout
            if f1 is None or f2 is None:
                raise RuntimeError('Failed to open pigz stdout stream')

            recno = 0
            while True:
                recno += 1
                r1 = self._read_fastq_record(f1, self.filea, recno)
                r2 = self._read_fastq_record(f2, self.fileb, recno)
                if r1 is None and r2 is None:
                    break
                if r1 is None or r2 is None:
                    raise RuntimeError('FASTQ input files ended at different positions')
                header1, seq1, qual1 = r1
                header2, seq2, qual2 = r2

                name1 = self._fastq_name(header1)
                name2 = self._fastq_name(header2)
                if name1 != name2:
                    raise RuntimeError('FastQ input files not in same read order')

                self.buffer[name1] = (name1, seq1, seq2, qual1, qual2)

                if name1 == self.search_read or self.search_read is True:
                    self.read_found_or_buffer_full.set()

                if len(self.buffer) >= self.buffer_size:
                    self.read_found_or_buffer_full.set()
                    self.continue_reading.clear()
                    self.continue_reading.wait()

            rc1 = proc1.wait()
            rc2 = proc2.wait()
            if t1 is not None:
                t1.join()
            if t2 is not None:
                t2.join()
            if rc1 != 0:
                err = ''.join(stderr1 or [])
                raise RuntimeError(f"pigz failed for {self.filea} (exit code {rc1}):\n{err}")
            if rc2 != 0:
                err = ''.join(stderr2 or [])
                raise RuntimeError(f"pigz failed for {self.fileb} (exit code {rc2}):\n{err}")

            self.read_found_or_buffer_full.set()
            self.file_finished = True
        except Exception as e:
            self.error = e
            self.file_finished = True
            self.read_found_or_buffer_full.set()
            self.continue_reading.set()
        finally:
            try:
                if proc1 is not None and proc1.poll() is None:
                    proc1.terminate()
            except Exception:
                pass
            try:
                if proc2 is not None and proc2.poll() is None:
                    proc2.terminate()
            except Exception:
                pass

    def popRead(self, name):
        self._raise_if_error()
        while not name in self.buffer:
            self._raise_if_error()
            self.search_read = name
            self.read_found_or_buffer_full.clear()
            while True:
                self.read_found_or_buffer_full.wait(timeout=5)
                if name in self.buffer:
                    break
                elif self.file_finished:
                    self._raise_if_error()
                    raise RuntimeError(f"Read {name} not found while file has ended. Is BAM file reordered?")
                elif len(self.buffer) >= self.buffer_size: #buffer full, read not found
                    self.buffer_size *=2
                    sys.stderr.write(f"Order of reads has changed between BAM and FASTQ. Growing look-back buffer to {self.buffer_size}\n")
                    sys.stderr.flush()
                    self.continue_reading.set()
                    break
                else:    
                    self.continue_reading.set()
                    #wait a bit more
                    continue

        read = self.buffer[name]
        del self.buffer[name]
        
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()
        return read

    def retrieveRead(self):
        self._raise_if_error()
        if len(self.buffer) < self.buffer_size:
            self.continue_reading.set()

        if not self.buffer:
            if self.file_finished:
                self._raise_if_error()
                raise StopIteration()
                
            self.read_found_or_buffer_full.clear()
            self.search_read = True
            self.read_found_or_buffer_full.wait()
            if not self.buffer:
                assert self.file_finished
                self._raise_if_error()
                raise StopIteration()

        return self.buffer.popitem()[1]


class PairedFastQReaderSimple:
    def __init__(self, filea,fileb):
        self.filea = filea
        self.fileb = fileb
        self.splitchar = None

        self.file_finished = False
        self.search_read = None


    def retrieveRead(self):
        proc1 = None
        proc2 = None
        t1 = None
        t2 = None
        stderr1 = []
        stderr2 = []

        def _start(path, stderr_buf):
            proc = subprocess.Popen(
                ['pigz', '-dc', '--', path],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding='utf-8',
                errors='replace',
                bufsize=1,
            )

            def _drain():
                try:
                    for line in proc.stderr:
                        stderr_buf.append(line)
                        if len(stderr_buf) > 2000:
                            stderr_buf.pop(0)
                except Exception:
                    pass

            t = threading.Thread(target=_drain, daemon=True)
            t.start()
            return proc, t

        def _read_record(handle, path, recno):
            h = handle.readline()
            if not h:
                return None
            s = handle.readline()
            p = handle.readline()
            q = handle.readline()
            if not s or not p or not q:
                raise RuntimeError(f"FASTQ appears truncated/corrupted while reading {path} at record {recno}")
            h = h.strip()
            s = s.strip()
            p = p.strip()
            q = q.strip()
            if not h.startswith('@'):
                raise RuntimeError(f"FASTQ parsing error in {path}: expected header starting with '@' but got: {h[:200]}")
            if not p.startswith('+'):
                raise RuntimeError(f"FASTQ parsing error in {path}: expected '+' line but got: {p[:200]}")
            return (h, s, q)

        def _name(header):
            token = header.split()[0]
            if token.startswith('@'):
                token = token[1:]
            if token.endswith('/1') or token.endswith('/2'):
                token = token[:-2]
            return token

        try:
            proc1, t1 = _start(self.filea, stderr1)
            proc2, t2 = _start(self.fileb, stderr2)
            f1 = proc1.stdout
            f2 = proc2.stdout
            if f1 is None or f2 is None:
                raise RuntimeError('Failed to open pigz stdout stream')

            recno = 0
            while True:
                recno += 1
                r1 = _read_record(f1, self.filea, recno)
                r2 = _read_record(f2, self.fileb, recno)
                if r1 is None and r2 is None:
                    break
                if r1 is None or r2 is None:
                    raise RuntimeError('FASTQ input files ended at different positions')
                header1, seq1, qual1 = r1
                header2, seq2, qual2 = r2
                name1 = _name(header1)
                name2 = _name(header2)
                if name1 != name2:
                    raise RuntimeError('FastQ input files not in same read order')
                yield (name1, seq1, seq2, qual1, qual2)
        finally:
            rc1 = None
            rc2 = None
            try:
                if proc1 is not None:
                    rc1 = proc1.wait()
            except Exception:
                pass
            try:
                if proc2 is not None:
                    rc2 = proc2.wait()
            except Exception:
                pass
            try:
                if t1 is not None:
                    t1.join()
            except Exception:
                pass
            try:
                if t2 is not None:
                    t2.join()
            except Exception:
                pass
            if rc1 not in (0, None):
                raise RuntimeError(f"pigz failed for {self.filea} (exit code {rc1}):\n{''.join(stderr1)}")
            if rc2 not in (0, None):
                raise RuntimeError(f"pigz failed for {self.fileb} (exit code {rc2}):\n{''.join(stderr2)}")


