#!/usr/bin/env python3
"""Inspect or change this job's node-local ZSlurm resource lease.

This is deliberately a small, pipeline-local client. It uses only Python's
standard library and the socket, token, and job id injected by zslurm_chief.
The manager/chief remains the authority for validation, memory safety floors,
idempotent releases, and admission control.
"""

import argparse
import json
import math
import os
import socket
import sys
import uuid


PROTOCOL_VERSION = 1
MAX_REQUEST_BYTES = 64 * 1024

ENV_SOCKET = "ZSLURM_LEASE_SOCKET"
ENV_TOKEN = "ZSLURM_LEASE_TOKEN"
ENV_JOB_ID = "ZSLURM_JOB_ID"


def finite_float(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def send_request(socket_path, request, timeout_s=30.0):
    timeout_s = max(0.1, finite_float(timeout_s, "timeout_s"))
    payload = json.dumps(request, sort_keys=True).encode("utf-8") + b"\n"
    if len(payload) > MAX_REQUEST_BYTES:
        raise ValueError("request is too large")

    client = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        client.settimeout(timeout_s)
        client.connect(socket_path)
        client.sendall(payload)

        chunks = []
        total = 0
        while True:
            chunk = client.recv(4096)
            if not chunk:
                break
            chunks.append(chunk)
            total += len(chunk)
            if total > MAX_REQUEST_BYTES:
                raise RuntimeError("lease response is too large")
            if b"\n" in chunk:
                break
        if not chunks:
            raise RuntimeError("local zslurm chief returned an empty response")
        response = b"".join(chunks).split(b"\n", 1)[0]
        return json.loads(response.decode("utf-8"))
    finally:
        client.close()


def request_from_environment(action, socket_timeout_s=30.0, **fields):
    socket_path = os.environ.get(ENV_SOCKET)
    token = os.environ.get(ENV_TOKEN)
    job_id = os.environ.get(ENV_JOB_ID)
    if not socket_path or not token or not job_id:
        raise RuntimeError(
            "this process has no zslurm lease environment; expected "
            f"{ENV_SOCKET}, {ENV_TOKEN}, and {ENV_JOB_ID}"
        )
    request = {
        "version": PROTOCOL_VERSION,
        "request_id": str(uuid.uuid4()),
        "action": action,
        "job_id": job_id,
        "token": token,
    }
    request.update(fields)
    return send_request(socket_path, request, timeout_s=socket_timeout_s)


def parser():
    result = argparse.ArgumentParser(
        prog="zslurm_lease_client.py",
        description="Inspect or change this job's local ZSlurm CPU/memory lease.",
    )
    result.add_argument(
        "--json", action="store_true", help="print the complete JSON response"
    )
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="show the current and maximum lease")

    set_parser = commands.add_parser("set", help="set an absolute lease target")
    set_parser.add_argument("--cores", type=float, help="target held CPU cores")
    set_memory = set_parser.add_mutually_exclusive_group()
    set_memory.add_argument("--mem-mb", type=float, help="target held memory in MB")
    set_memory.add_argument("--mem-gb", type=float, help="target held memory in GiB")
    set_parser.add_argument(
        "--wait",
        type=float,
        default=3600.0,
        metavar="SECONDS",
        help="maximum time to wait for growth (default: 3600)",
    )
    set_parser.add_argument(
        "--phase",
        help="name the phase started by this request, even if resources are unchanged",
    )

    phase_parser = commands.add_parser(
        "phase", help="start a named phase without changing held resources"
    )
    phase_parser.add_argument("name", help="semantic phase name")

    release_parser = commands.add_parser(
        "release", help="atomically release resources relative to the current holding"
    )
    release_parser.add_argument("--cores", type=float, help="CPU cores to release")
    release_memory = release_parser.add_mutually_exclusive_group()
    release_memory.add_argument("--mem-mb", type=float, help="memory MB to release")
    release_memory.add_argument("--mem-gb", type=float, help="memory GiB to release")
    release_parser.add_argument(
        "--release-id",
        help="stable logical completion id; retries with this id are idempotent",
    )
    return result


def print_human(response):
    if response.get("ok"):
        if "released_cores" in response:
            print(
                "{status}: released {released_cores:g} cores and "
                "{released_mem_mb:g} MB; holding {held_cores:g}/{max_cores:g} "
                "cores and {held_mem_mb:g}/{max_mem_mb:g} MB "
                "(epoch {epoch})".format(**response)
            )
            return
        print(
            "{status}: {held_cores:g}/{max_cores:g} cores, "
            "{held_mem_mb:g}/{max_mem_mb:g} MB (epoch {epoch})".format(**response)
        )
        safety = response.get("safety") or {}
        if safety.get("adjusted"):
            print(
                "safety floor applied: "
                f"{safety.get('cpu_floor', 0):g} cores, "
                f"{safety.get('memory_floor_mb', 0):g} MB memory "
                f"(observed {safety.get('observed_memory_mb', 0):g} MB)",
                file=sys.stderr,
            )
    else:
        print(
            f"{response.get('status', 'error')}: "
            f"{response.get('message', 'lease request failed')}",
            file=sys.stderr,
        )


def main(argv=None):
    argument_parser = parser()
    args = argument_parser.parse_args(argv)
    try:
        if args.command == "status":
            response = request_from_environment("status", socket_timeout_s=10.0)
        elif args.command == "set":
            mem_mb = args.mem_mb
            if args.mem_gb is not None:
                mem_mb = args.mem_gb * 1024.0
            if args.cores is None and mem_mb is None and args.phase is None:
                argument_parser.error(
                    "set requires --cores, --mem-mb, --mem-gb, or --phase"
                )
            response = request_from_environment(
                "set",
                socket_timeout_s=max(10.0, args.wait + 5.0),
                cores=args.cores,
                mem_mb=mem_mb,
                timeout_s=args.wait,
                phase=args.phase,
            )
        elif args.command == "phase":
            response = request_from_environment(
                "phase", socket_timeout_s=10.0, phase=args.name
            )
        else:
            mem_mb = args.mem_mb
            if args.mem_gb is not None:
                mem_mb = args.mem_gb * 1024.0
            if args.cores is None and mem_mb is None:
                argument_parser.error(
                    "release requires --cores, --mem-mb, or --mem-gb"
                )
            response = request_from_environment(
                "release",
                socket_timeout_s=10.0,
                cores=args.cores or 0.0,
                mem_mb=mem_mb or 0.0,
                release_id=args.release_id,
            )
    except Exception as exc:
        print(f"zslurm_lease_client.py: {exc}", file=sys.stderr)
        return 5

    if args.json:
        print(json.dumps(response, sort_keys=True))
    else:
        print_human(response)
    return int(response.get("code", 0 if response.get("ok") else 5))


if __name__ == "__main__":
    sys.exit(main())
