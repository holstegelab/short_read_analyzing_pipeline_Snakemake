import json
import os
import socket
import subprocess
import threading
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
CLIENT = REPO / "scripts" / "zslurm_lease_client.py"


def _serve(socket_path: Path, requests: list[dict], count: int):
    ready = threading.Event()

    def target():
        completed_release_ids = set()
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as server:
            server.bind(str(socket_path))
            server.listen()
            ready.set()
            for _ in range(count):
                connection, _ = server.accept()
                with connection:
                    raw = b""
                    while b"\n" not in raw:
                        raw += connection.recv(4096)
                    request = json.loads(raw.split(b"\n", 1)[0])
                    requests.append(request)
                    action = request["action"]
                    if action == "status":
                        response = {
                            "ok": True,
                            "code": 0,
                            "status": "current",
                            "held_cores": 12.0,
                            "max_cores": 12.0,
                            "held_mem_mb": 12000.0,
                            "max_mem_mb": 12000.0,
                            "epoch": 0,
                        }
                    elif action == "set":
                        response = {
                            "ok": True,
                            "code": 0,
                            "status": "updated",
                            "held_cores": request["cores"],
                            "max_cores": 12.0,
                            "held_mem_mb": request["mem_mb"],
                            "max_mem_mb": 12000.0,
                            "epoch": 1,
                        }
                    else:
                        release_id = request["release_id"]
                        duplicate = release_id in completed_release_ids
                        completed_release_ids.add(release_id)
                        response = {
                            "ok": True,
                            "code": 0,
                            "status": "duplicate" if duplicate else "released",
                            "released_cores": 0.0 if duplicate else request["cores"],
                            "released_mem_mb": 0.0 if duplicate else request["mem_mb"],
                            "held_cores": 6.0,
                            "max_cores": 12.0,
                            "held_mem_mb": 6000.0,
                            "max_mem_mb": 12000.0,
                            "epoch": 2,
                        }
                    connection.sendall(
                        json.dumps(response, sort_keys=True).encode() + b"\n"
                    )

    thread = threading.Thread(target=target, daemon=True)
    thread.start()
    assert ready.wait(timeout=5)
    return thread


def _run(client_env, *args):
    return subprocess.run(
        [str(CLIENT), "--json", *args],
        check=True,
        capture_output=True,
        text=True,
        env=client_env,
    )


def test_pipeline_local_client_status_set_and_idempotent_release(tmp_path):
    assert CLIENT.is_file()
    assert os.access(CLIENT, os.X_OK)

    requests = []
    socket_path = tmp_path / "lease.sock"
    server = _serve(socket_path, requests, count=4)
    environment = os.environ.copy()
    environment.update(
        {
            "ZSLURM_LEASE_SOCKET": str(socket_path),
            "ZSLURM_LEASE_TOKEN": "test-token",
            "ZSLURM_JOB_ID": "123",
        }
    )

    status = json.loads(_run(environment, "status").stdout)
    updated = json.loads(
        _run(
            environment,
            "set",
            "--cores",
            "2.5",
            "--mem-gb",
            "4",
            "--wait",
            "7",
        ).stdout
    )
    first_release = json.loads(
        _run(
            environment,
            "release",
            "--cores",
            "1.5",
            "--mem-mb",
            "1000",
            "--release-id",
            "bam-qc:S1:coverage",
        ).stdout
    )
    duplicate_release = json.loads(
        _run(
            environment,
            "release",
            "--cores",
            "1.5",
            "--mem-mb",
            "1000",
            "--release-id",
            "bam-qc:S1:coverage",
        ).stdout
    )
    server.join(timeout=5)

    assert status["status"] == "current"
    assert updated["held_cores"] == 2.5
    assert updated["held_mem_mb"] == 4096.0
    assert first_release["status"] == "released"
    assert duplicate_release["status"] == "duplicate"
    assert [request["action"] for request in requests] == [
        "status",
        "set",
        "release",
        "release",
    ]
    assert all(request["version"] == 1 for request in requests)
    assert all(request["job_id"] == "123" for request in requests)
    assert requests[1]["timeout_s"] == 7.0
    assert requests[2]["release_id"] == requests[3]["release_id"]


def test_pipeline_local_client_fails_cleanly_without_lease_environment():
    environment = os.environ.copy()
    for name in ("ZSLURM_LEASE_SOCKET", "ZSLURM_LEASE_TOKEN", "ZSLURM_JOB_ID"):
        environment.pop(name, None)
    result = subprocess.run(
        [str(CLIENT), "--json", "status"],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert result.returncode == 5
    assert "has no zslurm lease environment" in result.stderr
