import hashlib
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
HOOK = ROOT / "envs" / "preprocess.post-deploy.sh"


def test_native_source_checksums_match_post_deploy_manifest():
    content = HOOK.read_text(encoding="utf-8")
    expected = dict(
        re.findall(r"^    \[([^]]+)\]=([0-9a-f]{64})$", content, re.MULTILINE)
    )
    assert expected
    for name, digest in expected.items():
        source = ROOT / "scripts" / name
        assert source.is_file()
        assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
