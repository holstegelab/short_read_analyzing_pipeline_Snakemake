"""Compare outputs (including temp/ensure) and resources to the frozen baseline."""
import ast
import json
from pathlib import Path
import re
import textwrap


REPO = Path(__file__).resolve().parents[1]
STAGES = {
    "Aligner.smk": ["adapter_removal", "external_adapter_fused", "align_reads_fused",
                    "kmer_sex_fused", "split_alignments_by_readgroup"],
    "Deepvariant.smk": ["deepvariant_phasing_fused"],
    "Stat.smk": ["bam_qc_fused"],
    "chrM_analysis.smk": ["chrm_extract_align_fused", "chrm_mutect_tail_fused"],
}


def contracts(root):
    result = {}
    for filename, rules in STAGES.items():
        text = (root / filename).read_text()
        for rule in rules:
            match = re.search(r"^( *)rule " + rule + r":.*$", text, re.M)
            assert match, rule
            indent = len(match[1])
            lines = [match[0]]
            for line in text[match.end() + 1:].splitlines():
                if line.strip() and not line.lstrip().startswith("#") and len(line) - len(line.lstrip()) <= indent:
                    break
                lines.append(line)
            block = "\n".join(line[indent:] if line.startswith(" " * indent) else line for line in lines)
            for section in ("input", "output", "resources"):
                found = re.search(r"^    " + section + r":(.*?)(?=^    [a-z_]+:|\Z)", block, re.M | re.S)
                assert found, (rule, section)
                expression = "contract(" + textwrap.dedent(found[1]) + "\n)"
                result[f"{filename}:{rule}:{section}"] = ast.dump(ast.parse(expression), include_attributes=False)
    return result


def test_surviving_stage_input_output_and_resource_contracts_are_unchanged():
    expected = json.loads((REPO / "tests/fixtures/fused_rule_contracts.json").read_text())
    actual = contracts(REPO)
    assert actual == {key: expected[key] for key in actual}
