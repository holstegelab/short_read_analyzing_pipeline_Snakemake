import ast
import csv
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def load_read_sexchrom():
    # Importing common.py loads the samplefiles of the working directory.
    source = (REPO / "common.py").read_text()
    node = next(
        n for n in ast.parse(source).body
        if isinstance(n, ast.FunctionDef) and n.name == "read_sexchrom"
    )
    namespace = {"csv": csv}
    exec(compile(ast.Module([node], []), "common.py", "exec"), namespace)
    return namespace["read_sexchrom"]


def test_read_sexchrom_selects_kmer_or_reported_column(tmp_path):
    table = tmp_path / "s.sex_chrom.tab"
    table.write_text("sample\treported\tsex\tyyratio\nA\tM\tF\t0.0\nB\tF\tF\t0.0\n")
    read_sexchrom = load_read_sexchrom()
    assert read_sexchrom(str(table)) == {"A": "F", "B": "F"}
    assert read_sexchrom(str(table), column="reported") == {"A": "M", "B": "F"}


def test_glnexus_inputs_use_the_configured_sex_column():
    source = (REPO / "GLnexus.smk").read_text()
    assert source.count("read_sexchrom(") == 4
    assert source.count("'.sex_chrom.tab'), column=sex_column)") == 4
    assert "config.get(\"glnexus_sex\", \"kmer\") == \"reported\"" in source
