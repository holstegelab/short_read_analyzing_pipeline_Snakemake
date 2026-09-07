from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parents[1]


def rule_block(text: str, rule_name: str) -> str:
    marker = f"rule {rule_name}:"
    start = text.index(marker)
    next_rule = text.find("\nrule ", start + len(marker))
    return text[start : next_rule if next_rule >= 0 else None]


def conda_dependencies(filename: str) -> set[str]:
    document = yaml.safe_load((REPO / "envs" / filename).read_text())
    dependencies = set()
    for dependency in document.get("dependencies", []):
        if not isinstance(dependency, str):
            continue
        dependencies.add(dependency.split("=", 1)[0].lower())
    return dependencies


def test_external_adapter_fusion_uses_complete_preprocess_environment():
    dependencies = conda_dependencies("preprocess.yaml")
    assert {"python", "samtools", "htslib", "adapterremoval", "pigz"} <= dependencies


def test_kmer_sex_fusion_has_runtime_and_analysis_modules():
    dependencies = conda_dependencies("kmc.yaml")
    # kmc and kmc_tools are installed by kmc.post-deploy.sh.
    assert {"python", "numpy", "pandas", "scipy", "pyyaml"} <= dependencies
    post_deploy = (REPO / "envs" / "kmc.post-deploy.sh").read_text()
    assert "KMC_BUILD_JOBS=${KMC_BUILD_JOBS:-8}" in post_deploy
    assert 'make -j"${KMC_BUILD_JOBS}"' in post_deploy
    assert "KMC_COMMIT=751ef36a3c1ccc6dda664f529ad218dc51d76f55" in post_deploy
    assert (
        "Params.n_threads = Params.n_readers + Params.n_splitters;"
        in post_deploy
    )
    assert "uint64 ReadPart(FILE* f, uchar* part" in post_deploy
    assert "max_read_attempts = 5" in post_deploy
    assert "if (!ferror(f))" in post_deploy
    assert "clearerr(f);" in post_deploy
    assert post_deploy.count("ReadPart(") == 4
    assert 'cp bin/* "${CONDA_PREFIX}/bin/"' in post_deploy


def test_kmc_reserves_average_cpu_without_reducing_tool_parallelism():
    aligner = (REPO / "Aligner.smk").read_text()
    runner = (REPO / "scripts" / "run_fused_kmer_sex.py").read_text()

    assert "KMC_RESERVED_CORES = 2" in aligner
    assert "rule kmer_reads:" not in aligner
    assert 'n="1.6"' in rule_block(aligner, "kmer_sex_fused")
    assert "use_threads=KMC_RESERVED_CORES" in rule_block(
        aligner, "kmer_sex_fused"
    )
    assert "--kmc-threads {resources.use_threads}" in rule_block(
        aligner, "kmer_sex_fused"
    )
    for option in ('"-sf12"', '"-sp12"', '"-sr1"'):
        assert option in runner


def test_alignment_fusion_combines_aligner_and_bam_tooling():
    dependencies = conda_dependencies("align_fused.yaml")
    assert {"python", "dragmap", "samtools", "htslib", "setuptools"} <= dependencies


def test_deepvariant_phasing_and_chrm_tail_share_complete_vcf_environment():
    dependencies = conda_dependencies("vcf_handling.yaml")
    assert {
        "python",
        "numpy",
        "cyvcf2",
        "whatshap",
        "bcftools",
        "samtools",
        "gatk4",
    } <= dependencies


def test_parallel_bam_qc_fusion_contains_every_qc_executable():
    dependencies = conda_dependencies("qc_fused.yaml")
    assert {
        "python",
        "numpy",
        "pypy",
        "pypy3.9",
        "samtools",
        "gatk4",
        "verifybamid2",
        "mosdepth",
    } <= dependencies


def test_production_fusions_are_unconditional_and_have_no_predecessors():
    stages = {
        "Aligner.smk": {
            "external_adapter_fused": ["external_alignments_to_fastq"],
            "align_reads_fused": ["align_reads", "merge_bam_alignment_dechimer", "sort_bam_alignment"],
            "kmer_sex_fused": ["kmer_reads", "get_validated_sex"],
        },
        "Deepvariant.smk": {"deepvariant_phasing_fused": ["deepvariant", "DVWhatshapPhasingMerge"]},
        "Stat.smk": {"bam_qc_fused": ["coverage", "verifybamid", "hs_stats", "artifacts_and_oxog_metrics", "samtools_stats", "bamstats_all_and_exome"]},
        "chrM_analysis.smk": {
            "chrm_extract_align_fused": ["extract_chrM_reads", "extract_NUMTs_reads", "align_chrM_and_NUMTs"],
            "chrm_mutect_tail_fused": ["mutect_calls_both", "merge_and_filter_both", "mutect_bp_resolution_both"],
        },
    }
    for filename, rules in stages.items():
        source = (REPO / filename).read_text()
        assert "FUSE_" not in source
        assert "ruleorder:" not in source
        for fused, predecessors in rules.items():
            assert f"\nrule {fused}:" in source
            for old in predecessors:
                assert f"rule {old}:" not in source


def test_adapter_routes_use_one_shared_implementation():
    source = (REPO / "Aligner.smk").read_text()
    assert "--identify-adapters" not in source
    assert "sample=NATIVE_FASTQ_SAMPLE_PATTERN" in rule_block(source, "adapter_removal")
    assert "sample=EXTERNAL_ALIGNMENT_SAMPLE_PATTERN" in rule_block(source, "external_adapter_fused")
    for runner in ("run_fastq_adapter.py", "run_fused_external_adapter.py"):
        assert "from adapter_processing import prepare_adapters" in (REPO / "scripts" / runner).read_text()


def test_runtime_helpers_are_not_owned_by_an_alignment_runner():
    for runner in (REPO / "scripts").glob("run_fused_*.py"):
        source = runner.read_text()
        assert "from run_fused_alignment import" not in source
        assert "from pipeline_runtime import" in source
