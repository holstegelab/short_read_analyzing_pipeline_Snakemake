from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parents[1]


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
    assert aligner.count("n=str(KMC_RESERVED_CORES)") == 2
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


def test_all_production_fusions_are_enabled_by_default():
    defaults = {
        "Aligner.smk": (
            "config.get('fuse_external_adapter', True)",
            "config.get('fuse_kmer_sex', True)",
            "config.get('fuse_alignment_phases', True)",
        ),
        "Deepvariant.smk": (
            "config.get('fuse_deepvariant_phasing', True)",
        ),
        "Stat.smk": ("config.get('fuse_bam_qc', True)",),
        "chrM_analysis.smk": (
            "config.get('fuse_chrm_extract_align', True)",
            "config.get('fuse_chrm_mutect_tail', True)",
        ),
    }
    for filename, expected_fragments in defaults.items():
        source = (REPO / filename).read_text()
        for fragment in expected_fragments:
            assert fragment in source
