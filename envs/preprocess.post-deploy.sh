#!/usr/bin/env bash
set -euo pipefail

# Snakemake copies this file next to the hashed Conda environment before it is
# executed.  SHORT_READ_PIPELINE_ROOT is therefore exported by the root
# Snakefile instead of deriving the checkout from BASH_SOURCE.
: "${CONDA_PREFIX:?CONDA_PREFIX is required by the preprocess post-deploy script}"
# Detect whether the script is running in an Apptainer container
[ "$(dirname "$CONDA_PREFIX")" = "/conda-envs" ] && export SHORT_READ_PIPELINE_ROOT=/tmp/bind

: "${SHORT_READ_PIPELINE_ROOT:?Run environment creation through the pipeline Snakefile}"
conda install -y conda-forge::flock conda-forge::make conda-forge::coreutils conda-forge::gcc==13.4.0

readonly source_dir="${SHORT_READ_PIPELINE_ROOT}/scripts"
readonly lock_file="${source_dir}/.native-build.lock"

required_sources=(
    Makefile
    setup.py
    setup_hts.py
    fastcheck.c
    fastcheck_hts.c
    bam_merge.c
    bam_dechimer.c
    fix_bam_rg_pairs.c
)

# These checksums make changes to native sources explicit deployment inputs.
# Update the matching value whenever one of the sources changes; because this
# file participates in Snakemake's environment hash, that creates a new Conda
# environment and reruns this hook.
declare -A expected_sha256=(
    [Makefile]=c59a8063b021080bf6b9e06601017fade504d2f0c0a54f58b92e3867a4d8c5fa
    [setup.py]=938f74e798df312c6f915c15f6d8124538cd0db0acbf3365a8e1f53af6cbce19
    [setup_hts.py]=e89025be23766f07dd2f8d69aaf9af8c96eeecccf65e5e91480928a1eb836b6c
    [fastcheck.c]=d1ba52105c7980fe76946d98bb29dcdeb69f1178c023f21c7ebd9e05417c09f2
    [fastcheck_hts.c]=d758ce95e99936611caaeb75cb7aa4c45cb0d6ae60f1f74b6d358158cfd78e9e
    [bam_merge.c]=bf0f35ae473a20cc5eee00a2a806d4710dba57fa0dc98c0cf67012ee0c91713a
    [bam_dechimer.c]=fd063b98ea535aa01f291a36aa39ce02ebd11840e3441856ef3048fbbe4ab588
    [fix_bam_rg_pairs.c]=ac7e7b9c24146a9f7a4f8177fd347be1af37670f9563d78eb5865344c63c0379
)

for name in "${required_sources[@]}"; do
    path="${source_dir}/${name}"
    [[ -f "${path}" ]] || { echo "Missing native source: ${path}" >&2; exit 1; }
    actual=$(sha256sum "${path}" | awk '{print $1}')
    if [[ "${actual}" != "${expected_sha256[${name}]}" ]]; then
        echo "Native source checksum mismatch for ${path}" >&2
        echo "expected ${expected_sha256[${name}]}; found ${actual}" >&2
        exit 1
    fi
done

command -v flock >/dev/null || { echo "flock is required to serialize native builds" >&2; exit 1; }
exec 9>"${lock_file}"
flock 9

build_root=$(mktemp -d "${TMPDIR:-/tmp}/short-read-native.XXXXXX")
trap 'rm -rf -- "${build_root}"' EXIT
for name in "${required_sources[@]}"; do
    cp -- "${source_dir}/${name}" "${build_root}/${name}"
done

echo "[preprocess.post-deploy] Building native tools with $(${CONDA_PREFIX}/bin/python -V 2>&1)" >&2
make -C "${build_root}" all-hts \
    PYTHON="${CONDA_PREFIX}/bin/python" \
    HTSLIB_PREFIX="${CONDA_PREFIX}"

install_atomic() {
    local source=$1
    local destination=$2
    local mode=$3
    local temporary
    temporary=$(mktemp "${destination}.install.XXXXXX")
    install -m "${mode}" "${source}" "${temporary}"
    mv -f -- "${temporary}" "${destination}"
}

for name in bam_merge bam_dechimer fix_bam_rg_pairs; do
    install_atomic "${build_root}/${name}" "${source_dir}/${name}" 0755
done

extension_suffix=$("${CONDA_PREFIX}/bin/python" -c \
    'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
fastcheck_module="${build_root}/fastcheck${extension_suffix}"
fastcheck_hts_module="${build_root}/fastcheck_hts${extension_suffix}"
[[ -f "${fastcheck_module}" ]] || { echo "Missing ${fastcheck_module}" >&2; exit 1; }
[[ -f "${fastcheck_hts_module}" ]] || { echo "Missing ${fastcheck_hts_module}" >&2; exit 1; }
install_atomic "${fastcheck_module}" \
    "${source_dir}/$(basename "${fastcheck_module}")" 0755
install_atomic "${fastcheck_hts_module}" \
    "${source_dir}/$(basename "${fastcheck_hts_module}")" 0755

PYTHONPATH="${source_dir}" "${CONDA_PREFIX}/bin/python" -c \
    'import fastcheck, fastcheck_hts; print("native preprocessing imports OK")'
echo "[preprocess.post-deploy] Native tools installed in ${source_dir}" >&2

conda remove -y coreutils make flock gcc
