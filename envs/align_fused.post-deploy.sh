#!/usr/bin/env bash
set -euo pipefail

# The Bioconda dragmap 1.3.0 release binary is unsafe: its release CPPFLAGS define NDEBUG,
# while upstream performs posix_memalign() inside BOOST_ASSERT().  The call is
# consequently optimized away and DRAGMAP segfaults as soon as AVX2
# Smith-Waterman is used.
readonly DRAGMAP_REPOSITORY='https://github.com/populationgenomics/DRAGMAP.git'
readonly DRAGMAP_REVISION='4f98e00e2aedc85e27ea6c118cf7b16663036c14'
readonly EXPECTED_VERSION='1.3.0-tokenizer-next-fix'
readonly BUILD_JOBS=8

: "${CONDA_PREFIX:?CONDA_PREFIX is required by the DRAGMAP post-deploy script}"

readonly SYSTEM_CC="${CONDA_PREFIX}/bin/gcc"
readonly SYSTEM_CXX="${CONDA_PREFIX}/bin/g++"

conda install -y conda-forge::gcc conda-forge::gxx conda-forge::git conda-forge::coreutils conda-forge::grep conda-forge::make conda-forge::binutils
  
[[ -x "${SYSTEM_CC}" ]] || { echo "Missing compiler: ${SYSTEM_CC}" >&2; exit 1; }
[[ -x "${SYSTEM_CXX}" ]] || { echo "Missing compiler: ${SYSTEM_CXX}" >&2; exit 1; }

software_dir="${CONDA_PREFIX}/software"
source_dir="${software_dir}/DRAGMAP"
binary="${CONDA_PREFIX}/bin/dragen-os"

mkdir -p "${software_dir}"
if [[ -e "${source_dir}" && ! -d "${source_dir}/.git" ]]; then
    echo "Refusing to overwrite non-git path: ${source_dir}" >&2
    exit 1
fi

if [[ ! -d "${source_dir}/.git" ]]; then
    git clone --filter=blob:none --no-checkout \
        "${DRAGMAP_REPOSITORY}" "${source_dir}"
fi

git -C "${source_dir}" fetch --depth 1 origin "${DRAGMAP_REVISION}"
git -C "${source_dir}" checkout --detach "${DRAGMAP_REVISION}"

actual_revision=$(git -C "${source_dir}" rev-parse HEAD)
if [[ "${actual_revision}" != "${DRAGMAP_REVISION}" ]]; then
    echo "Unexpected DRAGMAP revision: ${actual_revision}" >&2
    exit 1
fi

# GCC 14 correctly detects that a 20-byte buffer cannot hold all decimal
# uint64_t values plus the terminating NUL.  Older compilers accepted this
# upstream code, but the project promotes warnings to errors.  Widen only that
# diagnostic buffer; this is outside the alignment path and also removes the
# one-byte overflow instead of suppressing the compiler warning.
hash_table_source="${source_dir}/thirdparty/dragen/src/common/hash_generation/hash_table.c"
if grep -Fq 'binStr[20], valStr[20], pctStr[20]' "${hash_table_source}"; then
    git -C "${source_dir}" apply <<'PATCH'
diff --git a/thirdparty/dragen/src/common/hash_generation/hash_table.c b/thirdparty/dragen/src/common/hash_generation/hash_table.c
--- a/thirdparty/dragen/src/common/hash_generation/hash_table.c
+++ b/thirdparty/dragen/src/common/hash_generation/hash_table.c
@@ -107,3 +107,3 @@ void printHistogram(FILE* file, uint64_t* hist, int bins, int indent, int lastPlus)
 {
-  char     binStr[20], valStr[20], pctStr[20], binLine[120] = "", valLine[120] = "", pctLine[120] = "";
+  char     binStr[20], valStr[32], pctStr[20], binLine[120] = "", valLine[120] = "", pctLine[120] = "";
   int      i, x, len = 0, first = 1, max, last = 0;
PATCH
elif ! grep -Fq 'binStr[20], valStr[32], pctStr[20]' "${hash_table_source}"; then
    echo "Unexpected DRAGMAP histogram buffer declaration" >&2
    exit 1
fi


git -C "${source_dir}" apply <<'PATCH'
diff --git a/thirdparty/dragen/src/common/hash_generation/gen_hash_table.c b/thirdparty/dragen/src/common/hash_generation/gen_hash_table.c
index cdca3df..5a55699 100644
--- a/thirdparty/dragen/src/common/hash_generation/gen_hash_table.c
+++ b/thirdparty/dragen/src/common/hash_generation/gen_hash_table.c
@@ -249,7 +249,7 @@ void setDefaultHashParams(hashTableConfig_t* defConfig, const char* destDir, Has
     free(dir);
   }
 
-  defConfig->hostVersion = (char*)getHostVersion(0);
+  defConfig->hostVersion = (char*)getHostVersion();
 }
 
 //-------------------------------------------------------------------------------swhitmore
PATCH


(
    cd "${source_dir}"

    # conda compiler activation injects CPPFLAGS=-DNDEBUG.  DRAGMAP 1.3.0 has
    # allocation calls inside BOOST_ASSERT(), so retaining NDEBUG produces a
    # deterministic NULL write in ssw_init_avx2().  Use the system compiler
    # flags, matching the known-good build used by this pipeline since 2024.
    unset CC CXX CPPFLAGS CFLAGS CXXFLAGS LDFLAGS
    unset BOOST_ROOT BOOST_INCLUDEDIR BOOST_LIBRARYDIR
    CC="${SYSTEM_CC} -Wno-unused-but-set-variable" CXX="${SYSTEM_CXX} -include cstdint -Wno-unused-but-set-variable -Wno-nonnull" \
        HAS_GTEST=0 make -j "${BUILD_JOBS}"
)

install -m 0755 "${source_dir}/build/release/dragen-os" "${binary}"

actual_version=$("${binary}" --version)
if [[ "${actual_version}" != "${EXPECTED_VERSION}" ]]; then
    echo "Unexpected DRAGMAP version: ${actual_version}" >&2
    exit 1
fi

dynamic_symbols=$(nm -D "${binary}")
if [[ "${dynamic_symbols}" != *' posix_memalign@'* ]]; then
    echo "Unsafe DRAGMAP build: posix_memalign is absent from ${binary}" >&2
    exit 1
fi

printf 'repository=%s\nrevision=%s\nlocal_patch=%s\nversion=%s\n' \
    "${DRAGMAP_REPOSITORY}" "${actual_revision}" \
    'hash-table-uint64-buffer-20-to-32' "${actual_version}" \
    > "${CONDA_PREFIX}/dragen-os.build-info"

#remove depencies only needed at compile time. This does not benefit Docker, but apptainer containers do, since they are not layered
conda remove -y gcc gxx git coreutils grep
