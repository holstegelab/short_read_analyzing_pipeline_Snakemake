#!/usr/bin/env bash
set -euo pipefail

GATK_ROOT="${GATK_CNV_ROOT:-/gpfs/work3/0/qtholstg/hg38_res_v2/software/gatk_4.4}"

can_import_gcnvkernel() {
    python -c 'import gcnvkernel' >/dev/null 2>&1
}

install_gatk_python_package() {
    local package_path
    local package_candidates=()

    for package_path in \
        "${GATK_ROOT}/build/gatkPythonPackageArchive.zip" \
        "${GATK_ROOT}/build/distributions/gatkPythonPackageArchive.zip"
    do
        if [[ -e "${package_path}" ]]; then
            package_candidates+=("${package_path}")
        fi
    done

    if [[ -d "${GATK_ROOT}/build" ]]; then
        while IFS= read -r package_path; do
            package_candidates+=("${package_path}")
        done < <(
            find "${GATK_ROOT}/build" -maxdepth 5 -type f \
                \( -iname "*python*package*.zip" -o -iname "*PythonPackageArchive*.zip" \) \
                -print 2>/dev/null | sort
        )
    fi

    for package_path in "${package_candidates[@]}"; do
        if python -m pip install --no-deps "${package_path}" && can_import_gcnvkernel; then
            return 0
        fi
    done

    return 1
}

add_gatk_source_to_pythonpath() {
    local gcnv_init
    local gcnv_parent
    local site_packages

    gcnv_init=$(
        find "${GATK_ROOT}/src/main/python" -path "*/gcnvkernel/__init__.py" \
            -print -quit 2>/dev/null || true
    )
    if [[ -z "${gcnv_init}" ]]; then
        return 1
    fi

    gcnv_parent="$(dirname "$(dirname "${gcnv_init}")")"
    site_packages=$(
        python - <<'PY'
import site
import sysconfig

paths = site.getsitepackages()
print(paths[0] if paths else sysconfig.get_paths()["purelib"])
PY
    )
    printf '%s\n' "${gcnv_parent}" > "${site_packages}/gatk_gcnvkernel.pth"
}

if ! can_import_gcnvkernel; then
    if ! install_gatk_python_package; then
        add_gatk_source_to_pythonpath || true
    fi
fi

if ! can_import_gcnvkernel; then
    echo "Could not make gcnvkernel importable from ${GATK_ROOT}." >&2
    echo "Set GATK_CNV_ROOT to the local GATK 4.4 checkout before creating this env if the path differs." >&2
    exit 1
fi
