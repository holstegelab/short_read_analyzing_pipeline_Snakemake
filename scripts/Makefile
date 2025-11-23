# Build the fastcheck C extension in-place and the bam_merge binary

PYTHON ?= python3
CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra
LDFLAGS ?= -lz
HTSLIB_PREFIX ?= $(CONDA_PREFIX)
HTSLIB_INCLUDE ?= $(HTSLIB_PREFIX)/include
HTSLIB_LIB ?= $(HTSLIB_PREFIX)/lib
ENV_NAME ?= preprocess
ENV_YAML ?= ../envs/preprocess.yaml

.PHONY: all all-hts build build-hts clean preprocess-env build-preprocess build-preprocess-hts bam_merge bam_dechimer fix_bam_rg_pairs

all: build bam_merge bam_dechimer fix_bam_rg_pairs

# Build both Python extensions and bam_merge
all-hts: build build-hts bam_merge fix_bam_rg_pairs

build:
	$(PYTHON) setup.py build_ext --inplace

# Build the htslib-based extension in-place
build-hts:
	$(PYTHON) setup_hts.py build_ext --inplace

# Ensure the conda environment defined in ../envs/preprocess.yaml exists (python=3.10)
preprocess-env:
	conda env update -n $(ENV_NAME) -f $(ENV_YAML) --prune || conda env create -n $(ENV_NAME) -f $(ENV_YAML)

# Build the extension using the preprocess conda environment (python 3.10 ABI)
build-preprocess: preprocess-env
	conda run -n $(ENV_NAME) $(PYTHON) setup.py build_ext --inplace

# Build the htslib-based extension using the preprocess conda environment
build-preprocess-hts: preprocess-env
	conda run -n $(ENV_NAME) $(PYTHON) setup_hts.py build_ext --inplace

bam_merge: bam_merge.c
	$(CC) $(CFLAGS) -o bam_merge bam_merge.c $(LDFLAGS)

bam_dechimer: bam_dechimer.c
	$(CC) $(CFLAGS) -o bam_dechimer bam_dechimer.c $(LDFLAGS)

# Build BAM/CRAM RG/flags fixer (requires htslib)
fix_bam_rg_pairs: fix_bam_rg_pairs.c
	$(CC) $(CFLAGS) -I$(HTSLIB_INCLUDE) -o fix_bam_rg_pairs fix_bam_rg_pairs.c -L$(HTSLIB_LIB) -Wl,-rpath,$(HTSLIB_LIB) -lhts $(LDFLAGS)

clean:
	rm -rf build/ __pycache__
	rm -f fastcheck.*.so fastcheck.*.pyd fastcheck_hts.*.so fastcheck_hts.*.pyd bam_merge bam_dechimer fix_bam_rg_pairs
