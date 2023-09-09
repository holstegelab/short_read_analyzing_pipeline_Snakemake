#!/usr/bin/env bash
set -e



mkdir -p ${CONDA_PREFIX}/software

git clone --recurse-submodules https://github.com/refresh-bio/KMC ${CONDA_PREFIX}/software/kmc
cd ${CONDA_PREFIX}/software/kmc

echo '--- Makefile	2023-09-06 12:06:55.211683000 +0200
+++ Makefile.2	2023-09-06 12:14:27.316297000 +0200
@@ -1,4 +1,4 @@
-all: kmc kmc_dump kmc_tools py_kmc_api
+all: kmc kmc_dump kmc_tools
 
 UNAME_S := $(shell uname -s)
 UNAME_M := $(shell uname -m)
@@ -68,7 +68,7 @@
 
 
 CFLAGS	= -Wall -O3 -fsigned-char $(CPU_FLAGS) $(STATIC_CFLAGS) -std=c++14
-CLINK	= -lm $(STATIC_LFLAGS) -O3 -std=c++14
+CLINK	= -lm -lpthread -O3 -std=c++14
 PY_KMC_API_CFLAGS = $(PY_FLAGS) -Wall -shared -std=c++14 -O3
 
 KMC_CLI_OBJS = \
' > kmc_make.patch

patch < kmc_make.patch
make -j32
make

cp ${CONDA_PREFIX}/software/kmc/bin/* ${CONDA_PREFIX}/bin


