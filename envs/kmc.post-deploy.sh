#!/usr/bin/env bash
set -e



mkdir -p ${CONDA_PREFIX}/software

KMC_COMMIT=751ef36a3c1ccc6dda664f529ad218dc51d76f55
KMC_BUILD_JOBS=${KMC_BUILD_JOBS:-32}
git clone --recurse-submodules https://github.com/refresh-bio/KMC ${CONDA_PREFIX}/software/kmc
cd ${CONDA_PREFIX}/software/kmc
git checkout --detach ${KMC_COMMIT}
git submodule update --init --recursive

echo '--- Makefile
+++ Makefile.2
@@ -1,4 +1,4 @@
-all: kmc kmc_dump kmc_tools py_kmc_api
+all: kmc kmc_dump kmc_tools
 
 dummy := $(shell git submodule update --init --recursive)
 
@@ -62,8 +62,8 @@ else
 		STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread	
 	else
 		CPU_FLAGS = -m64
-		STATIC_CFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
-		STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
+		STATIC_CFLAGS = -lpthread 
+		STATIC_LFLAGS = -lpthread
 	endif
 	PY_FLAGS = -fPIC
 endif
@@ -151,11 +151,11 @@ $(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(KMC_DUMP_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(
 	$(CC) $(CFLAGS) -I 3rd_party/cloudflare -c $< -o $@
 
 $(KMC_MAIN_DIR)/raduls_sse2.o: $(KMC_MAIN_DIR)/raduls_sse2.cpp
-	$(CC) $(CFLAGS) -msse2 -c $< -o $@
+	$(CC) $(CFLAGS) -msse2 -mno-sse4 -mno-avx -mno-avx2 -c $< -o $@
 $(KMC_MAIN_DIR)/raduls_sse41.o: $(KMC_MAIN_DIR)/raduls_sse41.cpp
-	$(CC) $(CFLAGS) -msse4.1 -c $< -o $@
+	$(CC) $(CFLAGS) -msse4.1 -mno-avx -mno-avx2 -c $< -o $@
 $(KMC_MAIN_DIR)/raduls_avx.o: $(KMC_MAIN_DIR)/raduls_avx.cpp
-	$(CC) $(CFLAGS) -mavx -c $< -o $@
+	$(CC) $(CFLAGS) -mavx -mno-avx2 -c $< -o $@
 $(KMC_MAIN_DIR)/raduls_avx2.o: $(KMC_MAIN_DIR)/raduls_avx2.cpp
 	$(CC) $(CFLAGS) -mavx2 -c $< -o $@
 
@@ -172 +172 @@
-kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE) $(LIB_ZLIB)
+kmc: $(RADULS_OBJS) $(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(LIB_ZLIB)' > kmc_make.patch

echo '--- kmc_core/kmc.h
+++ kmc_core/kmc.h
@@ -164,0 +165,3 @@
+		// The total Stage 1 thread budget is also used to initialize memory in
+		// Stage 2. Keep it in sync when readers and splitters are set explicitly.
+		Params.n_threads = Params.n_readers + Params.n_splitters;' >> kmc_make.patch


patch --batch -p0 < kmc_make.patch
make -j${KMC_BUILD_JOBS}
make

cp ${CONDA_PREFIX}/software/kmc/bin/* ${CONDA_PREFIX}/bin
