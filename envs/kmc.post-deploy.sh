#!/usr/bin/env bash
set -e



mkdir -p ${CONDA_PREFIX}/software

git clone --recurse-submodules https://github.com/refresh-bio/KMC ${CONDA_PREFIX}/software/kmc
cd ${CONDA_PREFIX}/software/kmc

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
 
@@ -169,7 +169,7 @@ $(LIB_KMC_CORE): $(KMC_CORE_OBJS) $(RADULS_OBJS) $(KMC_API_OBJS) $(KFF_OBJS)
 	-mkdir -p $(OUT_BIN_DIR)
 	ar rcs $@ $^
 
-kmc: $(KMC_CLI_OBJS) $(LIB_KMC_CORE) $(LIB_ZLIB)
+kmc: $(RADULS_OBJS) $(KMC_CLI_OBJS) $(KMC_CORE_OBJS) $(KMC_API_OBJS) $(KFF_OBJS) $(LIB_ZLIB)
 	-mkdir -p $(OUT_BIN_DIR)
 	$(CC) $(CLINK) -o $(OUT_BIN_DIR)/$@ $^' > kmc_make.patch


patch < kmc_make.patch
make -j32
make

cp ${CONDA_PREFIX}/software/kmc/bin/* ${CONDA_PREFIX}/bin


