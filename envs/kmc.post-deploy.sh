#!/usr/bin/env bash
set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is required by the KMC post-deploy script}"

conda install -y conda-forge::git conda-forge::patch conda-forge::make conda-forge::binutils conda-forge::gxx==13.4.0 conda-forge::gcc==13.4.0
KMC_COMMIT=751ef36a3c1ccc6dda664f529ad218dc51d76f55
KMC_BUILD_JOBS=${KMC_BUILD_JOBS:-8}

# Build outside of CONDA_PREFIX. Snakemake considers an environment incomplete
# until this script finishes, so another Snakemake invocation can remove and
# recreate that shared prefix concurrently. Building inside it made active
# compilers lose their source/output directories. A unique local build tree also
# avoids putting the compilation's metadata load on GPFS.
KMC_BUILD_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/kmc-build.XXXXXX")
trap 'rm -rf -- "$KMC_BUILD_ROOT"' EXIT
KMC_SOURCE_DIR=${KMC_BUILD_ROOT}/kmc

git clone --recurse-submodules https://github.com/refresh-bio/KMC "$KMC_SOURCE_DIR"
cd "$KMC_SOURCE_DIR"
git checkout --detach "$KMC_COMMIT"
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
+		STATIC_CFLAGS = -lpthread -Wno-deprecated-declarations
+		STATIC_LFLAGS = -lpthread -Wno-deprecated-declarations
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

echo '--- kmc_core/binary_reader.h
+++ kmc_core/binary_reader.h
@@ -21,0 +22,3 @@
+#include <cerrno>
+#include <chrono>
+#include <cstring>
@@ -22,0 +26 @@
+#include <thread>
@@ -78,0 +83,39 @@
+	uint64 ReadPart(FILE* f, uchar* part, const string& file_name)
+	{
+		constexpr uint32 max_read_attempts = 5;
+		constexpr uint32 initial_retry_delay_ms = 100;
+
+		for (uint32 attempt = 1; attempt <= max_read_attempts; ++attempt)
+		{
+			errno = 0;
+			uint64 readed = fread(part, 1, part_size, f);
+			if (!ferror(f))
+				return readed;
+
+			const int read_errno = errno;
+			clearerr(f);
+
+			// fread may return valid bytes and set the error indicator at the
+			// same time. Process those bytes now; the next call resumes from the
+			// current file position with a cleared error indicator.
+			if (readed)
+				return readed;
+
+			if (attempt < max_read_attempts)
+			{
+				const uint32 retry_delay_ms = initial_retry_delay_ms << (attempt - 1);
+				std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
+				continue;
+			}
+
+			std::ostringstream ostr;
+			ostr << "Error while reading file: " << file_name
+				 << " after " << max_read_attempts << " attempts";
+			if (read_errno)
+				ostr << ": " << std::strerror(read_errno);
+			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
+		}
+
+		return 0;
+	}
+
@@ -417 +460 @@
-		vector<tuple<FILE*, CBinaryPackQueue*, CompressionType>> files;
+		vector<tuple<FILE*, CBinaryPackQueue*, CompressionType, string>> files;
@@ -430 +473 @@
-			files.push_back(make_tuple(f, q, mode));
+			files.push_back(make_tuple(f, q, mode, file_name));
@@ -432 +475 @@
-			uint64 readed = fread(part, 1, part_size, f);
+			uint64 readed = ReadPart(f, part, file_name);
@@ -453 +496 @@
-				uint64 readed = fread(part, 1, part_size, get<0>(f));
+				uint64 readed = ReadPart(get<0>(f), part, get<3>(f));
@@ -467,0 +511 @@
+						get<3>(f) = file_name;
@@ -469 +513 @@
-						readed = fread(part, 1, part_size, get<0>(f));
+						readed = ReadPart(get<0>(f), part, get<3>(f));' >> kmc_make.patch


patch --batch --ignore-whitespace -p0 < kmc_make.patch
make -j"${KMC_BUILD_JOBS}"

cp bin/* "${CONDA_PREFIX}/bin/"

conda remove -y git make patch gcc gxx binutils
