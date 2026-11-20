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
# A stats reader stops once the global sample quota is exhausted. The shared
# binary producer then stops ALL input streams, not just that reader's stream.
# Publish cancellation on every queue before completing any of them, so another
# gzip reader cannot mistake this intentional cutoff for a truncated file.
# Full counting uses fresh, uncancelled queues: real EOF/CRC failures stay fatal.
printf '%b\n' '--- kmc_core/queues.h
+++ kmc_core/queues.h
@@ -109,6 +109,11 @@
 \t\tstop = true;
 \t\tcv_push.notify_all();
 \t}
+\tbool is_cancelled()
+\t{
+\t\tlock_guard<mutex> lck(mtx);
+\t\treturn stop;
+\t}
 };
\040
\040
--- kmc_core/binary_reader.h
+++ kmc_core/binary_reader.h
@@ -504,6 +504,13 @@
 \t\t}
\040
+\t\t// A rejected push means the statistics sample quota was reached.
+\t\t// Mark EVERY stream cancelled before waking readers with completion.
+\t\t// Otherwise a different gzip reader sees a false, corrupt-file EOF.
+\t\tif (forced_to_finish)
+\t\t\tfor (auto& e : binary_pack_queues)
+\t\t\t\te->ignore_rest();
+
 \t\t//user may specify more fastq_readers than input files
 \t\tfor (auto& e : binary_pack_queues)
 \t\t\te->mark_completed();
 \t}
--- kmc_core/fastq_reader.cpp
+++ kmc_core/fastq_reader.cpp
@@ -1041,8 +1041,11 @@
-\t\t\t\t\tif (ret != Z_STREAM_END) {
+\t\t\t\t\t// Statistics sampling can intentionally stop the producer
+\t\t\t\t\t// mid-gzip-stream. Only explicit queue cancellation permits
+\t\t\t\t\t// this; physical truncation in sampling/counting stays fatal.
+\t\t\t\t\tif (ret != Z_STREAM_END && !binary_pack_queue->is_cancelled()) {
 \t\t\t\t\t\tstd::ostringstream ostr;
 \t\t\t\t\t\tostr << "Unexpected end of gzip file";
 \t\t\t\t\t\tCCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
 \t\t\t\t\t}
\040
 \t\t\t\t\treturn ret_val;
 \t\t\t\t}' >> kmc_make.patch

patch --batch --fuzz=0 --ignore-whitespace -p0 < kmc_make.patch
make -j"${KMC_BUILD_JOBS}"

cp bin/* "${CONDA_PREFIX}/bin/"

conda remove -y git make patch gcc gxx binutils
