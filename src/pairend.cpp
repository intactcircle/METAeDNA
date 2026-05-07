// [[Rcpp::plugins(cpp11)]]
//#include <Rcpp.h>
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>         // For FILE*, fopen, fgets, fread, fclose, fwrite
#include <cstdlib>        // For getenv, exit
#include <cstring>        // For strlen, memset, memcpy
#include <iostream>       // For CERR
#include <string>         // For std::string
#include <vector>         // For std::vector
#include <algorithm>      // For std::min, std::reverse, std::transform, std::max
#include <zlib.h>         // For gzFile, gzopen, gzgets, gzwrite, gzclose
#include <cmath>          // For pow
#include <array>          // For std::array
#include <chrono>
#include <ctime>
#include <iomanip>  // for put_time
#include <fstream>
#include <bitset>
//#include <emmintrin.h>  // SSE2
//#include <immintrin.h>  // popcnt intrinsics
#include "utils.hpp"
#ifdef _OPENMP
#include <omp.h> 
#else
#endif
#ifdef USE_RCPP_OUTPUT
#include <Rcpp.h>
using namespace Rcpp;
#else
#endif
// ------------------------------
// Rcpp Exported Function: pairend
// ------------------------------
//' Merge paired-end reads from FASTQ files
//'
//' This function merges paired-end reads from two FASTQ files by identifying
//' overlapping regions between read 1 and the reverse-complement of read 2.
//' Overlap length, mismatch allowance, and overlap scoring can be customized.
//'
//' @param fastq_1 (character) path to the FASTQ file containing read 1.
//'   The file may be plain text or gzip-compressed.
//' @param fastq_2 (character) path to the FASTQ file containing read 2.
//'   The file may be plain text or gzip-compressed.
//' @param out_file (character) path to the output FASTQ file containing merged
//'   paired-end reads. If \code{compress_output = TRUE} or the file name ends
//'   with \code{.gz}, the output is gzip-compressed.
//'
//' @param compress_level (integer, 0–9) gzip compression level passed to
//'   \code{gzopen()} when writing compressed output.
//'
//' @param min_overlap (integer) minimum overlap length required for merging
//'   two reads. Overlaps shorter than this threshold are not used.
//'
//' @param max_mismatches (integer) maximum number of mismatched bases allowed
//'   within the overlapping region between the two reads.
//'
//' @param compress_output (logical) if \code{TRUE}, write the output FASTQ file
//'   in gzip format regardless of the suffix of \code{out_file}.
//'
//'
//' @param lambda (double) penalty weight applied to mismatches during overlap
//'   scoring. This parameter controls the mismatch contribution to the overlap
//'   score returned by the algorithm.
//'
//' @param score_threshold (double) minimum overlap score required for accepting
//'   a merged read. Merged reads with scores below this threshold are discarded.
//'
//' @param n_threads (integer) number of parallel threads used for merging.
//'   When OpenMP is available, multi-threading accelerates overlap detection
//'   and merging within each processed batch.
//'
//' @param batch_size (integer) number of read pairs processed per batch.
//'   Larger batch sizes reduce I/O overhead but increase memory usage.
//'
//' @return
//' Invisibly returns \code{NULL}. The merged reads are written to
//' \code{out_file}, with each merged record containing updated sequence,
//' quality string, and a JSON block appended to the identifier summarizing
//' overlap size, mismatch count, and overlap score.
//' @export
// [[Rcpp::export]]
void pairend(std::string fastq_1, std::string fastq_2, std::string out_file,
    size_t compress_level = 6, size_t min_overlap = 10, size_t max_mismatches = 5,
    bool compress_output = false, bool enable_oes_test = true,
    double alpha = 1.0, double beta = -1.0,
    int n_threads = 4, int batch_size = 1000000)
{
    int max_threads = detect_system_threads();
    if (n_threads > max_threads) n_threads = max_threads;
    auto start_time = std::chrono::system_clock::now();
    

    if (fastq_1 == fastq_2) {
        CERR << "[Error] Check your input, The paired-end sequencing files cannot be the same."<< std::endl;
        return;
    }
    bool is_gz1 = isGzFile(fastq_1), is_gz2 = isGzFile(fastq_2);
    size_t record_count = 0, valid_count = 0;
    FILE* fp1 = nullptr; gzFile gz1 = nullptr;
    FILE* fp2 = nullptr; gzFile gz2 = nullptr;
    auto now = std::chrono::system_clock::now();
    show_time("Start streaming reads...");
    show_time("Using " + std::to_string(n_threads) + " thread(s)");
    //
	// Basic statistical data ------------------START
    if (is_gz1) gz1 = gzopen(expandTilde(fastq_1).c_str(), "rb");
    else        fp1 = fopen(expandTilde(fastq_1).c_str(), "rb");
    int i = 0;
    base_probability base_pb;
    std::vector<FastqRecord> rec_v;
    while (i < batch_size) {
        FastqRecord tempr;
        readRec(is_gz1, fp1, gz1, tempr, true);
        rec_v.push_back(tempr);
        i++;
    }
    base_pb = calculate_base_probabilities(rec_v);
    // Basic statistical data ------------------END
	//Close connection and reopen to reset file pointer
    if (is_gz1) gzclose(gz1);
    else        fclose(fp1);
    if (is_gz2) gzclose(gz2);
    else        fclose(fp2);
    if (is_gz1) gz1 = gzopen(expandTilde(fastq_1).c_str(), "rb");
    else        fp1 = fopen(expandTilde(fastq_1).c_str(), "rb");
    if (is_gz2) gz2 = gzopen(expandTilde(fastq_2).c_str(), "rb");
    else        fp2 = fopen(expandTilde(fastq_2).c_str(), "rb");

    std::string lower_out = expandTilde(out_file);
    out_file = expandTilde(out_file);
    std::transform(lower_out.begin(), lower_out.end(), lower_out.begin(), ::tolower);
    if (lower_out.size() >= 3 && lower_out.substr(lower_out.size() - 3) == ".gz")
        compress_output = true;
    FILE* outFp = nullptr; gzFile gzOut = nullptr;
    if (compress_output) {
        std::string mode = "wb" + std::to_string(compress_level);
        gzOut = gzopen(out_file.c_str(), mode.c_str());
    }
    else {
        outFp = fopen(out_file.c_str(), "w");
    }

    const unsigned buffer_size = 512 * 1024;
    char* buffer1 = new char[buffer_size];
    char* buffer2 = new char[buffer_size];

    // prepare batch containers
    const size_t BATCH = batch_size;
    struct Pair { FastqRecord r1, r2; };
    std::vector<Pair> batch;
    batch.reserve(BATCH);
    std::vector<std::string> outBuf;
    std::vector<bool> keep;

    //std::vector<std::string> merged_ids, merged_seqs, merged_plus, merged_quals;
    std::vector<double> ovl_score;
    std::vector<size_t> mismatch_num;
    std::vector<size_t> overlap_length;

#ifdef USE_RCPP_OUTPUT
    size_t intrrupt_count = 0;
#endif // check interrupt
    while (true) {
        
        batch.clear();
        for (size_t i = 0; i < BATCH; ++i) {
#ifdef USE_RCPP_OUTPUT
            intrrupt_count++;
            if (intrrupt_count % 100000 == 0) {
                Rcpp::checkUserInterrupt();
            }
#endif // check interrupt
            FastqRecord rec1, rec2;
            bool ok1 = readRec(is_gz1, fp1, gz1, rec1, true);
            bool ok2 = readRec(is_gz2, fp2, gz2, rec2, true);
            //----------------------------------------------------
            // seqs to upper -------------------------------------
            //----------------------------------------------------

            std::transform(rec1.seq.begin(), rec1.seq.end(), rec1.seq.begin(),
                [](unsigned char c) { return std::toupper(c); });
            std::transform(rec2.seq.begin(), rec2.seq.end(), rec2.seq.begin(),
                [](unsigned char c) { return std::toupper(c); });

            /*
* 
            bool ok1 = is_gz1
                ? readRecordGZ(gz1, rec1, buffer1, buffer_size,true)
                : readRecordPlain(fp1, rec1, buffer1, buffer_size, true);
            bool ok2 = is_gz2
                ? readRecordGZ(gz2, rec2, buffer2, buffer_size, true)
                : readRecordPlain(fp2, rec2, buffer2, buffer_size, true);
            */
            if (!ok1 || !ok2) break;
            batch.push_back({ rec1, rec2 });
        }
        if (batch.empty()) break;

        record_count += batch.size();
        
        //-----------------
            now = std::chrono::system_clock::now();
            std::ostringstream temp_string;
            temp_string << "Processed:" << record_count << " records(streaming)...";
            show_time(temp_string.str(), 2);
        //------------------



#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)

        {
				// local buffers to avoid write conflicts
            std::vector<std::string> local_buf;
            local_buf.reserve(batch.size() / n_threads + 100);
#else
            
            std::vector<std::string> local_buf;
			local_buf.reserve(batch.size() + 100);
#endif

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
			for (int i = 0; i < (int)batch.size(); ++i) {
				auto& rec1 = batch[i].r1;
				auto& rec2 = batch[i].r2;
				size_t mismatch_count = 0;
				size_t ovl = 0;
				double score = 0.0;
                double p_value = 0;
				FastqRecord merged_rec;
				find_best_overlap(rec1, rec2, base_pb, min_overlap,max_mismatches, mismatch_count, ovl, score, merged_rec, p_value, alpha, beta, enable_oes_test);

				std::ostringstream id_js;
				std::string clean_id = merged_rec.id;        // copy
				std::replace(clean_id.begin(), clean_id.end(), ' ', '_');   // replace spaces with underscores
				id_js << clean_id
					<< ";{\"pairend\":"
					<< "{\"mismatch\":" << mismatch_count
					<< ",\"overlap\":" << ovl
					<< ",\"ovl_score\":" << score
					<< "}}";
				std::string mergedID = id_js.str();
				std::ostringstream ss;
				ss << mergedID << '\n'
					<< merged_rec.seq << '\n'
					<< merged_rec.plus << '\n'
					<< merged_rec.qual << '\n';
				if (enable_oes_test ) {
					if (p_value <0.01&&mismatch_count <= max_mismatches && min_overlap <= ovl) {
                        local_buf.push_back(ss.str());
                    }
                }
                else {
					if (mismatch_count <= max_mismatches && min_overlap <= ovl) {
                        local_buf.push_back(ss.str());
                    }
                }
			}
            
			// write local buffer to output
#ifdef _OPENMP
#pragma omp critical
            {
#endif
                for (auto& s : local_buf) {
                    writeBuf(s, compress_output, outFp, gzOut);
                    ++valid_count;
                }
#ifdef _OPENMP
            }

} // end parallel
#endif


    }

    if (!is_gz1) fclose(fp1); else gzclose(gz1);
    if (!is_gz2) fclose(fp2); else gzclose(gz2);
    if (compress_output) gzclose(gzOut); else fclose(outFp);
    delete[] buffer1; delete[] buffer2;

    now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = now - start_time;
    std::ostringstream temp_string;
    COUT << std::endl;
    temp_string << "Processed: [Valid:" << valid_count << "] / [Total:" << record_count << "]";
    show_time(temp_string.str());
    show_time("Done.");
    show_time("",1, elapsed_seconds.count());

}
