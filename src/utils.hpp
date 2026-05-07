#pragma once
#include <cstdio>
#include <zlib.h>
#include <cstdlib>        // For getenv, exit
#include <cstring>        // For strlen, memset, memcpy
#include <iostream>       // For CERR
#include <string>         // For std::string
#include <vector>         // For std::vector
#include <algorithm>      // For std::min, std::reverse, std::transform, std::max
#include <cmath>          // For pow lgamma
#include <array>          // For std::array
#include <chrono>
#include <ctime>
#include <iomanip>  // for put_time
#include <fstream>
#include <bitset>
#include <limits>
//#include <emmintrin.h>  // SSE2
//#include <immintrin.h>  // popcnt intrinsics (可选)
//-----------------------------------------------------------
// Unified output stream selector.
// If compiled with -DUSE_RCPP_OUTPUT or defined before include,
// output will use Rcpp::Rcout / Rcpp::Rcerr.
// Otherwise defaults to COUT / CERR.
//-----------------------------------------------------------
#ifdef _OPENMP
#include <omp.h> 
#else
#endif
#ifdef USE_RCPP_OUTPUT
#include <Rcpp.h>
using namespace Rcpp;
#define COUT Rcpp::Rcout
#define CERR Rcpp::Rcerr
#else
#define COUT std::cout
#define CERR std::cerr
#endif
#if (defined(_WIN32) || defined(_WIN64) || defined(__linux__)) && \
    (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86))
#define METAEDNA_X86 1
#else
#define METAEDNA_X86 0
#endif
#if METAEDNA_X86
#include <immintrin.h>
#endif

struct base_probability
{
    double A_sqr, C_sqr, G_sqr, T_sqr, N_sqr;
    double A, C, G, T, N;
    int max_length;
};

// bool isFastqFile1(const std::string& filename); deprecated....

//using namespace Rcpp;
struct FastqRecord {
    std::string id;    // Header line (record identifier)
    std::string seq;   // Nucleotide sequence
    std::string plus;  // Separator line (usually a '+' sign)
    std::string qual;  // Quality score string
};
bool isFastqFile(const std::string& filename);
std::string expandTilde(const std::string& path);
bool isGzFile(const std::string& file_path, unsigned int buf_size = 1024);
//bool readNextRecordPlain(FILE* fp, FastqRecord& rec, char* buffer, unsigned int buf_size);
//bool readNextRecordGz(gzFile fp, FastqRecord& rec, char* buffer, unsigned int buf_size);
bool readRecordPlain(FILE* fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq = false);
bool readRecordGZ(gzFile fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq = false);
void writeBuf(const std::string& buf, bool is_gz_out, FILE* out_plain, gzFile out_gz);
bool readRec(bool is_gz_in, FILE* fp_plain, gzFile fp_gz, FastqRecord& rec, bool fastq = false, int buffer_size = 1024*1024);

 const std::array<unsigned char, 256>& get_comp_table();
// fast_reverse_complement computes the reverse complement of the input sequence.
// It pre-allocates the destination string to the proper size and iterates over the input from end to start.
// For each character, it uses the lookup table to get the complement if available.
 void fast_reverse_complement(const std::string& src, std::string& dest);

// fast_reverse performs a simple string reversal, used for quality strings.
 void fast_reverse(const std::string& src, std::string& dest);
std::vector<double> build_quality_table();
extern const std::vector<double> QUALITY_WEIGHT_TABLE;
size_t find_overlap_with_mismatches(const std::string& s1,
    const std::string& s2,
    bool& is_offset,
    size_t min_overlap,
    size_t max_mismatches,
    size_t& mismatch_count,
    size_t mismatch_positions[1000]);
//
double calculate_overlap_score(const std::string& s1, const std::string& q1,
    const std::string& s2, const std::string& q2,
    size_t overlap_start, size_t overlap_length, double lambda = 0.5);

// Structure representing a single FASTQ record (only FASTQ format is handled).
// trim function


void trimEnd(char* s);
void show_time(const std::string& s, size_t mode = 0, double secs = 0);
int detect_system_threads();

base_probability calculate_base_probabilities(const std::vector <FastqRecord>& rec_v);
void find_best_overlap(
    const FastqRecord& r1,
    const FastqRecord& r2,
    const base_probability& base_pb,
    const size_t& min_overlap,
    const size_t max_mismatches,
    size_t& mismatch_count,
    size_t& overlap_length,
    double& overlap_score,
    FastqRecord& merged_rec,
    double& p_value,
    double alpha, double beta,
    bool enable_oes_test = false,
    bool MAP = true
);
inline bool  CPU_SUPPORTS_AVX2() {
#if METAEDNA_X86 && (defined(__GNUC__) || defined(__clang__))
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}

