#pragma once
#include <cstdio>
#include <zlib.h>
#include <cstdlib>        // For getenv, exit
#include <cstring>        // For strlen, memset, memcpy
#include <iostream>       // For CERR
#include <string>         // For std::string
#include <string_view>
#include <charconv>
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
#include <climits>
#include <stdexcept>
#include <sstream>
#include <filesystem>
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

struct FastqRecord {
    std::string id;    // Header line (record identifier)
    std::string seq;   // Nucleotide sequence
    std::string plus;  // Separator line (usually a '+' sign)
    std::string qual;  // Quality score string
};

// Read size= values or quoted JSON keys followed by a colon and a number.
inline bool header_count(std::string_view header, std::string_view key,
                         size_t& value, size_t& position) {
    while ((position = header.find(key, position)) != std::string_view::npos) {
        position += key.size();
        if (!key.empty() && key.back() == '"') {
            position = header.find_first_not_of(" \t\r\n\f\v", position);
            if (position == std::string_view::npos) return false;
            if (header[position] != ':') continue;
            position = header.find_first_not_of(" \t\r\n\f\v", position + 1);
            if (position == std::string_view::npos) return false;
        }
        const char* first = header.data() + position;
        const auto parsed = std::from_chars(first, header.data() + header.size(), value);
        if (parsed.ec == std::errc::result_out_of_range)
            throw std::overflow_error("Header count exceeds the supported integer range.");
        if (parsed.ec == std::errc()) {
            position = static_cast<size_t>(parsed.ptr - header.data());
            return true;
        }
    }
    return false;
}

// Return the biological read identifier without allocating a new string.
// FASTA/FASTQ markers, Illumina fields after whitespace, and conventional
// paired-end suffixes (/1 and /2) are ignored.
inline std::string_view paired_read_id(const std::string& header) {
    size_t start = (!header.empty() &&
                    (header.front() == '@' || header.front() == '>')) ? 1 : 0;
    size_t end = header.find_first_of(" \t", start);
    if (end == std::string::npos) end = header.size();

    size_t length = end - start;
    if (length > 2 && header[start + length - 2] == '/' &&
        (header[start + length - 1] == '1' ||
         header[start + length - 1] == '2')) {
        length -= 2;
    }
    return std::string_view(header.data() + start, length);
}

bool isFastqFile(const std::string& filename);
std::string expandTilde(const std::string& path);
bool isGzFile(const std::string& file_path, unsigned int buf_size = 1024);
bool readRecordPlain(FILE* fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq = false);
bool readRecordGZ(gzFile fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq = false);
void writeBuf(const std::string& buf, bool is_gz_out, FILE* out_plain, gzFile out_gz);
void closeFile(FILE*& file);
void closeGzip(gzFile& file);
void check_output_path(const std::string& output, const std::vector<std::string>& inputs);
void check_processing_parameters(int threads, double batch_size, double compression);

// Close open handles during normal return, an I/O error, or an R interrupt.
struct FileGuard {
    FILE*& plain;
    gzFile& gzip;
    FileGuard(FILE*& p, gzFile& g) : plain(p), gzip(g) {}
    ~FileGuard() { if (plain) std::fclose(plain); if (gzip) gzclose(gzip); }
};
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
base_probability calculate_base_probabilities(
    const std::vector<FastqRecord>& rec_v1,
    const std::vector<FastqRecord>& rec_v2,
    bool reverse_complement_r2 = true
);
bool find_best_overlap(
    const FastqRecord& r1,
    const FastqRecord& r2,
    const base_probability& base_pb,
    size_t min_overlap,
    size_t max_mismatches,
    size_t& mismatch_count,
    size_t& overlap_length,
    double& overlap_score,
    FastqRecord& merged_rec,
    double& p_value,
    double alpha, double beta
);
inline bool  CPU_SUPPORTS_AVX2() {
#if METAEDNA_X86 && (defined(__GNUC__) || defined(__clang__))
    return __builtin_cpu_supports("avx2");
#else
    return false;
#endif
}
