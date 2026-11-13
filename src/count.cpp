#include <string>
#include <vector>
#include <chrono>
#include "utils.hpp"

// -------------------------------------------------------------------
// helper: extract all integers after patterns like:
//   size=123;
//   "total_count":123
//   "count":123
//
// returns sum of all detected values
// -------------------------------------------------------------------
static size_t extract_size_fields(const std::string& header) {
    size_t sum = 0;
    size_t position = 0, value = 0;
    while (header_count(header, "size=", value, position)) sum += value;
    return sum;
}

// -------------------------------------------------------------------
//' Count FASTA or FASTQ records and stored read abundances
//'
//' Count records and sum the \code{size=} fields in their headers. Missing
//' size fields contribute zero. JSON counts are not added to this sum.
//'
//' @param input_file (character) path to the input FASTA or FASTQ file. The
//'   file may be plain text or gzip-compressed (suffix \code{.gz}).
//'
//' @return
//' Invisibly returns \code{NULL}. The console output includes:
//' \itemize{
//'   \item total number of sequence records,
//'   \item sum of all \code{size=} values (if present),
//'   \item time elapsed during processing.
//' }
//' @export
// [[Rcpp::export]]
void count_records(std::string input_file) {
    show_time("Counting records and extracting size fields...");
    input_file = expandTilde(input_file);
    if (!std::filesystem::exists(input_file)) {
        throw std::runtime_error("Input file does not exist: " + input_file);
    }
    auto start_time = std::chrono::system_clock::now();

    bool is_fastq = isFastqFile(input_file);
    bool is_gz = isGzFile(input_file);

    FILE* fp_plain = nullptr;
    gzFile fp_gz = nullptr;
    FileGuard guard(fp_plain, fp_gz);

    if (is_gz) {
        fp_gz = gzopen(input_file.c_str(), "rb");
        if (!fp_gz) {
            throw std::runtime_error("Cannot open input file: " + input_file);
        }
    }
    else {
        fp_plain = fopen(input_file.c_str(), "rb");
        if (!fp_plain) {
            throw std::runtime_error("Cannot open input file: " + input_file);
        }
    }


    size_t record_count = 0;
    size_t total_size_from_headers = 0;

    FastqRecord rec;

    while (readRec(is_gz, fp_plain, fp_gz, rec, is_fastq)) {
        record_count++;

        // Extract size values from header (rec.id)
        total_size_from_headers += extract_size_fields(rec.id);
#ifdef USE_RCPP_OUTPUT
        if (record_count % 10000 == 0) Rcpp::checkUserInterrupt();
#endif
    }

    closeFile(fp_plain);
    closeGzip(fp_gz);

    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    show_time("Total records: " + std::to_string(record_count));
/*
    show_time("Sum of size fields: " + std::to_string(total_size_from_headers));
    show_time("", 1, elapsed.count());

*/
}
