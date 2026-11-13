// [[Rcpp::plugins(cpp17)]]
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include "utils.hpp"
// Only for Illumina

static void close_fastq_input(FILE*& plain, gzFile& gzip) {
    if (plain) closeFile(plain);
    if (gzip) closeGzip(gzip);
}

//' Reorder R2 FASTQ records to match R1
//'
//' Reorder the records in a paired-end R2 FASTQ file so that their identifiers
//' occur in the same order as the corresponding records in R1. Both files
//' must contain the same reads; only their order may differ.
//'
//' @param fastq_1 Character string giving the path to the R1 FASTQ file. Its
//'   record order is used as the reference. Plain and gzip-compressed files are
//'   supported.
//' @param fastq_2 Character string giving the path to the R2 FASTQ file that
//'   will be reordered. Plain and gzip-compressed files are supported.
//' @param output_fastq Character string giving the output path for the reordered
//'   R2 records. The output must not overwrite either input file.
//' @param compress_output Logical value indicating whether the output should be
//'   gzip-compressed. A file name ending in `.gz` also enables compression.
//' @param compress_level Integer from 0 to 9 giving the gzip compression level.
//'   This argument is ignored for uncompressed output.
//'
//' @details
//' Only the first field of each FASTQ identifier is compared. Text following a
//' space or tab is ignored, as are conventional `/1` and `/2` mate suffixes.
//' For example, `@read42/1 1:N:0:1` and `@read42/2 2:N:0:1` are treated as a
//' pair.
//'
//' R2 records are indexed in memory and then written in R1 order. Memory use
//' therefore increases with the size of the R2 file. Both files must contain
//' exactly the same unique identifiers. An error is raised for invalid FASTQ
//' records, duplicate or empty identifiers, missing mates, extra R2 records,
//' unreadable files, or an output path that refers to an input file.
//'
//' The sequence, quality string, and complete original R2 header are preserved;
//' only the order of the R2 records changes. R1 is never modified.
//'
//' @return Invisibly returns `NULL`. Reordered R2 records are written to
//'   `output_fastq`.
//'
//' @examples
//' \donttest{
//' r1 <- system.file("extdata", "trnL_40k_1.fq.gz", package = "METAeDNA")
//' r2 <- system.file("extdata", "trnL_40k_2.fq.gz", package = "METAeDNA")
//' output <- tempfile(fileext = ".fastq.gz")
//' reorder_fastq(r1, r2, output)
//' validate_paired_fastq(r1, output)
//' }
//'
//' @export
// [[Rcpp::export]]
void reorder_fastq(std::string fastq_1, std::string fastq_2,
                   std::string output_fastq,
                   bool compress_output = false,
                   int compress_level = 6) {
    check_processing_parameters(1, 1, compress_level);
    check_output_path(output_fastq, {fastq_1, fastq_2});
    fastq_1 = expandTilde(fastq_1);
    fastq_2 = expandTilde(fastq_2);
    output_fastq = expandTilde(output_fastq);

    if (fastq_1 == fastq_2) {
        CERR << "[Error] fastq_1 and fastq_2 must be different files." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    if (output_fastq == fastq_1 || output_fastq == fastq_2) {
        CERR << "[Error] output_fastq must not overwrite an input file." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    if (compress_level < 0 || compress_level > 9) {
        CERR << "[Error] compress_level must be between 0 and 9." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }

    const bool gz1 = isGzFile(fastq_1);
    const bool gz2 = isGzFile(fastq_2);
    FILE* p1 = gz1 ? nullptr : fopen(fastq_1.c_str(), "rb");
    FILE* p2 = gz2 ? nullptr : fopen(fastq_2.c_str(), "rb");
    gzFile z1 = gz1 ? gzopen(fastq_1.c_str(), "rb") : nullptr;
    gzFile z2 = gz2 ? gzopen(fastq_2.c_str(), "rb") : nullptr;
    FileGuard guard1(p1, z1), guard2(p2, z2);
    if ((gz1 ? z1 == nullptr : p1 == nullptr) ||
        (gz2 ? z2 == nullptr : p2 == nullptr)) {
        close_fastq_input(p1, z1);
        close_fastq_input(p2, z2);
        CERR << "[Error] Cannot open one or both FASTQ input files." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }

    std::unordered_map<std::string, FastqRecord> r2_by_id;
    FastqRecord record;
    while (readRec(gz2, p2, z2, record, true)) {
        if (record.id.empty() || record.id[0] != '@' ||
            record.plus.empty() || record.plus[0] != '+' ||
            record.seq.size() != record.qual.size()) {
            close_fastq_input(p1, z1);
            close_fastq_input(p2, z2);
            CERR << "[Error] Invalid FASTQ record in fastq_2." << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
        std::string id(paired_read_id(record.id));
        if (id.empty() || !r2_by_id.emplace(id, std::move(record)).second) {
            close_fastq_input(p1, z1);
            close_fastq_input(p2, z2);
            CERR << "[Error] Duplicate or empty read identifier in fastq_2: "
                 << id << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
    }
    close_fastq_input(p2, z2);
    p2 = nullptr;
    z2 = nullptr;

    std::string lower_output = output_fastq;
    std::transform(lower_output.begin(), lower_output.end(), lower_output.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    const bool gz_out = compress_output ||
        (lower_output.size() >= 3 && lower_output.substr(lower_output.size() - 3) == ".gz");
    FILE* out_plain = gz_out ? nullptr : fopen(output_fastq.c_str(), "wb");
    std::string mode = "wb" + std::to_string(compress_level);
    gzFile out_gz = gz_out ? gzopen(output_fastq.c_str(), mode.c_str()) : nullptr;
    FileGuard output_guard(out_plain, out_gz);
    if (gz_out ? out_gz == nullptr : out_plain == nullptr) {
        close_fastq_input(p1, z1);
        CERR << "[Error] Cannot open output FASTQ file: " << output_fastq
             << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }

    std::unordered_map<std::string, bool> seen_r1;
    size_t written = 0;
    bool failed = false;
    while (readRec(gz1, p1, z1, record, true)) {
            if (record.id.empty() || record.id[0] != '@' ||
                record.plus.empty() || record.plus[0] != '+' ||
                record.seq.size() != record.qual.size())
            {
                CERR << "[Error] Invalid FASTQ record in fastq_1." << std::endl;
                failed = true;
                break;
            }
            std::string id(paired_read_id(record.id));
            if (id.empty() || !seen_r1.emplace(id, true).second)
            {
                CERR << "[Error] Duplicate or empty read identifier in fastq_1: "
                     << id << std::endl;
                failed = true;
                break;
            }
            auto match = r2_by_id.find(id);
            if (match == r2_by_id.end())
            {
                CERR << "[Error] R2 is missing read identifier: " << id
                     << std::endl;
                failed = true;
                break;
            }
            const FastqRecord& r2 = match->second;
            writeBuf(r2.id + "\n" + r2.seq + "\n" + r2.plus + "\n" + r2.qual + "\n",
                     gz_out, out_plain, out_gz);
            r2_by_id.erase(match);
            ++written;
            
#ifdef USE_RCPP_OUTPUT
        if (written % 100000 == 0) {
            Rcpp::checkUserInterrupt();
        }
#endif
    }
    if (!failed && !r2_by_id.empty()) failed = true;

    close_fastq_input(p1, z1);
    if (out_plain) closeFile(out_plain);
    if (out_gz) closeGzip(out_gz);
    if (failed) throw std::runtime_error("R1 and R2 identifiers or record counts differ; output is incomplete.");
}
