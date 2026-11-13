// [[Rcpp::plugins(cpp17)]]
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

// ------------------------------
// Rcpp Exported Function: pairend
// ------------------------------
//' Merge paired-end reads from FASTQ files
//'
//' \code{pairend()} merges paired-end reads from two FASTQ files by identifying
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
//'   within the overlapping region between the two reads. This is an additional
//'   METAeDNA filter applied after the PEAR OES p-value test.
//'
//' @param compress_output (logical) if \code{TRUE}, write the output FASTQ file
//'   in gzip format regardless of the suffix of \code{out_file}.
//'
//'
//' @param alpha (double) weight assigned to evidence supporting a base match in
//'   overlap scoring.
//'
//' @param beta (double) weight assigned to evidence against a base match in
//'   overlap scoring.
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
             size_t compress_level = 6, size_t min_overlap = 15, size_t max_mismatches = 3,
             bool compress_output = false, double alpha = 1.0, double beta = -1.0,
             int n_threads = 4, int batch_size = 1000000)
{
    check_processing_parameters(n_threads, batch_size, compress_level);
    if (min_overlap < 1 || min_overlap > INT_MAX || max_mismatches > INT_MAX ||
        !std::isfinite(alpha) || !std::isfinite(beta) || alpha <= beta)
        throw std::invalid_argument("Require positive min_overlap, non-negative max_mismatches, and finite alpha > beta.");
    fastq_1 = expandTilde(fastq_1);
    fastq_2 = expandTilde(fastq_2);
    out_file = expandTilde(out_file);
    if (std::filesystem::weakly_canonical(fastq_1) == std::filesystem::weakly_canonical(fastq_2))
        throw std::invalid_argument("R1 and R2 must be different files.");
    check_output_path(out_file, {fastq_1, fastq_2});
    n_threads = std::min(n_threads, detect_system_threads());
    const auto start = std::chrono::steady_clock::now();
    const bool gz1 = isGzFile(fastq_1), gz2 = isGzFile(fastq_2);
    FILE *p1 = nullptr, *p2 = nullptr, *out = nullptr;
    gzFile z1 = nullptr, z2 = nullptr, zout = nullptr;
    FileGuard guard1(p1, z1), guard2(p2, z2), guard_out(out, zout);
    if (gz1)
        z1 = gzopen(fastq_1.c_str(), "rb");
    else
        p1 = fopen(fastq_1.c_str(), "rb");
    if (gz2)
        z2 = gzopen(fastq_2.c_str(), "rb");
    else
        p2 = fopen(fastq_2.c_str(), "rb");
    if ((gz1 ? !z1 : !p1) || (gz2 ? !z2 : !p2))
        throw std::runtime_error("Cannot open paired FASTQ inputs.");
    std::string lower = out_file;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    compress_output = compress_output || (lower.size() >= 3 && lower.substr(lower.size() - 3) == ".gz");
    std::vector<FastqRecord> r1, r2;
    size_t total = 0, retained = 0, mismatches = 0, overlap_bases = 0;
    bool seq_upper_sign = true;
    bool seq_checked = false;
    show_time("Merging paired reads using " + std::to_string(n_threads) + " thread(s)");
    if (compress_output)
        zout = gzopen(out_file.c_str(), ("wb" + std::to_string(compress_level)).c_str());
    else
        out = fopen(out_file.c_str(), "wb");
    if (compress_output ? !zout : !out)
        throw std::runtime_error("Cannot open merged FASTQ output.");
    for (;;)
    {
        r1.clear();
        r2.clear();
        for (int i = 0; i < batch_size; ++i)
        {
            FastqRecord a, b;
            bool ok1 = readRec(gz1, p1, z1, a, true);
            bool ok2 = readRec(gz2, p2, z2, b, true);
            if (ok1 != ok2)
                throw std::runtime_error("R1 and R2 FASTQ files have different record counts.");
            if (!ok1)
                break;
                /* Pending
                
            if (paired_read_id(a.id) != paired_read_id(b.id))
                throw std::runtime_error("R1 and R2 identifiers are out of sync: " + a.id + " / " + b.id);
                */
            if (!seq_checked) {
                seq_upper_sign =
                    std::any_of(a.seq.begin(), a.seq.end(),
                                [](unsigned char c) { return std::islower(c); }) ||
                    std::any_of(b.seq.begin(), b.seq.end(),
                                [](unsigned char c) { return std::islower(c); });
                seq_checked = true;
            }
            if (seq_upper_sign) {
                std::transform(a.seq.begin(), a.seq.end(), a.seq.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
                std::transform(b.seq.begin(), b.seq.end(), b.seq.begin(),
                               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            }
            r1.push_back(std::move(a));
            r2.push_back(std::move(b));
        }
        if (r1.empty())
            break;
        // Estimate the empirical background from the read pairs in the
        // current processing batch, with R2 counted in merged orientation.
        const base_probability probabilities =
            calculate_base_probabilities(r1, r2, true);
        std::vector<std::string> records(r1.size());
        std::vector<size_t> mismatch(r1.size()), overlaps(r1.size());
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static)
#endif
        for (int i = 0; i < static_cast<int>(r1.size()); ++i)
        {
            FastqRecord merged;
            double score = 0, p_value = 1;
            const bool accepted = find_best_overlap(
                r1[i], r2[i], probabilities, min_overlap, max_mismatches,
                mismatch[i], overlaps[i], score, merged, p_value, alpha, beta);
            if (accepted)
            {
                std::string id = r1[i].id;
                std::replace(id.begin(), id.end(), ' ', '_');

                std::string text;
                // make readable and fast
                text.reserve(
                    id.size() +
                    merged.seq.size() +
                    merged.qual.size() +
                    128);

                text += id;
                text += ";{\"pairend\":{\"mismatch\":";
                text += std::to_string(mismatch[i]);
                text += ",\"overlap\":";
                text += std::to_string(overlaps[i]);
                text += ",\"ovl_score\":";
                text += std::to_string(score);
                text += "}}\n";
                text += merged.seq;
                text += "\n+\n";
                text += merged.qual;
                text += '\n';

                records[i] = std::move(text);
                //                std::ostringstream text;
                //               text << merged.id << ";{\"pairend\":{\"mismatch\":" << mismatch[i]
                //                  << ",\"overlap\":" << overlaps[i] << ",\"ovl_score\":" << score
                //                << "}}\n"
                //              << merged.seq << "\n+\n"
                //            << merged.qual << '\n';
                //     records[i] = text.str();
            }
        }
        for (size_t i = 0; i < records.size(); ++i)
        {
            if (records[i].empty())
                continue;
            writeBuf(records[i], compress_output, out, zout);
            ++retained;
            mismatches += mismatch[i];
            overlap_bases += overlaps[i];
        }
        total += r1.size();
        show_time("Processed " + std::to_string(total) + " pairs", 2);
#ifdef USE_RCPP_OUTPUT
        Rcpp::checkUserInterrupt();
#endif
    }
    if (!total)
        throw std::runtime_error("Paired FASTQ inputs contain no records.");
    closeFile(out);
    closeGzip(zout);
    COUT << "\n[Total pairs]: " << total << " [Merged pairs]: " << retained
         << " [Overlap mismatches]: " << mismatches << " [Overlap bases]: " << overlap_bases
         << " [Mismatch rate]: " << (overlap_bases ? double(mismatches) / overlap_bases : 0) << std::endl;
    show_time("", 1, std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
}
