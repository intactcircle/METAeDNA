#include <algorithm>
#include <cctype>
#include <chrono>
#include <ctime>
#include <fstream>
#include <zlib.h>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
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
// 解析信息结构
struct UniqueInfo {
    std::string prefix;
    std::string pairend_json;
    size_t total_count = 0;
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>> exp_samples;
};

// ------------------------------
// Rcpp Exported Function: unique
// ------------------------------
//' Collapse duplicate sequences in a FASTA/FASTQ file
//'
//' This function identifies all identical sequences in an input FASTA or FASTQ
//' file and collapses them into unique sequence entries. For each unique
//' sequence, the function aggregates the total read count and the per-PCR
//' (experiment–sample) counts derived from the input header metadata. The
//' resulting dereplicated FASTA file contains one record per unique sequence,
//' with an updated `"stat"` JSON block summarizing the merged counts and the
//' sequence length. If paired-end metadata was present (in a `"pairend"`
//' block), it is preserved in the output.
//'
//' @param input_fasta (character) path to the input FASTA or FASTQ file. The
//'   file may be plain text or gzip-compressed (suffix \code{.gz}).
//'
//' @param output_fasta (character) path to the dereplicated output FASTA file.
//'   If \code{compress_output = TRUE} or if the file name ends with \code{.gz},
//'   the file is written in gzip-compressed format.
//'
//' @param n_threads (integer) number of threads used for processing. When
//'   OpenMP support is available, multiple threads accelerate sequence hashing
//'   and merging.
//'
//' @param compress_output (logical) whether to write the output file in gzip
//'   format regardless of the suffix of \code{output_fasta}.
//'
//' @param compress_level (integer, between 0 and 9) gzip compression level used
//'   when writing compressed output.
//'
//' @param batch_size (integer) number of records to read and process per batch.
//'   Larger values reduce I/O overhead but require more memory.
//'
//' @details
//' Each input record is parsed into:
//' \itemize{
//'   \item the sequence identifier (prefix before any metadata block),
//'   \item optional \code{"pairend"} metadata block,
//'   \item optional experiment and sample identifiers contained in
//'         \code{"stat"} or previous pipeline modules.
//' }
//' The function maintains per-sequence counters:
//' \itemize{
//'   \item \code{total_count}: total occurrences of the same sequence,
//'   \item per-PCR counts: \code{exp → sample → count}.
//' }
//' After all batches are processed, unique sequences are written to the output
//' FASTA file. Each header contains:
//' \itemize{
//'   \item \code{size=}: total count of the sequence,
//'   \item a `"stat"` JSON block summarizing per-PCR counts,
//'   \item preserved `"pairend"` metadata when present.
//' }
//'
//' @return
//' Invisibly returns \code{NULL}. The dereplicated sequences are written to
//' \code{output_fasta}. Each output record includes updated counts, preserved
//' metadata, and a single representative sequence.
//' @export
// [[Rcpp::export]]
void dereplicate(std::string input_fasta,
    std::string output_fasta,
    int n_threads = 4,
    bool compress_output=false,
    int compress_level=6,
    size_t batch_size = 100000)
{
    int max_threads = detect_system_threads();
    if (n_threads > max_threads) n_threads = max_threads;
	input_fasta = expandTilde(input_fasta);
    output_fasta = expandTilde(output_fasta);
    auto start_time = std::chrono::system_clock::now();

    show_time("Merging duplicates...");
    
    bool is_gz_in = isGzFile(input_fasta);
    FILE* fp_plain = nullptr;
    gzFile fp_gz = nullptr;
    if (is_gz_in) fp_gz = gzopen(input_fasta.c_str(), "rb");
    else fp_plain = fopen(input_fasta.c_str(), "rb");

    bool is_gz_out = false;
    {
        auto lo = output_fasta;
        transform(lo.begin(), lo.end(), lo.begin(), ::tolower);
        if ((lo.size() >= 3 && lo.substr(lo.size() - 3) == ".gz") || compress_output) {
            is_gz_out = true;
            compress_output = true;
        }
    }
    bool fastq_file = isFastqFile(input_fasta);
    // bool fastq_file = isFastqFile1(input_fasta);
    FILE* out_plain = nullptr;
    gzFile out_gz = nullptr;
    if (is_gz_out) {
        std::string mode = "wb" + std::to_string(compress_level);
        out_gz = gzopen(output_fasta.c_str(), mode.c_str());
    }
    else {
        out_plain = fopen(output_fasta.c_str(), "w");
    }
    show_time("Using " + std::to_string(n_threads) + " thread(s)");
    const size_t BATCH = batch_size;
    std::vector<FastqRecord> batch;
    batch.reserve(BATCH);

    std::unordered_map<std::string, UniqueInfo> merged;
#ifdef _OPENMP // _OPENMP
    omp_set_num_threads(n_threads);
#endif
#ifdef USE_RCPP_OUTPUT
	size_t intrrupt_count = 0;
#endif // check interrupt
    while (true) {
        batch.clear();
        FastqRecord rec;
        for (size_t i = 0; i < BATCH; ++i) {
            bool ok = readRec(is_gz_in, fp_plain, fp_gz, rec, fastq_file);
            if (!ok) break;
            batch.push_back(rec);
#ifdef USE_RCPP_OUTPUT
            intrrupt_count++;
            if (intrrupt_count % 100000 == 0) {
                Rcpp::checkUserInterrupt();
            }
#endif // check interrupt
        }
        if (batch.empty()) break;

        std::unordered_map<std::string, UniqueInfo> local_map;

#ifdef _OPENMP
#pragma omp parallel
        {
#endif
            std::unordered_map<std::string, UniqueInfo> local_private;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (size_t i = 0; i < batch.size(); ++i) {
                const auto& r = batch[i];
                std::string seq_copy = r.seq;

                std::string hdr = r.id.substr(1);
                size_t pos = hdr.find(';');
                std::string prefix = pos == std::string::npos ? hdr : hdr.substr(0, pos);
                std::string json = pos == std::string::npos ? "" : hdr.substr(pos + 1);

                // extract pairend information
                std::string pairend_json;
                size_t pe_s = json.find("\"pairend\":");
                if (pe_s != std::string::npos) {
                    size_t pe_e = json.find("}", pe_s);
                    if (pe_e != std::string::npos)
                        pairend_json = "{" + json.substr(pe_s, pe_e - pe_s + 1) + ",";
                }
                // get experiment values
                std::string exp, sample;
                size_t ex_p = json.find("\"exp\":\"");
                if (ex_p != std::string::npos) {
                    ex_p += 7;
                    size_t ex_e = json.find('"', ex_p);
                    exp = json.substr(ex_p, ex_e - ex_p);
                }
                // sample position
                size_t sm_p = json.find("\"sample\":\"");
                if (sm_p != std::string::npos) {
                    sm_p += 10;
                    size_t sm_e = json.find('"', sm_p);
                    sample = json.substr(sm_p, sm_e - sm_p);
                }

                if (exp.empty()) exp = "none";
                if (sample.empty()) sample = "no_sample";
                // create unordered_map
                auto& ui = local_private[seq_copy];
                if (ui.total_count == 0) {
                    ui.prefix = prefix;
                    ui.pairend_json = pairend_json;
                }
                ui.total_count++;
                ui.exp_samples[exp][sample]++;
            }
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                for (auto& kv : local_private) {
                    auto it = local_map.find(kv.first);
                    if (it == local_map.end()) // if seqs do not present in current container, move it in.
                        local_map.emplace(std::move(kv));
                    else { // if seqs presented in container, sum their count.
                        auto& tgt = it->second;
                        auto& src = kv.second;
                        tgt.total_count += src.total_count;
                        for (auto& e : src.exp_samples) { 
                            auto& dst_e = tgt.exp_samples[e.first]; // experiment
                            for (auto& s : e.second) // tags
                                dst_e[s.first] += s.second;
                        }
                    }
                }
            }
#ifdef _OPENMP
        }
#endif
        // merge all results
        for (auto& kv : local_map) {
            auto itG = merged.find(kv.first);
            if (itG == merged.end())
                merged.emplace(std::move(kv));
            else {
                auto& g = itG->second;
                auto& src = kv.second;
                g.total_count += src.total_count;
                for (auto& e : src.exp_samples) {
                    auto& dst_e = g.exp_samples[e.first];
                    for (auto& s : e.second)
                        dst_e[s.first] += s.second;
                }
            }
        }
    }

    // output 
    for (const auto& kv : merged) {
        const auto& seq = kv.first;
        const auto& ui = kv.second;

        std::ostringstream exp_block;
        exp_block << "{\"exp\":{";
        bool first_exp = true;
        for (const auto& e : ui.exp_samples) {
            if (!first_exp) exp_block << ",";
            first_exp = false;
            exp_block << "\"" << e.first << "\":{\"samples\":{";
            bool first_s = true;
            for (const auto& s : e.second) {
                if (!first_s) exp_block << ",";
                first_s = false;
                exp_block << "\"" << s.first << "\":" << s.second;
            }
            exp_block << "}}";
        }
        exp_block << "},\"total_count\":" << ui.total_count
            << ",\"seq_length\":" << seq.size() << "}";

        std::ostringstream full_js;
        full_js << ui.pairend_json << "\"stat\":" << exp_block.str() << "}";

        std::string header = ">" + ui.prefix +
            ";size=" + std::to_string(ui.total_count) +
            ";" + full_js.str() + "\n";
        writeBuf(header, is_gz_out, out_plain, out_gz);
        writeBuf(seq + "\n", is_gz_out, out_plain, out_gz);
    }

    if (fp_plain) fclose(fp_plain);
    if (fp_gz) gzclose(fp_gz);
    if (out_plain) fclose(out_plain);
    if (out_gz) gzclose(out_gz);

    auto end_time = std::chrono::system_clock::now();
   
    show_time("Done.");
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    show_time("", 1, elapsed_seconds.count());
}
