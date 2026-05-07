// denoise.cpp
// PCR punctual error detection using per-PCR directed labeling and cascade depth.
//
// Method summary:
// - Parse per-PCR counts from headers produced by "unique" step:
//     "stat":{"exp":{"EXP1":{"samples":{"S1":n,...}}, "EXP2":...}, ...}
//   Use composite key "EXP|SAMPLE" to avoid cross-experiment name collisions.
// - Build a global undirected neighbor list connecting sequences at exactly one edit:
//     either one substitution OR a single ±1 bp indel.
// - For each PCR independently, direct edges high_count -> low_count IF (low/high) < threshold.
//   Then assign status per PCR:
//     outdegree > 0                => "internal"
//     outdegree == 0 && indegree > 0 => "head"
//     indegree == 0 && outdegree == 0 => "singleton"
// - Compute cascade depth per PCR for internal nodes as minimal steps to a head
//   following decreasing abundance.
// - Append JSON block {"error_detect":{...}} to each FASTA header and write sequences.
//   Summary tallies across PCRs are also written.

#define _CRT_SECURE_NO_WARNINGS
#include "utils.hpp"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>
#include <regex>
#include "seq_compare.hpp"
#ifdef _OPENMP
#include <omp.h> 
#else
#endif

#ifdef USE_RCPP_OUTPUT
#include <Rcpp.h>
using namespace Rcpp;
#else
#endif
// -------------------------------------------------------------------
// Core data structures
// -------------------------------------------------------------------

struct Node {
    std::string seq;
    std::vector<uint64_t> encode;
    std::string raw_header; // full header line including '>'
    size_t total_count = 0;
    size_t seq_length;
    // Exp-Sample
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>> exp_sample_count;

    // Status
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> status; // "head" | "internal" | "singleton" | "absent"
    std::unordered_map<std::string, std::unordered_map<std::string, int>>   depth;  // min steps to nearest head; head/singleton=0; absent=-1

    // indeg and outdeg equals to 0: no neighbors, singleton
    // outdeg equals to 0 and intdeg > 0: head
    // outdeg > 0: internal
    int indeg = 0; // from others 
    int outdeg = 0; // to others

    // summary across PCRs where present
    int sum_head = 0, sum_internal = 0, sum_singleton = 0;
    // indexing
    size_t index;
};
struct Edge { size_t u, v; }; // directed one-edit neighbor pair u (child) --> v (parent)


// ------------
// -------------------------------------------------------------------
// Header parsing utilities
// -------------------------------------------------------------------

// Parse per-PCR counts from "unique" JSON header.
// Expected layout snippet:
//   ..."stat":{"exp":{"EXP1":{"samples":{"S1":n,"S2":m}}, "EXP2":{...}}, "total_count":...}...
// Fills dst["EXP|SMP"] = count.
static void parse_per_pcr_counts(const std::string& hdr,
    std::unordered_map<std::string, std::unordered_map<std::string, size_t>>& dst) {
    dst.clear();

    size_t stat = hdr.find("\"stat\"");
    if (stat == std::string::npos) return;
    size_t exp_root = hdr.find("\"exp\"", stat);
    if (exp_root == std::string::npos) return;
    
    size_t p = exp_root + 7;
    // Find exp end
    size_t depth = 1; size_t i = hdr.find('{', exp_root) + 1;
    while (i < hdr.size() && depth != 0) {
        if (hdr[i] == '{') { ++depth; }
        else if (hdr[i] == '}') { --depth; }
        
        ++i;
    }
    size_t exp_end = i;
    
    while (1) {
        // Find next "EXP_NAME"
        size_t q1 = hdr.find('\"', p);
        // q1 = hdr.find('\"', q1 );
        if (q1 == std::string::npos) break;
        size_t q2 = hdr.find('\"', q1 + 1);
       
        if (q2 == std::string::npos) break;
        std::string exp_name = hdr.substr(q1 + 1, q2 - (q1 + 1));

        // Find samples block under this EXP
        size_t samp = hdr.find("\"samples\"", q2);
        if (samp == std::string::npos) { p = q2 + 1; continue; }
        size_t lb = hdr.find('{', samp);
        if (lb == std::string::npos) break;

        // Match the { ... } of samples
        int depth = 1; size_t i = lb + 1;
        while (i < hdr.size() && depth > 0) {
            if (hdr[i] == '{') ++depth;
            else if (hdr[i] == '}') --depth;
            ++i;
        }
        if (depth != 0) break;
        size_t rb = i - 1;
        std::string block = hdr.substr(lb + 1, rb - (lb + 1));

        // Parse "Sample": number
        size_t t = 0;
        while (1) {
            size_t s1 = block.find('\"', t);
            if (s1 == std::string::npos) break;
            size_t s2 = block.find('\"', s1 + 1);
            if (s2 == std::string::npos) break;
            std::string smp = block.substr(s1 + 1, s2 - (s1 + 1));
            size_t colon = block.find(':', s2);
            if (colon == std::string::npos) { t = s2 + 1; continue; }
            size_t pos = colon + 1; while (pos < block.size() && isspace((unsigned char)block[pos])) ++pos;
            size_t val = 0; bool ok = false;
            // Safety extraction
            while (pos < block.size() && isdigit((unsigned char)block[pos])) { val = val * 10 + (block[pos] - '0'); ok = true; ++pos; }
            if (ok) dst[exp_name][smp] += val;
            t = pos;
        }
        p = hdr.find(',', i) + 1;
        if (p >= exp_end) { break; }
        std::string test = hdr.substr(p,5);
    }
}
//static const std::regex re_size(R"(size=\d+;)");
//static const std::regex re_total(R"("total_count":\d+)");
// Parse overall total_count from header if present, else sum per-PCR counts, fallback to 1.
static size_t parse_total_or_sum(const std::string& hdr,
    const std::unordered_map<std::string, std::unordered_map<std::string, size_t>>& pcr) {
    size_t p = hdr.find("\"total_count\"");
    if (p != std::string::npos) {
        p = hdr.find(':', p);
        if (p != std::string::npos) {
            ++p; while (p < hdr.size() && isspace((unsigned char)hdr[p])) ++p;
            size_t v = 0; bool ok = false;
            while (p < hdr.size() && isdigit((unsigned char)hdr[p])) { v = v * 10 + (hdr[p] - '0'); ok = true; ++p; }
            if (ok) return v;
        }
    }
    size_t sum = 0;
    for (auto& kv : pcr) {
        for (auto& skv : kv.second) {
            sum += skv.second;
        }
        
    }

    return (sum > 0) ? sum : 1u;
}

// -------------------------------------------------------------------
// One-edit neighbor tests
// -------------------------------------------------------------------

// Exactly one indel: |len diff| == 1 and contents otherwise align with one gap.
static inline bool is_indel(const std::string& a, const std::string& b, size_t max_indel)
{
    int la = (int)a.size(), lb = (int)b.size();
    size_t diff = abs(lb - la);

    // legnt diff <= max_indel
    if (diff > max_indel) return false;

    const std::string& s = (la <= lb ? a : b);  // shorter
    const std::string& t = (la <= lb ? b : a);  // longer

    size_t Ls = s.size(), Lt = t.size();

    size_t i = 0, j = 0;
    size_t used = 0;

    while (i < Ls && j < Lt) {
        if (s[i] == t[j]) {
            ++i; ++j;
        }
        else {
            ++j;
            if (++used > diff) return false;
        }
    }

    // 
    used += (Lt - j);

    return used == diff;   //
}
// -------------------------------------------------------------------
// Bit operations
// -------------------------------------------------------------------
static inline bool is_hamming_bit(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b, const size_t max_diff, bool enable_avx2){
    if (a.size() != b.size()) return false;
    size_t nword = a.size();
    return hamming_k(a.data(), b.data(), nword, max_diff, enable_avx2);
}
//--------------------------------------------------------------------
// Get count 
//--------------------------------------------------------------------
size_t get_count(size_t idx,
    const std::vector<Node>& nodes,
    const std::string& exp_key,
    const std::string& sample_key)
{

    const auto& mp = nodes[idx].exp_sample_count;

    auto itE = mp.find(exp_key);
    if (itE == mp.end()) return 0;

    auto itS = itE->second.find(sample_key);
    if (itS == itE->second.end()) return 0;

    return itS->second;
}

//--------------------------------------------------------------------
// Calculate real memory usage of nodes including dynamic allocations.
//--------------------------------------------------------------------
size_t calc_nodes_memory(const std::vector<Node>& nodes) {

    size_t total_bytes = 0;

    // vector header（24 字节左右，不重要）
    total_bytes += sizeof(nodes);

    // Node 本体占用（capacity * sizeof(Node)）
    total_bytes += nodes.capacity() * sizeof(Node);

    // 遍历每一个 Node
    for (const auto& nd : nodes) {

        // -------- raw_header --------
        total_bytes += nd.raw_header.capacity();

        // -------- seq --------
        total_bytes += nd.seq.capacity();

        // =========================================================
        // exp_sample_count: unordered_map<string, unordered_map<string, size_t>>
        // =========================================================
        total_bytes += nd.exp_sample_count.bucket_count() * sizeof(void*);  // bucket array

        for (const auto& kv : nd.exp_sample_count) {
            const std::string& exp_key = kv.first;
            const auto& inner_map = kv.second;

            // key string dynamic space
            total_bytes += exp_key.capacity();

            // inner map bucket array
            total_bytes += inner_map.bucket_count() * sizeof(void*);

            for (const auto& kv2 : inner_map) {
                total_bytes += kv2.first.capacity();  // sample name
                // kv2.second 是 size_t，不占额外动态内存
            }
        }

        // =========================================================
        // status: unordered_map<string, unordered_map<string, string>>
        // =========================================================
        total_bytes += nd.status.bucket_count() * sizeof(void*);

        for (const auto& kv : nd.status) {
            const std::string& exp_key = kv.first;
            const auto& inner_map = kv.second;

            total_bytes += exp_key.capacity();
            total_bytes += inner_map.bucket_count() * sizeof(void*);

            for (const auto& kv2 : inner_map) {
                total_bytes += kv2.first.capacity();   // sample name
                total_bytes += kv2.second.capacity();  // status string
            }
        }

        // =========================================================
        // depth: unordered_map<string, unordered_map<string, size_t>>
        // =========================================================
        total_bytes += nd.depth.bucket_count() * sizeof(void*);

        for (const auto& kv : nd.depth) {
            const std::string& exp_key = kv.first;
            const auto& inner_map = kv.second;

            total_bytes += exp_key.capacity();
            total_bytes += inner_map.bucket_count() * sizeof(void*);

            for (const auto& kv2 : inner_map) {
                total_bytes += kv2.first.capacity(); // sample name
                // kv2.second 是 size_t，占用固定，不算动态空间
            }
        }
    }

    // 打印
    double mb = static_cast<double>(total_bytes) / 1024.0 / 1024.0;
    double gb = mb / 1024.0;
   
    COUT << "Memory usage = " << total_bytes << " bytes"
        << " (" << mb << std::fixed << std::setprecision(3) <<" MB, " <<  std::fixed << std::setprecision(3) << gb << " GB)" << std::endl;

    return total_bytes;
}

// -------------------------------------------------------------------
// Build global directed neighbor list once to reuse across PCRs.
// -------------------------------------------------------------------
static std::vector<Edge> build_neighbors(const std::vector<Node>& nodes, const std::string & exp_key,const std::string &sample_key,
    std::unordered_map< std::string,std::unordered_map<std::string, std::vector<size_t>>> & nodes_index, double threshold,
    size_t max_diff, size_t indel,int n_threads, bool enable_avx2)
{
    int max_threads = detect_system_threads();
	if (n_threads > max_threads) n_threads = max_threads;

    
    const auto& temp_index = nodes_index.at(exp_key).at(sample_key);


    // Bucket by length.
    std::unordered_map<size_t, std::vector<size_t>> buckets;
    buckets.reserve(temp_index.size() * 2);
    for (size_t i : temp_index) {
        buckets[nodes[i].seq.size()].push_back(i);
    }

    std::vector<size_t> lengths;
    lengths.reserve(buckets.size());
    for (auto& kv : buckets) lengths.push_back(kv.first);
    sort(lengths.begin(), lengths.end());

    std::vector<Edge> edges;
    edges.reserve(temp_index.size() * 8);
    //
    for (size_t L : lengths) {
        const auto& current_length_index = buckets[L];// length id
        

        std::vector<size_t> current_length_indel_index;
        current_length_indel_index.reserve(1024); // 

        for (size_t diff = 1; diff <= indel; ++diff) {
            size_t target_length = L + diff;

            auto it = buckets.find(target_length);
            if (it != buckets.end()) {
                // 
                current_length_indel_index.insert(
                    current_length_indel_index.end(),
                    it->second.begin(),
                    it->second.end()
                );
            }
        }
        //if (current_length_index.size() < 2 && current_length_indel_index.size() <2) continue;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(dynamic)
#endif
        // Equal-length: one substitution
        for (size_t i = 0; i < current_length_index.size(); ++i) {
            size_t u = current_length_index[i];
            size_t u_count = get_count(u, nodes, exp_key, sample_key);
            std::vector<Edge> local;
            
                for (size_t j = i + 1; j < current_length_index.size(); ++j) {

                   
                    size_t v = current_length_index[j]; 
                    

                    size_t v_count = get_count(v, nodes, exp_key, sample_key);
                    double ratio = static_cast<double>(std::min(u_count, v_count)) / static_cast<double>(std::max(u_count, v_count));

                    if (ratio > threshold) {
                        continue;
                    }

                        uint64_t pre_u = nodes[u].encode[0];
                        uint64_t pre_v = nodes[v].encode[0];
                        uint64_t diff = pre_u ^ pre_v;
                        diff = diff | (diff >> 1) | (diff >> 2) | (diff >> 3);
                        diff = diff & 0x1111111111111111ULL;
                        if (__builtin_popcountll(diff) > max_diff) {
                            continue;
                        }
                    
                    


                    if (is_hamming_bit(nodes[u].encode, nodes[v].encode, max_diff, enable_avx2)) {
                        if (u_count < v_count) {
                            local.push_back({ nodes[u].index, nodes[v].index });
                        }
                        else {
                            local.push_back({ nodes[v].index, nodes[u].index });
                        }
                    }
                }
            
            
                // One bp difference
                for (size_t ide = 0; ide < current_length_indel_index.size(); ++ide) {
                    size_t v = current_length_indel_index[ide];

                    size_t v_count = get_count(v, nodes, exp_key, sample_key);
                    double ratio = static_cast<double>(std::min(u_count, v_count)) / static_cast<double>(std::max(u_count, v_count));
                    if (ratio > threshold) {
                        continue;
                    }

                    if (is_indel(nodes[u].seq, nodes[v].seq, indel))
                    {

                        if (u_count < v_count) {
                            local.push_back({ nodes[u].index, nodes[v].index });
                        }
                        else {
                            local.push_back({ nodes[v].index, nodes[u].index });
                        }

                    }
                }
            
            if (!local.empty()) { // 
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    edges.insert(edges.end(), local.begin(), local.end());
                }
            }
        }
    }
    return edges;
}

// -------------------------------------------------------------------
// Label for a single PCR key and compute depths.
// Direction rule: high_count -> low_count only if (low/high) < threshold.
// -------------------------------------------------------------------
static void label_one_pcr(const std::vector<Edge>& directed,
    std::vector<Node>& nodes,
    const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<size_t>>> & nodes_index,
    const std::string& exp_key,
    const std::string& sample_key,
    const double threshold)
{

	// used index
    std::vector<size_t> used_idx = nodes_index.at(exp_key).at(sample_key);

    // Reset degrees only for used nodes
    for (size_t idx : used_idx) {
        nodes[idx].indeg = 0;
        nodes[idx].outdeg = 0;
    }

    // Direct edges per PCR counts 

    for (const auto& e : directed) {
        if(e.u == e.v) continue;
		nodes[e.u].outdeg++; // 指向高丰度节点
        nodes[e.v].indeg++; //被低丰度节点指向
    }

    // Status and init depth
    for (size_t i : used_idx) {
        auto itE = nodes[i].exp_sample_count.find(exp_key);
        if (itE == nodes[i].exp_sample_count.end()) continue;

        auto itS = itE->second.find(sample_key);
       
        if (itS == itE->second.end() || itS->second == 0) {
            continue;
        }
      
        if (nodes[i].indeg == 0 && nodes[i].outdeg == 0) {
            nodes[i].status[exp_key][sample_key] = "singleton";
            nodes[i].depth[exp_key][sample_key] = 0;
        }
        else if (nodes[i].outdeg == 0 && nodes[i].indeg > 0) {
            nodes[i].status[exp_key][sample_key] = "head";
            nodes[i].depth[exp_key][sample_key] = 0;
        }
        else {
            nodes[i].status[exp_key][sample_key] = "internal";
            nodes[i].depth[exp_key][sample_key] = -1; // will be computed
        }
        
    }

    
}
// -------------------------------------------------------------------
// remove internal key function 
// -------------------------------------------------------------------
void remove_internal_key(Node& nd) {
    for (auto ep_it = nd.status.begin(); ep_it != nd.status.end(); ) {
        const std::string& ep = ep_it->first;

        for (auto sm_it = ep_it->second.begin(); sm_it != ep_it->second.end(); ) {
            if (sm_it->second == "internal") {
                const std::string& sm = sm_it->first;
				// remove count and status key
                auto ec = nd.exp_sample_count.find(ep);
                if (ec != nd.exp_sample_count.end()) ec->second.erase(sm);

                sm_it = ep_it->second.erase(sm_it);   // 
            }
            else {
                ++sm_it;
            }
        }
		// remove all empty exp key
        bool empty_status = ep_it->second.empty();
        bool empty_counts = (nd.exp_sample_count.find(ep) == nd.exp_sample_count.end()) ||
            nd.exp_sample_count[ep].empty();

        if (empty_status && empty_counts) {
            // remove useless exp_sample_count
            auto ec = nd.exp_sample_count.find(ep);
            if (ec != nd.exp_sample_count.end()) nd.exp_sample_count.erase(ec);
            ep_it = nd.status.erase(ep_it);
        }
        else {
            ++ep_it;
        }
    }
}

// ------------------------------
// Rcpp Exported Function: denoise
// ------------------------------
//' Denoise FASTA sequences using abundance and local sequence similarity
//'
//'
//' For each sequence, a `"denoise"` JSON block is appended to the header
//' summarizing per-PCR status and counts, and the overall
//' \code{size=} and \code{"total_count"} fields are updated to match the
//' retained counts. Output can be written either uncompressed or gzip-compressed.
//'
//' @param input_fasta (character) string; path to the input FASTA file. The file
//'   may be plain text or gzip-compressed (suffix \code{.gz}).
//' @param output_fasta (character )string; path to the output FASTA file. If the
//'   file name ends with \code{.gz} or \code{compress_output = TRUE}, the
//'   output is gzip-compressed.
//' @param sequence_min_L (integer); minimum sequence length to retain. Sequences
//'   shorter than this threshold are discarded before denoising.
//' @param sequence_max_L (integer); maximum sequence length to retain. Sequences
//'   longer than this threshold are discarded before denoising.
//' @param min_count (integer) minimum total read count (summed over all PCRs)
//'   required to keep a sequence. Sequences with total counts below this
//'   threshold are removed before graph construction.
//' @param threshold (double) scalar in (0,1); abundance ratio threshold used
//'   to direct edges in the neighbor graph. For a pair of neighboring
//'   sequences with counts \code{high} and \code{low}, an edge
//'   \code{high -> low} is created only if \code{low/high < threshold}. These
//'   directed edges are then used to classify \code{"head"}, \code{"internal"},
//'   and \code{"singleton"} sequences per PCR.
//' @param n_threads (integer) number of parallel threads used for neighbor
//'   search and labeling. When OpenMP is available, values greater than 1
//'   enable multi-threaded computation; otherwise this argument is ignored.
//' @param max_diff  (integer >= 0) maximum number of substitutions allowed
//'   between two sequences for them to be considered neighbors. Setting
//'   \code{max_diff = 1} restricts neighbors to Hamming distance 1 (plus any
//'   allowed indels via \code{max_indel}).
//' @param max_indel (integer >= 0) maximum absolute indel size (in base pairs)
//'   allowed between two sequences for them to be considered neighbors.
//'   \code{max_indel = 0} restricts neighbors to equal length; \code{max_indel = 1}
//'   additionally allows single-base insertions/deletions.
//' @param keep_internal (logical) if \code{FALSE}, sequences that are labeled
//'   \code{"internal"} in all PCRs are removed from the output. For sequences
//'   that are kept, any per-PCR entries with status \code{"internal"} may be
//'   dropped from the per-PCR statistics.
//' @param keep_previous_stat (logical) if \code{TRUE}, the existing JSON block
//'   in the header (typically produced by previous steps such as \code{unique})
//'   is preserved and the new `"denoise"` block is appended to it. If
//'   \code{FALSE}, the old statistics block is discarded and the header is
//'   rebuilt to contain only the sequence identifier and the new `"denoise"`
//'   block (with updated \code{size=} and \code{"total_count"} fields).
//' @param compress_output (logical) if \code{TRUE}, the output FASTA is written
//'   as gzip regardless of the file name suffix.
//' @param compress_level (integer,between 0 and 9) gzip compression level passed
//'   to \code{gzopen()} when \code{compress_output} is enabled or the output
//'   file name ends with \code{.gz}. Higher values give better compression at
//'   the cost of slower writing.
//'
//' @return
//' Invisibly returns \code{NULL}. The main result is the FASTA file written to
//' \code{output_fasta}, where each header line has been updated to include a
//' `"denoise"` JSON block with per-PCR counts and status labels, together with
//' summary tallies of \code{"head"}, \code{"internal"}, and \code{"singleton"}
//' occurrences across all PCRs.
//' @details
//' For each PCR independently, sequences are first filtered by length and total
//' abundance thresholds. A local neighbor graph is then constructed by
//' connecting sequences that differ by at most a user-defined number of
//' substitutions (\code{max_diff}) and/or short indels (\code{max_indel}). Edges
//' are directed from higher-abundance sequences to lower-abundance sequences
//' when their abundance ratio satisfies the specified \code{threshold}. This
//' directed graph is used to classify sequences into three categories per PCR:
//' \code{"head"} (no outgoing edges), \code{"internal"} (receiving edges from
//' higher-abundance neighbors), and \code{"singleton"} (no neighbors). This algorithm was described by Boyer et. al.,(2016)
//'
//' @section Citations: Boyer, F., Mercier, C., Bonin, A., Le Bras, Y., Taberlet, P. and Coissac, E. (2016), obitools: a unix-inspired software package for DNA metabarcoding. Mol Ecol Resour, 16: 176-182. \url{https://doi.org/10.1111/1755-0998.12428}
//' @export
// [[Rcpp::export]]
void denoise(std::string input_fasta,
    std::string output_fasta,
    size_t sequence_min_L = 1,
    size_t sequence_max_L = 1000,
    size_t min_count = 1,
    double threshold = 0.05,
    int n_threads = 32,
    size_t max_diff = 1,
    size_t max_indel = 0,
    bool keep_internal = true,
    bool keep_previous_stat = false,
    bool compress_output = false,
	bool enable_avx2 = true,
    bool debug = false,
    int compress_level = 6)
{
    auto start_time= std::chrono::system_clock::now();
    show_time("Denoising...");
    input_fasta = expandTilde(input_fasta);
    output_fasta = expandTilde(output_fasta);
    // Detect input type and open.
    bool is_gz_in = isGzFile(input_fasta);
    FILE* fp_plain = nullptr;
    gzFile fp_gz = nullptr;
    if (is_gz_in) fp_gz = gzopen(expandTilde(input_fasta).c_str(), "rb");
    else          fp_plain = fopen(expandTilde(input_fasta).c_str(), "rb");

    // Output type: .gz or compress_output => gzip.
    
	enable_avx2 = CPU_SUPPORTS_AVX2()  && enable_avx2;
    enable_avx2 ? show_time("AVX2 enable.") : show_time("AVX2 disable.");
    bool is_gz_out = false;
    {
       
        auto lo = output_fasta;
        std::transform(lo.begin(), lo.end(), lo.begin(), ::tolower);
        if ((lo.size() >= 3 && lo.substr(lo.size() - 3) == ".gz") || compress_output) {
            is_gz_out = true;

        }
    }
    FILE* out_plain = nullptr;
    gzFile out_gz = nullptr;
    if (is_gz_out) {
        std::string mode = "wb" + std::to_string(compress_level);
        out_gz = gzopen(output_fasta.c_str(), mode.c_str());
    }
    else {
        out_plain = fopen(output_fasta.c_str(), "w");
    }


    // Read all records into nodes.
    std::vector<Node> nodes;
    nodes.reserve(1 << 16);
    std::unordered_set<std::string> exp_name_set;
    std::unordered_set <std::string> sample_set;
    FastqRecord rec;
    show_time("Loading sequences...");
    while (true) {
        if (!readRec(is_gz_in, fp_plain, fp_gz, rec)) break;
        Node nd;
        nd.raw_header = rec.id;
        nd.seq = rec.seq;
		nd.seq_length = nd.seq.size();
        parse_per_pcr_counts(nd.raw_header, nd.exp_sample_count);
        for (auto& kv : nd.exp_sample_count) {
            exp_name_set.insert(kv.first);
             for (auto& skv : kv.second) {
                 sample_set.insert(skv.first);
             }
        }
        nd.total_count = parse_total_or_sum(nd.raw_header, nd.exp_sample_count);
        // Filter out not target length sequences
        if (sequence_min_L > nd.seq_length || sequence_max_L < nd.seq_length) continue;
        // Filter out low count sequences
        if (nd.total_count < min_count) continue;
        //
		nd.encode = encode4bit(nd.seq);
        nodes.push_back(std::move(nd));
    }
    if (fp_plain) fclose(fp_plain);
    if (fp_gz)    gzclose(fp_gz);


    std::vector<std::string> exp_name(exp_name_set.begin(), exp_name_set.end());

	std::vector<std::string> sample(sample_set.begin(), sample_set.end());
    calc_nodes_memory(nodes);
    // Stable output order
    show_time("Sorting sequences...");
    sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        if (a.total_count != b.total_count) return a.total_count > b.total_count;
        return a.seq.size() < b.seq.size();
        });

    // nodes_index：exp  sample vector<node_index>
    std::unordered_map<
        std::string,
        std::unordered_map<std::string, std::vector<size_t>>
    > nodes_index;

    for (size_t N = 0; N < nodes.size(); ++N) {
        nodes[N].index = N; // indexing
        const auto& m = nodes[N].exp_sample_count;

        for (const auto& kv_exp : m) {
            const std::string& exp = kv_exp.first;   // experiment 
            const auto& sample_map = kv_exp.second;  // sample

            for (const auto& kv_smp : sample_map) {
                const std::string& smp = kv_smp.first;
                size_t count = kv_smp.second;

                if (count > 0) {
                    nodes_index[exp][smp].push_back(N);
                }
            }
        }
    }
   
   
	// int count = 0;
    show_time("Labeling sequences...");
	// progress tracking
    const size_t total = exp_name.size() * sample.size();
    size_t progress = 0;
    std::string pecbar = std::string(50, '=');
    size_t bar_current_pos = 1;
	double percent = 0.0;
 
    std::ostringstream oss;
    //
#ifdef USE_RCPP_OUTPUT
    size_t intrrupt_count = 0;
#endif // check interrupt
    for(size_t e =0; e < exp_name.size();e++){
     
#ifdef USE_RCPP_OUTPUT
        intrrupt_count++;
        if (intrrupt_count % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
#endif // check interrupt
		std::string ep = exp_name[e];
       
        auto it_exp = nodes_index.find(ep);
        if (it_exp == nodes_index.end()) continue;
        for(size_t s =0; s < sample.size(); s ++){
            oss.str("");   // 
            oss.clear();
            std::vector<Edge> directed;
			std::string sm = sample[s];
            //
            

            auto it_smp = it_exp->second.find(sm);
            if (it_smp == it_exp->second.end()) continue;
            progress++;
            
            //
            percent = (double)progress / double(total);
           
                bar_current_pos = int(percent * 100) / 2;
				pecbar.std::string::replace(0,bar_current_pos, std::string(bar_current_pos, '-'));
                pecbar[bar_current_pos] = '>';
                
               
            oss << "[" + pecbar + "] " << std::fixed << std::setprecision(1) << percent * 100 << "% Labeling experienment [" << ep <<"] sample " << "[" << sm <<"]     ";
            //std::this_thread::sleep_for(std::chrono::milliseconds(500));
            show_time(oss.str(), 2);
            
          
            directed = build_neighbors(nodes, ep, sm, nodes_index, threshold,max_diff, max_indel, n_threads, enable_avx2);
		
            label_one_pcr(directed, nodes, nodes_index , ep, sm, threshold);

            
            
        }

    }
    oss.str("");   // 
    oss.clear();
    oss << "[" + std::string(50, '-') + "] " << std::fixed << std::setprecision(1) << 100.0 << "% Done.                                          ";
    show_time(oss.str(), 2);
    // Summaries freq
	COUT << std::endl;
    show_time("Summaries...");
    for (auto& nd : nodes) {
        int h = 0, in = 0, s = 0;
        for (auto& kv : nd.status) {
            for (auto& skv : kv.second) {
                if (skv.second == "head")      ++h;
                else if (skv.second == "internal")  ++in;
                else if (skv.second == "singleton") ++s;
            }
        }
        nd.sum_head = h; nd.sum_internal = in; nd.sum_singleton = s;
    }


    show_time("Writting sequences...");
    // Write with appended JSON block.
    for (auto& nd : nodes) {
		// Filter out internal-only sequences if requested
        if (!keep_internal) {
            if(nd.sum_head == 0 && nd.sum_singleton == 0) {
                continue; // skip internal-only sequences
			}
            remove_internal_key(nd);
        }
        // Build JSON
        std::ostringstream js;
        js << "\"denoise\":{"
            << "\"threshold\":" << std::fixed << std::setprecision(3) << threshold << ","
            << "\"total_count\":" << nd.total_count << ","
            << "\"summary\":{\"head\":" << nd.sum_head
            << ",\"internal\":" << nd.sum_internal
            << ",\"singleton\":" << nd.sum_singleton << "},"
            << "\"by_sample\":{";
        bool first_exp = true;
        size_t new_size = 0;
        for (const auto &ep : exp_name) {
			auto exp_status = nd.exp_sample_count.find(ep);
			if  (exp_status == nd.exp_sample_count.end()) continue;
            // Decide filter if label_only is off.
            //
            if (!first_exp) {
                js << ",";
            }
            js << "\"" << ep << "\":{";
            bool first_sample = true;
            first_exp = false;
            for (const auto& sm : sample) {
				auto sample_status = exp_status->second.find(sm);
                if (sample_status == exp_status->second.end() || sample_status->second == 0)  continue; 
                new_size += sample_status->second;
                if (!first_sample) {
					js << ",";
                }
                first_sample = false;
                js << "\"" << sm << "\":{" << "\"count\":" << sample_status->second << ",\"status\":\"" << nd.status.at(ep).at(sm)  << "\"}"; //<< "\",\"depth\":" << nd.depth.at(ep).at(sm)
            }
            js << "}";
        }
        js << "}}}";
		if (new_size < 10) continue; // skip low count seqs
        // Header + sequence
        std::string hdr = nd.raw_header;
        // Update count if remove internal
        // 1. size=xxx;

        // 
// -------- Replace size=xxx; -----------
        {
            size_t p = hdr.find("size=");
            if (p != std::string::npos) {
                size_t q = hdr.find(";", p);
                if (q != std::string::npos) {
                    hdr.replace(p, q - p + 1, "size=" + std::to_string(new_size) + ";");
                }
            }
		
            p = hdr.find("\"total_count\":");
            if (p != std::string::npos) {
                size_t q = hdr.find(",", p);
                if (q != std::string::npos) {
                    hdr.replace(p, q - p + 1, "\"total_count\":" + std::to_string(new_size) + ",");
                }
            }
        }
        //------------------------------
        
        if (!hdr.empty() && hdr.back() == '\n') hdr.pop_back();
        if (keep_previous_stat) {
            hdr.erase(hdr.rfind('}'));
            writeBuf(hdr + "," + js.str() + "\n", is_gz_out, out_plain, out_gz);
            writeBuf(nd.seq + "\n", is_gz_out, out_plain, out_gz);
        }
        else {
            hdr.erase(hdr.find_first_of("{"));
            writeBuf(hdr + "{" + js.str() + "\n", is_gz_out, out_plain, out_gz);
            writeBuf(nd.seq + "\n", is_gz_out, out_plain, out_gz);
        }

    }

    if (out_plain) fclose(out_plain);
    if (out_gz)    gzclose(out_gz);

    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    show_time("Done.");
    show_time("", 1, elapsed.count());
}
