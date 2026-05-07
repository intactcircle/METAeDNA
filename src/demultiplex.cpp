
#define _CRT_SECURE_NO_WARNINGS
// [[Rcpp::plugins(cpp11)]]
//#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdio>         // For FILE*, fopen, fgets, fread, fclose, fwrite
#include <cstdlib>        // For getenv, exit
#include <cstring>        // For strlen, memset, memcpy
#include <algorithm>      // For std::min, std::reverse, std::transform, std::max
#include <zlib.h>         // For gzFile, gzopen, gzgets, gzwrite, gzclose
#include <cmath>          // For pow
#include <array>          // For std::array
#include <chrono>
#include <ctime>
#include <iomanip>  // for put_time
#include <bitset>
#include <set>
#include <thread>
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
//


// ---------------------------------------
// FUNCTION: read primer & tag
// ---------------------------------------
struct SampleInfo {
    std::vector<std::string> exp, sample, ftags,rtags, forward_primer, reverse_primer;
     std::set<std::pair<std::string, std::string >> unique_primer_pair;
     std::vector<size_t> fwd_pos, rev_pos;
     std::vector<bool> try_reverse;
};
struct hitInfo {
    std::string fwd, rev, ftag, rtag, exp, sample;
    int f_mismatch, r_mismatch;
};
static void readPrimerFile(
    const std::string& filename,
    SampleInfo &dem_infos
) {
    COUT << "[Supported COMMA and TAB delimiter.]" << std::endl;
    std::ifstream file(filename);
    if (!file.is_open()) {
        CERR << "Could not open: " << filename << std::endl;
        return;
    }
    show_time("Table use for demultiplexing was opened...");
    
    std::string line;
    char delimiter = '\t';  // tab delim was default

    // reading header
    if (!std::getline(file, line)) {
        CERR << "[ERROR] Table format was illegal!" << std::endl;
        return;
    }

    // auto detect delim
    if (line.find(',') != std::string::npos) {
        delimiter = ',';
        COUT << "Table delimiter [COMMA]." << std::endl;
    }else {
        COUT << "Table delimiter [TAB]." << std::endl;
    }

    
    // using istringstream to parse header
    std::istringstream headerStream(line);
    std::string header;
    std::vector<std::string> columnNames;
    while (std::getline(headerStream, header, delimiter)) {
        // remove useless characters
        header.erase(0, header.find_first_not_of(" \t\r\n\"#"));
        header.erase(header.find_last_not_of(" \t\r\n\"") + 1);
        columnNames.push_back(header);
    }

    // build index map, seperate column names to index its own
    std::unordered_map<std::string, size_t> colIndex;
    for (size_t i = 0; i < columnNames.size(); ++i) {
        colIndex[columnNames[i]] = i;
    }

    // print debug info
    COUT << "Table colnames: ";
    for (const auto& col : columnNames) {
        COUT << "[" << col << "] ";
    }
    std::vector < std::pair < std::string, std::string >> primer_pair;
    COUT << std::endl;
    // read data line
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;  // skip empty
        }

        std::istringstream lineStream(line);
        std::string cell;
        std::vector<std::string> row;

        while (std::getline(lineStream, cell, delimiter)) {
            cell.erase(0, cell.find_first_not_of(" \t\r\n\""));
            cell.erase(cell.find_last_not_of(" \t\r\n\"") + 1);
            row.push_back(cell);
        }


		// check index and push back
        if (!row.empty()) {
            if (colIndex.count("exp") && colIndex["exp"] < row.size()) {
                dem_infos.exp.push_back(row[colIndex["exp"]]);
            }
            if (colIndex.count("sample") && colIndex["sample"] < row.size()) {
                dem_infos.sample.push_back(row[colIndex["sample"]]);
            }
            if (colIndex.count("tags") && colIndex["tags"] < row.size()) {
                std::string temp_tag = row[colIndex["tags"]];
                for (auto& x : temp_tag) {
                    x = static_cast<char>(std::toupper(static_cast<unsigned char>(x)));
                }
                size_t colon_pos = temp_tag.find(':');
                if (colon_pos != std::string::npos&&colon_pos > 0 && colon_pos + 1 < temp_tag.size()) {
                    dem_infos.ftags.push_back(temp_tag.substr(0, colon_pos));
                    dem_infos.rtags.push_back(temp_tag.substr(colon_pos + 1));
                }
                else {
                    dem_infos.ftags.push_back(temp_tag);
                    dem_infos.rtags.push_back(temp_tag);
        
                }

            }
            if (colIndex.count("forward_primer") && colIndex["forward_primer"] < row.size()) {
                auto& primer = row[colIndex["forward_primer"]];
                std::transform(primer.begin(), primer.end(), primer.begin(),
                    [](unsigned char c) { return std::toupper(c); });
                dem_infos.forward_primer.push_back(row[colIndex["forward_primer"]]);
                
            }
            if (colIndex.count("reverse_primer") && colIndex["reverse_primer"] < row.size()) {
                auto& primer = row[colIndex["reverse_primer"]];
                std::transform(primer.begin(), primer.end(), primer.begin(),
                    [](unsigned char c) { return std::toupper(c); });
                dem_infos.reverse_primer.push_back(row[colIndex["reverse_primer"]]);
               
            }
            primer_pair.push_back(std::make_pair(row[colIndex["forward_primer"]], row[colIndex["reverse_primer"]]));
        }
        //find unique primer pairs
        std::set<std::pair<std::string, std::string>> unique_pairs(
            primer_pair.begin(), primer_pair.end());
        dem_infos.unique_primer_pair = unique_pairs;


        
    }
    COUT << "Table contains " + std::to_string(dem_infos.sample.size()) + " entries." << std::endl;
    file.close();
}

// ---------------------------------------
// FUNCTION: find primer and tag
// ---------------------------------------
std::vector<size_t> find_base_exact(const std::string& pattern, const std::string& text)
{
    size_t psl = pattern.size();
    size_t tsl = text.size();
    if (tsl < psl)
    {
        COUT << "Pattern is longer than text, please use a shorter pattern." << std::endl;
        return {};
    }
    if (psl > 64)
    {
        COUT << "Pattern is too long, please use a shorter pattern." << std::endl;
        return {};
    }
    std::vector<size_t> positions;
    // cout << pattern.size() << endl;
    long long unsigned bitmask[256];
    std::fill(bitmask, bitmask + 256, ~0ULL);
    for (size_t i = 0; i < psl; i++)
    {
        bitmask[(unsigned char)pattern[i]] &= ~(1ULL << i);
        // cout << bitset<64>(bitmask[(unsigned char)pattern[i]]) << endl;
        if (pattern[i] != 'A' && pattern[i] != 'C' && pattern[i] != 'G' && pattern[i] != 'T')

        {
            switch (pattern[i])
            {
            case 'W': // W = A or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'K': // K = G or T
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);

                break;
            case 'M': // M = A or C
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);

                break;
            case 'B': // B = C or G or T
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'D': // D = A or G or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'H': // H = A or C or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'V': // V = A or C or G
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            case 'N': // N = A or C or G or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'R': // R = A or G
                bitmask['A'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            case 'Y': // Y = C or T
                bitmask['C'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'S': // S = C or G
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            default:
                break;
            }
        }
    }

    long long unsigned D = ~0ULL;
    long long unsigned matchbit = (1ULL << (psl - 1));
    // cout << "test: "<< bitset<64>(bitmask['D']) << endl;
    //cout << "initial state: " << bitset<64>(D) << endl;
    for (size_t s = 0; s < tsl; s++)
    {
        //cout << "state: " << bitset<64>(D) << endl;
        //cout << "CHAR: " << text[s] << " bitmask: " << bitset<64>(bitmask[(unsigned char)text[s]]) << endl;
        D = (D << 1) | bitmask[(unsigned char)text[s]];
        if ((D & matchbit) == 0)
        {
            // cout << "found at " << s - pattern.size() + 1 << endl;
            positions.push_back(s - psl + 1);
        }
    }
    return positions;
}


std::vector<size_t> find_base_mismatch(const std::string& pattern, const std::string& text, const int& max_mismatch, int& mis_num)   // 
{
    size_t psl = pattern.size();
    size_t tsl = text.size();
    if (tsl < psl)
    {
        COUT << "Pattern is longer than text." << std::endl;
        return {};
    }
    if (psl >= 64)
    {
        COUT << "Pattern must shorter than 64 bp." << std::endl;
        return {};
    }

    // cout << pattern.size() << endl;
    long long unsigned bitmask[256];
    std::fill(bitmask, bitmask + 256, ~0ULL);
    for (size_t i = 0; i < psl; i++)
    {
        bitmask[(unsigned char)pattern[i]] &= ~(1ULL << i);
        if (pattern[i] != 'A' && pattern[i] != 'C' && pattern[i] != 'G' && pattern[i] != 'T')

        {
            switch (pattern[i])
            {
            case 'W': // W = A or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'K': // K = G or T
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);

                break;
            case 'M': // M = A or C
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);

                break;
            case 'B': // B = C or G or T
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'D': // D = A or G or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'H': // H = A or C or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'V': // V = A or C or G
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            case 'N': // N = A or C or G or T
                bitmask['A'] &= ~(1ULL << i);
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'R': // R = A or G
                bitmask['A'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            case 'Y': // Y = C or T
                bitmask['C'] &= ~(1ULL << i);
                bitmask['T'] &= ~(1ULL << i);
                break;
            case 'S': // S = C or G
                bitmask['C'] &= ~(1ULL << i);
                bitmask['G'] &= ~(1ULL << i);
                break;
            default:
                break;
            }
        }
    }
    //////////////

    // Bitap 状态
    std::vector<uint64_t> R(max_mismatch + 1, ~1ULL);
    std::vector<size_t> positions;
    int min_used = max_mismatch + 1;

    // 扫描 text
    for (size_t s = 0; s < tsl; ++s) {
        uint64_t prev = R[0];
        // d = 0
        R[0] = (R[0] | bitmask[(unsigned char)text[s]]) << 1;
        // d > 0
        for (int d = 1; d <= max_mismatch; ++d) {
            uint64_t tmp = R[d];
            R[d] = ((R[d] | bitmask[(unsigned char)text[s]]) & prev) << 1;
            prev = tmp;
        }
        // 检测匹配
        if ((R[max_mismatch] & (1ULL << psl)) == 0) {
            positions.push_back(s - psl + 1);
            // 找最小的 d
            for (int d = 0; d <= max_mismatch; ++d) {
                if ((R[d] & (1ULL << psl)) == 0) {
                    min_used = std::min(min_used, d);
                    break;
                }
            }
        }
    }

    // set -1
    mis_num = positions.empty() ? -1 : min_used;
    return positions;
}

//---------------------------
//---------------------------
std::string extract_tagged_seq(
    FastqRecord& rec,
    const SampleInfo& dem_infos,
    const bool& with_tag,
    int& primer_mismatch,
    hitInfo& hit_info
)
{
    bool try_reverse;
    hit_info = {};
    std::string temp_result;
    std::string temp_match_rev_primer, temp_match_fwd_primer, temp_match_ftag, temp_match_rtag;
    // matching primer 
    struct matchedInfo {
        std::vector <std::string> forward, reverse;
        std::vector <size_t> fwd_pos, rev_pos;
        std::vector <bool> try_reverse;
    };
    // Store match data
    matchedInfo minfo;

    for (const auto& p : dem_infos.unique_primer_pair) {
        std::string temp_rev_primer;
        try_reverse = false;
        std::vector<size_t> fwd_pos, rev_pos;
        std::string tmp_raw_fwd = p.first;
        std::string tmp_raw_rev = p.second;
        std::string rc;

        // 1. raw direction
        fwd_pos = find_base_mismatch(tmp_raw_fwd, rec.seq, primer_mismatch, hit_info.f_mismatch);
        fast_reverse_complement(tmp_raw_rev, temp_rev_primer);
        rev_pos = find_base_mismatch(temp_rev_primer, rec.seq, primer_mismatch, hit_info.r_mismatch);

        // 2. try direction
        if (fwd_pos.empty() || rev_pos.empty()) {
            try_reverse = true;
            fast_reverse_complement(rec.seq, rc);
            fwd_pos = find_base_mismatch(tmp_raw_fwd, rc, primer_mismatch, hit_info.f_mismatch);
            fast_reverse_complement(tmp_raw_rev, temp_rev_primer);
            rev_pos = find_base_mismatch(temp_rev_primer, rc, primer_mismatch, hit_info.r_mismatch);

        }
        // keep it, very useful
        if (fwd_pos.empty() || rev_pos.empty()) {}
        else {
            minfo.forward.push_back(tmp_raw_fwd);
            minfo.reverse.push_back(tmp_raw_rev);
            minfo.try_reverse.push_back(try_reverse);
            // Do not consider nesting

            minfo.fwd_pos.push_back(fwd_pos[0]);
            minfo.rev_pos.push_back(rev_pos.back());

        }
    }
    if (minfo.reverse.empty()) { return ""; }
    temp_result.clear();
    SampleInfo filtered_dem_infos;
    // filter
    for (size_t ind = 0; ind < dem_infos.forward_primer.size(); ind++) {
        for (size_t m = 0; m < minfo.forward.size(); m++) {
            if ((dem_infos.forward_primer[ind] == minfo.forward[m]) && (dem_infos.reverse_primer[ind] == minfo.reverse[m])) {
                filtered_dem_infos.reverse_primer.push_back(dem_infos.reverse_primer[ind]);
                filtered_dem_infos.forward_primer.push_back(dem_infos.forward_primer[ind]);
                filtered_dem_infos.ftags.push_back(dem_infos.ftags[ind]);
                filtered_dem_infos.rtags.push_back(dem_infos.rtags[ind]);
                filtered_dem_infos.sample.push_back(dem_infos.sample[ind]);
                filtered_dem_infos.exp.push_back(dem_infos.exp[ind]);
                //
                filtered_dem_infos.fwd_pos.push_back(minfo.fwd_pos[m]);
                filtered_dem_infos.rev_pos.push_back(minfo.rev_pos[m]);
                filtered_dem_infos.try_reverse.push_back(minfo.try_reverse[m]);
            }
        }
    }
    //
    for (size_t ind = 0; ind < filtered_dem_infos.forward_primer.size(); ind++) {
        std::string temp_seq_data;
        if (filtered_dem_infos.try_reverse[ind]) {
            fast_reverse_complement(rec.seq, temp_seq_data);
        }
        else
        {
            temp_seq_data = rec.seq;
        }

        // 3. 匹配tag逻辑
        std::string temp_rev_revtag;
        size_t f_primer_size, r_primer_size;
        f_primer_size = filtered_dem_infos.forward_primer[ind].size();
        r_primer_size = filtered_dem_infos.reverse_primer[ind].size();
        //
        size_t fwd_tag_len = 0, rev_tag_len = 0;
        std::string tmp_raw_f_tag, tmp_raw_r_tag;

        // 4. 判断序列区间
        size_t f1 = filtered_dem_infos.fwd_pos[ind];
        size_t r1 = filtered_dem_infos.rev_pos[ind];



        std::string fwd_subseq, rev_subseq;
        if (f1 >= r1) {
			continue;
            /*
            size_t temp;
            temp = filtered_dem_infos.fwd_pos[ind];
            filtered_dem_infos.fwd_pos[ind] = filtered_dem_infos.rev_pos[ind];
            filtered_dem_infos.rev_pos[ind] = temp;
            size_t temp_f;
            temp_f = f1;
            f1 = r1;
            r1 = temp_f;
            temp_f = f_primer_size;
            f_primer_size = r_primer_size;
            r_primer_size = temp_f;
            */
        }

        //cout << rec.seq << endl;
        // 5. 不含tag的情况
        if (!with_tag) {

            // raw
            size_t start = f1 + f_primer_size;
            size_t len = r1 - f1 - f_primer_size;
            if (start + len > temp_seq_data.size()) continue;
            temp_result.assign(temp_seq_data.begin() + start, temp_seq_data.begin() + start + len);


            if (temp_result.size() != 0) {

                //hit_info.ftag = tmp_raw_f_tag;
                //hit_info.rtag = tmp_raw_r_tag;
                hit_info.fwd = filtered_dem_infos.forward_primer[ind];
                hit_info.rev = filtered_dem_infos.reverse_primer[ind];
                hit_info.sample = filtered_dem_infos.sample[ind];
                hit_info.exp = filtered_dem_infos.exp[ind];
                return temp_result;

            }

        }
        else {
            //
            fwd_tag_len = filtered_dem_infos.ftags[ind].size();
            rev_tag_len = filtered_dem_infos.rtags[ind].size();
            //
            fwd_subseq = temp_seq_data.substr(((int)f1 - (int)fwd_tag_len) >= 0 ? ((int)f1 - (int)fwd_tag_len) : 0, fwd_tag_len);
            rev_subseq = temp_seq_data.substr((r1 + r_primer_size), (r1 + r_primer_size + rev_tag_len) >= temp_seq_data.size() ? std::string::npos : rev_tag_len);

            //
            tmp_raw_f_tag = filtered_dem_infos.ftags[ind];
            tmp_raw_r_tag = filtered_dem_infos.rtags[ind];
            // 6. if left seq length shorter than tag length, continue
            if (rev_subseq.size() < rev_tag_len || fwd_subseq.size() < fwd_tag_len) continue;

            // 7.match tag
            std::vector<size_t> fwd_tag_pos, rev_tag_pos;


            fwd_tag_pos = find_base_exact(tmp_raw_f_tag, fwd_subseq);
            fast_reverse_complement(tmp_raw_r_tag, temp_rev_revtag);
            rev_tag_pos = find_base_exact(temp_rev_revtag, rev_subseq);
            // 8. 判断tag匹配边界
            if (!fwd_tag_pos.empty() && !rev_tag_pos.empty()) {
                if (fwd_tag_pos[0] == 0 && rev_tag_pos[0] == 0) {
                    size_t start = f1 + f_primer_size;
                    size_t len;
                    if (r1 >= (f1 + f_primer_size)) {
                        len = r1 - f1 - f_primer_size;
                    }
                    else { continue; }
                    if (start + len > temp_seq_data.size()) continue;
                    temp_result.assign(temp_seq_data.begin() + start, temp_seq_data.begin() + start + len);

                }

            }
            else {
                continue;
            }

            if (temp_result.size() != 0) {
                hit_info.ftag = tmp_raw_f_tag;
                hit_info.rtag = tmp_raw_r_tag;
                hit_info.fwd = filtered_dem_infos.forward_primer[ind];
                hit_info.rev = filtered_dem_infos.reverse_primer[ind];
                hit_info.sample = filtered_dem_infos.sample[ind];
                hit_info.exp = filtered_dem_infos.exp[ind];
                return temp_result;
            }
        }


    }
    return "";

}

// ------------------------------
// Rcpp Exported Function: demultiplex
// ------------------------------
//' Demultiplex FASTQ reads using primer and tag information
//'
//' This function assigns each FASTQ read to a sample based on matching forward
//' and reverse primers, as well as optional tag sequences. For each read, the
//' identified experiment, sample, tag, and primer mismatches are recorded in a
//' `"demultiplex"` JSON block appended to the FASTA header.
//'
//' @param fastq_file (character) path to the input FASTQ file. The file may be
//'   plain text or gzip-compressed (suffix \code{.gz}).
//'
//' @param primer_table (character) path to a CSV/TSV file containing primer and
//'   tag definitions. This table must include experiment names, sample names,
//'   forward and reverse primers, and optional tag sequences.
//'
//' @param identified_output (character) path to the output file that receives
//'   reads successfully assigned to a sample. If \code{compress_output = TRUE}
//'   or the file name ends with \code{.gz}, the file is gzip-compressed.
//'
//' @param unidentified_output (character) path to the output file for reads
//'   that cannot be matched to any sample based on the specified primers and
//'   tags. Compression behavior is identical to \code{identified_output}.
//'
//' @param n_threads (integer) number of parallel threads used for matching
//'   primers and tags. When OpenMP support is available, multiple threads
//'   accelerate demultiplexing.
//'
//' @param with_tag (logical) if \code{TRUE}, tag sequences must be present
//'   and matched in addition to primers. If \code{FALSE}, only primer matching
//'   is performed.
//'
//' @param compress_output (logical) if \code{TRUE}, write all output files in
//'   gzip format regardless of output file name suffix.
//'
//' @param compress_level (integer, between 0 and 9) gzip compression level
//'   applied when writing compressed output.
//'
//' @param primer_mismatch (integer) maximum number of mismatches allowed when
//'   matching primers. Tag matching does not allow mismatches.
//'
//' @return
//' Invisibly returns \code{NULL}. Identified reads are written to
//' \code{identified_output}, whereas unmatched reads are written to
//' \code{unidentified_output}. Each identified read receives a `"demultiplex"`
//' JSON block in the header containing experiment, sample, tag, primer matches,
//' and mismatch counts.
//' @export
// [[Rcpp::export]]
void demultiplex( std::string fastq_file,  std::string primer_table,
     std::string identified_output,  std::string unidentified_output,
    int n_threads = 4, bool with_tag = true,
    bool compress_output = false, int compress_level = 6,
    int primer_mismatch = 3)
{
    int max_threads = detect_system_threads();
    if (n_threads > max_threads) n_threads = max_threads;
    size_t batch_size = 1000000;
    size_t total_count = 0;
    size_t identified_count = 0;
    size_t unidentified_count = 0;
    auto start_time = std::chrono::system_clock::now();
    auto now = std::chrono::system_clock::now();
    fastq_file =  expandTilde(fastq_file);
    primer_table = expandTilde(primer_table);
    identified_output = expandTilde(identified_output);
    unidentified_output = expandTilde(unidentified_output);

    // ------------------------------
    // 1. Read demutiplex table
    // ------------------------------
    SampleInfo infos;
    readPrimerFile(primer_table, infos);

    // ------------------------------
    // 2. detect file type
    // ------------------------------
    bool is_gz_in = isGzFile(fastq_file);
    FILE* fp_plain = nullptr;
    gzFile fp_gz = nullptr;
    if (is_gz_in) {

        fp_gz = gzopen(fastq_file.c_str(), "rb");
    }
    else {
        fp_plain = fopen(fastq_file.c_str(), "rb");
    }
    /*
    
    auto readRec = [&](FastqRecord& r)->bool {
        static const unsigned BUF_SZ = 512 * 1024;
        static thread_local char* buffer = new char[BUF_SZ];
        return is_gz_in
            ? readNextRecordGz(fp_gz, r, buffer, BUF_SZ)
            : readNextRecordPlain(fp_plain, r, buffer, BUF_SZ);
        };
     */
    // ------------------------------
    // 3. output file type
    // ------------------------------
    bool is_gz_out = false;
    {
        auto lo = identified_output;
        std::transform(lo.begin(), lo.end(), lo.begin(), ::tolower);
        if ((lo.size() >= 3 && lo.substr(lo.size() - 3) == ".gz")||compress_output) {
            is_gz_out = true;
            compress_output = true;
        }
    }
    FILE* out_plain_id = nullptr;
    FILE* out_plain_un = nullptr;
    gzFile out_gz_id = nullptr;
    gzFile out_gz_un = nullptr;
    if (is_gz_out) {
        std::string mode = "wb" + std::to_string(compress_level);
        out_gz_id = gzopen(identified_output.c_str(), mode.c_str());
        out_gz_un = gzopen(unidentified_output.c_str(), mode.c_str());
    }
    else {
        out_plain_id = fopen(identified_output.c_str(), "wb");
        out_plain_un = fopen(unidentified_output.c_str(), "wb");
    }
    /*
    auto writeBuf = [&](const std::string& buf, bool identified) {
        if (is_gz_out) {
            gzFile out = identified ? out_gz_id : out_gz_un;
            gzwrite(out, buf.data(), static_cast<unsigned int>(buf.size()));
        }
        else {
            FILE* out = identified ? out_plain_id : out_plain_un;
            fwrite(buf.data(), 1, buf.size(), out);
        }
        };
    */
    // ------------------------------
    // 4. read in batch and then demultiplex
    // ------------------------------
    
    const size_t BATCH = batch_size;
    std::vector<FastqRecord> batch;
    batch.reserve(BATCH);
    show_time("Demultiplexing...");
    show_time("Using " + std::to_string(n_threads) + " thread(s)");
#ifdef USE_RCPP_OUTPUT
    size_t intrrupt_count = 0;
#endif // check interrupt
    while (true) {
        batch.clear();
        FastqRecord rec;
        for (size_t i = 0; i < BATCH; ++i) {
            if (!readRec(is_gz_in, fp_plain, fp_gz,rec,true)) break;
           
            batch.push_back(rec);
        }
        if (batch.empty()) break;
        total_count += batch.size();


        now = std::chrono::system_clock::now();
        // std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        std::ostringstream temp_string;
        temp_string << "Processed:" << total_count << " records(demultiplexing)...";
        show_time(temp_string.str(), 2);
        // Setting openmp

        std::vector<std::vector<std::string>> thrId(n_threads), thrUn(n_threads);
#ifdef USE_RCPP_OUTPUT
        intrrupt_count++;
        if (intrrupt_count % 100000 == 0) {
            Rcpp::checkUserInterrupt();
        }
#endif // check interrupt
#ifdef _OPENMP
#pragma omp parallel num_threads(n_threads)
        {
#endif
         
            std::string local;
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int i = 0; i < (int)batch.size(); ++i) {
                auto& r = batch[i];
                bool matched = true;
                hitInfo hit_info;
               
                std::string sub_seq = extract_tagged_seq(r, infos, with_tag, 
                    primer_mismatch, hit_info);

                local.clear();
                local.reserve(sub_seq.size() + 256);
                    if (sub_seq.empty()) matched = false;
                    if (matched) {
                       

                        local = ">" + r.id.substr(1, r.id.size() - 1)
                            + "{\"demultiplex\":{\"sample\":\"" + hit_info.sample
                            + "\",\"exp\":\"" + hit_info.exp
                            + "\",\"tag\":\"" + (with_tag ? hit_info.ftag + "+" + hit_info.rtag : "none")
                            + "\",\"fwd\":\"" + hit_info.fwd
                            + "\",\"rev\":\"" + hit_info.rev
                            + "\",\"fwd_mismatch\":\"" + std::to_string(hit_info.f_mismatch)
                            + "\",\"rev_mismatch\":\"" + std::to_string(hit_info.r_mismatch)
                            + "\"}\n"
                            + sub_seq + "\n";
                    }
                    else { 
                       
                        local = ">" + r.id.substr(1, r.id.size() - 1)
                            +  "\n"
                            + r.seq + "\n";
                    }

                    // 
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
                    {
                        if (matched) { ++identified_count; }
                        else { ++unidentified_count; }
                        //writeBuf(local, matched);
                        writeBuf(local, is_gz_out,
                            matched ? out_plain_id : out_plain_un,
                            matched ? out_gz_id : out_gz_un);
                    }
                       
                    
                }
#ifdef _OPENMP
    }
#endif
        }

    
    
    // ------------------------------
    // 5. 清理
    // ------------------------------
    if (fp_plain)     fclose(fp_plain);
    if (fp_gz)        gzclose(fp_gz);
    if (out_plain_id) fclose(out_plain_id);
    if (out_plain_un) fclose(out_plain_un);
    if (out_gz_id)    gzclose(out_gz_id);
    if (out_gz_un)    gzclose(out_gz_un);

    std::ostringstream temp_string;
    temp_string << std::endl << "[Identified]: " + std::to_string(identified_count) + "   [Unidentified]: " + std::to_string(unidentified_count);

    show_time(temp_string.str(), 0);

    auto end_time = std::chrono::system_clock::now();
    show_time("Done.", 0);
    
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    show_time("", 1, elapsed_seconds.count());
}
