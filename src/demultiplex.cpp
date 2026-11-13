
#define _CRT_SECURE_NO_WARNINGS
// [[Rcpp::plugins(cpp17)]]
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cstdio>         // For FILE*, fopen, fgets, fread, fclose, fwrite
#include <cstdlib>        // For getenv, exit
#include <cctype>
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




// ---------------------------------------
// FUNCTION: read primer & tag
// ---------------------------------------
struct SampleInfo {
    std::vector<std::string> exp, sample, ftags,rtags, forward_primer, reverse_primer;
     std::set<std::pair<std::string, std::string >> unique_primer_pair;
     std::vector<size_t> fwd_pos, rev_pos;
     std::vector<bool> try_reverse;
     std::vector<int> f_mismatch, r_mismatch;
};
struct hitInfo {
    std::string fwd, rev, ftag, rtag, exp, sample, direction;
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
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    show_time("Table use for demultiplexing was opened...");
    
    std::string line;
    char delimiter = '\t';  // tab delim was default

    // reading header
    if (!std::getline(file, line)) {
        CERR << "[ERROR] Table format was illegal!" << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
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
    std::vector<std::string> displayedColumnNames;
    while (std::getline(headerStream, header, delimiter)) {
        header.erase(0, header.find_first_not_of(" \t\r\n\""));
        header.erase(header.find_last_not_of(" \t\r\n\"") + 1);
        displayedColumnNames.push_back(header);

        std::string key;
        for (unsigned char c : header) {
            if (std::isalnum(c)) key.push_back(static_cast<char>(std::tolower(c)));
        }
        if (!key.empty() && key.back() == 's') key.pop_back();

        if (key == "exp" || key == "xexp" || key.find("experiment") != std::string::npos)
            columnNames.push_back("exp");
        else if (key.find("primer") != std::string::npos &&
                 (key.find("forward") != std::string::npos || key.find("fwd") != std::string::npos))
            columnNames.push_back("forward_primer");
        else if (key.find("primer") != std::string::npos &&
                 (key.find("reverse") != std::string::npos || key.find("rev") != std::string::npos))
            columnNames.push_back("reverse_primer");
        else if (key.find("tag") != std::string::npos)
            columnNames.push_back("reaction_tag");
        else if (key.find("reaction") != std::string::npos || key.find("sample") != std::string::npos)
            columnNames.push_back("reaction_id");
        else
            columnNames.push_back(key);
    }

    // build index map, seperate column names to index its own
    std::unordered_map<std::string, size_t> colIndex;
    for (size_t i = 0; i < columnNames.size(); ++i) {
        colIndex[columnNames[i]] = i;
    }
    const std::vector<std::string> required_cols = {
        "exp", "reaction_id", "forward_primer", "reverse_primer"
    };
    for (const auto& col : required_cols) {
        if (!colIndex.count(col)) {
            CERR << "[ERROR] Primer table missing required column: " << col << std::endl; // Check that the required columns exist
            return;
        }
    }

    // print debug info
    COUT << "Table colnames: ";
    for (const auto& col : displayedColumnNames) {
        COUT << "[" << col << "] ";
    }
    std::vector < std::pair < std::string, std::string >> primer_pair;
    COUT << std::endl;
    // read data line
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) continue;

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
            if (colIndex["exp"] >= row.size() ||
                colIndex["reaction_id"] >= row.size() ||
                colIndex["forward_primer"] >= row.size() ||
                colIndex["reverse_primer"] >= row.size()) {
                CERR << "[Warning] Skipping incomplete primer table row: " << line << std::endl;
                continue;
            }

            dem_infos.exp.push_back(row[colIndex["exp"]]);
            dem_infos.sample.push_back(row[colIndex["reaction_id"]]);

            if (colIndex.count("reaction_tag") && colIndex["reaction_tag"] < row.size()) {
                std::string temp_tag = row[colIndex["reaction_tag"]];
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
            else {
                dem_infos.ftags.push_back("");
                dem_infos.rtags.push_back("");
            }

            auto& fwd_primer = row[colIndex["forward_primer"]];
            std::transform(fwd_primer.begin(), fwd_primer.end(), fwd_primer.begin(),
                [](unsigned char c) { return std::toupper(c); });
            dem_infos.forward_primer.push_back(fwd_primer);

            auto& rev_primer = row[colIndex["reverse_primer"]];
            std::transform(rev_primer.begin(), rev_primer.end(), rev_primer.begin(),
                [](unsigned char c) { return std::toupper(c); });
            dem_infos.reverse_primer.push_back(rev_primer);

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
        COUT << "Pattern is longer than the target text." << std::endl;
        return {};
    }
    if (psl > 64)
    {
        COUT << "Pattern exceeds the supported length." << std::endl;
        return {};
    }
    std::vector<size_t> positions;
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
        std::vector <int> f_mismatch, r_mismatch;
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
            minfo.f_mismatch.push_back(hit_info.f_mismatch);
            minfo.r_mismatch.push_back(hit_info.r_mismatch);

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
                filtered_dem_infos.f_mismatch.push_back(minfo.f_mismatch[m]);
                filtered_dem_infos.r_mismatch.push_back(minfo.r_mismatch[m]);
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
                hit_info.direction = filtered_dem_infos.try_reverse[ind] ? "reverse" : "forward";
                hit_info.f_mismatch = filtered_dem_infos.f_mismatch[ind];
                hit_info.r_mismatch = filtered_dem_infos.r_mismatch[ind];
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
                hit_info.direction = filtered_dem_infos.try_reverse[ind] ? "reverse" : "forward";
                hit_info.f_mismatch = filtered_dem_infos.f_mismatch[ind];
                hit_info.r_mismatch = filtered_dem_infos.r_mismatch[ind];
                return temp_result;
            }
        }


    }
    return "";

}

static std::string append_demultiplex_header(const std::string& raw_id, const hitInfo& hit_info, bool with_tag) {
    std::string clean_id = raw_id;
    if (!clean_id.empty() && clean_id[0] == '@') {
        clean_id.erase(0, 1);
    }

    std::ostringstream block;
    block << "\"demultiplex\":{\"reaction_id\":\"" << hit_info.sample
        << "\",\"exp\":\"" << hit_info.exp
        << "\",\"direction\":\"" << hit_info.direction
        << "\",\"reaction_tag\":\"" << (with_tag ? hit_info.ftag + "+" + hit_info.rtag : "none")
        << "\",\"fwd\":\"" << hit_info.fwd
        << "\",\"rev\":\"" << hit_info.rev
        << "\",\"fwd_mismatch\":\"" << hit_info.f_mismatch
        << "\",\"rev_mismatch\":\"" << hit_info.r_mismatch
        << "\"}";

    size_t json_start = clean_id.find('{');
    if (json_start != std::string::npos && !clean_id.empty() && clean_id.back() == '}') {
        clean_id.pop_back();
        return ">" + clean_id + "," + block.str() + "}\n";
    }

    return ">" + clean_id + ";{" + block.str() + "}\n";
}

// ------------------------------
// Rcpp Exported Function: demultiplex
// ------------------------------
//' Demultiplex FASTQ reads using primer and tag information
//'
//' \code{demultiplex()} assigns each FASTQ read to a reaction based on matching forward
//' and reverse primers, as well as optional tag sequences. For each read, the
//' identified experiment, reaction ID, reaction tag, and primer mismatches are recorded in a
//' `"demultiplex"` JSON block appended to the FASTA header.
//'
//' Reaction tags are short oligonucleotide labels attached to primers before
//' PCR products are pooled for library preparation (tag-based multiplexing).
//' The supplied reaction table uses paired 8-bp tags. Tag length is determined
//' by the table; the function does not require 8-bp tags. Tags identify PCR
//' reactions, not different target-gene markers or sequencing-library indices.
//'
//' @param fastq_file (character) path to the input FASTQ file. The file may be
//'   plain text or gzip-compressed (suffix \code{.gz}).
//'
//' @param reaction_table (character) path to a CSV/TSV file containing primer and
//'   tag definitions. The columns are \code{experiment},
//'   \code{reaction_id}, \code{reaction_tag}, \code{forward_primer}, and
//'   \code{reverse_primer}. A paired tag is written as
//'   \code{forward_tag:reverse_tag}; a single tag is used for both directions.
//'
//' @param identified_output (character) path to the output file that receives
//'   reads successfully assigned to a reaction. If \code{compress_output = TRUE}
//'   or the file name ends with \code{.gz}, the file is gzip-compressed.
//'
//' @param unidentified_output (character) path to the output file for reads
//'   that cannot be matched to any reaction based on the specified primers and
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
void demultiplex( std::string fastq_file,  std::string reaction_table,
     std::string identified_output,  std::string unidentified_output,
    int n_threads = 4, bool with_tag = true,
    bool compress_output = false, int compress_level = 6,
    int primer_mismatch = 3)
{
    check_processing_parameters(n_threads, 1000000, compress_level);
    if (primer_mismatch < 0 || primer_mismatch > 63) throw std::invalid_argument("primer_mismatch must be between 0 and 63.");
    check_output_path(identified_output, {fastq_file, reaction_table, unidentified_output});
    check_output_path(unidentified_output, {fastq_file, reaction_table});
    int max_threads = detect_system_threads();
    if (n_threads > max_threads) n_threads = max_threads;
    size_t batch_size = 1000000;
    size_t total_count = 0;
    size_t identified_count = 0;
    size_t unidentified_count = 0;
    auto start_time = std::chrono::system_clock::now();
    auto now = std::chrono::system_clock::now();
    fastq_file =  expandTilde(fastq_file);
    reaction_table = expandTilde(reaction_table);
    identified_output = expandTilde(identified_output);
    unidentified_output = expandTilde(unidentified_output);

    // ------------------------------
    // 1. Read demutiplex table
    // ------------------------------
    SampleInfo infos;
    readPrimerFile(reaction_table, infos);
    if (infos.forward_primer.empty() || infos.reverse_primer.empty()) {
        CERR << "[Error] No valid primer definitions were loaded." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    if (with_tag) {
        for (size_t i = 0; i < infos.ftags.size(); ++i) {
            if (infos.ftags[i].empty() || infos.rtags[i].empty()) {
                CERR << "[Error] with_tag = TRUE requires a non-empty reaction_tag column for every primer-table row." << std::endl;
                throw std::runtime_error("Sequence processing failed; see the preceding error.");
            }
        }
    }

    // A read must map to exactly one table row. With tags enabled, repeated
    // primers are valid, but the complete primer/tag key must be unique. With
    // tags disabled, the primer pair itself is the assignment key.
    std::unordered_map<std::string, size_t> assignment_keys;
    for (size_t i = 0; i < infos.sample.size(); ++i) {
        std::string key = infos.forward_primer[i] + "\t" + infos.reverse_primer[i];
        if (with_tag) {
            key += "\t" + infos.ftags[i] + "\t" + infos.rtags[i];
        }
        auto inserted = assignment_keys.emplace(key, i);
        if (!inserted.second) {
            size_t previous = inserted.first->second;
            CERR << "[Error] Duplicate demultiplex assignment key for samples '"
                 << infos.sample[previous] << "' and '" << infos.sample[i] << "': "
                 << (with_tag ? "primer pair and tag pair are identical."
                              : "primer pair is identical while with_tag = FALSE.")
                 << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
    }

    // ------------------------------
    // 2. detect file type
    // ------------------------------
    bool is_gz_in = isGzFile(fastq_file);
    FILE* fp_plain = nullptr;
    gzFile fp_gz = nullptr;
    FileGuard input_guard(fp_plain, fp_gz);
    if (is_gz_in) {

        fp_gz = gzopen(fastq_file.c_str(), "rb");
    }
    else {
        fp_plain = fopen(fastq_file.c_str(), "rb");
    }
    if ((is_gz_in && fp_gz == nullptr) || (!is_gz_in && fp_plain == nullptr)) {
        CERR << "[Error] Cannot open input file: " << fastq_file << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
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
    std::string un_suffix = unidentified_output;
    std::transform(un_suffix.begin(), un_suffix.end(), un_suffix.begin(), ::tolower);
    const bool is_gz_un = compress_output || (un_suffix.size() >= 3 && un_suffix.substr(un_suffix.size() - 3) == ".gz");
    {
        auto lo = identified_output;
        std::transform(lo.begin(), lo.end(), lo.begin(), ::tolower);
        if ((lo.size() >= 3 && lo.substr(lo.size() - 3) == ".gz")||compress_output) {
            is_gz_out = true;
        }
    }
    FILE* out_plain_id = nullptr;
    FILE* out_plain_un = nullptr;
    gzFile out_gz_id = nullptr;
    gzFile out_gz_un = nullptr;
    FileGuard identified_guard(out_plain_id, out_gz_id), unidentified_guard(out_plain_un, out_gz_un);
    if (is_gz_out) {
        std::string mode = "wb" + std::to_string(compress_level);
        out_gz_id = gzopen(identified_output.c_str(), mode.c_str());
    }
    else {
        out_plain_id = fopen(identified_output.c_str(), "wb");
    }
    if (is_gz_un) {
        std::string mode = "wb" + std::to_string(compress_level);
        out_gz_un = gzopen(unidentified_output.c_str(), mode.c_str());
    } else {
        out_plain_un = fopen(unidentified_output.c_str(), "wb");
    }
    if ((is_gz_out ? !out_gz_id : !out_plain_id) ||
        (is_gz_un ? !out_gz_un : !out_plain_un)) {
        CERR << "[Error] Cannot open one or more output files." << std::endl;
        if (fp_plain) closeFile(fp_plain);
        if (fp_gz) closeGzip(fp_gz);
        if (out_plain_id) closeFile(out_plain_id);
        if (out_plain_un) closeFile(out_plain_un);
        if (out_gz_id) closeGzip(out_gz_id);
        if (out_gz_un) closeGzip(out_gz_un);
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
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

        std::vector<std::string> output_records(batch.size());
        std::vector<unsigned char> assigned(batch.size());
#ifdef USE_RCPP_OUTPUT
        Rcpp::checkUserInterrupt();
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
                        local = append_demultiplex_header(r.id, hit_info, with_tag)
                            + sub_seq + "\n";
                    }
                    else { 
                       
                        local = ">" + r.id.substr(1, r.id.size() - 1)
                            +  "\n"
                            + r.seq + "\n";
                    }

                    output_records[i] = std::move(local);
                    assigned[i] = matched;
                       
                    
                }
#ifdef _OPENMP
    }
#endif
        for (size_t i = 0; i < output_records.size(); ++i) {
            if (assigned[i]) ++identified_count; else ++unidentified_count;
            writeBuf(output_records[i], assigned[i] ? is_gz_out : is_gz_un, assigned[i] ? out_plain_id : out_plain_un,
                     assigned[i] ? out_gz_id : out_gz_un);
        }
        }

    
    
    // ------------------------------
    // 5. 清理
    // ------------------------------
    if (fp_plain)     closeFile(fp_plain);
    if (fp_gz)        closeGzip(fp_gz);
    if (out_plain_id) closeFile(out_plain_id);
    if (out_plain_un) closeFile(out_plain_un);
    if (out_gz_id)    closeGzip(out_gz_id);
    if (out_gz_un)    closeGzip(out_gz_un);

    std::ostringstream temp_string;
    temp_string << std::endl << "[Identified]: " + std::to_string(identified_count) + "   [Unidentified]: " + std::to_string(unidentified_count);

    show_time(temp_string.str(), 0);

    auto end_time = std::chrono::system_clock::now();
    show_time("Done.", 0);
    
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    show_time("", 1, elapsed_seconds.count());
}

//' Split tagged paired-end reads by PCR reaction
//'
//' Read paired FASTQ files in batches, assign each pair by its forward and
//' reverse tags, remove the matched tags and primers, and write one FASTQ pair
//' per reaction. This function only splits reads; it does not run DADA2.
//' Reaction tags are oligonucleotide labels used for tag-based multiplexing
//' before library pooling. The example tags are 8 bp long; other lengths are
//' accepted as specified in `reaction_table`.
//'
//' @param fastq_1 Path to the R1 FASTQ file.
//' @param fastq_2 Path to the R2 FASTQ file.
//' @param reaction_table CSV or tab-delimited reaction table containing
//'   `experiment`, `reaction_id`, `reaction_tag`, `forward_primer`, and
//'   `reverse_primer`. Write paired tags as `forward_tag:reverse_tag`.
//' @param output_dir Empty directory in which the split FASTQ files are written.
//' @param n_threads Number of threads used to match and trim each batch.
//' @param batch_size Number of paired reads processed in each batch.
//' @param compress_output Write gzip-compressed FASTQ files.
//' @param compress_level Gzip compression level from 0 to 9.
//' @param with_primer Remove the primers following the tags.
//'
//' @return `NULL`, invisibly. Split FASTQ files are written to `output_dir`.
//'
//' @details Tags are matched exactly. Primers may contain IUPAC bases but are
//'   matched without mismatches. A read pair is assigned only when R1 and R2
//'   match the same row in `reaction_table`.
//'
//' @examples
//' \dontrun{
//' barcode_split(
//'   fastq_1 = "run_R1.fastq.gz",
//'   fastq_2 = "run_R2.fastq.gz",
//'   reaction_table = "reactions.csv",
//'   output_dir = "split_fastq",
//'   n_threads = 8
//' )
//' }
//' @export
// [[Rcpp::export]]
void barcode_split(
    std::string fastq_1, std::string fastq_2, std::string reaction_table,
    std::string output_dir,
    int n_threads = 4, int batch_size = 100000, bool compress_output = true,
    int compress_level = 6, bool with_primer = true) {
    check_processing_parameters(n_threads, batch_size, compress_level);
    fastq_1 = expandTilde(fastq_1); fastq_2 = expandTilde(fastq_2);
    if (std::filesystem::weakly_canonical(fastq_1) == std::filesystem::weakly_canonical(fastq_2))
        throw std::invalid_argument("R1 and R2 must be different files.");
    output_dir = expandTilde(output_dir);
    if (std::filesystem::exists(output_dir) && !std::filesystem::is_empty(output_dir))
        throw std::invalid_argument("output_dir must be empty.");
    if (!std::filesystem::exists(output_dir) &&
        !std::filesystem::create_directories(output_dir))
        throw std::runtime_error("Cannot create output_dir.");
    const auto start_time = std::chrono::system_clock::now();
    const int max_threads = detect_system_threads();
    if (n_threads < 1) n_threads = 1;
    if (n_threads > max_threads) n_threads = max_threads;
    SampleInfo infos;
    readPrimerFile(expandTilde(reaction_table), infos);
    const size_t ns = infos.sample.size();
    if (ns == 0) {
        CERR << "[Error] No valid sample definitions were loaded." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    if (batch_size < 1) {
        CERR << "[Error] batch_size must be positive." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    show_time("Using " + std::to_string(n_threads) + " thread(s)");
    std::set<std::pair<std::string, std::string>> tag_pairs;
    std::set<std::pair<size_t, size_t>> tag_length_pairs;
    std::unordered_map<std::string, size_t> sample_by_tag;
    for (size_t i = 0; i < ns; ++i) {
        if (infos.ftags[i].empty() || infos.rtags[i].empty()) {
            CERR << "[Error] Every sample needs a non-empty paired tag."
                 << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
        if (!tag_pairs.insert({infos.ftags[i], infos.rtags[i]}).second) {
            CERR << "[Error] Tag pairs must be unique." << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
        tag_length_pairs.insert({infos.ftags[i].size(), infos.rtags[i].size()});
        sample_by_tag[infos.ftags[i] + "\t" + infos.rtags[i]] = i;
    }
    const std::vector<std::pair<size_t, size_t>> tag_lengths(
        tag_length_pairs.begin(), tag_length_pairs.end());

    std::set<std::string> unique_samples(infos.sample.begin(), infos.sample.end());
    const bool use_experiment = unique_samples.size() != ns;
    std::set<std::string> unique_file_names;
    std::vector<std::string> forward_output(ns);
    std::vector<std::string> reverse_output(ns);
    const std::string suffix = compress_output ? ".fastq.gz" : ".fastq";
    for (size_t i = 0; i < ns; ++i) {
        std::string file_name = use_experiment
            ? infos.exp[i] + "_" + infos.sample[i]
            : infos.sample[i];
        std::transform(file_name.begin(), file_name.end(), file_name.begin(),
            [](unsigned char c) {
                return (std::isalnum(c) || c == '_' || c == '-' || c == '.')
                    ? static_cast<char>(c) : '_';
            });
        if (!unique_file_names.insert(file_name).second) {
            CERR << "[Error] Sample names produce duplicate output file names."
                 << std::endl;
            throw std::runtime_error("Sequence processing failed; see the preceding error.");
        }
        forward_output[i] = output_dir + "/" + file_name + "_R1" + suffix;
        reverse_output[i] = output_dir + "/" + file_name + "_R2" + suffix;
    }

    const bool gz1 = isGzFile(fastq_1), gz2 = isGzFile(fastq_2);
    FILE *p1 = gz1 ? nullptr : fopen(fastq_1.c_str(), "rb");
    FILE *p2 = gz2 ? nullptr : fopen(fastq_2.c_str(), "rb");
    gzFile z1 = gz1 ? gzopen(fastq_1.c_str(), "rb") : nullptr;
    gzFile z2 = gz2 ? gzopen(fastq_2.c_str(), "rb") : nullptr;
    FileGuard guard1(p1, z1), guard2(p2, z2);
    if ((gz1 ? z1 == nullptr : p1 == nullptr) || (gz2 ? z2 == nullptr : p2 == nullptr)) {
        if (p1) closeFile(p1);
         if (p2) closeFile(p2);
        if (z1) closeGzip(z1);
         if (z2) closeGzip(z2);
        CERR << "[Error] Cannot open paired FASTQ input." << std::endl;
        throw std::runtime_error("Sequence processing failed; see the preceding error.");
    }
    std::vector<bool> written(ns, false);
    size_t total_count = 0;
    size_t assigned_count = 0;
    size_t unassigned_count = 0;
    bool failed = false;

    auto record_text = [](const FastqRecord& r) {
        return r.id + "\n" + r.seq + "\n" + r.plus + "\n" + r.qual + "\n";
    };
    // Write one batch of demultiplexed reads. The two vectors contain one
    // accumulated FASTQ string per sample: bf[i] for R1 and br[i] for R2.
    // `[&]` captures output paths, compression settings, and `written` by
    // reference, so the lambda updates the same state used by the main loop.
    auto flush = [&](std::vector<std::string>& bf,
                     std::vector<std::string>& br) -> bool {
        for (size_t i = 0; i < ns; ++i) {
            // No reads from this batch were assigned to sample i.
            if (bf[i].empty()) continue;

            // forward_output and reverse_output have the same sample order as
            // bf and br, so index i identifies both files for this sample.
            const std::string& pf = forward_output[i];
            const std::string& pr = reverse_output[i];
            if (compress_output) {
                // The first batch creates the gzip file (wb); later batches
                // append another valid gzip member (ab). Readers transparently
                // decode the concatenated members as a single FASTQ stream.
                std::string mode = std::string(written[i] ? "ab" : "wb") + std::to_string(compress_level);
                gzFile of = gzopen(pf.c_str(), mode.c_str());
                gzFile orv = gzopen(pr.c_str(), mode.c_str());
                if (!of || !orv) {
                    if (of) closeGzip(of);
                    if (orv) closeGzip(orv);
                    CERR << "[Error] Cannot open output files: "
                         << pf << " and " << pr << std::endl;
                    return false;
                }
                FILE *plain_f = nullptr, *plain_r = nullptr;
                FileGuard guard_f(plain_f, of), guard_r(plain_r, orv);
                writeBuf(bf[i], true, nullptr, of);
                writeBuf(br[i], true, nullptr, orv);
                closeGzip(of); closeGzip(orv);
            } else {
                // Plain FASTQ files use the equivalent binary write/append
                // modes. Binary mode prevents newline conversion on Windows.
                FILE* of = fopen(pf.c_str(), written[i] ? "ab" : "wb");
                FILE* orv = fopen(pr.c_str(), written[i] ? "ab" : "wb");
                if (!of || !orv) {
                    if (of) closeFile(of);
                     if (orv) closeFile(orv);
                    CERR << "[Error] Cannot open output files: "
                         << pf << " and " << pr << std::endl;
                    return false;
                }
                gzFile gzip_f = nullptr, gzip_r = nullptr;
                FileGuard guard_f(of, gzip_f), guard_r(orv, gzip_r);
                writeBuf(bf[i], false, of, nullptr);
                writeBuf(br[i], false, orv, nullptr);
                closeFile(of); closeFile(orv);
            }

            // Remember that subsequent batches must append, then release the
            // batch strings before the next group of reads is processed.
            written[i] = true;
            bf[i].clear(); br[i].clear();
        }
        return true;
    };

    // Process consecutive batches until both input files reach EOF.
    // `for (;;)` has no loop condition because EOF is detected by readRec()
    // inside the loop. A short final batch is written before the loop ends.
    for (;;) {
            // R1 and R2 are read together into equally sized vectors. Input
            // remains serial because gzip streams are not safe for concurrent
            // reads; matching and trimming below operate in parallel.
            std::vector<FastqRecord> batch_r1;
            std::vector<FastqRecord> batch_r2;
            batch_r1.reserve(static_cast<size_t>(batch_size));
            batch_r2.reserve(static_cast<size_t>(batch_size));
            bool reached_eof = false;
            for (int i = 0; i < batch_size; ++i) {
                FastqRecord r1, r2;
                const bool ok1 = readRec(gz1, p1, z1, r1, true);
                const bool ok2 = readRec(gz2, p2, z2, r2, true);
                if (!ok1 && !ok2) break;
                if (ok1 != ok2) {
                    CERR << std::endl
                         << "[Error] Paired FASTQ files have different record counts."
                         << std::endl;
                    failed = true;
                    break;
                }
                if (paired_read_id(r1.id) != paired_read_id(r2.id)) {
                    CERR << std::endl
                         << "[Error] R1 and R2 identifiers are out of sync. "
                         << "Check the FASTQ files." << std::endl;
                    failed = true;
                    break;
                }
                batch_r1.push_back(std::move(r1));
                batch_r2.push_back(std::move(r2));
            }
            if (failed) break;
            if (batch_r1.empty()) break;
            reached_eof = batch_r1.size() < static_cast<size_t>(batch_size);
            total_count += batch_r1.size();

            // Each index is owned by one OpenMP iteration. Threads therefore
            // modify separate records and result slots without locks.
            std::vector<size_t> matched_sample(batch_r1.size(), ns);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(n_threads)
#endif
            for (long long read_index = 0;
                 read_index < static_cast<long long>(batch_r1.size());
                 ++read_index) {
                FastqRecord& r1 = batch_r1[static_cast<size_t>(read_index)];
                FastqRecord& r2 = batch_r2[static_cast<size_t>(read_index)];
                std::transform(r1.seq.begin(), r1.seq.end(), r1.seq.begin(), ::toupper);
                std::transform(r2.seq.begin(), r2.seq.end(), r2.seq.begin(), ::toupper);
                // Identify the sample from the leading tag pair. Only the
                // distinct tag-length combinations are tested, so each read
                // is processed once regardless of the number of samples.
                size_t hit = ns;
                for (const auto& tag_length : tag_lengths) {
                    const size_t forward_tag_length = tag_length.first;
                    const size_t reverse_tag_length = tag_length.second;
                    if (r1.seq.size() < forward_tag_length ||
                        r2.seq.size() < reverse_tag_length) continue;
                    const std::string key =
                        r1.seq.substr(0, forward_tag_length) + "\t" +
                        r2.seq.substr(0, reverse_tag_length);
                    auto sample_it = sample_by_tag.find(key);
                    if (sample_it == sample_by_tag.end()) continue;

                    const size_t candidate = sample_it->second;
                    if (with_primer) {
                        const std::string forward_insert = r1.seq.substr(forward_tag_length);
                        const std::string reverse_insert = r2.seq.substr(reverse_tag_length);
                        auto forward_match = find_base_exact(
                            infos.forward_primer[candidate], forward_insert);
                        auto reverse_match = find_base_exact(
                            infos.reverse_primer[candidate], reverse_insert);
                        if (forward_match.empty() || forward_match[0] != 0 ||
                            reverse_match.empty() || reverse_match[0] != 0) continue;
                    }
                    if (hit != ns) {
                        hit = ns;
                        break;
                    }
                    hit = candidate;
                }
                if (hit == ns) {
                    continue;
                }
                size_t nf = infos.ftags[hit].size() + (with_primer ? infos.forward_primer[hit].size() : 0);
                size_t nr = infos.rtags[hit].size() + (with_primer ? infos.reverse_primer[hit].size() : 0);
                if (nf >= r1.seq.size() || nr >= r2.seq.size()) {
                    continue;
                }
                r1.seq.erase(0, nf); r1.qual.erase(0, nf);
                r2.seq.erase(0, nr); r2.qual.erase(0, nr);
                matched_sample[static_cast<size_t>(read_index)] = hit;
            }

            // Merge parallel results in original read order, then write once
            // per sample for this batch.
            std::vector<std::string> bf(ns), br(ns);
            for (size_t i = 0; i < matched_sample.size(); ++i) {
                const size_t hit = matched_sample[i];
                if (hit == ns) {
                    ++unassigned_count;
                    continue;
                }
                bf[hit] += record_text(batch_r1[i]);
                br[hit] += record_text(batch_r2[i]);
                ++assigned_count;
            }
            if (!flush(bf, br)) {
                failed = true;
                break;
            }
            show_time("Processed " + std::to_string(total_count) +
                      " paired reads...", 2);
            if (reached_eof) break;
#ifdef USE_RCPP_OUTPUT
            Rcpp::checkUserInterrupt();
#endif // check interrupt
    }
    if (p1) closeFile(p1);
    if (p2) closeFile(p2);
    if (z1) closeGzip(z1);
    if (z2) closeGzip(z2);
    if (failed) throw std::runtime_error("Paired FASTQ splitting failed; output is incomplete.");

    COUT << std::endl
         << "[Total pairs]: " << total_count
         << "   [Assigned]: " << assigned_count
         << "   [Unassigned]: " << unassigned_count
         << std::endl;
    const auto end_time = std::chrono::system_clock::now();
    const std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    show_time("", 1, elapsed_seconds.count());
}
