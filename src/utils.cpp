// [[Rcpp::plugins(cpp17)]]
#define _CRT_SECURE_NO_WARNINGS
#include "utils.hpp"
#include <string>
#include <vector>
//using namespace Rcpp;
//-------------------------------------------------------------------
// Expand "~" to the user's HOME directory if present in the path. 
//-------------------------------------------------------------------
std::string expandTilde(const std::string& path) {
    if (!path.empty() && path[0] == '~') {
        const char* home = getenv("HOME");
        if (home)
            return std::string(home) + path.substr(1);
    }
    return path;
}

// Check whether a file is gzipped by reading the first 2 bytes and comparing to the gzip signature.
bool isGzFile(const std::string& file_path, unsigned int buf_size) {
    std::string realPath = expandTilde(file_path);
    FILE* fp = fopen(realPath.c_str(), "rb");  // open file in binary mode
    if (!fp) return false;
    unsigned char header[2];
    size_t n = fread(header, 1, 2, fp);  // read first 2 bytes
    fclose(fp);
    return (n == 2 && header[0] == 0x1F && header[1] == 0x8B);
}
//-------------------------------------------------------------------
// Check whether a file is in FASTQ format by examining the first character.
//-------------------------------------------------------------------

bool isFastqFile(const std::string& filename) {
    bool gz = isGzFile(filename);
    char buffer[4096];

    if (gz) {
        gzFile f = gzopen(filename.c_str(), "rb");
        if (!f) return false;
        while (gzgets(f, buffer, sizeof(buffer))) {
            if (buffer[0] == '\n' || buffer[0] == '\r') continue;
            gzclose(f);
            return buffer[0] == '@';
        }
        gzclose(f);
        return false;
    }
    else {
        FILE* f = fopen(filename.c_str(), "rb");
        if (!f) return false;
        while (fgets(buffer, sizeof(buffer), f)) {
            if (buffer[0] == '\n' || buffer[0] == '\r') continue;
            fclose(f);
            return buffer[0] == '@';
        }
        fclose(f);
        return false;
    }
}
// store base probability----------------------------------


//-------------------------------------------------------------------
// Remove newline and carriage return characters from the end of a C string 
//-------------------------------------------------------------------
void trimEnd(char* s) {
    size_t len = strlen(s);
    while (len > 0 && (s[len - 1] == '\n' || s[len - 1] == '\r'))
        s[--len] = '\0';
}

//-------------------------------------------------------------------
// Enhanced records reading function [text file]
//-------------------------------------------------------------------
bool readRecordPlain(FILE* fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq) {
    if (fastq) {
        // Read header line
        if (!fgets(buffer, buf_size, fp)) return false;
        trimEnd(buffer);
        rec.id = buffer;
        // Read sequence line
        if (!fgets(buffer, buf_size, fp)) return false;
        trimEnd(buffer);
        rec.seq = buffer;
        // Plus
        if (!fgets(buffer, buf_size, fp))  return false;
        trimEnd(buffer);
        rec.plus = buffer;
        // Quality
        if (!fgets(buffer, buf_size, fp))  return false;
        trimEnd(buffer);
        rec.qual = buffer;
        return true;
    }
    else {
    // Read id
        while (fgets(buffer, buf_size, fp)) {
            if (buffer[0] == '>') {
                trimEnd(buffer);
                rec.id = buffer;
                break;
            }
        }
        if (rec.id.empty()) { return false; }
    // Read sequence
        std::ostringstream temp_seq;
        while (fgets(buffer, buf_size, fp)) {
            size_t read_len = strlen(buffer);
            if (buffer[0] == '>') {
                fseek(fp, -static_cast<long>(read_len), SEEK_CUR);
                break;
            }
            trimEnd(buffer);
            temp_seq << buffer;
        }
        rec.seq = temp_seq.str();
        if (rec.seq.empty()) { return false; }
    }
    return true;
   
}

//-------------------------------------------------------------------
// Enhanced records reading function [gz file]
//-------------------------------------------------------------------

bool readRecordGZ(gzFile fp, FastqRecord& rec, char* buffer, int buf_size, bool fastq) {
    if (fastq) {
        // Read header line 
        if (!gzgets(fp, buffer, buf_size)) return false;
        trimEnd(buffer);
        rec.id = buffer;

        // Read sequence line
        if (!gzgets(fp, buffer, buf_size)) return false;
        trimEnd(buffer);
        rec.seq = buffer;

        // Plus
        if (!gzgets(fp, buffer, buf_size)) return false;
        trimEnd(buffer);
        rec.plus = buffer;
        // Quality
        if (!gzgets(fp, buffer, buf_size)) return false;
        trimEnd(buffer);
        rec.qual = buffer;
        return true;
    }
    else {
        while (gzgets(fp, buffer, buf_size)) {
            if (buffer[0] == '>') {
                trimEnd(buffer);
                rec.id = buffer;
                break;
            }
        }
        if (rec.id.empty()) { return false; }
        std::ostringstream temp_seq;
        while (gzgets(fp, buffer, buf_size)) {
            size_t read_length = strlen(buffer);
            if (buffer[0] == '>') {
                gzseek(fp, -static_cast<long>(read_length), SEEK_CUR);
                break;
            }
            trimEnd(buffer);
            temp_seq << buffer;
        }
        rec.seq = temp_seq.str();
        if (rec.seq.empty()) { return false; }
    }
    return true;
}

//-------------------------------------------------------------------
// Read function, auto detect file format
//-------------------------------------------------------------------
/**
* @param is_gz_in Is it a gz file or not.
* @param fp_plain Pointer, for text file.
* @param rec FastqRecord structrue.
* @param fastq Is it fastq.
* @param buffer_size The buffer during reading. If reading a genomic file, make it larger.
*/
bool readRec(bool is_gz_in, FILE* fp_plain, gzFile fp_gz, FastqRecord& rec, bool fastq, int buffer_size)
{

    static thread_local std::vector<char> buffer(buffer_size);

    if (is_gz_in)
        return readRecordGZ(fp_gz, rec, buffer.data(), buffer_size, fastq);
    else
        return readRecordPlain(fp_plain, rec, buffer.data(), buffer_size, fastq);
}

//-------------------------------------------------------------------
//  Write function, auto detect format
//-------------------------------------------------------------------
/**
* @param buf The buffer that needs to be written.
* @param is_gz_out If true, write our the gz file.
* @param r FastqRecord structrue.
* @param fastq Is it fastq.
* @param buffer_size The buffer during reading. If reading a genomic file, make it larger.
*/
void writeBuf(const std::string& buf, bool is_gz_out, FILE* out_plain, gzFile out_gz)
{
    const char* data = buf.data();
    size_t total = buf.size();
    const size_t CHUNK = std::numeric_limits<unsigned int>::max();

    while (total > 0) {
        unsigned int nbytes = static_cast<unsigned int>(std::min(total, CHUNK));
        if (is_gz_out)
            gzwrite(out_gz, data, nbytes);
        else
            fwrite(data, 1, nbytes, out_plain);
        data += nbytes;
        total -= nbytes;
    }
}

//-------------------------------------------------------------------
// Build a quality score lookup table for Phred Q values from 0 to 93.
// The table is used to compute weights during overlap scoring.
// Values in table represent the probability of a base being correct, calculated as 1 - 10^(-Q/10).
//-------------------------------------------------------------------
std::vector<double> build_quality_table() {
    std::vector<double> table(94);
    for (int q = 0; q <= 93; ++q)
        table[q] = 1.0 - pow(10.0, -q / 10.0); 
    return table; // Correct probability
}


// Global quality weight table (constant after initialization)
extern const std::vector<double> QUALITY_WEIGHT_TABLE = build_quality_table();

//-------------------------------------------------------------------
// Overlap Detection and Scoring Functions
//-------------------------------------------------------------------
// find_overlap_with_mismatches scans for an overlap between the end of s1 and the beginning of s2.
// It iterates from the maximum possible overlap down to a given minimum.
// For each candidate overlap length L, it counts mismatches and records their positions (up to 60 mismatches).
// Returns the overlap length if the mismatch count is within the allowed limit.
//-------------------------------------------------------------------


inline double cal_pr_diff(
    char x, char y,
    double correct_x, double correct_y,
    double error_x, double error_y,
    const base_probability& pb
) {
    auto get = [&](char b) -> double {
        switch (b) {
        case 'A': return pb.A;
        case 'C': return pb.C;
        case 'G': return pb.G;
        case 'T': return pb.T;
        default: return 0.0;
        }
    };

    double px = get(x);
    double py = get(y);

    // Σ_{b≠Xi}
    double sum_not_x = 1.0 - px;

    // Σ_{b≠Yi}
    double sum_not_y = 1.0 - py;

    // Σ_{b≠Xi,Yi}
    double sum_not_xy = 1.0 - px - py;

    // Σ P_b^2 over b≠Xi,Yi
    double sum_sq_not_xy = 0.0;
    for (char b : {'A', 'C', 'G', 'T'}) {
        if (b != x && b != y) {
            double p = get(b);
            sum_sq_not_xy += p * p;
        }
    }

    double term1 = correct_y * error_x * (py / sum_not_x);
    double term2 = correct_x * error_y * (px / sum_not_y);
    double term3 = error_x * error_y *
        (sum_sq_not_xy / (sum_not_xy * sum_not_xy));

    return term1 + term2 + term3;
}


// ------------------------------
// ------------------------------
double oes_test(double q, int omega, double oes, double alpha, double beta, int max_length) {
	if (omega > max_length) return 0.0;
	if (q <= 0.0) return 1.0;
	if (q >= 1.0) return 0.0;
	if (alpha == beta) return 1.0;

	const double one_minus_q = 1.0 - q;
	const double ratio = q / one_minus_q;
	long double p = 1.0;

	for (int i = omega; i <= max_length; ++i) {
		int l_c = static_cast<int>(std::ceil((oes - beta * i) / (alpha - beta))) - 1;

		if (l_c < 0) {
			return 1.0;
		}
		if (l_c >= i) {
			continue;
		}

		// term = C(i,0) * q^0 * (1-q)^i
		long double term = std::pow(one_minus_q, static_cast<double>(i));
		long double sigma = term;

		// 
		// T(k+1) = T(k) * ((i-k)/(k+1)) * (q/(1-q))
		for (int k = 0; k < l_c; ++k) {
			term *= (static_cast<double>(i - k) / static_cast<double>(k + 1)) * ratio;
			sigma += term;
		}

		p *= sigma;
		if (p <= 0.0) {
			return 1.0;
		}
	}
    double return_p;
    return_p = (double)(1.0 - p * p) < 0 ? 0: (double)(1.0 - p * p);
    return return_p;
}

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
    bool enable_oes_test,
    bool MAP
) {
    int r1_size = static_cast<int>(r1.seq.size());
    int r2_size = static_cast<int>(r2.seq.size());

    std::string r2_reverse;
    std::string r2_q_reverse = r2.qual;
    fast_reverse_complement(r2.seq, r2_reverse);
    std::reverse(r2_q_reverse.begin(), r2_q_reverse.end());

    std::string_view r1_seq_v(r1.seq);
    std::string_view r2_seq_v(r2_reverse);
    std::string_view r1_q_v(r1.qual);
    std::string_view r2_q_v(r2_q_reverse);

    merged_rec.seq.clear();
    merged_rec.qual.clear();
    merged_rec.id = r1.id;
    merged_rec.plus = r1.plus;
    /*
    // Non-OES mode: use the original overlap-by-mismatch strategy.
    if (!enable_oes_test) {
        bool is_offset = false;
        size_t mismatch_positions[1000];
        mismatch_count = 0;
        p_value = 0.0;

        size_t ovl = find_overlap_with_mismatches(
            r1.seq, r2_reverse, is_offset,
            min_overlap, max_mismatches,
            mismatch_count, mismatch_positions
        );

        overlap_length = ovl;
        overlap_score = 0.0;

        if (ovl > 0) {
            std::string mergedSeq, mergedQual;
            double score = 0.0;

            if (mismatch_count > 0) {
                std::string prefixSeq, prefixQl, overlapSeq, overlapQl;

                if (is_offset) {
                    prefixSeq = r2_reverse.substr(0, r2_reverse.size() - ovl);
                    prefixQl = r2_q_reverse.substr(0, r2_q_reverse.size() - ovl);
                    overlapSeq = r2_reverse.substr(r2_reverse.size() - ovl, ovl);
                    overlapQl = r2_q_reverse.substr(r2_q_reverse.size() - ovl, ovl);

                    for (size_t p = 0; p < mismatch_count; ++p) {
                        char b2 = r1.seq[mismatch_positions[p]];
                        char q1 = overlapQl[mismatch_positions[p]];
                        char q2 = r1.qual[mismatch_positions[p]];
                        if (q2 > q1) {
                            overlapSeq[mismatch_positions[p]] = b2;
                            overlapQl[mismatch_positions[p]] = q2;
                        }
                    }

                    score = calculate_overlap_score(
                        r2_reverse, r2_q_reverse,
                        r1.seq, r1.qual,
                        r2_reverse.size() - ovl,
                        ovl
                    );

                    mergedSeq = prefixSeq + overlapSeq + r1.seq.substr(ovl);
                    mergedQual = prefixQl + overlapQl + r1.qual.substr(ovl);
                }
                else {
                    prefixSeq = r1.seq.substr(0, r1.seq.size() - ovl);
                    prefixQl = r1.qual.substr(0, r1.qual.size() - ovl);
                    overlapSeq = r1.seq.substr(r1.seq.size() - ovl, ovl);
                    overlapQl = r1.qual.substr(r1.qual.size() - ovl, ovl);

                    for (size_t p = 0; p < mismatch_count; ++p) {
                        char b2 = r2_reverse[mismatch_positions[p]];
                        char q1 = overlapQl[mismatch_positions[p]];
                        char q2 = r2_q_reverse[mismatch_positions[p]];
                        if (q2 > q1) {
                            overlapSeq[mismatch_positions[p]] = b2;
                            overlapQl[mismatch_positions[p]] = q2;
                        }
                    }

                    score = calculate_overlap_score(
                        r1.seq, r1.qual,
                        r2_reverse, r2_q_reverse,
                        r1.seq.size() - ovl,
                        ovl
                    );

                    mergedSeq = prefixSeq + overlapSeq + r2_reverse.substr(ovl);
                    mergedQual = prefixQl + overlapQl + r2_q_reverse.substr(ovl);
                }
            }
            else {
                if (is_offset) {
                    mergedSeq = r2_reverse + r1.seq.substr(ovl);
                    mergedQual = r2_q_reverse + r1.qual.substr(ovl);
                }
                else {
                    mergedSeq = r1.seq + r2_reverse.substr(ovl);
                    mergedQual = r1.qual + r2_q_reverse.substr(ovl);
                }
                score = 1.0;
            }

            merged_rec.seq = std::move(mergedSeq);
            merged_rec.qual = std::move(mergedQual);
            overlap_score = score;
            overlap_length = ovl;
        }
        else {
            merged_rec.seq.reserve(r1.seq.size() + r2_reverse.size());
            merged_rec.qual.reserve(r1.qual.size() + r2_q_reverse.size());
            merged_rec.seq += r1.seq;
            merged_rec.seq += r2_reverse;
            merged_rec.qual += r1.qual;
            merged_rec.qual += r2_q_reverse;

            overlap_score = 0.0;
            overlap_length = 0;
            mismatch_count = 0;
        }

        return;
    }
    */
    // OES mode: use the PEAR-style scoring and OES filtering.
    double q = base_pb.A * base_pb.A +
        base_pb.C * base_pb.C +
        base_pb.G * base_pb.G +
        base_pb.T * base_pb.T;

    bool found_best = false;
    double best_score = std::numeric_limits<double>::lowest();
    double best_oes_score = 0.0;
    size_t best_mismatch = 0;
    size_t best_length = 0;
    int best_offset = 0;
    size_t best_r1_start = 0;
    size_t best_r2_start = 0;
    std::string best_overlap_region;
    std::string best_overlap_region_q;

    // ofs is the start position of r2 relative to r1.
    int min_offset = static_cast<int>(min_overlap) - r2_size;
    int max_offset = r1_size - static_cast<int>(min_overlap);

    for (int ofs = min_offset; ofs <= max_offset; ++ofs) {
        int ovl_begin = std::max(0, ofs);
        int ovl_end = std::min(r1_size, ofs + r2_size);
        int ovl_len_int = ovl_end - ovl_begin;

        if (ovl_len_int < static_cast<int>(min_overlap)) {
            continue;
        }

        size_t current_overlap_length = static_cast<size_t>(ovl_len_int);
        size_t r1_start = static_cast<size_t>(ovl_begin);
        size_t r2_start = static_cast<size_t>(ovl_begin - ofs);

        std::string_view r1_ovl = r1_seq_v.substr(r1_start, current_overlap_length);
        std::string_view r1_ovl_q = r1_q_v.substr(r1_start, current_overlap_length);
        std::string_view r2_ovl = r2_seq_v.substr(r2_start, current_overlap_length);
        std::string_view r2_ovl_q = r2_q_v.substr(r2_start, current_overlap_length);

        std::string overlap_region;
        std::string overlap_region_q;
        overlap_region.reserve(current_overlap_length);
        overlap_region_q.reserve(current_overlap_length);

        size_t current_mismatch = 0;
        double ovl_score = 0.0;
        double oes_score = 0.0;

        // Important: do not stop early on mismatch threshold here.
        // OES mode should evaluate the full overlap, otherwise the score is incomplete.
        for (size_t o = 0; o < current_overlap_length; ++o) {
            const char b1 = r1_ovl[o];
            const char b2 = r2_ovl[o];
            const char q1 = r1_ovl_q[o];
            const char q2 = r2_ovl_q[o];

            if (b1 != b2) {
                ++current_mismatch;
            }

            double Pr_same = 0.0;
            double Pr_diff = 0.0;
            double correct_x = QUALITY_WEIGHT_TABLE[q1 - 33];
            double correct_y = QUALITY_WEIGHT_TABLE[q2 - 33];
            double error_x = 1.0 - correct_x;
            double error_y = 1.0 - correct_y;

            if ((b1 == b2) && (b1 != 'N')) {
                overlap_region.push_back(b2);
                overlap_region_q.push_back(correct_x < correct_y ? q2 : q1);

                switch (b1) {
                case 'A': Pr_same = correct_x * correct_y + error_x * error_y * base_pb.A_sqr; break;
                case 'C': Pr_same = correct_x * correct_y + error_x * error_y * base_pb.C_sqr; break;
                case 'T': Pr_same = correct_x * correct_y + error_x * error_y * base_pb.T_sqr; break;
                case 'G': Pr_same = correct_x * correct_y + error_x * error_y * base_pb.G_sqr; break;
                default:  Pr_same = correct_x * correct_y; break;
                }

                oes_score += Pr_same * alpha + (1.0 - Pr_same) * beta;
                ovl_score += Pr_same * alpha;
            }
            else if (b1 == 'N' || b2 == 'N') {
                if (b1 != b2) {
                    overlap_region.push_back(b1 == 'N' ? b2 : b1);
                    overlap_region_q.push_back(b1 == 'N' ? q2 : q1);

                    Pr_diff = 1.0 - q;
                    oes_score += Pr_diff * alpha + (1.0 - Pr_diff) * beta;
                    ovl_score += (1.0 - Pr_diff) * beta;
                }
                else {
                    Pr_same = q;
                    overlap_region.push_back('N');
                    overlap_region_q.push_back(correct_x < correct_y ? q1 : q2);

                    oes_score += Pr_same * alpha + (1.0 - Pr_same) * beta;
                    ovl_score += Pr_same * alpha;
                }
            }
            else {
                overlap_region.push_back(correct_x >= correct_y ? b1 : b2);
                overlap_region_q.push_back(correct_x >= correct_y ? q1 : q2);

                Pr_diff = cal_pr_diff(b1, b2, correct_x, correct_y, error_x, error_y, base_pb);
                oes_score += Pr_diff * alpha + (1.0 - Pr_diff) * beta;
                ovl_score += (1.0 - Pr_diff) * beta;
            }
        }

        if (!found_best ||
            ovl_score > best_score ||
            (ovl_score == best_score && current_overlap_length > best_length) ||
            (ovl_score == best_score && current_overlap_length == best_length && current_mismatch < best_mismatch)) {
            found_best = true;
            best_score = ovl_score;
            best_oes_score = oes_score;
            best_mismatch = current_mismatch;
            best_length = current_overlap_length;
            best_offset = ofs;
            best_r1_start = r1_start;
            best_r2_start = r2_start;
            best_overlap_region = overlap_region;
            best_overlap_region_q = overlap_region_q;
        }
    }

    if (!found_best) {
        p_value = 1.0;
        overlap_score = 0.0;
        overlap_length = 0;
        mismatch_count = 0;
        merged_rec.seq.reserve(r1_seq_v.size() + r2_seq_v.size());
        merged_rec.qual.reserve(r1_q_v.size() + r2_q_v.size());
        merged_rec.seq += r1_seq_v;
        merged_rec.seq += r2_seq_v;
        merged_rec.qual += r1_q_v;
        merged_rec.qual += r2_q_v;
        return;
    }

    p_value = oes_test(
        q,
        MAP ? static_cast<int>(best_length) : static_cast<int>(min_overlap),
        best_oes_score,
        alpha,
        beta,
        base_pb.max_length
    );

    if (p_value < 0.0) p_value = 0.0;
    if (p_value > 1.0) p_value = 1.0;

    if (p_value >= 0.01) {
        overlap_score = 0.0;
        overlap_length = 0;
        mismatch_count = 0;

        merged_rec.seq.reserve(r1_seq_v.size() + r2_seq_v.size());
        merged_rec.qual.reserve(r1_q_v.size() + r2_q_v.size());
        merged_rec.seq += r1_seq_v;
        merged_rec.seq += r2_seq_v;
        merged_rec.qual += r1_q_v;
        merged_rec.qual += r2_q_v;
        return;
    }

    merged_rec.seq.reserve(r1_seq_v.size() + r2_seq_v.size() - best_length);
    merged_rec.qual.reserve(r1_q_v.size() + r2_q_v.size() - best_length);

    // Left non-overlapping segment.
    if (best_offset >= 0) {
        merged_rec.seq += r1_seq_v.substr(0, best_r1_start);
        merged_rec.qual += r1_q_v.substr(0, best_r1_start);
    }
    else {
        merged_rec.seq += r2_seq_v.substr(0, best_r2_start);
        merged_rec.qual += r2_q_v.substr(0, best_r2_start);
    }

    // Merged overlap segment.
    merged_rec.seq += best_overlap_region;
    merged_rec.qual += best_overlap_region_q;

    // Right non-overlapping tail.
    int r1_end_pos = r1_size;
    int r2_end_pos = best_offset + r2_size;
    size_t r1_tail_start = best_r1_start + best_length;
    size_t r2_tail_start = best_r2_start + best_length;

    if (r1_end_pos >= r2_end_pos) {
        merged_rec.seq += r1_seq_v.substr(r1_tail_start);
        merged_rec.qual += r1_q_v.substr(r1_tail_start);
    }
    else {
        merged_rec.seq += r2_seq_v.substr(r2_tail_start);
        merged_rec.qual += r2_q_v.substr(r2_tail_start);
    }

    overlap_score = best_score;
    overlap_length = best_length;
    mismatch_count = best_mismatch;
}


//-------------------------------------------------------------------
// Score calculation Function
//-------------------------------------------------------------------
double calculate_overlap_score(const std::string& s1, const std::string& q1,
    const std::string& s2, const std::string& q2,
    size_t overlap_start, size_t overlap_length, double lambda) {
    double total_weight = 0.0, penalty_weight = 0.0;
    // Iterate over the overlap region.
    for (size_t i = 0; i < overlap_length; i++) {
        int Q1 = q1[overlap_start + i] - 33; // Convert quality char to integer.
        int Q2 = q2[i] - 33;
        // Lookup weights from the precomputed quality table.
        double w1 = (Q1 >= 0 && Q1 < 94) ? QUALITY_WEIGHT_TABLE[Q1] : 0.0;
        double w2 = (Q2 >= 0 && Q2 < 94) ? QUALITY_WEIGHT_TABLE[Q2] : 0.0;
        double w = (w1 + w2) * 0.5;  // Average the two quality weights.
        total_weight += w;
        // If the bases do not match, add penalty weighted by lambda.
        if (s1[overlap_start + i] != s2[i]) {
            if (lambda == 0) {
                penalty_weight = w;
            }
            else { penalty_weight += lambda * w; }
        }
            
    }
    // Return the normalized overlap score.
    return total_weight > 0.0 ? (1.0 - penalty_weight / total_weight) : 0.0;
}



//------------------------------------------------------------------
// Show time function
//------------------------------------------------------------------
/*
 * @param s    Message text to display.
 * @param mode Display mode selector:
 *              - 0 : Log message with current time.
 *              - 1 : Show elapsed time summary.
 *              - 2 : Inline updating progress line.
 * @param secs Elapsed seconds (used only when mode = 1).
*/
void show_time(const std::string& s, size_t mode, double secs) {
    // now time
    auto now = std::chrono::system_clock::now();
    std::time_t now_t = std::chrono::system_clock::to_time_t(now);
    std::tm* local_time = std::localtime(&now_t);
    if (mode == 0) {
        // formatting time
        COUT << '[' << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << "] ";
        COUT << s << std::endl;
    }
    else if (mode == 1) {
        COUT << std::fixed << std::setprecision(2);
        COUT << "Elapsed time: ";
        if (secs < 60.0) {
            COUT << secs << " seconds";
        }
        else if (secs < 3600.0) {
            double mins = secs / 60.0;
            COUT << mins << " minutes";
        }
        else if (secs < 86400.0) {
            double hours = secs / 3600.0;
            COUT << hours << " hours";
        }
        else {
            double days = secs / 86400.0;
            COUT << days << " days";
        }
        COUT << std::endl << "----------------------------------------------------------------" << std::endl;
    }
    else if (mode == 2) {
        COUT << "\r" << '[' << std::put_time(local_time, "%Y-%m-%d %H:%M:%S") << "] " <<  s << std::flush;
    }
}
// ------------------------------
// Optimized Reverse Complement Functions
// ------------------------------
// get_comp_table returns a constant lookup table mapping nucleotide characters to their complements.
// This is defined using std::array and a lambda to initialize it once.
const std::array<unsigned char, 256>& get_comp_table() {
    static const std::array<unsigned char, 256> table = []() {
        std::array<unsigned char, 256> t{};

        // 标准碱基互补
        t['A'] = 'T'; t['T'] = 'A'; t['C'] = 'G'; t['G'] = 'C';
        t['a'] = 't'; t['t'] = 'a'; t['c'] = 'g'; t['g'] = 'c';

        // 简并碱基互补（大写）
        t['R'] = 'Y'; t['Y'] = 'R';
        t['S'] = 'S'; t['W'] = 'W';
        t['K'] = 'M'; t['M'] = 'K';
        t['B'] = 'V'; t['D'] = 'H';
        t['H'] = 'D'; t['V'] = 'B';
        t['N'] = 'N';

        // 简并碱基互补（小写）
        t['r'] = 'y'; t['y'] = 'r';
        t['s'] = 's'; t['w'] = 'w';
        t['k'] = 'm'; t['m'] = 'k';
        t['b'] = 'v'; t['d'] = 'h';
        t['h'] = 'd'; t['v'] = 'b';
        t['n'] = 'n';

        return t;
        }();
    return table;
}
// fast_reverse_complement computes the reverse complement of the input sequence.
// It pre-allocates the destination string to the proper size and iterates over the input from end to start.
// For each character, it uses the lookup table to get the complement if available.
void fast_reverse_complement(const std::string& src, std::string& dest) {
    // Resize destination to match source length.
    dest.resize(src.size());
    // Retrieve the prebuilt complement lookup table.
    const auto& comp_table = get_comp_table();
    size_t n = src.size();
    // Loop over each character in reverse order.
    for (size_t i = 0; i < n; i++) {
        unsigned char c = src[n - 1 - i];  // Read character from the end.
        // Use lookup table; if no complement is available (table[c] is 0), use the original character.
        dest[i] = comp_table[c] ? comp_table[c] : c;
    }
}
//----------------------------------------------------------------
int detect_system_threads() {
    int system_threads = 1;

#ifdef _OPENMP
    // 
#pragma omp parallel
    {
#pragma omp single
        system_threads = omp_get_num_threads();
    }
#endif

    return system_threads;
}
//---------------------------------------------------------------
// Base probability


 base_probability calculate_base_probabilities(const std::vector <FastqRecord>& rec_v) {
    base_probability base_pb{0.0, 0.0, 0.0, 0.0, 0.0};
    double total_length = 0;
    size_t freq[256] = { 0 };
    int max_length = 0;
    for (const FastqRecord &rec:  rec_v) {
        int temp_length = 0;
        std::string_view sv(rec.seq);
        for (unsigned char c : sv) {
            temp_length++;
            ++freq[c];
            total_length++;
        }
		max_length = temp_length > max_length ? temp_length : max_length;
    }
    
    base_pb.A = static_cast<double>(freq['A'] + freq['a']) / total_length;
    base_pb.T = static_cast<double>(freq['T'] + freq['t']) / total_length;
    base_pb.C = static_cast<double>(freq['C'] + freq['c']) / total_length;
    base_pb.G = static_cast<double>(freq['G'] + freq['g']) / total_length;
    base_pb.N = static_cast<double>(freq['N'] + freq['n']) / total_length;
	
    // ------
    base_pb.A_sqr = (base_pb.A * base_pb.A) / ((base_pb.T + base_pb.C + base_pb.G) * (base_pb.T + base_pb.C + base_pb.G));
    base_pb.T_sqr = (base_pb.T * base_pb.T) / ((base_pb.A + base_pb.C + base_pb.G) * (base_pb.A + base_pb.C + base_pb.G));
    base_pb.C_sqr = (base_pb.C * base_pb.C) / ((base_pb.A + base_pb.T + base_pb.G) * (base_pb.A + base_pb.T + base_pb.G));
    base_pb.G_sqr = (base_pb.G * base_pb.G) / ((base_pb.A + base_pb.C + base_pb.T) * (base_pb.A + base_pb.C + base_pb.T));
    //
    base_pb.max_length = max_length;
	return base_pb;
}