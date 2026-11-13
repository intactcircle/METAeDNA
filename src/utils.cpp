// [[Rcpp::plugins(cpp17)]]
#define _CRT_SECURE_NO_WARNINGS
#include "utils.hpp"
#include <string>
#include <vector>
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
    rec = FastqRecord{};
    if (fastq) {
        // Read header line
        if (!fgets(buffer, buf_size, fp)) {
            if (ferror(fp)) throw std::runtime_error("Failed to read FASTQ.");
            return false;
        }
        trimEnd(buffer);
        rec.id = buffer;
        // Read sequence line
        if (!fgets(buffer, buf_size, fp)) throw std::runtime_error("Truncated FASTQ sequence: " + rec.id);
        trimEnd(buffer);
        rec.seq = buffer;
        // Plus
        if (!fgets(buffer, buf_size, fp)) throw std::runtime_error("Truncated FASTQ separator: " + rec.id);
        trimEnd(buffer);
        rec.plus = buffer;
        // Quality
        if (!fgets(buffer, buf_size, fp)) throw std::runtime_error("Truncated FASTQ quality: " + rec.id);
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
    rec = FastqRecord{};
    if (fastq) {
        // Read header line 
        if (!gzgets(fp, buffer, buf_size)) {
            if (!gzeof(fp)) throw std::runtime_error("Failed to read gzip FASTQ.");
            return false;
        }
        trimEnd(buffer);
        rec.id = buffer;

        // Read sequence line
        if (!gzgets(fp, buffer, buf_size)) throw std::runtime_error("Truncated gzip FASTQ sequence: " + rec.id);
        trimEnd(buffer);
        rec.seq = buffer;

        // Plus
        if (!gzgets(fp, buffer, buf_size)) throw std::runtime_error("Truncated gzip FASTQ separator: " + rec.id);
        trimEnd(buffer);
        rec.plus = buffer;
        // Quality
        if (!gzgets(fp, buffer, buf_size)) throw std::runtime_error("Truncated gzip FASTQ quality: " + rec.id);
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

    if (buffer_size < 2) throw std::invalid_argument("buffer_size must be at least 2.");
    static thread_local std::vector<char> buffer;
    buffer.resize(buffer_size);
    bool found = is_gz_in
        ? readRecordGZ(fp_gz, rec, buffer.data(), buffer_size, fastq)
        : readRecordPlain(fp_plain, rec, buffer.data(), buffer_size, fastq);
    if (is_gz_in) {
        int error = Z_OK;
        gzerror(fp_gz, &error);
        if (error != Z_OK && error != Z_STREAM_END) throw std::runtime_error("Corrupt or truncated gzip input.");
    } else if (ferror(fp_plain)) throw std::runtime_error("Input read failed.");
    if (found && fastq) {
        if (rec.id.empty() || rec.id[0] != '@' || rec.plus.empty() || rec.plus[0] != '+' ||
            rec.seq.empty() || rec.seq.size() != rec.qual.size()) {
            throw std::runtime_error("Invalid FASTQ record or line exceeding buffer size: " + rec.id);
        }
        for (unsigned char q : rec.qual) if (q < 33 || q > 126)
            throw std::runtime_error("FASTQ quality must use Phred+33 encoding: " + rec.id);
    }
    return found;
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
    const size_t CHUNK = std::numeric_limits<int>::max();

    while (total > 0) {
        unsigned int nbytes = static_cast<unsigned int>(std::min(total, CHUNK));
        if (is_gz_out) {
            if (gzwrite(out_gz, data, nbytes) != static_cast<int>(nbytes))
                throw std::runtime_error("Failed to write gzip output; check free disk space.");
        } else if (fwrite(data, 1, nbytes, out_plain) != nbytes) {
            throw std::runtime_error("Failed to write output; check free disk space.");
        }
        data += nbytes;
        total -= nbytes;
    }
}

void closeFile(FILE*& file) {
    if (!file) return;
    FILE* handle = file;
    file = nullptr;
    if (std::fclose(handle) != 0) throw std::runtime_error("Failed to close file; output may be incomplete.");
}

void closeGzip(gzFile& file) {
    if (!file) return;
    gzFile handle = file;
    file = nullptr;
    if (gzclose(handle) != Z_OK) throw std::runtime_error("Failed to close gzip file; output may be incomplete.");
}

void check_output_path(const std::string& output, const std::vector<std::string>& inputs) {
    if (output.empty()) throw std::invalid_argument("Output path cannot be empty.");
    const auto path = std::filesystem::weakly_canonical(expandTilde(output));
    for (const auto& input : inputs) {
        if (path == std::filesystem::weakly_canonical(expandTilde(input)))
            throw std::invalid_argument("Output must not overwrite an input file.");
    }
}

void check_processing_parameters(int threads, double batch_size, double compression) {
    if (threads < 1) throw std::invalid_argument("n_threads must be positive.");
    if (!std::isfinite(batch_size) || batch_size < 1 || batch_size > std::numeric_limits<int>::max() || batch_size != std::floor(batch_size))
        throw std::invalid_argument("batch_size must be a positive integer no larger than INT_MAX.");
    if (!std::isfinite(compression) || compression < 0 || compression > 9 || compression != std::floor(compression))
        throw std::invalid_argument("compress_level must be an integer from 0 to 9.");
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


size_t find_overlap_with_mismatches(
    const std::string& s1,
    const std::string& s2,
    bool& is_offset,
    size_t min_overlap,
    size_t max_mismatches,
    size_t& mismatch_count,
    size_t mismatch_positions[1000]
) {
    size_t max_possible = std::min(s1.size(), s2.size());
    mismatch_count = 0;
    is_offset = false;
    if (max_possible < min_overlap) {
        return 0;
    }

    auto count_mismatches = [&](const std::string& left,
                                size_t left_start,
                                const std::string& right,
                                size_t right_start,
                                size_t overlap,
                                size_t& mismatches) -> bool {
        mismatches = 0;
        for (size_t i = 0; i < overlap; ++i) {
            if (left[left_start + i] != right[right_start + i]) {
                if (mismatches < 1000) {
                    mismatch_positions[mismatches] = i;
                }
                ++mismatches;
                if (mismatches > max_mismatches) {
                    return false;
                }
            }
        }
        return true;
    };

    for (size_t L = max_possible; L >= min_overlap; --L) {
        size_t mismatches = 0;

        // R1 suffix overlaps R2 prefix.
        if (count_mismatches(s1, s1.size() - L, s2, 0, L, mismatches)) {
            mismatch_count = mismatches;
            is_offset = false;
            return L;
        }

        // R2 suffix overlaps R1 prefix.
        if (count_mismatches(s2, s2.size() - L, s1, 0, L, mismatches)) {
            mismatch_count = mismatches;
            is_offset = true;
            return L;
        }

        if (L == min_overlap) {
            break;
        }
    }

    mismatch_count = 0;
    is_offset = false;
    return 0;
}

static const std::array<unsigned char, 256> BASE_PROBABILITY_INDEX = [] {
    std::array<unsigned char, 256> index{};
    index.fill(4);
    index['A'] = 0; index['C'] = 1; index['G'] = 2; index['T'] = 3;
    return index;
}();

inline double cal_pr_diff(
    char x, char y,
    double correct_x, double correct_y,
    double error_x, double error_y,
    const base_probability& pb
) {
    const double frequencies[4] = {pb.A, pb.C, pb.G, pb.T};
    const unsigned char ix = BASE_PROBABILITY_INDEX[static_cast<unsigned char>(x)];
    const unsigned char iy = BASE_PROBABILITY_INDEX[static_cast<unsigned char>(y)];
    const double px = frequencies[ix];
    const double py = frequencies[iy];

    // Σ_{b≠Xi}
    double sum_not_x = 1.0 - px;

    // Σ_{b≠Yi}
    double sum_not_y = 1.0 - py;

    // Σ_{b≠Xi,Yi}
    double sum_not_xy = 1.0 - px - py;

    // Σ P_b^2 over b≠Xi,Yi
    const double all_squared = pb.A * pb.A + pb.C * pb.C +
                               pb.G * pb.G + pb.T * pb.T;
    const double sum_sq_not_xy = all_squared - px * px - py * py;

    double term1 = sum_not_x > 0.0 ? correct_y * error_x * (py / sum_not_x) : 0.0;
    double term2 = sum_not_y > 0.0 ? correct_x * error_y * (px / sum_not_y) : 0.0;
    double term3 = sum_not_xy > 0.0 ? error_x * error_y *
        (sum_sq_not_xy / (sum_not_xy * sum_not_xy)) : 0.0;

    return term1 + term2 + term3;
}


// ------------------------------
// ------------------------------
double oes_test(double q, int omega, double oes, double alpha, double beta, int max_length) {
	if (omega > max_length) return 0.0;
	if (q <= 0.0) return 1.0;
	// PEAR limits the lookup-table background probability to 0.49. Apply the
	// same bound before evaluating the equivalent tail probability here.
	q = std::min(q, 0.49);
	if (alpha == beta) return 1.0;

	const double one_minus_q = 1.0 - q;
	const double ratio = q / one_minus_q;
	long double p = 1.0;
    //
    long double one_minus_q_power =
        std::pow(one_minus_q, static_cast<double>(omega));

	for (int i = omega; i <= max_length; ++i) {
		int l_c = static_cast<int>(std::ceil((oes - beta * i) / (alpha - beta))) - 1;

		if (l_c < 0) {
			return 1.0;
		}
		if (l_c >= i) {
			continue;
		}

		// term = C(i,0) * q^0 * (1-q)^i
        long double term = one_minus_q_power;  
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
        one_minus_q_power *= one_minus_q; //next i (1-q)^(i+1)
	}
    double return_p;
    return_p = (double)(1.0 - p * p) < 0 ? 0: (double)(1.0 - p * p);
    return return_p;
}

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
) {
    const int r1_size = static_cast<int>(r1.seq.size());
    const int r2_size = static_cast<int>(r2.seq.size());

    std::string r2_reverse;
    std::string r2_q_reverse = r2.qual;
    fast_reverse_complement(r2.seq, r2_reverse);
    std::reverse(r2_q_reverse.begin(), r2_q_reverse.end());

    const std::string_view r1_seq(r1.seq), r2_seq(r2_reverse);
    const std::string_view r1_qual(r1.qual), r2_qual(r2_q_reverse);

    // Lookup tables avoid switches and repeated Phred conversions in the
    // candidate-scoring loop.
    const std::array<double, 5> base_probability_squared = {
        base_pb.A_sqr, base_pb.C_sqr, base_pb.G_sqr, base_pb.T_sqr, 0.0625
    };

    const double q = base_pb.A * base_pb.A +
    base_pb.C * base_pb.C +
    base_pb.G * base_pb.G +
    base_pb.T * base_pb.T;

    bool found = false;
    double best_score = std::numeric_limits<double>::lowest();
    double best_oes_score = 0.0;
    size_t best_mismatch = 0, best_length = 0;
    int best_offset = 0;

    const int min_offset = static_cast<int>(min_overlap) - r2_size;
    const int max_offset = r1_size - static_cast<int>(min_overlap);

    // Phase 1: score candidates. No consensus sequence is allocated here.
    for (int offset = min_offset; offset <= max_offset; ++offset) {
        const int overlap_begin = std::max(0, offset);
        const int overlap_end = std::min(r1_size, offset + r2_size);
        const int overlap_size = overlap_end - overlap_begin;
        if (overlap_size < static_cast<int>(min_overlap)) continue;

        const size_t length = static_cast<size_t>(overlap_size);
        const size_t r1_start = static_cast<size_t>(overlap_begin);
        const size_t r2_start = static_cast<size_t>(overlap_begin - offset);
        size_t candidate_mismatch = 0;
        double candidate_score = 0.0;
        double candidate_oes = 0.0;

        for (size_t i = 0; i < length; ++i) {
            const unsigned char b1 = static_cast<unsigned char>(r1_seq[r1_start + i]);
            const unsigned char b2 = static_cast<unsigned char>(r2_seq[r2_start + i]);
            const unsigned char q1 = static_cast<unsigned char>(r1_qual[r1_start + i]);
            const unsigned char q2 = static_cast<unsigned char>(r2_qual[r2_start + i]);
            const double correct1 = QUALITY_WEIGHT_TABLE[q1 - 33];
            const double correct2 = QUALITY_WEIGHT_TABLE[q2 - 33];
            const double error1 = 1.0 - correct1;
            const double error2 = 1.0 - correct2;

            const bool same = b1 == b2;
            candidate_mismatch += !same;
            double evidence;
            const unsigned char base_index1 = BASE_PROBABILITY_INDEX[b1];
            const unsigned char base_index2 = BASE_PROBABILITY_INDEX[b2];
            if (base_index1 >= 4 || base_index2 >= 4) {
                // PEAR treats every position containing N as background
                // evidence, regardless of whether both calls are N.
                evidence = q;
            } else if (same) {
                evidence = correct1 * correct2 + error1 * error2 *
                           base_probability_squared[base_index1];
            } else {
                evidence = cal_pr_diff(static_cast<char>(b1), static_cast<char>(b2),
                                       correct1, correct2, error1, error2, base_pb);
            }
            // PEAR uses AS to choose the overlap: an observed match receives
            // its true-match contribution, whereas an observed mismatch
            // receives its true-mismatch contribution.
            candidate_score += same
                ? evidence * alpha
                : (1.0 - evidence) * beta;

            // OES includes both possible true states and is used only for the
            // statistical test of the selected overlap.
            candidate_oes += evidence * alpha + (1.0 - evidence) * beta;
        }

        if (!found || candidate_score > best_score ||
            (candidate_score == best_score && length > best_length) ||
            (candidate_score == best_score && length == best_length &&
             candidate_mismatch < best_mismatch)) {
            found = true;
            best_score = candidate_score;
            best_oes_score = candidate_oes;
            best_mismatch = candidate_mismatch;
            best_length = length;
            best_offset = offset;
        }
    }

    mismatch_count = found ? best_mismatch : 0;
    overlap_length = found ? best_length : 0;
    overlap_score = found ? best_score : 0.0;
    const int maximum_overlap = std::max(r1_size, r2_size);
    p_value = found ? oes_test(q, static_cast<int>(min_overlap), best_oes_score,
                               alpha, beta, maximum_overlap) : 1.0;
    p_value = std::max(0.0, std::min(1.0, p_value));

    // Failed reads leave merged_rec empty and perform no consensus work.
    if (!found || p_value >= 0.01 || best_mismatch > max_mismatches) {
        merged_rec = FastqRecord{};
        return false;
    }

    // Phase 2: construct the consensus once, for the accepted best overlap.
    const size_t r1_start = static_cast<size_t>(std::max(0, best_offset));
    const size_t r2_start = static_cast<size_t>(std::max(0, -best_offset));
    merged_rec.seq.clear();
    merged_rec.qual.clear();
    merged_rec.seq.reserve(r1_seq.size() + r2_seq.size() - best_length);
    merged_rec.qual.reserve(r1_qual.size() + r2_qual.size() - best_length);

    if (best_offset >= 0) {
        merged_rec.seq.append(r1_seq.data(), r1_start);
        merged_rec.qual.append(r1_qual.data(), r1_start);
    } else {
        merged_rec.seq.append(r2_seq.data(), r2_start);
        merged_rec.qual.append(r2_qual.data(), r2_start);
    }

    for (size_t i = 0; i < best_length; ++i) {
        const char b1 = r1_seq[r1_start + i], b2 = r2_seq[r2_start + i];
        const char q1 = r1_qual[r1_start + i], q2 = r2_qual[r2_start + i];
        char consensus_base, consensus_quality;
        const bool called1 = BASE_PROBABILITY_INDEX[static_cast<unsigned char>(b1)] < 4;
        const bool called2 = BASE_PROBABILITY_INDEX[static_cast<unsigned char>(b2)] < 4;
        if (!called1 && !called2) {
            consensus_base = 'N';
            consensus_quality = std::min(q1, q2);
        } else if (!called1) {
            consensus_base = b2;
            consensus_quality = q2;
        } else if (!called2) {
            consensus_base = b1;
            consensus_quality = q1;
        } else if (b1 == b2) {
            consensus_base = b1;
            const int combined_q = (q1 - 33) + (q2 - 33);
            consensus_quality = static_cast<char>(33 + std::min(combined_q, 93));
        } else {
            // PEAR keeps the second call when qualities are tied.
            const bool choose_first = q1 > q2;
            consensus_base = choose_first ? b1 : b2;
            consensus_quality = choose_first ? q1 : q2;
        }
        merged_rec.seq.push_back(consensus_base);
        merged_rec.qual.push_back(consensus_quality);
    }

    const size_t r1_tail = r1_start + best_length;
    const size_t r2_tail = r2_start + best_length;
    if (r1_size >= best_offset + r2_size) {
        merged_rec.seq.append(r1_seq.data() + r1_tail, r1_seq.size() - r1_tail);
        merged_rec.qual.append(r1_qual.data() + r1_tail, r1_qual.size() - r1_tail);
    } else {
        merged_rec.seq.append(r2_seq.data() + r2_tail, r2_seq.size() - r2_tail);
        merged_rec.qual.append(r2_qual.data() + r2_tail, r2_qual.size() - r2_tail);
    }
    return true;
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


static void add_base_counts(
    const FastqRecord& rec,
    size_t freq[256],
    double& total_length,
    int& max_length,
    bool reverse_complement
) {
    int temp_length = 0;
    const auto& comp_table = get_comp_table();
    for (unsigned char c : rec.seq) {
        ++temp_length;
        unsigned char counted = reverse_complement && comp_table[c] ? comp_table[c] : c;
        ++freq[counted];
        const unsigned char upper = static_cast<unsigned char>(std::toupper(counted));
        if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
            ++total_length;
        }
    }
    max_length = temp_length > max_length ? temp_length : max_length;
}

base_probability calculate_base_probabilities(const std::vector <FastqRecord>& rec_v) {
    base_probability base_pb{0.0, 0.0, 0.0, 0.0, 0.0};
    double total_length = 0;
    size_t freq[256] = { 0 };
    int max_length = 0;
    for (const FastqRecord &rec:  rec_v) {
        add_base_counts(rec, freq, total_length, max_length, false);
    }
    if (total_length <= 0.0) {
        base_pb.A = base_pb.T = base_pb.C = base_pb.G = 0.25;
        base_pb.N = 0.0;
        base_pb.A_sqr = base_pb.T_sqr = base_pb.C_sqr = base_pb.G_sqr = 1.0 / 9.0;
        base_pb.max_length = 0;
        return base_pb;
    }
    
    base_pb.A = static_cast<double>(freq['A'] + freq['a']) / total_length;
    base_pb.T = static_cast<double>(freq['T'] + freq['t']) / total_length;
    base_pb.C = static_cast<double>(freq['C'] + freq['c']) / total_length;
    base_pb.G = static_cast<double>(freq['G'] + freq['g']) / total_length;
    const double n_count = static_cast<double>(freq['N'] + freq['n']);
    base_pb.N = n_count / (total_length + n_count);
	
    // ------
    double denom_A = (base_pb.T + base_pb.C + base_pb.G) * (base_pb.T + base_pb.C + base_pb.G);
    double denom_T = (base_pb.A + base_pb.C + base_pb.G) * (base_pb.A + base_pb.C + base_pb.G);
    double denom_C = (base_pb.A + base_pb.T + base_pb.G) * (base_pb.A + base_pb.T + base_pb.G);
    double denom_G = (base_pb.A + base_pb.C + base_pb.T) * (base_pb.A + base_pb.C + base_pb.T);
    base_pb.A_sqr = denom_A > 0.0 ?
        (base_pb.T * base_pb.T + base_pb.C * base_pb.C + base_pb.G * base_pb.G) / denom_A : 0.0;
    base_pb.T_sqr = denom_T > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.C * base_pb.C + base_pb.G * base_pb.G) / denom_T : 0.0;
    base_pb.C_sqr = denom_C > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.T * base_pb.T + base_pb.G * base_pb.G) / denom_C : 0.0;
    base_pb.G_sqr = denom_G > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.C * base_pb.C + base_pb.T * base_pb.T) / denom_G : 0.0;
    //
    base_pb.max_length = max_length;
	return base_pb;
}

base_probability calculate_base_probabilities(
    const std::vector<FastqRecord>& rec_v1,
    const std::vector<FastqRecord>& rec_v2,
    bool reverse_complement_r2
) {
    std::vector<FastqRecord> rec_v;
    rec_v.reserve(rec_v1.size() + rec_v2.size());
    rec_v.insert(rec_v.end(), rec_v1.begin(), rec_v1.end());

    if (!reverse_complement_r2) {
        rec_v.insert(rec_v.end(), rec_v2.begin(), rec_v2.end());
        return calculate_base_probabilities(rec_v);
    }

    base_probability base_pb{0.0, 0.0, 0.0, 0.0, 0.0};
    double total_length = 0;
    size_t freq[256] = { 0 };
    int max_length = 0;

    for (const FastqRecord& rec : rec_v1) {
        add_base_counts(rec, freq, total_length, max_length, false);
    }
    // R2 is compared after reverse-complementing. Reversal does not change
    // counts, but complementation swaps A/T and C/G, so count R2 in merged
    // orientation before estimating the OES background match probability.
    for (const FastqRecord& rec : rec_v2) {
        add_base_counts(rec, freq, total_length, max_length, true);
    }

    if (total_length <= 0.0) {
        base_pb.A = base_pb.T = base_pb.C = base_pb.G = 0.25;
        base_pb.N = 0.0;
        base_pb.A_sqr = base_pb.T_sqr = base_pb.C_sqr = base_pb.G_sqr = 1.0 / 9.0;
        base_pb.max_length = 0;
        return base_pb;
    }

    base_pb.A = static_cast<double>(freq['A'] + freq['a']) / total_length;
    base_pb.T = static_cast<double>(freq['T'] + freq['t']) / total_length;
    base_pb.C = static_cast<double>(freq['C'] + freq['c']) / total_length;
    base_pb.G = static_cast<double>(freq['G'] + freq['g']) / total_length;
    const double n_count = static_cast<double>(freq['N'] + freq['n']);
    base_pb.N = n_count / (total_length + n_count);

    double denom_A = (base_pb.T + base_pb.C + base_pb.G) * (base_pb.T + base_pb.C + base_pb.G);
    double denom_T = (base_pb.A + base_pb.C + base_pb.G) * (base_pb.A + base_pb.C + base_pb.G);
    double denom_C = (base_pb.A + base_pb.T + base_pb.G) * (base_pb.A + base_pb.T + base_pb.G);
    double denom_G = (base_pb.A + base_pb.C + base_pb.T) * (base_pb.A + base_pb.C + base_pb.T);
    base_pb.A_sqr = denom_A > 0.0 ?
        (base_pb.T * base_pb.T + base_pb.C * base_pb.C + base_pb.G * base_pb.G) / denom_A : 0.0;
    base_pb.T_sqr = denom_T > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.C * base_pb.C + base_pb.G * base_pb.G) / denom_T : 0.0;
    base_pb.C_sqr = denom_C > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.T * base_pb.T + base_pb.G * base_pb.G) / denom_C : 0.0;
    base_pb.G_sqr = denom_G > 0.0 ?
        (base_pb.A * base_pb.A + base_pb.C * base_pb.C + base_pb.T * base_pb.T) / denom_G : 0.0;
    base_pb.max_length = max_length;
    return base_pb;
}
