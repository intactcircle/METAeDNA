#include "utils.hpp"
#include <exception>

struct FastqProfile {
    std::vector<std::array<size_t, 5>> bases;
    std::vector<std::array<size_t, 94>> qualities;
    size_t reads = 0, total_bases = 0, q20 = 0, q30 = 0;
    size_t overlap_bases = 0, mismatches = 0, annotated_reads = 0;
};

FastqProfile profile_fastq_core(const std::string& path, size_t max_records,
                               int n_threads = 4, int batch_size = 1000000,
                               bool debug = false) {
    check_processing_parameters(n_threads, batch_size, 0);
    const std::string input = expandTilde(path);
    const bool compressed = isGzFile(input);
    FILE* plain = compressed ? nullptr : fopen(input.c_str(), "rb");
    gzFile gzip = compressed ? gzopen(input.c_str(), "rb") : nullptr;
    FileGuard guard(plain, gzip);
    if (compressed ? !gzip : !plain) throw std::runtime_error("Cannot open FASTQ: " + input);
    if (debug) COUT << "[fastq_profile] Input opened." << std::endl;
    FastqProfile result;
    int workers = 1;
#ifdef _OPENMP
    workers = n_threads;
#endif
    std::vector<FastqProfile> partial(workers);
    std::vector<FastqRecord> records;
    size_t reads = 0;
    for (;;) {
        records.clear();
        FastqRecord record;
        while (records.size() < static_cast<size_t>(batch_size) &&
               (!max_records || reads + records.size() < max_records) &&
               readRec(compressed, plain, gzip, record, true)) {
            records.push_back(std::move(record));
#ifdef USE_RCPP_OUTPUT
            if (records.size() % 10000 == 0) Rcpp::checkUserInterrupt();
#endif
        }
        if (records.empty()) break;
        if (debug) COUT << "[fastq_profile] Read " << records.size() << " records; starting statistics." << std::endl;
        reads += records.size();
        std::vector<std::exception_ptr> errors(workers);
        // Thread-local histograms avoid locks in the per-base loop.
#ifdef _OPENMP
#pragma omp parallel num_threads(workers)
#endif
        {
            int worker = 0;
#ifdef _OPENMP
            worker = omp_get_thread_num();
#endif
            FastqProfile& local = partial[worker];
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (int index = 0; index < static_cast<int>(records.size()); ++index) {
                try {
                    const FastqRecord& record = records[index];
                    ++local.reads;
                    local.total_bases += record.seq.size();
                    if (record.seq.size() > local.bases.size()) {
                        local.bases.resize(record.seq.size());
                        local.qualities.resize(record.seq.size());
                    }
                    for (size_t i = 0; i < record.seq.size(); ++i) {
                        char base = std::toupper(static_cast<unsigned char>(record.seq[i]));
                        size_t index = base == 'A' ? 0 : base == 'C' ? 1 : base == 'G' ? 2 : base == 'T' ? 3 : 4;
                        ++local.bases[i][index];
                        unsigned int quality = static_cast<unsigned char>(record.qual[i]) - 33;
                        ++local.qualities[i][quality];
                        local.q20 += quality >= 20; local.q30 += quality >= 30;
                    }
                    if (record.id.find("\"pairend\"") != std::string::npos) {
                        size_t a = 0, b = 0, overlap_position = 0, mismatch_position = 0;
                        if (header_count(record.id, "\"overlap\"", a, overlap_position) &&
                            header_count(record.id, "\"mismatch\"", b, mismatch_position)) {
                            local.overlap_bases += a;
                            local.mismatches += b;
                            ++local.annotated_reads;
                        }
                    }
                } catch (...) {
                    errors[worker] = std::current_exception();
                }
            }
        }
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
        if (debug) COUT << "[fastq_profile] Batch statistics complete." << std::endl;
#ifdef USE_RCPP_OUTPUT
        Rcpp::checkUserInterrupt();
#endif
    }
    for (const FastqProfile& local : partial) {
        result.reads += local.reads;
        result.total_bases += local.total_bases;
        result.q20 += local.q20;
        result.q30 += local.q30;
        result.overlap_bases += local.overlap_bases;
        result.mismatches += local.mismatches;
        result.annotated_reads += local.annotated_reads;
        result.bases.resize(std::max(result.bases.size(), local.bases.size()));
        result.qualities.resize(result.bases.size());
        for (size_t i = 0; i < local.bases.size(); ++i) {
            for (size_t b = 0; b < 5; ++b) result.bases[i][b] += local.bases[i][b];
            for (size_t q = 0; q < 94; ++q) result.qualities[i][q] += local.qualities[i][q];
        }
    }
    return result;
}

#ifdef USE_RCPP_OUTPUT
//' Validate paired FASTQ files
//'
//' Read R1 and R2 together and check that every FASTQ record has a corresponding
//' mate before read merging or paired demultiplexing.
//'
//' @param forward Character string giving the path to the R1 FASTQ file. Plain
//'   text and gzip-compressed files are supported.
//' @param reverse Character string giving the path to the R2 FASTQ file. Plain
//'   text and gzip-compressed files are supported.
//' @param n_threads Number of threads. Defaults to 4, as in `pairend()`.
//' @param batch_size Number of read pairs per batch. Defaults to 1000000,
//'   as in `pairend()`.
//'
//' @details
//' Both files are read from beginning to end. For every pair, the function
//' checks the FASTQ structure, sequence and quality lengths, Phred+33 quality
//' range, and read identifier. Only the first field of each identifier is
//' compared; text following a space or tab is ignored. Conventional `/1` and
//' `/2` mate suffixes are also ignored. For example, `@read42/1` and
//' `@read42/2` are treated as the same read identifier.
//'
//' An error is raised when the files are identical, cannot be opened, contain
//' invalid FASTQ records, contain different numbers of records, or have mates
//' in different orders. Empty input files are not accepted. The function only
//' validates the files and does not modify or create sequence data.
//' Files are read in batches using the shared FASTQ reader. With multiple
//' threads, the two files are read concurrently and identifiers are checked
//' in parallel. Without OpenMP, processing is sequential.
//'
//' @return A numeric value giving the number of validated read pairs.
//'
//' @examples
//' forward <- system.file("extdata", "trnL_40k_1.fq.gz", package = "METAeDNA")
//' reverse <- system.file("extdata", "trnL_40k_2.fq.gz", package = "METAeDNA")
//' validate_paired_fastq(forward, reverse)
//'
//' @export
// [[Rcpp::export]]
double validate_paired_fastq(std::string forward, std::string reverse,
                             int n_threads = 4, int batch_size = 1000000) {
    check_processing_parameters(n_threads, batch_size, 0);
    forward = expandTilde(forward); reverse = expandTilde(reverse);
    if (std::filesystem::weakly_canonical(forward) == std::filesystem::weakly_canonical(reverse))
        throw std::invalid_argument("R1 and R2 must be different files.");
    const bool gz1 = isGzFile(forward), gz2 = isGzFile(reverse);
    FILE *p1 = gz1 ? nullptr : fopen(forward.c_str(), "rb"), *p2 = gz2 ? nullptr : fopen(reverse.c_str(), "rb");
    gzFile z1 = gz1 ? gzopen(forward.c_str(), "rb") : nullptr, z2 = gz2 ? gzopen(reverse.c_str(), "rb") : nullptr;
    FileGuard guard1(p1, z1), guard2(p2, z2);
    if ((gz1 ? !z1 : !p1) || (gz2 ? !z2 : !p2)) throw std::runtime_error("Cannot open paired FASTQ inputs.");
    size_t count = 0;
    std::vector<FastqRecord> records[2];
    for (;;) {
        std::exception_ptr errors[2];
        // Each worker owns one file handle; shared readers validate FASTQ records.
#ifdef _OPENMP
#pragma omp parallel for num_threads(std::min(n_threads, 2))
#endif
        for (int mate = 0; mate < 2; ++mate) {
            records[mate].clear();
            try {
                FastqRecord record;
                while (records[mate].size() < static_cast<size_t>(batch_size) &&
                       readRec(mate == 0 ? gz1 : gz2, mate == 0 ? p1 : p2,
                               mate == 0 ? z1 : z2, record, true)) {
                    records[mate].push_back(std::move(record));
                }
            } catch (...) {
                errors[mate] = std::current_exception();
            }
        }
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
        if (records[0].size() != records[1].size())
            throw std::runtime_error("R1 and R2 FASTQ files have different record counts.");
        const int size = static_cast<int>(records[0].size());
        if (!size) break;
        int mismatched = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) reduction(|:mismatched)
#endif
        for (int i = 0; i < size; ++i) {
            mismatched |= paired_read_id(records[0][i].id) != paired_read_id(records[1][i].id);
        }
        if (mismatched) throw std::runtime_error("R1 and R2 identifiers are out of sync.");
        count += size;
        Rcpp::checkUserInterrupt();
    }
    if (!count) throw std::runtime_error("Paired FASTQ inputs contain no reads.");
    return double(count);
}

//' Summarize FASTQ quality and base composition
//'
//' Inspect reads from the beginning of a FASTQ file and summarize sequencing
//' quality and nucleotide composition at each read position. The result can be
//' used to examine quality loss toward the read end, positional nucleotide
//' bias, and the number of reads that reach each position.
//'
//' @param path Character string giving the path to a plain or gzip-compressed
//'   FASTQ file. Quality scores must use Phred+33 encoding.
//' @param max_records Non-negative integer giving the maximum number of records
//'   read from the beginning of the file. The default is 100000. Use 0 to scan
//'   the complete file. This is a prefix preview rather than a random sample.
//' @param n_threads Number of threads used for statistics. Defaults to 4.
//' @param batch_size Number of reads per batch. Defaults to 1000000.
//' @param debug Print processing stages to help locate slow operations.
//'
//' @details
//' The shared FASTQ reader reads one batch at a time; quality and nucleotide
//' counts are computed in parallel. Only the current batch is retained.
//' Without OpenMP, statistics are computed sequentially.
//' A, C, G, and T are reported
//' separately. N and any other ambiguity codes are combined as `N_other`.
//'
//' If the identifiers contain the `pairend` metadata written by `pairend()`,
//' the function also summarizes overlap lengths and mismatch counts. These
//' fields describe only the merged reads present in the inspected file. The
//' resulting `mismatch_rate` is not an error rate measured against a reference
//' sequence.
//'
//' @return A list with two data frames:
//' \describe{
//'   \item{`positions`}{One row per read position, containing:
//'     \itemize{
//'       \item `position`: base position counted from the beginning of the
//'         stored read, starting at 1. Reads are not realigned or reoriented.
//'       \item `reads`: number of inspected reads containing a base at this
//'         position. Shorter reads do not contribute to later positions.
//'       \item `mean_quality`: arithmetic mean of the Phred scores at this
//'         position, after subtracting the Phred+33 encoding offset. This is
//'         a quality score, not a mean error probability.
//'       \item `q10`: 10th percentile of the Phred scores at this position.
//'       \item `q25` and `q75`: 25th and 75th percentiles, used as the lower
//'         and upper edges of quality boxes.
//'       \item `median`: 50th percentile of the Phred scores at this position.
//'       \item `q90`: 90th percentile of the Phred scores at this position.
//'         Each percentile is the smallest observed score whose cumulative
//'         count reaches the requested fraction; no interpolation is used.
//'       \item `A`: fraction of bases called A at this position.
//'       \item `C`: fraction of bases called C at this position.
//'       \item `G`: fraction of bases called G at this position.
//'       \item `T`: fraction of bases called T at this position.
//'       \item `N_other`: fraction of N bases and all other ambiguity codes
//'         at this position. Base counting is case-insensitive. All five
//'         fractions use `reads` as the denominator, range from 0 to 1,
//'         and sum to 1 at each position; they are not percentages.
//'     }}
//'   \item{`summary`}{One row containing:
//'     \itemize{
//'       \item `filename`: input filename without its directory.
//'       \item `file_type`: conventional nucleotide base calls in FASTQ format.
//'       \item `encoding`: required Sanger / Illumina 1.8+ (Phred+33) encoding;
//'         this is the assumed encoding, not automatic encoding detection.
//'       \item `total_megabases`: inspected bases divided by 1000000.
//'       \item `min_length`, `max_length`, `mean_length`, and `median_length`:
//'         shortest, longest, mean, and median inspected read lengths in bases.
//'       \item `gc_percent`: 100 times the number of G and C calls divided by
//'         all inspected bases, including ambiguous bases.
//'       \item `reads` and `bases`: inspected read and base counts;
//'       \item `q20_fraction` and `q30_fraction`: fractions of bases with
//'         Phred scores of at least 20 and 30;
//'       \item `annotated_reads`: reads containing valid `pairend` metadata;
//'       \item `overlap_bases` and `overlap_mismatches`: corresponding totals;
//'       \item `mismatch_rate`: `overlap_mismatches / overlap_bases`, or `NA`
//'         when overlap metadata is unavailable.
//'     }}
//' }
//' All summary counts and length statistics describe only the inspected
//' records. Use `max_records = 0` to summarize the whole file. Empty files
//' have zero counts and missing length and GC statistics.
//'
//' @examples
//' fastq <- system.file("extdata", "trnL_40k_1.fq.gz", package = "METAeDNA")
//' result <- fastq_profile(fastq, max_records = 1000)
//' head(result$positions)
//' result$summary
//' # Boxes: 25th--75th percentiles; whiskers: 10th--90th percentiles.
//' ggplot2::ggplot(result$positions, ggplot2::aes(position, mean_quality)) +
//'   ggplot2::geom_boxplot(ggplot2::aes(group = position, ymin = q10,
//'     lower = q25, middle = median, upper = q75, ymax = q90),
//'     stat = "identity", width = 0.7, fill = "#A6CEE3") +
//'   ggplot2::geom_line(colour = "#3B4992") +
//'   ggplot2::labs(x = "Base position", y = "Phred quality") +
//'   ggplot2::theme_classic(base_size = 10)
//'
//' \donttest{
//' merged <- tempfile(fileext = ".fastq.gz")
//' reverse <- system.file("extdata", "trnL_40k_2.fq.gz", package = "METAeDNA")
//' pairend(fastq, reverse, merged, n_threads = 1)
//' fastq_profile(merged, max_records = 1000)$summary
//' }
//'
//' @export
// [[Rcpp::export]]
List fastq_profile(std::string path, int max_records = 100000,
                   int n_threads = 4, int batch_size = 1000000, bool debug = false) {
    if (debug) COUT << "[fastq_profile] Entered C++." << std::endl;
    if (max_records < 0) throw std::invalid_argument("max_records must be non-negative.");
    FastqProfile profile = profile_fastq_core(path, max_records, n_threads, batch_size, debug);
    if (debug) COUT << "[fastq_profile] Building R result." << std::endl;
    size_t n = profile.bases.size();
    IntegerVector position(n);
    NumericVector depth(n), mean_quality(n), q10(n), q25(n), median(n), q75(n), q90(n);
    NumericMatrix frequencies(n, 5);
    size_t gc = 0;
    for (size_t i = 0; i < n; ++i) {
        position[i] = i + 1;
        double total = 0, weighted = 0;
        for (size_t q = 0; q < 94; ++q) { total += profile.qualities[i][q]; weighted += q * profile.qualities[i][q]; }
        depth[i] = total;
        gc += profile.bases[i][1] + profile.bases[i][2];
        mean_quality[i] = weighted / total;
        size_t cumulative = 0;
        bool low_found = false, lower_found = false, median_found = false,
             upper_found = false, high_found = false;
        for (size_t q = 0; q < 94; ++q) {
            cumulative += profile.qualities[i][q];
            if (!low_found && cumulative >= total * 0.1) { q10[i] = q; low_found = true; }
            if (!lower_found && cumulative >= total * 0.25) { q25[i] = q; lower_found = true; }
            if (!median_found && cumulative >= total * 0.5) { median[i] = q; median_found = true; }
            if (!upper_found && cumulative >= total * 0.75) { q75[i] = q; upper_found = true; }
            if (!high_found && cumulative >= total * 0.9) { q90[i] = q; high_found = true; }
        }
        for (size_t b = 0; b < 5; ++b) frequencies(i, b) = profile.bases[i][b] / total;
    }
    // Position depth also gives read-length counts, without storing lengths.
    double min_length = NA_REAL, median_length = NA_REAL;
    size_t cumulative_lengths = 0;
    double middle_sum = 0;
    const size_t middle1 = (profile.reads + 1) / 2, middle2 = profile.reads / 2 + 1;
    for (size_t i = 0; i < n; ++i) {
        const size_t count = static_cast<size_t>(depth[i] - (i + 1 < n ? depth[i + 1] : 0));
        if (count && Rcpp::NumericVector::is_na(min_length)) min_length = i + 1;
        if (cumulative_lengths < middle1 && cumulative_lengths + count >= middle1) middle_sum += i + 1;
        if (cumulative_lengths < middle2 && cumulative_lengths + count >= middle2) middle_sum += i + 1;
        cumulative_lengths += count;
    }
    if (profile.reads) median_length = middle_sum / 2;
    return List::create(
        Named("positions") = DataFrame::create(Named("position") = position, Named("reads") = depth,
            Named("mean_quality") = mean_quality, Named("q10") = q10, Named("q25") = q25,
            Named("median") = median, Named("q75") = q75, Named("q90") = q90,
            Named("A") = frequencies(_, 0), Named("C") = frequencies(_, 1), Named("G") = frequencies(_, 2),
            Named("T") = frequencies(_, 3), Named("N_other") = frequencies(_, 4)),
        Named("summary") = DataFrame::create(Named("reads") = double(profile.reads), Named("bases") = double(profile.total_bases),
            Named("filename") = std::filesystem::path(expandTilde(path)).filename().string(),
            Named("file_type") = "Conventional base calls",
            Named("encoding") = "Sanger / Illumina 1.8+ (Phred+33, assumed)",
            Named("total_megabases") = double(profile.total_bases) / 1000000,
            Named("min_length") = min_length, Named("max_length") = profile.reads ? double(n) : NA_REAL,
            Named("mean_length") = profile.reads ? double(profile.total_bases) / profile.reads : NA_REAL,
            Named("median_length") = median_length,
            Named("gc_percent") = profile.total_bases ? 100.0 * gc / profile.total_bases : NA_REAL,
            Named("q20_fraction") = profile.total_bases ? double(profile.q20) / profile.total_bases : NA_REAL,
            Named("q30_fraction") = profile.total_bases ? double(profile.q30) / profile.total_bases : NA_REAL,
            Named("annotated_reads") = double(profile.annotated_reads),
            Named("overlap_bases") = double(profile.overlap_bases), Named("overlap_mismatches") = double(profile.mismatches),
            Named("mismatch_rate") = profile.overlap_bases ? double(profile.mismatches) / profile.overlap_bases : NA_REAL));
}
#endif
