#define _CRT_SECURE_NO_WARNINGS
#include <cerrno>         // For errno
#include <climits>        // For INT_MAX
#include <cstdlib>        // For strtol, strtod
#include <functional>     // For std::function
#include <stdexcept>      // For runtime_error
#include <string>         // For std::string
#include <vector>         // For std::vector
#include <zlib.h>         // For gzFile, gzopen, gzgets, gzclose
#include "utils.hpp"

namespace {

bool read_gzip_line(gzFile input, std::string& line) {
    line.clear();
    char buffer[65536];

    while (true) {
        char* result = gzgets(input, buffer, sizeof(buffer));
        if (result == nullptr) {
            if (line.empty()) return false;
            break;
        }

        line.append(buffer);
        if (!line.empty() && line.back() == '\n') break;
    }

    if (!line.empty() && line.back() == '\n') line.pop_back();
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return true;
}

void split_accession_row(
    const std::string& line,
    std::string& accession,
    std::string& accession_version,
    int& taxid,
    double& gi,
    std::size_t line_number) {
    const std::size_t first = line.find('\t');
    const std::size_t second =
        first == std::string::npos ? std::string::npos : line.find('\t', first + 1);
    const std::size_t third =
        second == std::string::npos ? std::string::npos : line.find('\t', second + 1);

    if (first == std::string::npos || second == std::string::npos ||
        third == std::string::npos || line.find('\t', third + 1) != std::string::npos) {
        throw std::runtime_error(
            "Invalid accession2taxid row at line " + std::to_string(line_number) +
            ": expected four tab-separated fields.");
    }

    accession.assign(line, 0, first);
    accession_version.assign(line, first + 1, second - first - 1);
    const std::string taxid_text = line.substr(second + 1, third - second - 1);
    const std::string gi_text = line.substr(third + 1);

    if (accession.empty() || accession_version.empty() ||
        taxid_text.empty() || gi_text.empty()) {
        throw std::runtime_error(
            "Invalid accession2taxid row at line " + std::to_string(line_number) +
            ": fields must not be empty.");
    }

    char* end = nullptr;
    errno = 0;
    const long parsed_taxid = std::strtol(taxid_text.c_str(), &end, 10);
    if (errno != 0 || end == taxid_text.c_str() || *end != '\0' ||
        parsed_taxid < 0 || parsed_taxid > INT_MAX) {
        throw std::runtime_error(
            "Invalid taxid at line " + std::to_string(line_number) + ".");
    }

    end = nullptr;
    errno = 0;
    const double parsed_gi = std::strtod(gi_text.c_str(), &end);
    if (errno != 0 || end == gi_text.c_str() || *end != '\0' || parsed_gi < 0) {
        throw std::runtime_error(
            "Invalid gi at line " + std::to_string(line_number) + ".");
    }

    taxid = static_cast<int>(parsed_taxid);
    gi = parsed_gi;
}

using BatchCallback = std::function<void(
    const std::vector<std::string>& accessions,
    const std::vector<std::string>& accession_versions,
    const std::vector<int>& taxids,
    const std::vector<double>& gis)>;

double stream_accession2taxid_core(
    const std::string& input_file,
    int chunk_size,
    const BatchCallback& callback) {
    if (chunk_size < 1) {
        throw std::invalid_argument("chunk_size must be a positive integer.");
    }

    std::string real_path = expandTilde(input_file);
    gzFile input = gzopen(real_path.c_str(), "rb");
    if (input == nullptr) {
        throw std::runtime_error("Cannot open gzip file: " + input_file);
    }

    std::vector<std::string> accessions;
    std::vector<std::string> accession_versions;
    std::vector<int> taxids;
    std::vector<double> gis;
    accessions.reserve(chunk_size);
    accession_versions.reserve(chunk_size);
    taxids.reserve(chunk_size);
    gis.reserve(chunk_size);

    std::string line;
    std::size_t line_number = 0;
    double total_records = 0;

    try {
        if (!read_gzip_line(input, line)) {
            throw std::runtime_error("The accession2taxid file is empty.");
        }
        ++line_number;
        if (line != "accession\taccession.version\ttaxid\tgi") {
            throw std::runtime_error(
                "Invalid accession2taxid header in: " + input_file);
        }

        while (read_gzip_line(input, line)) {
            ++line_number;
            if (line.empty()) continue;

            std::string accession;
            std::string accession_version;
            int taxid = 0;
            double gi = 0;
            split_accession_row(
                line, accession, accession_version, taxid, gi, line_number);

            accessions.push_back(std::move(accession));
            accession_versions.push_back(std::move(accession_version));
            taxids.push_back(taxid);
            gis.push_back(gi);

            if (accessions.size() == static_cast<std::size_t>(chunk_size)) {
                callback(accessions, accession_versions, taxids, gis);
                total_records += static_cast<double>(accessions.size());
                accessions.clear();
                accession_versions.clear();
                taxids.clear();
                gis.clear();
            }
        }

        int zlib_error = Z_OK;
        const char* zlib_message = gzerror(input, &zlib_error);
        if (zlib_error != Z_OK && zlib_error != Z_STREAM_END) {
            throw std::runtime_error(
                "Failed while reading gzip file: " +
                std::string(zlib_message == nullptr ? "unknown zlib error" : zlib_message));
        }

        if (!accessions.empty()) {
            callback(accessions, accession_versions, taxids, gis);
            total_records += static_cast<double>(accessions.size());
        }
    } catch (...) {
        gzclose(input);
        throw;
    }

    if (gzclose(input) != Z_OK) {
        throw std::runtime_error("Failed to close gzip file: " + input_file);
    }
    return total_records;
}

}  // namespace

#ifdef USE_RCPP_OUTPUT
// ------------------------------
// Rcpp Exported Function: stream_accession2taxid_gz
// ------------------------------
// Stream an NCBI accession2taxid gzip file into an R callback.
//
// The parser retains at most chunk_size records before invoking callback.
// [[Rcpp::export]]
double stream_accession2taxid_gz(
    const std::string& input_file,
    int chunk_size,
    Function callback) {
    BatchCallback emit_batch = [&callback](
        const std::vector<std::string>& accessions,
        const std::vector<std::string>& accession_versions,
        const std::vector<int>& taxids,
        const std::vector<double>& gis) {
        callback(DataFrame::create(
            Named("accession") = wrap(accessions),
            Named("accession.version") = wrap(accession_versions),
            Named("taxid") = wrap(taxids),
            Named("gi") = wrap(gis)
        ));
        checkUserInterrupt();
    };

    return stream_accession2taxid_core(input_file, chunk_size, emit_batch);
}
#endif
