/**
 * mzgaf.hpp: Functions to parse output from minigraph --write-mz which outputs GAF-like records that contain minimizer hits
 *          mapping the query to the target.  Rely on gafkluge to read the normal GAF lines, and add here some functions
 *          to read the minimizer lines. 
 */


#pragma once
// copied from https://github.com/vgteam/libvgio/blob/master/include/vg/io/gafkluge.hpp
#include <cassert>
#include "gafkluge.hpp"

using namespace std;

namespace gafkluge {

/** from minigraph -S --write-mz
 *      segName segLen  nMinis  seqDiv  segStart segEnd qStart  qEnd    k-mer   segOffsets                      qOffsets
 *	>s79	124	15	0.0096	7	 118	193184	193299	19	8,9,3,8,9,11,9,4,2,7,2,9,3,8	8,9,3,8,9,15,9,4,2,7,2,9,3,8
 */
struct MzGafRecord {
    string query_name;  // parsed from parent gaf record
    int64_t query_length; // parsed from parent gaf record
    string target_name;
    bool is_reverse;
    int64_t target_length;
    int64_t num_minimizers;
    double seq_div;
    int64_t target_start;
    int64_t target_end;
    int64_t query_start;
    int64_t query_end;
    int64_t kmer_size;
    vector<int32_t> target_mz_offsets;
    vector<int32_t> query_mz_offsets;
};

inline double string_to_float(const std::string& s) {
    return s == missing_string ? (double)missing_int : std::stod(s);
}

/**
 * Parse comma-separated string of minimizer offsets into a vector of numbers
 * returning their total to be used as sanity check
 */
inline int64_t parse_minimizers(const std::string& buffer, vector<int32_t> offsets) {
    int64_t span = 0;
    int64_t i = 0;
    for (int64_t j = 0; j < buffer.size(); ++j) {
        if (buffer[j] == ',') {
            assert(j > i);
            offsets.push_back(std::stoi(buffer.substr(i, j-1)));
            span += offsets.back();
            i = j + 1;
        }
    }
    assert(i <= buffer.length() - 1);
    offsets.push_back(std::stoi(buffer.substr(i)));
    span += offsets.back();
    return span;
}

/**
 * Parse a mz GAF line into a struct
 */
inline void parse_mzgaf_record(const std::string& gaf_line, MzGafRecord& gaf_record) {

    std::istringstream in(gaf_line);
    std::string buffer;
    
    int col = 1;
    
    std::function<void(void)> scan_column = [&]() {
        getline(in, buffer, '\t');
        if (!in || buffer.empty()) {
            throw std::runtime_error("Error parsing GAF column " + std::to_string(col));
        }
        ++col;
    };
    
    scan_column();
    assert(buffer == "*");
    
    scan_column();
    assert(buffer[0] == '<' || buffer[0] == '>');
    gaf_record.target_name = buffer.substr(1);
    gaf_record.is_reverse = buffer[0] == '<';

    scan_column();
    gaf_record.target_length = string_to_int(buffer);

    scan_column();
    gaf_record.num_minimizers = string_to_int(buffer);

    scan_column();
    gaf_record.seq_div = string_to_float(buffer);
    
    scan_column();
    gaf_record.target_start = string_to_int(buffer);
    gaf_record.target_end = string_to_int(buffer);

    scan_column();
    gaf_record.query_start = string_to_int(buffer);
    gaf_record.query_end = string_to_int(buffer);

    scan_column();
    gaf_record.target_mz_offsets.clear();
    if (gaf_record.num_minimizers > 0) {
        gaf_record.target_mz_offsets.reserve(gaf_record.num_minimizers);
    }
    int64_t span = parse_minimizers(buffer, gaf_record.target_mz_offsets);
    assert(gaf_record.target_mz_offsets.size() == gaf_record.num_minimizers);
    assert(span + gaf_record.kmer_size == gaf_record.target_length);

    scan_column();
    gaf_record.query_mz_offsets.clear();
    if (gaf_record.num_minimizers > 0) {
        gaf_record.query_mz_offsets.reserve(gaf_record.num_minimizers);
    }
    span = parse_minimizers(buffer, gaf_record.query_mz_offsets);
    assert(gaf_record.query_mz_offsets.size() == gaf_record.num_minimizers);
    assert(span + gaf_record.kmer_size == gaf_record.query_length);
}

/**
 *  Visit every mz GAF record in a file, and do a callback.  The regulare GAF records 
 *  are ignored except for the query name (though they are parsed, and could maybe get 
 *  used for filtering or validation)
 */
inline void scan_mzgaf(istream& in_stream, function<void(MzGafRecord& gaf_record)> visit_fn) {
    string line_buffer;
    GafRecord gaf_record;
    MzGafRecord mz_record;
    while (getline(in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            assert(!mz_record.query_name.empty());
            parse_mzgaf_record(line_buffer, mz_record);
            visit_fn(mz_record);
        } else {
            parse_gaf_record(line_buffer, gaf_record);
            mz_record.query_name = gaf_record.query_name;
            mz_record.query_length = gaf_record.query_length;
        }
    }
}

}
