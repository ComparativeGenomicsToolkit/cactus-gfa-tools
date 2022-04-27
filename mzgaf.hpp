/**
 * mzgaf.hpp: Functions to parse output from minigraph --write-mz which outputs GAF-like records that contain minimizer hits
 *          mapping the query to the target.  Rely on gafkluge to read the normal GAF lines, and add here some functions
 *          to read the minimizer lines. 
 */


#pragma once
// copied from https://github.com/vgteam/libvgio/blob/master/include/vg/io/gafkluge.hpp
#include <cassert>
#include <list>
#include "gafkluge.hpp"

using namespace std;

namespace gafkluge {

/** from minigraph -S --write-mz
 *      segName segLen  nMinis  seqDiv  segStart segEnd qStart  qEnd    k-mer   segOffsets                      qOffsets
 *	>s79	124	15	0.0096	7	 118	193184	193299	19	8,9,3,8,9,11,9,4,2,7,2,9,3,8	8,9,3,8,9,15,9,4,2,7,2,9,3,8
 */
struct MzGafRecord {
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
    vector<pair<char, int64_t>> cigar;
};

inline double string_to_float(const std::string& s) {
    return s == missing_string ? (double)missing_int : std::stod(s);
}

/**
 * Parse comma-separated string of minimizer offsets into a vector of numbers
 * returning their total to be used as sanity check
 */
inline int64_t parse_minimizers(const std::string& buffer, vector<int32_t>& offsets) {
    int64_t span = 0;
    int64_t i = 0;
    for (int64_t j = 0; j < buffer.size(); ++j) {
        if (buffer[j] == ',') {
            assert(j > i);
            offsets.push_back(std::stoi(buffer.substr(i, j-i)));
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
    gaf_record.target_length = std::stol(buffer);

    scan_column();
    gaf_record.num_minimizers = std::stol(buffer);

    if (gaf_record.num_minimizers == 0) {
        // it turns out the remaining columns are optional
        gaf_record.target_start = missing_int;
        gaf_record.target_end = missing_int;
        gaf_record.query_start = missing_int;
        gaf_record.query_end = missing_int;
        gaf_record.kmer_size = missing_int;
        return;
    }
        
    scan_column();
    gaf_record.seq_div = string_to_float(buffer);
    
    scan_column();
    gaf_record.target_start = std::stol(buffer);
    scan_column();
    gaf_record.target_end = std::stol(buffer);

    scan_column();
    gaf_record.query_start = std::stol(buffer);
    scan_column();
    gaf_record.query_end = std::stol(buffer);

    gaf_record.target_mz_offsets.clear();
    gaf_record.query_mz_offsets.clear();

    if (in.peek() == EOF) {
        // we now support -S without --write-mz (but with -c).  in this case, we're at the end
        gaf_record.kmer_size = missing_int;
        return;
    }
    
    scan_column();
    gaf_record.kmer_size = std::stol(buffer);
    scan_column();
    if (gaf_record.num_minimizers > 0) {
        gaf_record.target_mz_offsets.reserve(gaf_record.num_minimizers);
    }
    int64_t span = parse_minimizers(buffer, gaf_record.target_mz_offsets);
    assert(gaf_record.target_mz_offsets.size() + 1 == gaf_record.num_minimizers);
    assert(span + gaf_record.kmer_size == gaf_record.target_end - gaf_record.target_start);

    scan_column();
    if (gaf_record.num_minimizers > 0) {
        gaf_record.query_mz_offsets.reserve(gaf_record.num_minimizers);
    }
    span = parse_minimizers(buffer, gaf_record.query_mz_offsets);
    assert(gaf_record.query_mz_offsets.size() + 1 == gaf_record.num_minimizers);
    assert(span + gaf_record.kmer_size == gaf_record.query_end - gaf_record.query_start);
}

/**
 *  Visit every mz GAF record in a file, and do a callback.  The regular GAF records 
 *  are ignored except for the query name (though they are parsed, and could maybe get 
 *  used for filtering or validation)
 */
inline void scan_mzgaf(istream& in_stream, function<void(MzGafRecord& gaf_record, GafRecord& parent_record)> visit_fn,
                       function<void(GafRecord& parent_record)> parent_fn = nullptr) {
    string line_buffer;
    GafRecord gaf_record;
    list<pair<char, size_t>> gaf_cigar;
    list<pair<char, size_t>>::iterator gaf_cig_it;
    MzGafRecord mz_record;
    while (getline(in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            assert(!gaf_record.query_name.empty());
            parse_mzgaf_record(line_buffer, mz_record);
            // scan the cigar to consume the step
            if (!gaf_cigar.empty()) {
                int64_t cig_target_len = 0;
                auto gaf_cig_it2 = gaf_cig_it;
                int64_t target_length = mz_record.target_end - mz_record.target_start;
                for (; gaf_cig_it2 != gaf_cigar.end() && cig_target_len < target_length; ++gaf_cig_it2) {
                    if (gaf_cig_it2->first == 'X' || gaf_cig_it2->first == '=' || gaf_cig_it2->first == 'D'
                        || gaf_cig_it2->first == 'M' || gaf_cig_it2->first == 'N') {
                        cig_target_len += gaf_cig_it2->second;
                    }
                }
                if (cig_target_len < target_length) {
                    cerr << "[mzgaf] error: Ran out of cigar to consume GAF step\n" << gaf_record << endl << line_buffer << endl;
                    exit(1);
                }
                // if our cigar overhangs the step we just chop it
                // maintaining that the step exactly spans range [it, it2)
                if (cig_target_len > target_length) {
                    int64_t delta = cig_target_len - target_length;
                    // insert a new (empty) step before it2
                    gaf_cigar.insert(gaf_cig_it2, pair<char, size_t>());
                    --gaf_cig_it2;

                    // set the new step's length to delta
                    auto gaf_cig_it3 = gaf_cig_it2;
                    --gaf_cig_it3;
                    gaf_cig_it2->first = gaf_cig_it3->first;
                    gaf_cig_it2->second = delta;
                    // and subtract it off the previous step
                    assert(gaf_cig_it3->second > delta);
                    gaf_cig_it3->second -= delta;
                }
                mz_record.cigar.clear();
                for (; gaf_cig_it != gaf_cig_it2; ++gaf_cig_it) {
                    mz_record.cigar.push_back(*gaf_cig_it);
                }
            }            
            visit_fn(mz_record, gaf_record);
        } else {
            parse_gaf_record(line_buffer, gaf_record);
            // parse the cigar if it's there
            gaf_cigar.clear();
            for_each_cg(gaf_record, [&](const char& c, const size_t& s) {
                    gaf_cigar.push_back(make_pair(c, s));
                });
            gaf_cig_it = gaf_cigar.begin();
            if (parent_fn) {
                parent_fn(gaf_record);
            }
        }
    }
}

}
