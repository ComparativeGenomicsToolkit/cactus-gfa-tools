/**
 * mzgaf.hpp: Functions to parse output from minigraph --write-mz which outputs GAF-like records that contain minimizer hits
 *          mapping the query to the target.  Rely on gafkluge to read the normal GAF lines, and add here some functions
 *          to read the minimizer lines. 
 */


#pragma once
#include <cassert>
#include <unordered_map>
#include <limits>
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

    scan_column();
    gaf_record.kmer_size = std::stol(buffer);
    scan_column();
    gaf_record.target_mz_offsets.clear();
    if (gaf_record.num_minimizers > 0) {
        gaf_record.target_mz_offsets.reserve(gaf_record.num_minimizers);
    }
    int64_t span = parse_minimizers(buffer, gaf_record.target_mz_offsets);
    assert(gaf_record.target_mz_offsets.size() + 1 == gaf_record.num_minimizers);
    assert(span + gaf_record.kmer_size == gaf_record.target_end - gaf_record.target_start);

    scan_column();
    gaf_record.query_mz_offsets.clear();
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
                       function<void(GafRecord& parent_record)> parent_fn = nullptr,
                       const unordered_map<string, tuple<string, string, size_t>>* fa_header_table = nullptr) {
    string line_buffer;
    GafRecord gaf_record;
    MzGafRecord mz_record;
    size_t step_no = 0;
    size_t target_pos = 0;
    bool single_contig = false; // just a single contig like "chr1"
    size_t step_offset = 0; // where we are in the step's stable sequence
    MzGafRecord mz_copy;
    while (getline(in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            assert(!gaf_record.query_name.empty());
            parse_mzgaf_record(line_buffer, mz_record);
            if (fa_header_table) {
                // override node coordinates with stable coordinate from parent's path step
                // this only works under the assumption that there is one step per * line
                assert(step_no < gaf_record.path.size());
                auto table_it = fa_header_table->find(gaf_record.path[step_no].name);
                if (table_it == fa_header_table->end()) {
                    throw runtime_error("Unable to find contig " + gaf_record.path[step_no].name + " in header table");
                }
                mz_copy.target_name = mz_record.target_name;
                mz_copy.target_length = mz_record.target_length;
                mz_copy.target_start = mz_record.target_start;
                mz_copy.target_end = mz_record.target_end;
                mz_record.target_name = "id=" + get<1>(table_it->second) + "|" + get<0>(table_it->second);
                mz_record.target_length = get<2>(table_it->second);
                mz_record.target_start += gaf_record.path[step_no].start + step_offset;
                mz_record.target_end = mz_record.target_start + (mz_copy.target_end - mz_copy.target_start);
                assert(mz_record.target_start <= mz_record.target_end);
                assert(mz_record.target_start >= 0);
                assert(mz_record.target_start < mz_record.target_length);
                assert(mz_record.target_end <= mz_record.target_length);
            }
            visit_fn(mz_record, gaf_record);
            if (fa_header_table) {
                // put everything back to the way it was before overriding
                mz_record.target_name = mz_copy.target_name;
                mz_record.target_length = mz_copy.target_length;
                mz_record.target_start = mz_copy.target_start;
                mz_record.target_end = mz_copy.target_end;                
            }
            // we move forward in the path by the length (second column of the mz line)
            step_offset += mz_record.target_length;
#ifdef debug
            cerr << "added " << mz_record.target_length << " to get " << step_offset;            
            cerr << " our current contig size is " << gaf_record.path[step_no].end << " with step no " << step_no << endl;
#endif
            assert(mz_record.target_end - mz_record.target_start <=  mz_record.target_length);
            if (gaf_record.path[step_no].start + step_offset == gaf_record.path[step_no].end) {
                // if we're off the step, advance
                ++step_no;
                if (step_no < gaf_record.path.size()) {
                    step_offset = 0;
                } else {
                    step_offset = numeric_limits<size_t>::max();
                }
            } else if (!single_contig) {
                assert (step_offset < gaf_record.path[step_no].end);
            }
        } else {
            step_offset = 0;
            if (step_no > 0) {
                assert(step_no == gaf_record.path.size());
            }
            parse_gaf_record(line_buffer, gaf_record);
            target_pos = gaf_record.path_start;
            single_contig = gaf_record.path.size() == 1 && gaf_record.path[0].is_stable && gaf_record.path[0].is_interval == false;
            if (parent_fn) {
                parent_fn(gaf_record);
            }
            step_no = 0;
        }
    }
}

}
