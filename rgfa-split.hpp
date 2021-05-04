/**
 * rgfa-split.hpp: Get a mapping from minigraph rGFA node ids to reference contigs, and use to partition paf
 */

#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include "IntervalTree.h"

using namespace std;

/**
 * Returns a map from node-id (string, ex S1) to reference contig id (int32_t)
 * as well as a map of reference contig id to reference contig name (ex chr2) (vector<string>)
 */
pair<unordered_map<int64_t, int64_t>, vector<string>> rgfa2contig(const string& gfa_path);

/*
 * Load table (2 column tsv) of above data
 */
pair<unordered_map<int64_t, int64_t>, vector<string>> load_contig_map(const string& contgs_path);

/*
 * Load bed file into coverage map to keep track of masked intervals
 */
unordered_map<string, int64_t> load_query_mask_stats(const string& bed_path);

/*
 * Reorganize mapping so all unselected reference contigs map to an "other" category of given name
 */
void set_other_contig(unordered_map<int64_t, int64_t>& contig_map,
                      vector<string>& contigs,
                      function<bool(const string&)> visit_contig,
                      const string& other_name);

/**
 * Use contigs identified above to split PAF
 */
void paf_split(const string& input_paf_path,
               function<int64_t(const string&)> name_to_refid,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_prefix,
               double min_query_coverage,
               double min_small_query_coverage,
               int64_t small_coverage_threshold,
               double min_query_uniqueness,
               int64_t min_query_chunk,
               bool allow_softclip,
               int64_t ambiguous_id,
               const string& reference_prefix,
               const unordered_map<string, int64_t>& mask_stats,
               int64_t max_gap_as_match,
               int64_t min_mapq,
               ostream& log_stream); 

/**
 * Use the contigs to split the GFA
 */
void gfa_split(const string& rgfa_path,
               const unordered_map<int64_t, int64_t>& contig_map,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_gfa_prefix);

/**
 * Map from minigraph string ID to numeric ID assuming S<ID> naming convention
 */
inline int64_t node_id(const string& rgfa_id) {
    // cactus may have added some nonsense -- search past it:
    size_t offset = rgfa_id.find('s') + 1;
    return stol(rgfa_id.substr(offset));
}

/**
 * Count up bases of small indels bookended by matches, which we will apply to the coverage
 */
int64_t count_small_gap_bases(const vector<string>& toks, int64_t max_gap_as_match);

/**
 * Keep track of PAF coverage by remembering intervals (generalizes previous logic that
 * just counted bases.  The value here is pair<int64_t, int64_t> = <aligned bases, reference contig>
 */
typedef IntervalTree<int64_t, pair<int64_t, int64_t>> CoverageIntervalTree;
typedef CoverageIntervalTree::interval CoverageInterval;
inline ostream& operator<<(ostream& os, pair<int64_t, int64_t> v) {
    os << "(" << v.first << ", " << v.second << ")";
    return os;
}
void scan_coverage_intervals(CoverageIntervalTree& interval_tree, int64_t padding, function<void(int64_t, int64_t, pair<int64_t, int64_t>)> fn);

/**
 * Smooth out the individual paf interval mappings from a given query to try to come up with some kind of consensus assignment
 * 
 */
void smooth_query_intervals(const string& query_name, int64_t query_length, int64_t masked_bases,
                            vector<CoverageInterval> & intervals, double min_coverage, double min_uniqueness,
                            int64_t min_chunk, const vector<string>& ref_contigs, bool allow_softclip, ostream& log_stream);

/**
 * Rename the query contig to a sub-fragment in order to reflect the fact that it will be cut in the output
 */
void apply_paf_query_offsets(vector<string>& paf_toks, int64_t query_fragment_start, int64_t query_fragment_end); 


// TODO: make consisten with vg's naming scheme (requires change in grpahmap-split that would postprocess faidx output)

/**
 * Turn chr1:10-100 into <"chr1", 9, 99>
 */
tuple<string, int64_t, int64_t>  parse_faidx_subpath(const string& name);


/**
 * Turn <"chr1", 9, 99> into chr1:10-100
 */
string make_faidx_subpath(const string& name, int64_t start, int64_t end);
