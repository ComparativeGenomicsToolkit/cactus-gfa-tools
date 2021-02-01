/**
 * rgfa-split.hpp: Get a mapping from minigraph rGFA node ids to reference contigs, and use to partition paf
 */

#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

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
               const unordered_map<int64_t, int64_t>& contig_map,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_prefix,
               const string& minigraph_prefix, // this is the cactus unique id prefix (ex id=0|)
               double min_query_coverage,
               double min_small_query_coverage,
               int64_t small_coverage_threshold,
               double min_query_uniqueness,
               int64_t ambiguous_id,
               const string& reference_prefix,
               const unordered_map<string, int64_t>& mask_stats); 

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



