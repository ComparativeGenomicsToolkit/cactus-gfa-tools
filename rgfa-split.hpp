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

/**
 * Use contigs identified above to split PAF
 */
void paf_split(istream& input_paf_stream,
               const unordered_map<int64_t, int64_t>& contig_map,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_prefix);

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
    return stol(rgfa_id.substr(1));
}



