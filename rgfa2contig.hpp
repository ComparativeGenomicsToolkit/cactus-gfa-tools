/**
 * rgfa2contig.hpp: Get a mapping from minigraph rGFA node ids to reference contigs
 */


#pragma once
#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

/**
 * Returns a map from node-id (string, ex S1) to reference contig id (int32_t)
 * as well as a map of reference contig id to reference contig name (ex chr2) (vector<string>)
 */
pair<unordered_map<int64_t, int64_t>, vector<string>> rgfa2contig(const string& gfa_path);

/**
 * Map from minigraph string ID to numeric ID assuming S<ID> naming convention
 */
inline int64_t node_id(const string& rgfa_id) {
    return stol(rgfa_id.substr(1));
}

