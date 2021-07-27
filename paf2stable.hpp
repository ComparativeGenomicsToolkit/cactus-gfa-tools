/**
 * paf2stable.hpp: Use transitivity to remove all targets from a PAF file, and replace with queries
 *                 Used to drop minigraph node sequences out of an assembly-to-minigraph mapping file
 */


#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <functional>
#include <set>
#include "IntervalTree.h"

using namespace std;

// maps an interval in a minigraph node sequence to <QUERY_ID, QUERY_START, IS_REVERSED>
// note that QUERY_END can inferred from the interval's length
typedef IntervalTree<int64_t, tuple<int64_t, int64_t, bool>> StableIntervalTree;
typedef StableIntervalTree::interval StableInterval;

void update_stable_mapping_info(const vector<string>& paf_toks,
                                unordered_map<string, int64_t>& query_name_to_id,
                                vector<pair<string, int64_t>>& query_id_to_info,
                                unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals);


// build an interval tree for each target out of the mapped intervals
// the input map is emptied in the process
//
// UPDATE: no longer using interval trees as they are too slow to query.  now just keep vectors in sorted order!
//
void create_interval_trees(unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals);

// clip out an interval against a set, updating the clipped_intervals list
void clip_interval(const StableInterval& interval,
                   const set<int64_t>& cut_points, vector<StableInterval>& clipped_intervals);

// apply the interval map to a PAF line to rewrite it as a series of query-to-query mappings
// output PAF line(s) prtined to cout and return number of lines
size_t paf_to_stable(const vector<string>& paf_toks,
                     const vector<pair<string, int64_t>>& query_id_to_info,
                     const unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals);

// write a new paf line (in stable coordinates) for a given sub-interval of a match block
// return number of lines written
void make_paf_line_for_interval(const vector<string>& paf_toks,
                                const vector<pair<string, int64_t>>& query_id_to_info,
                                const StableInterval& overlapping_interval,
                                int64_t query_pos);

// for debugging
inline ostream& operator<<(ostream& os, const tuple<int64_t, int64_t, bool>& triple) {
    os << "<" << get<0>(triple) << ", " << get<1>(triple) << ", " << get<2>(triple) << ">";
    return os;
}
