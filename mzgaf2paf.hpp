/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */


#pragma once
#include "mzgaf.hpp"
#include <iostream>
#include <unordered_map>

/** offset -> <number of minimzers, number of mappings> of given target sequence
 */
typedef vector<pair<uint16_t, uint16_t>> MZCount;

/** map a target sequence name to a pair of the number of times its been mapped to, along with
 * the number of times each of its minimizers has been mapped to
 */
typedef unordered_map<string, MZCount> MZMap;

/**
 * Conevert mzgaf to paf (returns total length of all matches)
 */
size_t mzgaf2paf(const gafkluge::MzGafRecord& gaf_record,
               const gafkluge::GafRecord& parent_record,
               ostream& paf_stream,
               int64_t min_match_length,
               int64_t min_gap,
               MZMap& mz_map,
               double universal_filter,
               const string& target_prefix = "");

/* update the counts for one mapping of query to target
 */
void update_mz_map(const gafkluge::MzGafRecord& gaf_record,
                   const gafkluge::GafRecord& parent_record,
                   MZMap& mz_map,
                   bool node_based_universal);

/* combine the two maps, adding all elements of map1 into map2 *and removing* them from map1
 */
void combine_mz_maps(MZMap& map1, MZMap& map2, bool reset_multiple_counts_to_0);
