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

/* array of 2-bit values based on vector<bool>
 * (too lazy to import bit vector library)
 */
class TwoBitVec {
public:
    TwoBitVec() = default;
    ~TwoBitVec() = default;
    void resize(size_t n) {
        v1.resize(n, false);
        v2.resize(n, false);
    }
    void set(size_t pos, size_t val) {
        assert(val < 4);
        v1[pos] = val & 1u;
        v2[pos] = val & 2u;
    }
    uint8_t get(size_t pos) {
        uint8_t val = 0;
        val |= 1u * v1[pos];
        val |= 2u * v2[pos];
        return val;
    }
    void increment(size_t pos) {
        uint8_t v = get(pos) + 1;
        if (v < 4) {
            set(pos, v);
        }
    }
    void clear() {
        v1.clear();
        v2.clear();
    }
    size_t size() {
        return v1.size();
    }
protected:
    vector<bool> v1;
    vector<bool> v2;
};

typedef unordered_map<string, TwoBitVec> QueryCoverage;

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
                 QueryCoverage& query_coverage,
                 int64_t min_overlap_len,
                 const string& target_prefix = "");

/* update the counts for one mapping of query to target.  if the any of the filters
 * fail, the coverage count is updated, but none of the minimizer counts.
 */
void update_mz_map(const gafkluge::MzGafRecord& gaf_record,
                   const gafkluge::GafRecord& parent_record,
                   MZMap& mz_map,
                   int64_t min_mapq,
                   int64_t min_block_len,
                   int64_t min_node_len,
                   bool node_based_universal);

/* combine the two maps, adding all elements of map1 into map2 *and removing* them from map1
 */
void combine_mz_maps(MZMap& map1, MZMap& map2, bool reset_multiple_counts_to_0);

/*
 * Update the query coverage map.  It stores if a query base has been covered 0 1 or 2 times. 
 */
void update_query_coverage(const gafkluge::GafRecord& parent_record,
                           QueryCoverage& query_coverage);
