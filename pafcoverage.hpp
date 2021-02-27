/**
 * pafcoverage.hpp: Get some stats from cigar pafs in order to see what kind of anchors they make for Cactus
 *                 (not done within mzgaf2paf so it can be used on paf's from other sources too)
 */


#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <functional>

using namespace std;

// map a sequence name to all its covered bases (not caring about depth, so just using bools)
typedef unordered_map<string, vector<bool>> CoverageMap;

/** update bases covered in the query sequence of a paf line
 */
void update_coverage_map(const string& paf_line, CoverageMap& coverage_map);

/** print some stats
 */
void print_coverage_summary(const CoverageMap& coverage_map, ostream& out);

/** print bed of coverage gaps
 */
void print_coverage_gaps_as_bed(const CoverageMap& coverage_map, ostream& out, int64_t min_gap_length);

// some parsing functions more or less copied from vg
vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems);

// parse cigar
void for_each_cg(const string& cg_tok, function<void(const string&, const string&)> fn);
