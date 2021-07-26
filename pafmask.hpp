#pragma once

#include <fstream>
#include <iostream>
#include <sstream>

#include "pafcoverage.hpp"
#include "rgfa-split.hpp"

//#define debug

using namespace std;

// make a map of sequence name -> interval tree from a bed file
// intervals get merged up if they overlap or are within padding of each other
unordered_map<string, CoverageIntervalTree> load_bed(istream& bed_stream, int64_t padding);

// cut interval_b out of interval_a, any peices of a that remain get added to out_fragments
// the number of such pieces is 0 (a in b), 1 (b overlaps one end of a) or 2 (b in a)
void interval_subtract(CoverageInterval& interval_a, CoverageInterval& interval_b,
                              vector<CoverageInterval>& out_fragments);

// subtract all masked intervals from a paf line and output what's left
size_t mask_paf_line(const string& paf_line, int64_t min_length, const unordered_map<string,
                          CoverageIntervalTree>& ref_to_intervals, bool validate);

// output paf line(s) corresponding to given sub-interval of the original paf line
void clip_paf(const vector<string>& toks, const string& query_name, int64_t query_length, int64_t query_start, int64_t query_end,
                     CoverageInterval& interval, int64_t min_length, bool validate);

// make sure every homology in the fragment_paf corresponds to a homology in toks
void validate_paf(const vector<string>& toks, const string& fragment_paf);
