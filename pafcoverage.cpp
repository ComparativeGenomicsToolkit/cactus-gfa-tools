#include <functional>
#include <cassert>
#include <iostream>
#include "pafcoverage.hpp"

using namespace std;

void update_coverage_map(const string& paf_line, CoverageMap& coverage_map) {

    // split into array of tokens
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);

    if (toks.size() < 12) {
        throw runtime_error("too few tokens in PAF line: " + paf_line);
    }

    string& query_name = toks[0];
    int64_t query_length = stol(toks[1]);    
    
    vector<bool>& query_coverage = coverage_map[query_name];
    if (query_coverage.empty()) {
        query_coverage.resize(query_length, false);
    } 
    assert(query_coverage.size() == query_length);

    for (int i = 12; i < toks.size(); ++i) {
        if (toks[i].substr(0, 5) == "cg:Z:") {
            // note: query coordinates always on forward strand in paf
            int64_t query_pos = stol(toks[2]);
            for_each_cg(toks[i], [&](const string& val, const string& cat) {
                    int64_t len = stol(val);
                    if (cat == "M" || cat == "=" || cat == "X") {
                        for (int64_t j = 0; j < len; ++j) {
                            query_coverage[query_pos + j] = true;
                        }
                    }
                    if (cat == "M" || cat == "=" || cat == "X" || cat == "I") {
                        query_pos += len;
                    }
                });
        }
    }
}

void print_coverage_summary(const CoverageMap& coverage_map, ostream& out) {
    out << "query-name" << "\t" << "pct-coverage" << "\t" << "max-gap" << "\t" << "avg-gap" << endl;
    out << "----------" << "\t" << "------------" << "\t" << "-------" << "\t" << "-------" << endl;
    for (const auto& kv : coverage_map) {
        const string& query_name = kv.first;
        const vector<bool>& coverage = kv.second;
        int64_t last_covered = -1;
        int64_t coverage_count = 0;
        vector<int64_t> gaps;

        for (int64_t i = 0; i < coverage.size(); ++i) {
            if (coverage[i] == true) {
                ++coverage_count;
                if (i - last_covered > 1) {
                    gaps.push_back(i - last_covered - 1);
                }
                last_covered = i;
            }
        }

        int64_t i = coverage.size();
        if (i - last_covered > 1) {
            gaps.push_back(i - last_covered - 1);
        }

        int64_t max_gap = 0;
        int64_t total_gap = 0;
        for (int i = 0; i < gaps.size(); ++i) {
            max_gap = max(gaps[i], max_gap);
            total_gap += gaps[i];
        }

        out << query_name << "\t"
            << ((float)coverage_count / coverage.size()) << "\t"
            << max_gap << "\t"
            << (!gaps.empty() ? (total_gap / gaps.size()) : 0) << endl;
    }    
}

void print_coverage_gaps_as_bed(const CoverageMap& coverage_map, ostream& out, int64_t min_gap_length) {
    for (const auto& kv : coverage_map) {
        const string& query_name = kv.first;
        const vector<bool>& coverage = kv.second;
        int64_t last_covered = -1;

        for (int64_t i = 0; i < coverage.size(); ++i) {
            if (coverage[i] == true) {
                if (i - last_covered > min_gap_length) {
                    out << query_name << "\t" << (last_covered + 1) << "\t" << i << "\t" << "pafcoverage-m" << min_gap_length << endl;
                }
                last_covered = i;
            }
        }

        int64_t i = coverage.size();
        if (i - last_covered > min_gap_length) {
            out << query_name << "\t" << (last_covered + 1) << "\t" << i << "\t" << "pafcoverage-m" << min_gap_length << endl;
        }
    }
}
