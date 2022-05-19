#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <map>
using namespace std;

struct PafLine {
    string query_name;
    int64_t query_len;
    int64_t query_start;
    int64_t query_end;
    char strand;
    string target_name;
    int64_t target_len;
    int64_t target_start;
    int64_t target_end;
    int64_t num_matching;
    int64_t num_bases;
    int64_t mapq;
    string cigar;

    // Map a tag name to its type and value
    // ex: "de:f:0.2183" in the GAF would appear as opt_fields["de"] = ("f", "0.2183")
    // note: cigar not stored here, but rather in cigar string above
    std::map<std::string, std::pair<std::string, std::string>>  opt_fields;
};

inline vector<string> split_delims(const string &s, const string& delims, vector<string> &elems) {
    size_t start = string::npos;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) != string::npos) {
            if (start != string::npos && i > start) {
                elems.push_back(s.substr(start, i - start));
            }
            start = string::npos;
        } else if (start == string::npos) {
            start = i;
        }
    }
    if (start != string::npos && start < s.size()) {
        elems.push_back(s.substr(start, s.size() - start));
    }
    return elems;
}

inline PafLine parse_paf_line(const string& paf_line) {
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);
    assert(toks.size() > 12);

    PafLine paf;
    paf.query_name = toks[0];
    paf.query_len = stol(toks[1]);
    paf.query_start = stol(toks[2]);
    paf.query_end = stol(toks[3]);
    assert(toks[4] == "+" || toks[4] == "-");
    paf.strand = toks[4][0];
    paf.target_name = toks[5];
    paf.target_len = stol(toks[6]);
    paf.target_start = stol(toks[7]);
    paf.target_end = stol(toks[8]);
    paf.num_matching = stol(toks[9]);
    paf.num_bases = stol(toks[10]);
    paf.mapq = stol(toks[11]);

    for (size_t i = 12; i < toks.size(); ++i) {
        if (toks[i].compare(0, 3, "cg:Z:") == 0) {
            paf.cigar = toks[i].substr(5);
        } else {
            vector<string> tag_toks;
            split_delims(toks[i], ":", tag_toks);
            assert(tag_toks.size() == 3);
            paf.opt_fields[tag_toks[0]] = make_pair(tag_toks[1], tag_toks[2]);
        }
    }

    return paf;
}

inline ostream& operator<<(ostream& os, const PafLine& paf) {
    os << paf.query_name << "\t" << paf.query_len << "\t" << paf.query_start << "\t" << paf.query_end << "\t"
       << string(1, paf.strand) << "\t"
       << paf.target_name << "\t" << paf.target_len << "\t" << paf.target_start << "\t" << paf.target_end << "\t"
       << paf.num_matching << "\t" << paf.num_bases << "\t" << paf.mapq;
    if (!paf.cigar.empty()) {
        os << "\tcg:Z:" << paf.cigar;        
    }
    for (const auto& kv : paf.opt_fields) {
        os << "\t" << kv.first << ":" << kv.second.first << ":" << kv.second.second;
    }
    return os;
}

inline void for_each_cg(const string& cg_tok, function<void(const string&, const string&)> fn) {
    size_t next;
    for (size_t co = 5; co != string::npos; co = next) {
        next = cg_tok.find_first_of("M=XDI", co + 1);
        if (next != string::npos) {
            fn(cg_tok.substr(co, next - co), cg_tok.substr(next, 1));
            ++next;
        }
    }
}
