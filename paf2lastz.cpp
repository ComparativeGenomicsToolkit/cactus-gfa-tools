#include <sstream>
#include <vector>
#include <functional>
#include <cassert>
#include <iostream>
#include "paf2lastz.hpp"

using namespace std;

// some parsing functions more or less copied from vg
static std::vector<std::string> &split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
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

static void for_each_cg(const string& cg_tok, std::function<void(const std::string&, const std::string&)> fn) {
    size_t next;
    for (size_t co = 5; co != std::string::npos; co = next) {
        next = cg_tok.find_first_of("MDI", co + 1);
        if (next != std::string::npos) {
            fn(cg_tok.substr(co, next - co), cg_tok.substr(next, 1));
            ++next;
        }
    }
}

pair<std::string, bool> paf2lastz(const std::string& paf_line, bool use_mapq) {

    // split into array of tokens
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);

    if (toks.size() < 12) {
        throw runtime_error("too few tokens in PAF line: " + paf_line);
    }

    // use PAF mapq as the lastz score
    string lz_score = toks[11];
    if (!use_mapq) {
        // or take it from the optional score field
        lz_score = "0";
        for (int i = 12; i < toks.size(); ++i) {
            if (toks[i].substr(0, 5) == "AS:i:") {
                lz_score = toks[i].substr(5);
                break;
            }
        }
    }

    assert(toks[4] == "+" || toks[4] == "-");
    // hat-tip: @Robin-Rounthwaite
    // https://github.com/Robin-Rounthwaite/reference-based-cactus-aligner/blob/master/src/paf_to_lastz.py#L49-L71
    if (toks[4] == "-") {
        std::swap(toks[2], toks[3]);
    }

    string lastz_line = "cigar: " +
        toks[0] + " " +
        toks[2] + " " +
        toks[3] + " " +
        toks[4] + " " +
        toks[5] + " " +
        toks[7] + " " +
        toks[8] + " " +
        "+ " +
        lz_score;

    // do the cigar
    bool found_cigar = false;
    bool is_secondary = false;
    for (int i = 12; i < toks.size(); ++i) {
        if (toks[i].substr(0, 5) == "cg:Z:") {
            found_cigar = true;
            for_each_cg(toks[i], [&](const string& val, const string& cat) {
                    lastz_line += " " + cat + " " + val;
                });
        } else if (toks[i].substr(0, 5) == "tp:A:") {
            is_secondary = toks[i].length() == 6 && toks[i][5] == 'S';
        }
    }

    if (!found_cigar) {
        cerr << "Warning: cg tag not found on PAF line: " << paf_line << endl;
    }

    return make_pair(lastz_line, is_secondary);
}
    
