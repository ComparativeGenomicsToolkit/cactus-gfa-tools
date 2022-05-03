#include <sstream>
#include <vector>
#include <functional>
#include <cassert>
#include <iostream>
#include "paf2lastz.hpp"
#include "paf.hpp"

using namespace std;

pair<std::string, bool> paf2lastz(const std::string& paf_line, bool use_mapq) {

    // split into array of tokens
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);

    // handle empty line
    if (toks.size() == 0) {
        return make_pair("", false);
    }
    
    if (toks.size() < 12) {
        throw runtime_error("[paf2lastz] error: too few tokens in PAF line: " + paf_line);
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
                    lastz_line += " " + (cat == "X" || cat == "=" ? "M" : cat) + " " + val;
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
    
