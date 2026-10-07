#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <map>
#include <list>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include "gafkluge.hpp"
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


// ---- gaf2paf's split of a GAF record into PAF lines -----------------------------------------
// gaf2paf turns each step of a record's path into one PAF line.  The code lives here, rather than
// in gaf2paf_main.cpp, because gaffilter -x has to know exactly which lines gaf2paf will print for a
// record (it judges what the stock filter chain finally anchors), and the only way to be sure of
// that is to run the same code.  gaf2paf_main.cpp prints what these functions produce.
namespace gaf2paf_split {

typedef pair<char, int64_t> Cig;
typedef list<Cig> Cigar;

inline bool consumes_query(const Cig& c) {
    return c.first == 'M' || c.first == 'I' || c.first == 'S' || c.first == '=' || c.first == 'X';
}

inline bool consumes_target(const Cig& c) {
    return  c.first == 'M' || c.first == 'D' || c.first == 'N' || c.first == '=' || c.first == 'X';
}

// cut the cigar at pos removing cut_len. cut_len is inserted as a new element after pos
inline void cigar_cut(Cigar& cigar, Cigar::iterator pos, int64_t cut_len) {
    assert (cut_len > 0);
    int64_t remainder = pos->second - cut_len;
    assert(remainder > 0);
    auto pos2 = pos;
    ++pos2;
    Cig new_item = make_pair(pos->first, cut_len);
    cigar.insert(pos2, new_item);
    pos->second = remainder;
};

// get the next "target_len" bases worth of cigar, starting at pos, clipping at the end if necessary
inline pair<Cigar::iterator, Cigar::iterator> cigar_next_by_target(Cigar& cigar, Cigar::iterator pos, int64_t target_len) {
    int64_t cur_len = 0;
    auto pos2 = pos;
    for (pos2 = pos; pos2 != cigar.end() && cur_len < target_len; ++pos2) {
        if (consumes_target(*pos2)) {
            cur_len += pos2->second;
        }
    }
    if (cur_len != target_len) {
        assert(cur_len > target_len);
        int64_t cut_len = cur_len - target_len;
        --pos2;
        cur_len -= pos2->second;
        cigar_cut(cigar, pos2, cut_len);
        cur_len += pos2->second;
        ++pos2;
    }
    assert(cur_len == target_len);
    return make_pair(pos, pos2);
}

inline void flip_gaf(gafkluge::GafRecord& gaf_record, const unordered_map<string, int64_t>& len_map) {
    // flip strand
    gaf_record.strand = gaf_record.strand == '+' ? '-' : '+';
    // flip the cigar
    Cigar cigar;
    gafkluge::for_each_cg(gaf_record, [&](const char& c, const int64_t& s) {
            cigar.push_back(make_pair(c, s));
        });
    cigar.reverse();
    assert(!cigar.empty());
    stringstream flipped_cigar;
    for (const auto& cig : cigar) {
        flipped_cigar << cig.second << cig.first;
    }
    gaf_record.opt_fields["cg"].second = flipped_cigar.str();
    // flip the path
    std::reverse(gaf_record.path.begin(), gaf_record.path.end());
    // flip the path offsets. do to this, we first measure its total (target) length
    int64_t path_target_len = 0;
    for (auto& step : gaf_record.path) {
        step.is_reverse = !step.is_reverse;
        //assert(step.is_stable);
        int64_t step_len = -1;
        // if the step is just a chromosome, we shimmy it into an interval so we treat consistently
        if (!step.is_interval) {
            if (!len_map.count(step.name)) {
                cerr << "[gaf2paf] error: unable to find " << step.name << " in lengths map" << endl;
                exit(1);
            }            
            step_len = len_map.at(step.name);
        } else {
            step_len = step.end - step.start;
        }           
        path_target_len += step_len;
    }
    int64_t rev_start = path_target_len - gaf_record.path_end;
    int64_t rev_end = path_target_len - gaf_record.path_start;
    gaf_record.path_start = rev_start;
    gaf_record.path_end = rev_end;
}

/* split a GAF record (forward strand: flip_gaf a '-' one first) into one PAF line per path step,
   calling fn on each line gaf2paf prints, with the line's cigar (in target order).  The PafLine
   carries columns 1-12 only: gaf2paf adds the parent's tags (see print_paf_line) */
inline void for_each_paf_line(const gafkluge::GafRecord& gaf_record, const unordered_map<string, int64_t>& len_map,
                              function<void(const PafLine&, const string&)> fn) {

    assert(gaf_record.strand == '+');
    
    // load up the cg_cigar
    Cigar cigar;
    gafkluge::for_each_cg(gaf_record, [&](const char& c, const int64_t& s) {
            cigar.push_back(make_pair(c, s));
        });

    // make a template output paf record
    PafLine paf_record;
    paf_record.query_name = gaf_record.query_name;
    paf_record.query_len = gaf_record.query_length;
    paf_record.mapq = gaf_record.mapq;

    int64_t path_len = gaf_record.path_end - gaf_record.path_start;
    Cigar::iterator cigar_pos = cigar.begin();
    
    int64_t query_base_count = 0; // keep track of bases in query
    int64_t target_base_count = 0; // and target
 
    // for every GAF step
    for (int64_t step_idx = 0; step_idx < gaf_record.path.size(); ++ step_idx) {
        auto step = gaf_record.path[step_idx];
        
        //assert(step.is_stable);
        //  get the length from the lookup
        if (!len_map.count(step.name)) {
            cerr << "[gaf2paf] error: unable to find " << step.name << " in lengths map" << endl;
            exit(1);
        }
        paf_record.target_name = step.name;
        paf_record.target_len = len_map.at(step.name);

        // if the step is just a chromosome, we shimmy it into an interval so we treat consistently
        if (!step.is_interval) {
            step.start = 0;
            step.end = paf_record.target_len;
            //assert(gaf_record.path.size() == 1);
        }
        // the path offsets affect the first and last step. end_offset here is the distance cut from the end of the last step
        int64_t start_offset = step_idx == 0 ? gaf_record.path_start : 0;
        int64_t end_offset = step_idx == gaf_record.path.size() - 1 ? target_base_count + (step.end - step.start) - path_len - start_offset: 0;
        assert(start_offset >= 0 && end_offset >= 0);

        // gobble up the step's worth of target bases from the cigar
        // we use target, because that's the only measure we have -- it's embedded in the step
        auto cig_range = cigar_next_by_target(cigar, cigar_pos, (step.end - end_offset) - (step.start + start_offset));

        if (step.is_reverse) {
            std::swap(start_offset, end_offset);
            std::reverse(cig_range.first, cig_range.second);
            paf_record.strand = '-';
        } else {
            paf_record.strand = '+';
        }

        // turn the cigar into a string
        stringstream cig_string;        
        int64_t cig_query_bases = 0;
        int64_t cig_target_bases = 0;
        paf_record.num_matching = 0;
        paf_record.num_bases  = 0;
        // todo: may need to reverse it!!
        for (auto i = cig_range.first; i != cig_range.second; ++i) {
            if (consumes_query(*i)) {
                cig_query_bases += i->second;
            }
            if (consumes_target(*i)) {
                cig_target_bases += i->second;
            }
            if (i->first == 'M' || i->first == '=') {
                paf_record.num_matching += i->second;
            }
            paf_record.num_bases += i->second;
            cig_string << i->second << i->first;
        }

        // make a new paf record
        paf_record.query_start = gaf_record.query_start + query_base_count;
        paf_record.query_end = paf_record.query_start + cig_query_bases;
        paf_record.target_start = step.start + start_offset;
        paf_record.target_end = step.end - end_offset;
        assert((step.end - end_offset) - (step.start + start_offset) == cig_target_bases);
        assert(paf_record.target_end - paf_record.target_start == cig_target_bases);
        assert(paf_record.query_end - paf_record.query_start == cig_query_bases);

        // we can delete over some nodes.  this leaves obnoxious 0-length paf lines
        // that don't serve a purpose (and all have overlapping query coordinates which
        // can raise flags when debugging gaffilter)... so don't print them
        if (paf_record.num_matching > 0) {
            fn(paf_record, cig_string.str());
        }
        
        // advance our counters
        query_base_count += cig_query_bases;
        target_base_count += cig_target_bases;
        cigar_pos = cig_range.second;
    }
}

// the identity gaf2paf writes as gi:f: (col 10 / col 11 of the parent GAF record, to 3 places)
inline double parent_identity(const gafkluge::GafRecord& gaf_record) {
    double identity = 0;
    if (gaf_record.block_length > 0) {
        identity = (double)gaf_record.matches / (double)gaf_record.block_length;
        identity = std::floor(identity * 1000 + 0.5 )/1000;
    }
    return identity;
}

/* print one line as gaf2paf does: columns 1-12, the parent record's tp/rc, then gm/gl/gi, which
   describe the parent GAF record (for downstream filtering), then the line's own cigar */
inline void print_paf_line(ostream& os, const gafkluge::GafRecord& gaf_record, const PafLine& paf_record,
                           const string& cigar) {
    // output the record
    os << paf_record;

    // todo: are there other optional tags we want to preserve? most would need to be recomputed to be
    // valid on alignment subregion
    if (gaf_record.opt_fields.count("tp")) {
        const auto& tp = gaf_record.opt_fields.at("tp");
        os << "\ttp:" << tp.first << ":" << tp.second;
    }

    if (gaf_record.opt_fields.count("rc")) {
        const auto& rc = gaf_record.opt_fields.at("rc");
        os << "\trc:" << rc.first << ":" << rc.second;
    }

    // throw in some information about the parent gaf (for downstream filtering)
    // gm: number of matches in the gaf
    os << "\tgm:i:" << gaf_record.matches;
    // gl: block length in the gaf
    os << "\tgl:i:" << gaf_record.block_length;
    // gi: the identity (col 10 / 11) in the gaf
    os << "\tgi:f:" << parent_identity(gaf_record);

    // output the cigar last
    os << "\tcg:Z:" << cigar << "\n";
}

typedef pair<int64_t, int64_t> KeptInterval;

/* parse the kq:Z: tag gaffilter -x writes on a record it cut: the query intervals the record keeps,
   0-based half-open, as s1-e1,s2-e2,...  Returns false if the tag is malformed */
inline bool parse_kept_intervals(const string& tag, vector<KeptInterval>& out) {
    out.clear();
    size_t pos = 0;
    while (pos < tag.length()) {
        size_t comma = tag.find(',', pos);
        string tok = tag.substr(pos, comma == string::npos ? string::npos : comma - pos);
        size_t dash = tok.find('-');
        if (dash == string::npos || dash == 0 || dash + 1 >= tok.length()) {
            return false;
        }
        try {
            size_t used1 = 0, used2 = 0;
            int64_t s = stoll(tok.substr(0, dash), &used1);
            int64_t e = stoll(tok.substr(dash + 1), &used2);
            if (used1 != dash || used2 != tok.length() - dash - 1 || e < s ||
                (!out.empty() && s < out.back().second)) {
                return false;
            }
            out.push_back(make_pair(s, e));
        } catch (...) {
            return false;
        }
        if (comma == string::npos) {
            break;
        }
        pos = comma + 1;
    }
    return !out.empty();
}

/* cut one PAF line (as for_each_paf_line makes it, with its cigar in target order) down to the
   query intervals keep (sorted, disjoint), calling fn once per piece that survives with the piece's
   columns 1-12 and cigar.  A line that lies inside one kept interval is passed through untouched.
   Otherwise every op is clipped to each interval: a deletion is kept only strictly inside it, the
   indels a cut leaves dangling at either end are stripped, and a piece with no '='/'M' column left
   makes no line.  On a '-' line the query runs backwards along the cigar */
inline void cut_paf_line(const PafLine& line, const string& cigar, const vector<KeptInterval>& keep,
                         function<void(const PafLine&, const string&)> fn) {
    int64_t qs = line.query_start, qe = line.query_end;
    vector<KeptInterval> ivs;
    for (const auto& k : keep) {
        if (k.second > qs && k.first < qe) {
            ivs.push_back(make_pair(std::max(k.first, qs), std::min(k.second, qe)));
        }
    }
    if (ivs.empty()) {
        return;
    }
    if (ivs.size() == 1 && ivs[0].first == qs && ivs[0].second == qe) {
        fn(line, cigar);
        return;
    }
    vector<pair<int64_t, char>> ops;
    for (size_t co = 0; co < cigar.length();) {
        size_t next = cigar.find_first_of("MIDNSHPX=", co);
        assert(next != string::npos);
        ops.push_back(make_pair((int64_t)stoll(cigar.substr(co, next - co)), cigar[next]));
        co = next + 1;
    }
    bool rev = line.strand == '-';
    struct Sel {
        int64_t n; char op; int64_t qlo, qhi, tlo, thi;
    };
    for (const auto& iv : ivs) {
        int64_t a = iv.first, b = iv.second;
        int64_t q = rev ? qe : qs;
        int64_t tp = line.target_start;
        vector<Sel> sel;
        for (const auto& o : ops) {
            int64_t n = o.first;
            char op = o.second;
            if (op == '=' || op == 'X' || op == 'M' || op == 'I') {
                int64_t lo, hi;
                if (!rev) {
                    lo = std::max(q, a);
                    hi = std::min(q + n, b);
                } else {
                    lo = std::max(q - n, a);
                    hi = std::min(q, b);
                }
                if (hi > lo) {
                    if (op == 'I') {
                        sel.push_back({hi - lo, 'I', lo, hi, tp, tp});
                    } else if (!rev) {
                        sel.push_back({hi - lo, op, lo, hi, tp + (lo - q), tp + (hi - q)});
                    } else {
                        sel.push_back({hi - lo, op, lo, hi, tp + (q - hi), tp + (q - lo)});
                    }
                }
                q += rev ? -n : n;
                if (op != 'I') {
                    tp += n;
                }
            } else {
                // a deletion is target only, between query bases q-1 and q
                if (a < q && q < b) {
                    sel.push_back({n, 'D', q, q, tp, tp + n});
                }
                tp += n;
            }
        }
        size_t first = 0, last = sel.size();
        while (first < last && (sel[first].op == 'I' || sel[first].op == 'D')) ++first;
        while (last > first && (sel[last - 1].op == 'I' || sel[last - 1].op == 'D')) --last;
        bool any_match = false;
        for (size_t i = first; i < last; ++i) {
            if (sel[i].op == '=' || sel[i].op == 'M') {
                any_match = true;
            }
        }
        if (!any_match) {
            continue;
        }
        PafLine piece = line;
        bool have_q = false, have_t = false;
        vector<pair<int64_t, char>> merged;
        for (size_t i = first; i < last; ++i) {
            const Sel& x = sel[i];
            if (x.op != 'D') {
                piece.query_start = have_q ? std::min(piece.query_start, x.qlo) : x.qlo;
                piece.query_end = have_q ? std::max(piece.query_end, x.qhi) : x.qhi;
                have_q = true;
            }
            if (x.op != 'I') {
                piece.target_start = have_t ? std::min(piece.target_start, x.tlo) : x.tlo;
                piece.target_end = have_t ? std::max(piece.target_end, x.thi) : x.thi;
                have_t = true;
            }
            if (!merged.empty() && merged.back().second == x.op) {
                merged.back().first += x.n;
            } else {
                merged.push_back(make_pair(x.n, x.op));
            }
        }
        piece.num_matching = 0;
        piece.num_bases = 0;
        stringstream cg;
        for (const auto& m : merged) {
            if (m.second == '=' || m.second == 'M') {
                piece.num_matching += m.first;
            }
            piece.num_bases += m.first;
            cg << m.first << m.second;
        }
        fn(piece, cg.str());
    }
}

} // namespace gaf2paf_split
