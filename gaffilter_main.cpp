/*
  Filter GAF records according to query overlap
 */

#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <list>
#include <cassert>

#include "gafkluge.hpp"
#include "paf.hpp"
#include "IntervalTree.h"
#include "pafcoverage.hpp"

//#define debug

using namespace std;
using namespace gafkluge;

typedef IntervalTree<int64_t, const GafRecord*> GafIntervalTree;
typedef GafIntervalTree::interval GafInterval;

// TODO: we have two ways of filtering below, dominates and dominates_mzgaf()
// need to figure out which one is better (all signs point to the old one, which,
// if I remember, is Gospel From Heng.  If that's the case we can just forget the new one
// if not here, then at least hide the options in cactus config.

// test if one record "dominates" another, using primary/secondary, mapq, block length in that order
static bool dominates(const GafRecord& gaf1, const GafRecord& gaf2, double ratio) {
    bool primary1 = !gaf1.opt_fields.count("tp") || gaf1.opt_fields.at("tp").second == "P";
    bool primary2 = !gaf2.opt_fields.count("tp") || gaf2.opt_fields.at("tp").second == "P";
    // empty interval can't dominate
    if (gaf1.query_start >= gaf1.query_end) {
        return false;
    } else if (gaf2.query_start >= gaf2.query_end) {
        return true;
    }
    if (primary1 && !primary2) {
        return true;
    } else if (primary2 && !primary1) {
        return false;
    }
    if ((double)gaf1.mapq / ((double)gaf2.mapq + 0.000001) >= ratio) {
        return true;
    } else if ((double)gaf2.mapq / ((double)gaf1.mapq + 0.000001) >= ratio) {
        return false;
    }
    if ((double)gaf1.block_length / ((double)gaf2.block_length + 0.000001) >= ratio) {
        return true;
    } else if ((double)gaf2.block_length / ((double)gaf1.block_length + 0.000001) >= ratio) {
        return false;
    }    
    return false;
}

// test if one record "dominates" another, using same logic as mzgaf2paf 
static bool dominates_mzgaf2paf(const GafRecord& gaf1, const GafRecord& gaf2, int64_t query_overlap_threshold) {
    return (gaf1.block_length >= query_overlap_threshold && gaf2.block_length < query_overlap_threshold ||
            gaf1.block_length < query_overlap_threshold && gaf2.block_length < query_overlap_threshold);
}

// measure the overlap
static int64_t overlap_size(const GafRecord& gaf1, const GafRecord& gaf2) {
    int64_t ostart = std::max(gaf1.query_start, gaf2.query_start);
    int64_t oend = std::min(gaf1.query_end, gaf2.query_end);
    assert(oend >= ostart);
    return oend - ostart;
}


// ---- subtractive ("trim") mode -------------------------------------------------------------
// The default action on losing an overlap is to delete the whole record, which also throws away
// the part of it nothing ever contested.  In trim mode the record instead gives up only the
// contested span.  Everything below exists to make that cut exactly, or refuse to make it.

typedef vector<pair<char, int64_t>> CigarVec;
typedef pair<int64_t, int64_t> QueryInterval;

static inline bool cig_query(char c) { return c == '=' || c == 'X' || c == 'M' || c == 'I'; }
static inline bool cig_path(char c)  { return c == '=' || c == 'X' || c == 'M' || c == 'D'; }
static inline bool cig_aligned(char c) { return c == '=' || c == 'X' || c == 'M'; }

static vector<QueryInterval> merge_intervals(vector<QueryInterval> ivs) {
    std::sort(ivs.begin(), ivs.end());
    vector<QueryInterval> out;
    for (const auto& iv : ivs) {
        if (!out.empty() && iv.first <= out.back().second) {
            out.back().second = std::max(out.back().second, iv.second);
        } else {
            out.push_back(iv);
        }
    }
    return out;
}

// the parts of [start, end) that no interval of cut covers.  cut must be merged and sorted
static vector<QueryInterval> subtract_intervals(int64_t start, int64_t end, const vector<QueryInterval>& cut) {
    vector<QueryInterval> out;
    int64_t cur = start;
    for (const auto& c : cut) {
        if (c.second <= cur) continue;
        if (c.first >= end) break;
        if (c.first > cur) out.push_back(make_pair(cur, std::min(c.first, end)));
        cur = std::max(cur, c.second);
        if (cur >= end) break;
    }
    if (cur < end) out.push_back(make_pair(cur, end));
    return out;
}

// every interval of ivs, with cut removed from it
static vector<QueryInterval> remove_interval(const vector<QueryInterval>& ivs, const QueryInterval& cut) {
    vector<QueryInterval> out;
    for (const auto& iv : ivs) {
        if (cut.second <= iv.first || cut.first >= iv.second) {
            out.push_back(iv);
            continue;
        }
        if (iv.first < cut.first) out.push_back(make_pair(iv.first, cut.first));
        if (cut.second < iv.second) out.push_back(make_pair(cut.second, iv.second));
    }
    return out;
}

// drop the first n query bases from cig, returning the path bases that went with them
static int64_t cigar_cut_front(CigarVec& cig, int64_t n) {
    int64_t path_used = 0;
    size_t i = 0;
    while (i < cig.size() && n > 0) {
        char c = cig[i].first;
        int64_t len = cig[i].second;
        if (!cig_query(c)) {
            // a path-only op inside the removed prefix goes with it
            path_used += len;
            ++i;
            continue;
        }
        int64_t take = std::min(n, len);
        n -= take;
        if (cig_path(c)) path_used += take;
        if (take == len) {
            ++i;
        } else {
            cig[i].second = len - take;
            break;
        }
    }
    cig.erase(cig.begin(), cig.begin() + i);
    return path_used;
}

// drop the last n query bases from cig, returning the path bases that went with them
static int64_t cigar_cut_back(CigarVec& cig, int64_t n) {
    int64_t path_used = 0;
    while (!cig.empty() && n > 0) {
        char c = cig.back().first;
        int64_t len = cig.back().second;
        if (!cig_query(c)) {
            path_used += len;
            cig.pop_back();
            continue;
        }
        int64_t take = std::min(n, len);
        n -= take;
        if (cig_path(c)) path_used += take;
        if (take == len) {
            cig.pop_back();
        } else {
            cig.back().second = len - take;
            break;
        }
    }
    return path_used;
}

// An alignment has to begin and end on an aligned column, so a leading or trailing indel left by
// the cut is dangling and comes off too.  Which end of the query that moves depends on the strand,
// because the cigar always runs along the path while a '-' record's query runs the other way.
static void cigar_strip_front(CigarVec& cig, GafRecord& rec) {
    while (!cig.empty() && !cig_aligned(cig.front().first)) {
        char c = cig.front().first;
        int64_t len = cig.front().second;
        if (cig_path(c)) rec.path_start += len;
        if (cig_query(c)) {
            if (rec.strand == '-') rec.query_end -= len;
            else rec.query_start += len;
        }
        cig.erase(cig.begin());
    }
}

static void cigar_strip_back(CigarVec& cig, GafRecord& rec) {
    while (!cig.empty() && !cig_aligned(cig.back().first)) {
        char c = cig.back().first;
        int64_t len = cig.back().second;
        if (cig_path(c)) rec.path_end -= len;
        if (cig_query(c)) {
            if (rec.strand == '-') rec.query_start += len;
            else rec.query_end -= len;
        }
        cig.pop_back();
    }
}

// Rebase the path onto the steps the cut actually left.  minigraph never emits a path whose
// leading nodes the alignment does not enter, and gaf2paf relies on that: it reads path_start as
// an offset into the FIRST step and path_end as one into the last.  Advancing those offsets past
// whole steps without dropping the steps produces a record that is arithmetically fine and that
// gaf2paf then misreads.  An unstable GAF names bare nodes, whose lengths only the -l map knows;
// without it the record cannot be rebased and is left for the caller to delete instead.
static bool rebase_path(GafRecord& rec, const unordered_map<string, int64_t>& node_lengths) {
    vector<int64_t> len(rec.path.size(), -1);
    int64_t known = 0;
    int64_t unknown_count = 0;
    size_t unknown_idx = 0;
    for (size_t i = 0; i < rec.path.size(); ++i) {
        const GafStep& step = rec.path[i];
        if (step.is_interval) {
            len[i] = step.end - step.start;
        } else {
            unordered_map<string, int64_t>::const_iterator it = node_lengths.find(step.name);
            if (it != node_lengths.end()) {
                len[i] = it->second;
            }
        }
        if (len[i] < 0) {
            ++unknown_count;
            unknown_idx = i;
        } else {
            known += len[i];
        }
    }
    // A step that names a whole stable sequence carries no interval, and on a stable GAF there is
    // no lengths file to ask.  The path length column is the sum over all steps, so on a ONE-STEP
    // path the missing length is pinned by arithmetic -- which is the common case (10,532 of
    // 25,763 records on HPRC chr20) and is safe because a one-step path has nothing to drop.
    // The inference is refused on a multi-step path: there it would satisfy the total check below
    // by construction, disabling the only cross-check this function has, and a wrong length would
    // then drop the wrong steps and leave the record naming a node it does not align to.
    if (unknown_count == 1 && rec.path.size() == 1) {
        len[unknown_idx] = rec.path_length - known;
        if (len[unknown_idx] < 0) {
            return false;
        }
        known += len[unknown_idx];
    } else if (unknown_count > 0) {
        return false;
    }
    int64_t total = known;
    if (total != rec.path_length) {
        return false;                     // our model of the path disagrees with the record
    }
    size_t first = 0;
    int64_t head = 0;
    while (first < rec.path.size() && head + len[first] <= rec.path_start) {
        head += len[first];
        ++first;
    }
    size_t last = rec.path.size();
    int64_t tail = total;
    while (last > first && tail - len[last - 1] >= rec.path_end) {
        tail -= len[last - 1];
        --last;
    }
    if (last <= first) {
        return false;
    }
    rec.path = vector<GafStep>(rec.path.begin() + first, rec.path.begin() + last);
    rec.path_start -= head;
    rec.path_end -= head;
    rec.path_length = tail - head;
    return true;
}

// Cut rec down to the query sub-interval [new_qs, new_qe), rewriting its coordinates and cigar.
// Returns false, with rec left unusable, if the cigar cannot be reconciled with the columns it is
// supposed to describe or if nothing survives -- the caller then deletes the record whole, which
// is what would have happened anyway.  Better a record lost than a record that lies about where
// it aligns.
static bool trim_gaf_record(GafRecord& rec, int64_t new_qs, int64_t new_qe,
                            const unordered_map<string, int64_t>& node_lengths) {
    if (new_qe <= new_qs || !rec.opt_fields.count("cg")) {
        return false;
    }
    CigarVec cig;
    for_each_cg(rec, [&](const char& c, const size_t& l) {
            cig.push_back(make_pair(c, (int64_t)l));
        });
    if (cig.empty()) {
        return false;
    }
    // only safe to cut a cigar that agrees with the record it belongs to.  this also rejects M
    // cigars, where '=' cannot be counted and so matches cannot be recomputed after the cut
    int64_t q = 0, p = 0, m = 0, bl = 0;
    for (const auto& e : cig) {
        if (cig_query(e.first)) q += e.second;
        if (cig_path(e.first)) p += e.second;
        if (e.first == '=') m += e.second;
        bl += e.second;
    }
    if (q != rec.query_end - rec.query_start || p != rec.path_end - rec.path_start ||
        m != rec.matches || bl != rec.block_length) {
        return false;
    }
    if (new_qs <= rec.query_start && new_qe >= rec.query_end) {
        return true;                          // nothing to cut
    }

    // the cigar runs along the path, so for a '-' record the low query end is at its BACK
    bool rev = rec.strand == '-';
    int64_t lo_cut = std::max((int64_t)0, new_qs - rec.query_start);
    int64_t hi_cut = std::max((int64_t)0, rec.query_end - new_qe);
    int64_t front_cut = rev ? hi_cut : lo_cut;
    int64_t back_cut = rev ? lo_cut : hi_cut;

    rec.query_start = new_qs;
    rec.query_end = new_qe;
    rec.path_start += cigar_cut_front(cig, front_cut);
    rec.path_end -= cigar_cut_back(cig, back_cut);
    if (front_cut) cigar_strip_front(cig, rec);
    if (back_cut) cigar_strip_back(cig, rec);
    if (cig.empty() || rec.query_end <= rec.query_start || rec.path_end <= rec.path_start) {
        return false;
    }

    // recompute everything the cut invalidated
    int64_t nm = 0;
    q = p = m = bl = 0;
    stringstream cg;
    for (const auto& e : cig) {
        if (cig_query(e.first)) q += e.second;
        if (cig_path(e.first)) p += e.second;
        if (e.first == '=') m += e.second;
        else nm += e.second;
        bl += e.second;
        cg << e.second << e.first;
    }
    if (q != rec.query_end - rec.query_start || p != rec.path_end - rec.path_start) {
        return false;
    }
    rec.matches = m;
    rec.block_length = bl;
    rec.opt_fields["cg"] = make_pair("Z", cg.str());
    rec.opt_fields["NM"] = make_pair("i", std::to_string(nm));
    if (rec.opt_fields.count("gi")) {
        stringstream ss;
        ss << (double)m / (double)bl;
        rec.opt_fields["gi"] = make_pair("f", ss.str());
    }
    // these describe the untrimmed alignment and nothing here can cut them.  dv is minigraph's
    // own chain-derived divergence estimate, not (block_length - matches)/block_length, so it
    // cannot be recomputed from the cigar either -- dropping it beats redefining it in place.
    rec.opt_fields.erase("cs");
    rec.opt_fields.erase("ds");
    rec.opt_fields.erase("dv");
    return rebase_path(rec, node_lengths);
}

// The order dominates() applies, as a total order, for the -R fallback: primary before
// secondary, then MAPQ, then block length, then input position so the result cannot depend on it.
static bool record_precedes(const GafRecord& a, const GafRecord& b) {
    bool pa = !a.opt_fields.count("tp") || a.opt_fields.at("tp").second == "P";
    bool pb = !b.opt_fields.count("tp") || b.opt_fields.at("tp").second == "P";
    if (pa != pb) return pa;
    if (a.mapq != b.mapq) return a.mapq > b.mapq;
    return a.block_length > b.block_length;
}

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <gaf> > output.gaf" << endl
         << "Filter GAF record if its query interval overlaps another query interval and\n"
         << "  1) the record is secondary and the overlapping record is primary or\n"
         << "  2) the record's MAPQ is lower than {ratio, see -r} times the overlapping record's MAPQ or\n"
         << "  3) the record's block length is less than {ratio, see -r} times larger than the overlapping record's block length (and its MAPQ isn't higher)" << endl
         << "  Also: the -o option can be used to mimic mzgaf2paf's query overlap filter" << endl
         << endl
         << "options: " << endl
         << "    -r, --ratio N                   If two query blocks overlap, and one is Nx bigger than the other, the bigger one is kept (otherwise both deleted) [0]" << endl
         << "    -m, --min-overlap N             Ignore overlaps that consitute <N% of the length [0]" << endl
         << "    -o, --min-overlap-length N      If >= 2 query regions with size >= N overlap, ignore the query region.  If 1 query region with size >= N overlaps any regions of size <= N, ignore the smaller ones only. Works separate to -r/-m but can be used in conjunction with them to combine the two filters (0 = disable) [0]" << endl
         << "    -q, --min-mapq N                Don't let an interval with MAPQ < N cause something to be filtered out" << endl
         << "    -b, --min-block-length N        Don't let an interval with block length < N cause something to be filtered out" << endl
         << "    -i, --min-identity N            Don't let an interval with identity < N cause something to be filtered out" << endl       
         << "    -t, --trim                      Instead of deleting a record that loses an overlap, cut the contested span out of it and keep the rest (GAF input only)" << endl
         << "    -g, --trim-min-gap N            With -t, only a hole longer than N, with alignment still on both sides, is worth closing at all (see -R); shorter ones do not split a path downstream [10000]" << endl
         << "    -l, --node-lengths FILE         Node lengths (as written by gaf2unstable -o). Needed by -t on an unstable GAF, whose path names carry no interval" << endl
         << "    -e, --trim-edge N               With -t, also cut N bases beyond each side of a contested span. The bases abutting an overlap are the least trustworthy part of the alignment, and cactus keeps unaligned stretches shorter than its own clip threshold anyway [5000]" << endl
         << "    -R, --rescue-weak               With -t, close such a hole by giving the span to the best claimant (primary, then MAPQ, then block length). Off by default, because no such choice can meet the bar -r sets: leaving a hole clips sequence out, but a wrong placement puts a wrong alignment in" << endl
         << "    -p, --paf                       Input is PAF, not GAF" << endl;
}    

int main(int argc, char** argv) {
    // the result goes to standard output, and a failed write there is
    // otherwise reported to nobody
    check_stdout_at_exit();


    double ratio = 0.;
    double min_overlap_pct = 0.;
    int64_t min_overlap_len = 0;
    int64_t min_block_len = 0;
    int64_t min_mapq = 0;
    double min_identity = 0;
    bool trim_mode = false;
    bool rescue_weak = false;
    int64_t trim_min_gap = 10000;
    int64_t trim_edge = 5000;
    string node_lengths_path;

    int c;
    bool is_paf = false;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"ratio", required_argument, 0, 'r'},
            {"min-overlap", required_argument, 0, 'm'},
            {"min-overlap-length", required_argument, 0, 'o'},
            {"min-block-length", required_argument, 0, 'b'},
            {"min-mapq", required_argument, 0, 'q'},
            {"min-identity", required_argument, 0, 'i'},
            {"trim", no_argument, 0, 't'},
            {"trim-min-gap", required_argument, 0, 'g'},
            {"rescue-weak", no_argument, 0, 'R'},
            {"trim-edge", required_argument, 0, 'e'},
            {"node-lengths", required_argument, 0, 'l'},
            {"paf", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:r:m:po:b:q:i:tg:l:Re:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            ratio = stof(optarg);
            break;
        case 'm':
            min_overlap_pct = stof(optarg);
            break;
        case 'o':
            min_overlap_len = std::stol(optarg);
            break;            
        case 'p':
            is_paf = true;
            break;
        case 'b':
            min_block_len = std::stol(optarg);
            break;
        case 'i':
            min_identity = std::stof(optarg);
            break;
        case 'q':
            min_mapq = std::stol(optarg);
            break;
        case 't':
            trim_mode = true;
            break;
        case 'g':
            trim_min_gap = std::stol(optarg);
            break;
        case 'l':
            node_lengths_path = optarg;
            break;
        case 'R':
            rescue_weak = true;
            break;
        case 'e':
            trim_edge = std::stol(optarg);
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }

    if (trim_mode && is_paf) {
        cerr << "[gaffilter] error: -t/--trim needs the cigar and coordinates of a GAF record, "
             << "and cannot be used with -p" << endl;
        return 1;
    }

    if ((ratio == 0) && (min_overlap_len == 0)) {
        cerr << "[gaffilter] error: at least one of -r or -o must be used to specify filter" << endl;
        return 1;
    }

    // Parse the positional argument
    if (optind >= argc) {
        cerr << "[gaffilter] error: too few arguments" << endl;
        help(argv);
        return 1;
    }
    
    string gaf_path = argv[optind++];
    
    // open the gaf file
    ifstream in_file;
    istream* in_stream;
    if (gaf_path == "-") {
        in_stream = &cin;
    } else {
        in_file.open(gaf_path);
        if (!in_file) {
            cerr << "[gaffilter] error: unable to open input: " << gaf_path << endl;
            return 1;
        }
        in_stream = &in_file;
    }

    // -t needs to know how long each path step is in order to drop the ones the cut left behind.
    // Read AFTER the input below, not here: the usual caller is
    //   gaf2unstable ... -o lengths.tsv | gaffilter - -l lengths.tsv
    // and both sides of that pipe start at once, so the file does not exist yet.  gaf2unstable
    // writes it in full before it emits its first GAF line, so by the time this process has read
    // its input to EOF the file is complete.
    unordered_map<string, int64_t> node_lengths;

    // shimmy in paf support post hoc (at the cost of storing a dummy gaf record list in memory!)
    vector<PafLine> paf_records;
    function<string(const GafRecord&)> print_record = [&](const GafRecord& gaf_record) {
        stringstream ss;
        if (is_paf) {
            // hack alert: hijack path_length with offset in paf_records
            ss << paf_records[gaf_record.path_length];
        } else {
            ss << gaf_record;
        }
        return ss.str();
    };

    // just load the gaf into memory
    vector<GafRecord> gaf_records;
    string line_buffer;
    while (getline(*in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            // skip -S stuff
            continue;
        }        
        GafRecord gaf_record;
        if (is_paf) {
            PafLine paf_record = parse_paf_line(line_buffer);
            paf_records.push_back(paf_record);
            // just copy what we (might) need
            gaf_record.query_name = paf_record.query_name;
            gaf_record.query_length = paf_record.query_len;
            gaf_record.query_start = paf_record.query_start;
            gaf_record.query_end = paf_record.query_end;
            gaf_record.strand = paf_record.strand;
            gaf_record.mapq = paf_record.mapq;
            if (paf_record.opt_fields.count("gl")) {
                gaf_record.block_length = stol(paf_record.opt_fields.at("gl").second);
            } else {
                gaf_record.block_length = paf_record.num_bases;
            }
            if (paf_record.opt_fields.count("gm")) {
                gaf_record.matches = stol(paf_record.opt_fields.at("gm").second);
            } else {
                gaf_record.matches = paf_record.num_matching;
            }
            if (paf_record.opt_fields.count("tp")) {
                gaf_record.opt_fields["tp"] = paf_record.opt_fields.at("tp");
            }
            if (paf_record.opt_fields.count("rc")) {
                gaf_record.opt_fields["rc"] = paf_record.opt_fields.at("rc");
            }
            // hack alert: hijack path_length with offset in paf_records
            gaf_record.path_length = paf_records.size() - 1;
        } else {
            parse_gaf_record(line_buffer, gaf_record);
        }
        gaf_records.push_back(gaf_record);
    }
    cerr << "[gaffilter]: Loaded " << gaf_records.size() << (is_paf ? " PAF" : " GAF") << " records" << endl;

    // now that the input is exhausted, the upstream writer of the lengths file has finished
    if (!node_lengths_path.empty()) {
        ifstream len_file(node_lengths_path);
        if (!len_file) {
            cerr << "[gaffilter] error: unable to open node lengths: " << node_lengths_path << endl;
            return 1;
        }
        // line-oriented and first-two-fields, because a .fai has five columns and is exactly
        // what the sibling tools document for -l.  Reading with >> would silently take columns
        // 3 and 4 of line 1 as the next name/length pair and then stop, leaving a map that is
        // not empty but is garbage -- which disabled the trim without saying so.
        string line;
        int64_t nline = 0;
        while (getline(len_file, line)) {
            ++nline;
            if (line.empty() || line[0] == '#') {
                continue;
            }
            stringstream ss(line);
            string name, len_tok;
            if (!(ss >> name >> len_tok)) {
                cerr << "[gaffilter] error: " << node_lengths_path << ":" << nline
                     << ": expected a name and a length" << endl;
                return 1;
            }
            try {
                node_lengths[name] = stol(len_tok);
            } catch (...) {
                cerr << "[gaffilter] error: " << node_lengths_path << ":" << nline
                     << ": second field is not a length: " << len_tok << endl;
                return 1;
            }
        }
        cerr << "[gaffilter]: Loaded " << node_lengths.size() << " node lengths" << endl;
    }

    // make an interval tree for each query sequence
    unordered_map<string, vector<GafInterval>> gaf_intervals;
    for (const auto& gaf_record : gaf_records) {
        int64_t end_point = gaf_record.query_end;
        if (end_point > gaf_record.query_start) {
            // interval tree expects closed coordinates.  but it also expects the end point >= start point
            --end_point;
        }
        GafInterval gaf_interval(gaf_record.query_start, gaf_record.query_end - 1, &gaf_record);
        gaf_intervals[gaf_record.query_name].push_back(gaf_interval);
    }
    unordered_map<string, GafIntervalTree*> gaf_trees;
    for (const auto& qi : gaf_intervals) {
        gaf_trees[qi.first] = new GafIntervalTree(qi.second);
    }
    gaf_intervals.clear();
    cerr << "[gaffilter]: Constructed interval trees" << endl;


    int64_t filter_count = 0;
    int64_t filter_len_count = 0;
    int64_t trim_count = 0;
    int64_t trim_len_count = 0;
    int64_t rescue_count = 0;
    int64_t rescue_declined = 0;
    int64_t holes_opened = 0;
    int64_t holes_opened_bp = 0;
    int64_t trim_fail_count = 0;

    // in trim mode a record that loses an overlap is not deleted: it gives up only the span it
    // lost, collected here.  without -t contested stays empty and the original keep/drop applies.
    vector<char> keep(gaf_records.size(), 1);
    // sized only under -t: an empty vector per record is 24 bytes, which on a production-scale
    // PAF is over a hundred megabytes paid by runs that never asked to trim
    vector<vector<QueryInterval> > contested(trim_mode ? gaf_records.size() : 0);

    // simple algorithm:
    // for each record, scan its overlaps and flag it if it finds anything
    // that overlaps that isn't ratio X smaller.
    // this is a really inefficient in worst-case (where everything overlaps) but that's not at all what we expect
    for (int64_t i = 0; i < gaf_records.size(); ++i) {
        int64_t end_point = gaf_records[i].query_end;
        if (end_point > gaf_records[i].query_start) {
            // interval tree expects closed coordinates.  but it also expects the end point >= start point
            --end_point;
        }
        vector<GafInterval> overlapping;
        string ref_contig;
        if (gaf_records[i].opt_fields.count("rc")) {
            ref_contig = gaf_records[i].opt_fields.at("rc").second;
        }
        gaf_trees[gaf_records[i].query_name]->visit_overlapping(gaf_records[i].query_start, end_point, [&](const GafInterval& interval) {
                // filter self alignments
                double identity = interval.value->matches ? (double)interval.value->block_length / (double)interval.value->matches  : 0;
                assert(identity >= 0);
                if (interval.value->opt_fields.count("gi")) {
                    identity = min(identity, (double)stof(interval.value->opt_fields.at("gi").second));
                }
                if (interval.value != &gaf_records[i] &&
                    //and mapq/block length failing alignments
                    interval.value->mapq >= min_mapq && (interval.value->query_length <= min_block_len || interval.value->block_length >= min_block_len) &&
                    // and identity failing alignments
                    identity >= min_identity) {
                    string overlap_contig;
                    if (interval.value->opt_fields.count("rc")) {
                        overlap_contig = interval.value->opt_fields.at("rc").second;
                    }
                    // also ignore things that map to different contigs
                    if (ref_contig == overlap_contig || ref_contig.empty() || overlap_contig.empty()) {
                        int64_t overlap_bases = overlap_size(gaf_records[i], *interval.value);
                        // filter overlaps that are too small to matter (via min_overlap_pct)
                        if (gaf_records[i].block_length == 0 ||
                            (double)overlap_bases / (double)gaf_records[i].block_length >= min_overlap_pct) {
                            overlapping.push_back(interval);
                        }
                    }
                }
            });
        for (const auto& ogi : overlapping) {
            bool is_dominant = true;
            if (ratio) {
                is_dominant = dominates(gaf_records[i], *ogi.value, ratio);
            }
            if (is_dominant && min_overlap_len) {
                is_dominant = dominates_mzgaf2paf(gaf_records[i], *ogi.value, min_overlap_len);
            }
            if (!is_dominant) {
                if (!trim_mode) {
                    keep[i] = 0;
                    break;
                }
                contested[i].push_back(make_pair(std::max(gaf_records[i].query_start, ogi.value->query_start),
                                                 std::min(gaf_records[i].query_end, ogi.value->query_end)));
            }
        }
        if (trim_mode) {
            contested[i] = merge_intervals(contested[i]);
            // Widen each contested span by -e before anything downstream looks at it, so the
            // rescue and the emit loop agree on what is actually given up.  The bases butting up
            // against an overlap are where the alignment is least certain, and cactus will carry
            // an unaligned stretch shorter than its clip threshold regardless.  Clipping to the
            // record's own span makes the outer side a no-op, so in practice only the interior
            // borders move; a record shorter than the widening is given up whole.
            if (trim_edge > 0 && !contested[i].empty()) {
                vector<QueryInterval> widened;
                for (const auto& c : contested[i]) {
                    widened.push_back(make_pair(
                        std::max(gaf_records[i].query_start, c.first - trim_edge),
                        std::min(gaf_records[i].query_end, c.second + trim_edge)));
                }
                contested[i] = merge_intervals(widened);
            }
        }
#ifdef debug
        if (!keep[i] || !contested[i].empty()) {
            cerr << "\nfiltering record " << i << " (" << &gaf_records[i] << ") because it doesn't dominate its "
                 << overlapping.size() << " overlaps\n  " << print_record(gaf_records[i]) << endl;
            int64_t ocount = 0;
            for (const auto& ogi : overlapping) {
                cerr << "overlap " << ocount++ << " (" << ogi.value << "):\n  " << print_record(*ogi.value) << endl;
            }
        }
#endif
    }

    // The filter may trim, but it must not perforate.  Each cut above is justified on its own, but
    // a contested span given up by every record that claims it leaves an unaligned hole with
    // alignment still standing on both sides -- a breakpoint the assembly does not have.  The -m
    // guard cannot see this: it judges one pair at a time, while a hole is a property of what the
    // whole query has left.  Where a hole would open, the span is handed whole to one claimant
    // instead, so the bases stay placed exactly once rather than zero times.  Ties go to the
    // largest record and then to the earliest, so the result does not depend on input order.
    // cactus-graphmap-join --clip only breaks a path on unaligned stretches longer than its
    // threshold, so -g leaves shorter holes alone.  Iterated, because closing one hole adds
    // coverage that can turn a neighbouring end trim into an interior one.
    if (trim_mode && trim_min_gap >= 0) {
        unordered_map<string, vector<int64_t> > by_query;
        for (int64_t i = 0; i < (int64_t)gaf_records.size(); ++i) {
            if (!contested[i].empty()) {
                by_query[gaf_records[i].query_name];
            }
        }
        for (int64_t i = 0; i < (int64_t)gaf_records.size(); ++i) {
            if (by_query.count(gaf_records[i].query_name)) {
                by_query[gaf_records[i].query_name].push_back(i);
            }
        }
        for (auto& q : by_query) {
            // coverage BEFORE any cut, so a gap can be told apart from one that was always there.
            // Only a gap inside this is something the trim opened and therefore ours to answer for.
            vector<QueryInterval> pre;
            for (int64_t i : q.second) {
                if (gaf_records[i].query_end > gaf_records[i].query_start) {
                    pre.push_back(make_pair(gaf_records[i].query_start, gaf_records[i].query_end));
                }
            }
            pre = merge_intervals(pre);
            for (int round = 0; round < 32; ++round) {
                vector<QueryInterval> cov;
                for (int64_t i : q.second) {
                    if (gaf_records[i].query_end <= gaf_records[i].query_start) {
                        continue;
                    }
                    for (const auto& f : subtract_intervals(gaf_records[i].query_start,
                                                            gaf_records[i].query_end, contested[i])) {
                        cov.push_back(f);
                    }
                }
                cov = merge_intervals(cov);
                bool changed = false;
                for (size_t k = 0; k + 1 < cov.size(); ++k) {
                    QueryInterval gap(cov[k].second, cov[k + 1].first);
                    if (gap.second - gap.first <= trim_min_gap) {
                        continue;
                    }
                    // A record may claim the gap if its own alignment spans it.  Testing instead
                    // that ONE contested interval contains the gap missed any hole formed by two
                    // adjacent seams contested by different records, and left it open while
                    // reporting nothing.  A record that spans an uncovered gap must have given all
                    // of it up, so spanning is the right and sufficient test.
                    vector<int64_t> claim;
                    for (int64_t i : q.second) {
                        if (!contested[i].empty() &&
                            gaf_records[i].query_start <= gap.first &&
                            gap.second <= gaf_records[i].query_end) {
                            claim.push_back(i);
                        }
                    }
                    bool was_covered = false;
                    for (const auto& pv : pre) {
                        if (pv.first <= gap.first && gap.second <= pv.second) {
                            was_covered = true;
                            break;
                        }
                    }
                    if (!was_covered) {
                        continue;             // a gap that was already there; not ours to close
                    }
                    if (claim.empty()) {
                        // the gap is the union of two adjacent seams contested by different
                        // records, so no single record spans it and none can close it alone
                        ++holes_opened;
                        holes_opened_bp += gap.second - gap.first;
                        continue;
                    }
                    // No claimant can ever be the one the filter would pick.  A hole exists only
                    // if EVERY record spanning it gave it up, and a record gives a span up only to
                    // a competitor it does not dominate -- a competitor which therefore spans the
                    // hole too and is itself a claimant.  So for any hole, no claimant dominates
                    // all the others; gating the rescue on dominates() is exactly "never rescue".
                    // (Measured before the argument was noticed: 0 rescues pass that gate on HPRC
                    // chr9, chr15 and chr20.)  Closing the hole therefore always means choosing on
                    // evidence the filter itself rejects, so it is opt-in.  A hole clips sequence
                    // out of the graph; a wrong placement puts a wrong alignment into it, and the
                    // second is the failure this tool exists to prevent.
                    if (!rescue_weak) {
                        ++rescue_declined;
                        holes_opened_bp += gap.second - gap.first;
                        continue;
                    }
                    // -R: order by the filter's own precedence, not raw block length, so a
                    // secondary or a lower-MAPQ record cannot take the span from a primary
                    int64_t best = -1;
                    for (int64_t i : claim) {
                        if (best < 0 || record_precedes(gaf_records[i], gaf_records[best])) {
                            best = i;
                        }
                    }
                    contested[best] = remove_interval(contested[best], gap);
                    ++rescue_count;
                    changed = true;
                }
                if (!changed) {
                    break;
                }
            }
        }
    }

    for (int64_t i = 0; i < (int64_t)gaf_records.size(); ++i) {
        // a survivor with nothing to give up goes out untouched.  this is decided on the verdict
        // and not on whether the record has any query span left, because a record with an empty
        // query interval yields no fragments and the old code still printed it
        if (keep[i] && (!trim_mode || contested[i].empty())) {
            cout << print_record(gaf_records[i]) << "\n";
            continue;
        }
        int64_t emitted = 0;
        bool cut_failed = false;
        if (keep[i]) {
            for (const auto& f : subtract_intervals(gaf_records[i].query_start,
                                                    gaf_records[i].query_end, contested[i])) {
                GafRecord cut = gaf_records[i];
                if (!trim_gaf_record(cut, f.first, f.second, node_lengths)) {
                    cut_failed = true;
                    continue;
                }
                cout << cut << "\n";
                ++emitted;
                trim_len_count += cut.block_length;
            }
        }
        if (cut_failed && !emitted) {
            ++trim_fail_count;
        }
        if (emitted) {
            ++trim_count;
        } else {
            ++filter_count;
            if (is_paf) {
                filter_len_count += paf_records[i].num_bases;
            } else {
                filter_len_count += gaf_records[i].block_length;
            }
        }
    }

    for (auto qt : gaf_trees) {
        delete qt.second;
    }
    
    cerr << "[gaffilter]: filtered " << filter_count << " / " << gaf_records.size() << ". total block lengths filtered: " << filter_len_count << endl;
    if (trim_mode) {
        cerr << "[gaffilter]: trimmed " << trim_count << " records, keeping " << trim_len_count
             << " block length that whole-record deletion would have dropped. rescued "
             << rescue_count << " contested spans that would have left a hole" << endl;
        if (rescue_declined || holes_opened) {
            cerr << "[gaffilter]: left " << (rescue_declined + holes_opened) << " hole(s), "
                 << holes_opened_bp << " bp, in coverage the input had: " << rescue_declined
                 << " ambiguous"
                 << (rescue_weak ? "" : " (-R would close them, on evidence this filter rejects)")
                 << ", " << holes_opened << " spanned by no single record" << endl;
        }
        if (trim_fail_count) {
            cerr << "[gaffilter]: warning: " << trim_fail_count << " record(s) could not be cut "
                 << "(cigar or path inconsistent with the record"
                 << (node_lengths.empty() ? ", or node lengths needed -- see -l" : "")
                 << ") and were deleted whole" << endl;
        }
    }
    return 0;
}
