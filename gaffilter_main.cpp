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
#include <map>
#include <set>
#include <tuple>
#include <cmath>
#include <unordered_set>

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

// which test of dominates() decided a contest (gaffilter -x needs to know, not only who won)
enum DomClause { DOM_EMPTY, DOM_PRIMARY, DOM_SECONDARY, DOM_MAPQ, DOM_LENGTH, DOM_TIE };

// dominates(), on the values it reads, also saying which test decided
static bool dominance(int64_t query_start1, int64_t query_end1, bool primary1, double mapq1, double block_length1,
                      int64_t query_start2, int64_t query_end2, bool primary2, double mapq2, double block_length2,
                      double ratio, DomClause& clause) {
    // empty interval can't dominate
    if (query_start1 >= query_end1) {
        clause = DOM_EMPTY;
        return false;
    } else if (query_start2 >= query_end2) {
        clause = DOM_EMPTY;
        return true;
    }
    if (primary1 && !primary2) {
        clause = DOM_PRIMARY;
        return true;
    } else if (primary2 && !primary1) {
        clause = DOM_SECONDARY;
        return false;
    }
    clause = DOM_MAPQ;
    if (mapq1 / (mapq2 + 0.000001) >= ratio) {
        return true;
    } else if (mapq2 / (mapq1 + 0.000001) >= ratio) {
        return false;
    }
    clause = DOM_LENGTH;
    if (block_length1 / (block_length2 + 0.000001) >= ratio) {
        return true;
    } else if (block_length2 / (block_length1 + 0.000001) >= ratio) {
        return false;
    }
    clause = DOM_TIE;
    return false;
}

static inline bool is_primary(const GafRecord& gaf) {
    return !gaf.opt_fields.count("tp") || gaf.opt_fields.at("tp").second == "P";
}

// test if one record "dominates" another, using primary/secondary, mapq, block length in that order
static bool dominates(const GafRecord& gaf1, const GafRecord& gaf2, double ratio) {
    DomClause clause;
    return dominance(gaf1.query_start, gaf1.query_end, is_primary(gaf1), (double)gaf1.mapq, (double)gaf1.block_length,
                     gaf2.query_start, gaf2.query_end, is_primary(gaf2), (double)gaf2.mapq, (double)gaf2.block_length,
                     ratio, clause);
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


// ---- intervals and paths shared with exact mode (-x) ---------------------------------------

typedef pair<int64_t, int64_t> QueryInterval;

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

// Rebase the path onto the steps its alignment actually enters.  minigraph never emits a path
// whose leading nodes the alignment does not enter, and gaf2paf relies on that: it reads
// path_start as an offset into the FIRST step and path_end as one into the last.  A record with
// offsets past whole steps is arithmetically fine but gaf2paf misreads it.  An unstable GAF names
// bare nodes, whose lengths only node_lengths knows; without them the record cannot be rebased and
// false is returned.
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

// make an interval tree for each query sequence
static unordered_map<string, GafIntervalTree*> build_query_trees(const vector<GafRecord>& gaf_records) {
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
    return gaf_trees;
}

// The overlap filter itself: keep[i] is cleared for a record that fails to dominate an overlap.
// -x runs this as the stock baseline it is measured against.
static void overlap_filter(const vector<GafRecord>& gaf_records, unordered_map<string, GafIntervalTree*>& gaf_trees,
                           double ratio, double min_overlap_pct, int64_t min_overlap_len, int64_t min_mapq,
                           int64_t min_block_len, double min_identity, vector<char>& keep,
                           function<string(const GafRecord&)>& print_record) {
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
                keep[i] = 0;
                break;
            }
        }
#ifdef debug
        if (!keep[i]) {
            cerr << "\nfiltering record " << i << " (" << &gaf_records[i] << ") because it doesn't dominate its "
                 << overlapping.size() << " overlaps\n  " << print_record(gaf_records[i]) << endl;
            int64_t ocount = 0;
            for (const auto& ogi : overlapping) {
                cerr << "overlap " << ocount++ << " (" << ogi.value << "):\n  " << print_record(*ogi.value) << endl;
            }
        }
#endif
    }

}

// ---- exact mode (-x) -------------------------------------------------------------------------
// The stock filter deletes a record whole when it loses an overlap.  -x instead resolves every
// elementary query segment on its own (it replaces the old -t, which cut the whole contested span,
// widened, out of every loser): a record keeps a
// segment if it beats (dominates(), with --exact-ratio) every compatible record that also claims it.
// A record that would have been deleted whole by the stock rule (it loses to a competitor overlapping
// -m of its block) is "demoted": it only fills segments no compatible non-demoted record claims, and
// what it keeps (its "remainders") must be long and similar enough.  Whatever the result anchors that
// the stock filter chain would not ("new" sequence) is gated against the reference the haplotype
// already covers, and the new structure it forms -- every off-backbone excursion along the contig --
// is tested: an excursion that cannot be a real two-sided rearrangement keeps its sequence, but its
// junctions to the backbone are broken by a query gap wide enough that cactus cannot rejoin them.
// The reference implementation is the overlap-pilot emulator (emu3j/emu3s); this is a port of it.
//
// Output: whole records as today; a cut record is printed whole, with the query intervals it keeps
// in a kq:Z: tag, which gaf2paf applies line by line.  The record's own block length and matches
// thus reach the PAF's gl/gm unchanged, which the block-length line filter downstream relies on.
namespace exact {

typedef pair<int64_t, int64_t> IV;
typedef vector<IV> IVs;

struct Params {
    // the stock stage (A0), exactly as cactus runs it: parsed as the stock options are (stof)
    double stock_ratio = 0;
    double stock_min_overlap = 0;
    double stock_min_identity = 0;
    int64_t min_mapq = 0;
    int64_t min_block = 0;
    // the arm's own thresholds, parsed as doubles (the emulator's are)
    double stock_ratio_d = 0;       // -r, for the guard's "stock would tie" test
    double ratio = 0;               // --exact-ratio (default -r)
    double min_overlap = 0;         // -m as the demotion threshold
    double min_identity = 0;        // -i as eligibility (matches/block)
    // A0's -p stage, as filter_paf runs it (0 = off)
    double paf_ratio = 0;
    double paf_min_overlap = 0;
    string nodes_path;
    string ref_event;
    string a0_paf_path;
    string only_chrom;
    string plan_path;
    string log_path;
    string junctions_path;
    string summary_path;
    double rem_floor = 0.95;
    int64_t rmin = 10000;
    double id_floor = 0.98;
    bool r1ref = true;
    bool onetoone = true;
    int64_t gate_chunk = 10000;
    int64_t gate_minrun = 20000;
    int64_t gap = 31000;
    int64_t clip = 10000;
    int64_t r2_tol = 100000;
    int64_t r2_slack = 1000000;
    int64_t lin_maxindel = 25000000;
    int64_t lin_ovtol = 50000;
    int64_t junction_min = 50000;
    bool guard = false;
    double guard_trip = 0.01;
    int64_t guard_minrun = 20000;
    int64_t guard_chunk = 10000;
    int64_t guard_minchunk = 2000;
    int64_t guard_svlen = 50;
};

// hard-coded, as in the emulator
static const int64_t NEW_MIN = 1000;            // a0new sub-interval minimum
static const int64_t BACKBONE_QOVERLAP = 1000;  // backbone query overlap allowance
static const double COLLIDE_FRAC = 0.5;         // R4: chunk collides at this fraction
static const double R1REF_COVER = 0.5;          // R1': judged when this much of a span sits on reference nodes
static const int MAX_ITERATIONS = 60;

// ---- intervals (half-open; merge joins overlapping and touching intervals)

static IVs merge(IVs v) {
    return merge_intervals(v);
}

static int64_t ivlen(const IVs& v) {
    int64_t t = 0;
    for (const IV& x : merge(v)) {
        t += x.second - x.first;
    }
    return t;
}

static IVs subtract(const IVs& ivs, const IVs& cut_in) {
    IVs cut = merge(cut_in);
    IVs out;
    for (const IV& iv : merge(ivs)) {
        int64_t cur = iv.first;
        for (const IV& c : cut) {
            if (c.second <= cur) continue;
            if (c.first >= iv.second) break;
            if (c.first > cur) out.push_back(make_pair(cur, c.first));
            cur = std::max(cur, c.second);
            if (cur >= iv.second) break;
        }
        if (cur < iv.second) out.push_back(make_pair(cur, iv.second));
    }
    return out;
}

static IVs intersect(const IVs& a_in, const IVs& b_in) {
    IVs a = merge(a_in), b = merge(b_in);
    IVs out;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        int64_t s = std::max(a[i].first, b[j].first), e = std::min(a[i].second, b[j].second);
        if (e > s) out.push_back(make_pair(s, e));
        if (a[i].second < b[j].second) ++i; else ++j;
    }
    return out;
}

// Python's str(round(x, places)), which the emulator's logs print
static string pyround(double x, int places) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.*f", places, x);
    string s = buf;
    if (s.find('.') != string::npos) {
        while (s.back() == '0') s.pop_back();
        if (s.back() == '.') s.push_back('0');
    }
    if (s == "-0.0") s = "0.0";
    return s;
}

static string fixed4(double x) {
    char buf[64];
    snprintf(buf, sizeof(buf), "%.4f", x);
    return buf;
}

// ---- chromosome names (reference contigs), interned

struct Chroms {
    vector<string> names;
    unordered_map<string, int> ids;
    int id(const string& name) {
        auto it = ids.find(name);
        if (it != ids.end()) return it->second;
        ids[name] = names.size();
        names.push_back(name);
        return names.size() - 1;
    }
    int find(const string& name) const {
        auto it = ids.find(name);
        return it == ids.end() ? -1 : it->second;
    }
};

// the reference contig a rank-0 node's SN names, as rgfa-split names it (and so as rc:Z: does):
// cactus's id=EVENT| prefix stripped, or failing that a PanSN SAMPLE#HAP# prefix
static string ref_chrom_name(const string& sn) {
    if (sn.compare(0, 3, "id=") == 0) {
        size_t p = sn.find('|', 3);
        if (p != string::npos) return sn.substr(p + 1);
    }
    size_t h1 = sn.find('#');
    if (h1 != string::npos) {
        size_t h2 = sn.find('#', h1 + 1);
        if (h2 != string::npos) return sn.substr(h2 + 1);
    }
    return sn;
}

struct Node {
    int64_t len;
    int chrom;          // reference contig of a rank-0 node, -1 otherwise
    int64_t so;
};

// ---- records

struct Place {
    bool ok = false;
    int chrom = -1;
    int64_t lo = 0, hi = 0;
    char ori = '+';
    int64_t entry = 0, exit = 0;
    bool interp = false;
};

struct Rec {
    int64_t idx;                // 0-based input position
    GafRecord* g;
    string path_text;
    int64_t qlen, qs, qe, ps, pe, matches, block, mapq;
    char strand;
    bool primary;
    string rc;
    double gi;                  // matches / block
    bool eligible = false;
    bool keep0 = false;         // the stock GAF stage keeps it
    int group = -1;             // chromosome group (rc, or its placement without one)
    IVs a0fin;                  // query intervals the stock chain's final PAF anchors with it
    // path steps
    vector<int64_t> soff, slen, sso, sref;
    vector<char> srev;
    vector<int> schrom;
    // cigar, query-increasing order, built on first use
    bool cig = false;
    vector<int64_t> QS, QE, PA, PB, N;
    vector<char> K, EQ;         // K: 0 aligned (=XM), 1 I, 2 D.  EQ: the op is '='
    map<IV, Place> place_cache;
};

// the canonical order of records: every tie-break uses it, never the input order
static int key_cmp(const Rec* a, const Rec* b) {
    if (a->g->query_name != b->g->query_name) return a->g->query_name < b->g->query_name ? -1 : 1;
    if (a->qs != b->qs) return a->qs < b->qs ? -1 : 1;
    if (a->qe != b->qe) return a->qe < b->qe ? -1 : 1;
    if (a->strand != b->strand) return a->strand < b->strand ? -1 : 1;
    if (a->path_text != b->path_text) return a->path_text < b->path_text ? -1 : 1;
    if (a->block != b->block) return a->block < b->block ? -1 : 1;
    if (a->matches != b->matches) return a->matches < b->matches ? -1 : 1;
    if (a->mapq != b->mapq) return a->mapq < b->mapq ? -1 : 1;
    if (a->primary != b->primary) return a->primary ? -1 : 1;
    return 0;
}

// Records equal in all of the key (a duplicated GAF line) are told apart by input position.  That
// cannot change the output: they tie each other wherever they meet, as they do in the stock filter,
// which deletes both, so the order only decides which idx the logs name
static bool key_less(const Rec* a, const Rec* b) {
    int c = key_cmp(a, b);
    return c != 0 ? c < 0 : a->idx < b->idx;
}

static bool key_equal(const Rec* a, const Rec* b) {
    return key_cmp(a, b) == 0;
}

static bool dom(const Rec* a, const Rec* b, double ratio, DomClause& clause) {
    return dominance(a->qs, a->qe, a->primary, (double)a->mapq, (double)a->block,
                     b->qs, b->qe, b->primary, (double)b->mapq, (double)b->block, ratio, clause);
}

static inline bool compat(const Rec* a, const Rec* b) {
    return a->rc == b->rc || a->rc.empty() || b->rc.empty();
}

static inline int64_t overlap(const Rec* a, const Rec* b) {
    return std::max((int64_t)0, std::min(a->qe, b->qe) - std::max(a->qs, b->qs));
}

static void build_cigar(Rec& r) {
    if (r.cig) return;
    r.cig = true;
    vector<pair<char, int64_t>> ops;
    for_each_cg(*r.g, [&](const char& c, const size_t& l) {
            ops.push_back(make_pair(c, (int64_t)l));
        });
    size_t m = ops.size();
    r.QS.resize(m); r.QE.resize(m); r.PA.resize(m); r.PB.resize(m); r.N.resize(m); r.K.resize(m); r.EQ.resize(m);
    int64_t p = r.ps;
    int64_t q = r.strand == '-' ? r.qe : r.qs;
    for (size_t i = 0; i < m; ++i) {
        char c = ops[i].first;
        int64_t n = ops[i].second;
        bool al = c == '=' || c == 'X' || c == 'M';
        bool ins = c == 'I';
        bool del = c == 'D';
        if (!al && !ins && !del) {
            cerr << "[gaffilter] error: -x cannot use cigar op " << c << " (record " << r.idx << ")" << endl;
            exit(1);
        }
        int64_t qn = (al || ins) ? n : 0;
        int64_t pn = (al || del) ? n : 0;
        r.PA[i] = p;
        r.PB[i] = p + pn;
        p += pn;
        if (r.strand == '-') {
            r.QE[i] = q;
            r.QS[i] = q - qn;
            q -= qn;
        } else {
            r.QS[i] = q;
            r.QE[i] = q + qn;
            q += qn;
        }
        r.N[i] = n;
        r.K[i] = al ? 0 : (ins ? 1 : 2);
        r.EQ[i] = c == '=';
    }
    if (r.strand == '-') {
        std::reverse(r.QS.begin(), r.QS.end()); std::reverse(r.QE.begin(), r.QE.end());
        std::reverse(r.PA.begin(), r.PA.end()); std::reverse(r.PB.begin(), r.PB.end());
        std::reverse(r.N.begin(), r.N.end()); std::reverse(r.K.begin(), r.K.end());
        std::reverse(r.EQ.begin(), r.EQ.end());
    }
}

// ops of r (query order) with QE >= s and QS < e
static inline pair<size_t, size_t> op_range(const Rec& r, int64_t s, int64_t e) {
    size_t i0 = std::lower_bound(r.QE.begin(), r.QE.end(), s) - r.QE.begin();
    size_t i1 = std::lower_bound(r.QS.begin(), r.QS.end(), e) - r.QS.begin();
    return make_pair(i0, i1);
}

// '=' columns over all columns of r's query [s,e) (a deletion counts where s <= its position < e).
// svlen > 0 leaves out indels of at least that size (the guard's identity)
static double local_ident(Rec& r, int64_t s, int64_t e, int64_t svlen = 0) {
    if (e <= s) return 0.0;
    build_cigar(r);
    pair<size_t, size_t> rg = op_range(r, s, e);
    int64_t eq = 0, cols = 0;
    for (size_t i = rg.first; i < rg.second; ++i) {
        int64_t ov;
        if (r.K[i] == 2) {
            ov = (r.QS[i] >= s && r.QS[i] < e) ? r.N[i] : 0;
        } else {
            ov = std::max((int64_t)0, std::min(r.QE[i], e) - std::max(r.QS[i], s));
        }
        if (svlen > 0 && r.K[i] != 0 && r.N[i] >= svlen) {
            ov = 0;
        }
        cols += ov;
        if (r.EQ[i]) eq += ov;
    }
    return cols ? (double)eq / (double)cols : 0.0;
}

// reference bp of the path before path offset p
static int64_t ref_at(const Rec& r, int64_t p) {
    int64_t i = (int64_t)(std::upper_bound(r.soff.begin(), r.soff.end(), p) - r.soff.begin()) - 1;
    i = std::max((int64_t)0, std::min(i, (int64_t)r.soff.size() - 1));
    int64_t inside = r.schrom[i] >= 0 ? std::max((int64_t)0, std::min(p - r.soff[i], r.slen[i])) : 0;
    return r.sref[i] + inside;
}

// identity of r's query [s,e) over its aligned columns on reference nodes only -> (identity, columns)
static pair<double, int64_t> ref_only_ident(Rec& r, int64_t s, int64_t e) {
    if (e <= s) return make_pair(0.0, (int64_t)0);
    build_cigar(r);
    pair<size_t, size_t> rg = op_range(r, s, e);
    int64_t eq = 0, al = 0;
    bool rev = r.strand == '-';
    for (size_t i = rg.first; i < rg.second; ++i) {
        if (r.K[i] != 0) continue;
        int64_t qa = std::max(r.QS[i], s), qb = std::min(r.QE[i], e);
        if (qb <= qa) continue;
        int64_t pa, pb;
        if (!rev) {
            pa = r.PA[i] + (qa - r.QS[i]);
            pb = r.PA[i] + (qb - r.QS[i]);
        } else {
            pa = r.PB[i] - (qb - r.QS[i]);
            pb = r.PB[i] - (qa - r.QS[i]);
        }
        int64_t refal = ref_at(r, pb) - ref_at(r, pa);
        al += refal;
        if (r.EQ[i]) eq += refal;
    }
    return make_pair(al > 0 ? (double)eq / (double)al : 0.0, al);
}

// path extent [p0,p1) of the aligned columns of r inside query [s,e)
static bool q_to_path(Rec& r, int64_t s, int64_t e, int64_t& p0, int64_t& p1) {
    s = std::max(s, r.qs);
    e = std::min(e, r.qe);
    if (e <= s) return false;
    if (s == r.qs && e == r.qe) {
        p0 = r.ps;
        p1 = r.pe;
        return true;
    }
    build_cigar(r);
    size_t i0 = std::upper_bound(r.QE.begin(), r.QE.end(), s) - r.QE.begin();
    size_t i1 = std::lower_bound(r.QS.begin(), r.QS.end(), e) - r.QS.begin();
    bool rev = r.strand == '-';
    bool found = false;
    int64_t lo = 0, hi = 0;
    auto visit = [&](size_t i) -> bool {
        if (r.K[i] != 0) return false;
        int64_t a = std::max(r.QS[i], s), b = std::min(r.QE[i], e);
        if (b <= a) return false;
        int64_t x0, x1;
        if (!rev) {
            x0 = r.PA[i] + (a - r.QS[i]);
            x1 = r.PA[i] + (b - r.QS[i]);
        } else {
            x0 = r.PB[i] - (b - r.QS[i]);
            x1 = r.PB[i] - (a - r.QS[i]);
        }
        lo = found ? std::min(lo, x0) : x0;
        hi = found ? std::max(hi, x1) : x1;
        found = true;
        return true;
    };
    for (size_t i = i0; i < i1; ++i) {
        if (visit(i)) break;
    }
    for (size_t i = i1; i > i0; --i) {
        if (visit(i - 1)) break;
    }
    if (!found) return false;
    p0 = lo;
    p1 = hi;
    return true;
}

static inline size_t first_step(const Rec& r, int64_t p0) {
    int64_t j = (int64_t)(std::upper_bound(r.soff.begin(), r.soff.end(), p0) - r.soff.begin()) - 1;
    return (size_t)std::max((int64_t)0, j);
}

// the reference intervals of contig c that r's query [s,e) sits on (whole steps: indels inside a
// step are not split out)
static IVs ref_ivs(Rec& r, int64_t s, int64_t e, int c) {
    IVs out;
    int64_t p0, p1;
    if (!q_to_path(r, s, e, p0, p1)) return out;
    for (size_t j = first_step(r, p0); j < r.soff.size() && r.soff[j] < p1; ++j) {
        int64_t off = r.soff[j], L = r.slen[j];
        int64_t a = std::max(p0, off) - off, b = std::min(p1, off + L) - off;
        if (b > a && r.schrom[j] == c) {
            if (!r.srev[j]) out.push_back(make_pair(r.sso[j] + a, r.sso[j] + b));
            else out.push_back(make_pair(r.sso[j] + L - b, r.sso[j] + L - a));
        }
    }
    return merge(out);
}

// reference placement of r's query [s,e): dominant contig, its extent there, orientation, and the
// reference positions where the query enters and leaves it.  Falls back on the nearest reference
// steps along the path when the span sits on non-reference nodes only
static Place place(Rec& r, int64_t s, int64_t e, const Chroms& chroms) {
    IV key(s, e);
    auto cached = r.place_cache.find(key);
    if (cached != r.place_cache.end()) return cached->second;
    Place res;
    int64_t p0, p1;
    if (q_to_path(r, s, e, p0, p1)) {
        struct Seg { int chrom; int64_t rl, rh; bool fwd; };
        vector<Seg> segs;
        size_t j0 = first_step(r, p0);
        for (size_t j = j0; j < r.soff.size(); ++j) {
            int64_t off = r.soff[j], L = r.slen[j];
            if (off >= p1) break;
            int64_t a = std::max(p0, off) - off, b = std::min(p1, off + L) - off;
            if (b <= a || r.schrom[j] < 0) continue;
            Seg sg;
            sg.chrom = r.schrom[j];
            if (!r.srev[j]) { sg.rl = r.sso[j] + a; sg.rh = r.sso[j] + b; }
            else { sg.rl = r.sso[j] + L - b; sg.rh = r.sso[j] + L - a; }
            sg.fwd = (r.strand == '+') != (bool)r.srev[j];
            segs.push_back(sg);
        }
        bool interp = false;
        if (segs.empty()) {
            int64_t jl = -1, jr = -1;
            for (int64_t j = (int64_t)j0; j >= 0; --j) {
                if (r.schrom[j] >= 0 && r.soff[j] + r.slen[j] <= p0) { jl = j; break; }
            }
            for (size_t j = j0; j < r.soff.size(); ++j) {
                if (r.schrom[j] >= 0 && r.soff[j] >= p1) { jr = j; break; }
            }
            for (int side = 0; side < 2; ++side) {
                int64_t j = side == 0 ? jl : jr;
                if (j < 0) continue;
                bool rev = r.srev[j];
                int64_t x;
                // the reference base adjacent to the non-reference run
                if (side == 0) x = !rev ? r.sso[j] + r.slen[j] : r.sso[j];
                else x = !rev ? r.sso[j] : r.sso[j] + r.slen[j];
                segs.push_back({r.schrom[j], x, x, (r.strand == '+') != rev});
            }
            interp = true;
        }
        if (!segs.empty()) {
            // the contig with the most reference bp, ties to the larger name
            map<int, int64_t> byc;
            for (const Seg& sg : segs) byc[sg.chrom] += std::max((int64_t)1, sg.rh - sg.rl);
            int best = -1;
            for (const auto& kv : byc) {
                if (best < 0 || kv.second > byc[best] ||
                    (kv.second == byc[best] && chroms.names[kv.first] > chroms.names[best])) {
                    best = kv.first;
                }
            }
            vector<Seg> mine;
            int64_t fb = 0, rb = 0;
            for (const Seg& sg : segs) {
                if (sg.chrom != best) continue;
                mine.push_back(sg);
                (sg.fwd ? fb : rb) += std::max((int64_t)1, sg.rh - sg.rl);
            }
            res.ok = true;
            res.chrom = best;
            res.lo = mine[0].rl;
            res.hi = mine[0].rh;
            for (const Seg& sg : mine) {
                res.lo = std::min(res.lo, sg.rl);
                res.hi = std::max(res.hi, sg.rh);
            }
            res.ori = fb >= rb ? '+' : '-';
            res.interp = interp;
            // the path runs along the query on a '+' record, against it on a '-' one; along the
            // query a forward segment runs lo->hi and a reverse one hi->lo
            const Seg& fq = r.strand == '+' ? mine.front() : mine.back();
            const Seg& lq = r.strand == '+' ? mine.back() : mine.front();
            res.entry = fq.fwd ? fq.rl : fq.rh;
            res.exit = lq.fwd ? lq.rh : lq.rl;
        }
    }
    r.place_cache[key] = res;
    return res;
}

static string place_str(const Place& p, const Chroms& chroms) {
    if (!p.ok) return ".";
    stringstream ss;
    ss << chroms.names[p.chrom] << ":" << p.lo << "-" << p.hi << p.ori << (p.interp ? "~" : "");
    return ss.str();
}

// n roughly equal chunks of about size bp (round half to even, as the emulator's Python does)
static IVs chunks(int64_t s, int64_t e, int64_t size) {
    IVs out;
    int64_t n = std::max((int64_t)1, (int64_t)std::nearbyint((double)(e - s) / (double)size));
    double step = (double)(e - s) / (double)n;
    for (int64_t i = 0; i < n; ++i) {
        int64_t a = s + (int64_t)std::nearbyint((double)i * step);
        int64_t b = i < n - 1 ? s + (int64_t)std::nearbyint((double)(i + 1) * step) : e;
        if (b > a) out.push_back(make_pair(a, b));
    }
    return out;
}

// ---- one query contig

// one kept piece of a record, as the structure test sees it
struct Entry {
    int ri;                     // index into the contig's E
    int64_t s, e;
    Place pl;
    double li;
    bool rem;                   // a remainder: a piece of a demoted record
    IVs a0new;                  // >= NEW_MIN bp the stock chain does not anchor with this record
    bool is_new() const { return rem || !a0new.empty(); }
};

enum Sided { S_TWO, S_ISOLATED, S_ONE, S_ALONE, S_END_JOINED, S_END_ISOLATED };
static const char* sided_name[] = {"two", "isolated", "one", "alone", "end-joined", "end-isolated"};

struct Excursion {
    vector<int> idxs;
    int L = -1, R = -1;
    bool has_fp = false;
    int64_t flo = 0, fhi = 0;
    bool bounded = false;
    bool has_gl = false, has_gr = false;
    int64_t gl = 0, gr = 0;
    int64_t ig = 0;             // widest query gap between the excursion's own pieces
    Sided sided;
};

struct Shared {
    const Params& P;
    Chroms& chroms;
    vector<string>& log_rows;
    vector<string>& junction_rows;
    // summary
    int64_t isolations = 0, gate_drops = 0, guard_applied = 0, guard_vetoes = 0, capped = 0;
    map<string, int64_t> hole_bp;
    Shared(const Params& p, Chroms& c, vector<string>& l, vector<string>& j) : P(p), chroms(c), log_rows(l), junction_rows(j) {}
};

class Contig {
public:
    Contig(Shared& sh, const string& q, vector<Rec*>& eligible, function<IVs(int)> hco_fn) :
        sh(sh), P(sh.P), q(q), E(eligible), hco_fn(hco_fn) {}

    // final kept pieces, per E index, and whether it was demoted
    vector<IVs> pieces;
    vector<char> dem;
    int iterations = 0;
    bool capped = false;

    void run();

private:
    Shared& sh;
    const Params& P;
    const string& q;
    vector<Rec*>& E;
    function<IVs(int)> hco_fn;
    vector<IVs> excluded;
    vector<vector<pair<IV, string>>> why;          // drop reasons, for the hole summary
    map<pair<int, int>, IVs> veto;
    set<pair<int, int>> veto_used;
    vector<int> groups;
    map<int, set<vector<int64_t>>> a0sig;
    map<int, IVs> hco_cache;

    bool beats(int i, int j) const {
        DomClause c;
        return dom(E[i], E[j], P.ratio, c);
    }
    // is piece p new at one end (its end if end, else its start)?  A remainder is new throughout;
    // otherwise some of its a0new lies within --gap of that end.  The junction log's test, and
    // isolation's for a backbone neighbour
    bool side_new(const Entry& p, bool end) const {
        if (p.rem) return true;
        for (const IV& x : p.a0new) {
            if (end ? (x.second >= p.e - P.gap && x.first < p.e) : (x.first <= p.s + P.gap && x.second > p.s)) return true;
        }
        return false;
    }
    void log_row(const string& kind, int64_t a, int64_t b, const string& reason, int it, int64_t idx,
                 const string& pl, const string& li) {
        stringstream ss;
        ss << kind << "\t" << q << "\t" << a << "\t" << b << "\t" << reason << "\t" << it << "\t" << idx;
        if (!pl.empty()) ss << "\t" << pl << "\t" << li;
        sh.log_rows.push_back(ss.str());
    }
    void drop(const Entry& p, int64_t a, int64_t b, const string& reason, int it) {
        excluded[p.ri].push_back(make_pair(a, b));
        string kind = reason.substr(0, reason.find(':'));
        why[p.ri].push_back(make_pair(make_pair(a, b), kind));
        log_row(kind, a, b, reason, it, E[p.ri]->idx, place_str(p.pl, sh.chroms), pyround(p.li, 4));
    }
    const IVs& hco(int c) {
        auto it = hco_cache.find(c);
        if (it == hco_cache.end()) it = hco_cache.insert(make_pair(c, hco_fn(c))).first;
        return it->second;
    }
    void init_guard();
    vector<IVs> resolve(bool track_veto);
    vector<Entry> piece_table(int c, const vector<IVs>& pcs);
    bool gates(int c, vector<Entry>& pl, int it);
    vector<int> backbone(const vector<Entry>& pl, int c) const;
    vector<Excursion> excursions(const vector<Entry>& pl, const vector<int>& on, bool newgap) const;
    Excursion judge(const vector<Entry>& pl, const vector<int>& grp, int L, int R, bool newgap) const;
    bool colinear(const Entry& a, const Entry& b) const;
    vector<int64_t> sig(const Excursion& ex, const vector<Entry>& pl) const;
    // (action, reason): action 0 pass, 1 isolate, 2 exclude
    pair<int, string> reality(const Excursion& ex, const vector<Entry>& pl) const;
    bool structure(int c, const vector<Entry>& pl, int it, vector<Excursion>& exs_out);
    void junctions(int c, const vector<Entry>& pl, const vector<Excursion>& exs);
    void holes();
};

void Contig::init_guard() {
    if (!P.guard) return;
    for (int w = 0; w < (int)E.size(); ++w) {
        for (int l = 0; l < (int)E.size(); ++l) {
            if (w == l || !compat(E[w], E[l]) || overlap(E[w], E[l]) == 0) continue;
            // in scope: the arm's ratio decides the contest by MAPQ or length, the stock one ties
            DomClause cw, cs;
            if (!dom(E[w], E[l], P.ratio, cw) || (cw != DOM_MAPQ && cw != DOM_LENGTH)) continue;
            if (dom(E[w], E[l], P.stock_ratio_d, cs) || cs != DOM_TIE) continue;
            Rec& W = *E[w];
            Rec& Lr = *E[l];
            int64_t S = std::max(W.qs, Lr.qs), En = std::min(W.qe, Lr.qe);
            struct Ch { int64_t a, b; double m; };
            vector<Ch> ch;
            for (int64_t a = S; a < En; a += P.guard_chunk) {
                int64_t b = std::min(a + P.guard_chunk, En);
                if (b - a >= P.guard_minchunk) {
                    ch.push_back({a, b, local_ident(Lr, a, b, P.guard_svlen) - local_ident(W, a, b, P.guard_svlen)});
                }
            }
            size_t n = ch.size();
            vector<char> v(n, 0);
            // trip: runs of margin >= trip spanning >= minrun
            for (size_t i = 0; i < n;) {
                if (ch[i].m >= P.guard_trip) {
                    size_t j = i;
                    while (j + 1 < n && ch[j + 1].m >= P.guard_trip) ++j;
                    if (ch[j].b - ch[i].a >= P.guard_minrun) {
                        for (size_t k = i; k <= j; ++k) v[k] = 1;
                    }
                    i = j + 1;
                } else {
                    ++i;
                }
            }
            // extend right, then left, while the loser is still better
            for (size_t k = 1; k < n; ++k) {
                if (v[k - 1] && !v[k] && ch[k].m > 0) v[k] = 1;
            }
            for (int64_t k = (int64_t)n - 2; k >= 0; --k) {
                if (v[k + 1] && !v[k] && ch[k].m > 0) v[k] = 1;
            }
            IVs runs;
            for (size_t k = 0; k < n; ++k) {
                if (!v[k]) continue;
                if (!runs.empty() && runs.back().second == ch[k].a) runs.back().second = ch[k].b;
                else runs.push_back(make_pair(ch[k].a, ch[k].b));
            }
            // merge runs closer than the isolation gap
            IVs merged;
            for (const IV& x : runs) {
                if (!merged.empty() && x.first - merged.back().second < P.gap) merged.back().second = x.second;
                else merged.push_back(x);
            }
            if (merged.empty()) continue;
            // winner islands shorter than the gap, at a span end where the winner record ends
            if (W.qs == S && merged.front().first - S > 0 && merged.front().first - S < P.gap) merged.front().first = S;
            if (W.qe == En && En - merged.back().second > 0 && En - merged.back().second < P.gap) merged.back().second = En;
            veto[make_pair(w, l)] = merged;
            for (const IV& x : merged) {
                double mx = 0;
                bool any = false;
                for (const Ch& c : ch) {
                    if (c.a < x.second && c.b > x.first) {
                        mx = any ? std::max(mx, c.m) : c.m;
                        any = true;
                    }
                }
                stringstream reason;
                reason << "guard:" << W.idx << ">" << Lr.idx << ":" << fixed4(mx);
                log_row("guard", x.first, x.second, reason.str(), 0, W.idx, ".", "0");
                ++sh.guard_vetoes;
            }
        }
    }
}

// every elementary segment goes to the claimant that beats every compatible other claimant, if one does
vector<IVs> Contig::resolve(bool track_veto) {
    vector<int64_t> bps;
    for (Rec* r : E) {
        bps.push_back(r->qs);
        bps.push_back(r->qe);
    }
    for (const IVs& x : excluded) {
        for (const IV& iv : x) { bps.push_back(iv.first); bps.push_back(iv.second); }
    }
    for (const auto& kv : veto) {
        for (const IV& iv : kv.second) { bps.push_back(iv.first); bps.push_back(iv.second); }
    }
    std::sort(bps.begin(), bps.end());
    bps.erase(std::unique(bps.begin(), bps.end()), bps.end());
    vector<IVs> keep(E.size());
    if (track_veto) veto_used.clear();
    size_t next = 0;
    vector<int> active, C, C2;
    for (size_t k = 0; k + 1 < bps.size(); ++k) {
        int64_t s = bps[k], e = bps[k + 1];
        while (next < E.size() && E[next]->qs <= s) active.push_back(next++);
        // E is in key order, which within one query is query-start order, so active stays in it
        active.erase(std::remove_if(active.begin(), active.end(), [&](int i) { return E[i]->qe < e; }), active.end());
        C.clear();
        for (int i : active) {
            bool ex = false;
            for (const IV& iv : excluded[i]) {
                if (iv.first <= s && e <= iv.second) { ex = true; break; }
            }
            if (!ex) C.push_back(i);
        }
        // tier: a demoted record yields the segment to any compatible non-demoted claimant
        C2.clear();
        for (int i : C) {
            bool shadowed = false;
            if (dem[i]) {
                for (int o : C) {
                    if (!dem[o] && compat(E[i], E[o])) { shadowed = true; break; }
                }
            }
            if (!shadowed) C2.push_back(i);
        }
        for (int r : C2) {
            bool ok = true;
            for (int o : C2) {
                if (o == r || !compat(E[r], E[o])) continue;
                if (!beats(r, o)) { ok = false; break; }
            }
            if (ok && !veto.empty()) {
                // the guard: only contests as decided, after the tier (a tier-shadowed loser is none)
                for (int o : C2) {
                    if (o == r || !compat(E[r], E[o])) continue;
                    auto vt = veto.find(make_pair(r, o));
                    if (vt == veto.end()) continue;
                    bool hit = false;
                    for (const IV& iv : vt->second) {
                        if (iv.first <= s && e <= iv.second) { hit = true; break; }
                    }
                    if (hit) {
                        ok = false;
                        if (track_veto) veto_used.insert(make_pair(r, o));
                        break;
                    }
                }
            }
            if (ok) keep[r].push_back(make_pair(s, e));
        }
    }
    for (IVs& x : keep) x = merge(x);
    return keep;
}

vector<Entry> Contig::piece_table(int c, const vector<IVs>& pcs) {
    vector<Entry> out;
    for (int i = 0; i < (int)E.size(); ++i) {
        Rec& r = *E[i];
        if (r.group != c) continue;
        for (const IV& iv : pcs[i]) {
            Entry p;
            p.ri = i;
            p.s = iv.first;
            p.e = iv.second;
            p.pl = place(r, p.s, p.e, sh.chroms);
            p.rem = dem[i];
            p.li = dem[i] ? local_ident(r, p.s, p.e) : r.gi;
            if (!dem[i]) {
                for (const IV& x : subtract({iv}, r.a0fin)) {
                    if (x.second - x.first >= NEW_MIN) p.a0new.push_back(x);
                }
            }
            out.push_back(p);
        }
    }
    std::stable_sort(out.begin(), out.end(), [](const Entry& a, const Entry& b) {
            if (a.s != b.s) return a.s < b.s;
            if (a.e != b.e) return a.e < b.e;
            return a.ri < b.ri;
        });
    return out;
}

// R1' and R4 on every new piece: identity to the reference copy it is placed on, and one-to-one
bool Contig::gates(int c, vector<Entry>& pl, int it) {
    vector<IVs> cov_of(pl.size());
    for (size_t k = 0; k < pl.size(); ++k) {
        cov_of[k] = ref_ivs(*E[pl[k].ri], pl[k].s, pl[k].e, c);
    }
    struct Drop { size_t k; int64_t a, b; string why; };
    vector<Drop> drops;
    for (size_t k = 0; k < pl.size(); ++k) {
        const Entry& p = pl[k];
        if (!p.is_new()) continue;
        Rec& r = *E[p.ri];
        IVs spans = p.rem ? IVs{make_pair(p.s, p.e)} : merge(p.a0new);
        IVs other_v = hco(c);
        for (size_t j = 0; j < pl.size(); ++j) {
            if (j != k) other_v.insert(other_v.end(), cov_of[j].begin(), cov_of[j].end());
        }
        IVs other = merge(other_v);
        for (const IV& span : spans) {
            int64_t a0 = span.first, b0 = span.second;
            if (P.r1ref) {
                pair<double, int64_t> rio = ref_only_ident(r, a0, b0);
                if ((double)rio.second >= R1REF_COVER * (double)(b0 - a0) && rio.first < P.id_floor) {
                    drops.push_back({k, a0, b0, "r1ref:" + fixed4(rio.first)});
                    continue;
                }
            }
            if (!P.onetoone) continue;
            IVs bad;
            vector<double> col_of;
            for (const IV& ch : chunks(a0, b0, P.gate_chunk)) {
                IVs civ = ref_ivs(r, ch.first, ch.second, c);
                int64_t L = ivlen(civ);
                double col = L ? (double)ivlen(intersect(civ, other)) / (double)L : 0.0;
                if (col >= COLLIDE_FRAC) {
                    bad.push_back(ch);
                    col_of.push_back(col);
                }
            }
            IVs runs = merge(bad);
            for (const IV& run : runs) {
                if (run.second - run.first >= P.gate_minrun) {
                    // the first colliding chunk inside the run gives the logged fraction
                    double inf = 0;
                    for (size_t x = 0; x < bad.size(); ++x) {
                        if (bad[x].first >= run.first && bad[x].second <= run.second) { inf = col_of[x]; break; }
                    }
                    drops.push_back({k, run.first, run.second, "collide:" + pyround(inf, 3)});
                }
            }
        }
    }
    for (const Drop& d : drops) {
        drop(pl[d.k], d.a, d.b, d.why, it);
        ++sh.gate_drops;
    }
    return !drops.empty();
}

bool Contig::colinear(const Entry& a, const Entry& b) const {
    if (!a.pl.ok || !b.pl.ok || a.pl.chrom != b.pl.chrom || a.pl.ori != b.pl.ori) return false;
    int64_t dq = b.s - a.e;
    int64_t dr = a.pl.ori == '+' ? b.pl.lo - a.pl.hi : a.pl.lo - b.pl.hi;
    if (dr < -P.lin_ovtol) return false;
    return std::abs(dr - dq) <= P.lin_maxindel;
}

// the max-weight (query bp) colinear chain of placed pieces on contig c, in table order
vector<int> Contig::backbone(const vector<Entry>& pl, int c) const {
    size_t n = pl.size();
    vector<int> on;
    if (n == 0) return on;
    vector<int64_t> best(n, 0);
    vector<int64_t> prev(n, -1);
    for (size_t i = 0; i < n; ++i) {
        int64_t w = pl[i].e - pl[i].s;
        if (!pl[i].pl.ok || pl[i].pl.chrom != c) {
            best[i] = -1;
            continue;
        }
        best[i] = w;
        for (size_t j = 0; j < i; ++j) {
            if (best[j] < 0 || pl[j].e > pl[i].s + BACKBONE_QOVERLAP) continue;
            if (colinear(pl[j], pl[i]) && best[j] + w > best[i]) {
                best[i] = best[j] + w;
                prev[i] = j;
            }
        }
    }
    size_t top = 0;
    for (size_t i = 1; i < n; ++i) {
        if (best[i] > best[top]) top = i;
    }
    if (best[top] < 0) return on;
    for (int64_t i = top; i >= 0; i = prev[i]) on.push_back(i);
    std::sort(on.begin(), on.end());
    return on;
}

Excursion Contig::judge(const vector<Entry>& pl, const vector<int>& grp, int L, int R, bool newgap) const {
    Excursion d;
    d.idxs = grp;
    d.L = L;
    d.R = R;
    const Place* pL = L >= 0 && pl[L].pl.ok ? &pl[L].pl : nullptr;
    const Place* pR = R >= 0 && pl[R].pl.ok ? &pl[R].pl : nullptr;
    vector<int64_t> bpts;
    if (pL) bpts.push_back(pL->exit);
    if (pR) bpts.push_back(pR->entry);
    int bchrom = pL ? pL->chrom : (pR ? pR->chrom : -1);
    bool any_pl = false, same_chrom = true;
    int64_t flo = 0, fhi = 0;
    for (int k : grp) {
        const Place& p = pl[k].pl;
        if (!p.ok) continue;
        flo = any_pl ? std::min(flo, p.lo) : p.lo;
        fhi = any_pl ? std::max(fhi, p.hi) : p.hi;
        any_pl = true;
        if (p.chrom != bchrom) same_chrom = false;
    }
    if (any_pl && same_chrom) {
        d.has_fp = true;
        d.flo = flo;
        d.fhi = fhi;
    }
    if (!bpts.empty() && d.has_fp) {
        int64_t wlo = *std::min_element(bpts.begin(), bpts.end());
        int64_t whi = *std::max_element(bpts.begin(), bpts.end());
        d.bounded = d.flo >= wlo - P.r2_tol && d.fhi <= whi + P.r2_tol && (whi - wlo) <= (d.fhi - d.flo) + P.r2_slack;
    }
    // the widest query gap inside the excursion: unanchored query between its pieces, which the
    // clip downstream removes
    int64_t reach = pl[grp.front()].e;
    for (size_t i = 1; i < grp.size(); ++i) {
        d.ig = std::max(d.ig, pl[grp[i]].s - reach);
        reach = std::max(reach, pl[grp[i]].e);
    }
    // sidedness: is each junction path-continuous?  Up to the clip threshold always; with a new
    // piece on either side, anything under the isolation gap counts
    if (L >= 0) { d.has_gl = true; d.gl = pl[grp.front()].s - pl[L].e; }
    if (R >= 0) { d.has_gr = true; d.gr = pl[R].s - pl[grp.back()].e; }
    bool nl = newgap && L >= 0 && (pl[L].is_new() || pl[grp.front()].is_new());
    bool nr = newgap && R >= 0 && (pl[R].is_new() || pl[grp.back()].is_new());
    bool cl = d.has_gl && (d.gl <= P.clip || (nl && d.gl < P.gap));
    bool cr = d.has_gr && (d.gr <= P.clip || (nr && d.gr < P.gap));
    if (L >= 0 && R >= 0) d.sided = (cl && cr) ? S_TWO : ((!cl && !cr) ? S_ISOLATED : S_ONE);
    else if (L < 0 && R < 0) d.sided = S_ALONE;
    else d.sided = (cl || cr) ? S_END_JOINED : S_END_ISOLATED;
    return d;
}

vector<Excursion> Contig::excursions(const vector<Entry>& pl, const vector<int>& on, bool newgap) const {
    vector<Excursion> out;
    vector<char> is_on(pl.size(), 0);
    for (int k : on) is_on[k] = 1;
    int n = pl.size();
    int lastL = -1;
    for (int i = 0; i < n;) {
        if (is_on[i]) { lastL = i; ++i; continue; }
        int j = i;
        while (j + 1 < n && !is_on[j + 1]) ++j;
        int R = j + 1 < n ? j + 1 : -1;
        vector<int> grp;
        for (int k = i; k <= j; ++k) grp.push_back(k);
        out.push_back(judge(pl, grp, lastL, R, newgap));
        i = j + 1;
    }
    return out;
}

vector<int64_t> Contig::sig(const Excursion& ex, const vector<Entry>& pl) const {
    vector<int64_t> v;
    v.push_back(ex.idxs.size());
    for (int k : ex.idxs) { v.push_back(pl[k].ri); v.push_back(pl[k].s); v.push_back(pl[k].e); }
    v.push_back(ex.sided);
    for (int k : {ex.L, ex.R}) {
        if (k < 0) { v.push_back(-1); v.push_back(-1); v.push_back(-1); }
        else { v.push_back(pl[k].ri); v.push_back(pl[k].s); v.push_back(pl[k].e); }
    }
    return v;
}

pair<int, string> Contig::reality(const Excursion& ex, const vector<Entry>& pl) const {
    for (int k : ex.idxs) {
        if (pl[k].rem && pl[k].li < P.id_floor) return make_pair(2, "R1:ident" + fixed4(pl[k].li));
    }
    if (ex.sided == S_ISOLATED || ex.sided == S_END_ISOLATED || ex.sided == S_ALONE) return make_pair(0, string("isolated"));
    // two-sided only until the clip removes the unanchored query inside it: that leaves two halves
    // joined on one side each, which no later stage severs
    if (ex.sided == S_TWO && ex.ig >= P.gap) return make_pair(1, string("R2:gapped-two"));
    if (ex.sided == S_TWO && ex.bounded && ex.has_fp) return make_pair(0, string("bounded-two"));
    if (ex.sided == S_TWO) return make_pair(1, string("R2:unbounded-two"));
    if (ex.sided == S_ONE) return make_pair(1, string("R2:one-sided"));
    return make_pair(1, string("R2:") + sided_name[ex.sided]);
}

bool Contig::structure(int c, const vector<Entry>& pl, int it, vector<Excursion>& exs) {
    vector<int> on = backbone(pl, c);
    exs = excursions(pl, on, true);
    bool changed = false;
    const set<vector<int64_t>>& sigs = a0sig[c];
    for (const Excursion& ex : exs) {
        vector<int> nw;
        for (int k : ex.idxs) if (pl[k].is_new()) nw.push_back(k);
        bool nbnew = (ex.L >= 0 && pl[ex.L].is_new()) || (ex.R >= 0 && pl[ex.R].is_new());
        bool chg = !sigs.count(sig(ex, pl));
        if (nw.empty() && !nbnew && !chg) continue;
        pair<int, string> act = reality(ex, pl);
        if (act.first == 0) continue;
        if (act.first == 2) {
            // only the new sequence under the identity floor goes (piecewise); the excursion is
            // tested again on the next round.  A remainder is new throughout and goes whole.  A piece
            // of a record the stock rule keeps loses only the stretches the stock chain does not
            // anchor, never its stock-anchored rest (its identity is the whole record's)
            for (int k : nw) {
                string reason = "reality:" + act.second;
                if (pl[k].li >= P.id_floor) {
                    log_row("r1-spared", pl[k].s, pl[k].e, reason, it, E[pl[k].ri]->idx,
                            place_str(pl[k].pl, sh.chroms), pyround(pl[k].li, 4));
                } else if (pl[k].rem) {
                    drop(pl[k], pl[k].s, pl[k].e, reason, it);
                } else {
                    for (const IV& x : pl[k].a0new) drop(pl[k], x.first, x.second, reason, it);
                }
            }
            if (!nw.empty()) changed = true;
            continue;
        }
        // isolate: break each path-continuous junction, trimming its new side back to the gap.  The
        // excursion's own entry is trimmed if any of it is new.  Otherwise the backbone neighbour is,
        // but only if it is new at the end that faces the junction: a backbone record that is new
        // somewhere far from the junction (for instance a stretch the stock line filter drops) must
        // not lose stock-anchored sequence at it.  If neither is new there, the junction is the
        // stock chain's own, and is logged as a0-junction
        struct Todo { int inner, outer; bool inner_start; int64_t gap; };
        vector<Todo> todo;
        if (ex.L >= 0 && ex.has_gl && ex.gl < P.gap) todo.push_back({ex.idxs.front(), ex.L, true, ex.gl});
        if (ex.R >= 0 && ex.has_gr && ex.gr < P.gap) todo.push_back({ex.idxs.back(), ex.R, false, ex.gr});
        for (const Todo& t : todo) {
            int tp;
            bool start;
            if (pl[t.inner].is_new()) {
                tp = t.inner;
                start = t.inner_start;
            } else if (side_new(pl[t.outer], t.inner_start)) {
                tp = t.outer;
                start = !t.inner_start;
            } else {
                log_row("a0-junction", pl[t.inner].s, pl[t.inner].e, act.second, it, E[pl[t.inner].ri]->idx, ".", "0");
                continue;
            }
            int64_t cut = P.gap - std::max((int64_t)0, t.gap);
            int64_t a, b;
            if (pl[tp].e - pl[tp].s - cut < P.rmin) { a = pl[tp].s; b = pl[tp].e; }
            else if (start) { a = pl[tp].s; b = pl[tp].s + cut; }
            else { a = pl[tp].e - cut; b = pl[tp].e; }
            drop(pl[tp], a, b, "isolate:" + act.second, it);
            ++sh.isolations;
            changed = true;
        }
    }
    return changed;
}

// review log: every junction under the gap that joins pieces placed apart, and what it is
void Contig::junctions(int c, const vector<Entry>& pl, const vector<Excursion>& exs) {
    set<pair<int, int>> exempt;
    for (const Excursion& ex : exs) {
        pair<int, string> act = reality(ex, pl);
        if (act.first == 0 && act.second == "bounded-two") {
            vector<int> ks;
            if (ex.L >= 0) ks.push_back(ex.L);
            ks.insert(ks.end(), ex.idxs.begin(), ex.idxs.end());
            if (ex.R >= 0) ks.push_back(ex.R);
            for (size_t i = 0; i + 1 < ks.size(); ++i) exempt.insert(make_pair(ks[i], ks[i + 1]));
        }
    }
    for (size_t k = 0; k + 1 < pl.size(); ++k) {
        const Entry& p = pl[k];
        const Entry& qn = pl[k + 1];
        int64_t gap = qn.s - p.e;
        if (gap >= P.gap) continue;
        bool np = side_new(p, true), nq = side_new(qn, false);
        if (!p.pl.ok || !qn.pl.ok) continue;
        string kind;
        bool has_dr = false;
        int64_t dr = 0;
        if (p.pl.chrom != qn.pl.chrom) kind = "cross";
        else if (p.pl.ori != qn.pl.ori) kind = "inv";
        else {
            has_dr = true;
            dr = p.pl.ori == '+' ? qn.pl.entry - p.pl.exit : p.pl.exit - qn.pl.entry;
            kind = dr >= P.junction_min ? "fjump" : (dr <= -P.junction_min ? "back" : "local");
        }
        if (kind == "local") continue;
        bool ex = exempt.count(make_pair((int)k, (int)k + 1));
        string act;
        if (!np && !nq) act = "a0-structure";
        else if (ex) act = "exempt:bounded-two";
        else if (kind == "inv") act = "inv-not-tested";
        else {
            // what a junction rule would cut: the new side (the shorter if both are)
            const Entry* tp;
            bool start;
            if (np && nq) {
                if (p.e - p.s <= qn.e - qn.s) { tp = &p; start = false; } else { tp = &qn; start = true; }
            } else if (np) { tp = &p; start = false; }
            else { tp = &qn; start = true; }
            int64_t cut = P.gap - std::max((int64_t)0, gap);
            int64_t a, b;
            if (tp->e - tp->s - cut < P.rmin) { a = tp->s; b = tp->e; }
            else if (start) { a = tp->s; b = tp->s + cut; }
            else { a = tp->e - cut; b = tp->e; }
            stringstream ss;
            ss << "trim:" << E[tp->ri]->idx << ":" << a << "-" << b << ":"
               << ((a == tp->s && b == tp->e) ? "whole" : (start ? "start" : "end"));
            act = ss.str();
        }
        stringstream row;
        row << q << "\t" << sh.chroms.names[c] << "\t" << E[p.ri]->idx << "\t" << p.s << "\t" << p.e << "\t"
            << E[qn.ri]->idx << "\t" << qn.s << "\t" << qn.e << "\t" << (int)np << "\t" << (int)nq << "\t" << gap << "\t"
            << kind << "\t";
        if (has_dr) row << dr; else row << ".";
        row << "\t" << sh.chroms.names[p.pl.chrom] << "\t" << p.pl.exit << "\t" << sh.chroms.names[qn.pl.chrom] << "\t"
            << qn.pl.entry << "\t" << (int)ex << "\t" << act;
        sh.junction_rows.push_back(row.str());
    }
}

// query bp eligible records claim but no piece keeps, by cause (summary only)
void Contig::holes() {
    IVs kept;
    vector<int64_t> bps;
    for (int i = 0; i < (int)E.size(); ++i) {
        bps.push_back(E[i]->qs);
        bps.push_back(E[i]->qe);
        for (const IV& x : excluded[i]) { bps.push_back(x.first); bps.push_back(x.second); }
        for (const IV& x : pieces[i]) { kept.push_back(x); bps.push_back(x.first); bps.push_back(x.second); }
    }
    kept = merge(kept);
    std::sort(bps.begin(), bps.end());
    bps.erase(std::unique(bps.begin(), bps.end()), bps.end());
    size_t kp = 0;
    for (size_t k = 0; k + 1 < bps.size(); ++k) {
        int64_t s = bps[k], e = bps[k + 1];
        while (kp < kept.size() && kept[kp].second <= s) ++kp;
        if (kp < kept.size() && kept[kp].first < e) continue;
        vector<int> C, act;
        for (int i = 0; i < (int)E.size(); ++i) {
            if (E[i]->qs <= s && E[i]->qe >= e) C.push_back(i);
        }
        if (C.empty()) continue;
        for (int i : C) {
            bool ex = false;
            for (const IV& x : excluded[i]) if (x.first <= s && e <= x.second) ex = true;
            if (!ex) act.push_back(i);
        }
        string reason;
        if (act.empty()) {
            set<string> kinds;
            for (int i : C) for (const auto& w : why[i]) if (w.first.first <= s && e <= w.first.second) kinds.insert(w.second);
            for (const string& x : kinds) reason += (reason.empty() ? "" : "+") + x;
            if (reason.empty()) reason = "excluded";
        } else {
            bool all_dem = true, any_dem = false;
            for (int i : act) { all_dem = all_dem && dem[i]; any_dem = any_dem || dem[i]; }
            reason = (!all_dem && any_dem) ? "tie(tier1)" : (act.size() > 1 ? "tie" : "lone-loser");
        }
        sh.hole_bp[reason] += e - s;
    }
}

void Contig::run() {
    excluded.assign(E.size(), IVs());
    why.assign(E.size(), vector<pair<IV, string>>());
    // which groups are gated and tested
    set<int> gs;
    for (Rec* r : E) {
        if (r->group >= 0) gs.insert(r->group);
    }
    for (int c : gs) {
        if (P.only_chrom.empty() || sh.chroms.names[c] == P.only_chrom) groups.push_back(c);
    }
    std::sort(groups.begin(), groups.end(), [&](int a, int b) { return sh.chroms.names[a] < sh.chroms.names[b]; });
    init_guard();
    // the stock outcome's own structure on each group: an excursion with the same pieces and
    // sidedness is not re-tested
    for (int c : groups) {
        vector<Entry> pl0;
        for (int i = 0; i < (int)E.size(); ++i) {
            Rec& r = *E[i];
            if (r.group != c) continue;
            IVs runs;
            for (const IV& x : r.a0fin) {
                if (!runs.empty() && x.first - runs.back().second <= P.clip) runs.back().second = x.second;
                else runs.push_back(x);
            }
            for (const IV& x : runs) {
                Entry p;
                p.ri = i; p.s = x.first; p.e = x.second; p.pl = place(r, p.s, p.e, sh.chroms); p.li = r.gi; p.rem = false;
                pl0.push_back(p);
            }
        }
        std::stable_sort(pl0.begin(), pl0.end(), [](const Entry& a, const Entry& b) {
                if (a.s != b.s) return a.s < b.s;
                if (a.e != b.e) return a.e < b.e;
                return a.ri < b.ri;
            });
        vector<int> on = backbone(pl0, c);
        set<vector<int64_t>>& sigs = a0sig[c];
        for (const Excursion& ex : excursions(pl0, on, false)) sigs.insert(sig(ex, pl0));
    }
    // demoted: loses to a compatible eligible competitor overlapping -m of its block, so the stock
    // rule would delete it whole.  Fixed for the contig (beats() does not change between rounds)
    dem.assign(E.size(), 0);
    for (int i = 0; i < (int)E.size(); ++i) {
        for (int o = 0; o < (int)E.size(); ++o) {
            if (E[o]->qs >= E[i]->qe) break;
            if (o == i || !compat(E[i], E[o])) continue;
            int64_t ov = overlap(E[i], E[o]);
            if (ov == 0 || !(E[i]->block == 0 || (double)ov / (double)E[i]->block >= P.min_overlap)) continue;
            if (!beats(i, o)) { dem[i] = 1; break; }
        }
    }
    map<int, vector<Entry>> PL;
    map<int, vector<Excursion>> EXS;
    int it = 0;
    while (true) {
        ++it;
        if (it > MAX_ITERATIONS) {
            cerr << "[gaffilter] warning: -x found no fixpoint on " << q << " in " << MAX_ITERATIONS
                 << " rounds: its records get the stock decision" << endl;
            capped = true;
            ++sh.capped;
            for (int i = 0; i < (int)E.size(); ++i) {
                pieces[i].clear();
                if (E[i]->keep0) pieces[i].push_back(make_pair(E[i]->qs, E[i]->qe));
            }
            it = MAX_ITERATIONS;
            break;
        }
        pieces = resolve(true);
        bool changed = false;
        // remainders must be long and similar enough
        for (int i = 0; i < (int)E.size(); ++i) {
            if (!dem[i]) continue;
            Rec& r = *E[i];
            for (const IV& x : pieces[i]) {
                if (x.second - x.first < P.rmin) {
                    excluded[i].push_back(x);
                    why[i].push_back(make_pair(x, string("rmin")));
                    log_row("rmin", x.first, x.second, "rmin", it, r.idx, "", "");
                    changed = true;
                    continue;
                }
                double li = local_ident(r, x.first, x.second);
                if (li < P.rem_floor) {
                    excluded[i].push_back(x);
                    why[i].push_back(make_pair(x, string("floor")));
                    log_row("floor", x.first, x.second, "floor:" + fixed4(li), it, r.idx, "", "");
                    changed = true;
                }
            }
        }
        if (changed) continue;
        PL.clear();
        for (int c : groups) PL[c] = piece_table(c, pieces);
        for (int c : groups) changed = gates(c, PL[c], it) || changed;
        if (changed) continue;
        EXS.clear();
        for (int c : groups) changed = structure(c, PL[c], it, EXS[c]) || changed;
        if (!changed) break;
    }
    iterations = it;
    if (!capped) {
        for (int c : groups) junctions(c, PL[c], EXS[c]);
        for (const auto& u : veto_used) {
            for (const IV& x : veto[u]) {
                stringstream reason;
                reason << "guard:" << E[u.first]->idx << ">" << E[u.second]->idx;
                log_row("guard-applied", x.first, x.second, reason.str(), iterations, E[u.first]->idx, ".", "0");
            }
            ++sh.guard_applied;
        }
    }
    holes();
}

// ---- the stock chain's final anchors (A0), per record

struct A0Line {
    int rec;            // index into the record table
    int64_t qs, qe, qlen, mapq, gl, gm;
    bool primary;
    const string* rc;
    const string* query;
};

// gaf2paf every record the stock GAF stage keeps, cactus's python line filter, then its -p stage
static void compute_a0fin(vector<Rec>& recs, const Params& P, const unordered_map<string, int64_t>& len_map) {
    vector<A0Line> lines;
    string ref_prefix = P.ref_event.empty() ? string() : "id=" + P.ref_event + "|";
    for (int i = 0; i < (int)recs.size(); ++i) {
        Rec& r = recs[i];
        if (!r.keep0) continue;
        GafRecord work;
        const GafRecord* g = r.g;
        // A GAF reused against a graph extended since it was mapped (cactus --inGAF) can have path
        // offsets that reach past its first or last step.  gaf2paf cannot split such a record (it
        // asserts), so cactus drops those steps (trim_unstable_gaf) between gaffilter and gaf2paf.
        // Do the same here, or the split below aborts on exactly the records the reuse path makes.
        if (!g->path.empty() && !(g->path_start < len_map.at(g->path.front().name) &&
                                  g->path_end > g->path_length - len_map.at(g->path.back().name))) {
            work = *g;
            if (!rebase_path(work, len_map)) {
                cerr << "[gaffilter] error: -x: record " << r.idx << " (" << g->query_name << ":" << g->query_start
                     << "-" << g->query_end << ") has path offsets outside its end steps that its node lengths "
                     << "cannot account for" << endl;
                exit(1);
            }
            g = &work;
        }
        if (g->strand == '-') {
            if (g != &work) {
                work = *g;
            }
            gaf2paf_split::flip_gaf(work, len_map);
            g = &work;
        }
        double gi = gaf2paf_split::parent_identity(*g);
        bool is_ref = !ref_prefix.empty() && g->query_name.compare(0, ref_prefix.length(), ref_prefix) == 0;
        gaf2paf_split::for_each_paf_line(*g, len_map, [&](const PafLine& pl, const string& cigar) {
                // filter_paf's line filter (the reference's own lines are exempt)
                double ident = std::min((double)pl.num_matching / ((double)pl.num_bases + 0.00000001), gi);
                bool pass = is_ref || (pl.mapq >= P.min_mapq &&
                                       (pl.query_len <= P.min_block || g->block_length >= P.min_block) &&
                                       ident >= P.min_identity);
                if (pass) {
                    lines.push_back({i, pl.query_start, pl.query_end, pl.query_len, pl.mapq, g->block_length,
                                g->matches, r.primary, &r.rc, &r.g->query_name});
                }
            });
    }
    vector<char> keep(lines.size(), 1);
    if (P.paf_ratio > 0) {
        // gaffilter -p on the lines: each carries its parent's block length, matches, tp and rc
        unordered_map<string, vector<int>> byq;
        for (int i = 0; i < (int)lines.size(); ++i) byq[*lines[i].query].push_back(i);
        for (auto& qv : byq) {
            vector<int>& idx = qv.second;
            std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return lines[a].qs < lines[b].qs; });
            vector<int64_t> starts;
            int64_t maxlen = 0;
            for (int i : idx) {
                starts.push_back(lines[i].qs);
                maxlen = std::max(maxlen, lines[i].qe - lines[i].qs);
            }
            for (int i : idx) {
                const A0Line& a = lines[i];
                size_t lo = std::lower_bound(starts.begin(), starts.end(), a.qs - maxlen) - starts.begin();
                size_t hi = std::lower_bound(starts.begin(), starts.end(), a.qe) - starts.begin();
                for (size_t x = lo; x < hi; ++x) {
                    int j = idx[x];
                    if (j == i) continue;
                    const A0Line& b = lines[j];
                    // closed-interval overlap, as the interval tree tests it
                    if (!(b.qs <= a.qe - 1 && a.qs <= b.qe - 1)) continue;
                    double identity = b.gm ? (double)b.gl / (double)b.gm : 0;
                    if (!(b.mapq >= P.min_mapq && (b.qlen <= P.min_block || b.gl >= P.min_block) &&
                          identity >= P.stock_min_identity)) continue;
                    if (!(*a.rc == *b.rc || a.rc->empty() || b.rc->empty())) continue;
                    int64_t ov = std::min(a.qe, b.qe) - std::max(a.qs, b.qs);
                    if (!(a.gl == 0 || (double)ov / (double)a.gl >= P.paf_min_overlap)) continue;
                    DomClause cl;
                    if (!dominance(a.qs, a.qe, a.primary, (double)a.mapq, (double)a.gl,
                                   b.qs, b.qe, b.primary, (double)b.mapq, (double)b.gl, P.paf_ratio, cl)) {
                        keep[i] = 0;
                        break;
                    }
                }
            }
        }
    }
    vector<IVs> fin(recs.size());
    for (size_t i = 0; i < lines.size(); ++i) {
        if (keep[i]) fin[lines[i].rec].push_back(make_pair(lines[i].qs, lines[i].qe));
    }
    for (size_t i = 0; i < recs.size(); ++i) recs[i].a0fin = merge(fin[i]);
}

// --a0-paf (testing): the stock chain's final PAF, read instead of computed.  A line belongs to the
// one record with its query, gm and gl whose query span holds it
static void read_a0fin(vector<Rec>& recs, const string& path) {
    map<tuple<string, int64_t, int64_t>, vector<int>> bykey;
    for (int i = 0; i < (int)recs.size(); ++i) {
        bykey[make_tuple(recs[i].g->query_name, recs[i].matches, recs[i].block)].push_back(i);
    }
    ifstream in(path);
    if (!in) {
        cerr << "[gaffilter] error: unable to open " << path << endl;
        exit(1);
    }
    vector<IVs> fin(recs.size());
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        PafLine pl = parse_paf_line(line);
        if (!pl.opt_fields.count("gm") || !pl.opt_fields.count("gl")) continue;
        auto it = bykey.find(make_tuple(pl.query_name, (int64_t)stol(pl.opt_fields["gm"].second),
                                        (int64_t)stol(pl.opt_fields["gl"].second)));
        if (it == bykey.end()) continue;
        int found = -1, n = 0;
        for (int i : it->second) {
            if (recs[i].qs <= pl.query_start && pl.query_end <= recs[i].qe) { found = i; ++n; }
        }
        if (n > 1) {
            cerr << "[gaffilter] error: ambiguous parent for --a0-paf line: " << line.substr(0, 200) << endl;
            exit(1);
        }
        if (n == 1) fin[found].push_back(make_pair(pl.query_start, pl.query_end));
    }
    for (size_t i = 0; i < recs.size(); ++i) recs[i].a0fin = merge(fin[i]);
}

static bool write_lines(const string& path, const string& header, const vector<string>& rows) {
    ofstream out(path);
    if (!out) {
        cerr << "[gaffilter] error: unable to open " << path << endl;
        return false;
    }
    if (!header.empty()) out << header << "\n";
    for (const string& r : rows) out << r << "\n";
    out.close();
    if (!out) {
        cerr << "[gaffilter] error: failed to write " << path << endl;
        return false;
    }
    return true;
}

static int run(vector<GafRecord>& gaf_records, unordered_map<string, GafIntervalTree*>& gaf_trees, const Params& P,
               function<string(const GafRecord&)>& print_record) {

    // the stock GAF stage, unchanged, as the baseline
    vector<char> keep0(gaf_records.size(), 1);
    overlap_filter(gaf_records, gaf_trees, P.stock_ratio, P.stock_min_overlap, 0, P.min_mapq, P.min_block,
                   P.stock_min_identity, keep0, print_record);

    // only the nodes the records name.  Read now, after the input: like -l, it is written by the
    // gaf2unstable at the other end of the pipe, in full before its first GAF line
    unordered_set<string> named;
    for (const GafRecord& g : gaf_records) {
        for (const GafStep& st : g.path) named.insert(st.name);
    }
    Chroms chroms;
    unordered_map<string, Node> nodes;
    unordered_map<string, int64_t> len_map;
    {
        ifstream nf(P.nodes_path);
        if (!nf) {
            cerr << "[gaffilter] error: unable to open node table: " << P.nodes_path << endl;
            return 1;
        }
        string line;
        int64_t nline = 0;
        while (getline(nf, line)) {
            ++nline;
            if (line.empty() || line[0] == '#') continue;
            vector<string> toks;
            split_delims(line, "\t", toks);
            if (toks.size() < 2) {
                cerr << "[gaffilter] error: " << P.nodes_path << ":" << nline << ": expected name, length, SN, SO, SR" << endl;
                return 1;
            }
            if (!named.count(toks[0])) continue;
            Node nd;
            nd.len = stol(toks[1]);
            nd.so = toks.size() > 3 && !toks[3].empty() ? stol(toks[3]) : 0;
            bool ref = toks.size() > 4 && toks[4] == "0";
            nd.chrom = ref ? chroms.id(ref_chrom_name(toks[2])) : -1;
            nodes[toks[0]] = nd;
            len_map[toks[0]] = nd.len;
        }
        cerr << "[gaffilter]: Loaded " << nodes.size() << " nodes from " << P.nodes_path << endl;
    }

    // the records
    vector<Rec> recs(gaf_records.size());
    for (size_t i = 0; i < gaf_records.size(); ++i) {
        GafRecord& g = gaf_records[i];
        Rec& r = recs[i];
        r.idx = i;
        r.g = &g;
        if (!g.opt_fields.count("cg")) {
            cerr << "[gaffilter] error: -x needs the cg cigar of every record (minigraph -c)" << endl;
            return 1;
        }
        if (g.opt_fields.count("kq")) {
            cerr << "[gaffilter] error: -x input already carries kq:Z: (the output of -x?)" << endl;
            return 1;
        }
        stringstream pt;
        for (const GafStep& st : g.path) pt << st;
        r.path_text = pt.str();
        r.qlen = g.query_length; r.qs = g.query_start; r.qe = g.query_end;
        r.ps = g.path_start; r.pe = g.path_end;
        r.matches = g.matches; r.block = g.block_length;
        r.mapq = g.mapq == missing_int ? 255 : g.mapq;
        r.strand = g.strand;
        r.primary = is_primary(g);
        r.rc = g.opt_fields.count("rc") ? g.opt_fields.at("rc").second : string();
        r.gi = r.block ? (double)r.matches / (double)r.block : 0.0;
        // eligible to compete and survive: the line filter's own record-level test, with the
        // identity read the right way round (the stock competitor test reads block/matches).  The
        // line filter reads the identity from gaf2paf's gi:f:, rounded to 3 places, so it is
        // rounded here too.  Tested unrounded, a record with matches/block in [i - 0.0005, i) (i the
        // -i value) would pass through uncontested as ineligible, and its lines be kept downstream
        r.eligible = r.mapq >= P.min_mapq && (r.qlen <= P.min_block || r.block >= P.min_block) &&
            gaf2paf_split::parent_identity(g) >= P.min_identity;
        r.keep0 = keep0[i];
        int64_t off = 0, refbp = 0;
        for (const GafStep& st : g.path) {
            auto nt = nodes.find(st.name);
            if (st.is_stable || nt == nodes.end()) {
                cerr << "[gaffilter] error: -x: record " << i << " steps through " << st.name
                     << ", which is not in the node table (the GAF must be gaf2unstable's output, and the table its -n)" << endl;
                return 1;
            }
            r.soff.push_back(off); r.slen.push_back(nt->second.len); r.sso.push_back(nt->second.so);
            r.srev.push_back(st.is_reverse); r.schrom.push_back(nt->second.chrom); r.sref.push_back(refbp);
            if (nt->second.chrom >= 0) refbp += nt->second.len;
            off += nt->second.len;
        }
        if (!r.rc.empty()) {
            r.group = chroms.id(r.rc);
        } else {
            Place pl = place(r, r.qs, r.qe, chroms);
            r.group = pl.ok ? pl.chrom : -1;
        }
    }

    // what the stock chain finally anchors, per record: "new" is judged against it
    if (!P.a0_paf_path.empty()) {
        read_a0fin(recs, P.a0_paf_path);
    } else {
        compute_a0fin(recs, P, len_map);
    }

    // per query, and per haplotype (the query name up to its first |)
    map<string, vector<int>> byq;
    for (int i = 0; i < (int)recs.size(); ++i) byq[recs[i].g->query_name].push_back(i);
    map<string, vector<string>> hap_queries;
    for (const auto& qv : byq) hap_queries[qv.first.substr(0, qv.first.find('|'))].push_back(qv.first);
    // the reference the stock outcome covers with each contig, per chromosome group
    map<pair<string, int>, IVs> hapcov;
    auto hapcov_of = [&](const string& q2, int c) -> const IVs& {
        auto key = make_pair(q2, c);
        auto it = hapcov.find(key);
        if (it != hapcov.end()) return it->second;
        IVs v;
        for (int i : byq[q2]) {
            Rec& r = recs[i];
            if (r.keep0 && r.eligible && (r.rc.empty() || r.rc == chroms.names[c])) {
                IVs x = ref_ivs(r, r.qs, r.qe, c);
                v.insert(v.end(), x.begin(), x.end());
            }
        }
        return hapcov[key] = merge(v);
    };

    vector<string> log_rows, junction_rows;
    Shared sh(P, chroms, log_rows, junction_rows);
    int64_t duplicates = 0;
    // per record: kept pieces (whole, deleted, or cut) and status
    vector<int> plan_kind(recs.size(), 0);         // 0 whole, 1 cut, 2 deleted
    vector<IVs> plan_pieces(recs.size());
    vector<string> status(recs.size());
    string ref_prefix = P.ref_event.empty() ? string() : "id=" + P.ref_event + "|";
    for (auto& qv : byq) {
        const string& q = qv.first;
        if (!ref_prefix.empty() && q.compare(0, ref_prefix.length(), ref_prefix) == 0) {
            // the reference's own records pass under the stock rule
            for (int i : qv.second) plan_kind[i] = recs[i].keep0 ? 0 : 2;
            continue;
        }
        vector<Rec*> E;
        for (int i : qv.second) {
            if (recs[i].eligible) {
                E.push_back(&recs[i]);
            } else {
                // passes through whole: the line filter downstream drops it
                plan_kind[i] = 0;
                status[i] = "ineligible";
            }
        }
        std::sort(E.begin(), E.end(), key_less);
        for (size_t i = 1; i < E.size(); ++i) {
            if (key_equal(E[i - 1], E[i])) {
                // not fatal: the stock filter takes such input in its stride (the copies delete each
                // other), and failing a whole mapping job over it would cost more than it tells
                cerr << "[gaffilter] warning: -x: two records of " << q << " are identical (input records "
                     << E[i - 1]->idx << " and " << E[i]->idx << ")" << endl;
                ++duplicates;
            }
        }
        string hap = q.substr(0, q.find('|'));
        auto hco_fn = [&, q, hap](int c) {
            IVs v;
            for (const string& q2 : hap_queries[hap]) {
                if (q2 == q) continue;
                const IVs& x = hapcov_of(q2, c);
                v.insert(v.end(), x.begin(), x.end());
            }
            return merge(v);
        };
        Contig contig(sh, q, E, hco_fn);
        contig.pieces.assign(E.size(), IVs());
        contig.run();
        for (size_t k = 0; k < E.size(); ++k) {
            int i = E[k]->idx;
            const IVs& pv = contig.pieces[k];
            status[i] = contig.dem[k] ? "demoted" : "ok";
            if (pv.size() == 1 && pv[0].first == recs[i].qs && pv[0].second == recs[i].qe) {
                plan_kind[i] = 0;
            } else if (pv.empty()) {
                plan_kind[i] = 2;
            } else {
                plan_kind[i] = 1;
                plan_pieces[i] = pv;
            }
        }
    }

    // the GAF: whole records as they are, cut ones whole with the kept query intervals in kq:Z:
    int64_t n_whole = 0, n_cut = 0, n_del = 0, kept_bp = 0, stock_bp = 0;
    for (size_t i = 0; i < recs.size(); ++i) {
        if (recs[i].keep0) stock_bp += recs[i].qe - recs[i].qs;
        if (plan_kind[i] == 2) {
            ++n_del;
            continue;
        }
        if (plan_kind[i] == 0) {
            ++n_whole;
            kept_bp += recs[i].qe - recs[i].qs;
            cout << gaf_records[i] << "\n";
            continue;
        }
        ++n_cut;
        stringstream kq;
        for (size_t k = 0; k < plan_pieces[i].size(); ++k) {
            kq << (k ? "," : "") << plan_pieces[i][k].first << "-" << plan_pieces[i][k].second;
            kept_bp += plan_pieces[i][k].second - plan_pieces[i][k].first;
        }
        gaf_records[i].opt_fields["kq"] = make_pair("Z", kq.str());
        cout << gaf_records[i] << "\n";
        gaf_records[i].opt_fields.erase("kq");
    }

    bool ok = true;
    if (!P.plan_path.empty()) {
        vector<string> rows;
        for (size_t i = 0; i < recs.size(); ++i) {
            const Rec& r = recs[i];
            stringstream row;
            int64_t kb = plan_kind[i] == 0 ? r.qe - r.qs : ivlen(plan_pieces[i]);
            row << r.idx << "\t" << r.g->query_name << "\t" << r.qs << "\t" << r.qe << "\t" << r.block << "\t" << r.mapq
                << "\t" << status[i] << "\t" << kb << "\t";
            if (plan_kind[i] == 0) row << "whole";
            else if (plan_kind[i] == 2) row << "deleted";
            else for (size_t k = 0; k < plan_pieces[i].size(); ++k) row << (k ? "," : "") << plan_pieces[i][k].first << "-" << plan_pieces[i][k].second;
            rows.push_back(row.str());
        }
        ok = write_lines(P.plan_path, "#idx\tquery\tqs\tqe\tblock\tmapq\tstatus\tkept_bp\tpieces", rows) && ok;
    }
    if (!P.log_path.empty()) {
        ok = write_lines(P.log_path, "", log_rows) && ok;
    }
    if (!P.junctions_path.empty()) {
        ok = write_lines(P.junctions_path, "", junction_rows) && ok;
    }

    // the summary: on stderr, and in --exact-summary for a caller whose stderr is shared with a pipe
    vector<string> summary;
    {
        stringstream ss;
        ss << "-x: " << n_whole << " records whole, " << n_cut << " cut, " << n_del << " deleted, of "
           << recs.size() << ". kept " << kept_bp << " query bp, against " << stock_bp << " for the stock filter ("
           << (kept_bp >= stock_bp ? "+" : "") << (kept_bp - stock_bp) << ")";
        summary.push_back(ss.str());
    }
    {
        stringstream ss;
        ss << "-x: holes (bp claimed by an eligible record, kept by none):";
        if (sh.hole_bp.empty()) ss << " none";
        for (const auto& h : sh.hole_bp) ss << " " << h.first << "=" << h.second;
        summary.push_back(ss.str());
    }
    {
        stringstream ss;
        ss << "-x: " << sh.isolations << " junction trims (isolations), " << sh.gate_drops
           << " gate drops (r1ref/collide), " << junction_rows.size() << " junctions logged";
        if (P.guard) ss << ", guard vetoes " << sh.guard_vetoes << " computed, " << sh.guard_applied << " applied";
        summary.push_back(ss.str());
    }
    if (duplicates) {
        stringstream ss;
        ss << "-x: warning: " << duplicates << " record(s) duplicate another record of their contig";
        summary.push_back(ss.str());
    }
    if (sh.capped) {
        stringstream ss;
        ss << "-x: warning: " << sh.capped << " contig(s) found no fixpoint and got the stock decision";
        summary.push_back(ss.str());
    }
    for (const string& line : summary) {
        cerr << "[gaffilter]: " << line << endl;
    }
    if (!P.summary_path.empty()) {
        ok = write_lines(P.summary_path, "", summary) && ok;
    }
    return ok ? 0 : 1;
}

} // namespace exact


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
         << "    -p, --paf                       Input is PAF, not GAF" << endl
         << endl
         << "exact mode (GAF input only; not with -p or -o; replaces the removed trim mode -t and its -e, -g, -l, -Q, -R):" << endl
         << "    -x, --exact                     Resolve each query segment on its own instead of deleting whole records: a record keeps the segments" << endl
         << "                                    where it beats (as -r does) every compatible record also claiming them.  A record the stock rule deletes" << endl
         << "                                    (it loses to a competitor overlapping -m of its block) only fills segments nobody else claims, and only" << endl
         << "                                    in long, similar enough pieces.  New sequence (what the stock filter chain would not anchor) is gated" << endl
         << "                                    against reference the haplotype already covers, and off-backbone excursions that cannot be real" << endl
         << "                                    two-sided rearrangements have their junctions broken by a query gap.  -q/-b/-i select the records" << endl
         << "                                    that take part (-i as matches/block).  A cut record is printed whole with the query intervals it keeps" << endl
         << "                                    in kq:Z:, which gaf2paf applies" << endl
         << "    --exact-nodes FILE              Node table written by gaf2unstable -n (required with -x)" << endl
         << "    --exact-ratio X                 Dominance ratio for -x's contests and demotion [-r]" << endl
         << "    --exact-ref NAME                Records of queries id=NAME|* get the stock rule, and are exempt from the line filter [none]" << endl
         << "    --a0-paf-ratio X                Ratio of the stock -p stage that follows (PAFOverlapFilterRatio; 0 = none) [0]" << endl
         << "    --a0-paf-min-overlap X          -m of that -p stage [0]" << endl
         << "    --min-remainder N               Shortest piece a demoted record may keep [10000]" << endl
         << "    --min-remainder-ident X         Lowest identity of a piece a demoted record keeps [0.95]" << endl
         << "    --id-floor X                    Lowest identity of a new piece (path identity, and over its reference-node columns) [0.98]" << endl
         << "    --skip-r1ref                    Do not test new pieces' identity over their reference-node columns" << endl
         << "    --skip-onetoone                 Do not drop new sequence placed on reference its haplotype already covers" << endl
         << "    --gate-chunk N                  Chunk size of the one-to-one test [10000]" << endl
         << "    --gate-minrun N                 Drop failing chunks in runs of at least this [20000]" << endl
         << "    --gap N                         Query gap that breaks a junction; a junction touching a new piece joins under it, and a" << endl
         << "                                    two-sided excursion with a query gap this wide inside it has its junctions broken [31000]" << endl
         << "    --junction-clip N               Query gap under which a junction between stock pieces joins (cactus's clip) [10000]" << endl
         << "    --bound-tol N                   A bounded excursion's reference footprint lies within the backbone window +- this [100000]" << endl
         << "    --bound-slack N                 ...and the window exceeds the footprint by at most this [1000000]" << endl
         << "    --lin-maxindel N                Largest implied indel between colinear backbone pieces [25000000]" << endl
         << "    --lin-ovtol N                   Reference step-back allowed between colinear backbone pieces [50000]" << endl
         << "    --junction-min N                Reference jump that makes a logged junction a jump [50000]" << endl
         << "    --guard                         Identity guard: where --exact-ratio < -r decides a contest -r would tie, cancel the win over runs" << endl
         << "                                    where the loser is the more similar (--guard-trip 0.01 --guard-minrun 20000 --guard-chunk 10000" << endl
         << "                                    --guard-minchunk 2000 --guard-svlen 50)" << endl
         << "    --exact-plan FILE               Write the per-record plan (kept query pieces) to FILE" << endl
         << "    --exact-log FILE                Write every drop, isolation and guard decision to FILE" << endl
         << "    --junctions FILE                Write the junctions between pieces placed apart (review log) to FILE" << endl
         << "    --exact-summary FILE            Also write the summary printed on stderr to FILE" << endl
         << "    --a0-paf FILE                   (testing) Read the stock chain's final PAF instead of computing it" << endl
         << "    --only-chrom NAME               (testing) Gate and test only this chromosome group" << endl;
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
    // the trim mode (-t and its -e, -g, -l, -Q, -R) was removed; -x replaces it
    string removed_option;

    bool exact_mode = false;
    exact::Params xp;
    string ratio_arg, min_overlap_arg, min_identity_arg, exact_ratio_arg;
    // long-only options of -x.  None starts with a letter an existing long option starts with, so
    // every abbreviation that resolved before still resolves the same way
    enum { X_NODES = 256, X_RATIO, X_REF, X_PAF_RATIO, X_PAF_MIN_OVERLAP, X_RMIN, X_REM_FLOOR, X_ID_FLOOR,
           X_NO_R1REF, X_NO_ONETOONE, X_GATE_CHUNK, X_GATE_MINRUN, X_GAP, X_CLIP, X_R2_TOL, X_R2_SLACK,
           X_LIN_MAXINDEL, X_LIN_OVTOL, X_JUNCTION_MIN, X_GUARD, X_GUARD_TRIP, X_GUARD_MINRUN, X_GUARD_CHUNK,
           X_GUARD_MINCHUNK, X_GUARD_SVLEN, X_PLAN, X_LOG, X_JUNCTIONS, X_SUMMARY, X_A0_PAF, X_ONLY_CHROM };

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
            // the removed trim mode's options, still recognised so they fail with a pointer to -x
            // rather than as unknown options
            {"trim", no_argument, 0, 't'},
            {"trim-min-gap", required_argument, 0, 'g'},
            {"close-holes", no_argument, 0, 'R'},
            {"trim-edge", required_argument, 0, 'e'},
            {"trim-min-mapq", required_argument, 0, 'Q'},
            {"node-lengths", required_argument, 0, 'l'},
            {"paf", no_argument, 0, 'p'},
            {"exact", no_argument, 0, 'x'},
            {"exact-nodes", required_argument, 0, X_NODES},
            {"exact-ratio", required_argument, 0, X_RATIO},
            {"exact-ref", required_argument, 0, X_REF},
            {"a0-paf-ratio", required_argument, 0, X_PAF_RATIO},
            {"a0-paf-min-overlap", required_argument, 0, X_PAF_MIN_OVERLAP},
            {"min-remainder", required_argument, 0, X_RMIN},
            {"min-remainder-ident", required_argument, 0, X_REM_FLOOR},
            {"id-floor", required_argument, 0, X_ID_FLOOR},
            {"skip-r1ref", no_argument, 0, X_NO_R1REF},
            {"skip-onetoone", no_argument, 0, X_NO_ONETOONE},
            {"gate-chunk", required_argument, 0, X_GATE_CHUNK},
            {"gate-minrun", required_argument, 0, X_GATE_MINRUN},
            {"gap", required_argument, 0, X_GAP},
            {"junction-clip", required_argument, 0, X_CLIP},
            {"bound-tol", required_argument, 0, X_R2_TOL},
            {"bound-slack", required_argument, 0, X_R2_SLACK},
            {"lin-maxindel", required_argument, 0, X_LIN_MAXINDEL},
            {"lin-ovtol", required_argument, 0, X_LIN_OVTOL},
            {"junction-min", required_argument, 0, X_JUNCTION_MIN},
            {"guard", no_argument, 0, X_GUARD},
            {"guard-trip", required_argument, 0, X_GUARD_TRIP},
            {"guard-minrun", required_argument, 0, X_GUARD_MINRUN},
            {"guard-chunk", required_argument, 0, X_GUARD_CHUNK},
            {"guard-minchunk", required_argument, 0, X_GUARD_MINCHUNK},
            {"guard-svlen", required_argument, 0, X_GUARD_SVLEN},
            {"exact-plan", required_argument, 0, X_PLAN},
            {"exact-log", required_argument, 0, X_LOG},
            {"junctions", required_argument, 0, X_JUNCTIONS},
            {"exact-summary", required_argument, 0, X_SUMMARY},
            {"a0-paf", required_argument, 0, X_A0_PAF},
            {"only-chrom", required_argument, 0, X_ONLY_CHROM},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:r:m:po:b:q:i:tg:l:Re:Q:x",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            ratio = stof(optarg);
            ratio_arg = optarg;
            break;
        case 'm':
            min_overlap_pct = stof(optarg);
            min_overlap_arg = optarg;
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
            min_identity_arg = optarg;
            break;
        case 'q':
            min_mapq = std::stol(optarg);
            break;
        case 't':
        case 'g':
        case 'l':
        case 'R':
        case 'e':
        case 'Q':
            if (removed_option.empty()) {
                removed_option = long_options[option_index].val == c ?
                    string("--") + long_options[option_index].name : string("-") + (char)c;
            }
            break;
        case 'x':
            exact_mode = true;
            break;
        // -x's thresholds are parsed as doubles (the emulator it reproduces uses Python floats)
        case X_NODES: xp.nodes_path = optarg; break;
        case X_RATIO: exact_ratio_arg = optarg; break;
        case X_REF: xp.ref_event = optarg; break;
        // these two are the stock -p stage's -r and -m, which it parses as floats
        case X_PAF_RATIO: xp.paf_ratio = stof(optarg); break;
        case X_PAF_MIN_OVERLAP: xp.paf_min_overlap = stof(optarg); break;
        case X_RMIN: xp.rmin = stol(optarg); break;
        case X_REM_FLOOR: xp.rem_floor = stod(optarg); break;
        case X_ID_FLOOR: xp.id_floor = stod(optarg); break;
        case X_NO_R1REF: xp.r1ref = false; break;
        case X_NO_ONETOONE: xp.onetoone = false; break;
        case X_GATE_CHUNK: xp.gate_chunk = stol(optarg); break;
        case X_GATE_MINRUN: xp.gate_minrun = stol(optarg); break;
        case X_GAP: xp.gap = stol(optarg); break;
        case X_CLIP: xp.clip = stol(optarg); break;
        case X_R2_TOL: xp.r2_tol = stol(optarg); break;
        case X_R2_SLACK: xp.r2_slack = stol(optarg); break;
        case X_LIN_MAXINDEL: xp.lin_maxindel = stol(optarg); break;
        case X_LIN_OVTOL: xp.lin_ovtol = stol(optarg); break;
        case X_JUNCTION_MIN: xp.junction_min = stol(optarg); break;
        case X_GUARD: xp.guard = true; break;
        case X_GUARD_TRIP: xp.guard_trip = stod(optarg); break;
        case X_GUARD_MINRUN: xp.guard_minrun = stol(optarg); break;
        case X_GUARD_CHUNK: xp.guard_chunk = stol(optarg); break;
        case X_GUARD_MINCHUNK: xp.guard_minchunk = stol(optarg); break;
        case X_GUARD_SVLEN: xp.guard_svlen = stol(optarg); break;
        case X_PLAN: xp.plan_path = optarg; break;
        case X_LOG: xp.log_path = optarg; break;
        case X_JUNCTIONS: xp.junctions_path = optarg; break;
        case X_SUMMARY: xp.summary_path = optarg; break;
        case X_A0_PAF: xp.a0_paf_path = optarg; break;
        case X_ONLY_CHROM: xp.only_chrom = optarg; break;
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

    if (!removed_option.empty()) {
        cerr << "[gaffilter] error: " << removed_option << ": the trim mode (-t/--trim, with -e, -g, -l, -Q and -R) "
             << "has been removed. Use -x/--exact instead, which cuts a record that loses an overlap down to "
             << "exactly the contested segment (see the exact mode options)" << endl;
        return 1;
    }

    if (exact_mode) {
        if (is_paf || min_overlap_len) {
            cerr << "[gaffilter] error: -x/--exact cannot be used with -p or -o" << endl;
            return 1;
        }
        if (ratio == 0) {
            cerr << "[gaffilter] error: -x/--exact needs -r" << endl;
            return 1;
        }
        if (xp.nodes_path.empty()) {
            cerr << "[gaffilter] error: -x/--exact needs --exact-nodes (gaf2unstable -n)" << endl;
            return 1;
        }
        xp.stock_ratio = ratio;
        xp.min_mapq = min_mapq;
        xp.min_block = min_block_len;
        xp.stock_min_overlap = min_overlap_pct;
        xp.stock_min_identity = min_identity;
        xp.stock_ratio_d = stod(ratio_arg);
        xp.ratio = exact_ratio_arg.empty() ? xp.stock_ratio_d : stod(exact_ratio_arg);
        xp.min_overlap = min_overlap_arg.empty() ? 0. : stod(min_overlap_arg);
        xp.min_identity = min_identity_arg.empty() ? 0. : stod(min_identity_arg);
        // a ratio of 0 lets every record beat every other (double placement); a chunk size of 0
        // divides by zero (--gate-chunk) or never advances (--guard-chunk)
        if (!(xp.ratio > 0) || xp.gate_chunk <= 0 || xp.guard_chunk <= 0) {
            cerr << "[gaffilter] error: -x needs --exact-ratio > 0, --gate-chunk > 0 and --guard-chunk > 0" << endl;
            return 1;
        }
    } else if (!exact_ratio_arg.empty() || !xp.nodes_path.empty() || xp.guard || !xp.plan_path.empty()) {
        cerr << "[gaffilter] error: --exact-* options and --guard need -x" << endl;
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

    // make an interval tree for each query sequence
    unordered_map<string, GafIntervalTree*> gaf_trees = build_query_trees(gaf_records);
    cerr << "[gaffilter]: Constructed interval trees" << endl;

    if (exact_mode) {
        int ret = exact::run(gaf_records, gaf_trees, xp, print_record);
        for (auto qt : gaf_trees) {
            delete qt.second;
        }
        return ret;
    }


    int64_t filter_count = 0;
    int64_t filter_len_count = 0;

    vector<char> keep(gaf_records.size(), 1);
    overlap_filter(gaf_records, gaf_trees, ratio, min_overlap_pct, min_overlap_len, min_mapq, min_block_len,
                   min_identity, keep, print_record);

    for (int64_t i = 0; i < (int64_t)gaf_records.size(); ++i) {
        if (keep[i]) {
            cout << print_record(gaf_records[i]) << "\n";
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
    return 0;
}
