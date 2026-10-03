/*
  zip_alt.cpp -- spec step 5: alt-vs-alt (v2).  See zip_alt.hpp for the definitions.

  Per site:
    1. residual excursions (kinds F, I, BK; no node used by the reference pass), with their walks
       dep, alts, arr; an index handle -> (representative candidate, position) over the handles of
       creator excursions whose node occurs once in the walk;
    2. per round: every member candidate M (query >= b) meets every representative candidate X of
       lower rank that shares >= 2 handles with it (through the index); consecutive shared handles
       (u, v) bound a branch of M (non-empty, >= b, no node on X, no node used); it is a member of
       group (u, v) in this round if X is the group's representative of the round (the next distinct,
       unused branch of at least ceil(b*i) in the group's order).  Identical member branches merge
       (their excursions' owners);
    3. per candidate: the sizes (no row if they fail), pair-too-big, shares-nodes, not-private,
       in-series, the site cap, then one Aligner::align_jobs call for the round, and the confident
       chains go to SitePlanner::accept_alt in a total order;
    4. near-parallel diagnostics, aligned in one more align_jobs call.
  Every order is total and nothing depends on the thread count.
*/
#include "zip_alt.hpp"

#include <chrono>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace zip {

void AltStats::add(const AltStats& o) {
    residual += o.residual;
    pairs += o.pairs;
    groups += o.groups;
    candidates += o.candidates;
    shares_nodes += o.shares_nodes;
    not_private += o.not_private;
    in_series += o.in_series;
    pair_too_big += o.pair_too_big;
    capped += o.capped;
    aligned += o.aligned;
    confident += o.confident;
    zipped += o.zipped;
    refused += o.refused;
    reverted += o.reverted;
    zipped_bp += o.zipped_bp;
    near_parallel += o.near_parallel;
}

namespace {

// near-parallel: window starts and ends within this distance
const int64_t NEAR_BP = 200;

// a residual excursion
struct RExc {
    const Unit* unit = nullptr;
    std::vector<Handle> walk;           // dep, alts, arr (canonical orientation)
    std::vector<int64_t> pre;           // pre[k]: bp of walk[0..k)
    std::vector<uint8_t> once;          // the node of walk[k] occurs once in the walk
    std::vector<NodeId> nodes;          // the walk's nodes, sorted, unique
    int32_t rank = INT32_MAX;           // the highest SR of its alt nodes (for a creator excursion: its run's SR)
    int32_t rank_rep = INT32_MAX;       // = rank for a creator excursion; INT32_MAX: never a representative
    const Owner* rep_owner = nullptr;   // the lowest creator owner (tie-breaks)
    const Owner* first_owner = nullptr; // the lowest owner
    bool has(NodeId n) const { return std::binary_search(nodes.begin(), nodes.end(), n); }
    int64_t bp(uint32_t a, uint32_t b) const { return pre[b] - pre[a]; }   // walk[a..b)
};

// a representative candidate of a group (u, v): excursion x, with u at pu and v at pv
struct RepCand {
    uint32_t x = NONE;
    uint32_t pu = 0, pv = 0;
    int64_t bp = 0;                     // its branch walk[pu+1..pv)
};

struct Group {
    bool listed = false;
    std::vector<RepCand> order;         // representative order: rank, then the longest branch, SN, SO, name
    size_t next = 0;                    // the next round starts looking here
    int round = 0;                      // the round `cur` belongs to
    size_t cur = SIZE_MAX;              // this round's representative (index into order), SIZE_MAX: none
};

// a candidate: a member branch at (u, v) against the round's representative
struct CKey {
    Handle u = 0, v = 0;
    std::vector<Handle> mb;
    bool operator<(const CKey& o) const {
        if (u != o.u) return u < o.u;
        if (v != o.v) return v < o.v;
        return mb < o.mb;
    }
};
struct ACand {
    Handle u = 0, v = 0;
    std::vector<Handle> mb;             // member branch, walk orientation
    uint32_t x = NONE;                  // representative excursion
    uint32_t pu = 0, pv = 0;            // the representative branch is rex[x].walk[pu+1..pv)
    std::vector<uint32_t> ms;           // member excursions that walk this branch at (u, v)
};

// a branch of member M relative to representative candidate X (round-independent): consecutive shared
// handles at M[k1], M[k2] and X[p1], X[p2]
struct BEntry {
    uint32_t x = NONE, m = NONE;
    uint32_t k1 = 0, k2 = 0, p1 = 0, p2 = 0;
};

inline int64_t iabs64(int64_t x) { return x < 0 ? -x : x; }

std::string handles_str(const Graph& g, const std::vector<Handle>& hs) {
    std::string s;
    for (size_t i = 0; i < hs.size(); ++i) {
        if (i) s.push_back(',');
        s += g.handle_str(hs[i]);
    }
    return s;
}

// a branch for keys and the report: every handle, or the first three and the last two
std::string branch_str(const Graph& g, const std::vector<Handle>& b) {
    std::string s;
    const size_t n = b.size();
    for (size_t i = 0; i < n; ++i) {
        if (n > 8 && i == 3) {
            s += strf(",..%zu more..", n - 5);
            i = n - 3;
            continue;
        }
        if (i) s.push_back(',');
        s += g.handle_str(b[i]);
    }
    return s;
}

class AltSite {
public:
    AltSite(const Graph& g, const SiteData& sd, SitePlanner& planner, const Options& opt, AltStats& st)
        : g_(g), sd_(sd), s_(*sd.site), planner_(planner), opt_(opt), st_(st), label_(s_.label(g)) {}

    void run(Aligner& aligner, std::deque<ReportRow>& rows);

private:
    const Graph& g_;
    const SiteData& sd_;
    const Site& s_;
    SitePlanner& planner_;
    const Options& opt_;
    AltStats& st_;
    const std::string label_;
    std::vector<RExc> rex;
    std::unordered_map<Handle, std::vector<std::pair<uint32_t, uint32_t>>> index;   // creator walks: handle -> (x, pos), by (rank, x)
    std::vector<BEntry> entries;                                                      // every branch of every pair, once per site
    bool scanned = false;
    std::unordered_map<NodeId, std::vector<uint32_t>> node_units;                     // alt node -> units walking it
    std::map<std::pair<Handle, Handle>, Group> groups;
    int64_t budget = 0;
    std::string dump;
    // scratch
    std::vector<uint32_t> hstamp, ustamp;
    uint32_t hround = 0, uround = 0;

    void collect();
    bool used_any(const std::vector<Handle>& hs) const {
        for (Handle h : hs)
            if (planner_.node_used(handle_node(h))) return true;
        return false;
    }
    bool rep_less(const RepCand& a, const RepCand& b) const;
    uint32_t round_rep(Handle u, Handle v, int round);
    void pair_branches(uint32_t xi, uint32_t mi, const std::vector<std::pair<uint32_t, uint32_t>>& C);
    void scan_pairs();
    std::map<CKey, ACand> round_candidates(int round);
    bool shares_nodes(const std::vector<Handle>& mb, const std::vector<Handle>& rb);
    bool is_private(const std::vector<Handle>& branch, Handle u, Handle v) const;
    bool in_series(const std::vector<Handle>& mb, const std::vector<Handle>& rb);
    ReportRow base_row() const;
    // the planner key of the site's alt row k (AltCandidate::key, ZippedPiece::unit, the edit dump's
    // unit column): past the site's unit ids, so it never equals a reference-pass unit id
    size_t key_of(size_t k) const { return sd_.units.size() + k; }
    void near_parallel(Aligner& aligner, std::deque<ReportRow>& rows);
};

void AltSite::collect() {
    // residual excursions: kinds F, I, BK; no node (anchors included) zipped or used as a target
    for (const Unit& u : sd_.units) {
        if (u.kind == Kind::J) continue;
        if (planner_.node_used(handle_node(u.exc.dep)) || planner_.node_used(handle_node(u.exc.arr)) || used_any(u.exc.alts)) continue;
        RExc r;
        r.unit = &u;
        r.walk.reserve(u.exc.alts.size() + 2);
        r.walk.push_back(u.exc.dep);
        r.walk.insert(r.walk.end(), u.exc.alts.begin(), u.exc.alts.end());
        r.walk.push_back(u.exc.arr);
        r.pre.assign(r.walk.size() + 1, 0);
        for (size_t k = 0; k < r.walk.size(); ++k) r.pre[k + 1] = r.pre[k] + g_.len(handle_node(r.walk[k]));
        for (Handle h : r.walk) r.nodes.push_back(handle_node(h));
        std::sort(r.nodes.begin(), r.nodes.end());
        std::vector<NodeId> twice;
        for (size_t k = 1; k < r.nodes.size(); ++k)
            if (r.nodes[k] == r.nodes[k - 1]) twice.push_back(r.nodes[k]);
        r.nodes.erase(std::unique(r.nodes.begin(), r.nodes.end()), r.nodes.end());
        r.once.assign(r.walk.size(), 1);
        for (size_t k = 0; k < r.walk.size(); ++k)
            if (std::binary_search(twice.begin(), twice.end(), handle_node(r.walk[k]))) r.once[k] = 0;
        // rank: when the path came into existence, the highest SR of its alt nodes.  For a creator
        // excursion that is its run's SR (reused nodes are older); for an observed (GAF) excursion it
        // is the insertion, not the contigs that happen to walk it
        r.rank = 0;
        for (Handle h : u.exc.alts) r.rank = std::max(r.rank, g_.nodes[handle_node(h)].sr);
        for (const Owner& o : u.exc.owners) {          // owners are sorted by owner_less (SR first)
            if (o.sr < 0) continue;
            if (!r.first_owner) r.first_owner = &o;
            if (o.src == Src::CREATOR && !r.rep_owner) r.rep_owner = &o;
        }
        if (!r.first_owner && !u.exc.owners.empty()) r.first_owner = &u.exc.owners[0];
        if (r.rep_owner) r.rank_rep = r.rank;
        rex.push_back(std::move(r));
    }
    st_.residual += rex.size();
    // the index over representative candidates (x in increasing order, so every list is sorted by x)
    for (uint32_t x = 0; x < (uint32_t)rex.size(); ++x) {
        if (!rex[x].rep_owner) continue;
        for (uint32_t k = 0; k < (uint32_t)rex[x].walk.size(); ++k)
            if (rex[x].once[k]) index[rex[x].walk[k]].push_back(std::make_pair(x, k));
    }
    // each list by (rank, x): a member scans only the prefix of lower-rank candidates
    for (auto& kv : index)
        std::sort(kv.second.begin(), kv.second.end(), [&](const std::pair<uint32_t, uint32_t>& a, const std::pair<uint32_t, uint32_t>& b) {
            if (rex[a.first].rank_rep != rex[b.first].rank_rep) return rex[a.first].rank_rep < rex[b.first].rank_rep;
            return a.first < b.first;
        });
    // every unit of the site (residual or not) through each alt node: the disjointness check
    for (const Unit& u : sd_.units) {
        NodeId last = NONE;
        for (Handle h : u.exc.alts) {
            NodeId n = handle_node(h);
            if (n == last) continue;
            std::vector<uint32_t>& v = node_units[n];
            if (v.empty() || v.back() != u.id) v.push_back(u.id);
            last = n;
        }
    }
    ustamp.assign(sd_.units.size(), 0);
    hstamp.assign(sd_.sg.n_handles(), 0);
}

bool AltSite::rep_less(const RepCand& a, const RepCand& b) const {
    const RExc& X = rex[a.x];
    const RExc& Y = rex[b.x];
    if (X.rank_rep != Y.rank_rep) return X.rank_rep < Y.rank_rep;
    if (a.bp != b.bp) return a.bp > b.bp;
    const Owner& ox = *X.rep_owner;
    const Owner& oy = *Y.rep_owner;
    if (ox.contig != oy.contig) return ox.contig < oy.contig;
    if (ox.so != oy.so) return ox.so < oy.so;
    if (ox.first != oy.first) return ox.first < oy.first;
    return X.unit->id < Y.unit->id;
}

// The representative of group (u, v) in this round: the next branch in the group's order (distinct
// branches of at least ceil(b*i), by rank) whose nodes (and the bounds) are unused.  Members that did
// not merge regroup under it.
uint32_t AltSite::round_rep(Handle u, Handle v, int round) {
    Group& G = groups[std::make_pair(u, v)];
    if (!G.listed) {
        G.listed = true;
        auto iu = index.find(u), iv = index.find(v);
        if (iu != index.end() && iv != index.end()) {
            std::unordered_map<uint32_t, uint32_t> pu;
            for (const auto& e : iu->second) pu[e.first] = e.second;
            for (const auto& e : iv->second) {
                auto it = pu.find(e.first);
                if (it == pu.end() || it->second >= e.second) continue;
                RepCand rc;
                rc.x = e.first;
                rc.pu = it->second;
                rc.pv = e.second;
                rc.bp = rex[rc.x].bp(rc.pu + 1, rc.pv);
                G.order.push_back(rc);
            }
            std::sort(G.order.begin(), G.order.end(), [&](const RepCand& a, const RepCand& b) { return rep_less(a, b); });
            // one entry per distinct branch (several excursions can share it): a later round takes
            // the next branch, not the same branch of another excursion.  A branch shorter than
            // ceil(b*i) can never be the target of a unit, so it is no round's representative.
            std::set<std::vector<Handle>> seen;
            std::vector<RepCand> uniq;
            for (const RepCand& rc : G.order) {
                const RExc& X = rex[rc.x];
                if (rc.bp < opt_.min_window()) continue;
                if (seen.insert(std::vector<Handle>(X.walk.begin() + rc.pu + 1, X.walk.begin() + rc.pv)).second) uniq.push_back(rc);
            }
            G.order.swap(uniq);
        }
    }
    if (G.round != round) {
        G.round = round;
        G.cur = SIZE_MAX;
        if (!planner_.node_used(handle_node(u)) && !planner_.node_used(handle_node(v))) {
            for (size_t i = G.next; i < G.order.size(); ++i) {
                const RepCand& rc = G.order[i];
                const RExc& X = rex[rc.x];
                bool intact = true;
                for (uint32_t k = rc.pu + 1; k < rc.pv && intact; ++k)
                    if (planner_.node_used(handle_node(X.walk[k]))) intact = false;
                if (intact) { G.cur = i; break; }
            }
        }
        G.next = G.cur == SIZE_MAX ? G.order.size() : G.cur + 1;
    }
    return G.cur == SIZE_MAX ? NONE : G.order[G.cur].x;
}

// The branches of M relative to X, independent of the round: consecutive shared handles (u, v) in the
// same order in both walks, with a non-empty branch of M of >= b bp, none of its nodes on X.  C holds
// (position in M, position in X) of their shared handles, in M's order.
void AltSite::pair_branches(uint32_t xi, uint32_t mi, const std::vector<std::pair<uint32_t, uint32_t>>& C) {
    const RExc& X = rex[xi];
    const RExc& M = rex[mi];
    for (size_t j = 0; j + 1 < C.size(); ++j) {
        const uint32_t k1 = C[j].first, k2 = C[j + 1].first, p1 = C[j].second, p2 = C[j + 1].second;
        if (k2 <= k1 + 1 || p2 <= p1) continue;            // empty, or not in the same order
        if (M.bp(k1 + 1, k2) < opt_.b) continue;
        bool on = false;
        for (uint32_t k = k1 + 1; k < k2 && !on; ++k)
            if (X.has(handle_node(M.walk[k]))) on = true;   // none of them on R
        if (on) continue;
        BEntry e;
        e.x = xi;
        e.m = mi;
        e.k1 = k1;
        e.k2 = k2;
        e.p1 = p1;
        e.p2 = p2;
        entries.push_back(e);
    }
}

// Every (representative candidate, member) pair sharing >= 2 handles, once per site: members in
// unit order, candidates in index order, so the entries (and every round's candidates) come out in
// one total order.
void AltSite::scan_pairs() {
    scanned = true;
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> hits(rex.size());
    std::vector<uint32_t> touched;
    for (uint32_t mi = 0; mi < (uint32_t)rex.size(); ++mi) {
        const RExc& M = rex[mi];
        if (M.unit->query_bp < opt_.b) continue;
        touched.clear();
        for (uint32_t k = 0; k < (uint32_t)M.walk.size(); ++k) {
            if (!M.once[k]) continue;
            auto it = index.find(M.walk[k]);
            if (it == index.end()) continue;
            for (const auto& xp : it->second) {
                if (!(rex[xp.first].rank_rep < M.rank)) break;     // lower rank only: a prefix of the list
                if (xp.first == mi) continue;
                std::vector<std::pair<uint32_t, uint32_t>>& hv = hits[xp.first];
                if (hv.empty()) touched.push_back(xp.first);
                hv.push_back(std::make_pair(k, xp.second));
            }
        }
        std::sort(touched.begin(), touched.end());
        for (uint32_t x : touched) {
            std::vector<std::pair<uint32_t, uint32_t>>& C = hits[x];
            if (C.size() >= 2) {
                ++st_.pairs;
                pair_branches(x, mi, C);
            }
            C.clear();
        }
    }
}

// This round's candidates: the branches whose nodes are unused, at a group whose representative of
// the round is their X.  Identical member branches merge.
std::map<CKey, ACand> AltSite::round_candidates(int round) {
    if (!scanned) scan_pairs();
    std::map<CKey, ACand> out;
    for (const BEntry& e : entries) {
        const RExc& M = rex[e.m];
        bool used = false;
        for (uint32_t k = e.k1 + 1; k < e.k2 && !used; ++k)
            if (planner_.node_used(handle_node(M.walk[k]))) used = true;   // merged in an earlier round
        if (used) continue;
        const Handle u = M.walk[e.k1], v = M.walk[e.k2];
        if (round_rep(u, v, round) != e.x) continue;
        CKey key;
        key.u = u;
        key.v = v;
        key.mb.assign(M.walk.begin() + e.k1 + 1, M.walk.begin() + e.k2);
        ACand& c = out[key];
        if (c.ms.empty()) {
            c.u = u;
            c.v = v;
            c.mb = key.mb;
            c.x = e.x;
            c.pu = e.p1;
            c.pv = e.p2;
        }
        c.ms.push_back(e.m);
    }
    return out;
}

// Disjoint: no excursion of the site walks both a member-branch node and a representative-branch
// node (so no image can reach the representative twice); and nothing is used yet.
bool AltSite::shares_nodes(const std::vector<Handle>& mb, const std::vector<Handle>& rb) {
    if (used_any(mb) || used_any(rb)) return true;
    std::unordered_set<NodeId> rset;
    for (Handle h : rb) rset.insert(handle_node(h));
    for (Handle h : mb)
        if (rset.count(handle_node(h))) return true;
    if (++uround == 0) { std::fill(ustamp.begin(), ustamp.end(), 0); uround = 1; }
    for (Handle h : mb) {
        auto it = node_units.find(handle_node(h));
        if (it == node_units.end()) continue;
        for (uint32_t ui : it->second) {
            if (ustamp[ui] == uround) continue;
            ustamp[ui] = uround;
            for (Handle w : sd_.units[ui].exc.alts)
                if (rset.count(handle_node(w))) return true;
        }
    }
    return false;
}

// Private: a side-aware flood from the branch that never crosses u or v reaches no reference node,
// and touches u only on u's outgoing side and v only on v's incoming side.
bool AltSite::is_private(const std::vector<Handle>& branch, Handle u, Handle v) const {
    const NodeId nu = handle_node(u), nv = handle_node(v);
    const Side out_u = exit_side(u), in_v = entry_side(v);
    std::unordered_set<NodeId> seen;
    std::vector<NodeId> stk;
    for (Handle h : branch)
        if (seen.insert(handle_node(h)).second) stk.push_back(handle_node(h));
    while (!stk.empty()) {
        NodeId x = stk.back();
        stk.pop_back();
        for (int o = 0; o < 2; ++o)
            for (const Edge& e : g_.out(make_handle(x, o == 1))) {
                const NodeId y = handle_node(e.to);
                const Side t = entry_side(e.to);
                if (y == nu) {
                    if (t != out_u) return false;
                    continue;
                }
                if (y == nv) {
                    if (t != in_v) return false;
                    continue;
                }
                if (g_.is_ref(y) || !s_.is_alt(y)) return false;
                if (seen.insert(y).second) stk.push_back(y);
            }
    }
    return true;
}

// In series: an alt-only path joins the two branches, in either orientation; u and v may be crossed.
// Walks from both orientations of every member-branch node over the site's alt handles.
bool AltSite::in_series(const std::vector<Handle>& mb, const std::vector<Handle>& rb) {
    const SiteGraph& sg = sd_.sg;
    std::unordered_set<uint32_t> target;
    for (Handle h : rb) {
        uint32_t l = sg.local_node(handle_node(h));
        if (l != NONE) target.insert(l);
    }
    if (++hround == 0) { std::fill(hstamp.begin(), hstamp.end(), 0); hround = 1; }
    std::vector<uint32_t> stk;
    for (Handle h : mb) {
        uint32_t l = sg.local_node(handle_node(h));
        if (l == NONE) continue;
        for (uint32_t lh : {2 * l, 2 * l + 1})
            if (hstamp[lh] != hround) { hstamp[lh] = hround; stk.push_back(lh); }
    }
    while (!stk.empty()) {
        uint32_t h = stk.back();
        stk.pop_back();
        if (target.count(h >> 1)) return true;
        for (uint32_t k = sg.succ_off[h]; k < sg.succ_off[h + 1]; ++k) {
            uint32_t w = sg.succ[k];
            if (hstamp[w] != hround) { hstamp[w] = hround; stk.push_back(w); }
        }
    }
    return false;
}

ReportRow AltSite::base_row() const {
    ReportRow r;
    r.site = label_;
    r.sn = g_.sn_names[(size_t)s_.sn];
    r.lo = s_.lo;
    r.hi = s_.hi;
    return r;
}

struct Pending {                        // an eligible candidate waiting for the cap and the aligner
    size_t row = 0;
    const ACand* c = nullptr;
    std::vector<Handle> rb;
    int64_t mbp = 0;
    uint32_t order = 0;
};

void AltSite::run(Aligner& aligner, std::deque<ReportRow>& rows) {
    const auto t0 = std::chrono::steady_clock::now();
    AlignTimingScope timing;            // this pass's aligner calls (they also count for the site)
    const size_t rows0 = rows.size();
    const uint64_t pairs0 = st_.pairs, aligned0 = st_.aligned;
    collect();
    budget = opt_.max_site_query;
    const int64_t minw = opt_.min_window();
    if (rex.size() >= 2) {
        for (int round = 1; round <= opt_.alt_rounds; ++round) {
            std::map<CKey, ACand> cands = round_candidates(round);
            if (cands.empty()) break;
            std::vector<Pending> pend;
            std::set<std::pair<Handle, Handle>> grp;
            uint32_t order = 0;
            for (const auto& kv : cands) {
                const ACand& c = kv.second;
                const RExc& X = rex[c.x];
                int64_t mbp = 0;
                for (Handle h : c.mb) mbp += g_.len(handle_node(h));
                const int64_t tbp = X.bp(c.pu + 1, c.pv);
                if (mbp < opt_.b || tbp < minw) continue;   // not a unit: no row
                std::vector<Handle> rb(X.walk.begin() + c.pu + 1, X.walk.begin() + c.pv);
                // the row
                ReportRow row = base_row();
                row.pass = "alt";
                row.round = round;
                std::vector<Owner> ow;
                Src src = Src::WITNESS;
                for (uint32_t m : c.ms) {
                    const Excursion& e = rex[m].unit->exc;
                    ow.insert(ow.end(), e.owners.begin(), e.owners.end());
                    if ((int)e.src < (int)src) src = e.src;
                }
                std::sort(ow.begin(), ow.end(), owner_less);
                ow.erase(std::unique(ow.begin(), ow.end(), owner_equal), ow.end());
                row.source = src_name(src);
                row.owners = owners_str(g_, ow);
                row.kind = kind_name(rex[c.ms[0]].unit->kind);
                row.anchors = g_.handle_str(c.u) + ">" + g_.handle_str(c.v);
                row.window = strf("0-%lld", (long long)tbp);
                row.query_bp = mbp;
                row.representative = owners_str(g_, std::vector<Owner>(1, *X.rep_owner)) + " " + branch_str(g_, rb);
                row.unit = rex[c.ms[0]].unit->id;
                ++st_.candidates;
                grp.insert(std::make_pair(c.u, c.v));
                std::string why;
                if (mbp > opt_.max_pair || tbp > opt_.max_pair) why = "pair-too-big";
                else if (shares_nodes(c.mb, rb)) why = "shares-nodes";
                else if (!is_private(c.mb, c.u, c.v)) why = "not-private";
                else if (in_series(c.mb, rb)) why = "in-series";
                row.outcome = why.empty() ? "candidate" : why;
                if (!opt_.dump_dir.empty())
                    dump += strf("cand\t%zu\t%d\t%s\t%s\t%s\t%lld\t%s\t%s\t%lld\t%s\t%s\n", key_of(rows.size()), round, g_.handle_str(c.u).c_str(),
                                 g_.handle_str(c.v).c_str(),
                                 row.owners.c_str(), (long long)mbp, handles_str(g_, c.mb).c_str(),
                                 owners_str(g_, std::vector<Owner>(1, *X.rep_owner)).c_str(), (long long)tbp, handles_str(g_, rb).c_str(),
                                 row.outcome.c_str());
                rows.push_back(std::move(row));
                if (!why.empty()) continue;
                Pending p;
                p.row = rows.size() - 1;
                p.c = &c;
                p.rb = std::move(rb);
                p.mbp = mbp;
                p.order = order++;
                pend.push_back(std::move(p));
            }
            st_.groups += grp.size();
            // the site cap, by member bp (descending), then candidate order
            {
                std::vector<size_t> by(pend.size());
                for (size_t i = 0; i < by.size(); ++i) by[i] = i;
                std::sort(by.begin(), by.end(), [&](size_t a, size_t b) {
                    if (pend[a].mbp != pend[b].mbp) return pend[a].mbp > pend[b].mbp;
                    return a < b;
                });
                std::vector<uint8_t> capped(pend.size(), 0);
                bool over = false;
                for (size_t i : by) {
                    if (over || pend[i].mbp > budget) { over = true; capped[i] = 1; continue; }
                    budget -= pend[i].mbp;
                }
                std::vector<Pending> keep;
                for (size_t i = 0; i < pend.size(); ++i) {
                    if (capped[i]) rows[pend[i].row].outcome = "site:capped";
                    else keep.push_back(std::move(pend[i]));
                }
                pend.swap(keep);
            }
            if (pend.empty()) continue;
            // align the round's members against their representatives (alt targets in walk orientation)
            std::vector<AlignJob> jobs(pend.size());
            for (size_t i = 0; i < pend.size(); ++i) {
                const ACand& c = *pend[i].c;
                AlignJob& j = jobs[i];
                j.key = strf("%s:alt:%d:%s>%s>%s", label_.c_str(), round, g_.handle_str(c.u).c_str(), branch_str(g_, c.mb).c_str(),
                             g_.handle_str(c.v).c_str());
                j.query = c.mb;
                j.target.ref = false;
                j.target.walk = pend[i].rb;
                // every target base is a representative node: no feasibility cut
            }
            std::vector<AlignResult> res;
            aligner.align_jobs(g_, jobs, res);
            st_.aligned += jobs.size();
            std::vector<size_t> conf;
            for (size_t i = 0; i < pend.size(); ++i) {
                fill_align_row(rows[pend[i].row], res[i], opt_.G);
                if (res[i].outcome.empty() && !res[i].chain.parts.empty()) conf.push_back(i);
            }
            st_.confident += conf.size();
            // acceptance order: support (complete walks), aligned bp, sum '=', candidate order
            auto weight = [&](size_t i) {
                uint64_t w = 0;
                for (uint32_t m : pend[i].c->ms) w += rex[m].unit->exc.weight;
                return w;
            };
            std::sort(conf.begin(), conf.end(), [&](size_t a, size_t b) {
                if (sd_.complete_walks) {
                    uint64_t wa = weight(a), wb = weight(b);
                    if (wa != wb) return wa > wb;
                }
                const ChainResult& x = res[a].chain;
                const ChainResult& y = res[b].chain;
                if (x.aligned_q != y.aligned_q) return x.aligned_q > y.aligned_q;
                if (x.score != y.score) return x.score > y.score;
                return pend[a].order < pend[b].order;
            });
            std::vector<AltCandidate> acs;
            for (size_t i : conf) {
                AltCandidate a;
                a.query = pend[i].c->mb;
                a.u = pend[i].c->u;
                a.v = pend[i].c->v;
                a.target = pend[i].rb;
                a.res = &res[i];
                a.row = &rows[pend[i].row];
                a.key = (uint32_t)key_of(pend[i].row);
                a.round = round;
                acs.push_back(std::move(a));
            }
            planner_.accept_alt(acs);
        }
    }
    near_parallel(aligner, rows);
    // wall time less the waits for an aligner slot (-j), which other sites' work decides
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double wait = timing.timing().waiting();
    if (secs - wait > 5 || log_level() > 1)
        ZLOG("site %s: alt-vs-alt %.1f s, of which %.1f s waiting for an aligner slot (%zu residual excursion(s), %llu pair(s), "
             "%zu branch(es), %zu row(s), %llu aligned; minimap2 running %.1f s, CPU %.1f s)",
             label_.c_str(), secs, wait, rex.size(), (unsigned long long)(st_.pairs - pairs0), entries.size(), rows.size() - rows0,
             (unsigned long long)(st_.aligned - aligned0), timing.timing().running(), timing.timing().mm2_cpu);
    if (!opt_.dump_dir.empty() && !dump.empty()) {
        std::string dir = opt_.dump_dir + "/alt";
        if (!make_dirs(dir)) fail(EXIT_IO, "cannot create " + dir);
        AtomicWriter w(dir + "/" + label_ + ".tsv");
        w.write("#type\tkey\tround\tu\tv\tmember_owners\tmember_bp\tmember_branch\trepresentative\ttarget_bp\trepresentative_branch\t"
                "check\n");
        w.write(dump);
        w.commit();
    }
}

// Residual excursions with different anchors whose windows lie within NEAR_BP of each other's: one
// row per other anchor pair of a cluster, against the cluster's representative.  Never zipped.
void AltSite::near_parallel(Aligner& aligner, std::deque<ReportRow>& rows) {
    std::vector<uint32_t> big;
    for (uint32_t i = 0; i < (uint32_t)rex.size(); ++i)
        if (rex[i].unit->query_bp >= opt_.b) big.push_back(i);
    if (big.size() < 2) return;
    std::sort(big.begin(), big.end(), [&](uint32_t a, uint32_t b) {
        if (rex[a].unit->wlo != rex[b].unit->wlo) return rex[a].unit->wlo < rex[b].unit->wlo;
        return a < b;
    });
    std::vector<uint32_t> parent(big.size());
    for (uint32_t i = 0; i < (uint32_t)big.size(); ++i) parent[i] = i;
    auto find = [&](uint32_t i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    };
    for (size_t i = 0; i < big.size(); ++i)
        for (size_t j = i + 1; j < big.size() && rex[big[j]].unit->wlo - rex[big[i]].unit->wlo <= NEAR_BP; ++j) {
            const Unit& a = *rex[big[i]].unit;
            const Unit& b = *rex[big[j]].unit;
            if (iabs64(a.whi - b.whi) > NEAR_BP) continue;
            uint32_t x = find((uint32_t)i), y = find((uint32_t)j);
            if (x != y) parent[std::max(x, y)] = std::min(x, y);
        }
    std::map<uint32_t, std::vector<uint32_t>> clusters;   // root -> rex indices
    for (uint32_t i = 0; i < (uint32_t)big.size(); ++i) clusters[find(i)].push_back(big[i]);
    // excursion order: rank (creator rank for representatives), then the longest, then SN, SO, first node, unit
    auto exc_less = [&](uint32_t a, uint32_t b, bool rep) {
        const RExc& X = rex[a];
        const RExc& Y = rex[b];
        int32_t rx = rep ? X.rank_rep : X.rank, ry = rep ? Y.rank_rep : Y.rank;
        if (rx != ry) return rx < ry;
        if (X.unit->query_bp != Y.unit->query_bp) return X.unit->query_bp > Y.unit->query_bp;
        const Owner* ox = rep ? X.rep_owner : X.first_owner;
        const Owner* oy = rep ? Y.rep_owner : Y.first_owner;
        if (ox && oy) {
            if (ox->contig != oy->contig) return ox->contig < oy->contig;
            if (ox->so != oy->so) return ox->so < oy->so;
            if (ox->first != oy->first) return ox->first < oy->first;
        }
        return X.unit->id < Y.unit->id;
    };
    struct NearPair { uint32_t rep, mem; };
    std::vector<NearPair> pairs;
    for (const auto& kv : clusters) {
        const std::vector<uint32_t>& cl = kv.second;
        std::map<std::pair<Handle, Handle>, std::vector<uint32_t>> sub;
        for (uint32_t i : cl) sub[std::make_pair(rex[i].walk.front(), rex[i].walk.back())].push_back(i);
        if (sub.size() < 2) continue;
        uint32_t rep = NONE;
        for (uint32_t i : cl)
            if (rex[i].rep_owner && (rep == NONE || exc_less(i, rep, true))) rep = i;
        if (rep == NONE) continue;
        const std::pair<Handle, Handle> rk(rex[rep].walk.front(), rex[rep].walk.back());
        for (const auto& sv : sub) {
            if (sv.first == rk) continue;
            uint32_t m = sv.second[0];
            for (uint32_t i : sv.second)
                if (exc_less(i, m, false)) m = i;
            bool share = false;
            for (NodeId n : rex[m].nodes)
                if (!g_.is_ref(n) && rex[rep].has(n)) share = true;
            if (share) continue;
            pairs.push_back(NearPair{rep, m});
        }
    }
    if (pairs.empty()) return;
    std::sort(pairs.begin(), pairs.end(), [&](const NearPair& a, const NearPair& b) {
        if (rex[a.mem].unit->id != rex[b.mem].unit->id) return rex[a.mem].unit->id < rex[b.mem].unit->id;
        return rex[a.rep].unit->id < rex[b.rep].unit->id;
    });
    std::vector<AlignJob> jobs;
    std::vector<size_t> jrow;
    for (const NearPair& p : pairs) {
        const Unit& M = *rex[p.mem].unit;
        const Unit& R = *rex[p.rep].unit;
        ReportRow row = base_row();
        row.pass = "diag";
        row.round = 0;
        row.source = src_name(M.exc.src);
        row.owners = owners_str(g_, M.exc.owners);
        row.kind = kind_name(M.kind);
        row.anchors = g_.handle_str(M.exc.dep) + ">" + g_.handle_str(M.exc.arr);
        row.window = strf("%lld-%lld", (long long)M.wlo, (long long)M.whi);
        row.query_bp = M.query_bp;
        row.unit = M.id;
        const bool overlap = std::max(M.wlo, R.wlo) <= std::min(M.whi, R.whi);
        const bool priv = is_private(M.exc.alts, M.exc.dep, M.exc.arr);
        row.outcome = std::string("near-parallel:") + (overlap ? "overlap" : "disjoint") + ":" + (priv ? "private" : "entangled");
        row.representative = owners_str(g_, std::vector<Owner>(1, *rex[p.rep].rep_owner)) + " " + g_.handle_str(R.exc.dep) + ">" +
                             branch_str(g_, R.exc.alts) + ">" + g_.handle_str(R.exc.arr);
        if (!opt_.dump_dir.empty())
            dump += strf("near\t%zu\t0\t%s\t%s\t%s\t%lld\t%s\t%s\t%lld\t%s\t%s\n", key_of(rows.size()), g_.handle_str(M.exc.dep).c_str(),
                         g_.handle_str(M.exc.arr).c_str(),
                         row.owners.c_str(), (long long)M.query_bp, handles_str(g_, M.exc.alts).c_str(),
                         owners_str(g_, std::vector<Owner>(1, *rex[p.rep].rep_owner)).c_str(), (long long)R.query_bp,
                         excursion_str(g_, R.exc).c_str(), row.outcome.c_str());
        rows.push_back(std::move(row));
        // the alignment of the two alt runs (a diagnostic; within the site's alt-vs-alt cap)
        if (M.query_bp > opt_.max_pair || R.query_bp > opt_.max_pair || M.query_bp > budget) continue;
        budget -= M.query_bp;
        AlignJob j;
        j.key = strf("%s:near:%u", label_.c_str(), M.id);
        j.query = M.exc.alts;
        j.target.ref = false;
        j.target.walk = R.exc.alts;
        jobs.push_back(std::move(j));
        jrow.push_back(rows.size() - 1);
    }
    if (jobs.empty()) return;
    std::vector<AlignResult> res;
    aligner.align_jobs(g_, jobs, res);
    for (size_t i = 0; i < jobs.size(); ++i) {
        ReportRow& row = rows[jrow[i]];
        std::string keep = row.outcome;
        fill_align_row(row, res[i], opt_.G);
        row.outcome = keep;
    }
}

} // namespace

void alt_pass(const Graph& g, const SiteData& sd, Aligner* aligner, SitePlanner& planner, const Options& opt, std::deque<ReportRow>& rows,
              AltStats& st) {
    if (!aligner) return;
    AltSite a(g, sd, planner, opt, st);
    a.run(*aligner, rows);
}

void alt_tally(const std::deque<ReportRow>& rows, AltStats& st) {
    for (const ReportRow& r : rows) {
        if (r.pass == "diag") {
            if (r.outcome.compare(0, 14, "near-parallel:") == 0) ++st.near_parallel;
            continue;
        }
        if (r.pass != "alt") continue;
        const std::string& o = r.outcome;
        if (o == "zipped") {
            ++st.zipped;
            if (r.kept_bp > 0) st.zipped_bp += r.kept_bp;
        } else if (o == "(D)" || o == "trimmed-below-b") {
            ++st.refused;
        } else if (o.compare(0, 9, "reverted:") == 0) {
            ++st.reverted;
        } else if (o == "shares-nodes") {
            ++st.shares_nodes;
        } else if (o == "not-private") {
            ++st.not_private;
        } else if (o == "in-series") {
            ++st.in_series;
        } else if (o == "pair-too-big") {
            ++st.pair_too_big;
        } else if (o == "site:capped") {
            ++st.capped;
        }
    }
}

} // namespace zip
