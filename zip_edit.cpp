/*
  zip_edit.cpp -- spec step 6: pieces, the left-snap projection, snapping, unaligned stretches,
  acceptance (A)-(D), cuts, the splitter, the pinch, the transaction with validators V1-V7,
  --check, id allocation and applying the plans.

  How a site is edited
    Every decision is made on the input graph.  An accepted piece gives its node interval a fate:
    "zipped onto [ta, tb) of target space S with rel r".  The edited site is then a pure function
    of (input site, fates, block-end cut points) and is built by build_model():
      1. ends: every site link end, and both sides of every fate end inside a node, is an End
         (node, offset, side); an End on a zipped piece is translated (the pinch):
           rel '+': p.L -> L side of the target piece starting at ta;  p.R -> R side of the piece ending at tb
           rel '-': p.L -> R side of the target piece ending at tb;    p.R -> L side of the piece starting at ta
         read in the space's walk orientation and converted to node sides;
      2. a translated link equal to the target's own adjacency inside a node is dropped (no cut);
      3. cuts: fate ends; every remaining link end inside a node (the pieces that carry a
         translated link); and P of every block end of the accepted chains;
      4. the splitter: pieces inherit SN/SR, SO = parent SO + offset; links on a node's L side go
         to its first piece, on its R side to its last (side-aware, both ends of a self-loop);
         consecutive kept pieces of a node are joined by a new ++ link carrying the node's SR;
      5. links are deduplicated under a canonical key (a link equals its reverse) with priority:
         input links, then new split links, then translated links; a translated image equal to the
         target's own adjacency or to an existing link is dropped;
      6. zipped pieces are deleted.  Nothing is wired to chain ends or snarl boundaries; an end
         that resolves to no live piece is an internal error.
    The validators read the model: V1 reference walk, V2 every unit's image exists link by link and
    spells the predicted sequence (and, with --check, every anchored path at sites with <= 5,000),
    V3 no piece twice in an image, V4 no new forward-strand SCC, V5 locality / windows / rank-0 back
    edges, V6 every kept piece reaches the reference, V7 pieces tile their parents (SN/SR/SO) and,
    per chain over its kept pieces, removed <= target / i + G.

  Acceptance (spec 6.7) is greedy in a total order; a piece failing (A), (B) or (C) is dropped and
  stays alt (it is then attached at its projected points by the translation); a chain is refused
  by trimmed-below-b or (D).  With observed walks (v3) a new piece also fails (G), the whole-record
  rule, when a GAF record that walks its node reads its target elsewhere (WalkSource::reads_target;
  the GAF audit's reading rule).  The transaction (6.8) validates the whole edit and, on failure,
  re-plans incrementally in acceptance order, dropping each chain whose addition fails.
  A chain's report kept_bp is the bp of its new pieces (compatible pieces are counted by the chain
  that zipped them first): the zipped rows sum to the bp the edit removes (checked).

  P (the left-snap projection) is zip_align's project_left_snap, shared with the feasibility cut.
*/
#include "zip_edit.hpp"

#include <chrono>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace zip {

namespace {

inline int64_t iabs64(int64_t x) { return x < 0 ? -x : x; }

// Debug hook for test Z29 (not a CLI option): RGFA_ZIP_DEBUG_DROP_LINK="<site label>:<unit id>[,...]"
// leaves out the first translated link (canonical order) that the named unit's chain creates,
// whenever that chain is in the edit.  The transaction must then revert that chain alone.
const std::vector<std::pair<std::string, uint32_t>>& debug_drop_links() {
    static std::once_flag once;
    static std::vector<std::pair<std::string, uint32_t>> v;
    std::call_once(once, []() {
        const char* e = getenv("RGFA_ZIP_DEBUG_DROP_LINK");
        if (!e || !*e) return;
        std::string s(e);
        size_t p = 0;
        while (p <= s.size()) {
            size_t q = s.find(',', p);
            if (q == std::string::npos) q = s.size();
            std::string item = s.substr(p, q - p);
            size_t c = item.rfind(':');
            int64_t u = 0;
            if (c == std::string::npos || !parse_int64(item.substr(c + 1), u) || u < 0)
                fail(EXIT_INPUT, "RGFA_ZIP_DEBUG_DROP_LINK: expected <site label>:<unit id>[,...], got '" + item + "'");
            v.push_back(std::make_pair(item.substr(0, c), (uint32_t)u));
            ZLOG("debug: RGFA_ZIP_DEBUG_DROP_LINK drops one translated link of unit %u at site %s", (uint32_t)u, item.substr(0, c).c_str());
            p = q + 1;
        }
    });
    return v;
}

// ================================================================ target spaces

// An oriented walk of nodes carrying target coordinates.  Element k (handle walk[k]) covers
// [c[k], c[k+1]).  The first and last elements are sentinels (A and B for the reference; u and v
// for an alt branch): targets lie in [c[1], c[n-1]].
struct TSpace {
    bool ref = true;
    std::vector<Handle> walk;
    std::vector<int64_t> c;
    int64_t lo() const { return c[1]; }
    int64_t hi() const { return c[c.size() - 2]; }
    // element k with c[k] < T <= c[k+1]
    size_t elem_ending(int64_t T) const {
        return (size_t)(std::lower_bound(c.begin() + 1, c.end(), T) - c.begin()) - 1;
    }
    // element k with c[k] <= T < c[k+1]
    size_t elem_starting(int64_t T) const {
        return (size_t)(std::upper_bound(c.begin(), c.end(), T) - c.begin()) - 1;
    }
};

// ================================================================ the site context

struct Ctx {
    const Graph& g;
    const Site& s;
    const SiteData& sd;
    const Options& opt;
    std::vector<NodeId> nodes;                    // site nodes (A, B, interior), sorted
    std::unordered_map<NodeId, uint32_t> index;   // node -> position in nodes
    std::vector<uint32_t> links;                  // input links with an end on a site side, sorted
    std::vector<uint8_t> is_ref, orient_rev;      // per site node; orient_rev: forward-strand orientation '-'
    uint32_t iA = NONE, iB = NONE;
    std::vector<TSpace> spaces;                   // [0] = the reference

    Ctx(const Graph& g_, const SiteData& sd_, const Options& opt_) : g(g_), s(*sd_.site), sd(sd_), opt(opt_) {}

    uint32_t idx(NodeId n) const {
        auto it = index.find(n);
        return it == index.end() ? NONE : it->second;
    }
    int64_t len(uint32_t i) const { return g.len(nodes[i]); }
};

void build_ctx(Ctx& cx) {
    const Graph& g = cx.g;
    const Site& s = cx.s;
    cx.nodes = s.interior;
    cx.nodes.push_back(s.A);
    cx.nodes.push_back(s.B);
    std::sort(cx.nodes.begin(), cx.nodes.end());
    cx.nodes.erase(std::unique(cx.nodes.begin(), cx.nodes.end()), cx.nodes.end());
    cx.index.reserve(cx.nodes.size() * 2);
    for (uint32_t i = 0; i < (uint32_t)cx.nodes.size(); ++i) cx.index[cx.nodes[i]] = i;
    cx.iA = cx.idx(s.A);
    cx.iB = cx.idx(s.B);
    // site links: every link of an interior node, plus the links on A.R and B.L
    for (NodeId n : s.interior)
        for (int o = 0; o < 2; ++o)
            for (const Edge& e : g.out(make_handle(n, o == 1))) cx.links.push_back(e.link);
    for (const Edge& e : g.out(make_handle(s.A, false))) cx.links.push_back(e.link);   // A.R
    for (const Edge& e : g.out(make_handle(s.B, true))) cx.links.push_back(e.link);    // B.L
    std::sort(cx.links.begin(), cx.links.end());
    cx.links.erase(std::unique(cx.links.begin(), cx.links.end()), cx.links.end());
    cx.is_ref.assign(cx.nodes.size(), 0);
    cx.orient_rev.assign(cx.nodes.size(), 0);
    for (uint32_t i = 0; i < (uint32_t)cx.nodes.size(); ++i) cx.is_ref[i] = g.is_ref(cx.nodes[i]) ? 1 : 0;
    // forward-strand orientation of alt nodes: as the node's creator excursion walks it (canonical
    // frame: anchors '+'), else as the first unit through it walks it, else '+'
    {
        std::vector<uint8_t> set(cx.nodes.size(), 0);
        for (int pass = 0; pass < 2; ++pass)
            for (const Unit& u : cx.sd.units) {
                for (Handle h : u.exc.alts) {
                    uint32_t i = cx.idx(handle_node(h));
                    if (i == NONE || set[i] || cx.is_ref[i]) continue;
                    if (pass == 0) {
                        bool own = false;
                        const std::string& sn = g.sn_names[g.nodes[handle_node(h)].sn];
                        for (const Owner& o : u.exc.owners)
                            if (o.contig == sn) own = true;
                        if (!own) continue;
                    }
                    cx.orient_rev[i] = handle_rev(h) ? 1 : 0;
                    set[i] = 1;
                }
            }
    }
    // the reference space: A+, backbone, B+
    TSpace S;
    S.ref = true;
    S.walk.push_back(make_handle(s.A, false));
    S.c.push_back(g.start(s.A));
    for (NodeId b : s.backbone) {
        S.walk.push_back(make_handle(b, false));
        S.c.push_back(g.start(b));
    }
    S.walk.push_back(make_handle(s.B, false));
    S.c.push_back(g.start(s.B));
    S.c.push_back(g.end(s.B));
    cx.spaces.push_back(std::move(S));
}

// ================================================================ fates, ends

// A zipped node interval: [a, b) of a node lands on [ta, tb) of target space `space` with rel.
struct Fate {
    int64_t a = 0, b = 0;
    uint32_t space = 0;
    int64_t ta = 0, tb = 0;
    char rel = '+';
    uint32_t chain = NONE;       // candidate index
};
typedef std::map<uint32_t, std::vector<Fate>> FateMap;   // site node index -> fates sorted by a

// An end of a link: the R side of the piece of node i ending at x (right), or the L side of the
// piece starting at x.
struct End {
    uint32_t i = NONE;
    int64_t x = 0;
    bool right = false;
};

End side_end(const Ctx& cx, Side sd) {
    uint32_t i = cx.idx(side_node(sd));
    End e;
    e.i = i;
    e.right = side_is_right(sd);
    e.x = e.right ? (i == NONE ? 0 : cx.len(i)) : 0;
    return e;
}

// walk-R side at T: the R side (in walk orientation) of the space piece ending at T
End walk_R(const Ctx& cx, const TSpace& S, int64_t T) {
    size_t k = S.elem_ending(T);
    Handle h = S.walk[k];
    End e;
    e.i = cx.idx(handle_node(h));
    int64_t d = T - S.c[k];
    if (!handle_rev(h)) { e.x = d; e.right = true; }
    else { e.x = cx.len(e.i) - d; e.right = false; }
    return e;
}

// walk-L side at T: the L side (in walk orientation) of the space piece starting at T
End walk_L(const Ctx& cx, const TSpace& S, int64_t T) {
    size_t k = S.elem_starting(T);
    Handle h = S.walk[k];
    End e;
    e.i = cx.idx(handle_node(h));
    int64_t d = T - S.c[k];
    if (!handle_rev(h)) { e.x = d; e.right = false; }
    else { e.x = cx.len(e.i) - d; e.right = true; }
    return e;
}

// the pinch: where a side of the zipped piece [f.a, f.b) goes
End translate(const Ctx& cx, const Fate& f, bool right) {
    const TSpace& S = cx.spaces[f.space];
    if (!right) return f.rel == '+' ? walk_L(cx, S, f.ta) : walk_R(cx, S, f.tb);
    return f.rel == '+' ? walk_R(cx, S, f.tb) : walk_L(cx, S, f.ta);
}

const Fate* fate_ending(const std::vector<Fate>& v, int64_t x) {
    for (const Fate& f : v)
        if (f.b == x) return &f;
    return nullptr;
}
const Fate* fate_starting(const std::vector<Fate>& v, int64_t x) {
    for (const Fate& f : v)
        if (f.a == x) return &f;
    return nullptr;
}

// ================================================================ the edited site (model)

struct MPiece {
    uint32_t i = NONE;           // site node index
    int64_t off = 0, len = 0;
    bool kept = true;
    uint32_t chain = NONE;       // zipped: the chain whose fate it is
};

struct MLink {
    uint32_t sa = 0, sb = 0;     // piece sides 2*piece + right; for class 0, sa is the input's first end
    int32_t sr = -1;
    uint32_t src = NONE;         // input link, or NONE (split link)
    uint32_t split_i = NONE;     // split link: site node index and offset
    int64_t split_x = 0;
    uint8_t cls = 0;             // 0 input, 1 new split (kept-kept), 2 translated
    uint32_t chain = NONE;       // translated: smallest chain among its translated ends
};

struct Model {
    std::vector<uint32_t> first;                 // per site node: first piece (size N + 1)
    std::vector<MPiece> pieces;
    std::vector<MLink> links;
    std::unordered_set<uint64_t> keys;
    std::vector<uint8_t> changed;                // per site node: cut or zipped
    FateMap fates;
    std::string error;                           // internal inconsistency (never expected)

    uint32_t npieces(uint32_t i) const { return first[i + 1] - first[i]; }
    uint32_t piece_ending(uint32_t i, int64_t x) const {
        uint32_t a = first[i], b = first[i + 1];
        while (a < b) {                          // first piece with end >= x
            uint32_t m = (a + b) / 2;
            if (pieces[m].off + pieces[m].len < x) a = m + 1; else b = m;
        }
        if (a < first[i + 1] && pieces[a].off + pieces[a].len == x) return a;
        return NONE;
    }
    uint32_t piece_starting(uint32_t i, int64_t x) const {
        uint32_t a = first[i], b = first[i + 1];
        while (a < b) {
            uint32_t m = (a + b) / 2;
            if (pieces[m].off < x) a = m + 1; else b = m;
        }
        if (a < first[i + 1] && pieces[a].off == x) return a;
        return NONE;
    }
    static uint64_t key(uint32_t sa, uint32_t sb) {
        uint32_t lo = std::min(sa, sb), hi = std::max(sa, sb);
        return ((uint64_t)lo << 32) | hi;
    }
    bool has_link(uint32_t sa, uint32_t sb) const { return keys.count(key(sa, sb)) != 0; }
};

struct RawLink {
    End a, b;
    int32_t sr = -1;
    uint32_t src = NONE;
    uint32_t split_i = NONE;
    int64_t split_x = 0;
    uint8_t cls = 0;
    uint32_t chain = NONE;
};

// Build the edited site from fates and extra (block-end) cut points per space.  drop_chain
// (debug, Z29): the first translated link of that chain is left out.
void build_model(const Ctx& cx, const FateMap& fates, const std::vector<std::set<int64_t>>& space_cuts, uint32_t drop_chain,
                 Model& m) {
    const Graph& g = cx.g;
    const uint32_t N = (uint32_t)cx.nodes.size();
    m = Model();
    m.fates = fates;
    std::unordered_map<uint32_t, std::vector<int64_t>> cuts;
    for (const auto& kv : fates) {
        int64_t L = cx.len(kv.first);
        for (const Fate& f : kv.second) {
            if (f.a > 0) cuts[kv.first].push_back(f.a);
            if (f.b < L) cuts[kv.first].push_back(f.b);
        }
    }
    auto tr = [&](const End& e, uint32_t& chain) -> End {
        if (e.i == NONE) return e;
        auto it = fates.find(e.i);
        if (it == fates.end()) return e;
        const Fate* f = e.right ? fate_ending(it->second, e.x) : fate_starting(it->second, e.x);
        if (!f) return e;
        chain = std::min(chain, f->chain);
        return translate(cx, *f, e.right);
    };
    std::vector<RawLink> raw;
    raw.reserve(cx.links.size() + 16);
    for (uint32_t li : cx.links) {
        const Link& l = g.links[li];
        RawLink r;
        uint32_t ch = NONE;
        r.a = tr(side_end(cx, l.a), ch);
        r.b = tr(side_end(cx, l.b), ch);
        r.sr = l.sr;
        r.src = li;
        r.cls = ch == NONE ? 0 : 2;
        r.chain = ch;
        if (r.a.i == NONE || r.b.i == NONE) { m.error = strf("link %u leaves the site", li); return; }
        raw.push_back(r);
    }
    for (const auto& kv : fates) {
        std::vector<int64_t> xs;
        int64_t L = cx.len(kv.first);
        for (const Fate& f : kv.second) {
            if (f.a > 0) xs.push_back(f.a);
            if (f.b < L) xs.push_back(f.b);
        }
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
        for (int64_t x : xs) {
            RawLink r;
            uint32_t ch = NONE;
            End er, el;
            er.i = el.i = kv.first;
            er.x = el.x = x;
            er.right = true;
            el.right = false;
            r.a = tr(er, ch);
            r.b = tr(el, ch);
            r.sr = g.nodes[cx.nodes[kv.first]].sr;
            r.split_i = kv.first;
            r.split_x = x;
            r.cls = 2;
            r.chain = ch;
            raw.push_back(r);
        }
    }
    // own adjacency inside a node: dropped, no cut
    auto own_adj = [&](const RawLink& r) {
        return r.a.i == r.b.i && r.a.x == r.b.x && r.a.right != r.b.right && r.a.x > 0 && r.a.x < cx.len(r.a.i);
    };
    for (const RawLink& r : raw) {
        if (own_adj(r)) continue;
        if (r.a.x > 0 && r.a.x < cx.len(r.a.i)) cuts[r.a.i].push_back(r.a.x);
        if (r.b.x > 0 && r.b.x < cx.len(r.b.i)) cuts[r.b.i].push_back(r.b.x);
    }
    // P of every block end
    for (size_t sp = 0; sp < space_cuts.size() && sp < cx.spaces.size(); ++sp) {
        const TSpace& S = cx.spaces[sp];
        for (int64_t T : space_cuts[sp]) {
            if (T <= S.lo() || T >= S.hi()) continue;
            size_t k = S.elem_starting(T);
            if (T == S.c[k]) continue;
            Handle h = S.walk[k];
            uint32_t i = cx.idx(handle_node(h));
            int64_t d = T - S.c[k];
            cuts[i].push_back(handle_rev(h) ? cx.len(i) - d : d);
        }
    }
    // pieces
    m.first.assign(N + 1, 0);
    m.changed.assign(N, 0);
    for (uint32_t i = 0; i < N; ++i) {
        m.first[i] = (uint32_t)m.pieces.size();
        int64_t L = cx.len(i);
        std::vector<int64_t> xs;
        auto it = cuts.find(i);
        if (it != cuts.end()) {
            xs = it->second;
            std::sort(xs.begin(), xs.end());
            xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
        }
        int64_t p = 0;
        auto fit = fates.find(i);
        for (size_t k = 0; k <= xs.size(); ++k) {
            int64_t q = k < xs.size() ? xs[k] : L;
            if (q <= p || q > L) continue;
            MPiece mp;
            mp.i = i;
            mp.off = p;
            mp.len = q - p;
            if (fit != fates.end()) {
                for (const Fate& f : fit->second) {
                    if (f.a == p && f.b == q) { mp.kept = false; mp.chain = f.chain; break; }
                    if (f.a < q && p < f.b)
                        m.error = strf("a cut at %lld..%lld falls inside a zipped piece of %s", (long long)p, (long long)q,
                                       g.name(cx.nodes[i]).c_str());
                }
            }
            m.pieces.push_back(mp);
            p = q;
        }
        if (m.pieces.size() - m.first[i] > 1 || fit != fates.end()) m.changed[i] = 1;
    }
    m.first[N] = (uint32_t)m.pieces.size();
    if (!m.error.empty()) return;
    auto resolve = [&](const End& e, uint32_t& side) -> bool {
        uint32_t p = e.right ? m.piece_ending(e.i, e.x) : m.piece_starting(e.i, e.x);
        if (p == NONE || !m.pieces[p].kept) return false;
        side = 2 * p + (e.right ? 1u : 0u);
        return true;
    };
    // links in priority order: input (no translated end), new split links, translated
    std::vector<MLink> cand;
    cand.reserve(raw.size() + 16);
    std::vector<MLink> translated;
    for (const RawLink& r : raw) {
        if (own_adj(r)) continue;
        MLink l;
        if (!resolve(r.a, l.sa) || !resolve(r.b, l.sb)) {
            const End& e = r.a;
            m.error = strf("a link end resolves to no live piece (%s %lld %s)", g.name(cx.nodes[e.i]).c_str(), (long long)e.x,
                           e.right ? "R" : "L");
            return;
        }
        l.sr = r.sr;
        l.src = r.src;
        l.split_i = r.split_i;
        l.split_x = r.split_x;
        l.cls = r.cls;
        l.chain = r.chain;
        if (r.cls == 0) cand.push_back(l);
        else translated.push_back(l);
    }
    for (uint32_t i = 0; i < N; ++i) {
        if (fates.count(i)) continue;
        for (uint32_t p = m.first[i]; p + 1 < m.first[i + 1]; ++p) {
            MLink l;
            l.sa = 2 * p + 1;
            l.sb = 2 * (p + 1);
            l.sr = g.nodes[cx.nodes[i]].sr;
            l.split_i = i;
            l.split_x = m.pieces[p + 1].off;
            l.cls = 1;
            cand.push_back(l);
        }
    }
    std::stable_sort(translated.begin(), translated.end(), [&](const MLink& x, const MLink& y) {
        if (x.src != y.src) return x.src < y.src;              // translated input links first (NONE last)
        if (x.split_i != y.split_i) return cx.nodes[x.split_i] < cx.nodes[y.split_i];
        return x.split_x < y.split_x;
    });
    bool dropped_debug = false;
    for (MLink& l : translated) {
        if (!dropped_debug && drop_chain != NONE && l.chain == drop_chain) { dropped_debug = true; continue; }
        cand.push_back(l);
    }
    m.links.reserve(cand.size());
    for (MLink& l : cand) {
        if (!m.keys.insert(Model::key(l.sa, l.sb)).second) continue;
        m.links.push_back(l);
    }
}

// ================================================================ forward-strand SCCs

// Non-trivial SCCs of the forward-strand graph of the model (reference '+', alt pieces as their
// node's creator excursion walks it, inversion-type links ignored), each as the sorted set of the
// site nodes it involves.
std::vector<std::vector<uint32_t>> forward_sccs(const Ctx& cx, const Model& m) {
    const uint32_t P = (uint32_t)m.pieces.size();
    std::vector<uint32_t> off(P + 1, 0);
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    edges.reserve(m.links.size());
    std::vector<uint8_t> selfloop(P, 0);
    for (const MLink& l : m.links) {
        uint32_t pa = l.sa >> 1, pb = l.sb >> 1;
        bool ra = (l.sa & 1u) != 0, rb = (l.sb & 1u) != 0;
        bool fa = cx.orient_rev[m.pieces[pa].i] != 0, fb = cx.orient_rev[m.pieces[pb].i] != 0;
        // the exit side of a '+' piece is R, its entry side L ('-': the reverse)
        bool ab = ra == !fa && rb == fb;
        bool ba = rb == !fb && ra == fa;
        if (ab) edges.push_back(std::make_pair(pa, pb));
        if (ba && !(ab && pa == pb)) edges.push_back(std::make_pair(pb, pa));
        if ((ab || ba) && pa == pb) selfloop[pa] = 1;
    }
    for (const auto& e : edges) off[e.first + 1]++;
    for (uint32_t i = 0; i < P; ++i) off[i + 1] += off[i];
    std::vector<uint32_t> adj(edges.size());
    {
        std::vector<uint32_t> fill(off.begin(), off.end() - 1);
        for (const auto& e : edges) adj[fill[e.first]++] = e.second;
    }
    std::vector<uint32_t> index(P, NONE), low(P, 0);
    std::vector<uint8_t> onst(P, 0);
    std::vector<uint32_t> st;
    struct Fr { uint32_t v, ei; };
    std::vector<Fr> cs;
    uint32_t idx = 0;
    std::vector<std::vector<uint32_t>> out;
    for (uint32_t s0 = 0; s0 < P; ++s0) {
        if (index[s0] != NONE || !m.pieces[s0].kept) continue;
        index[s0] = low[s0] = idx++;
        st.push_back(s0);
        onst[s0] = 1;
        cs.push_back(Fr{s0, off[s0]});
        while (!cs.empty()) {
            uint32_t v = cs.back().v;
            if (cs.back().ei < off[v + 1]) {
                uint32_t w = adj[cs.back().ei++];
                if (index[w] == NONE) {
                    index[w] = low[w] = idx++;
                    st.push_back(w);
                    onst[w] = 1;
                    cs.push_back(Fr{w, off[w]});
                } else if (onst[w]) {
                    low[v] = std::min(low[v], index[w]);
                }
            } else {
                if (low[v] == index[v]) {
                    std::vector<uint32_t> members;
                    while (true) {
                        uint32_t w = st.back();
                        st.pop_back();
                        onst[w] = 0;
                        members.push_back(w);
                        if (w == v) break;
                    }
                    if (members.size() > 1 || selfloop[members[0]]) {
                        std::vector<uint32_t> parents;
                        for (uint32_t p : members) parents.push_back(m.pieces[p].i);
                        std::sort(parents.begin(), parents.end());
                        parents.erase(std::unique(parents.begin(), parents.end()), parents.end());
                        out.push_back(std::move(parents));
                    }
                }
                cs.pop_back();
                if (!cs.empty()) {
                    uint32_t u = cs.back().v;
                    low[u] = std::min(low[u], low[v]);
                }
            }
        }
    }
    std::sort(out.begin(), out.end());
    return out;
}

// an SCC of `after` whose node set is not inside one SCC of `before`
bool new_scc(const std::vector<std::vector<uint32_t>>& after, const std::vector<std::vector<uint32_t>>& before) {
    for (const auto& a : after) {
        bool inside = false;
        for (const auto& b : before)
            if (std::includes(b.begin(), b.end(), a.begin(), a.end())) { inside = true; break; }
        if (!inside) return true;
    }
    return false;
}

// rank-0 back edges: a link joining the R side of reference piece r to the L side of reference
// piece l with start(l) < end(r)
uint64_t ref_back_edges(const Ctx& cx, const Model& m) {
    uint64_t n = 0;
    for (const MLink& l : m.links) {
        uint32_t pa = l.sa >> 1, pb = l.sb >> 1;
        bool ra = (l.sa & 1u) != 0, rb = (l.sb & 1u) != 0;
        const MPiece& A = m.pieces[pa];
        const MPiece& B = m.pieces[pb];
        if (!cx.is_ref[A.i] || !cx.is_ref[B.i] || ra == rb) continue;
        const MPiece& R = ra ? A : B;
        const MPiece& L = ra ? B : A;
        int64_t r_end = cx.g.start(cx.nodes[R.i]) + R.off + R.len;
        int64_t l_start = cx.g.start(cx.nodes[L.i]) + L.off;
        if (l_start < r_end) ++n;
    }
    return n;
}

// ================================================================ chains, pieces

enum PState : uint8_t { P_NEW = 0, P_COMPAT = 1, P_DROP = 2 };

struct XPiece {
    uint32_t part = 0, wpos = 0;
    uint32_t i = NONE;           // site node index of the query node
    bool wrev = false;           // walked '-'
    int64_t a = 0, b = 0;        // node-forward
    int64_t qa = 0, qb = 0;      // walk
    int64_t ta = 0, tb = 0;
    char rel = '+';
    char strand = '+';           // its part's strand
    int64_t aligned = 0;         // '=' and 'X' columns
    uint8_t state = P_NEW;
    std::string why;             // drop rule
    Fate img;                    // P_COMPAT: the accepted image
};

struct Cand {
    bool alt = false;
    uint32_t unit = NONE;        // reference: unit id; alt: AltCandidate::key
    uint32_t space = 0;
    const Excursion* exc = nullptr;
    std::vector<Handle> walk;    // query walk (alt handles, walk orientation)
    std::vector<int64_t> off;    // walk offsets (size walk + 1)
    int64_t wlo = 0, whi = 0;    // window in space coordinates
    std::vector<ChainPart> parts;
    std::vector<XPiece> base;    // pieces after trimming and snapping (acceptance starts from these)
    int64_t score = 0, order_bp = 0;
    int64_t order_w = 0;         // the order key: kept bp weighted by node support
    uint32_t weight = 0;
    ReportRow* row = nullptr;
    bool drop_link = false;      // debug
    // outcome of the last acceptance
    std::string outcome;
    int64_t kept_bp = -1, new_bp = 0;
    std::vector<std::string> dropped;
    std::vector<XPiece> final_pieces;
};

struct State {
    FateMap fates;
    std::vector<std::set<int64_t>> space_cuts;
    std::vector<std::vector<uint32_t>> sccs;     // forward SCCs of the current edit
    std::vector<uint32_t> accepted;              // candidate indices, acceptance order
};

// '=' and 'X' columns whose query base lies in walk interval [qa, qb) of part p
int64_t aligned_cols(const ChainPart& p, int64_t qa, int64_t qb) {
    int64_t c0, c1;
    if (p.strand == '+') { c0 = qa - p.qs; c1 = qb - p.qs; }
    else { c0 = p.qe - qb; c1 = p.qe - qa; }
    int64_t q = 0, n = 0;
    for (const CigarOp& op : p.ops) {
        if (q >= c1) break;
        if (op.op == '=' || op.op == 'X' || op.op == 'M') {
            int64_t a = std::max(q, c0), b = std::min(q + (int64_t)op.len, c1);
            if (b > a) n += b - a;
            q += op.len;
        } else if (op.op == 'I') {
            q += op.len;
        }
    }
    return n;
}

// target point at a piece's walk start / walk end
int64_t wstart_t(const XPiece& p) { return p.strand == '+' ? p.ta : p.tb; }
int64_t wend_t(const XPiece& p) { return p.strand == '+' ? p.tb : p.ta; }
void set_wstart_t(XPiece& p, int64_t t) { if (p.strand == '+') p.ta = t; else p.tb = t; }
void set_wend_t(XPiece& p, int64_t t) { if (p.strand == '+') p.tb = t; else p.ta = t; }

std::string piece_label(const Ctx& cx, uint32_t i, int64_t a, int64_t b) {
    return strf("%s[%lld-%lld)", cx.g.name(cx.nodes[i]).c_str(), (long long)a, (long long)b);
}

// pieces of part j of c over walk interval [x0, x1): spec 6.1 and 6.2's projection
void part_pieces(const Ctx& cx, const Cand& c, uint32_t j, int64_t x0, int64_t x1, std::vector<XPiece>& out) {
    const ChainPart& p = c.parts[j];
    if (x1 <= x0) return;
    size_t k = (size_t)(std::upper_bound(c.off.begin(), c.off.end(), x0) - c.off.begin());
    k = k == 0 ? 0 : k - 1;
    for (; k < c.walk.size() && c.off[k] < x1; ++k) {
        int64_t L = c.off[k + 1] - c.off[k];
        int64_t qa = std::max(c.off[k], x0), qb = std::min(c.off[k] + L, x1);
        if (qb <= qa) continue;
        XPiece x;
        x.part = j;
        x.wpos = (uint32_t)k;
        x.i = cx.idx(handle_node(c.walk[k]));
        x.wrev = handle_rev(c.walk[k]);
        x.qa = qa;
        x.qb = qb;
        if (!x.wrev) { x.a = qa - c.off[k]; x.b = qb - c.off[k]; }
        else { x.a = L - (qb - c.off[k]); x.b = L - (qa - c.off[k]); }
        x.strand = p.strand;
        if (p.strand == '+') { x.ta = project_left_snap(p, qa); x.tb = project_left_snap(p, qb); }
        else { x.ta = project_left_snap(p, qb); x.tb = project_left_snap(p, qa); }
        x.rel = (x.wrev == (p.strand == '-')) ? '+' : '-';
        x.aligned = aligned_cols(p, qa, qb);
        out.push_back(x);
    }
}

// smallest x in [lo, hi] with pred(P(x)), pred false then true along x (P is monotone); hi + 1 if none
template <typename F>
int64_t first_x(const ChainPart& p, int64_t lo, int64_t hi, F pred) {
    int64_t a = lo, b = hi + 1;
    while (a < b) {
        int64_t m = a + (b - a) / 2;
        if (pred(project_left_snap(p, m))) b = m; else a = m + 1;
    }
    return a;
}

// largest x in [lo, hi] with pred(P(x)), pred true then false along x; lo - 1 if none
template <typename F>
int64_t last_x(const ChainPart& p, int64_t lo, int64_t hi, F pred) {
    int64_t a = lo, b = hi + 1;
    while (a < b) {
        int64_t m = a + (b - a) / 2;
        if (pred(project_left_snap(p, m))) a = m + 1; else b = m;
    }
    return a - 1;
}

// Spec 6.1-6.3 for one chain: pieces (node-forward interval, target, rel) of every part, with the
// query/target overlaps that the chain's overlap tolerance G allows between parts trimmed from the
// later part, free block ends snapped, and folds inside the chain dropped.
//
// Overlap trimming.  The query overlap goes first: a part starts after every earlier part's
// (trimmed) query end.  Then its target must not overlap the target of any earlier part.  The
// chain DP lets two parts overlap by up to G at either end of the later one: at its query start
// when it continues a segment (same strand, collinear), but at its query END when it switches
// strand above the previous segment as '-' (frame F: A rc(B), a '+' record running into the
// breakpoint's microhomology) or below it as '+' (frame R), or when a part of a '-' segment
// touches the segment's floor.  So the overlap is cut from whichever end of the later part's
// target it lies at, by moving the part's query start or end (P is monotone).  A part left with
// nothing is a fold: its pieces are reported as dropped "(fold)" rather than lost silently.
void prepare_chain(const Ctx& cx, Cand& c, const std::vector<int64_t>& boundaries) {
    const int64_t G = cx.opt.G;
    const uint32_t np = (uint32_t)c.parts.size();
    std::vector<int64_t> x0(np, 0), x1(np, 0), hlo(np, 0), hhi(np, 0);
    std::vector<uint8_t> live(np, 0);
    std::vector<XPiece> folded;          // pieces of parts that vanished (reported, never zipped)
    int64_t qend = INT64_MIN;            // largest query end of the earlier live parts
    for (uint32_t j = 0; j < np; ++j) {
        const ChainPart& p = c.parts[j];
        const bool fwd = p.strand == '+';
        int64_t a = std::max<int64_t>(p.qs, 0), b = p.qe;
        a = std::max(a, qend);                                    // query overlap
        if (a >= b) continue;                                     // all of its query is in earlier parts
        const int64_t qa = a;
        bool gone = false;
        {
            int64_t lo = project_left_snap(p, fwd ? a : b), hi = project_left_snap(p, fwd ? b : a);
            for (uint32_t k = 0; k < j && !gone; ++k) {
                if (!live[k] || hhi[k] <= lo || hlo[k] >= hi) continue;
                bool bottom = hlo[k] <= lo, top = hhi[k] >= hi;
                if (bottom && top) { gone = true; break; }    // the earlier part covers this one's target
                if (!bottom && !top) bottom = hhi[k] - lo <= hi - hlo[k];   // strictly inside: the cheaper end
                if (bottom) lo = hhi[k];
                else hi = hlo[k];
                if (lo >= hi) gone = true;
            }
            if (!gone) {
                // [lo, hi) back to the query: '+' targets grow with x, '-' targets shrink
                if (fwd) {
                    a = first_x(p, a, b, [&](int64_t t) { return t >= lo; });
                    b = last_x(p, a, b, [&](int64_t t) { return t <= hi; });
                } else {
                    a = first_x(p, a, b, [&](int64_t t) { return t <= hi; });
                    b = last_x(p, a, b, [&](int64_t t) { return t >= lo; });
                }
                gone = a >= b;
            }
        }
        if (gone) {
            part_pieces(cx, c, j, qa, p.qe, folded);   // its query beyond the earlier parts stays alt
            continue;
        }
        x0[j] = a;
        x1[j] = b;
        int64_t ta = project_left_snap(p, a), tb = project_left_snap(p, b);
        hlo[j] = std::min(ta, tb);
        hhi[j] = std::max(ta, tb);
        live[j] = 1;
        qend = std::max(qend, b);
    }
    for (XPiece& x : folded) {
        x.state = P_DROP;
        x.why = "(fold)";
    }
    std::vector<XPiece> pcs;
    std::vector<std::pair<size_t, size_t>> range(np, std::make_pair((size_t)0, (size_t)0));
    for (uint32_t j = 0; j < np; ++j) {
        size_t b0 = pcs.size();
        if (live[j]) part_pieces(cx, c, j, x0[j], x1[j], pcs);
        range[j] = std::make_pair(b0, pcs.size());
    }
    // a partial node piece shorter than G stays alt: differences shorter than G are absorbed
    for (XPiece& p : pcs)
        if (p.b - p.a < G && p.b - p.a < cx.len(p.i)) { p.state = P_DROP; p.why = "(short)"; }
    auto first_live = [&](uint32_t j) -> size_t {
        for (size_t k = range[j].first; k < range[j].second; ++k)
            if (pcs[k].state != P_DROP) return k;
        return SIZE_MAX;
    };
    auto last_live = [&](uint32_t j) -> size_t {
        for (size_t k = range[j].second; k-- > range[j].first;)
            if (pcs[k].state != P_DROP) return k;
        return SIZE_MAX;
    };
    auto hull_of = [&](uint32_t j, int64_t& lo, int64_t& hi) {
        lo = INT64_MAX;
        hi = INT64_MIN;
        for (size_t k = range[j].first; k < range[j].second; ++k) {
            if (pcs[k].state == P_DROP) continue;
            lo = std::min(lo, pcs[k].ta);
            hi = std::max(hi, pcs[k].tb);
        }
    };
    // snap free block ends (spec 6.2): within G of a window end or an existing node boundary of
    // the target space, and the move shorter than the end piece's target; never at a junction
    // shared with another part, never into another part's target
    for (uint32_t j = 0; j < np; ++j) {
        if (first_live(j) == SIZE_MAX) continue;
        for (int end = 0; end < 2; ++end) {
            size_t ek = end == 0 ? first_live(j) : last_live(j);
            int64_t T = end == 0 ? wstart_t(pcs[ek]) : wend_t(pcs[ek]);
            bool shared = false;
            for (uint32_t k = 0; k < np; ++k) {
                if (k == j || first_live(k) == SIZE_MAX) continue;
                if (wstart_t(pcs[first_live(k)]) == T || wend_t(pcs[last_live(k)]) == T) shared = true;
            }
            if (shared) continue;
            int64_t tlen = pcs[ek].tb - pcs[ek].ta;
            int64_t best = T, bestd = INT64_MAX;
            auto consider = [&](int64_t B) {
                if (B < c.wlo || B > c.whi) return;
                int64_t d = iabs64(B - T);
                if (d == 0 || d > G || d >= tlen) return;
                if (d < bestd || (d == bestd && B < best)) { bestd = d; best = B; }
            };
            consider(c.wlo);
            consider(c.whi);
            for (auto it = std::lower_bound(boundaries.begin(), boundaries.end(), T - G); it != boundaries.end() && *it <= T + G; ++it)
                consider(*it);
            if (bestd == INT64_MAX) continue;
            XPiece trial = pcs[ek];
            if (end == 0) set_wstart_t(trial, best); else set_wend_t(trial, best);
            if (trial.tb <= trial.ta) continue;
            int64_t lo = INT64_MAX, hi = INT64_MIN;
            for (size_t k = range[j].first; k < range[j].second; ++k) {
                if (pcs[k].state == P_DROP) continue;
                const XPiece& z = k == ek ? trial : pcs[k];
                lo = std::min(lo, z.ta);
                hi = std::max(hi, z.tb);
            }
            bool clash = false;
            for (uint32_t k = 0; k < np; ++k) {
                if (k == j || first_live(k) == SIZE_MAX) continue;
                int64_t l2, h2;
                hull_of(k, l2, h2);
                if (std::max(lo, l2) < std::min(hi, h2)) clash = true;
            }
            if (!clash) pcs[ek] = trial;
        }
    }
    // folds inside the chain: a piece whose target overlaps an earlier piece of another part
    for (size_t k = 0; k < pcs.size(); ++k) {
        if (pcs[k].state == P_DROP) continue;
        for (size_t k2 = 0; k2 < k; ++k2) {
            if (pcs[k2].part == pcs[k].part || pcs[k2].state == P_DROP) continue;
            if (std::max(pcs[k].ta, pcs[k2].ta) < std::min(pcs[k].tb, pcs[k2].tb)) {
                pcs[k].state = P_DROP;
                pcs[k].why = "(fold)";
                break;
            }
        }
    }
    // the parts that vanished, after every live piece (so that piece order still gives junctions)
    pcs.insert(pcs.end(), folded.begin(), folded.end());
    c.base = std::move(pcs);
    c.score = 0;
    for (const ChainPart& p : c.parts) c.score += p.n_eq;
}

// junction neighbours: consecutive pieces of one part, or of two parts whose junction lands on
// one target point
bool junction(const XPiece& x, const XPiece& y) {
    if (x.qb != y.qa) return false;
    if (x.part == y.part) return true;
    return wend_t(x) == wstart_t(y);
}

// ================================================================ images (V2, V3, --check)

struct ImgEl {
    bool run = false;
    uint32_t piece = NONE;       // !run
    bool rev = false;
    uint32_t space = 0;          // run: [t0, t1) of space, walked reversed if rev
    int64_t t0 = 0, t1 = 0;
};

void append_oriented(std::string& out, const std::string& s, int64_t off, int64_t len, bool rev) {
    if (!rev) { out.append(s, (size_t)off, (size_t)len); return; }
    for (int64_t k = off + len; k-- > off;) out.push_back(comp_base(s[(size_t)k]));
}

// pieces of the space covering [t0, t1) in space order ('+' run); false if a boundary is missing
bool space_pieces(const Ctx& cx, const Model& m, const TSpace& S, int64_t t0, int64_t t1, std::vector<std::pair<uint32_t, bool>>& out,
                  std::string& why) {
    if (t1 <= t0) return true;
    size_t k = S.elem_starting(t0);
    for (; k + 1 < S.c.size() && S.c[k] < t1; ++k) {
        int64_t s0 = std::max(t0, S.c[k]), s1 = std::min(t1, S.c[k + 1]);
        if (s1 <= s0) continue;
        Handle h = S.walk[k];
        uint32_t i = cx.idx(handle_node(h));
        int64_t L = cx.len(i);
        int64_t na, nb;
        if (!handle_rev(h)) { na = s0 - S.c[k]; nb = s1 - S.c[k]; }
        else { na = L - (s1 - S.c[k]); nb = L - (s0 - S.c[k]); }
        uint32_t pa = m.piece_starting(i, na), pb = m.piece_ending(i, nb);
        if (pa == NONE || pb == NONE || pb < pa) {
            why = strf("no piece boundary of %s at %lld/%lld (target %lld-%lld)", cx.g.name(cx.nodes[i]).c_str(), (long long)na, (long long)nb,
                       (long long)t0, (long long)t1);
            return false;
        }
        if (!handle_rev(h))
            for (uint32_t p = pa; p <= pb; ++p) out.push_back(std::make_pair(p, false));
        else
            for (uint32_t p = pb + 1; p-- > pa;) out.push_back(std::make_pair(p, true));
    }
    return true;
}

// the input sequence of [t0, t1) of a space
void space_seq(const Ctx& cx, const TSpace& S, int64_t t0, int64_t t1, std::string& out) {
    if (t1 <= t0) return;
    size_t k = S.elem_starting(t0);
    for (; k + 1 < S.c.size() && S.c[k] < t1; ++k) {
        int64_t s0 = std::max(t0, S.c[k]), s1 = std::min(t1, S.c[k + 1]);
        if (s1 <= s0) continue;
        Handle h = S.walk[k];
        const std::string& seq = cx.g.nodes[handle_node(h)].seq;
        int64_t L = (int64_t)seq.size();
        if (!handle_rev(h)) append_oriented(out, seq, s0 - S.c[k], s1 - s0, false);
        else append_oriented(out, seq, L - (s1 - S.c[k]), s1 - s0, true);
    }
}

// The image of an input walk (anchors included) in the model.  Returns "" when it exists link by
// link, spells the predicted sequence (kept pieces as they are, zipped pieces as their targets
// read from the input) and reaches no target piece twice; else "V2" or "V3" with a message.
std::string check_image(const Ctx& cx, const Model& m, const std::vector<Handle>& walk, std::string& msg) {
    std::vector<ImgEl> els;
    for (Handle h : walk) {
        uint32_t i = cx.idx(handle_node(h));
        if (i == NONE) { msg = "walk leaves the site at " + cx.g.handle_str(h); return "V2"; }
        bool orev = handle_rev(h);
        uint32_t n = m.npieces(i);
        auto fit = m.fates.find(i);
        for (uint32_t k = 0; k < n; ++k) {
            uint32_t p = orev ? m.first[i] + n - 1 - k : m.first[i] + k;
            const MPiece& mp = m.pieces[p];
            ImgEl e;
            if (mp.kept) {
                e.piece = p;
                e.rev = orev;
            } else {
                const Fate* f = fit == m.fates.end() ? nullptr : fate_starting(fit->second, mp.off);
                if (!f || f->b != mp.off + mp.len) { msg = "zipped piece without a fate"; return "V2"; }
                e.run = true;
                e.space = f->space;
                e.t0 = f->ta;
                e.t1 = f->tb;
                e.rev = (f->rel == '-') != orev;
                if (!els.empty() && els.back().run && els.back().space == e.space && els.back().rev == e.rev) {
                    ImgEl& b = els.back();
                    if (!e.rev && b.t1 == e.t0) { b.t1 = e.t1; continue; }
                    if (e.rev && b.t0 == e.t1) { b.t0 = e.t0; continue; }
                }
            }
            els.push_back(e);
        }
    }
    std::vector<std::pair<uint32_t, bool>> img;
    std::vector<uint8_t> from_run;
    std::string pred, spell;
    for (const ImgEl& e : els) {
        if (!e.run) {
            img.push_back(std::make_pair(e.piece, e.rev));
            from_run.push_back(0);
            const MPiece& mp = m.pieces[e.piece];
            append_oriented(pred, cx.g.nodes[cx.nodes[mp.i]].seq, mp.off, mp.len, e.rev);
            continue;
        }
        const TSpace& S = cx.spaces[e.space];
        std::vector<std::pair<uint32_t, bool>> sp;
        std::string why;
        if (!space_pieces(cx, m, S, e.t0, e.t1, sp, why)) { msg = why; return "V2"; }
        if (e.rev) {
            std::reverse(sp.begin(), sp.end());
            for (auto& x : sp) x.second = !x.second;
        }
        for (auto& x : sp) { img.push_back(x); from_run.push_back(1); }
        std::string t;
        space_seq(cx, S, e.t0, e.t1, t);
        if (e.rev) t = revcomp(t);
        pred += t;
    }
    for (size_t k = 0; k < img.size(); ++k) {
        const MPiece& mp = m.pieces[img[k].first];
        if (!mp.kept) { msg = "image visits a zipped piece"; return "V2"; }
        append_oriented(spell, cx.g.nodes[cx.nodes[mp.i]].seq, mp.off, mp.len, img[k].second);
        if (k == 0) continue;
        uint32_t p1 = img[k - 1].first, p2 = img[k].first;
        bool r1 = img[k - 1].second, r2 = img[k].second;
        uint32_t s1 = 2 * p1 + (r1 ? 0u : 1u);      // exit side
        uint32_t s2 = 2 * p2 + (r2 ? 1u : 0u);      // entry side
        if (!m.has_link(s1, s2)) {
            const MPiece& a = m.pieces[p1];
            const MPiece& b = m.pieces[p2];
            msg = strf("no link %s%s -> %s%s", piece_label(cx, a.i, a.off, a.off + a.len).c_str(), r1 ? "-" : "+",
                       piece_label(cx, b.i, b.off, b.off + b.len).c_str(), r2 ? "-" : "+");
            return "V2";
        }
    }
    if (pred != spell) {
        msg = strf("image spells %zu bp, the prediction %zu bp, and they differ", spell.size(), pred.size());
        return "V2";
    }
    // V3: a target piece the walk reaches through zipped runs appears once, and not otherwise
    std::unordered_map<uint32_t, int> runs, other;
    for (size_t k = 0; k < img.size(); ++k) (from_run[k] ? runs : other)[img[k].first]++;
    for (const auto& kv : runs) {
        if (kv.second > 1 || other.count(kv.first)) {
            const MPiece& a = m.pieces[kv.first];
            msg = strf("piece %s occurs twice in the image (copy collapse or fold)", piece_label(cx, a.i, a.off, a.off + a.len).c_str());
            return "V3";
        }
    }
    return "";
}

} // namespace

// ================================================================ the planner

struct SitePlanner::Impl {
    const Graph& g;
    const SiteData& sd;
    const Options& opt;
    Ctx cx;
    std::vector<int64_t> boundaries0;          // node boundaries of the reference space
    std::vector<Cand> cands;
    std::vector<uint32_t> seq;                 // candidates in acceptance order
    State st;
    bool have_base = false;
    std::vector<std::vector<uint32_t>> base_sccs;
    uint64_t base_back = 0;
    std::unordered_map<uint32_t, std::vector<uint64_t>> reach;   // sg local node -> reachable local nodes
    std::set<NodeId> used;
    std::map<NodeId, uint32_t> alt_target;     // v2: representative node -> the alt space it is the target of
    std::map<std::vector<Handle>, uint32_t> alt_space;   // v2: target walk (u, branch, v) -> its space
    uint32_t drop_chain_unit = NONE;           // debug drop-link (unit id)
    bool finished = false;
    // the planner's own time: inside accept_reference, accept_alt and finish only (the alt pass's
    // alignments between them, and their waits for an aligner slot, are not planning)
    double own_s = 0;
    uint64_t n_models = 0, n_validations = 0;

    Impl(const Graph& g_, const SiteData& sd_, const Options& opt_) : g(g_), sd(sd_), opt(opt_), cx(g_, sd_, opt_) {
        build_ctx(cx);
        const TSpace& S = cx.spaces[0];
        boundaries0.assign(S.c.begin() + 1, S.c.end() - 1);
        st.space_cuts.assign(1, std::set<int64_t>());
    }

    // The unedited site as a model (forward SCCs and back edges of the input).  It cannot fail for a
    // site that does not leak (every link of an interior node, and on A.R and B.L, stays in the
    // site); if it ever does, that one site is not edited -- every candidate is reverted:model and
    // counted as a reverted chain (--strict exits 3) -- instead of failing the whole run.
    bool broken = false;
    uint32_t broken_cands = 0;
    bool baseline() {
        if (have_base) return !broken;
        have_base = true;
        Model m;
        build_model(cx, FateMap(), std::vector<std::set<int64_t>>(), NONE, m);
        if (!m.error.empty()) {
            broken = true;
            ZLOG("warning: site %s: the unedited site model fails (%s): the site is not edited", cx.s.label(g).c_str(), m.error.c_str());
            return false;
        }
        base_sccs = forward_sccs(cx, m);
        base_back = ref_back_edges(cx, m);
        st.sccs = base_sccs;
        return true;
    }

    // nodes i and j (site node indices) reach each other over alt nodes, in either orientation
    bool can_reach(uint32_t i, uint32_t j) {
        if (i == j) return true;
        uint32_t li = sd.sg.local_node(cx.nodes[i]), lj = sd.sg.local_node(cx.nodes[j]);
        if (li == NONE || lj == NONE) return false;
        auto it = reach.find(li);
        if (it == reach.end()) {
            std::vector<uint64_t> bits((sd.sg.n_nodes() + 63) / 64, 0);
            std::vector<uint8_t> seen(sd.sg.n_handles(), 0);
            std::vector<uint32_t> stk{2 * li, 2 * li + 1};
            seen[2 * li] = seen[2 * li + 1] = 1;
            while (!stk.empty()) {
                uint32_t h = stk.back();
                stk.pop_back();
                bits[(h >> 1) / 64] |= 1ULL << ((h >> 1) % 64);
                for (uint32_t k = sd.sg.succ_off[h]; k < sd.sg.succ_off[h + 1]; ++k) {
                    uint32_t w = sd.sg.succ[k];
                    if (!seen[w]) { seen[w] = 1; stk.push_back(w); }
                }
            }
            it = reach.emplace(li, std::move(bits)).first;
        }
        return ((it->second[lj / 64] >> (lj % 64)) & 1ULL) != 0;
    }

    // (B): the target lies inside Allowed(node) (reference pass) and inside the window
    bool allowed_B(const Cand& c, const XPiece& p) const {
        if (p.tb <= p.ta || p.ta < c.wlo || p.tb > c.whi) return false;
        if (c.alt) return true;
        return sd.allowed_of(cx.nodes[p.i]).contains(p.ta, p.tb);
    }

    // (G), the whole-record rule (observed walks, v3; reference-pass pieces): a record of the walk
    // source that walks the piece's node also reads its target elsewhere -- the zip would collapse
    // that haplotype's second copy.  The GAF audit's reading rule, applied at acceptance.
    bool reads_G(const Cand& c, const XPiece& p) const {
        if (c.alt || !sd.walks) return false;
        return sd.walks->reads_target(*sd.site, cx.nodes[p.i], p.b - p.a, p.ta, p.tb);
    }

    // Spec 6.7 for one candidate against state `s`.  Returns "" when accepted (then `add` holds its
    // new fates, `cuts` its block-end points and `sccs_after` the forward SCCs with it), else the
    // chain-level refusal.
    std::string try_accept(uint32_t ci, const State& s, FateMap& add, std::set<int64_t>& cuts, std::vector<std::vector<uint32_t>>& sccs_after) {
        Cand& c = cands[ci];
        const int64_t G = opt.G;
        std::vector<XPiece> pcs = c.base;
        c.dropped.clear();
        auto drop = [&](XPiece& p, const char* why) {
            if (p.state == P_DROP) return;
            p.state = P_DROP;
            p.why = why;
        };
        // no piece is deleted without an aligned target base; (B) after snapping
        for (XPiece& p : pcs) {
            if (p.state == P_DROP) continue;
            if (p.aligned <= 0 || p.tb <= p.ta) drop(p, "(empty)");
            else if (!allowed_B(c, p)) drop(p, "(B)");
        }
        // (A) compatibility with the accepted fates of the node
        for (XPiece& p : pcs) {
            if (p.state != P_NEW) continue;
            auto it = s.fates.find(p.i);
            if (it == s.fates.end()) continue;
            std::vector<const Fate*> ov;
            for (const Fate& f : it->second)
                if (f.a < p.b && p.a < f.b) ov.push_back(&f);
            if (ov.empty()) continue;
            Fate run = *ov[0];
            bool ok = true;
            for (size_t k = 1; k < ov.size() && ok; ++k) {
                const Fate& f = *ov[k];
                bool contig = f.a == run.b && f.space == run.space && f.rel == run.rel && (run.rel == '+' ? f.ta == run.tb : f.tb == run.ta);
                if (!contig) { ok = false; break; }
                run.b = f.b;
                if (run.rel == '+') run.tb = f.tb; else run.ta = f.ta;
            }
            if (ok && run.space == c.space && run.rel == p.rel && iabs64(run.a - p.a) <= G && iabs64(run.b - p.b) <= G &&
                iabs64(run.ta - p.ta) <= G && iabs64(run.tb - p.tb) <= G) {
                p.state = P_COMPAT;
                p.img = run;
            } else {
                drop(p, "(A)");
            }
        }
        // the later chain's junctions next to a compatible piece are pinned to the accepted image;
        // a shift <= G is absorbed by the neighbouring piece, a larger one drops it
        auto img_wstart = [](const XPiece& p) -> int64_t { return ((p.rel == '+') != p.wrev) ? p.img.ta : p.img.tb; };
        auto img_wend = [](const XPiece& p) -> int64_t { return ((p.rel == '+') != p.wrev) ? p.img.tb : p.img.ta; };
        auto requery = [&](XPiece& q) {
            int64_t L = c.off[q.wpos + 1] - c.off[q.wpos];
            if (!q.wrev) { q.qa = c.off[q.wpos] + q.a; q.qb = c.off[q.wpos] + q.b; }
            else { q.qa = c.off[q.wpos] + (L - q.b); q.qb = c.off[q.wpos] + (L - q.a); }
        };
        for (size_t k = 0; k < pcs.size(); ++k) {
            if (pcs[k].state != P_COMPAT) continue;
            const XPiece& p = pcs[k];
            if (k > 0 && pcs[k - 1].state == P_NEW && junction(pcs[k - 1], p)) {
                XPiece& q = pcs[k - 1];
                int64_t T0 = img_wstart(p);
                if (iabs64(T0 - wend_t(q)) > G) drop(q, "(A)");
                else {
                    set_wend_t(q, T0);
                    if (q.i == p.i && q.wpos == p.wpos) {
                        if (!q.wrev) q.b = p.img.a; else q.a = p.img.b;
                        requery(q);
                    }
                    if (q.tb <= q.ta || q.b <= q.a || !allowed_B(c, q)) drop(q, "(A)");
                }
            }
            if (k + 1 < pcs.size() && pcs[k + 1].state == P_NEW && junction(p, pcs[k + 1])) {
                XPiece& q = pcs[k + 1];
                int64_t T0 = img_wend(p);
                if (iabs64(T0 - wstart_t(q)) > G) drop(q, "(A)");
                else {
                    set_wstart_t(q, T0);
                    if (q.i == p.i && q.wpos == p.wpos) {
                        if (!q.wrev) q.a = p.img.b; else q.b = p.img.a;
                        requery(q);
                    }
                    if (q.tb <= q.ta || q.b <= q.a || !allowed_B(c, q)) drop(q, "(A)");
                }
            }
        }
        // (A) one fate per node interval: new pieces against accepted fates and earlier pieces
        for (size_t k = 0; k < pcs.size(); ++k) {
            XPiece& p = pcs[k];
            if (p.state != P_NEW) continue;
            auto it = s.fates.find(p.i);
            bool clash = false;
            if (it != s.fates.end())
                for (const Fate& f : it->second)
                    if (f.a < p.b && p.a < f.b) clash = true;
            for (size_t k2 = 0; k2 < k && !clash; ++k2) {
                const XPiece& o = pcs[k2];
                if (o.state == P_DROP || o.i != p.i) continue;
                int64_t oa = o.state == P_COMPAT ? o.img.a : o.a, ob = o.state == P_COMPAT ? o.img.b : o.b;
                if (oa < p.b && p.a < ob) clash = true;
            }
            if (clash) drop(p, "(A)");
        }
        // (G) on the new pieces' final targets (the pieces the edit would delete and the audit reads)
        for (XPiece& p : pcs)
            if (p.state == P_NEW && reads_G(c, p)) drop(p, "(G)");
        // (C) disjoint targets for pieces of different accepted chains whose nodes reach each other
        for (XPiece& p : pcs) {
            if (p.state != P_NEW) continue;
            bool bad = false;
            for (const auto& kv : s.fates) {
                for (const Fate& f : kv.second) {
                    if (f.space != c.space || f.chain == ci) continue;
                    if (std::max(f.ta, p.ta) >= std::min(f.tb, p.tb)) continue;
                    if (can_reach(p.i, kv.first)) { bad = true; break; }
                }
                if (bad) break;
            }
            if (bad) drop(p, "(C)");
        }
        // chain rules: >= b bp remain (new and compatible pieces); (D)
        int64_t remain = 0, newbp = 0;
        for (const XPiece& p : pcs) {
            if (p.state == P_NEW) { remain += p.qb - p.qa; newbp += p.qb - p.qa; }
            else if (p.state == P_COMPAT) remain += p.qb - p.qa;
        }
        for (const XPiece& p : pcs)
            if (p.state == P_DROP) c.dropped.push_back(piece_label(cx, p.i, p.a, p.b) + ":" + p.why);
        c.kept_bp = remain;
        c.new_bp = newbp;
        c.final_pieces = pcs;
        if (remain < opt.b) return "trimmed-below-b";
        add.clear();
        cuts.clear();
        std::map<uint32_t, std::pair<size_t, size_t>> live_end;    // part -> first and last live piece
        for (size_t k = 0; k < pcs.size(); ++k) {
            if (pcs[k].state == P_DROP) continue;
            auto it = live_end.find(pcs[k].part);
            if (it == live_end.end()) live_end[pcs[k].part] = std::make_pair(k, k);
            else it->second.second = k;
        }
        for (size_t k = 0; k < pcs.size(); ++k) {
            const XPiece& p = pcs[k];
            if (p.state != P_NEW) continue;
            Fate f;
            f.a = p.a;
            f.b = p.b;
            f.space = c.space;
            f.ta = p.ta;
            f.tb = p.tb;
            f.rel = p.rel;
            f.chain = ci;
            add[p.i].push_back(f);
            const std::pair<size_t, size_t>& le = live_end[p.part];
            if (k == le.first) cuts.insert(wstart_t(p));          // P of the block's ends
            if (k == le.second) cuts.insert(wend_t(p));
        }
        if (add.empty()) { sccs_after = s.sccs; return ""; }
        FateMap trial = s.fates;
        for (const auto& kv : add) {
            std::vector<Fate>& v = trial[kv.first];
            v.insert(v.end(), kv.second.begin(), kv.second.end());
            std::sort(v.begin(), v.end(), [](const Fate& x, const Fate& y) { return x.a < y.a; });
        }
        std::vector<std::set<int64_t>> sc = s.space_cuts;
        if (sc.size() < cx.spaces.size()) sc.resize(cx.spaces.size());
        sc[c.space].insert(cuts.begin(), cuts.end());
        Model m;
        build_model(cx, trial, sc, NONE, m);
        ++n_models;
        if (!m.error.empty()) {
            ZLOG("warning: site %s: internal edit error while checking (D): %s", cx.s.label(g).c_str(), m.error.c_str());
            return "reverted:model";
        }
        sccs_after = forward_sccs(cx, m);
        if (new_scc(sccs_after, s.sccs)) return "(D)";
        return "";
    }

    void commit(uint32_t ci, State& s, const FateMap& add, const std::set<int64_t>& cuts, std::vector<std::vector<uint32_t>>& sccs_after) {
        for (const auto& kv : add) {
            std::vector<Fate>& v = s.fates[kv.first];
            v.insert(v.end(), kv.second.begin(), kv.second.end());
            std::sort(v.begin(), v.end(), [](const Fate& x, const Fate& y) { return x.a < y.a; });
        }
        if (s.space_cuts.size() < cx.spaces.size()) s.space_cuts.resize(cx.spaces.size());
        s.space_cuts[cands[ci].space].insert(cuts.begin(), cuts.end());
        s.sccs.swap(sccs_after);
        s.accepted.push_back(ci);
    }

    // greedy acceptance of candidate ci against st
    void accept_one(uint32_t ci) {
        FateMap add;
        std::set<int64_t> cuts;
        std::vector<std::vector<uint32_t>> after;
        std::string r = try_accept(ci, st, add, cuts, after);
        Cand& c = cands[ci];
        if (r.empty()) {
            commit(ci, st, add, cuts, after);
            c.outcome = "zipped";
            for (const auto& kv : add) used.insert(cx.nodes[kv.first]);
            // the reference nodes a reference chain lands on are used as targets: alt-vs-alt never
            // takes them as a bound (spec section 5)
            if (!c.alt) {
                const TSpace& S = cx.spaces[0];
                for (const auto& kv : add)
                    for (const Fate& f : kv.second)
                        for (size_t k = S.elem_starting(f.ta); k + 1 < S.c.size() && S.c[k] < f.tb; ++k)
                            if (S.c[k + 1] > f.ta) used.insert(handle_node(S.walk[k]));
            }
        } else {
            c.outcome = r;
        }
    }

    uint32_t drop_index(const State& s) const {
        for (uint32_t ci : s.accepted)
            if (cands[ci].drop_link) return ci;
        return NONE;
    }

    // ---------------------------------------------------------------- validation (V1-V7)

    std::string validate(const Model& m, std::string& msg) {
        ++n_validations;
        if (!m.error.empty()) { msg = m.error; return "V5"; }
        // V1: the reference walk A+ ... B+ exists link by link and tiles start(A) .. end(B)
        {
            std::vector<uint32_t> walk;
            for (uint32_t p = m.first[cx.iA]; p < m.first[cx.iA + 1]; ++p) walk.push_back(p);
            for (NodeId b : cx.s.backbone) {
                uint32_t i = cx.idx(b);
                for (uint32_t p = m.first[i]; p < m.first[i + 1]; ++p) walk.push_back(p);
            }
            for (uint32_t p = m.first[cx.iB]; p < m.first[cx.iB + 1]; ++p) walk.push_back(p);
            int64_t pos = g.start(cx.s.A);
            for (size_t k = 0; k < walk.size(); ++k) {
                const MPiece& mp = m.pieces[walk[k]];
                if (!mp.kept) { msg = "a reference piece is deleted"; return "V1"; }
                int64_t so = g.start(cx.nodes[mp.i]) + mp.off;
                if (so != pos) { msg = strf("reference pieces do not tile at %lld", (long long)pos); return "V1"; }
                pos += mp.len;
                if (k > 0 && !m.has_link(2 * walk[k - 1] + 1, 2 * walk[k])) {
                    msg = strf("the reference walk has no link into %s", piece_label(cx, mp.i, mp.off, mp.off + mp.len).c_str());
                    return "V1";
                }
            }
            if (pos != g.end(cx.s.B)) { msg = "the reference walk does not reach the end of B"; return "V1"; }
        }
        // V2 and V3: the image of every unit whose walk the edit touches
        for (const Unit& u : sd.units) {
            std::vector<Handle> w;
            w.reserve(u.exc.alts.size() + 2);
            w.push_back(u.exc.dep);
            w.insert(w.end(), u.exc.alts.begin(), u.exc.alts.end());
            w.push_back(u.exc.arr);
            bool touched = false;
            for (Handle h : w) {
                uint32_t i = cx.idx(handle_node(h));
                if (i == NONE || m.changed[i]) { touched = true; break; }
            }
            if (!touched) continue;
            std::string r = check_image(cx, m, w, msg);
            if (!r.empty()) {
                msg = strf("unit %u (%s): %s", u.id, owners_str(g, u.exc.owners).c_str(), msg.c_str());
                return r;
            }
        }
        // --check: every anchored path at sites with <= 5000 of them
        if (opt.check) {
            std::string r = check_paths(m, msg);
            if (!r.empty()) return r;
        }
        // V4: no new forward-strand SCC; count <= the input's
        {
            std::vector<std::vector<uint32_t>> sccs = forward_sccs(cx, m);
            if (new_scc(sccs, base_sccs) || sccs.size() > base_sccs.size()) {
                msg = strf("%zu forward-strand SCC(s), the input site has %zu", sccs.size(), base_sccs.size());
                return "V4";
            }
        }
        // V5: links local, targets inside their windows, rank-0 back edges not increased
        {
            for (const MLink& l : m.links) {
                for (uint32_t s2 : {l.sa, l.sb}) {
                    const MPiece& mp = m.pieces[s2 >> 1];
                    bool right = (s2 & 1u) != 0;
                    if ((mp.i == cx.iA && !right) || (mp.i == cx.iB && right)) { msg = "a link leaves the site through A.L or B.R"; return "V5"; }
                }
            }
            for (const auto& kv : m.fates)
                for (const Fate& f : kv.second) {
                    const Cand& c = cands[f.chain];
                    if (f.ta < c.wlo || f.tb > c.whi || f.tb <= f.ta) {
                        msg = strf("target %lld-%lld of %s is outside its window %lld-%lld", (long long)f.ta, (long long)f.tb,
                                   piece_label(cx, kv.first, f.a, f.b).c_str(), (long long)c.wlo, (long long)c.whi);
                        return "V5";
                    }
                }
            uint64_t be = ref_back_edges(cx, m);
            if (be > base_back) {
                msg = strf("%llu rank-0 back edge(s), the input site has %llu", (unsigned long long)be, (unsigned long long)base_back);
                return "V5";
            }
        }
        // V6: every kept piece connects to the site's reference pieces
        {
            const uint32_t P = (uint32_t)m.pieces.size();
            std::vector<uint32_t> off(P + 1, 0);
            for (const MLink& l : m.links) { off[(l.sa >> 1) + 1]++; off[(l.sb >> 1) + 1]++; }
            for (uint32_t i = 0; i < P; ++i) off[i + 1] += off[i];
            std::vector<uint32_t> adj(off[P]);
            {
                std::vector<uint32_t> fill(off.begin(), off.end() - 1);
                for (const MLink& l : m.links) { adj[fill[l.sa >> 1]++] = l.sb >> 1; adj[fill[l.sb >> 1]++] = l.sa >> 1; }
            }
            std::vector<uint8_t> seen(P, 0);
            std::vector<uint32_t> stk;
            for (uint32_t p = 0; p < P; ++p)
                if (m.pieces[p].kept && cx.is_ref[m.pieces[p].i]) { seen[p] = 1; stk.push_back(p); }
            while (!stk.empty()) {
                uint32_t x = stk.back();
                stk.pop_back();
                for (uint32_t k = off[x]; k < off[x + 1]; ++k)
                    if (!seen[adj[k]]) { seen[adj[k]] = 1; stk.push_back(adj[k]); }
            }
            for (uint32_t p = 0; p < P; ++p)
                if (m.pieces[p].kept && !seen[p]) {
                    msg = "orphaned piece " + piece_label(cx, m.pieces[p].i, m.pieces[p].off, m.pieces[p].off + m.pieces[p].len);
                    return "V6";
                }
        }
        // V7: pieces tile their parents (they inherit SN and SR, SO = parent SO + offset);
        // removed <= target / i per chain
        {
            for (uint32_t i = 0; i < (uint32_t)cx.nodes.size(); ++i) {
                int64_t pos = 0;
                for (uint32_t p = m.first[i]; p < m.first[i + 1]; ++p) {
                    if (m.pieces[p].off != pos || m.pieces[p].len <= 0) { msg = "pieces do not tile " + g.name(cx.nodes[i]); return "V7"; }
                    pos += m.pieces[p].len;
                }
                if (pos != cx.len(i)) { msg = "pieces do not tile " + g.name(cx.nodes[i]); return "V7"; }
            }
            // a chain's mapping as a whole (its new and compatible pieces as it placed them, after
            // pinning): what it represents by the target is at most target / i, plus G for the
            // differences shorter than G that are absorbed.  (New pieces alone can be a few bp
            // left over next to compatible ones, with an absorbed indel.)
            std::set<uint32_t> chains;
            for (const auto& kv : m.fates)
                for (const Fate& f : kv.second) chains.insert(f.chain);
            for (uint32_t ci : chains) {
                int64_t removed = 0, target = 0;
                for (const XPiece& p : cands[ci].final_pieces) {
                    if (p.state == P_DROP) continue;
                    removed += p.b - p.a;
                    target += p.tb - p.ta;
                }
                if ((double)removed > (double)target / opt.ident + (double)opt.G + 1e-9) {
                    msg = strf("unit %u represents %lld bp by %lld bp of target (more than target / %.3f + G)", cands[ci].unit,
                               (long long)removed, (long long)target, opt.ident);
                    return "V7";
                }
            }
        }
        return "";
    }

    // --check: V2/V3 over every anchored alt-only path of the site (none repeats a handle);
    // skipped when the site has more than 5,000 of them
    std::string check_paths(const Model& m, std::string& msg) {
        const SiteGraph& sg = sd.sg;
        const uint32_t nh = (uint32_t)sg.n_handles();
        const size_t max_paths = 5000;
        const uint64_t max_steps = 5000000;
        std::vector<std::vector<Handle>> paths;
        uint64_t steps = 0;
        bool over = false;
        std::vector<uint8_t> onpath(nh, 0);
        std::vector<uint32_t> path;
        struct Fr { uint32_t h, k; };
        for (uint32_t h0 = 0; h0 < nh && !over; ++h0) {
            if (sg.dep_off[h0] == sg.dep_off[h0 + 1]) continue;
            auto emit = [&](uint32_t h) {
                for (uint32_t a = sg.arr_off[h]; a < sg.arr_off[h + 1] && !over; ++a)
                    for (uint32_t d = sg.dep_off[h0]; d < sg.dep_off[h0 + 1]; ++d) {
                        if (paths.size() >= max_paths) { over = true; return; }
                        std::vector<Handle> w;
                        w.push_back(sg.dep[d].ref);
                        for (uint32_t x : path) w.push_back(sg.global(x));
                        w.push_back(sg.arr[a].ref);
                        paths.push_back(std::move(w));
                    }
            };
            std::vector<Fr> stk;
            stk.push_back(Fr{h0, sg.succ_off[h0]});
            path.assign(1, h0);
            onpath[h0] = 1;
            emit(h0);
            while (!stk.empty() && !over) {
                if (++steps > max_steps) { over = true; break; }
                uint32_t h = stk.back().h;
                if (stk.back().k < sg.succ_off[h + 1]) {
                    uint32_t w = sg.succ[stk.back().k++];
                    if (onpath[w]) continue;
                    onpath[w] = 1;
                    path.push_back(w);
                    stk.push_back(Fr{w, sg.succ_off[w]});
                    emit(w);
                } else {
                    onpath[h] = 0;
                    path.pop_back();
                    stk.pop_back();
                }
            }
            for (uint32_t x : path) onpath[x] = 0;
            path.clear();
        }
        if (over) {
            ZDEBUG("--check: site %s has more than %zu anchored paths; the path check is skipped", cx.s.label(g).c_str(), max_paths);
            return "";
        }
        for (const std::vector<Handle>& w : paths) {
            bool touched = false;
            for (Handle h : w) {
                uint32_t i = cx.idx(handle_node(h));
                if (i == NONE || m.changed[i]) { touched = true; break; }
            }
            if (!touched) continue;
            std::string r = check_image(cx, m, w, msg);
            if (!r.empty()) {
                std::string ws;
                for (Handle h : w) ws += g.handle_str(h) + " ";
                msg = "--check path " + ws + ": " + msg;
                return r;
            }
        }
        return "";
    }
};

// ================================================================ SitePlanner

namespace {
// adds the time of a scope to `acc` (the planner's own time)
struct OwnClock {
    double& acc;
    double t0;
    explicit OwnClock(double& a) : acc(a), t0(steady_seconds()) {}
    ~OwnClock() { acc += steady_seconds() - t0; }
    OwnClock(const OwnClock&) = delete;
    OwnClock& operator=(const OwnClock&) = delete;
};
} // namespace

SitePlanner::SitePlanner(const Graph& g, const SiteData& sd, const Options& opt) : impl_(new Impl(g, sd, opt)) {}

SitePlanner::~SitePlanner() {}

void SitePlanner::accept_reference(std::vector<AlignResult>& res, std::vector<ReportRow>& rows) {
    Impl& I = *impl_;
    OwnClock own(I.own_s);
    const SiteData& sd = I.sd;
    if (rows.size() != sd.units.size()) fail(EXIT_INVARIANT, "accept_reference: one report row per unit is required");
    if (res.size() != sd.units.size()) res.assign(sd.units.size(), AlignResult());
    std::string label = sd.site->label(I.g);
    for (const std::pair<std::string, uint32_t>& d : debug_drop_links())
        if (d.first == label) I.drop_chain_unit = d.second;
    // node support: the number of excursion owners (creator runs, or observed contigs) whose
    // excursion passes through each node; the prototype's order weight (one unit per run)
    std::unordered_map<NodeId, uint32_t> support;
    for (const Unit& u : sd.units) {
        std::vector<NodeId> ns;
        for (Handle h : u.exc.alts) ns.push_back(handle_node(h));
        std::sort(ns.begin(), ns.end());
        ns.erase(std::unique(ns.begin(), ns.end()), ns.end());
        uint32_t w = std::max<uint32_t>(1, (uint32_t)u.exc.owners.size());
        for (NodeId n : ns) support[n] += w;
    }
    std::vector<uint32_t> order;
    for (const Unit& u : sd.units) {
        const AlignResult& r = res[u.id];
        if (!r.outcome.empty() || r.chain.parts.empty()) continue;
        bool inside = true;
        for (Handle h : u.exc.alts)
            if (I.cx.idx(handle_node(h)) == NONE) inside = false;
        if (!inside) continue;
        Cand c;
        c.unit = u.id;
        c.space = 0;
        c.exc = &u.exc;
        c.walk = u.exc.alts;
        c.off.assign(1, 0);
        for (Handle h : c.walk) c.off.push_back(c.off.back() + I.g.len(handle_node(h)));
        c.wlo = u.wlo;
        c.whi = u.whi;
        c.parts = r.chain.parts;
        c.weight = u.exc.weight;
        c.row = &rows[u.id];
        c.drop_link = u.id == I.drop_chain_unit;
        prepare_chain(I.cx, c, I.boundaries0);
        // kept bp for the order: pieces that pass (B) (and, with observed walks, (G)) after
        // snapping and have an aligned base, each weighted by its node's support
        c.order_bp = 0;
        c.order_w = 0;
        for (const XPiece& p : c.base)
            if (p.state == P_NEW && p.aligned > 0 && I.allowed_B(c, p) && !I.reads_G(c, p)) {
                c.order_bp += p.qb - p.qa;
                c.order_w += (p.qb - p.qa) * (int64_t)std::max<uint32_t>(1, support[I.cx.nodes[p.i]]);
            }
        order.push_back((uint32_t)I.cands.size());
        I.cands.push_back(std::move(c));
    }
    const bool complete = sd.complete_walks;
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        const Cand& x = I.cands[a];
        const Cand& y = I.cands[b];
        if (complete && x.weight != y.weight) return x.weight > y.weight;
        if (x.order_w != y.order_w) return x.order_w > y.order_w;
        if (x.score != y.score) return x.score > y.score;
        if (excursion_key_less(*x.exc, *y.exc)) return true;
        if (excursion_key_less(*y.exc, *x.exc)) return false;
        return x.unit < y.unit;
    });
    if (!order.empty() && !I.baseline()) {
        for (uint32_t ci : order) {
            I.seq.push_back(ci);
            I.cands[ci].outcome = "reverted:model";
            ++I.broken_cands;
        }
        return;
    }
    for (uint32_t ci : order) {
        I.seq.push_back(ci);
        I.accept_one(ci);
    }
}

void SitePlanner::accept_alt(std::vector<AltCandidate>& acs) {
    Impl& I = *impl_;
    if (acs.empty()) return;
    OwnClock own(I.own_s);
    if (!I.baseline()) {
        for (AltCandidate& a : acs) {
            if (!a.res || !a.res->outcome.empty() || a.res->chain.parts.empty()) continue;
            if (a.row) a.row->outcome = "reverted:model";
            ++I.broken_cands;
        }
        return;
    }
    for (AltCandidate& a : acs) {
        if (!a.res || !a.res->outcome.empty() || a.res->chain.parts.empty()) continue;
        // the target space: u, the representative branch, v
        TSpace S;
        S.ref = false;
        S.walk.push_back(a.u);
        S.c.push_back(-I.g.len(handle_node(a.u)));
        int64_t pos = 0;
        for (Handle h : a.target) {
            S.walk.push_back(h);
            S.c.push_back(pos);
            pos += I.g.len(handle_node(h));
        }
        S.walk.push_back(a.v);
        S.c.push_back(pos);
        S.c.push_back(pos + I.g.len(handle_node(a.v)));
        // the members of one group share their representative: an identical target walk (u, branch,
        // v) is one space, so (A) and (C) compare them like chains of one reference window
        auto sit = I.alt_space.find(S.walk);
        const uint32_t space = sit == I.alt_space.end() ? NONE : sit->second;
        // a node zipped or used as a target is never again a query, a target or a bound; the nodes
        // of a representative stay the target of its own group's further members
        auto free_target = [&](NodeId n) {
            if (!I.used.count(n)) return true;
            auto it = I.alt_target.find(n);
            return space != NONE && it != I.alt_target.end() && it->second == space;
        };
        bool ok = I.cx.idx(handle_node(a.u)) != NONE && I.cx.idx(handle_node(a.v)) != NONE && !a.target.empty() &&
                  !I.used.count(handle_node(a.u)) && !I.used.count(handle_node(a.v));
        for (Handle h : a.target)
            if (I.cx.idx(handle_node(h)) == NONE || !free_target(handle_node(h))) ok = false;
        for (Handle h : a.query)
            if (I.cx.idx(handle_node(h)) == NONE || I.used.count(handle_node(h))) ok = false;
        if (!ok) {
            if (a.row) a.row->outcome = "shares-nodes";
            continue;
        }
        Cand c;
        c.alt = true;
        c.unit = a.key;
        if (space == NONE) {
            c.space = (uint32_t)I.cx.spaces.size();
            I.alt_space[S.walk] = c.space;
            I.cx.spaces.push_back(std::move(S));
        } else {
            c.space = space;
        }
        const TSpace& CS = I.cx.spaces[c.space];
        std::vector<int64_t> bnd(CS.c.begin() + 1, CS.c.end() - 1);
        c.walk = a.query;
        c.off.assign(1, 0);
        for (Handle h : c.walk) c.off.push_back(c.off.back() + I.g.len(handle_node(h)));
        c.wlo = 0;
        c.whi = pos;
        c.parts = a.res->chain.parts;
        c.row = a.row;
        prepare_chain(I.cx, c, bnd);
        uint32_t ci = (uint32_t)I.cands.size();
        I.cands.push_back(std::move(c));
        I.seq.push_back(ci);
        I.accept_one(ci);
        if (I.cands[ci].outcome == "zipped")
            for (Handle h : a.target) {
                I.used.insert(handle_node(h));
                I.alt_target.insert(std::make_pair(handle_node(h), I.cands[ci].space));
            }
    }
}

bool SitePlanner::node_used(NodeId n) const { return impl_->used.count(n) != 0; }

SitePlan SitePlanner::finish(std::vector<ReportRow>& rows, std::vector<ReportRow>& extra_rows) {
    (void)rows;
    (void)extra_rows;
    Impl& I = *impl_;
    const double f0 = steady_seconds();
    const Graph& g = I.g;
    const Ctx& cx = I.cx;
    SitePlan plan;
    plan.site = I.sd.site->index;
    plan.sn = I.sd.site->sn;
    plan.lo = I.sd.site->lo;
    if (I.finished) fail(EXIT_INVARIANT, "SitePlanner::finish called twice");
    I.finished = true;
    std::string label = cx.s.label(g);
    Model m;
    if (!I.st.accepted.empty()) {
        build_model(cx, I.st.fates, I.st.space_cuts, I.drop_index(I.st), m);
        std::string msg;
        std::string bad = I.validate(m, msg);
        if (!bad.empty()) {
            ZLOG("site %s: the edit fails %s (%s); re-planning chain by chain", label.c_str(), bad.c_str(), msg.c_str());
            // re-plan incrementally in acceptance order: one bad chain costs only itself
            State s2;
            s2.space_cuts.assign(cx.spaces.size(), std::set<int64_t>());
            s2.sccs = I.base_sccs;
            for (uint32_t ci : I.seq) {
                Cand& c = I.cands[ci];
                FateMap add;
                std::set<int64_t> cuts;
                std::vector<std::vector<uint32_t>> after;
                std::string r = I.try_accept(ci, s2, add, cuts, after);
                if (!r.empty()) { c.outcome = r; continue; }
                State s3 = s2;
                I.commit(ci, s3, add, cuts, after);
                Model m3;
                build_model(cx, s3.fates, s3.space_cuts, I.drop_index(s3), m3);
                std::string msg3;
                std::string f3 = I.validate(m3, msg3);
                if (!f3.empty()) {
                    c.outcome = "reverted:" + f3;
                    ++plan.reverted;
                    ZLOG("site %s: the chain of unit %u is reverted by %s: %s", label.c_str(), c.unit, f3.c_str(), msg3.c_str());
                    continue;
                }
                s2 = std::move(s3);
                c.outcome = "zipped";
            }
            I.st = std::move(s2);
            I.used.clear();
            for (const auto& kv : I.st.fates) I.used.insert(cx.nodes[kv.first]);
            build_model(cx, I.st.fates, I.st.space_cuts, I.drop_index(I.st), m);
            std::string msg4;
            std::string f4 = I.st.accepted.empty() ? std::string() : I.validate(m, msg4);
            if (!f4.empty())
                fail(EXIT_INVARIANT, strf("site %s: the re-planned edit still fails %s: %s", label.c_str(), f4.c_str(), msg4.c_str()));
        }
    }
    // report rows (dropped: the align module's island-below-b entries stay first).  kept_bp is the
    // alt bp the chain removes from the graph: its new pieces.  A piece compatible with an earlier
    // chain's zip is counted by that chain, and a chain that is not zipped removes nothing (0), so
    // the column sums to the bp removed (plan.zipped_bp, checked below).  The chain's own total
    // (new and compatible pieces, the bp the >= b rule reads) stays in the --dump edit tables.
    for (Cand& c : I.cands) {
        if (!c.row) continue;
        ReportRow& r = *c.row;
        r.outcome = c.outcome;
        r.kept_bp = c.kept_bp < 0 ? -1 : (c.outcome == "zipped" ? c.new_bp : 0);
        std::string d = r.dropped;
        size_t shown = 0;
        for (const std::string& x : c.dropped) {
            if (shown == 50) { d += strf(";+%zu more", c.dropped.size() - 50); break; }
            if (!d.empty()) d.push_back(';');
            d += x;
            ++shown;
        }
        r.dropped = d;
    }
    // rule (G)'s refusals (observed walks), as the candidates last tried them
    for (const Cand& c : I.cands) {
        uint32_t n = 0;
        for (const XPiece& p : c.final_pieces)
            if (p.state == P_DROP && p.why == "(G)") { ++n; plan.g_bp += p.b - p.a; }
        plan.g_pieces += n;
        if (n) ++plan.g_chains;
    }
    // --dump: the site's decisions
    if (!I.opt.dump_dir.empty() && !I.cands.empty()) {
        std::string dir = I.opt.dump_dir + "/edit";
        if (!make_dirs(dir)) fail(EXIT_IO, "cannot create " + dir);
        AtomicWriter w(dir + "/" + label + ".tsv");
        w.write("#chain\tunit\tpass\toutcome\tchain_bp\tnew_bp\torder_bp\torder_w\tscore\tdropped\n");
        for (uint32_t ci : I.seq) {
            const Cand& c = I.cands[ci];
            w.write(strf("chain\t%u\t%s\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%zu\n", c.unit, c.alt ? "alt" : "ref", c.outcome.c_str(),
                         (long long)c.kept_bp, (long long)c.new_bp, (long long)c.order_bp, (long long)c.order_w, (long long)c.score,
                         c.dropped.size()));
        }
        w.write("#piece\tunit\tnode\ta\tb\tqa\tqb\tta\ttb\trel\tstate\n");
        for (const Cand& c : I.cands)
            for (const XPiece& p : c.final_pieces)
                w.write(strf("piece\t%u\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%c\t%s\n", c.unit, g.name(cx.nodes[p.i]).c_str(), (long long)p.a,
                             (long long)p.b, (long long)p.qa, (long long)p.qb, (long long)p.ta, (long long)p.tb, p.rel,
                             p.state == P_NEW ? "new" : p.state == P_COMPAT ? "compatible" : p.why.c_str()));
        w.commit();
    }
    // the plan
    {
        const double secs = I.own_s + (steady_seconds() - f0);
        if (secs > 5 || log_level() > 1)
            ZLOG("site %s: edit planning %.1f s (%zu candidate chain(s), %zu site node(s), %llu (D) model(s), %llu validation(s))", label.c_str(),
                 secs, I.cands.size(), cx.nodes.size(), (unsigned long long)I.n_models, (unsigned long long)I.n_validations);
    }
    plan.accepted = (uint32_t)I.st.accepted.size();
    plan.reverted += I.broken_cands;
    if (I.st.accepted.empty()) return plan;
    for (uint32_t i = 0; i < (uint32_t)cx.nodes.size(); ++i) {
        if (!m.changed[i]) continue;
        PlanNode pn;
        pn.node = cx.nodes[i];
        for (uint32_t p = m.first[i]; p < m.first[i + 1]; ++p) {
            pn.off.push_back(m.pieces[p].off);
            pn.len.push_back(m.pieces[p].len);
            pn.keep.push_back(m.pieces[p].kept ? 1 : 0);
            if (!m.pieces[p].kept) plan.zipped_bp += m.pieces[p].len;
        }
        if (pn.off.size() > 1) ++plan.cut_nodes;
        plan.nodes.push_back(std::move(pn));
    }
    for (uint32_t li : cx.links) {
        const Link& l = g.links[li];
        uint32_t a = cx.idx(side_node(l.a)), b = cx.idx(side_node(l.b));
        if ((a != NONE && m.changed[a]) || (b != NONE && m.changed[b])) plan.removed_links.push_back(li);
    }
    auto pend = [&](uint32_t side) {
        PlanEnd e;
        const MPiece& mp = m.pieces[side >> 1];
        e.node = cx.nodes[mp.i];
        e.off = mp.off;
        e.right = (side & 1u) != 0;
        return e;
    };
    for (const MLink& l : m.links) {
        const MPiece& pa = m.pieces[l.sa >> 1];
        const MPiece& pb = m.pieces[l.sb >> 1];
        if (l.cls == 0 && !m.changed[pa.i] && !m.changed[pb.i]) continue;
        PlanLink pl;
        pl.a = pend(l.sa);
        pl.b = pend(l.sb);
        pl.sr = l.sr;
        pl.src = l.src;
        pl.mode = l.cls;
        plan.links.push_back(pl);
    }
    plan.links_added = (uint32_t)plan.links.size();
    plan.links_removed = (uint32_t)plan.removed_links.size();
    int64_t new_sum = 0;
    for (uint32_t ci : I.st.accepted) {
        const Cand& c = I.cands[ci];
        new_sum += c.new_bp;
        if (c.new_bp > 0) ++plan.accepted_new;
        for (const XPiece& p : c.final_pieces) {
            if (p.state != P_NEW) continue;
            ZippedPiece z;
            z.unit = c.unit;
            z.alt = c.alt;
            z.node = cx.nodes[p.i];
            z.a = p.a;
            z.b = p.b;
            z.ta = p.ta;
            z.tb = p.tb;
            z.rel = p.rel;
            plan.zipped.push_back(z);
            plan.target_bp += p.tb - p.ta;
        }
    }
    // the report's kept_bp (new pieces of the zipped chains) must sum to the bp the edit removes
    if (new_sum != plan.zipped_bp)
        fail(EXIT_INVARIANT, strf("internal: site %s: the zipped chains' new pieces hold %lld bp, the edit removes %lld bp", label.c_str(),
                                  (long long)new_sum, (long long)plan.zipped_bp));
    return plan;
}

// ================================================================ apply

void apply_plans(Graph& g, std::vector<SitePlan>& plans, const Options& opt) {
    if (opt.id_base >= 0 && opt.id_base <= g.max_id)
        fail(EXIT_INPUT, strf("--id-base %lld is not above the largest input id s%lld", (long long)opt.id_base, (long long)g.max_id));
    const int64_t base = opt.id_base >= 0 ? opt.id_base : g.max_id + 1;
    int64_t next = base;
    std::vector<size_t> order(plans.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const SitePlan& x = plans[a];
        const SitePlan& y = plans[b];
        if (x.sn != y.sn) return x.sn < y.sn;
        if (x.lo != y.lo) return x.lo < y.lo;
        return x.site < y.site;
    });
    // ids in (SN, lo, original id, offset) order; (node, piece offset) -> new node
    std::unordered_map<NodeId, std::vector<std::pair<int64_t, NodeId>>> pieces;
    uint64_t n_new = 0, n_del = 0, n_links = 0, n_rm = 0;
    for (size_t oi : order) {
        SitePlan& p = plans[oi];
        std::sort(p.nodes.begin(), p.nodes.end(), [&](const PlanNode& a, const PlanNode& b) { return g.nodes[a.node].id < g.nodes[b.node].id; });
        for (const PlanNode& pn : p.nodes) {
            if (pieces.count(pn.node)) fail(EXIT_INVARIANT, "apply: " + g.name(pn.node) + " is edited by two sites");
            std::vector<std::pair<int64_t, NodeId>>& v = pieces[pn.node];
            for (size_t k = 0; k < pn.off.size(); ++k) {
                if (!pn.keep[k]) { v.push_back(std::make_pair(pn.off[k], NONE)); continue; }
                Node nd;
                {
                    const Node& parent = g.nodes[pn.node];     // not held across add_node
                    nd.id = next++;
                    nd.seq = parent.seq.substr((size_t)pn.off[k], (size_t)pn.len[k]);
                    nd.sn = parent.sn;
                    nd.sr = parent.sr;
                    nd.so = parent.so + pn.off[k];
                    nd.tags = make_node_tags(pn.len[k], g.sn_names[parent.sn], nd.so, parent.sr);
                }
                NodeId nid = g.add_node(std::move(nd));
                v.push_back(std::make_pair(pn.off[k], nid));
                ++n_new;
            }
            g.nodes[pn.node].deleted = true;
            ++n_del;
        }
    }
    auto resolve = [&](const PlanEnd& e) -> Side {
        auto it = pieces.find(e.node);
        NodeId n = e.node;
        if (it != pieces.end()) {
            n = NONE;
            for (const auto& pr : it->second)
                if (pr.first == e.off) n = pr.second;
            if (n == NONE) fail(EXIT_INVARIANT, strf("apply: no live piece of %s at %lld", g.name(e.node).c_str(), (long long)e.off));
        } else if (e.off != 0) {
            fail(EXIT_INVARIANT, strf("apply: %s is not cut but a link names its offset %lld", g.name(e.node).c_str(), (long long)e.off));
        }
        return e.right ? right_side(n) : left_side(n);
    };
    for (size_t oi : order) {
        SitePlan& p = plans[oi];
        for (uint32_t li : p.removed_links) {
            g.links[li].deleted = true;
            ++n_rm;
        }
        for (const PlanLink& pl : p.links) {
            Link l;
            Side a = resolve(pl.a), b = resolve(pl.b);
            if (pl.mode == 2) {
                // canonical reading of a translated link: the one with more '+' signs (a forward
                // link is written a+ b+), then the lower id first, then (one node) its R side first
                int pab = (side_is_right(a) ? 1 : 0) + (side_is_right(b) ? 0 : 1);
                int pba = (side_is_right(b) ? 1 : 0) + (side_is_right(a) ? 0 : 1);
                int64_t ia = g.nodes[side_node(a)].id, ib = g.nodes[side_node(b)].id;
                if (pba > pab || (pba == pab && (ia > ib || (ia == ib && !side_is_right(a) && side_is_right(b))))) std::swap(a, b);
            }
            l.a = a;
            l.b = b;
            l.sr = pl.sr;
            if (pl.src != NONE) {
                l.overlap = g.links[pl.src].overlap;
                l.tags = g.links[pl.src].tags;
            } else {
                l.overlap = "0M";
                l.tags = make_link_tags(pl.sr);
            }
            g.add_link(std::move(l));
            ++n_links;
        }
    }
    if (n_del || n_links)
        ZLOG("edit: %llu node(s) cut or zipped, %llu new piece(s) s%lld..s%lld; %llu link(s) replaced by %llu",
             (unsigned long long)n_del, (unsigned long long)n_new, (long long)base, (long long)(next - 1), (unsigned long long)n_rm,
             (unsigned long long)n_links);
}

} // namespace zip
