/*
  zip_tandem.cpp -- --tandem (see zip_tandem.hpp): tandem copies beside their insertion point become
  back links over the reference copy.
*/
#include "zip_tandem.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace zip {

void TandemStats::add(const TandemStats& o) {
    candidates += o.candidates;
    copies += o.copies;
    right += o.right;
    left += o.left;
    rotated += o.rotated;
    components += o.components;
    refused_mixed += o.refused_mixed;
    refused_edited += o.refused_edited;
    refused_unwalked += o.refused_unwalked;
    units += o.units;
    nodes += o.nodes;
    bp += o.bp;
    for (const auto& kv : o.why) why[kv.first] += kv.second;
    back_links += o.back_links;
    skip_links += o.skip_links;
    existing += o.existing;
    equivalent += o.equivalent;
    snapped += o.snapped;
    splits += o.splits;
}

double global_identity(const std::string& a, const std::string& b) {
    const size_t n = a.size(), m = b.size();
    if (n == 0 && m == 0) return 1.0;
    std::vector<uint32_t> prev(m + 1), cur(m + 1);
    for (size_t j = 0; j <= m; ++j) prev[j] = (uint32_t)j;
    for (size_t i = 1; i <= n; ++i) {
        cur[0] = (uint32_t)i;
        for (size_t j = 1; j <= m; ++j) {
            uint32_t sub = prev[j - 1] + (std::toupper((unsigned char)a[i - 1]) == std::toupper((unsigned char)b[j - 1]) ? 0u : 1u);
            cur[j] = std::min(sub, std::min(prev[j], cur[j - 1]) + 1u);
        }
        std::swap(prev, cur);
    }
    return 1.0 - (double)prev[m] / (double)std::max(n, m);
}

double banded_identity(const std::string& a, const std::string& b, int64_t w) {
    const int64_t n = (int64_t)a.size(), m = (int64_t)b.size();
    if (n == 0 && m == 0) return 1.0;
    if (std::llabs(n - m) > w) return 0.0;
    const uint32_t INF = UINT32_MAX / 2;
    // row i holds columns j in [i - w, i + w], at index j - i + w
    std::vector<uint32_t> prev((size_t)(2 * w + 1), INF), cur((size_t)(2 * w + 1), INF);
    for (int64_t j = 0; j <= std::min(m, w); ++j) prev[(size_t)(j + w)] = (uint32_t)j;
    for (int64_t i = 1; i <= n; ++i) {
        std::fill(cur.begin(), cur.end(), INF);
        for (int64_t j = std::max<int64_t>(0, i - w); j <= std::min(m, i + w); ++j) {
            const size_t k = (size_t)(j - i + w);
            uint32_t v = INF;
            if (j == 0) v = (uint32_t)i;
            else {
                // (i-1, j-1) is at the same index in prev; (i-1, j) one to the right; (i, j-1) one to the left
                if (prev[k] < INF) v = prev[k] + (std::toupper((unsigned char)a[(size_t)(i - 1)]) == std::toupper((unsigned char)b[(size_t)(j - 1)]) ? 0u : 1u);
                if (k + 1 < prev.size() && prev[k + 1] < INF) v = std::min(v, prev[k + 1] + 1u);
                if (k > 0 && cur[k - 1] < INF) v = std::min(v, cur[k - 1] + 1u);
            }
            cur[k] = v;
        }
        std::swap(prev, cur);
    }
    const uint32_t d = prev[(size_t)(m - n + w)];
    return d >= INF ? 0.0 : 1.0 - (double)d / (double)std::max(n, m);
}

namespace {

// a verified copy of one candidate unit
struct Copy {
    bool ok = false;
    std::string why;                   // refusal reason when !ok
    std::string shape;                 // right | left | rotated
    int64_t loop = 0;                  // bp of reference the loop runs over
    std::vector<std::pair<int64_t, int64_t>> back;   // (from, to), from > to: the back link
    std::vector<std::pair<int64_t, int64_t>> skip;   // (p1, p2): the window skip, when the haplotype lacks the window
    std::string extra;                 // the rotation's other block (head or tail), formatted like a chain block
    int64_t extra_q = 0;               // ... and its query bp
};

bool near(int64_t x, int64_t y, int64_t tol) { return x >= y - tol && x <= y + tol; }

// candidates whose query is more repetitive than this (repeat_fraction) are not aligned
const double TANDEM_MAX_SELF_REPEAT = 0.5;

// the target position of query position x on a '+' record (left-snap: a D run starting exactly at x
// is not consumed), qs <= x <= qe
int64_t project_record(const PafRecord& r, int64_t x) {
    int64_t q = r.qs, t = r.ts;
    for (const CigarOp& op : r.ops) {
        if (q >= x) break;
        if (op.op == 'D') { t += op.len; continue; }
        int64_t take = std::min<int64_t>(op.len, x - q);
        q += take;
        if (op.op != 'I') t += take;
    }
    return t;
}

// The '+' record (identity >= -i, aligned query >= --min-piece) that carries the rest of a rotated
// copy: the tail from about q to the query end, ending at p1 or p2 (tail), or the head from the
// query start to about q, starting at p1 or p2 (head).  It may run past q (a similar next unit);
// the caller clips it there.  -1 if none; phase gets the end it touches.
int rest_record(const AlignResult& r, const Options& opt, bool tail, int64_t q, int64_t qlen, int64_t p1, int64_t p2, int64_t& phase) {
    int best = -1;
    int64_t best_eq = -1;
    for (size_t k = 0; k < r.records.size(); ++k) {
        const PafRecord& x = r.records[k];
        if (x.strand != '+' || x.gci < opt.ident || x.qe - x.qs < opt.min_piece) continue;
        bool fits = tail ? (x.qs <= q + opt.G && x.qe >= qlen - opt.G && (near(x.te, p1, opt.G) || near(x.te, p2, opt.G)))
                         : (x.qs <= opt.G && x.qe >= q - opt.G && (near(x.ts, p2, opt.G) || near(x.ts, p1, opt.G)));
        if (fits && x.n_eq > best_eq) {
            best = (int)k;
            best_eq = x.n_eq;
            phase = tail ? (near(x.te, p1, opt.G) ? p1 : p2) : (near(x.ts, p2, opt.G) ? p2 : p1);
        }
    }
    return best;
}

// query bases of [qa, qb) that a '+' record aligns: '=' and 'X' columns, and I runs shorter than G
// (an I run >= G is unaligned query, as the chain's blocks count it)
int64_t aligned_in(const PafRecord& r, int64_t qa, int64_t qb, int64_t G) {
    int64_t q = r.qs, n = 0;
    for (const CigarOp& op : r.ops) {
        if (op.op == 'D') continue;
        const int64_t a = std::max(q, qa), b = std::min<int64_t>(q + op.len, qb);
        if (b > a && (op.op != 'I' || op.len < G)) n += b - a;
        q += op.len;
    }
    return n;
}

std::string block_str(int64_t qs, int64_t qe, int64_t ts, int64_t te, double id) {
    return strf("%lld-%lld:%lld-%lld:+:%.4f", (long long)qs, (long long)qe, (long long)ts, (long long)te, id);
}

// Step 3: the copy shape of one aligned candidate.  With p1 = end(dep) and p2 = start(arr), the
// haplotype is ref[..p1) + U + ref[p2..).  A copy touches the insertion point with its start (U
// starts at s, s = p1 or p2) or its end (U ends at e, e = p1 or p2).  The links that spell the
// haplotype through the reference copy:
//   right    U = ref[s, d), d > p2           the back link d -> p2
//   left     U = ref[c, e), c < p1           the back link p1 -> c
//   rotated  U = ref[s, d) + ref[c, e)       the back link d -> c
// plus, for an F unit (p1 < p2), the skip p1 -> p2 when s = p2 or e = p1: then the haplotype lacks
// the window, which is otherwise inside the copy.  Kind I (p1 = p2) never needs a skip.  At most
// max(2G, 5%) of the query may be unaligned (novel sequence inside a copy is not dropped).
Copy copy_shape(const Graph& g, int32_t sn, const Unit& u, const AlignResult& r, int64_t flo, int64_t fhi, const Options& opt) {
    Copy cp;
    const int64_t p1 = u.wlo, p2 = u.whi, Q = u.query_bp, G = opt.G;
    const bool F = p1 < p2;
    if (!r.outcome.empty()) { cp.why = r.outcome; return cp; }
    const ChainResult& ch = r.chain;
    if (ch.parts.empty()) { cp.why = "no-chain"; return cp; }
    if (ch.strands != "+" || ch.frame != 'F') { cp.why = "strand"; return cp; }
    const int64_t q0 = ch.parts.front().qs, q1 = ch.parts.back().qe;
    int64_t t0 = ch.parts.front().ts, t1 = ch.parts.back().te;
    for (const ChainPart& p : ch.parts) {
        t0 = std::min(t0, p.ts);
        t1 = std::max(t1, p.te);
    }
    // the phase a chain end touches: the nearer of p1 and p2 within G.  Both can be within G when the
    // window is <= 2G; a tie goes to the reading that keeps the window (start p1, end p2: no skip).
    auto phase = [&](int64_t t, int64_t keep, int64_t other) -> int64_t {
        const bool a = near(t, keep, G), b = near(t, other, G);
        if (a && b) return std::llabs(t - other) < std::llabs(t - keep) ? other : keep;
        return a ? keep : b ? other : -1;
    };
    const int64_t ps = q0 <= G ? phase(t0, p1, p2) : -1;
    const int64_t pe = q1 >= Q - G ? phase(t1, p2, p1) : -1;
    int64_t c = -1, d = -1, s = -1, e = -1;
    if (q0 <= G && q1 >= Q - G && ps >= 0) {
        cp.shape = "right";
        s = ps;
        d = t1;
        if (d <= p2) { cp.why = "in-window"; return cp; }
        cp.back.push_back(std::make_pair(d, p2));
        cp.loop = d - p2;
    } else if (q0 <= G && q1 >= Q - G && pe >= 0) {
        cp.shape = "left";
        e = pe;
        c = t0;
        if (c >= p1) { cp.why = "in-window"; return cp; }
        cp.back.push_back(std::make_pair(p1, c));
        cp.loop = p1 - c;
    } else if (ps >= 0 || pe >= 0) {
        cp.shape = "rotated";
        const std::string q = g.spell(u.exc.alts);
        if (ps >= 0) {
            // the chain is the head ref[s, d); the tail U[q1, Q) is ref[c, e)
            s = ps;
            d = t1;
            int k = rest_record(r, opt, true, q1, Q, p1, p2, e);
            if (k >= 0) {
                const PafRecord& x = r.records[(size_t)k];
                c = x.qs >= q1 ? x.ts : project_record(x, q1);
                cp.extra = block_str(std::max(q1, x.qs), x.qe, c, x.te, x.gci);
                cp.extra_q = aligned_in(x, std::max(q1, x.qs), x.qe, G);
            } else if (Q - q1 < opt.min_piece) {
                for (int64_t ph : {p1, p2}) {
                    if (ph - (Q - q1) < flo) continue;
                    double id = global_identity(q.substr((size_t)q1), g.ref_seq(sn, ph - (Q - q1), ph));
                    if (id >= opt.ident) {
                        e = ph;
                        c = ph - (Q - q1);
                        cp.extra = block_str(q1, Q, c, ph, id);
                        cp.extra_q = Q - q1;
                        break;
                    }
                }
            }
            if (c < 0) { cp.why = "rotation-tail"; return cp; }
        } else {
            // the chain is the tail ref[c, e); the head U[0, q0) is ref[s, d)
            e = pe;
            c = t0;
            int k = rest_record(r, opt, false, q0, Q, p1, p2, s);
            if (k >= 0) {
                const PafRecord& x = r.records[(size_t)k];
                d = x.qe <= q0 ? x.te : project_record(x, q0);
                cp.extra = block_str(x.qs, std::min(q0, x.qe), x.ts, d, x.gci);
                cp.extra_q = aligned_in(x, x.qs, std::min(q0, x.qe), G);
            } else if (q0 < opt.min_piece) {
                for (int64_t ph : {p2, p1}) {
                    if (ph + q0 > fhi) continue;
                    double id = global_identity(q.substr(0, (size_t)q0), g.ref_seq(sn, ph, ph + q0));
                    if (id >= opt.ident) {
                        s = ph;
                        d = ph + q0;
                        cp.extra = block_str(0, q0, ph, d, id);
                        cp.extra_q = q0;
                        break;
                    }
                }
            }
            if (d < 0) { cp.why = "rotation-head"; return cp; }
        }
        if (!(c < e && s < d && c < d)) { cp.why = "rotation-order"; return cp; }
        cp.back.push_back(std::make_pair(d, c));
        cp.loop = d - c;
    } else {
        cp.why = "not-adjacent";
        return cp;
    }
    if (F && (s == p2 || e == p1)) cp.skip.push_back(std::make_pair(p1, p2));
    const int64_t aligned = ch.aligned_q + cp.extra_q;
    if (Q - aligned > std::max<int64_t>(2 * G, Q / 20)) { cp.why = "novel"; return cp; }
    if (cp.loop < opt.b) { cp.why = "short-loop"; return cp; }
    // the copy is at least half of the loop: a few blocks spread over a long stretch are not a copy
    if (cp.loop > 2 * Q) { cp.why = "sparse"; return cp; }
    cp.ok = true;
    return cp;
}

} // namespace

void tandem_site(const Graph& g, const SiteData& sd, const SitePlan& plan, Aligner& aligner, const Options& opt,
                 std::vector<ReportRow>& unit_rows, std::vector<TandemEdit>& edits, TandemStats& st) {
    const Site& s = *sd.site;
    const std::string label = s.label(g);
    const std::vector<NodeId>& ref = g.ref_nodes[(size_t)s.sn];
    if (ref.empty()) return;
    const int64_t sn_lo = g.start(ref.front()), sn_hi = g.end(ref.back());
    // step 1: candidates and their flanks
    std::vector<uint32_t> cand;
    std::vector<AlignJob> jobs;
    std::vector<std::pair<int64_t, int64_t>> flank;
    for (const Unit& u : sd.units) {
        if (u.kind != Kind::I && u.kind != Kind::F) continue;
        if (u.outcome != "no-window" && u.outcome != "infeasible") continue;
        if (u.query_bp < opt.b || u.query_bp - u.window_bp() < opt.b) continue;
        ++st.candidates;
        ReportRow& row = unit_rows[u.id];
        if (handle_rev(u.exc.dep) || handle_rev(u.exc.arr)) {
            row.dropped = "tandem:orientation";
            ++st.why["orientation"];
            continue;
        }
        // one copy of a long unit has little internal repetition; a query that is mostly its own repeats
        // is a VNTR expansion or several units at once (not handled), and aligning it is the costly case
        const double rf = repeat_fraction(g.spell(u.exc.alts));
        if (rf > TANDEM_MAX_SELF_REPEAT) {
            row.dropped = strf("tandem:self-repeat:%.3f", rf);
            ++st.why["self-repeat"];
            continue;
        }
        const int64_t X = u.query_bp + u.query_bp / 10 + 1000;
        const int64_t flo = std::max(sn_lo, u.wlo - X), fhi = std::min(sn_hi, u.whi + X);
        if (u.query_bp > opt.max_pair || fhi - flo > opt.max_pair) {
            row.dropped = "tandem:pair-too-big";
            ++st.why["pair-too-big"];
            continue;
        }
        AlignJob j;
        j.key = strf("%s:tandem:%u", label.c_str(), u.id);
        j.query = u.exc.alts;
        j.target.ref = true;
        j.target.sn = s.sn;
        j.target.lo = flo;
        j.target.hi = fhi;
        j.keep_records = true;
        jobs.push_back(std::move(j));
        cand.push_back(u.id);
        flank.push_back(std::make_pair(flo, fhi));
    }
    if (jobs.empty()) return;
    // step 2
    std::vector<AlignResult> res;
    aligner.align_jobs(g, jobs, res);
    // step 3
    std::map<uint32_t, Copy> copies;
    for (size_t k = 0; k < cand.size(); ++k) {
        const Unit& u = sd.units[cand[k]];
        Copy cp = copy_shape(g, s.sn, u, res[k], flank[k].first, flank[k].second, opt);
        if (!cp.ok) {
            unit_rows[u.id].dropped = "tandem:" + cp.why;
            ++st.why[cp.why];
            continue;
        }
        ++st.copies;
        if (cp.shape == "right") ++st.right;
        else if (cp.shape == "left") ++st.left;
        else ++st.rotated;
        // the row's align fields come from the copy's chain; outcome and frame are set below
        ReportRow& row = unit_rows[u.id];
        const std::string outcome = row.outcome;
        fill_align_row(row, res[k], opt.G);
        if (!cp.extra.empty()) {
            // the rotation's other block, in query order
            row.blocks = cp.shape == "rotated" && res[k].chain.parts.front().qs <= opt.G ? row.blocks + ";" + cp.extra : cp.extra + ";" + row.blocks;
        }
        row.frame = "T";
        row.label = cp.shape;
        row.outcome = outcome;
        copies[u.id] = cp;
    }
    if (copies.empty()) return;
    // step 4: components
    // nodes the site's plan touches: cut or zipped nodes, and the ends of the links it adds or removes
    // (an uncut alt-vs-alt representative keeps the links of its member's kept pieces)
    std::unordered_set<NodeId> edited;
    for (const PlanNode& pn : plan.nodes) edited.insert(pn.node);
    for (const PlanLink& pl : plan.links) {
        edited.insert(pl.a.node);
        edited.insert(pl.b.node);
    }
    for (uint32_t li : plan.removed_links) {
        edited.insert(side_node(g.links[li].a));
        edited.insert(side_node(g.links[li].b));
    }
    std::unordered_map<NodeId, std::vector<uint32_t>> units_of;   // alt node -> units walking it
    for (const Unit& u : sd.units)
        for (Handle h : u.exc.alts) {
            std::vector<uint32_t>& v = units_of[handle_node(h)];
            if (v.empty() || v.back() != u.id) v.push_back(u.id);
        }
    std::set<NodeId> done;
    for (const auto& kv : copies) {
        const Unit& u0 = sd.units[kv.first];
        if (done.count(handle_node(u0.exc.alts.front()))) continue;
        // the component: alt nodes of the site reachable over alt-alt links (input links)
        std::vector<NodeId> comp;
        std::set<NodeId> seen;
        std::vector<NodeId> stack;
        for (Handle h : u0.exc.alts)
            if (seen.insert(handle_node(h)).second) stack.push_back(handle_node(h));
        while (!stack.empty()) {
            NodeId n = stack.back();
            stack.pop_back();
            comp.push_back(n);
            for (int o = 0; o < 2; ++o)
                for (const Edge& e : g.out(make_handle(n, o == 1))) {
                    NodeId m = handle_node(e.to);
                    if (s.is_alt(m) && seen.insert(m).second) stack.push_back(m);
                }
        }
        std::sort(comp.begin(), comp.end());
        done.insert(comp.begin(), comp.end());
        // every unit through the component, and whether each is a verified copy
        std::set<uint32_t> us;
        for (NodeId n : comp) {
            auto it = units_of.find(n);
            if (it != units_of.end()) us.insert(it->second.begin(), it->second.end());
        }
        // every unit through it is a verified copy, every node of it is on one (a node no walk covers,
        // e.g. a partial GAF record's or an unanchored run's, may carry novel sequence), and the plan
        // left it alone
        bool mixed = false, unwalked = false, touched = false;
        for (uint32_t id : us)
            if (!copies.count(id)) mixed = true;
        for (NodeId n : comp) {
            if (!units_of.count(n)) unwalked = true;
            if (edited.count(n)) touched = true;
        }
        if (mixed || unwalked || touched) {
            const std::string why = mixed ? "mixed" : unwalked ? "unwalked" : "edited";
            for (uint32_t id : us)
                if (copies.count(id)) unit_rows[id].outcome = "tandem-refused:" + why;
            ++(mixed ? st.refused_mixed : unwalked ? st.refused_unwalked : st.refused_edited);
            continue;
        }
        TandemEdit ed;
        ed.sn = s.sn;
        ed.remove = comp;
        int64_t bp = 0;
        for (NodeId n : comp) {
            bp += g.len(n);
            ed.sr = ed.sr < 0 ? g.nodes[n].sr : std::min(ed.sr, g.nodes[n].sr);
        }
        bool first = true;
        for (uint32_t id : us) {
            const Copy& cp = copies[id];
            for (const auto& bk : cp.back) {
                bool dup = false;
                for (const auto& b : ed.back)
                    if (near(b.first, bk.first, opt.G) && near(b.second, bk.second, opt.G)) dup = true;
                if (!dup) ed.back.push_back(bk);
            }
            for (const auto& sk : cp.skip)
                if (std::find(ed.skip.begin(), ed.skip.end(), sk) == ed.skip.end()) ed.skip.push_back(sk);
            ReportRow& row = unit_rows[id];
            row.outcome = opt.detect_only ? "would-tandem" : "tandem";
            row.kept_bp = first ? bp : 0;
            first = false;
            ++st.units;
        }
        ++st.components;
        st.nodes += comp.size();
        st.bp += bp;
        edits.push_back(std::move(ed));
    }
}

void apply_tandem(Graph& g, const std::vector<const std::vector<TandemEdit>*>& edits, const Options& opt, TandemStats& st) {
    size_t n_edits = 0;
    for (const auto* v : edits) n_edits += v->size();
    if (n_edits == 0) return;
    int64_t next = g.max_id + 1;
    for (const Node& n : g.nodes) next = std::max(next, n.id + 1);
    if (opt.id_base > next) next = opt.id_base;
    // delete the components and every live link with an end on them
    std::vector<uint8_t> gone(g.nodes.size(), 0);
    for (const auto* v : edits)
        for (const TandemEdit& ed : *v)
            for (NodeId n : ed.remove) {
                if (g.nodes[n].deleted) fail(EXIT_INVARIANT, "tandem: " + g.name(n) + " is already deleted");
                g.nodes[n].deleted = true;
                gone[n] = 1;
            }
    std::unordered_set<LinkKey, LinkKeyHash> keys;
    uint64_t n_rm = 0;
    for (Link& l : g.links) {
        if (l.deleted) continue;
        if (gone[side_node(l.a)] || gone[side_node(l.b)]) {
            l.deleted = true;
            ++n_rm;
            continue;
        }
        keys.insert(link_key(l.a, l.b));
    }
    // live reference nodes of each SN, by SO (the edit split some; reference pieces are never deleted)
    std::map<int32_t, std::vector<NodeId>> live;
    auto ref_of = [&](int32_t sn) -> std::vector<NodeId>& {
        auto it = live.find(sn);
        if (it != live.end()) return it->second;
        std::vector<NodeId>& v = live[sn];
        for (NodeId n = 0; n < (NodeId)g.nodes.size(); ++n)
            if (!g.nodes[n].deleted && g.nodes[n].sr == 0 && g.nodes[n].sn == sn) v.push_back(n);
        std::sort(v.begin(), v.end(), [&](NodeId a, NodeId b) { return g.start(a) < g.start(b); });
        return v;
    };
    // the live reference node covering pos (start <= pos < end), as an index into ref_of(sn)
    auto covering = [&](std::vector<NodeId>& v, int64_t pos) -> size_t {
        size_t k = (size_t)(std::upper_bound(v.begin(), v.end(), pos, [&](int64_t p, NodeId n) { return p < g.start(n); }) - v.begin());
        if (k == 0) fail(EXIT_INVARIANT, strf("tandem: no reference node at %lld", (long long)pos));
        return k - 1;
    };
    // a reference boundary at pos: snap to one within G (if allowed), or split the node there
    auto boundary = [&](int32_t sn, int64_t pos, bool snap_ok) -> int64_t {
        std::vector<NodeId>& v = ref_of(sn);
        size_t k = covering(v, pos);
        const NodeId n = v[k];
        if (g.start(n) == pos || g.end(n) == pos) return pos;
        const int64_t lo = g.start(n), hi = g.end(n);
        int64_t snap = pos - lo <= hi - pos ? lo : hi;
        if (snap_ok && std::llabs(snap - pos) <= opt.G) {
            ++st.snapped;
            return snap;
        }
        // split n at pos into two new pieces; the parent is deleted and its links move to the pieces
        Node a, b;
        {
            const Node& p = g.nodes[n];
            const int64_t off = pos - p.so;
            a.id = next++;
            a.seq = p.seq.substr(0, (size_t)off);
            a.sn = p.sn;
            a.sr = p.sr;
            a.so = p.so;
            a.tags = make_node_tags(off, g.sn_names[(size_t)p.sn], a.so, a.sr);
            b.id = next++;
            b.seq = p.seq.substr((size_t)off);
            b.sn = p.sn;
            b.sr = p.sr;
            b.so = pos;
            b.tags = make_node_tags((int64_t)b.seq.size(), g.sn_names[(size_t)p.sn], b.so, b.sr);
        }
        NodeId na = g.add_node(std::move(a));
        NodeId nb = g.add_node(std::move(b));
        g.nodes[n].deleted = true;
        for (Link& l : g.links) {
            if (l.deleted) continue;
            Side oa = l.a, ob = l.b;
            auto move = [&](Side x) -> Side {
                if (side_node(x) != n) return x;
                return side_is_right(x) ? right_side(nb) : left_side(na);
            };
            l.a = move(l.a);
            l.b = move(l.b);
            if (l.a != oa || l.b != ob) {
                keys.erase(link_key(oa, ob));
                keys.insert(link_key(l.a, l.b));
            }
        }
        Link l;
        l.a = right_side(na);
        l.b = left_side(nb);
        l.sr = 0;
        l.tags = make_link_tags(0);
        keys.insert(link_key(l.a, l.b));
        g.add_link(std::move(l));
        v[k] = na;
        v.insert(v.begin() + (long)k + 1, nb);
        ++st.splits;
        return pos;
    };
    // the R side of the live reference node ending at pos / the L side of the one starting at pos
    auto side_at = [&](int32_t sn, int64_t pos, bool right) -> Side {
        std::vector<NodeId>& v = ref_of(sn);
        size_t k = covering(v, right ? pos - 1 : pos);
        if (right ? g.end(v[k]) != pos : g.start(v[k]) != pos)
            fail(EXIT_INVARIANT, strf("tandem: no reference node %s at %lld", right ? "ends" : "starts", (long long)pos));
        return right ? right_side(v[k]) : left_side(v[k]);
    };
    auto add = [&](Side a, Side b, int32_t sr, uint64_t& counter) {
        if (keys.count(link_key(a, b))) { ++st.existing; return; }
        Link l;
        l.a = a;
        l.b = b;
        l.sr = sr;
        l.tags = make_link_tags(sr);
        keys.insert(link_key(a, b));
        g.add_link(std::move(l));
        ++counter;
    };
    const int64_t base = next;
    // Back links, shortest loop first.  A loop with the same period as one already added (within
    // max(G, 1%)), overlapping it by at least half the period, is the same route when the stretch
    // between their starts matches the stretch one period later (identity >= -i): copies of one array
    // unit placed at different phases.  It is not added, so each copy number keeps one route.
    struct Back { int32_t sn; int64_t from, to; int32_t sr; };
    std::vector<Back> backs;
    for (const auto* v : edits)
        for (const TandemEdit& ed : *v)
            for (const auto& bk : ed.back) backs.push_back(Back{ed.sn, bk.first, bk.second, ed.sr});
    std::stable_sort(backs.begin(), backs.end(), [](const Back& a, const Back& b) {
        if (a.sn != b.sn) return a.sn < b.sn;
        if (a.from - a.to != b.from - b.to) return a.from - a.to < b.from - b.to;
        return a.to < b.to;
    });
    std::vector<Back> kept;
    for (const Back& bk : backs) {
        const int64_t P = bk.from - bk.to;
        bool same = false;
        for (const Back& o : kept) {
            const int64_t p = o.from - o.to;
            if (o.sn != bk.sn || std::llabs(P - p) > std::max<int64_t>(opt.G, p / 100)) continue;
            const int64_t ov = std::min(bk.from, o.from) - std::max(bk.to, o.to);
            if (ov < p / 2) continue;
            const int64_t lo = std::min(bk.to, o.to), hi = std::max(bk.to, o.to);
            if (hi == lo) { same = true; break; }
            const std::string x = g.ref_seq(bk.sn, lo, hi), y = g.ref_seq(bk.sn, lo + p, hi + p);
            if (banded_identity(x, y, std::max<int64_t>(64, (hi - lo) / 50)) >= opt.ident) { same = true; break; }
        }
        if (same) { ++st.equivalent; continue; }
        kept.push_back(bk);
        // snapping must not shorten the loop to nothing (possible with a small -b)
        int64_t d = boundary(bk.sn, bk.from, true), c = boundary(bk.sn, bk.to, true);
        if (d <= c) {
            d = boundary(bk.sn, bk.from, false);
            c = boundary(bk.sn, bk.to, false);
        }
        add(side_at(bk.sn, d, true), side_at(bk.sn, c, false), bk.sr, st.back_links);
    }
    for (const auto* v : edits)
        for (const TandemEdit& ed : *v)
            for (const auto& sk : ed.skip) add(side_at(ed.sn, sk.first, true), side_at(ed.sn, sk.second, false), ed.sr, st.skip_links);
    ZLOG("tandem edit: %zu component(s): %llu alt node(s) and %llu link(s) removed; %llu back link(s) (%llu more were the same loop as one "
         "of them), %llu skip link(s) added (%llu already there); loop ends: %llu snapped within -G, %llu reference node(s) split%s",
         n_edits, (unsigned long long)st.nodes, (unsigned long long)n_rm, (unsigned long long)st.back_links, (unsigned long long)st.equivalent,
         (unsigned long long)st.skip_links, (unsigned long long)st.existing, (unsigned long long)st.snapped,
         (unsigned long long)st.splits, st.splits ? strf(" (new pieces s%lld..s%lld)", (long long)base, (long long)(next - 1)).c_str() : "");
}

} // namespace zip
