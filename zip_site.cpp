/*
  zip_site.cpp -- spec steps 1-3: snarls (real JSON), sites (side-aware interior, leaks), runs,
  creator excursions, witness fallback, canonical frame, deduplication, kinds, Allowed, pre-filter
  and caps.

  Replaces rgfa-collapse's snarl reader (403-437; it counted "name" substrings), its interior
  (480-499; it ignored link sides) and its traversal/component code (514-530, 647-692).
  The creator excursion is a port of the design prototype's rf_excursion and
  Allowed of all_paths_allowed, computed exactly over the SCC condensation instead of by
  worklist relaxation.
*/
#include "zip_site.hpp"

#include <numeric>
#include <queue>
#include <unordered_set>

namespace zip {

// ================================================================ names

const char* site_status_name(SiteStatus s) {
    switch (s) {
        case SiteStatus::OK: return "ok";
        case SiteStatus::NO_ALT: return "no-alt";
        case SiteStatus::LEAK: return "site:leak";
        case SiteStatus::REF_MISMATCH: return "site:ref-mismatch";
        case SiteStatus::TOO_BIG: return "site:too-big";
    }
    return "?";
}

const char* allowed_status_name(AllowedStatus s) {
    switch (s) {
        case AllowedStatus::NONE: return "none";
        case AllowedStatus::OK: return "ok";
        case AllowedStatus::EMPTY: return "empty";
        case AllowedStatus::BLOCKED: return "blocked";
    }
    return "?";
}

const char* src_name(Src s) {
    switch (s) {
        case Src::CREATOR: return "creator";
        case Src::GAF: return "gaf";
        case Src::WITNESS: return "witness";
    }
    return "?";
}

const char* kind_name(Kind k) {
    switch (k) {
        case Kind::F: return "F";
        case Kind::I: return "I";
        case Kind::BK: return "BK";
        case Kind::J: return "J";
    }
    return "?";
}

void WalkStats::add(const WalkStats& o) {
    runs += o.runs; creator_ok += o.creator_ok; fallback += o.fallback; witness_dag += o.witness_dag;
    witness_bubble += o.witness_bubble; unanchored += o.unanchored;
    fail_creator_links += o.fail_creator_links; fail_own_ambiguous += o.fail_own_ambiguous;
    fail_creator_edge += o.fail_creator_edge; fail_cycle += o.fail_cycle; fail_other += o.fail_other;
    excursions_raw += o.excursions_raw; excursions += o.excursions; merged += o.merged;
}

// ================================================================ JSON (one value per line)

namespace {

struct JVal {
    enum Type { NUL, BOOL, NUM, STR, ARR, OBJ } type = NUL;
    bool b = false;
    double num = 0;
    std::string str;
    std::vector<JVal> arr;
    std::vector<std::pair<std::string, JVal>> obj;
    const JVal* get(const char* k) const {
        for (const auto& kv : obj)
            if (kv.first == k) return &kv.second;
        return nullptr;
    }
};

class JsonParser {
public:
    JsonParser(const char* p, const char* e) : p_(p), e_(e) {}
    bool parse(JVal& v, std::string& err) {
        skip_ws();
        if (!value(v, 0)) { err = err_; return false; }
        skip_ws();
        if (p_ != e_) { err = "trailing characters"; return false; }
        return true;
    }
private:
    const char* p_;
    const char* e_;
    std::string err_;
    static bool digit(char c) { return c >= '0' && c <= '9'; }
    void skip_ws() {
        while (p_ < e_ && (*p_ == ' ' || *p_ == '\t' || *p_ == '\n' || *p_ == '\r')) ++p_;
    }
    bool bad(const char* m) {
        if (err_.empty()) err_ = m;
        return false;
    }
    bool value(JVal& v, int depth) {
        if (depth > 64) return bad("nesting too deep");
        skip_ws();
        if (p_ >= e_) return bad("unexpected end of line");
        char c = *p_;
        if (c == '{') return object(v, depth);
        if (c == '[') return array(v, depth);
        if (c == '"') { v.type = JVal::STR; return string(v.str); }
        if (c == 't') return literal("true", v, JVal::BOOL, true);
        if (c == 'f') return literal("false", v, JVal::BOOL, false);
        if (c == 'n') return literal("null", v, JVal::NUL, false);
        if (c == '-' || digit(c)) return number(v);
        return bad("unexpected character");
    }
    bool literal(const char* lit, JVal& v, JVal::Type t, bool b) {
        size_t n = strlen(lit);
        if ((size_t)(e_ - p_) < n || strncmp(p_, lit, n) != 0) return bad("bad literal");
        p_ += n;
        v.type = t;
        v.b = b;
        return true;
    }
    bool number(JVal& v) {
        const char* s = p_;
        if (p_ < e_ && *p_ == '-') ++p_;
        if (p_ >= e_) return bad("bad number");
        if (*p_ == '0') ++p_;
        else if (digit(*p_)) { while (p_ < e_ && digit(*p_)) ++p_; }
        else return bad("bad number");
        if (p_ < e_ && *p_ == '.') {
            ++p_;
            if (p_ >= e_ || !digit(*p_)) return bad("bad number");
            while (p_ < e_ && digit(*p_)) ++p_;
        }
        if (p_ < e_ && (*p_ == 'e' || *p_ == 'E')) {
            ++p_;
            if (p_ < e_ && (*p_ == '+' || *p_ == '-')) ++p_;
            if (p_ >= e_ || !digit(*p_)) return bad("bad number");
            while (p_ < e_ && digit(*p_)) ++p_;
        }
        v.type = JVal::NUM;
        v.num = strtod(std::string(s, p_).c_str(), nullptr);
        return true;
    }
    bool hex4(unsigned& u) {
        if (e_ - p_ < 4) return bad("bad \\u escape");
        u = 0;
        for (int i = 0; i < 4; ++i) {
            char c = *p_++;
            u <<= 4;
            if (c >= '0' && c <= '9') u |= (unsigned)(c - '0');
            else if (c >= 'a' && c <= 'f') u |= (unsigned)(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') u |= (unsigned)(c - 'A' + 10);
            else return bad("bad \\u escape");
        }
        return true;
    }
    static void utf8(std::string& out, unsigned cp) {
        if (cp < 0x80) out.push_back((char)cp);
        else if (cp < 0x800) { out.push_back((char)(0xC0 | (cp >> 6))); out.push_back((char)(0x80 | (cp & 0x3F))); }
        else if (cp < 0x10000) {
            out.push_back((char)(0xE0 | (cp >> 12))); out.push_back((char)(0x80 | ((cp >> 6) & 0x3F)));
            out.push_back((char)(0x80 | (cp & 0x3F)));
        } else {
            out.push_back((char)(0xF0 | (cp >> 18))); out.push_back((char)(0x80 | ((cp >> 12) & 0x3F)));
            out.push_back((char)(0x80 | ((cp >> 6) & 0x3F))); out.push_back((char)(0x80 | (cp & 0x3F)));
        }
    }
    bool string(std::string& out) {
        ++p_;   // opening quote
        while (true) {
            if (p_ >= e_) return bad("unterminated string");
            unsigned char c = (unsigned char)*p_++;
            if (c == '"') return true;
            if (c < 0x20) return bad("control character in string");
            if (c != '\\') { out.push_back((char)c); continue; }
            if (p_ >= e_) return bad("bad escape");
            char x = *p_++;
            switch (x) {
                case '"': out.push_back('"'); break;
                case '\\': out.push_back('\\'); break;
                case '/': out.push_back('/'); break;
                case 'b': out.push_back('\b'); break;
                case 'f': out.push_back('\f'); break;
                case 'n': out.push_back('\n'); break;
                case 'r': out.push_back('\r'); break;
                case 't': out.push_back('\t'); break;
                case 'u': {
                    unsigned u = 0;
                    if (!hex4(u)) return false;
                    if (u >= 0xD800 && u <= 0xDBFF) {
                        unsigned lo = 0;
                        if (e_ - p_ < 2 || p_[0] != '\\' || p_[1] != 'u') return bad("unpaired surrogate");
                        p_ += 2;
                        if (!hex4(lo) || lo < 0xDC00 || lo > 0xDFFF) return bad("unpaired surrogate");
                        u = 0x10000 + ((u - 0xD800) << 10) + (lo - 0xDC00);
                    } else if (u >= 0xDC00 && u <= 0xDFFF) {
                        return bad("unpaired surrogate");
                    }
                    utf8(out, u);
                    break;
                }
                default: return bad("bad escape");
            }
        }
    }
    bool array(JVal& v, int depth) {
        ++p_;
        v.type = JVal::ARR;
        skip_ws();
        if (p_ < e_ && *p_ == ']') { ++p_; return true; }
        while (true) {
            v.arr.emplace_back();
            if (!value(v.arr.back(), depth + 1)) return false;
            skip_ws();
            if (p_ >= e_) return bad("unterminated array");
            if (*p_ == ',') { ++p_; continue; }
            if (*p_ == ']') { ++p_; return true; }
            return bad("expected ',' or ']'");
        }
    }
    bool object(JVal& v, int depth) {
        ++p_;
        v.type = JVal::OBJ;
        skip_ws();
        if (p_ < e_ && *p_ == '}') { ++p_; return true; }
        while (true) {
            skip_ws();
            if (p_ >= e_ || *p_ != '"') return bad("expected a key");
            std::string k;
            if (!string(k)) return false;
            skip_ws();
            if (p_ >= e_ || *p_ != ':') return bad("expected ':'");
            ++p_;
            v.obj.emplace_back(std::move(k), JVal());
            if (!value(v.obj.back().second, depth + 1)) return false;
            skip_ws();
            if (p_ >= e_) return bad("unterminated object");
            if (*p_ == ',') { ++p_; continue; }
            if (*p_ == '}') { ++p_; return true; }
            return bad("expected ',' or '}'");
        }
    }
};

} // namespace

std::vector<Snarl> read_snarls(const std::string& path, const Graph& g) {
    static const char* hint = "snarls must come from: vg snarls -n -P <ref> <graph> | vg view -Rj -";
    LineReader in(path);
    std::vector<Snarl> out;
    std::string line;
    while (in.next(line)) {
        uint64_t ln = in.line_number();
        bool blank = true;
        for (char c : line)
            if (c != ' ' && c != '\t' && c != '\r') { blank = false; break; }
        if (blank) continue;
        JVal v;
        std::string err;
        JsonParser jp(line.data(), line.data() + line.size());
        if (!jp.parse(v, err))
            fail(EXIT_INPUT, strf("%s line %llu is not valid JSON (%s); %s", path.c_str(), (unsigned long long)ln, err.c_str(), hint));
        if (v.type != JVal::OBJ)
            fail(EXIT_INPUT, strf("%s line %llu is not a JSON object; %s", path.c_str(), (unsigned long long)ln, hint));
        Snarl s;
        s.line = ln;
        s.nested = v.get("parent") != nullptr;
        for (int k = 0; k < 2; ++k) {
            const char* key = k == 0 ? "start" : "end";
            const JVal* b = v.get(key);
            if (!b || b->type != JVal::OBJ)
                fail(EXIT_INPUT, strf("%s line %llu has no \"%s\" object; %s", path.c_str(), (unsigned long long)ln, key, hint));
            const JVal* nm = b->get("name");
            if (!nm || nm->type != JVal::STR)
                fail(EXIT_INPUT, strf("%s line %llu: \"%s\" has no \"name\" (run vg snarls with -n); %s", path.c_str(),
                                      (unsigned long long)ln, key, hint));
            bool back = false;
            const JVal* bw = b->get("backward");
            if (bw) {
                if (bw->type != JVal::BOOL)
                    fail(EXIT_INPUT, strf("%s line %llu: \"%s\".\"backward\" is not a boolean", path.c_str(), (unsigned long long)ln, key));
                back = bw->b;
            }
            NodeId n = g.find_name(nm->str);
            if (n == NONE)
                fail(EXIT_INPUT, strf("%s line %llu names %s, which is not a segment of the graph: the snarls belong to another graph",
                                      path.c_str(), (unsigned long long)ln, nm->str.c_str()));
            if (k == 0) { s.start = n; s.start_back = back; }
            else { s.end = n; s.end_back = back; }
        }
        out.push_back(s);
    }
    in.close();
    return out;
}

// ================================================================ sites

namespace {

struct FloodScratch {
    std::vector<uint32_t> stamp;
    uint32_t gen = 0;
    std::vector<NodeId> stack;
};

// Side-aware interior: a flood from both inner sides, A.R and B.L, that never enters A or B -- what
// vg calls the snarl's contents, including a tip or a fold-back that hangs off B.L only.  The site
// leaks if any inner side (A.R, B.L or an interior side) links to an outer side (A.L or B.R).
// Past `cap` nodes the site is site:too-big, but the flood goes on, so that the skip row carries
// the site's whole alt bp and every node counts as covered by a snarl (build_sites clears the
// interior of a skipped site); a leaking site stops there (its count is then a lower bound: a leak
// can reach the whole graph).
void flood_site(const Graph& g, Site& s, int64_t cap, FloodScratch& fs) {
    if (fs.stamp.size() != g.nodes.size()) { fs.stamp.assign(g.nodes.size(), 0); fs.gen = 0; }
    if (++fs.gen == 0) { std::fill(fs.stamp.begin(), fs.stamp.end(), 0); fs.gen = 1; }
    const uint32_t gen = fs.gen;
    std::vector<NodeId>& st = fs.stack;
    st.clear();
    s.interior.clear();
    bool leak = false, capped = false;
    int64_t alt_bp = 0;
    const Handle a_enter_left = make_handle(s.A, false);   // entering A through A.L
    const Handle b_enter_right = make_handle(s.B, true);   // entering B through B.R
    auto visit = [&](NodeId y) {
        if (fs.stamp[y] == gen) return;
        fs.stamp[y] = gen;
        st.push_back(y);
        s.interior.push_back(y);
        if (!g.is_ref(y)) alt_bp += g.len(y);
        if ((int64_t)s.interior.size() > cap) capped = true;
    };
    // the steps out of handle h: an outer side of A or B is a leak, an inner one is not followed
    auto step = [&](Handle h) {
        for (const Edge& e : g.out(h)) {
            NodeId y = handle_node(e.to);
            if (y == s.A || y == s.B) {
                if (e.to == a_enter_left || e.to == b_enter_right) leak = true;
                continue;
            }
            visit(y);
        }
    };
    step(make_handle(s.A, false));   // links on A.R
    step(make_handle(s.B, true));    // links on B.L
    while (!st.empty() && !(leak && capped)) {
        NodeId x = st.back();
        st.pop_back();
        step(make_handle(x, false));
        step(make_handle(x, true));
    }
    std::sort(s.interior.begin(), s.interior.end());
    s.n_interior = s.interior.size();
    s.alt_bp = alt_bp;
    if (leak) { s.status = SiteStatus::LEAK; return; }
    if (capped) { s.status = SiteStatus::TOO_BIG; return; }
    // backbone: exactly the SN nodes tiling [lo, hi)
    s.backbone.clear();
    s.alts.clear();
    bool mismatch = false;
    for (NodeId x : s.interior) {
        if (g.is_ref(x)) {
            if (g.nodes[x].sn != s.sn || g.start(x) < s.lo || g.end(x) > s.hi) mismatch = true;
            s.backbone.push_back(x);
        } else {
            s.alts.push_back(x);
        }
    }
    std::sort(s.backbone.begin(), s.backbone.end(), [&](NodeId a, NodeId b) { return g.start(a) < g.start(b); });
    if (!mismatch) {
        const std::vector<NodeId>& rv = g.ref_nodes[s.sn];
        auto it = std::lower_bound(rv.begin(), rv.end(), s.lo, [&](NodeId n, int64_t pos) { return g.start(n) < pos; });
        size_t k = 0;
        for (; it != rv.end() && g.start(*it) < s.hi; ++it, ++k)
            if (k >= s.backbone.size() || s.backbone[k] != *it) { mismatch = true; break; }
        if (k != s.backbone.size()) mismatch = true;
    }
    if (mismatch) { s.status = SiteStatus::REF_MISMATCH; return; }
    s.status = s.alts.empty() ? SiteStatus::NO_ALT : SiteStatus::OK;
}

} // namespace

std::vector<Site> build_sites(const Graph& g, const std::vector<Snarl>& snarls, const Options& opt, SiteStats& st) {
    st = SiteStats();
    st.snarls = snarls.size();
    struct Cand { NodeId A, B; uint64_t line; };
    std::vector<Cand> cand;
    for (const Snarl& sn : snarls) {
        if (sn.nested) { ++st.nested; continue; }
        if (sn.start == sn.end) { ++st.same_node; continue; }
        if (!g.is_ref(sn.start) || !g.is_ref(sn.end)) { ++st.nonref_boundary; continue; }
        if (g.nodes[sn.start].sn != g.nodes[sn.end].sn) { ++st.cross_sn; continue; }
        bool fwd = g.start(sn.start) < g.start(sn.end);
        NodeId A = fwd ? sn.start : sn.end, B = fwd ? sn.end : sn.start;
        bool consistent = fwd ? (!sn.start_back && !sn.end_back) : (sn.start_back && sn.end_back);
        if (!consistent) ++st.orientation_disagrees;
        cand.push_back(Cand{A, B, sn.line});
    }
    std::sort(cand.begin(), cand.end(), [&](const Cand& x, const Cand& y) {
        int32_t sx = g.nodes[x.A].sn, sy = g.nodes[y.A].sn;
        if (sx != sy) return sx < sy;
        if (g.end(x.A) != g.end(y.A)) return g.end(x.A) < g.end(y.A);
        if (x.A != y.A) return x.A < y.A;
        if (x.B != y.B) return x.B < y.B;
        return x.line < y.line;
    });
    std::vector<Site> sites;
    sites.reserve(cand.size());
    for (size_t i = 0; i < cand.size(); ++i) {
        if (i > 0 && cand[i].A == cand[i - 1].A && cand[i].B == cand[i - 1].B) { ++st.duplicate; continue; }
        Site s;
        s.index = (uint32_t)sites.size();
        s.A = cand[i].A;
        s.B = cand[i].B;
        s.sn = g.nodes[s.A].sn;
        s.lo = g.end(s.A);
        s.hi = g.start(s.B);
        s.snarl_line = cand[i].line;
        sites.push_back(std::move(s));
    }
    st.sites = sites.size();
    int nt = std::max(1, opt.threads);
    std::vector<FloodScratch> scratch((size_t)nt);
    parallel_for(sites.size(), nt, [&](size_t i, int t) { flood_site(g, sites[i], opt.max_site_nodes, scratch[(size_t)t]); });
    scratch.clear();
    // every alt node lies in some top-level snarl of a valid decomposition: count those that lie in
    // no site's interior (whatever the site's status), to catch a truncated or partial snarls file
    {
        std::vector<uint8_t> in_site(g.nodes.size(), 0);
        for (const Site& s : sites)
            for (NodeId x : s.interior) in_site[x] = 1;
        for (NodeId n = 0; n < (NodeId)g.n_input_nodes; ++n) {
            if (g.is_ref(n)) continue;
            ++st.alt_nodes;
            st.alt_bp += g.len(n);
            if (!in_site[n]) { ++st.uncovered_nodes; st.uncovered_bp += g.len(n); }
        }
    }
    for (Site& s : sites) {
        switch (s.status) {
            case SiteStatus::OK: ++st.ok; if (s.lo == s.hi) ++st.pure_insertion; break;
            case SiteStatus::NO_ALT: ++st.no_alt; break;
            case SiteStatus::LEAK: ++st.leak; break;
            case SiteStatus::REF_MISMATCH: ++st.ref_mismatch; break;
            case SiteStatus::TOO_BIG: ++st.too_big; break;
        }
        if (s.status != SiteStatus::OK && s.status != SiteStatus::NO_ALT) {
            std::vector<NodeId>().swap(s.interior);
            std::vector<NodeId>().swap(s.backbone);
            std::vector<NodeId>().swap(s.alts);
        }
    }
    // processed sites must be node-disjoint, and no boundary may lie inside another site
    {
        std::vector<uint32_t> owner(g.nodes.size(), NONE);
        for (const Site& s : sites) {
            if (s.status != SiteStatus::OK && s.status != SiteStatus::NO_ALT) continue;
            for (NodeId x : s.interior) {
                if (owner[x] != NONE)
                    fail(EXIT_INPUT, strf("snarls %s and %s overlap at %s: the snarls are not a decomposition of this graph",
                                          sites[owner[x]].label(g).c_str(), s.label(g).c_str(), g.name(x).c_str()));
                owner[x] = s.index;
            }
        }
        for (const Site& s : sites) {
            if (s.status != SiteStatus::OK && s.status != SiteStatus::NO_ALT) continue;
            for (NodeId x : {s.A, s.B})
                if (owner[x] != NONE)
                    fail(EXIT_INPUT, strf("boundary %s of snarl %s lies inside snarl %s: the snarls are not a decomposition of this graph",
                                          g.name(x).c_str(), s.label(g).c_str(), sites[owner[x]].label(g).c_str()));
        }
    }
    // One snarl that is not a snarl of this graph leaks into the sites that share its boundary
    // nodes (their inner sides reach the leaking node), so up to three leaking sites are one
    // problem and are tolerated (each is skipped with a site:leak row); beyond that, more than 1%
    // of sites leaking means the snarls belong to another graph.
    if (st.leak > 3 && st.leak * 100 > st.sites)
        fail(EXIT_INPUT, strf("%llu of %llu sites leak (an inner side links to A.L or B.R); more than 1%% means the snarls "
                              "belong to another graph",
                              (unsigned long long)st.leak, (unsigned long long)st.sites));
    return sites;
}

// ================================================================ the site's alt-handle graph

namespace {

// Tarjan's SCC, iterative.  Component ids come out in reverse topological order.
void tarjan(const SiteGraph& sg, std::vector<uint32_t>& comp, uint32_t& ncomp) {
    const uint32_t n = (uint32_t)sg.n_handles();
    std::vector<uint32_t> index(n, NONE), low(n, 0);
    std::vector<uint8_t> onstack(n, 0);
    std::vector<uint32_t> st;
    struct Frame { uint32_t v, ei; };
    std::vector<Frame> cs;
    uint32_t idx = 0;
    ncomp = 0;
    comp.assign(n, NONE);
    for (uint32_t s0 = 0; s0 < n; ++s0) {
        if (index[s0] != NONE) continue;
        index[s0] = low[s0] = idx++;
        st.push_back(s0);
        onstack[s0] = 1;
        cs.push_back(Frame{s0, sg.succ_off[s0]});
        while (!cs.empty()) {
            uint32_t v = cs.back().v;
            if (cs.back().ei < sg.succ_off[v + 1]) {
                uint32_t w = sg.succ[cs.back().ei++];
                if (index[w] == NONE) {
                    index[w] = low[w] = idx++;
                    st.push_back(w);
                    onstack[w] = 1;
                    cs.push_back(Frame{w, sg.succ_off[w]});
                } else if (onstack[w]) {
                    low[v] = std::min(low[v], index[w]);
                }
            } else {
                if (low[v] == index[v]) {
                    while (true) {
                        uint32_t w = st.back();
                        st.pop_back();
                        onstack[w] = 0;
                        comp[w] = ncomp;
                        if (w == v) break;
                    }
                    ++ncomp;
                }
                cs.pop_back();
                if (!cs.empty()) {
                    uint32_t u = cs.back().v;
                    low[u] = std::min(low[u], low[v]);
                }
            }
        }
    }
}

} // namespace

void build_site_graph(const Graph& g, const Site& s, SiteGraph& sg) {
    sg = SiteGraph();
    sg.alt = s.alts;
    sg.local.reserve(sg.alt.size() * 2);
    for (uint32_t l = 0; l < (uint32_t)sg.alt.size(); ++l) sg.local[sg.alt[l]] = l;
    const uint32_t nh = (uint32_t)sg.n_handles();
    sg.succ_off.assign(nh + 1, 0);
    sg.pred_off.assign(nh + 1, 0);
    sg.dep_off.assign(nh + 1, 0);
    sg.arr_off.assign(nh + 1, 0);
    const Handle a_enter_left = make_handle(s.A, false), b_enter_right = make_handle(s.B, true);
    const Handle b_leave_right = make_handle(s.B, false), a_leave_left = make_handle(s.A, true);
    for (uint32_t lh = 0; lh < nh; ++lh) {
        Handle h = sg.global(lh);
        for (const Edge& e : g.out(h)) {              // successors of h
            NodeId y = handle_node(e.to);
            uint32_t ly = sg.local_node(y);
            if (ly != NONE) { sg.succ.push_back(2 * ly + (handle_rev(e.to) ? 1u : 0u)); continue; }
            if (!s.is_anchor_node(g, y)) continue;     // cannot happen in a non-leaking site
            if (e.to == a_enter_left || e.to == b_enter_right) continue;   // leak sides (excluded by build_sites)
            Anchor a;
            a.ref = e.to;
            a.minus = handle_rev(e.to);
            a.coord = a.minus ? g.end(y) : g.start(y);   // (r,+) enters r.L at start(r); (r,-) enters r.R at end(r)
            sg.arr.push_back(a);
        }
        sg.succ_off[lh + 1] = (uint32_t)sg.succ.size();
        sg.arr_off[lh + 1] = (uint32_t)sg.arr.size();
        for (const Edge& e : g.out(flip(h))) {        // predecessors of h are flip(e.to)
            Handle q = flip(e.to);
            NodeId y = handle_node(q);
            uint32_t ly = sg.local_node(y);
            if (ly != NONE) { sg.pred.push_back(2 * ly + (handle_rev(q) ? 1u : 0u)); continue; }
            if (!s.is_anchor_node(g, y)) continue;
            if (q == b_leave_right || q == a_leave_left) continue;
            Anchor a;
            a.ref = q;
            a.minus = handle_rev(q);
            a.coord = a.minus ? g.start(y) : g.end(y);   // (r,+) leaves r.R at end(r); (r,-) leaves r.L at start(r)
            sg.dep.push_back(a);
        }
        sg.pred_off[lh + 1] = (uint32_t)sg.pred.size();
        sg.dep_off[lh + 1] = (uint32_t)sg.dep.size();
    }
    tarjan(sg, sg.scc, sg.n_scc);
}

// ================================================================ Allowed

void compute_allowed(const Site& s, const SiteGraph& sg, std::vector<AllowedIv>& out) {
    (void)s;
    const uint32_t nc = sg.n_scc;
    const uint32_t nh = (uint32_t)sg.n_handles();
    out.assign(sg.n_nodes(), AllowedIv());
    if (nh == 0) return;
    std::vector<uint32_t> coff(nc + 1, 0), ch(nh);
    for (uint32_t h = 0; h < nh; ++h) coff[sg.scc[h] + 1]++;
    for (uint32_t c = 0; c < nc; ++c) coff[c + 1] += coff[c];
    {
        std::vector<uint32_t> fill(coff.begin(), coff.end() - 1);
        for (uint32_t h = 0; h < nh; ++h) ch[fill[sg.scc[h]]++] = h;
    }
    // forward values from departures, backward values from arrivals (bit 1 = '+', bit 2 = '-')
    struct Fwd { int64_t dplus = INT64_MIN, dminus = INT64_MAX; uint8_t sign = 0; };
    struct Bwd { int64_t aplus = INT64_MAX, aminus = INT64_MIN; uint8_t sign = 0; };
    std::vector<Fwd> F(nc);
    std::vector<Bwd> B(nc);
    for (uint32_t h = 0; h < nh; ++h) {
        Fwd& f = F[sg.scc[h]];
        for (uint32_t i = sg.dep_off[h]; i < sg.dep_off[h + 1]; ++i) {
            const Anchor& a = sg.dep[i];
            if (a.minus) { f.dminus = std::min(f.dminus, a.coord); f.sign |= 2; }
            else { f.dplus = std::max(f.dplus, a.coord); f.sign |= 1; }
        }
        Bwd& b = B[sg.scc[h]];
        for (uint32_t i = sg.arr_off[h]; i < sg.arr_off[h + 1]; ++i) {
            const Anchor& a = sg.arr[i];
            if (a.minus) { b.aminus = std::max(b.aminus, a.coord); b.sign |= 2; }
            else { b.aplus = std::min(b.aplus, a.coord); b.sign |= 1; }
        }
    }
    // an edge c1 -> c2 between components has c2 < c1: sources first is decreasing id
    for (uint32_t c = nc; c-- > 0;) {
        const Fwd fc = F[c];
        if (!fc.sign) continue;
        for (uint32_t k = coff[c]; k < coff[c + 1]; ++k) {
            uint32_t h = ch[k];
            for (uint32_t i = sg.succ_off[h]; i < sg.succ_off[h + 1]; ++i) {
                uint32_t c2 = sg.scc[sg.succ[i]];
                if (c2 == c) continue;
                Fwd& f2 = F[c2];
                f2.dplus = std::max(f2.dplus, fc.dplus);
                f2.dminus = std::min(f2.dminus, fc.dminus);
                f2.sign |= fc.sign;
            }
        }
    }
    for (uint32_t c = 0; c < nc; ++c) {
        Bwd& bc = B[c];
        for (uint32_t k = coff[c]; k < coff[c + 1]; ++k) {
            uint32_t h = ch[k];
            for (uint32_t i = sg.succ_off[h]; i < sg.succ_off[h + 1]; ++i) {
                uint32_t c2 = sg.scc[sg.succ[i]];
                if (c2 == c) continue;
                const Bwd& b2 = B[c2];
                bc.aplus = std::min(bc.aplus, b2.aplus);
                bc.aminus = std::max(bc.aminus, b2.aminus);
                bc.sign |= b2.sign;
            }
        }
    }
    for (uint32_t l = 0; l < (uint32_t)sg.n_nodes(); ++l) {
        uint32_t c = sg.scc[2 * l];
        const Fwd& f = F[c];
        const Bwd& b = B[c];
        AllowedIv& a = out[l];
        if (!f.sign || !b.sign) { a.status = AllowedStatus::NONE; continue; }
        uint8_t sg2 = f.sign | b.sign;
        if (sg2 == 3) { a.status = AllowedStatus::BLOCKED; continue; }
        if (sg2 == 1) { a.lo = f.dplus; a.hi = b.aplus; a.minus = false; }
        else { a.lo = b.aminus; a.hi = f.dminus; a.minus = true; }
        a.status = a.hi > a.lo ? AllowedStatus::OK : AllowedStatus::EMPTY;
    }
}

// ================================================================ runs

Runs build_runs(const Graph& g) {
    Runs R;
    R.run_of.assign(g.nodes.size(), NONE);
    R.pos_in_run.assign(g.nodes.size(), NONE);
    std::vector<std::vector<NodeId>> by_sn(g.sn_names.size());
    for (NodeId n = 0; n < (NodeId)g.nodes.size(); ++n)
        if (!g.is_ref(n)) by_sn[g.nodes[n].sn].push_back(n);
    auto linked = [&](NodeId p, NodeId n) {
        Handle t = make_handle(n, false);
        for (const Edge& e : g.out(make_handle(p, false)))
            if (e.to == t) return true;
        return false;
    };
    std::vector<std::vector<NodeId>> runs;
    for (std::vector<NodeId>& v : by_sn) {
        if (v.empty()) continue;
        std::sort(v.begin(), v.end(), [&](NodeId a, NodeId b) { return g.start(a) != g.start(b) ? g.start(a) < g.start(b) : a < b; });
        std::vector<NodeId> cur{v[0]};
        for (size_t i = 1; i < v.size(); ++i) {
            NodeId p = cur.back(), n = v[i];
            if (g.end(p) == g.start(n) && linked(p, n)) cur.push_back(n);
            else { runs.push_back(cur); cur.assign(1, n); }
        }
        runs.push_back(cur);
    }
    std::sort(runs.begin(), runs.end(), [&](const std::vector<NodeId>& a, const std::vector<NodeId>& b) {
        const Node& x = g.nodes[a[0]];
        const Node& y = g.nodes[b[0]];
        if (x.sr != y.sr) return x.sr < y.sr;
        if (x.sn != y.sn) return x.sn < y.sn;
        if (x.so != y.so) return x.so < y.so;
        return a[0] < b[0];
    });
    for (uint32_t i = 0; i < (uint32_t)runs.size(); ++i)
        for (uint32_t k = 0; k < (uint32_t)runs[i].size(); ++k) {
            R.run_of[runs[i][k]] = i;
            R.pos_in_run[runs[i][k]] = k;
        }
    R.runs.swap(runs);
    return R;
}

std::vector<uint32_t> site_runs(const Runs& runs, const Site& s) {
    std::vector<uint32_t> v;
    for (NodeId n : s.alts)
        if (runs.run_of[n] != NONE) v.push_back(runs.run_of[n]);
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

// ================================================================ excursions

bool owner_less(const Owner& a, const Owner& b) {
    if (a.sr != b.sr) return a.sr < b.sr;
    if (a.contig != b.contig) return a.contig < b.contig;
    if (a.so != b.so) return a.so < b.so;
    if (a.first != b.first) return a.first < b.first;
    return (int)a.src < (int)b.src;
}

bool owner_equal(const Owner& a, const Owner& b) {
    return a.sr == b.sr && a.contig == b.contig && a.so == b.so && a.first == b.first && a.src == b.src;
}

bool canonicalize(Excursion& e) {
    bool dm = handle_rev(e.dep), am = handle_rev(e.arr);
    if (dm && am) {
        std::reverse(e.alts.begin(), e.alts.end());
        for (Handle& h : e.alts) h = flip(h);
        Handle d = flip(e.arr), a = flip(e.dep);
        e.dep = d;
        e.arr = a;
        return true;
    }
    if (dm != am) {
        // a junction can be walked either way; keep the orientation with the smaller key so that
        // both readings of one J excursion deduplicate
        Excursion r;
        r.dep = flip(e.arr);
        r.arr = flip(e.dep);
        r.alts.assign(e.alts.rbegin(), e.alts.rend());
        for (Handle& h : r.alts) h = flip(h);
        if (excursion_key_less(r, e)) {
            e.dep = r.dep;
            e.arr = r.arr;
            e.alts.swap(r.alts);
            return true;
        }
    }
    return false;
}

bool excursion_key_less(const Excursion& a, const Excursion& b) {
    if (a.dep != b.dep) return a.dep < b.dep;
    if (a.arr != b.arr) return a.arr < b.arr;
    return a.alts < b.alts;
}

bool excursion_key_equal(const Excursion& a, const Excursion& b) {
    return a.dep == b.dep && a.arr == b.arr && a.alts == b.alts;
}

std::string excursion_str(const Graph& g, const Excursion& e) {
    std::string s = g.handle_str(e.dep) + ">";
    for (size_t i = 0; i < e.alts.size(); ++i) {
        if (i) s.push_back(',');
        s += g.handle_str(e.alts[i]);
    }
    s += ">" + g.handle_str(e.arr);
    return s;
}

namespace {

Owner run_owner(const Graph& g, const Runs& R, uint32_t ri, Src src) {
    Owner o;
    NodeId f = R.runs[ri].front();
    o.contig = g.sn_names[g.nodes[f].sn];
    o.sr = g.nodes[f].sr;
    o.so = g.nodes[f].so;
    o.first = f;
    o.run = ri;
    o.src = src;
    return o;
}

enum CxResult { CX_OK = 0, CX_CREATOR_LINKS, CX_OWN_AMBIGUOUS, CX_CREATOR_EDGE, CX_CYCLE, CX_OTHER };

// Extend from h (walk orientation) to the first rank-0 handle, by the run's own links (SR == rk)
// or, at a reused node without one, along that node's own run to its creator link.
int extend_run(const Graph& g, const Runs& R, const Site& s, int32_t rk, Handle h, std::vector<Handle>& path, Handle& anchor) {
    std::unordered_set<Handle> seen;
    while (true) {
        NodeId n = handle_node(h);
        if (g.is_ref(n)) {
            if (!s.is_anchor_node(g, n)) return CX_OTHER;
            anchor = h;
            return CX_OK;
        }
        if (R.run_of[n] == NONE || !s.is_alt(n)) return CX_OTHER;
        if (!seen.insert(h).second) return CX_CYCLE;
        path.push_back(h);
        Handle own = NONE;
        int nown = 0;
        for (const Edge& e : g.out(h))
            if (g.link_sr(e.link) == rk) { ++nown; own = e.to; }
        if (nown == 1) { h = own; continue; }
        if (nown > 1) return CX_OWN_AMBIGUOUS;
        const std::vector<NodeId>& rr = R.runs[R.run_of[n]];
        uint32_t j = R.pos_in_run[n];
        int32_t rk2 = g.nodes[n].sr;
        Handle ce = NONE;
        int nce = 0;
        if (!handle_rev(h)) {
            if (j + 1 < rr.size()) { h = make_handle(rr[j + 1], false); continue; }
            for (const Edge& e : g.out(make_handle(rr.back(), false)))
                if (g.link_sr(e.link) == rk2) { ++nce; ce = e.to; }
        } else {
            if (j >= 1) { h = make_handle(rr[j - 1], true); continue; }
            for (const Edge& e : g.out(make_handle(rr.front(), true)))     // leaving rr[0].L
                if (g.link_sr(e.link) == rk2) { ++nce; ce = e.to; }
        }
        if (nce != 1) return CX_CREATOR_EDGE;
        h = ce;
    }
}

int creator_excursion(const Graph& g, const Runs& R, const Site& s, uint32_t ri, Excursion& out) {
    const std::vector<NodeId>& run = R.runs[ri];
    for (NodeId n : run)
        if (!s.is_alt(n)) return CX_OTHER;
    int32_t rk = g.nodes[run[0]].sr;
    Handle first = make_handle(run.front(), false), last = make_handle(run.back(), false);
    Handle ein = NONE, eout = NONE;
    int nin = 0, nout = 0;
    for (const Edge& e : g.out(flip(first)))              // predecessors of first are flip(e.to)
        if (g.link_sr(e.link) == rk) { ++nin; ein = flip(e.to); }
    for (const Edge& e : g.out(last))
        if (g.link_sr(e.link) == rk) { ++nout; eout = e.to; }
    if (nin != 1 || nout != 1) return CX_CREATOR_LINKS;
    std::vector<Handle> fw, bwr;
    Handle a1 = NONE, a0r = NONE;
    int r = extend_run(g, R, s, rk, eout, fw, a1);
    if (r != CX_OK) return r;
    r = extend_run(g, R, s, rk, flip(ein), bwr, a0r);
    if (r != CX_OK) return r;
    out.dep = flip(a0r);
    out.arr = a1;
    out.alts.clear();
    out.alts.reserve(bwr.size() + run.size() + fw.size());
    for (size_t i = bwr.size(); i-- > 0;) out.alts.push_back(flip(bwr[i]));
    for (NodeId n : run) out.alts.push_back(make_handle(n, false));
    out.alts.insert(out.alts.end(), fw.begin(), fw.end());
    return CX_OK;
}

// ---------------------------------------------------------------- witness fallback

class WitnessFinder {
public:
    // covered[l]: local node l lies on an excursion already chosen (weight 0)
    WitnessFinder(const Graph& g, const SiteData& sd, const std::vector<uint8_t>& covered)
        : g_(g), sg_(sd.sg) {
        const uint32_t nh = (uint32_t)sg_.n_handles();
        w_.assign(nh, 0);
        for (uint32_t lh = 0; lh < nh; ++lh) w_[lh] = covered[lh >> 1] ? 0 : g_.len(sg_.alt[lh >> 1]);
        // DFS topological order (reverse post-order); edges against it are back edges and are dropped
        std::vector<uint32_t> post;
        post.reserve(nh);
        std::vector<uint8_t> state(nh, 0);
        struct Frame { uint32_t v, ei; };
        std::vector<Frame> cs;
        for (uint32_t s0 = 0; s0 < nh; ++s0) {
            if (state[s0]) continue;
            state[s0] = 1;
            cs.push_back(Frame{s0, sg_.succ_off[s0]});
            while (!cs.empty()) {
                uint32_t v = cs.back().v;
                if (cs.back().ei < sg_.succ_off[v + 1]) {
                    uint32_t w = sg_.succ[cs.back().ei++];
                    if (!state[w]) { state[w] = 1; cs.push_back(Frame{w, sg_.succ_off[w]}); }
                } else {
                    post.push_back(v);
                    cs.pop_back();
                }
            }
        }
        rank_.assign(nh, 0);
        topo_.assign(post.rbegin(), post.rend());
        for (uint32_t k = 0; k < nh; ++k) rank_[topo_[k]] = k;
        F_.assign(nh, -1); B_.assign(nh, -1);
        Fprev_.assign(nh, NONE); Bnext_.assign(nh, NONE);
        Fanc_.assign(nh, NONE); Banc_.assign(nh, NONE);
        // pass 1: best prefix from a departure to h (including h)
        for (uint32_t k = 0; k < nh; ++k) {
            uint32_t h = topo_[k];
            int64_t best = -1;
            Key bkey{0, 0};
            uint32_t bprev = NONE, banc = NONE;
            for (uint32_t i = sg_.dep_off[h]; i < sg_.dep_off[h + 1]; ++i) {
                Key key{0, sg_.dep[i].ref};
                int64_t cand = w_[h];
                if (cand > best || (cand == best && key < bkey)) { best = cand; bkey = key; bprev = NONE; banc = i; }
            }
            for (uint32_t i = sg_.pred_off[h]; i < sg_.pred_off[h + 1]; ++i) {
                uint32_t q = sg_.pred[i];
                if (rank_[q] >= rank_[h] || F_[q] < 0) continue;
                int64_t cand = F_[q] + w_[h];
                Key key = alt_key(q);
                if (cand > best || (cand == best && key < bkey)) { best = cand; bkey = key; bprev = q; banc = NONE; }
            }
            F_[h] = best; Fprev_[h] = bprev; Fanc_[h] = banc;
        }
        // pass 2: best suffix from h (including h) to an arrival
        for (uint32_t k = nh; k-- > 0;) {
            uint32_t h = topo_[k];
            int64_t best = -1;
            Key bkey{0, 0};
            uint32_t bnext = NONE, banc = NONE;
            for (uint32_t i = sg_.arr_off[h]; i < sg_.arr_off[h + 1]; ++i) {
                Key key{0, sg_.arr[i].ref};
                int64_t cand = w_[h];
                if (cand > best || (cand == best && key < bkey)) { best = cand; bkey = key; bnext = NONE; banc = i; }
            }
            for (uint32_t i = sg_.succ_off[h]; i < sg_.succ_off[h + 1]; ++i) {
                uint32_t q = sg_.succ[i];
                if (rank_[q] <= rank_[h] || B_[q] < 0) continue;
                int64_t cand = B_[q] + w_[h];
                Key key = alt_key(q);
                if (cand > best || (cand == best && key < bkey)) { best = cand; bkey = key; bnext = q; banc = NONE; }
            }
            B_[h] = best; Bnext_[h] = bnext; Banc_[h] = banc;
        }
        dist_.assign(nh, INT64_MAX);
        par_.assign(nh, NONE);
    }

    // 0: no anchored path; 1: DAG witness; 2: shortest anchored bubble
    int find(const std::vector<NodeId>& run, Excursion& out) {
        std::vector<uint32_t> fwd, rev;
        for (NodeId n : run) {
            uint32_t l = sg_.local_node(n);
            if (l == NONE) return 0;
            fwd.push_back(2 * l);
        }
        for (size_t i = fwd.size(); i-- > 0;) rev.push_back(fwd[i] + 1);
        Cand cf, cr;
        bool okf = dag(fwd, cf), okr = dag(rev, cr);
        if (okf || okr) {
            emit((okf && (!okr || cf.total >= cr.total)) ? cf : cr, out);
            return 1;
        }
        if (bubble(fwd, cf) || bubble(rev, cf)) {
            emit(cf, out);
            return 2;
        }
        return 0;
    }

private:
    typedef std::pair<int64_t, uint64_t> Key;   // (SR, then name/handle); anchors (rank 0) sort first
    struct Cand {
        int64_t total = 0;
        Handle dep = 0, arr = 0;
        std::vector<uint32_t> path;   // local handles
    };
    const Graph& g_;
    const SiteGraph& sg_;
    std::vector<int64_t> w_;
    std::vector<uint32_t> rank_, topo_;
    std::vector<int64_t> F_, B_;
    std::vector<uint32_t> Fprev_, Bnext_, Fanc_, Banc_;
    std::vector<int64_t> dist_;
    std::vector<uint32_t> par_;

    Key alt_key(uint32_t lh) const { return Key(g_.nodes[sg_.alt[lh >> 1]].sr, lh); }

    static bool simple(const std::vector<uint32_t>& path) {
        std::vector<uint32_t> v(path);
        std::sort(v.begin(), v.end());
        return std::adjacent_find(v.begin(), v.end()) == v.end();
    }

    bool dag(const std::vector<uint32_t>& R, Cand& c) const {
        uint32_t h0 = R.front(), hk = R.back();
        if (F_[h0] < 0 || B_[hk] < 0) return false;
        int64_t total = F_[h0] + B_[hk];
        if (R.size() == 1) total -= w_[h0];
        else for (size_t i = 1; i + 1 < R.size(); ++i) total += w_[R[i]];
        std::vector<uint32_t> pre;
        uint32_t x = h0;
        while (Fprev_[x] != NONE) { x = Fprev_[x]; pre.push_back(x); }
        Handle dep = sg_.dep[Fanc_[x]].ref;
        std::vector<uint32_t> path(pre.rbegin(), pre.rend());
        path.insert(path.end(), R.begin(), R.end());
        uint32_t y = hk;
        while (Bnext_[y] != NONE) { y = Bnext_[y]; path.push_back(y); }
        Handle arr = sg_.arr[Banc_[y]].ref;
        if (!simple(path)) return false;
        c.total = total;
        c.dep = dep;
        c.arr = arr;
        c.path.swap(path);
        return true;
    }

    // shortest-bp half from `start` to an anchor, avoiding forbidden local nodes; backward = towards a departure
    bool half(uint32_t start, bool backward, const std::vector<uint8_t>& forbid, std::vector<uint32_t>& half_path, Handle& anchor) {
        typedef std::pair<int64_t, std::pair<Key, uint32_t>> QE;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        std::vector<uint32_t> touched;
        dist_[start] = 0;
        par_[start] = NONE;
        touched.push_back(start);
        pq.push(QE(0, std::make_pair(alt_key(start), start)));
        bool found = false;
        uint32_t end_h = NONE;
        while (!pq.empty()) {
            QE top = pq.top();
            pq.pop();
            uint32_t h = top.second.second;
            if (top.first > dist_[h]) continue;
            const uint32_t a0 = backward ? sg_.dep_off[h] : sg_.arr_off[h];
            const uint32_t a1 = backward ? sg_.dep_off[h + 1] : sg_.arr_off[h + 1];
            if (a0 < a1) {
                Handle best = NONE;
                for (uint32_t i = a0; i < a1; ++i) {
                    Handle r = backward ? sg_.dep[i].ref : sg_.arr[i].ref;
                    if (best == NONE || r < best) best = r;
                }
                anchor = best;
                end_h = h;
                found = true;
                break;
            }
            const uint32_t e0 = backward ? sg_.pred_off[h] : sg_.succ_off[h];
            const uint32_t e1 = backward ? sg_.pred_off[h + 1] : sg_.succ_off[h + 1];
            for (uint32_t i = e0; i < e1; ++i) {
                uint32_t q = backward ? sg_.pred[i] : sg_.succ[i];
                if (forbid[q >> 1]) continue;
                int64_t nd = top.first + g_.len(sg_.alt[q >> 1]);
                if (nd < dist_[q]) {
                    if (dist_[q] == INT64_MAX) touched.push_back(q);
                    dist_[q] = nd;
                    par_[q] = h;
                    pq.push(QE(nd, std::make_pair(alt_key(q), q)));
                }
            }
        }
        half_path.clear();
        if (found) {
            for (uint32_t x = end_h; x != start; x = par_[x]) half_path.push_back(x);
            // backward: x .. towards start is already walk order; forward: reverse into walk order
            if (!backward) std::reverse(half_path.begin(), half_path.end());
        }
        for (uint32_t t : touched) { dist_[t] = INT64_MAX; par_[t] = NONE; }
        return found;
    }

    bool bubble(const std::vector<uint32_t>& R, Cand& c) {
        std::vector<uint8_t> forbid(sg_.n_nodes(), 0);
        for (uint32_t lh : R) forbid[lh >> 1] = 1;
        std::vector<uint32_t> left, right;
        Handle dep = 0, arr = 0;
        if (!half(R.front(), true, forbid, left, dep)) return false;
        for (uint32_t lh : left) forbid[lh >> 1] = 1;
        if (!half(R.back(), false, forbid, right, arr)) return false;
        c.path = left;
        c.path.insert(c.path.end(), R.begin(), R.end());
        c.path.insert(c.path.end(), right.begin(), right.end());
        c.dep = dep;
        c.arr = arr;
        c.total = 0;
        return true;
    }

    void emit(const Cand& c, Excursion& out) const {
        out.dep = c.dep;
        out.arr = c.arr;
        out.alts.clear();
        for (uint32_t lh : c.path) out.alts.push_back(sg_.global(lh));
    }
};

void count_failure(WalkStats& st, int why) {
    switch (why) {
        case CX_CREATOR_LINKS: ++st.fail_creator_links; break;
        case CX_OWN_AMBIGUOUS: ++st.fail_own_ambiguous; break;
        case CX_CREATOR_EDGE: ++st.fail_creator_edge; break;
        case CX_CYCLE: ++st.fail_cycle; break;
        default: ++st.fail_other; break;
    }
}

void witness_runs(const Graph& g, const Runs& R, const SiteData& sd, const std::vector<uint32_t>& todo,
                  const std::vector<uint8_t>& covered, std::vector<Excursion>& out, WalkStats& st) {
    if (todo.empty()) return;
    WitnessFinder wf(g, sd, covered);
    for (uint32_t ri : todo) {
        Excursion e;
        int r = wf.find(R.runs[ri], e);
        if (r == 0) { ++st.unanchored; continue; }
        if (r == 1) ++st.witness_dag;
        else ++st.witness_bubble;
        e.owners.push_back(run_owner(g, R, ri, Src::WITNESS));
        e.src = Src::WITNESS;
        out.push_back(std::move(e));
    }
}

} // namespace

std::vector<Excursion> CreatorWalks::excursions(const SiteData& sd, WalkStats& st) const {
    const Site& s = *sd.site;
    std::vector<Excursion> out;
    std::vector<uint32_t> failed;
    std::vector<uint8_t> covered(sd.sg.n_nodes(), 0);
    for (uint32_t ri : site_runs(runs_, s)) {
        ++st.runs;
        Excursion e;
        int why = creator_excursion(g_, runs_, s, ri, e);
        if (why != CX_OK) {
            ++st.fallback;
            count_failure(st, why);
            failed.push_back(ri);
            continue;
        }
        ++st.creator_ok;
        for (Handle h : e.alts) {
            uint32_t l = sd.sg.local_node(handle_node(h));
            if (l != NONE) covered[l] = 1;
        }
        e.owners.push_back(run_owner(g_, runs_, ri, Src::CREATOR));
        e.src = Src::CREATOR;
        out.push_back(std::move(e));
    }
    witness_runs(g_, runs_, sd, failed, covered, out, st);
    return out;
}

std::vector<Excursion> WitnessWalks::excursions(const SiteData& sd, WalkStats& st) const {
    std::vector<Excursion> out;
    std::vector<uint32_t> todo = site_runs(runs_, *sd.site);
    st.runs += todo.size();
    st.fallback += todo.size();
    std::vector<uint8_t> covered(sd.sg.n_nodes(), 0);
    witness_runs(g_, runs_, sd, todo, covered, out, st);
    return out;
}

// ================================================================ units

Kind classify(const Graph& g, const Excursion& e, int64_t& wlo, int64_t& whi) {
    bool dm = handle_rev(e.dep), am = handle_rev(e.arr);
    if (dm != am) { wlo = whi = 0; return Kind::J; }
    if (!dm) { wlo = g.end(handle_node(e.dep)); whi = g.start(handle_node(e.arr)); }
    else { wlo = g.end(handle_node(e.arr)); whi = g.start(handle_node(e.dep)); }   // (not canonicalized)
    return wlo < whi ? Kind::F : (wlo == whi ? Kind::I : Kind::BK);
}

const AllowedIv& SiteData::allowed_of(NodeId n) const {
    static const AllowedIv none;
    uint32_t l = sg.local_node(n);
    return l == NONE ? none : allowed[l];
}

void prepare_site(const Graph& g, const Site& s, const WalkSource& src, const Options& opt, SiteData& sd) {
    sd.site = &s;
    sd.units.clear();
    sd.wstats = WalkStats();
    build_site_graph(g, s, sd.sg);
    compute_allowed(s, sd.sg, sd.allowed);
    std::vector<Excursion> ex = src.excursions(sd, sd.wstats);
    sd.complete_walks = src.complete();
    sd.walks = &src;
    if (sd.complete_walks) src.refine_allowed(sd, sd.allowed);
    sd.wstats.excursions_raw += ex.size();
    for (Excursion& e : ex) {
        canonicalize(e);
        std::sort(e.owners.begin(), e.owners.end(), owner_less);
    }
    // merge identical canonical excursions, keeping all owners (linear in total length after sorting)
    std::vector<uint32_t> idx(ex.size());
    std::iota(idx.begin(), idx.end(), 0u);
    std::sort(idx.begin(), idx.end(), [&](uint32_t a, uint32_t b) {
        if (excursion_key_less(ex[a], ex[b])) return true;
        if (excursion_key_less(ex[b], ex[a])) return false;
        return a < b;
    });
    std::vector<Excursion> merged;
    for (size_t k = 0; k < idx.size(); ++k) {
        Excursion& e = ex[idx[k]];
        if (!merged.empty() && excursion_key_equal(merged.back(), e)) {
            Excursion& m = merged.back();
            m.owners.insert(m.owners.end(), e.owners.begin(), e.owners.end());
            m.weight += e.weight;
            if ((int)e.src < (int)m.src) m.src = e.src;
            ++sd.wstats.merged;
            continue;
        }
        merged.push_back(std::move(e));
    }
    for (Excursion& m : merged) {
        std::sort(m.owners.begin(), m.owners.end(), owner_less);
        m.owners.erase(std::unique(m.owners.begin(), m.owners.end(), owner_equal), m.owners.end());
    }
    sd.wstats.excursions += merged.size();
    // units
    sd.units.reserve(merged.size());
    const int64_t minw = opt.min_window();
    for (Excursion& m : merged) {
        Unit u;
        u.exc = std::move(m);
        u.kind = classify(g, u.exc, u.wlo, u.whi);
        u.query_bp = 0;
        for (Handle h : u.exc.alts) u.query_bp += g.len(handle_node(h));
        u.feasible_bp = 0;
        if (u.kind == Kind::F) {
            for (Handle h : u.exc.alts) {
                const AllowedIv& a = sd.allowed_of(handle_node(h));
                if (a.status == AllowedStatus::OK && std::max(a.lo, u.wlo) < std::min(a.hi, u.whi)) u.feasible_bp += g.len(handle_node(h));
            }
        }
        sd.units.push_back(std::move(u));
    }
    std::sort(sd.units.begin(), sd.units.end(), [](const Unit& a, const Unit& b) {
        bool ea = a.exc.owners.empty(), eb = b.exc.owners.empty();
        if (ea != eb) return eb;
        if (!ea) {
            const Owner& x = a.exc.owners[0];
            const Owner& y = b.exc.owners[0];
            if (owner_less(x, y)) return true;
            if (owner_less(y, x)) return false;
        }
        return excursion_key_less(a.exc, b.exc);
    });
    std::vector<uint32_t> todo;
    for (uint32_t i = 0; i < (uint32_t)sd.units.size(); ++i) {
        Unit& u = sd.units[i];
        u.id = i;
        if (u.kind == Kind::F && u.query_bp >= opt.b && u.window_bp() >= minw) {
            u.ref_unit = true;
            if (u.query_bp > opt.max_pair || u.window_bp() > opt.max_pair) u.outcome = "pair-too-big";
            else if (u.feasible_bp < opt.b) u.outcome = "infeasible";
            else { u.outcome = ""; todo.push_back(i); }
        } else if (u.query_bp >= opt.b) {
            u.outcome = "no-window";
        } else {
            u.outcome = "small";
        }
    }
    // cap: at most --max-site-query of feasible query per site -- the query bp of the feasible units
    // that are aligned, which is what the aligner spends (the spec's "119 Mb feasible; 329 units
    // (50 Mb) under the cap" at chr1:2.65) -- taking units by feasible bp (descending), then owner
    // SR, SN, SO; the first unit that does not fit and every unit after it are site:capped
    std::sort(todo.begin(), todo.end(), [&](uint32_t a, uint32_t b) {
        const Unit& x = sd.units[a];
        const Unit& y = sd.units[b];
        if (x.feasible_bp != y.feasible_bp) return x.feasible_bp > y.feasible_bp;
        return a < b;   // unit ids already follow owner SR, SN, SO, then the key
    });
    int64_t used = 0;
    bool capped = false;
    for (uint32_t i : todo) {
        Unit& u = sd.units[i];
        if (capped || used + u.query_bp > opt.max_site_query) { capped = true; u.outcome = "site:capped"; continue; }
        used += u.query_bp;
    }
}

} // namespace zip
