/*
  zip_graph.cpp -- rGFA reading with input checks, CSR adjacency, emit, global placement assert.

  Reading keeps the S/L fields and raw tags, plus L-line SR, overlap and L1/L2; the parser reads
  gzip, checks every field, and is independent of line order.  The reference nodes of each SN are
  Graph::ref_nodes.  Placement is checked twice: on input (read_rgfa: reachability and one rank
  per SN, exit 2) and as a final assert (check_placement, exit 3).  Emit sorts nodes and links,
  keys links canonically, keeps overlaps and writes through AtomicWriter.
*/
#include "zip_graph.hpp"

#include <deque>
#include <unordered_set>

namespace zip {

// ---------------------------------------------------------------- small helpers

bool parse_node_name(const std::string& name, int64_t& id) {
    if (name.size() < 2 || name[0] != 's') return false;
    if (name.size() > 2 && name[1] == '0') return false;          // leading zero
    int64_t v = 0;
    for (size_t i = 1; i < name.size(); ++i) {
        char c = name[i];
        if (c < '0' || c > '9') return false;
        if (v > (INT64_MAX - (c - '0')) / 10) return false;
        v = v * 10 + (c - '0');
    }
    id = v;
    return true;
}

std::string make_node_tags(int64_t len, const std::string& sn, int64_t so, int32_t sr) {
    return strf("LN:i:%lld\tSN:Z:%s\tSO:i:%lld\tSR:i:%d", (long long)len, sn.c_str(), (long long)so, (int)sr);
}

std::string make_link_tags(int32_t sr) {
    return strf("SR:i:%d\tL1:i:0\tL2:i:0", (int)sr);
}

// FNV-1a over a byte range, chained
static inline uint64_t fnv1a(uint64_t h, const char* p, size_t n) {
    for (size_t i = 0; i < n; ++i) { h ^= (unsigned char)p[i]; h *= 1099511628211ULL; }
    return h;
}

// ---------------------------------------------------------------- Graph members

NodeId Graph::find_name(const std::string& name) const {
    int64_t id = 0;
    if (!parse_node_name(name, id)) return NONE;
    return find_id(id);
}

int32_t Graph::find_sn(const std::string& sn) const {
    auto it = std::lower_bound(sn_names.begin(), sn_names.end(), sn);
    if (it == sn_names.end() || *it != sn) return -1;
    return (int32_t)(it - sn_names.begin());
}

NodeId Graph::ref_node_at(int32_t sn, int64_t pos) const {
    if (sn < 0 || (size_t)sn >= ref_nodes.size()) return NONE;
    const std::vector<NodeId>& v = ref_nodes[sn];
    size_t lo = 0, hi = v.size();
    while (lo < hi) {                       // first node with so > pos
        size_t mid = (lo + hi) / 2;
        if (nodes[v[mid]].so <= pos) lo = mid + 1; else hi = mid;
    }
    if (lo == 0) return NONE;
    NodeId n = v[lo - 1];
    return pos < end(n) ? n : NONE;
}

std::string Graph::ref_seq(int32_t sn, int64_t lo, int64_t hi) const {
    std::string s;
    if (hi <= lo) return s;
    s.reserve((size_t)(hi - lo));
    int64_t p = lo;
    while (p < hi) {
        NodeId n = ref_node_at(sn, p);
        if (n == NONE) fail(EXIT_INVARIANT, strf("no reference node covers %s:%lld", sn_names[sn].c_str(), (long long)p));
        int64_t a = p - nodes[n].so, b = std::min(hi, end(n)) - nodes[n].so;
        s.append(nodes[n].seq, (size_t)a, (size_t)(b - a));
        p = nodes[n].so + b;
    }
    return s;
}

void Graph::append_handle_seq(std::string& out, Handle h) const {
    const std::string& s = nodes[handle_node(h)].seq;
    if (!handle_rev(h)) { out.append(s); return; }
    size_t base = out.size();
    out.resize(base + s.size());
    for (size_t i = 0; i < s.size(); ++i) out[base + s.size() - 1 - i] = comp_base(s[i]);
}

std::string Graph::spell(const std::vector<Handle>& walk) const {
    size_t n = 0;
    for (Handle h : walk) n += nodes[handle_node(h)].seq.size();
    std::string s;
    s.reserve(n);
    for (Handle h : walk) append_handle_seq(s, h);
    return s;
}

NodeId Graph::add_node(Node&& n) {
    NodeId id = (NodeId)nodes.size();
    id_index_[n.id] = id;
    nodes.push_back(std::move(n));
    return id;
}

uint32_t Graph::add_link(Link&& l) {
    uint32_t i = (uint32_t)links.size();
    links.push_back(std::move(l));
    return i;
}

void Graph::build_index() {
    id_index_.clear();
    id_index_.reserve(nodes.size() * 2);
    for (NodeId i = 0; i < (NodeId)nodes.size(); ++i) id_index_[nodes[i].id] = i;
}

void Graph::build_csr() {
    size_t nh = nodes.size() * 2;
    csr_off_.assign(nh + 1, 0);
    for (const Link& l : links) {
        if (l.deleted) continue;
        csr_off_[handle_leaving(l.a) + 1]++;
        if (l.a != l.b) csr_off_[handle_leaving(l.b) + 1]++;
    }
    for (size_t i = 0; i < nh; ++i) csr_off_[i + 1] += csr_off_[i];
    csr_.assign(csr_off_[nh], Edge{0, 0});
    std::vector<uint32_t> fill(csr_off_.begin(), csr_off_.end() - 1);
    for (uint32_t li = 0; li < (uint32_t)links.size(); ++li) {
        const Link& l = links[li];
        if (l.deleted) continue;
        csr_[fill[handle_leaving(l.a)]++] = Edge{handle_entering(l.b), li};
        if (l.a != l.b) csr_[fill[handle_leaving(l.b)]++] = Edge{handle_entering(l.a), li};
    }
    for (size_t h = 0; h < nh; ++h)
        std::sort(csr_.begin() + csr_off_[h], csr_.begin() + csr_off_[h + 1],
                  [](const Edge& x, const Edge& y) { return x.to != y.to ? x.to < y.to : x.link < y.link; });
}

// ---------------------------------------------------------------- reading

namespace {

struct Problems {
    std::vector<std::string> msgs;
    size_t count = 0;
    void add(const std::string& m) {
        ++count;
        if (msgs.size() < 12) msgs.push_back(m);
    }
    void raise(const std::string& path) const {
        if (!count) return;
        std::string s = strf("%s is not a valid rGFA for rgfa-zip (%zu problem(s)):", path.c_str(), count);
        for (const std::string& m : msgs) s += "\n  " + m;
        if (count > msgs.size()) s += strf("\n  ... and %zu more", count - msgs.size());
        fail(EXIT_INPUT, s);
    }
};

struct RawLink {
    int64_t from = 0, to = 0;
    bool from_fwd = true, to_fwd = true;
    int32_t sr = -1;
    std::string overlap, tags;
    int64_t l1 = -1, l2 = -1;
    uint64_t lineno = 0;
};

// value of the first tag `key` (e.g. "SR") of type `type` in a tab-joined tag string
bool find_int_tag(const std::string& tags, const char* key, int64_t& out) {
    size_t p = 0;
    while (p < tags.size()) {
        size_t q = tags.find('\t', p);
        if (q == std::string::npos) q = tags.size();
        if (q - p >= 5 && tags[p] == key[0] && tags[p + 1] == key[1] && tags[p + 2] == ':' && tags[p + 3] == 'i' &&
            tags[p + 4] == ':') {
            int64_t v = 0;
            if (parse_int64(tags.substr(p + 5, q - p - 5), v)) { out = v; return true; }
            return false;
        }
        p = q + 1;
    }
    return false;
}

} // namespace

void read_rgfa(const std::string& path, Graph& g) {
    LineReader in(path);
    Problems prob;
    std::string line;
    std::vector<std::string> f;
    std::vector<std::string> node_sn;
    std::vector<RawLink> raw;
    uint64_t n_comment = 0;
    while (in.next(line)) {
        if (line.empty()) continue;
        char t = line[0];
        uint64_t ln = in.line_number();
        if (t == 'H') { g.header_lines.push_back(line); continue; }
        if (t == '#') { ++n_comment; continue; }
        if (t == 'S') {
            if (line.size() < 2 || line[1] != '\t') { prob.add(strf("line %llu: malformed S line", (unsigned long long)ln)); continue; }
            size_t t2 = line.find('\t', 2);
            if (t2 == std::string::npos) { prob.add(strf("line %llu: S line without a sequence", (unsigned long long)ln)); continue; }
            std::string name = line.substr(2, t2 - 2);
            size_t t3 = line.find('\t', t2 + 1);
            Node n;
            n.seq = line.substr(t2 + 1, t3 == std::string::npos ? std::string::npos : t3 - t2 - 1);
            if (t3 != std::string::npos) n.tags = line.substr(t3 + 1);
            if (!parse_node_name(name, n.id)) {
                prob.add(strf("line %llu: segment name '%s' is not s<int>", (unsigned long long)ln, name.c_str()));
                continue;
            }
            if (n.seq.empty() || n.seq == "*") {
                prob.add(strf("line %llu: %s has no sequence", (unsigned long long)ln, name.c_str()));
                continue;
            }
            std::string sn;
            bool have_sn = false, have_so = false, have_sr = false, bad = false;
            int64_t lnv = -1;
            size_t p = 0;
            const std::string& tg = n.tags;
            while (p < tg.size()) {
                size_t q = tg.find('\t', p);
                if (q == std::string::npos) q = tg.size();
                std::string tag = tg.substr(p, q - p);
                p = q + 1;
                if (tag.size() < 5 || tag[2] != ':' || tag[4] != ':') {
                    prob.add(strf("line %llu: %s: malformed tag '%s'", (unsigned long long)ln, name.c_str(), tag.c_str()));
                    bad = true;
                    continue;
                }
                std::string key = tag.substr(0, 2), val = tag.substr(5);
                char type = tag[3];
                int64_t v = 0;
                if (key == "SN") {
                    if (have_sn || type != 'Z' || val.empty()) { prob.add(strf("line %llu: %s: bad or repeated SN tag", (unsigned long long)ln, name.c_str())); bad = true; }
                    have_sn = true; sn = val;
                } else if (key == "SO") {
                    if (have_so || type != 'i' || !parse_int64(val, v) || v < 0) { prob.add(strf("line %llu: %s: bad or repeated SO tag", (unsigned long long)ln, name.c_str())); bad = true; }
                    have_so = true; n.so = v;
                } else if (key == "SR") {
                    if (have_sr || type != 'i' || !parse_int64(val, v) || v < 0 || v > INT32_MAX) { prob.add(strf("line %llu: %s: bad or repeated SR tag", (unsigned long long)ln, name.c_str())); bad = true; }
                    have_sr = true; n.sr = (int32_t)v;
                } else if (key == "LN") {
                    if (type != 'i' || !parse_int64(val, v)) { prob.add(strf("line %llu: %s: bad LN tag", (unsigned long long)ln, name.c_str())); bad = true; }
                    lnv = v;
                }
            }
            if (!have_sn || !have_so || !have_sr) {
                prob.add(strf("line %llu: %s lacks %s%s%s (every S line needs SN, SO and SR)", (unsigned long long)ln, name.c_str(),
                              have_sn ? "" : "SN ", have_so ? "" : "SO ", have_sr ? "" : "SR"));
                continue;
            }
            if (lnv >= 0 && lnv != n.len()) {
                prob.add(strf("line %llu: %s: LN:i:%lld but the sequence has %lld bp", (unsigned long long)ln, name.c_str(), (long long)lnv, (long long)n.len()));
                continue;
            }
            if (bad) continue;
            g.nodes.push_back(std::move(n));
            node_sn.push_back(sn);
        } else if (t == 'L') {
            split_tabs(line, f);
            if (f.size() < 6) { prob.add(strf("line %llu: L line with %zu fields", (unsigned long long)ln, f.size())); continue; }
            RawLink r;
            r.lineno = ln;
            if (!parse_node_name(f[1], r.from) || !parse_node_name(f[3], r.to)) {
                prob.add(strf("line %llu: link end '%s' or '%s' is not s<int>", (unsigned long long)ln, f[1].c_str(), f[3].c_str()));
                continue;
            }
            if ((f[2] != "+" && f[2] != "-") || (f[4] != "+" && f[4] != "-")) {
                prob.add(strf("line %llu: link orientation must be + or -", (unsigned long long)ln));
                continue;
            }
            r.from_fwd = f[2] == "+";
            r.to_fwd = f[4] == "+";
            r.overlap = f[5];
            if (r.overlap != "0M" && r.overlap != "*") {
                prob.add(strf("line %llu: overlap '%s' (rgfa-zip needs blunt links: 0M or *)", (unsigned long long)ln, r.overlap.c_str()));
                continue;
            }
            for (size_t i = 6; i < f.size(); ++i) {
                if (i > 6) r.tags.push_back('\t');
                r.tags += f[i];
            }
            int64_t v = 0;
            if (find_int_tag(r.tags, "SR", v)) {
                if (v < 0 || v > INT32_MAX) { prob.add(strf("line %llu: bad SR tag", (unsigned long long)ln)); continue; }
                r.sr = (int32_t)v;
            }
            if (find_int_tag(r.tags, "L1", v)) r.l1 = v;
            if (find_int_tag(r.tags, "L2", v)) r.l2 = v;
            raw.push_back(std::move(r));
        } else {
            prob.add(strf("line %llu: unsupported GFA line type '%c' (rgfa-zip reads minigraph rGFA: H, S and L lines)",
                          (unsigned long long)ln, t));
        }
    }
    in.close();
    prob.raise(path);

    // ---- nodes: sort by id, unique ids
    {
        std::vector<size_t> perm(g.nodes.size());
        for (size_t i = 0; i < perm.size(); ++i) perm[i] = i;
        std::sort(perm.begin(), perm.end(), [&](size_t a, size_t b) { return g.nodes[a].id < g.nodes[b].id; });
        std::vector<Node> sorted;
        std::vector<std::string> sns;
        sorted.reserve(perm.size());
        sns.reserve(perm.size());
        for (size_t i : perm) { sorted.push_back(std::move(g.nodes[i])); sns.push_back(std::move(node_sn[i])); }
        g.nodes.swap(sorted);
        node_sn.swap(sns);
    }
    if (g.nodes.empty()) fail(EXIT_INPUT, strf("%s has no S lines", path.c_str()));
    if (g.nodes.size() >= (size_t)(UINT32_MAX / 2 - 16)) fail(EXIT_INPUT, "too many nodes");
    for (size_t i = 1; i < g.nodes.size(); ++i)
        if (g.nodes[i].id == g.nodes[i - 1].id) prob.add(strf("segment s%lld is defined twice", (long long)g.nodes[i].id));
    prob.raise(path);
    g.max_id = g.nodes.back().id;
    g.n_input_nodes = g.nodes.size();
    g.build_index();

    // ---- SN table, one rank per SN
    {
        std::vector<std::string> names = node_sn;
        std::sort(names.begin(), names.end());
        names.erase(std::unique(names.begin(), names.end()), names.end());
        g.sn_names = names;
        g.sn_rank.assign(names.size(), -1);
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            int32_t s = g.find_sn(node_sn[i]);
            g.nodes[i].sn = s;
            int32_t& rk = g.sn_rank[s];
            if (rk < 0) rk = g.nodes[i].sr;
            else if (rk != g.nodes[i].sr)
                prob.add(strf("SN %s has two ranks (%d and %d, e.g. at s%lld)", node_sn[i].c_str(), rk, g.nodes[i].sr,
                              (long long)g.nodes[i].id));
        }
        prob.raise(path);
    }

    // ---- rank-0 nodes tile their SN by SO
    g.ref_nodes.assign(g.sn_names.size(), std::vector<NodeId>());
    g.ref_hash.assign(g.sn_names.size(), 0);
    {
        bool any_ref = false;
        for (NodeId i = 0; i < (NodeId)g.nodes.size(); ++i)
            if (g.nodes[i].sr == 0) { g.ref_nodes[g.nodes[i].sn].push_back(i); any_ref = true; }
        if (!any_ref) fail(EXIT_INPUT, strf("%s has no rank-0 (SR:i:0) segment", path.c_str()));
        for (size_t s = 0; s < g.ref_nodes.size(); ++s) {
            std::vector<NodeId>& v = g.ref_nodes[s];
            if (v.empty()) continue;
            std::sort(v.begin(), v.end(), [&](NodeId a, NodeId b) {
                return g.nodes[a].so != g.nodes[b].so ? g.nodes[a].so < g.nodes[b].so : a < b;
            });
            if (g.nodes[v[0]].so != 0)
                ZLOG("warning: reference %s starts at SO %lld, not 0", g.sn_names[s].c_str(), (long long)g.nodes[v[0]].so);
            uint64_t hsh = 1469598103934665603ULL;
            for (NodeId x : v) hsh = fnv1a(hsh, g.nodes[x].seq.data(), g.nodes[x].seq.size());
            g.ref_hash[s] = hsh;
            for (size_t k = 1; k < v.size(); ++k) {
                if (g.nodes[v[k]].so != g.end(v[k - 1]))
                    prob.add(strf("rank-0 segments of %s do not tile it: %s ends at %lld, %s starts at %lld", g.sn_names[s].c_str(),
                                  g.name(v[k - 1]).c_str(), (long long)g.end(v[k - 1]), g.name(v[k]).c_str(),
                                  (long long)g.nodes[v[k]].so));
            }
        }
        prob.raise(path);
    }

    // ---- links: resolve, canonical order, duplicates
    {
        uint64_t no_sr = 0, bad_len = 0;
        g.links.clear();
        g.links.reserve(raw.size());
        for (RawLink& r : raw) {
            NodeId a = g.find_id(r.from), b = g.find_id(r.to);
            if (a == NONE || b == NONE) {
                prob.add(strf("line %llu: link end s%lld does not exist", (unsigned long long)r.lineno, (long long)(a == NONE ? r.from : r.to)));
                continue;
            }
            Link l;
            l.a = r.from_fwd ? right_side(a) : left_side(a);
            l.b = r.to_fwd ? left_side(b) : right_side(b);
            l.sr = r.sr;
            l.overlap = r.overlap;
            l.tags = std::move(r.tags);
            if ((r.l1 >= 0 && r.l1 != g.len(a)) || (r.l2 >= 0 && r.l2 != g.len(b))) ++bad_len;
            g.links.push_back(std::move(l));
        }
        raw.clear();
        raw.shrink_to_fit();
        prob.raise(path);
        // Duplicates (one link written twice, or in both readings) keep one copy: a copy with SR
        // first (a creator link must not lose its rank to an SR-less duplicate), the lowest SR,
        // then the reading, overlap and tags -- a total order on content, so the choice does not
        // depend on line order
        std::sort(g.links.begin(), g.links.end(), [](const Link& x, const Link& y) {
            LinkKey kx = link_key(x.a, x.b), ky = link_key(y.a, y.b);
            if (kx != ky) return kx < ky;
            if ((x.sr < 0) != (y.sr < 0)) return y.sr < 0;
            if (x.sr != y.sr) return x.sr < y.sr;
            if (x.a != y.a) return x.a < y.a;
            if (x.overlap != y.overlap) return x.overlap < y.overlap;
            return x.tags < y.tags;
        });
        size_t w = 0, dups = 0, sr_conflicts = 0;
        for (size_t i = 0; i < g.links.size(); ++i) {
            if (w > 0 && link_key(g.links[w - 1].a, g.links[w - 1].b) == link_key(g.links[i].a, g.links[i].b)) {
                ++dups;
                if (g.links[i].sr >= 0 && g.links[i].sr != g.links[w - 1].sr) ++sr_conflicts;
                continue;
            }
            if (w != i) g.links[w] = std::move(g.links[i]);
            ++w;
        }
        g.links.resize(w);
        g.n_input_links = g.links.size();
        for (const Link& l : g.links)
            if (l.sr < 0) ++no_sr;
        if (dups) ZLOG("warning: %zu duplicate link(s) in %s (a link equals its reverse); kept one of each, with SR where a copy has it",
                       dups, path.c_str());
        if (sr_conflicts)
            ZLOG("warning: %zu duplicate link(s) in %s carry different SR values; the lowest SR is kept", sr_conflicts, path.c_str());
        if (no_sr)
            ZLOG("warning: %llu L line(s) without SR in %s: their runs use the witness fallback", (unsigned long long)no_sr, path.c_str());
        if (bad_len) ZLOG("warning: %llu L line(s) with L1/L2 not matching the segment lengths; rewritten on output", (unsigned long long)bad_len);
    }
    g.build_csr();

    // ---- placeable: every node reachable from rank 0, ignoring orientation (as rgfa-split does)
    {
        std::vector<uint8_t> seen(g.nodes.size(), 0);
        std::vector<NodeId> st;
        for (NodeId i = 0; i < (NodeId)g.nodes.size(); ++i)
            if (g.nodes[i].sr == 0) { seen[i] = 1; st.push_back(i); }
        while (!st.empty()) {
            NodeId x = st.back();
            st.pop_back();
            for (int o = 0; o < 2; ++o)
                for (const Edge& e : g.out(make_handle(x, o == 1))) {
                    NodeId y = handle_node(e.to);
                    if (!seen[y]) { seen[y] = 1; st.push_back(y); }
                }
        }
        size_t bad = 0;
        NodeId first = NONE;
        for (NodeId i = 0; i < (NodeId)g.nodes.size(); ++i)
            if (!seen[i]) { ++bad; if (first == NONE) first = i; }
        if (bad)
            fail(EXIT_INPUT, strf("%s is not placeable: %zu segment(s) are not reachable from any rank-0 segment (first: %s)",
                                  path.c_str(), bad, g.name(first).c_str()));
    }
    ZLOG("read %s: %zu segments, %zu links, %zu SN (%zu with rank 0)%s", path.c_str(), g.nodes.size(), g.links.size(),
         g.sn_names.size(), (size_t)std::count_if(g.ref_nodes.begin(), g.ref_nodes.end(), [](const std::vector<NodeId>& v) { return !v.empty(); }),
         n_comment ? strf(", %llu comment line(s) dropped", (unsigned long long)n_comment).c_str() : "");
}

// ---------------------------------------------------------------- emit

namespace {

// rewrite L1:i:/L2:i: values from the current endpoint lengths
std::string emit_link_tags(const std::string& tags, int64_t l1, int64_t l2) {
    if (tags.find("L1:i:") == std::string::npos && tags.find("L2:i:") == std::string::npos) return tags;
    std::string out;
    out.reserve(tags.size() + 8);
    size_t p = 0;
    while (p <= tags.size()) {
        size_t q = tags.find('\t', p);
        if (q == std::string::npos) q = tags.size();
        if (p > 0) out.push_back('\t');
        if (q - p >= 5 && tags.compare(p, 5, "L1:i:") == 0) out += "L1:i:" + std::to_string(l1);
        else if (q - p >= 5 && tags.compare(p, 5, "L2:i:") == 0) out += "L2:i:" + std::to_string(l2);
        else out.append(tags, p, q - p);
        if (q == tags.size()) break;
        p = q + 1;
    }
    return out;
}

struct IdSide {
    int64_t id;
    int r;
    bool operator<(const IdSide& o) const { return id != o.id ? id < o.id : r < o.r; }
    bool operator==(const IdSide& o) const { return id == o.id && r == o.r; }
};

} // namespace

void write_rgfa(AtomicWriter& w, const Graph& g) {
    for (const std::string& h : g.header_lines) { w.write(h); w.put('\n'); }
    std::vector<NodeId> order;
    order.reserve(g.nodes.size());
    for (NodeId i = 0; i < (NodeId)g.nodes.size(); ++i)
        if (!g.nodes[i].deleted) order.push_back(i);
    std::sort(order.begin(), order.end(), [&](NodeId a, NodeId b) { return g.nodes[a].id < g.nodes[b].id; });
    std::string buf;
    for (NodeId n : order) {
        const Node& nd = g.nodes[n];
        buf.clear();
        buf += "S\ts";
        buf += std::to_string(nd.id);
        buf.push_back('\t');
        w.write(buf);
        w.write(nd.seq);
        if (!nd.tags.empty()) { w.put('\t'); w.write(nd.tags); }
        w.put('\n');
    }
    struct LK { IdSide lo, hi; uint32_t i; };
    std::vector<LK> lk;
    lk.reserve(g.links.size());
    for (uint32_t i = 0; i < (uint32_t)g.links.size(); ++i) {
        const Link& l = g.links[i];
        if (l.deleted) continue;
        IdSide a{g.nodes[side_node(l.a)].id, side_is_right(l.a) ? 1 : 0}, b{g.nodes[side_node(l.b)].id, side_is_right(l.b) ? 1 : 0};
        if (b < a) std::swap(a, b);
        lk.push_back(LK{a, b, i});
    }
    std::sort(lk.begin(), lk.end(), [](const LK& x, const LK& y) {
        if (!(x.lo == y.lo)) return x.lo < y.lo;
        if (!(x.hi == y.hi)) return x.hi < y.hi;
        return x.i < y.i;
    });
    for (const LK& k : lk) {
        const Link& l = g.links[k.i];
        NodeId a = side_node(l.a), b = side_node(l.b);
        buf.clear();
        buf += "L\ts";
        buf += std::to_string(g.nodes[a].id);
        buf += side_is_right(l.a) ? "\t+\ts" : "\t-\ts";
        buf += std::to_string(g.nodes[b].id);
        buf += side_is_right(l.b) ? "\t-\t" : "\t+\t";
        buf += l.overlap;
        if (!l.tags.empty()) {
            buf.push_back('\t');
            buf += emit_link_tags(l.tags, g.len(a), g.len(b));
        }
        buf.push_back('\n');
        w.write(buf);
    }
}

// ---------------------------------------------------------------- final assert

void check_placement(const Graph& g) {
    std::vector<std::string> errs;
    auto add = [&](const std::string& m) { if (errs.size() < 12) errs.push_back(m); };
    size_t n_err = 0;
    // unique ids
    {
        std::vector<int64_t> ids;
        ids.reserve(g.nodes.size());
        for (const Node& n : g.nodes) if (!n.deleted) ids.push_back(n.id);
        std::sort(ids.begin(), ids.end());
        for (size_t i = 1; i < ids.size(); ++i)
            if (ids[i] == ids[i - 1]) { ++n_err; add(strf("id s%lld is used twice", (long long)ids[i])); }
    }
    // dangling and duplicate links
    {
        std::unordered_set<LinkKey, LinkKeyHash> seen;
        seen.reserve(g.links.size() * 2);
        for (const Link& l : g.links) {
            if (l.deleted) continue;
            NodeId a = side_node(l.a), b = side_node(l.b);
            if (a >= g.nodes.size() || b >= g.nodes.size() || g.nodes[a].deleted || g.nodes[b].deleted) {
                ++n_err;
                add(strf("dangling link %s-%s", a < g.nodes.size() ? g.name(a).c_str() : "?", b < g.nodes.size() ? g.name(b).c_str() : "?"));
                continue;
            }
            if (!seen.insert(link_key(l.a, l.b)).second) { ++n_err; add(strf("duplicate link %s-%s", g.name(a).c_str(), g.name(b).c_str())); }
        }
    }
    // one rank per SN
    std::vector<int32_t> rank(g.sn_names.size(), -1);
    for (const Node& n : g.nodes) {
        if (n.deleted) continue;
        if (n.sn < 0 || (size_t)n.sn >= g.sn_names.size()) { ++n_err; add(strf("s%lld has no SN", (long long)n.id)); continue; }
        if (rank[n.sn] < 0) rank[n.sn] = n.sr;
        else if (rank[n.sn] != n.sr) { ++n_err; add(strf("SN %s carries two ranks", g.sn_names[n.sn].c_str())); }
    }
    // rank-0 tiling unchanged: same span, contiguous
    {
        std::vector<std::vector<NodeId>> ref(g.sn_names.size());
        for (NodeId i = 0; i < (NodeId)g.nodes.size(); ++i)
            if (!g.nodes[i].deleted && g.nodes[i].sr == 0 && g.nodes[i].sn >= 0) ref[g.nodes[i].sn].push_back(i);
        for (size_t s = 0; s < ref.size(); ++s) {
            std::vector<NodeId>& v = ref[s];
            const std::vector<NodeId>& in = g.ref_nodes[s];
            if (v.empty() && in.empty()) continue;
            if (v.empty() || in.empty()) { ++n_err; add(strf("reference %s gained or lost all its rank-0 segments", g.sn_names[s].c_str())); continue; }
            std::sort(v.begin(), v.end(), [&](NodeId a, NodeId b) { return g.nodes[a].so < g.nodes[b].so; });
            int64_t in_lo = g.nodes[in.front()].so, in_hi = g.nodes[in.back()].so + g.nodes[in.back()].len();
            if (g.nodes[v.front()].so != in_lo) { ++n_err; add(strf("reference %s now starts at %lld", g.sn_names[s].c_str(), (long long)g.nodes[v.front()].so)); }
            int64_t p = g.nodes[v.front()].so;
            for (NodeId x : v) {
                if (g.nodes[x].so != p) { ++n_err; add(strf("rank-0 segments of %s no longer tile it at %lld", g.sn_names[s].c_str(), (long long)p)); break; }
                p += g.nodes[x].len();
            }
            if (p != in_hi) { ++n_err; add(strf("reference %s now ends at %lld, not %lld", g.sn_names[s].c_str(), (long long)p, (long long)in_hi)); }
            uint64_t hsh = 1469598103934665603ULL;
            for (NodeId x : v) hsh = fnv1a(hsh, g.nodes[x].seq.data(), g.nodes[x].seq.size());
            if (hsh != g.ref_hash[s]) { ++n_err; add(strf("reference sequence of %s changed", g.sn_names[s].c_str())); }
        }
    }
    // every node reachable from rank 0
    {
        size_t N = g.nodes.size();
        std::vector<uint32_t> off(N + 1, 0);
        for (const Link& l : g.links) {
            if (l.deleted) continue;
            NodeId a = side_node(l.a), b = side_node(l.b);
            if (a >= N || b >= N) continue;
            off[a + 1]++;
            off[b + 1]++;
        }
        for (size_t i = 0; i < N; ++i) off[i + 1] += off[i];
        std::vector<NodeId> adj(off[N]);
        std::vector<uint32_t> fill(off.begin(), off.end() - 1);
        for (const Link& l : g.links) {
            if (l.deleted) continue;
            NodeId a = side_node(l.a), b = side_node(l.b);
            if (a >= N || b >= N) continue;
            adj[fill[a]++] = b;
            adj[fill[b]++] = a;
        }
        std::vector<uint8_t> seen(N, 0);
        std::vector<NodeId> st;
        for (NodeId i = 0; i < (NodeId)N; ++i)
            if (!g.nodes[i].deleted && g.nodes[i].sr == 0) { seen[i] = 1; st.push_back(i); }
        while (!st.empty()) {
            NodeId x = st.back();
            st.pop_back();
            for (uint32_t k = off[x]; k < off[x + 1]; ++k) {
                NodeId y = adj[k];
                if (!seen[y] && !g.nodes[y].deleted) { seen[y] = 1; st.push_back(y); }
            }
        }
        size_t bad = 0;
        for (NodeId i = 0; i < (NodeId)N; ++i)
            if (!g.nodes[i].deleted && !seen[i]) { ++bad; if (bad <= 3) add(strf("%s is not reachable from rank 0", g.name(i).c_str())); }
        n_err += bad;
    }
    if (n_err) {
        std::string s = strf("global placement assert failed (%zu problem(s)); refusing to emit:", n_err);
        for (const std::string& m : errs) s += "\n  " + m;
        fail(EXIT_INVARIANT, s);
    }
}

} // namespace zip
