/*
  zip_gaf.cpp -- GafWalks (v3): observed haplotype walks from a graphmap GAF, and the GAF audit.

  Layout
    names       GAF/graph name canonicalization (the design prototype's gaf_name_to_sn) and SN resolution
    GafReader   line input through `bgzip -@ t -dc` (BGZF) or `gzip -dc`, keeping only the head of
                each line (the 12 mandatory fields and the first 512 bytes of tags), so the multi-GB
                cg/ds tags are never copied
    GafIndex    projection of every record onto handles (the design prototype's path_to_handles, plus a
                link check), excursions (as the design prototype cut them, in canonical frame), owners and weights, the
                per-node excursion index, partial walk ends, and per-record reference coverage for
                the audit (as in the design prototype)
    GafWalks    the WalkSource: observed excursions of a site, creator owners, observed Allowed
                (as the prototype's iv_from_excs), and the creator-vs-observed diagnostic (the prototype's
                "creator reconstruction exactness")
    GafAudit    the prototype audit's window and revisit checks, per zipped piece
*/
#include "zip_gaf.hpp"
#include "zip_edit.hpp"

#include <chrono>
#include <numeric>
#include <unordered_map>

namespace zip {

// ================================================================ names

namespace {

// true while unwinding from an error (outputs and summaries are then skipped)
inline bool unwinding() {
#if __cplusplus >= 201703L
    return std::uncaught_exceptions() > 0;
#else
    return std::uncaught_exception();
#endif
}

const uint8_t KIND_X = 255;   // anchors on two different SNs (possible only with several rank-0 SNs)

// "id=HG00320.1|JBHIKI010000013.1" -> "HG00320#1#JBHIKI010000013.1"; "id=GRCh38|chr16" ->
// "GRCh38#0#chr16" (the design prototype's gaf_name_to_sn); "id=X" -> "X"; anything else is unchanged.
std::string canon_name(const std::string& x) {
    if (x.compare(0, 3, "id=") != 0) return x;
    std::string y = x.substr(3);
    size_t bar = y.find('|');
    if (bar == std::string::npos) return y;
    std::string samp = y.substr(0, bar), ctg = y.substr(bar + 1);
    size_t dot = samp.rfind('.');
    if (dot != std::string::npos) return samp.substr(0, dot) + "#" + samp.substr(dot + 1) + "#" + ctg;
    return samp + "#0#" + ctg;
}

// the haplotype of a canonical contig name: its first two '#' fields ("HG00320#1"), else the name
std::string hap_name(const std::string& c) {
    size_t p = c.find('#');
    if (p == std::string::npos) return c;
    size_t q = c.find('#', p + 1);
    if (q == std::string::npos) return c;
    return c.substr(0, q);
}

std::string last_field(const std::string& c) {
    size_t p = c.rfind('#');
    return p == std::string::npos ? c : c.substr(p + 1);
}

// the executable `name` on PATH, or ""
std::string find_on_path(const char* name) {
    const char* p = getenv("PATH");
    if (!p) return std::string();
    std::string path(p);
    size_t i = 0;
    while (i <= path.size()) {
        size_t j = path.find(':', i);
        if (j == std::string::npos) j = path.size();
        std::string dir = path.substr(i, j - i);
        if (dir.empty()) dir = ".";
        std::string f = dir + "/" + name;
        struct stat st;
        if (stat(f.c_str(), &st) == 0 && S_ISREG(st.st_mode) && access(f.c_str(), X_OK) == 0) return f;
        i = j + 1;
    }
    return std::string();
}

// a bgzip that decompresses BGZF with -@ threads: htslib's, version 1.4 or later
bool usable_bgzip(const std::string& path) {
    ProcResult r = run_process({path, "--version"}, std::string(), 1024);
    if (!r.ok()) return false;
    size_t p = r.out.find("(htslib) ");
    if (p == std::string::npos) return false;
    int major = 0, minor = 0;
    if (sscanf(r.out.c_str() + p + 9, "%d.%d", &major, &minor) != 2) return false;
    return major > 1 || (major == 1 && minor >= 4);
}

// ================================================================ GafReader

class GafReader {
public:
    GafReader(const std::string& path, int threads) : path_(path), buf_(4u << 20) {
        int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) fail(EXIT_INPUT, strf("cannot open GAF %s: %s", path.c_str(), strerror(errno)));
        unsigned char h[18];
        memset(h, 0, sizeof(h));
        ssize_t n = pread(fd, h, sizeof(h), 0);
        if (!(n >= 2 && h[0] == 0x1f && h[1] == 0x8b)) {
            fd_ = fd;
            method_ = "plain";
            return;
        }
        ::close(fd);
        // BGZF: gzip with FEXTRA whose first subfield is 'BC'; bgzip decompresses it in parallel
        bool bgzf = n >= 18 && (h[3] & 4) != 0 && h[12] == 'B' && h[13] == 'C';
        std::string bg = bgzf ? find_on_path("bgzip") : std::string();
        if (!bg.empty() && !usable_bgzip(bg)) bg.clear();
        std::vector<std::string> argv;
        if (!bg.empty()) {
            int t = std::max(1, std::min(threads, 16));
            argv = {bg, "-@", std::to_string(t), "-dc", "--", path};
            method_ = strf("bgzip -@ %d", t);
        } else {
            argv = {"gzip", "-dc", "--", path};
            method_ = "gzip";
        }
        int p[2];
        if (pipe2(p, O_CLOEXEC) != 0) fail(EXIT_INPUT, strf("pipe: %s", strerror(errno)));
        posix_spawn_file_actions_t fa;
        posix_spawn_file_actions_init(&fa);
        posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
        posix_spawn_file_actions_adddup2(&fa, p[1], 1);
        int rc = spawn_child(pid_, argv, &fa);
        posix_spawn_file_actions_destroy(&fa);
        ::close(p[1]);
        if (rc != 0) {
            ::close(p[0]);
            pid_ = -1;
            fail(EXIT_INPUT, strf("cannot run %s to read %s: %s", argv[0].c_str(), path.c_str(), strerror(rc)));
        }
        fd_ = p[0];
    }
    ~GafReader() {
        if (fd_ >= 0) ::close(fd_);
        if (pid_ > 0) {
            kill(pid_, SIGTERM);
            int st;
            reap_child(pid_, st, nullptr);
        }
    }
    GafReader(const GafReader&) = delete;
    GafReader& operator=(const GafReader&) = delete;

    // The next line's head: fields 1-12 in full, then at most 512 bytes of tags (the tail of a
    // longer line is skipped without copying).  False at end of input.
    bool next(std::string& head) {
        head.clear();
        int tabs = 0;
        size_t tag_bytes = 0;
        bool got = false;
        while (true) {
            if (pos_ >= len_ && !fill()) {
                if (!got) return false;
                ++lineno_;
                return true;
            }
            got = true;
            const char* s = buf_.data() + pos_;
            size_t n = len_ - pos_;
            const char* nl = (const char*)memchr(s, '\n', n);
            size_t seg = nl ? (size_t)(nl - s) : n;
            size_t i = 0;
            while (i < seg && tabs < 12) {
                const char* t = (const char*)memchr(s + i, '\t', seg - i);
                size_t upto = t ? (size_t)(t - s) + 1 : seg;
                head.append(s + i, upto - i);
                if (t) ++tabs;
                i = upto;
            }
            if (tabs >= 12 && i < seg && tag_bytes < 512) {
                size_t take = std::min(seg - i, (size_t)512 - tag_bytes);
                head.append(s + i, take);
                tag_bytes += take;
            }
            pos_ += seg + (nl ? 1 : 0);
            if (nl) {
                ++lineno_;
                if (!head.empty() && head.back() == '\r') head.pop_back();
                return true;
            }
        }
    }
    // close and check the decompressor's exit status
    void close() {
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
        if (pid_ > 0) {
            int st = 0;
            bool waited = reap_child(pid_, st, nullptr);
            pid_ = -1;
            if (!waited) fail(EXIT_INPUT, strf("cannot read the exit status of the decompressor of GAF %s", path_.c_str()));
            if (!(WIFEXITED(st) && WEXITSTATUS(st) == 0))
                fail(EXIT_INPUT, strf("decompressing GAF %s failed (corrupt or truncated?)", path_.c_str()));
        }
    }
    uint64_t line_number() const { return lineno_; }
    const std::string& method() const { return method_; }

private:
    bool fill() {
        if (eof_) return false;
        check_abort();   // once per buffer: a stop signal does not wait for a multi-GB GAF
        while (true) {
            ssize_t n = read(fd_, buf_.data(), buf_.size());
            if (n < 0) {
                if (errno == EINTR) continue;
                fail(EXIT_INPUT, strf("read error on GAF %s: %s", path_.c_str(), strerror(errno)));
            }
            pos_ = 0;
            len_ = (size_t)n;
            if (n == 0) { eof_ = true; return false; }
            return true;
        }
    }
    std::string path_, method_;
    int fd_ = -1;
    pid_t pid_ = -1;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
    bool eof_ = false;
    uint64_t lineno_ = 0;
};

inline bool parse_i64(const char* p, size_t n, int64_t& out) {
    if (n == 0 || n > 18) return false;
    int64_t v = 0;
    for (size_t i = 0; i < n; ++i) {
        if (p[i] < '0' || p[i] > '9') return false;
        v = v * 10 + (p[i] - '0');
    }
    out = v;
    return true;
}

// "name:lo-hi" (the last ':' starts the interval); false if the step is not of that form
bool parse_stable(const char* t, size_t n, size_t& name_len, int64_t& lo, int64_t& hi) {
    size_t c = n;
    while (c > 0 && t[c - 1] != ':') --c;
    if (c <= 1) return false;          // no ':' or an empty name
    size_t d = c;
    while (d < n && t[d] != '-') ++d;
    if (d >= n) return false;
    if (!parse_i64(t + c, d - c, lo) || !parse_i64(t + d + 1, n - d - 1, hi)) return false;
    name_len = c - 1;
    return true;
}

struct RefIv {
    int32_t sn;
    int64_t lo, hi;
    bool operator<(const RefIv& o) const { return sn != o.sn ? sn < o.sn : (lo != o.lo ? lo < o.lo : hi < o.hi); }
};

// sort and merge intervals that overlap or touch (the design prototype's merge_iv)
void merge_ivs(std::vector<RefIv>& v) {
    std::sort(v.begin(), v.end());
    size_t k = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        if (k > 0 && v[k - 1].sn == v[i].sn && v[i].lo <= v[k - 1].hi) v[k - 1].hi = std::max(v[k - 1].hi, v[i].hi);
        else v[k++] = v[i];
    }
    v.resize(k);
}

// longest overlap of [ta, tb) on `sn` with the sorted, merged intervals [b, e).  Merged intervals
// of one SN are disjoint, so their ends are sorted as well: the search starts at the first one
// ending after ta (those ending earlier overlap by <= 0), and the scan stops at the first one
// starting at or after tb.
int64_t max_overlap(const RefIv* b, const RefIv* e, int32_t sn, int64_t ta, int64_t tb) {
    const RefIv* it = std::partition_point(b, e, [&](const RefIv& x) { return x.sn < sn || (x.sn == sn && x.hi <= ta); });
    int64_t best = 0;
    for (; it != e && it->sn == sn && it->lo < tb; ++it) best = std::max(best, std::min(it->hi, tb) - std::max(it->lo, ta));
    return best;
}

} // namespace

// ================================================================ GafIndex

struct GafLoadStats {
    uint64_t lines = 0, comments = 0, records = 0, secondary = 0, unmapped = 0, foreign = 0, projected = 0;
    uint64_t failed = 0, f_name = 0, f_misaligned = 0, f_overrun = 0, f_gap = 0, f_empty = 0, f_no_link = 0, f_length = 0;
    uint64_t handles = 0, alt_handles = 0, occurrences = 0;
    uint64_t partial_walks = 0, partial_handles = 0, no_ref_records = 0;
    std::string first_failure;
};

class GafIndex {
public:
    GafIndex(const Graph& g, const std::string& path, const Options& opt);
    ~GafIndex();
    GafIndex(const GafIndex&) = delete;
    GafIndex& operator=(const GafIndex&) = delete;

    const Graph& g;
    std::string path;

    // contigs (GAF query names, canonical), haplotypes
    std::vector<std::string> contig;
    std::vector<uint32_t> contig_hap;
    std::vector<int32_t> contig_sn;          // the graph SN of the same name, -1 if none
    std::vector<uint32_t> contig_records;    // projected records (the design prototype's ctg_lines)
    std::vector<std::string> hap;
    std::vector<uint32_t> sn_contig;         // graph SN -> contig, NONE if absent from the GAF

    // distinct canonical observed excursions, in first-seen order
    std::vector<Handle> dep, arr;
    std::vector<uint32_t> alt_off;
    std::vector<Handle> alts;
    std::vector<uint32_t> own_off, own_contig;   // owners sorted by contig id (one per contig)
    std::vector<int64_t> own_qs;                 // smallest query start of that contig's records
    std::vector<uint32_t> weight;                // distinct haplotypes
    std::vector<uint8_t> kind;                   // Kind
    std::vector<int64_t> wlo, whi;

    // per node
    std::vector<uint32_t> node_exc_off, node_exc;   // excursions through each alt node (sorted, unique)
    std::vector<uint8_t> node_partial;              // touched by a partial walk end or an unprojectable record

    // per record (audit)
    std::vector<uint32_t> rec_contig;
    std::vector<uint32_t> rec_cov_off;
    std::vector<RefIv> rec_cov;                     // merged rank-0 coverage
    std::vector<uint32_t> node_rec_off, node_rec;   // records walking each alt node
    std::vector<uint32_t> ctg_cov_off;
    std::vector<RefIv> ctg_cov;

    GafLoadStats st;
    std::string method;
    double seconds = 0;

    size_t n_exc() const { return dep.size(); }
    bool anchored_in(const Site& s, uint32_t id) const {
        return s.is_anchor_node(g, handle_node(dep[id])) && s.is_anchor_node(g, handle_node(arr[id]));
    }
    // the excursion as a foundation Excursion with its GAF owners (src GAF)
    Excursion make_excursion(uint32_t id) const;
    bool key_equal(uint32_t id, const Excursion& e) const {
        if (dep[id] != e.dep || arr[id] != e.arr || alt_off[id + 1] - alt_off[id] != e.alts.size()) return false;
        return std::equal(e.alts.begin(), e.alts.end(), alts.begin() + alt_off[id]);
    }
    bool has_contig(uint32_t id, uint32_t c) const {
        return std::binary_search(own_contig.begin() + own_off[id], own_contig.begin() + own_off[id + 1], c);
    }

private:
    std::vector<std::vector<NodeId>> sn_nodes_;     // every SN: its nodes sorted by SO
    std::unordered_map<std::string, int32_t> canon_sn_, last_sn_, name_cache_;
    std::unique_ptr<AtomicWriter> dump_exc_, dump_ctg_;

    int32_t resolve_sn(const std::string& name);
    int project_interval(int32_t sn, int64_t lo, int64_t hi, bool rev, std::vector<Handle>& hs) const;
    int project_path(const char* p, size_t n, int64_t ps, int64_t pe, std::vector<Handle>& hs, uint32_t& resolved, bool& bare);
    void touch_path(const char* p, size_t n, int64_t ps, int64_t pe);
    void touch_interval(int32_t sn, int64_t lo, int64_t hi);
    void write_dumps(const std::string& dir);
public:
    // the finished --dump writers (gaf_excursions.tsv, gaf_contigs.tsv), committed by main
    std::vector<AtomicWriter*> dump_writers() const;
};

namespace {
enum ProjResult { P_OK = 0, P_NAME, P_MISALIGNED, P_OVERRUN, P_GAP, P_EMPTY, P_NO_LINK, P_LENGTH };
const char* proj_name(int r) {
    switch (r) {
        case P_NAME: return "unknown sequence or segment name";
        case P_MISALIGNED: return "step does not start on a node boundary";
        case P_OVERRUN: return "step does not end on a node boundary";
        case P_GAP: return "step crosses a gap between nodes";
        case P_EMPTY: return "empty step";
        case P_NO_LINK: return "consecutive steps are not linked in the graph";
        case P_LENGTH: return "path length differs from column 7";
        default: return "ok";
    }
}
} // namespace

int32_t GafIndex::resolve_sn(const std::string& name) {
    auto it = name_cache_.find(name);
    if (it != name_cache_.end()) return it->second;
    int32_t sn = g.find_sn(name);
    if (sn < 0) {
        auto c = canon_sn_.find(canon_name(name));
        if (c != canon_sn_.end()) sn = c->second;
    }
    if (sn < 0 && name.find('#') == std::string::npos && name.find('|') == std::string::npos) {
        auto c = last_sn_.find(name);
        if (c != last_sn_.end()) sn = c->second;
    }
    if (sn < -1) sn = -1;    // ambiguous canonical or last-field match
    name_cache_[name] = sn;
    return sn;
}

int GafIndex::project_interval(int32_t sn, int64_t lo, int64_t hi, bool rev, std::vector<Handle>& hs) const {
    if (hi <= lo) return P_EMPTY;
    const std::vector<NodeId>& v = sn_nodes_[(size_t)sn];
    auto it = std::upper_bound(v.begin(), v.end(), lo, [&](int64_t x, NodeId n) { return x < g.start(n); });
    if (it == v.begin()) return P_MISALIGNED;
    --it;
    if (g.start(*it) != lo) return P_MISALIGNED;
    size_t first = hs.size();
    int64_t pos = lo;
    while (pos < hi) {
        if (it == v.end()) return P_OVERRUN;
        if (g.start(*it) != pos) return P_GAP;
        hs.push_back(make_handle(*it, false));
        pos = g.end(*it);
        ++it;
    }
    if (pos != hi) return P_OVERRUN;
    if (rev) {
        std::reverse(hs.begin() + (std::ptrdiff_t)first, hs.end());
        for (size_t k = first; k < hs.size(); ++k) hs[k] = flip(hs[k]);
    }
    return P_OK;
}

// Project a GAF path onto handles.  `resolved` counts the steps whose name is known to this graph
// (0: the record belongs to another graph).  A bare sequence name (no orientation) is the '+'
// path over that sequence, of which [ps, pe) is aligned.
int GafIndex::project_path(const char* p, size_t n, int64_t ps, int64_t pe, std::vector<Handle>& hs, uint32_t& resolved, bool& bare) {
    hs.clear();
    resolved = 0;
    bare = false;
    if (n == 0) return P_EMPTY;
    if (p[0] != '>' && p[0] != '<') {
        bare = true;
        std::string name(p, n);
        int32_t sn = resolve_sn(name);
        if (sn >= 0) {
            ++resolved;
            if (pe <= ps) return P_EMPTY;
            const std::vector<NodeId>& v = sn_nodes_[(size_t)sn];
            auto it = std::upper_bound(v.begin(), v.end(), ps, [&](int64_t x, NodeId y) { return x < g.start(y); });
            if (it == v.begin() || g.end(*(it - 1)) <= ps) return P_GAP;
            --it;
            int64_t pos = g.start(*it);
            while (pos < pe) {
                if (it == v.end()) return P_OVERRUN;
                if (g.start(*it) != pos) return P_GAP;
                hs.push_back(make_handle(*it, false));
                pos = g.end(*it);
                ++it;
            }
            return P_OK;
        }
        NodeId x = g.find_name(name);
        if (x == NONE) return P_NAME;
        ++resolved;
        hs.push_back(make_handle(x, false));
        return P_OK;
    }
    int result = P_OK;
    size_t i = 0;
    while (i < n) {
        bool rev = p[i] == '<';
        size_t j = i + 1;
        while (j < n && p[j] != '>' && p[j] != '<') ++j;
        const char* t = p + i + 1;
        size_t tl = j - i - 1;
        i = j;
        if (tl == 0) { if (result == P_OK) result = P_EMPTY; continue; }
        // after a failure, names are still resolved, so that a record that touches this graph
        // anywhere counts as failed rather than as belonging to another graph
        size_t nl = 0;
        int64_t lo = 0, hi = 0;
        if (parse_stable(t, tl, nl, lo, hi)) {
            int32_t sn = resolve_sn(std::string(t, nl));
            if (sn >= 0) {
                ++resolved;
                if (result == P_OK) result = project_interval(sn, lo, hi, rev, hs);
                continue;
            }
        }
        NodeId x = g.find_name(std::string(t, tl));
        if (x == NONE) { if (result == P_OK) result = P_NAME; continue; }
        ++resolved;
        if (result == P_OK) hs.push_back(make_handle(x, rev));
    }
    return result;
}

void GafIndex::touch_interval(int32_t sn, int64_t lo, int64_t hi) {
    const std::vector<NodeId>& v = sn_nodes_[(size_t)sn];
    auto it = std::upper_bound(v.begin(), v.end(), lo, [&](int64_t x, NodeId n) { return x < g.start(n); });
    if (it != v.begin()) --it;
    for (; it != v.end() && g.start(*it) < hi; ++it)
        if (g.end(*it) > lo && !g.is_ref(*it)) node_partial[*it] = 1;
}

// best effort for a record that could not be projected: every alt node it may walk is "partial"
void GafIndex::touch_path(const char* p, size_t n, int64_t ps, int64_t pe) {
    if (n == 0) return;
    if (p[0] != '>' && p[0] != '<') {
        std::string name(p, n);
        int32_t sn = resolve_sn(name);
        if (sn >= 0) touch_interval(sn, ps, pe);
        else {
            NodeId x = g.find_name(name);
            if (x != NONE && !g.is_ref(x)) node_partial[x] = 1;
        }
        return;
    }
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && p[j] != '>' && p[j] != '<') ++j;
        const char* t = p + i + 1;
        size_t tl = j - i - 1;
        i = j;
        size_t nl = 0;
        int64_t lo = 0, hi = 0;
        if (tl > 0 && parse_stable(t, tl, nl, lo, hi)) {
            int32_t sn = resolve_sn(std::string(t, nl));
            if (sn >= 0) { touch_interval(sn, lo, hi); continue; }
        }
        NodeId x = g.find_name(std::string(t, tl));
        if (x != NONE && !g.is_ref(x)) node_partial[x] = 1;
    }
}

GafIndex::GafIndex(const Graph& g_, const std::string& path_, const Options& opt) : g(g_), path(path_) {
    auto t0 = std::chrono::steady_clock::now();
    long rss0 = self_peak_rss_kb();
    const size_t N = g.nodes.size();
    const size_t NS = g.sn_names.size();
    // every SN's nodes by SO, and the name maps (-2 marks an ambiguous canonical / last-field name)
    sn_nodes_.assign(NS, std::vector<NodeId>());
    for (NodeId n = 0; n < (NodeId)N; ++n)
        if (g.nodes[n].sn >= 0) sn_nodes_[(size_t)g.nodes[n].sn].push_back(n);
    for (std::vector<NodeId>& v : sn_nodes_)
        std::sort(v.begin(), v.end(), [&](NodeId a, NodeId b) { return g.start(a) != g.start(b) ? g.start(a) < g.start(b) : a < b; });
    for (size_t s = 0; s < NS; ++s) {
        std::string c = canon_name(g.sn_names[s]);
        auto ins = canon_sn_.insert(std::make_pair(c, (int32_t)s));
        if (!ins.second && ins.first->second != (int32_t)s) ins.first->second = -2;
        std::string lf = last_field(c);
        auto ins2 = last_sn_.insert(std::make_pair(lf, (int32_t)s));
        if (!ins2.second && ins2.first->second != (int32_t)s) ins2.first->second = -2;
    }
    node_partial.assign(N, 0);

    std::unordered_map<std::string, uint32_t> contig_id, hap_id;
    std::unordered_map<std::string, uint32_t> exc_id;
    struct Occ { uint32_t exc, contig; int64_t qs; };
    std::vector<Occ> occ;
    std::vector<std::pair<uint32_t, uint32_t>> node_rec_pairs;   // (alt node, record)
    alt_off.push_back(0);
    rec_cov_off.push_back(0);

    GafReader rd(path, opt.threads);
    method = rd.method();
    std::string head, key;
    std::vector<size_t> fb, fe;          // field begin/end offsets in head
    std::vector<Handle> hs;
    std::vector<RefIv> cov;
    Excursion e;
    while (rd.next(head)) {
        ++st.lines;
        if (head.empty()) continue;
        if (head[0] == '#') { ++st.comments; continue; }
        fb.clear();
        fe.clear();
        size_t q = 0;
        while (true) {
            size_t t = head.find('\t', q);
            fb.push_back(q);
            fe.push_back(t == std::string::npos ? head.size() : t);
            if (t == std::string::npos) break;
            q = t + 1;
        }
        const uint64_t ln = rd.line_number();
        if (fb.size() < 12)
            fail(EXIT_INPUT, strf("GAF %s line %llu has %zu tab-separated fields; a GAF line has at least 12", path.c_str(),
                                  (unsigned long long)ln, fb.size()));
        const char* path_p = head.data() + fb[5];
        const size_t path_n = fe[5] - fb[5];
        if (path_n == 1 && path_p[0] == '*') { ++st.unmapped; continue; }
        auto field_i64 = [&](int k, const char* what) {
            int64_t v = 0;
            if (!parse_i64(head.data() + fb[(size_t)k], fe[(size_t)k] - fb[(size_t)k], v))
                fail(EXIT_INPUT, strf("GAF %s line %llu: column %d (%s) is not a non-negative integer", path.c_str(),
                                      (unsigned long long)ln, k + 1, what));
            return v;
        };
        const int64_t qlen = field_i64(1, "query length"), qs = field_i64(2, "query start"), qe = field_i64(3, "query end");
        const int64_t plen = field_i64(6, "path length"), ps = field_i64(7, "path start"), pe = field_i64(8, "path end");
        field_i64(9, "residue matches");
        field_i64(10, "block length");
        field_i64(11, "mapping quality");
        const std::string strand = head.substr(fb[4], fe[4] - fb[4]);
        if (strand != "+" && strand != "-" && strand != "*")
            fail(EXIT_INPUT, strf("GAF %s line %llu: column 5 (strand) is '%s'", path.c_str(), (unsigned long long)ln, strand.c_str()));
        if (qs > qe || qe > qlen || ps > pe || pe > plen)
            fail(EXIT_INPUT, strf("GAF %s line %llu: inconsistent query or path interval", path.c_str(), (unsigned long long)ln));
        bool secondary = false;
        for (size_t k = 12; k < fb.size(); ++k)
            if (fe[k] - fb[k] == 6 && head.compare(fb[k], 5, "tp:A:") == 0 && head[fb[k] + 5] == 'S') secondary = true;
        if (secondary) { ++st.secondary; continue; }
        ++st.records;
        uint32_t resolved = 0;
        bool bare = false;
        int r = project_path(path_p, path_n, ps, pe, hs, resolved, bare);
        if (r == P_OK && !bare) {
            int64_t sum = 0;
            for (Handle h : hs) sum += g.len(handle_node(h));
            if (sum != plen) r = P_LENGTH;
        }
        if (r == P_OK) {
            for (size_t k = 1; k < hs.size() && r == P_OK; ++k) {
                bool linked = false;
                for (const Edge& ed : g.out(hs[k - 1]))
                    if (ed.to == hs[k]) { linked = true; break; }
                if (!linked) r = P_NO_LINK;
            }
        }
        if (r != P_OK) {
            if (resolved == 0) { ++st.foreign; continue; }
            ++st.failed;
            switch (r) {
                case P_NAME: ++st.f_name; break;
                case P_MISALIGNED: ++st.f_misaligned; break;
                case P_OVERRUN: ++st.f_overrun; break;
                case P_GAP: ++st.f_gap; break;
                case P_EMPTY: ++st.f_empty; break;
                case P_NO_LINK: ++st.f_no_link; break;
                default: ++st.f_length; break;
            }
            if (st.first_failure.empty())
                st.first_failure = strf("line %llu (%s): %s", (unsigned long long)ln,
                                        head.substr(fb[0], fe[0] - fb[0]).c_str(), proj_name(r));
            touch_path(path_p, path_n, ps, pe);
            continue;
        }
        ++st.projected;
        // contig and haplotype
        std::string cname = canon_name(head.substr(fb[0], fe[0] - fb[0]));
        uint32_t cid;
        auto ci = contig_id.find(cname);
        if (ci != contig_id.end()) cid = ci->second;
        else {
            cid = (uint32_t)contig.size();
            contig_id.emplace(cname, cid);
            std::string hn = hap_name(cname);
            auto hit = hap_id.find(hn);
            uint32_t hid;
            if (hit != hap_id.end()) hid = hit->second;
            else { hid = (uint32_t)hap.size(); hap_id.emplace(hn, hid); hap.push_back(hn); }
            contig.push_back(cname);
            contig_hap.push_back(hid);
            contig_sn.push_back(g.find_sn(cname));
            contig_records.push_back(0);
        }
        ++contig_records[cid];
        const uint32_t rec = (uint32_t)rec_contig.size();
        rec_contig.push_back(cid);
        // reference coverage, alt-node visits, excursions
        st.handles += hs.size();
        cov.clear();
        int64_t last_ref = -1;
        for (size_t k = 0; k < hs.size(); ++k) {
            NodeId x = handle_node(hs[k]);
            if (!g.is_ref(x)) {
                ++st.alt_handles;
                node_rec_pairs.push_back(std::make_pair(x, rec));
                continue;
            }
            cov.push_back(RefIv{g.nodes[x].sn, g.start(x), g.end(x)});
            if (last_ref < 0) {
                if (k > 0) {    // partial walk start: alt handles before the first rank-0 handle
                    ++st.partial_walks;
                    for (size_t m = 0; m < k; ++m) { node_partial[handle_node(hs[m])] = 1; ++st.partial_handles; }
                }
            } else if ((int64_t)k > last_ref + 1) {
                e.dep = hs[(size_t)last_ref];
                e.arr = hs[k];
                e.alts.assign(hs.begin() + last_ref + 1, hs.begin() + (std::ptrdiff_t)k);
                canonicalize(e);
                key.assign((const char*)&e.dep, sizeof(Handle));
                key.append((const char*)&e.arr, sizeof(Handle));
                key.append((const char*)e.alts.data(), e.alts.size() * sizeof(Handle));
                auto ins = exc_id.insert(std::make_pair(key, (uint32_t)dep.size()));
                if (ins.second) {
                    dep.push_back(e.dep);
                    arr.push_back(e.arr);
                    alts.insert(alts.end(), e.alts.begin(), e.alts.end());
                    alt_off.push_back((uint32_t)alts.size());
                }
                occ.push_back(Occ{ins.first->second, cid, qs});
                ++st.occurrences;
            }
            last_ref = (int64_t)k;
        }
        if (last_ref < 0) {
            if (!hs.empty()) {
                ++st.no_ref_records;
                ++st.partial_walks;
                for (Handle h : hs) { node_partial[handle_node(h)] = 1; ++st.partial_handles; }
            }
        } else if ((size_t)last_ref + 1 < hs.size()) {   // partial walk end
            ++st.partial_walks;
            for (size_t m = (size_t)last_ref + 1; m < hs.size(); ++m) { node_partial[handle_node(hs[m])] = 1; ++st.partial_handles; }
        }
        merge_ivs(cov);
        rec_cov.insert(rec_cov.end(), cov.begin(), cov.end());
        rec_cov_off.push_back((uint32_t)rec_cov.size());
    }
    rd.close();
    exc_id.clear();
    name_cache_.clear();

    // does the GAF belong to this graph?
    const uint64_t touching = st.projected + st.failed;
    if (st.records > 0 && touching == 0)
        fail(EXIT_INPUT, strf("GAF %s: none of its %llu records names a sequence or segment of this graph", path.c_str(),
                              (unsigned long long)st.records));
    if (st.failed > 1 && st.failed * 100 > touching)
        fail(EXIT_INPUT, strf("GAF %s: %llu of %llu records that touch this graph cannot be projected onto it (first: %s); "
                              "the GAF must come from mapping to this (unzipped) graph",
                              path.c_str(), (unsigned long long)st.failed, (unsigned long long)touching, st.first_failure.c_str()));
    if (st.records == 0) ZLOG("warning: GAF %s has no mapped records: no site has observed walks", path.c_str());

    // owners (one per contig, smallest query start) and weights (distinct haplotypes)
    const size_t E = dep.size();
    std::sort(occ.begin(), occ.end(), [](const Occ& a, const Occ& b) {
        if (a.exc != b.exc) return a.exc < b.exc;
        if (a.contig != b.contig) return a.contig < b.contig;
        return a.qs < b.qs;
    });
    own_off.assign(1, 0);
    weight.assign(E, 0);
    {
        std::vector<uint32_t> haps;
        size_t i = 0;
        for (uint32_t x = 0; x < (uint32_t)E; ++x) {
            haps.clear();
            while (i < occ.size() && occ[i].exc == x) {
                if (own_contig.size() == own_off.back() || own_contig.back() != occ[i].contig) {
                    own_contig.push_back(occ[i].contig);
                    own_qs.push_back(occ[i].qs);
                    haps.push_back(contig_hap[occ[i].contig]);
                }
                ++i;
            }
            own_off.push_back((uint32_t)own_contig.size());
            std::sort(haps.begin(), haps.end());
            weight[x] = (uint32_t)(std::unique(haps.begin(), haps.end()) - haps.begin());
        }
    }
    std::vector<Occ>().swap(occ);
    // kinds and windows
    kind.assign(E, 0);
    wlo.assign(E, 0);
    whi.assign(E, 0);
    for (uint32_t x = 0; x < (uint32_t)E; ++x) {
        if (g.nodes[handle_node(dep[x])].sn != g.nodes[handle_node(arr[x])].sn) { kind[x] = KIND_X; continue; }
        Excursion t;
        t.dep = dep[x];
        t.arr = arr[x];
        kind[x] = (uint8_t)classify(g, t, wlo[x], whi[x]);
    }
    // node -> excursions
    node_exc_off.assign(N + 1, 0);
    {
        std::vector<NodeId> nodes;
        std::vector<std::pair<uint32_t, uint32_t>> pairs;
        for (uint32_t x = 0; x < (uint32_t)E; ++x) {
            nodes.clear();
            for (uint32_t k = alt_off[x]; k < alt_off[x + 1]; ++k) nodes.push_back(handle_node(alts[k]));
            std::sort(nodes.begin(), nodes.end());
            nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
            for (NodeId n : nodes) pairs.push_back(std::make_pair(n, x));
        }
        std::sort(pairs.begin(), pairs.end());
        for (const auto& pr : pairs) node_exc_off[pr.first + 1]++;
        for (size_t n = 0; n < N; ++n) node_exc_off[n + 1] += node_exc_off[n];
        node_exc.resize(pairs.size());
        for (size_t k = 0; k < pairs.size(); ++k) node_exc[k] = pairs[k].second;
    }
    // node -> records
    std::sort(node_rec_pairs.begin(), node_rec_pairs.end());
    node_rec_pairs.erase(std::unique(node_rec_pairs.begin(), node_rec_pairs.end()), node_rec_pairs.end());
    node_rec_off.assign(N + 1, 0);
    for (const auto& pr : node_rec_pairs) node_rec_off[pr.first + 1]++;
    for (size_t n = 0; n < N; ++n) node_rec_off[n + 1] += node_rec_off[n];
    node_rec.resize(node_rec_pairs.size());
    for (size_t k = 0; k < node_rec_pairs.size(); ++k) node_rec[k] = node_rec_pairs[k].second;
    std::vector<std::pair<uint32_t, uint32_t>>().swap(node_rec_pairs);
    // contig coverage
    {
        std::vector<std::vector<uint32_t>> recs(contig.size());
        for (uint32_t r = 0; r < (uint32_t)rec_contig.size(); ++r) recs[rec_contig[r]].push_back(r);
        ctg_cov_off.assign(1, 0);
        std::vector<RefIv> v;
        for (size_t c = 0; c < contig.size(); ++c) {
            v.clear();
            for (uint32_t r : recs[c]) v.insert(v.end(), rec_cov.begin() + rec_cov_off[r], rec_cov.begin() + rec_cov_off[r + 1]);
            merge_ivs(v);
            ctg_cov.insert(ctg_cov.end(), v.begin(), v.end());
            ctg_cov_off.push_back((uint32_t)ctg_cov.size());
        }
    }
    sn_contig.assign(NS, NONE);
    for (uint32_t c = 0; c < (uint32_t)contig.size(); ++c)
        if (contig_sn[c] >= 0) sn_contig[(size_t)contig_sn[c]] = c;

    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    long rss1 = self_peak_rss_kb();
    ZLOG("GAF %s (read with %s): %llu lines; %llu records: %llu projected, %llu for other graphs, %llu failed "
         "(name %llu, misaligned %llu, overrun %llu, gap %llu, empty %llu, no link %llu, length %llu); %llu secondary and "
         "%llu unmapped skipped",
         path.c_str(), method.c_str(), (unsigned long long)st.lines, (unsigned long long)st.records, (unsigned long long)st.projected,
         (unsigned long long)st.foreign, (unsigned long long)st.failed, (unsigned long long)st.f_name, (unsigned long long)st.f_misaligned,
         (unsigned long long)st.f_overrun, (unsigned long long)st.f_gap, (unsigned long long)st.f_empty, (unsigned long long)st.f_no_link,
         (unsigned long long)st.f_length, (unsigned long long)st.secondary, (unsigned long long)st.unmapped);
    if (!st.first_failure.empty()) ZLOG("GAF first unprojectable record: %s", st.first_failure.c_str());
    ZLOG("GAF: %zu contigs, %zu haplotypes; %llu excursion occurrences, %zu distinct; %llu partial walk ends (%llu alt handles); "
         "%.1f s, peak RSS %.2f -> %.2f GB",
         contig.size(), hap.size(), (unsigned long long)st.occurrences, E, (unsigned long long)st.partial_walks,
         (unsigned long long)st.partial_handles, seconds, rss0 / 1048576.0, rss1 / 1048576.0);
    if (!opt.dump_dir.empty()) write_dumps(opt.dump_dir);
}

GafIndex::~GafIndex() {
    // the dumps are committed by main with the other outputs (gaf_dump_writers); uncommitted, the
    // writers remove their temporary files
}

std::vector<AtomicWriter*> GafIndex::dump_writers() const {
    std::vector<AtomicWriter*> v;
    if (dump_exc_) v.push_back(dump_exc_.get());
    if (dump_ctg_) v.push_back(dump_ctg_.get());
    return v;
}

Excursion GafIndex::make_excursion(uint32_t id) const {
    Excursion e;
    e.dep = dep[id];
    e.arr = arr[id];
    e.alts.assign(alts.begin() + alt_off[id], alts.begin() + alt_off[id + 1]);
    e.owners.reserve(own_off[id + 1] - own_off[id]);
    for (uint32_t k = own_off[id]; k < own_off[id + 1]; ++k) {
        uint32_t c = own_contig[k];
        Owner o;
        o.contig = contig[c];
        o.sr = contig_sn[c] >= 0 ? g.sn_rank[(size_t)contig_sn[c]] : -1;
        o.so = own_qs[k];
        o.first = NONE;
        o.run = NONE;
        o.src = Src::GAF;
        e.owners.push_back(std::move(o));
    }
    e.weight = weight[id];
    e.src = Src::GAF;
    return e;
}

void GafIndex::write_dumps(const std::string& dir) {
    if (!make_dirs(dir)) fail(EXIT_IO, "cannot create --dump directory " + dir);
    std::vector<uint32_t> order(n_exc());
    std::iota(order.begin(), order.end(), 0u);
    auto key_less = [&](uint32_t a, uint32_t b) {
        if (dep[a] != dep[b]) return dep[a] < dep[b];
        if (arr[a] != arr[b]) return arr[a] < arr[b];
        return std::lexicographical_compare(alts.begin() + alt_off[a], alts.begin() + alt_off[a + 1], alts.begin() + alt_off[b],
                                            alts.begin() + alt_off[b + 1]);
    };
    std::sort(order.begin(), order.end(), key_less);
    dump_exc_.reset(new AtomicWriter(dir + "/gaf_excursions.tsv"));
    dump_exc_->write("kind\twlo\twhi\tdep\tarr\talts\tweight\tn_contigs\tcontigs\n");
    std::string line;
    for (uint32_t x : order) {
        line = strf("%s\t%lld\t%lld\t%s\t%s\t", kind[x] == KIND_X ? "X" : kind_name((Kind)kind[x]), (long long)wlo[x], (long long)whi[x],
                    g.handle_str(dep[x]).c_str(), g.handle_str(arr[x]).c_str());
        for (uint32_t k = alt_off[x]; k < alt_off[x + 1]; ++k) {
            if (k > alt_off[x]) line.push_back(',');
            line += g.handle_str(alts[k]);
        }
        line += strf("\t%u\t%u\t", weight[x], own_off[x + 1] - own_off[x]);
        std::vector<std::string> names;
        for (uint32_t k = own_off[x]; k < own_off[x + 1]; ++k) names.push_back(contig[own_contig[k]]);
        std::sort(names.begin(), names.end());
        for (size_t k = 0; k < names.size(); ++k) {
            if (k) line.push_back(',');
            line += names[k];
        }
        line.push_back('\n');
        dump_exc_->write(line);
    }
    dump_exc_->finish();
    dump_ctg_.reset(new AtomicWriter(dir + "/gaf_contigs.tsv"));
    dump_ctg_->write("contig\thaplotype\trecords\tgraph_sn\n");
    std::vector<uint32_t> co(contig.size());
    std::iota(co.begin(), co.end(), 0u);
    std::sort(co.begin(), co.end(), [&](uint32_t a, uint32_t b) { return contig[a] < contig[b]; });
    for (uint32_t c : co)
        dump_ctg_->write(strf("%s\t%s\t%u\t%s\n", contig[c].c_str(), hap[contig_hap[c]].c_str(), contig_records[c],
                              contig_sn[c] >= 0 ? g.sn_names[(size_t)contig_sn[c]].c_str() : "."));
    dump_ctg_->finish();
}

std::shared_ptr<const GafIndex> load_gaf(const Graph& g, const std::string& path, const Options& opt) {
    return std::shared_ptr<const GafIndex>(new GafIndex(g, path, opt));
}

// ================================================================ GafWalks

namespace {

// observed excursions of a site: through its alt nodes, both anchors on its reference nodes
std::vector<uint32_t> site_observed(const GafIndex& gi, const Site& s) {
    std::vector<uint32_t> ids;
    for (NodeId n : s.alts)
        for (uint32_t k = gi.node_exc_off[n]; k < gi.node_exc_off[n + 1]; ++k) ids.push_back(gi.node_exc[k]);
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    size_t m = 0;
    for (uint32_t id : ids)
        if (gi.anchored_in(s, id)) ids[m++] = id;
    ids.resize(m);
    return ids;
}

struct Counter {
    std::atomic<uint64_t> v{0};
    void add(uint64_t x) { if (x) v.fetch_add(x, std::memory_order_relaxed); }
    unsigned long long get() const { return (unsigned long long)v.load(); }
};

} // namespace

struct GafWalks::Impl {
    const Graph& g;
    const Runs& runs;
    std::shared_ptr<const GafIndex> gi;
    CreatorWalks creator;
    int64_t big = 5000;
    // summary counters (site threads add to them)
    Counter sites, observed, creator_matched_exc, runs_total;
    Counter al_observed, al_same, al_widened, al_opened, al_narrowed, al_closed, al_partial, al_foreign, al_unobserved;
    Counter al_blocked, al_empty, al_ok, al_partial_relax, al_partial_relax_bp;
    // creator reconstruction exactness (as the design prototype measured it): all sites / sites with <= 6000 alt nodes (the
    // the prototype's --max-alts); exact walk, and anchors only (the spec's 175/199 and 168/191)
    Counter cr_cmp, cr_exact, cr_anch, cr_cmp_big, cr_exact_big, cr_anch_big, cr_absent, cr_not_own;
    Counter cr6_cmp, cr6_exact, cr6_anch, cr6_cmp_big, cr6_exact_big, cr6_anch_big;
    // observed excursions of the site that are not creator excursions, by count and by contigs
    // (the prototype's obs-not-in-rf; the spec's 13-16%)
    Counter ob_exc, ob_exc_ctg, ob_notin, ob_notin_ctg, ob6_exc, ob6_exc_ctg, ob6_notin, ob6_notin_ctg;
    // rule (G): piece checks, records scanned, pieces whose node's records read the target, time
    Counter g_checks, g_records, g_hits, g_nanos;

    Impl(const Graph& g_, const Runs& r_, std::shared_ptr<const GafIndex> gi_, const Options& opt)
        : g(g_), runs(r_), gi(std::move(gi_)), creator(g_, r_), big(opt.b) {}
};

GafWalks::GafWalks(const Graph& g, const Runs& runs, const std::string& path, const Options& opt)
    : impl_(new Impl(g, runs, load_gaf(g, path, opt), opt)) {}

GafWalks::GafWalks(const Graph& g, const Runs& runs, std::shared_ptr<const GafIndex> index, const Options& opt)
    : impl_(new Impl(g, runs, std::move(index), opt)) {}

GafWalks::~GafWalks() {
    if (unwinding() || !impl_) return;
    Impl& m = *impl_;
    if (m.sites.get() == 0) return;
    ZLOG("gaf walks: %llu sites, %llu observed excursions (%llu equal a creator excursion)", m.sites.get(), m.observed.get(),
         m.creator_matched_exc.get());
    ZLOG("gaf walks: Allowed of %llu alt nodes crossed by observed excursions: ok %llu, empty %llu, blocked %llu; versus strict: "
         "%llu same, %llu widened, %llu opened, %llu narrowed, %llu closed; kept strict: %llu (partial walk; %llu of them, %llu bp, "
         "would have been widened or opened), %llu (foreign excursion); %llu alt nodes unobserved",
         m.al_observed.get(), m.al_ok.get(), m.al_empty.get(), m.al_blocked.get(), m.al_same.get(), m.al_widened.get(),
         m.al_opened.get(), m.al_narrowed.get(), m.al_closed.get(), m.al_partial.get(), m.al_partial_relax.get(),
         m.al_partial_relax_bp.get(), m.al_foreign.get(), m.al_unobserved.get());
    ZLOG("gaf walks: runs whose creator contig walks the run's first node: %llu; the creator excursion equals one of those observed "
         "excursions for %llu (anchors equal: %llu); runs >= %lld bp: %llu of %llu (anchors %llu of %llu); in sites with <= 6000 alt "
         "nodes: %llu (%llu) of %llu, >= %lld bp: %llu (%llu) of %llu; creator contig absent from the GAF %llu, not on its own node %llu",
         m.cr_cmp.get(), m.cr_exact.get(), m.cr_anch.get(), (long long)m.big, m.cr_exact_big.get(), m.cr_cmp_big.get(),
         m.cr_anch_big.get(), m.cr_cmp_big.get(), m.cr6_exact.get(), m.cr6_anch.get(), m.cr6_cmp.get(), (long long)m.big,
         m.cr6_exact_big.get(), m.cr6_anch_big.get(), m.cr6_cmp_big.get(), m.cr_absent.get(), m.cr_not_own.get());
    auto pct = [](unsigned long long a, unsigned long long b) { return b ? 100.0 * (double)a / (double)b : 0.0; };
    ZLOG("gaf walks: observed excursions that are not creator excursions: %llu of %llu (%.1f%%), %.1f%% weighted by contigs; in sites "
         "with <= 6000 alt nodes %llu of %llu (%.1f%%), %.1f%% by contigs",
         m.ob_notin.get(), m.ob_exc.get(), pct(m.ob_notin.get(), m.ob_exc.get()), pct(m.ob_notin_ctg.get(), m.ob_exc_ctg.get()),
         m.ob6_notin.get(), m.ob6_exc.get(), pct(m.ob6_notin.get(), m.ob6_exc.get()), pct(m.ob6_notin_ctg.get(), m.ob6_exc_ctg.get()));
    ZLOG("gaf walks: rule (G): %llu piece check(s), %llu GAF record(s) scanned, %.3f s; %llu piece(s) whose node a record walks that "
         "also reads the target",
         m.g_checks.get(), m.g_records.get(), (double)m.g_nanos.get() / 1e9, m.g_hits.get());
}

bool GafWalks::complete() const { return true; }

std::shared_ptr<const GafIndex> GafWalks::index() const { return impl_->gi; }

std::vector<AtomicWriter*> gaf_dump_writers(const GafWalks& w) { return w.index()->dump_writers(); }

std::vector<Excursion> GafWalks::excursions(const SiteData& sd, WalkStats& st) const {
    Impl& m = *impl_;
    const GafIndex& gi = *m.gi;
    const Graph& g = m.g;
    const Site& s = *sd.site;
    std::vector<uint32_t> ids = site_observed(gi, s);

    // creator excursions of the site, canonical and ordered by key
    WalkStats cst;
    std::vector<Excursion> cre = m.creator.excursions(sd, cst);
    std::vector<uint32_t> cidx;
    for (uint32_t i = 0; i < (uint32_t)cre.size(); ++i) {
        if (cre[i].src != Src::CREATOR) continue;
        canonicalize(cre[i]);
        cidx.push_back(i);
    }
    std::sort(cidx.begin(), cidx.end(), [&](uint32_t a, uint32_t b) {
        if (excursion_key_less(cre[a], cre[b])) return true;
        if (excursion_key_less(cre[b], cre[a])) return false;
        return a < b;
    });

    std::vector<Excursion> out;
    out.reserve(ids.size());
    std::vector<uint32_t> matched_runs;
    uint64_t matched_exc = 0, ob_exc = 0, ob_exc_ctg = 0, ob_notin = 0, ob_notin_ctg = 0;
    for (uint32_t id : ids) {
        Excursion e = gi.make_excursion(id);
        auto lo = std::lower_bound(cidx.begin(), cidx.end(), e, [&](uint32_t c, const Excursion& x) { return excursion_key_less(cre[c], x); });
        bool any = false;
        for (auto it = lo; it != cidx.end() && excursion_key_equal(cre[*it], e); ++it) {
            for (const Owner& o : cre[*it].owners) {
                e.owners.push_back(o);
                if (o.run != NONE) matched_runs.push_back(o.run);
            }
            any = true;
        }
        const uint64_t nctg = gi.own_off[id + 1] - gi.own_off[id];
        ++ob_exc;
        ob_exc_ctg += nctg;
        if (any) { e.src = Src::CREATOR; ++matched_exc; }
        else { ++ob_notin; ob_notin_ctg += nctg; }
        out.push_back(std::move(e));
    }
    std::sort(matched_runs.begin(), matched_runs.end());
    matched_runs.erase(std::unique(matched_runs.begin(), matched_runs.end()), matched_runs.end());
    st.runs += cst.runs;
    st.creator_ok += matched_runs.size();

    // creator reconstruction exactness, as the design prototype measured it: the run's creator contig has an
    // observed excursion through the run's first node; is the creator excursion one of them?
    const bool small_site = s.alts.size() <= 6000;
    uint64_t cmp = 0, exact = 0, anch = 0, cmp_big = 0, exact_big = 0, anch_big = 0, absent = 0, not_own = 0;
    for (uint32_t i = 0; i < (uint32_t)cre.size(); ++i) {
        if (cre[i].owners.empty()) continue;
        const Owner& ow = cre[i].owners[0];
        if (ow.run == NONE) continue;
        const std::vector<NodeId>& r = m.runs.runs[ow.run];
        int64_t nbp = 0;
        for (NodeId x : r) nbp += g.len(x);
        uint32_t c = gi.sn_contig[(size_t)g.nodes[r[0]].sn];
        if (c == NONE || gi.contig_records[c] == 0) { ++absent; continue; }
        bool own = false, eq = false, an = false;
        const bool ok = cre[i].src == Src::CREATOR;
        for (uint32_t k = gi.node_exc_off[r[0]]; k < gi.node_exc_off[r[0] + 1]; ++k) {
            uint32_t id = gi.node_exc[k];
            if (!gi.has_contig(id, c)) continue;
            own = true;
            if (ok && gi.dep[id] == cre[i].dep && gi.arr[id] == cre[i].arr) {
                an = true;
                if (gi.key_equal(id, cre[i])) eq = true;
            }
        }
        if (!own) { ++not_own; continue; }
        ++cmp;
        exact += eq;
        anch += an;
        if (nbp >= m.big) { ++cmp_big; exact_big += eq; anch_big += an; }
    }
    m.cr_cmp.add(cmp); m.cr_exact.add(exact); m.cr_anch.add(anch); m.cr_cmp_big.add(cmp_big); m.cr_exact_big.add(exact_big);
    m.cr_anch_big.add(anch_big); m.cr_absent.add(absent); m.cr_not_own.add(not_own);
    m.ob_exc.add(ob_exc); m.ob_exc_ctg.add(ob_exc_ctg); m.ob_notin.add(ob_notin); m.ob_notin_ctg.add(ob_notin_ctg);
    if (small_site) {
        m.cr6_cmp.add(cmp); m.cr6_exact.add(exact); m.cr6_anch.add(anch); m.cr6_cmp_big.add(cmp_big);
        m.cr6_exact_big.add(exact_big); m.cr6_anch_big.add(anch_big);
        m.ob6_exc.add(ob_exc); m.ob6_exc_ctg.add(ob_exc_ctg); m.ob6_notin.add(ob_notin); m.ob6_notin_ctg.add(ob_notin_ctg);
    }
    m.sites.add(1);
    m.observed.add(out.size());
    m.creator_matched_exc.add(matched_exc);
    m.runs_total.add(cst.runs);
    return out;
}

void GafWalks::refine_allowed(const SiteData& sd, std::vector<AllowedIv>& allowed) const {
    Impl& m = *impl_;
    const GafIndex& gi = *m.gi;
    const Site& s = *sd.site;
    uint64_t n_obs = 0, same = 0, widened = 0, opened = 0, narrowed = 0, closed = 0, partial = 0, foreign = 0, unobserved = 0;
    uint64_t n_ok = 0, n_empty = 0, n_blocked = 0, partial_relax = 0, partial_relax_bp = 0;
    for (uint32_t l = 0; l < (uint32_t)sd.sg.n_nodes(); ++l) {
        const NodeId n = sd.sg.alt[l];
        const uint32_t k0 = gi.node_exc_off[n], k1 = gi.node_exc_off[n + 1];
        if (k0 == k1) { ++unobserved; continue; }
        bool frn = false, junction = false, insertion = false, plus = false, minus = false;
        int64_t lo = INT64_MIN, hi = INT64_MAX;
        for (uint32_t k = k0; k < k1 && !frn; ++k) {
            const uint32_t id = gi.node_exc[k];
            if (!gi.anchored_in(s, id)) { frn = true; break; }
            const Kind kd = (Kind)gi.kind[id];
            if (kd == Kind::J) { junction = true; continue; }
            if (kd != Kind::F) insertion = true;
            lo = std::max(lo, gi.wlo[id]);
            hi = std::min(hi, gi.whi[id]);
            for (uint32_t a = gi.alt_off[id]; a < gi.alt_off[id + 1]; ++a)
                if (handle_node(gi.alts[a]) == n) {
                    if (handle_rev(gi.alts[a])) minus = true;
                    else plus = true;
                }
        }
        if (frn) { ++foreign; continue; }
        AllowedIv iv;
        iv.observed = true;
        if (junction) {
            iv.status = AllowedStatus::BLOCKED;
        } else {
            iv.lo = lo;
            iv.hi = hi;
            iv.minus = minus && !plus;
            iv.status = (!insertion && hi > lo) ? AllowedStatus::OK : AllowedStatus::EMPTY;
        }
        const AllowedIv& strict = allowed[l];
        if (gi.node_partial[n]) {
            // a walk ends inside this node (or an unprojectable record may touch it): its window is
            // unknown, so the node keeps the strict value; count what that costs
            ++partial;
            if (iv.status == AllowedStatus::OK &&
                (strict.status != AllowedStatus::OK || iv.lo < strict.lo || iv.hi > strict.hi)) {
                ++partial_relax;
                partial_relax_bp += (uint64_t)m.g.len(n);
            }
            continue;
        }
        ++n_obs;
        const bool so = strict.status == AllowedStatus::OK, oo = iv.status == AllowedStatus::OK;
        if (so && oo) {
            if (iv.lo > strict.lo || iv.hi < strict.hi) ++narrowed;
            else if (iv.lo < strict.lo || iv.hi > strict.hi) ++widened;
            else ++same;
        } else if (!so && oo) ++opened;
        else if (so && !oo) ++closed;
        else ++same;
        if (iv.status == AllowedStatus::OK) ++n_ok;
        else if (iv.status == AllowedStatus::EMPTY) ++n_empty;
        else ++n_blocked;
        allowed[l] = iv;
    }
    m.al_observed.add(n_obs); m.al_same.add(same); m.al_widened.add(widened); m.al_opened.add(opened);
    m.al_narrowed.add(narrowed); m.al_closed.add(closed); m.al_partial.add(partial); m.al_foreign.add(foreign);
    m.al_unobserved.add(unobserved); m.al_ok.add(n_ok); m.al_empty.add(n_empty); m.al_blocked.add(n_blocked);
    m.al_partial_relax.add(partial_relax); m.al_partial_relax_bp.add(partial_relax_bp);
}

bool GafWalks::reads_target(const Site& s, NodeId n, int64_t piece_bp, int64_t ta, int64_t tb) const {
    Impl& m = *impl_;
    const GafIndex& gi = *m.gi;
    // the audit's reading rule: pieces >= GAF_MIN_READ_PIECE, overlap > GAF_READ_TOL (impossible
    // when the target itself is not longer than that)
    if (piece_bp < GAF_MIN_READ_PIECE || tb - ta <= GAF_READ_TOL || (size_t)n + 1 >= gi.node_rec_off.size()) return false;
    const auto t0 = std::chrono::steady_clock::now();
    bool hit = false;
    uint32_t k = gi.node_rec_off[n];
    for (; k < gi.node_rec_off[n + 1] && !hit; ++k) {
        const uint32_t r = gi.node_rec[k];
        hit = max_overlap(gi.rec_cov.data() + gi.rec_cov_off[r], gi.rec_cov.data() + gi.rec_cov_off[r + 1], s.sn, ta, tb) > GAF_READ_TOL;
    }
    m.g_checks.add(1);
    m.g_records.add(k - gi.node_rec_off[n]);
    if (hit) m.g_hits.add(1);
    m.g_nanos.add((uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t0).count());
    return hit;
}

// ================================================================ GafAudit

void GafAuditStats::add(const GafAuditStats& o) {
    units += o.units; pieces += o.pieces; bp += o.bp;
    viol_units += o.viol_units; viol_pieces += o.viol_pieces; viol_bp += o.viol_bp; ok_bp += o.ok_bp;
    unobserved_bp += o.unobserved_bp; max_bad_haps = std::max(max_bad_haps, o.max_bad_haps);
    read_units += o.read_units; read_lines += o.read_lines; read_bp += o.read_bp;
    ctg_read_units += o.ctg_read_units; ctg_read_bp += o.ctg_read_bp;
    alt_pieces += o.alt_pieces; alt_bp += o.alt_bp;
}

std::string GafAuditStats::summary() const {
    return strf("GAF audit: %llu chains, %llu pieces, %lld bp: %lld bp observed-safe, %lld bp unobserved, %lld bp in %llu chains "
                "violating an observed window (max %u haplotypes); %llu GAF lines in %llu chains (%lld bp) also read the target; "
                "%llu chains (%lld bp) whose contig reads the target in another line; %llu alt-vs-alt pieces (%lld bp) not "
                "audited; %s",
                (unsigned long long)units, (unsigned long long)pieces, (long long)bp, (long long)ok_bp, (long long)unobserved_bp,
                (long long)viol_bp, (unsigned long long)viol_units, max_bad_haps, (unsigned long long)read_lines,
                (unsigned long long)read_units, (long long)read_bp, (unsigned long long)ctg_read_units, (long long)ctg_read_bp,
                (unsigned long long)alt_pieces, (long long)alt_bp, clean() ? "PASS" : "FAIL");
}

GafAudit::GafAudit(const Graph& g, const std::string& path, const Options& opt) : g_(g), gi_(load_gaf(g, path, opt)) {}

GafAudit::GafAudit(const Graph& g, std::shared_ptr<const GafIndex> index) : g_(g), gi_(std::move(index)) {}

GafAudit::~GafAudit() {}

GafAuditStats GafAudit::audit_plan(const Site& s, const std::vector<ZippedPiece>& zipped, std::string* detail) const {
    std::vector<GafAuditPiece> pieces;
    uint64_t alt_pieces = 0;
    int64_t alt_bp = 0;
    for (const ZippedPiece& z : zipped) {
        if (z.alt) { ++alt_pieces; alt_bp += z.b - z.a; continue; }
        GafAuditPiece p;
        p.unit = z.unit;
        p.node = z.node;
        p.a = z.a;
        p.b = z.b;
        p.ta = z.ta;
        p.tb = z.tb;
        pieces.push_back(p);
    }
    GafAuditStats st = audit_site(s, pieces, detail);
    st.alt_pieces += alt_pieces;
    st.alt_bp += alt_bp;
    return st;
}

GafAuditStats GafAudit::audit_site(const Site& s, const std::vector<GafAuditPiece>& pieces, std::string* detail) const {
    const GafIndex& gi = *gi_;
    GafAuditStats tot;
    std::vector<size_t> order(pieces.size());
    std::iota(order.begin(), order.end(), (size_t)0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return pieces[a].unit < pieces[b].unit; });
    const std::string lab = s.label(g_);
    std::vector<uint32_t> badh, rd_ctg, other_ctg;
    size_t i = 0;
    while (i < order.size()) {
        const uint32_t unit = pieces[order[i]].unit;
        GafAuditStats u;
        u.units = 1;
        uint64_t max_lines = 0, max_ctgs = 0;
        std::string pdet;
        for (; i < order.size() && pieces[order[i]].unit == unit; ++i) {
            const GafAuditPiece& p = pieces[order[i]];
            const int64_t len = p.b - p.a;
            ++u.pieces;
            u.bp += len;
            // window rule (the prototype audit's window check)
            badh.clear();
            bool any = false, bad = false;
            std::string why;
            if (p.node < gi.node_exc_off.size() - 1) {
                for (uint32_t k = gi.node_exc_off[p.node]; k < gi.node_exc_off[p.node + 1]; ++k) {
                    const uint32_t id = gi.node_exc[k];
                    bool in_site = gi.anchored_in(s, id);
                    if (!in_site) {
                        // impossible for a non-leaking site (the interior is closed under links); count it as unsafe
                        bad = true;
                        if (why.empty()) why = "foreign excursion " + g_.handle_str(gi.dep[id]) + ">" + g_.handle_str(gi.arr[id]);
                    } else {
                        any = true;
                        const Kind kd = (Kind)gi.kind[id];
                        bool ok = kd == Kind::F && gi.wlo[id] - window_tol <= p.ta && p.tb <= gi.whi[id] + window_tol;
                        if (ok) continue;
                        bad = true;
                        if (why.empty())
                            why = strf("%s excursion %s>%s window [%lld,%lld)", kind_name(kd), g_.handle_str(gi.dep[id]).c_str(),
                                       g_.handle_str(gi.arr[id]).c_str(), (long long)gi.wlo[id], (long long)gi.whi[id]);
                    }
                    for (uint32_t o = gi.own_off[id]; o < gi.own_off[id + 1]; ++o) badh.push_back(gi.contig_hap[gi.own_contig[o]]);
                }
            }
            std::sort(badh.begin(), badh.end());
            uint32_t nbad = (uint32_t)(std::unique(badh.begin(), badh.end()) - badh.begin());
            if (bad) {
                ++u.viol_pieces;
                u.viol_bp += len;
                u.max_bad_haps = std::max(u.max_bad_haps, nbad);
                if (detail)
                    pdet += strf("GAFAUDIT-PIECE\t%s\t%u\t%s[%lld,%lld)\t->[%lld,%lld)\t%s\thaps=%u\n", lab.c_str(), unit,
                                 g_.name(p.node).c_str(), (long long)p.a, (long long)p.b, (long long)p.ta, (long long)p.tb,
                                 why.c_str(), nbad);
            } else if (any) {
                u.ok_bp += len;
            } else {
                u.unobserved_bp += len;
            }
            // reading rule (the prototype audit's revisit check): GAF lines through the node that also read the target
            if (len >= min_read_piece && p.node < gi.node_rec_off.size() - 1) {
                rd_ctg.clear();
                other_ctg.clear();
                uint64_t lines = 0;
                for (uint32_t k = gi.node_rec_off[p.node]; k < gi.node_rec_off[p.node + 1]; ++k) {
                    const uint32_t r = gi.node_rec[k];
                    const uint32_t c = gi.rec_contig[r];
                    if (max_overlap(gi.rec_cov.data() + gi.rec_cov_off[r], gi.rec_cov.data() + gi.rec_cov_off[r + 1], s.sn, p.ta, p.tb) > read_tol) {
                        ++lines;
                        rd_ctg.push_back(c);
                    } else {
                        other_ctg.push_back(c);
                    }
                }
                std::sort(rd_ctg.begin(), rd_ctg.end());
                rd_ctg.erase(std::unique(rd_ctg.begin(), rd_ctg.end()), rd_ctg.end());
                std::sort(other_ctg.begin(), other_ctg.end());
                other_ctg.erase(std::unique(other_ctg.begin(), other_ctg.end()), other_ctg.end());
                uint64_t ctgs = rd_ctg.size();
                for (uint32_t c : other_ctg) {
                    if (std::binary_search(rd_ctg.begin(), rd_ctg.end(), c)) continue;
                    if (max_overlap(gi.ctg_cov.data() + gi.ctg_cov_off[c], gi.ctg_cov.data() + gi.ctg_cov_off[c + 1], s.sn, p.ta, p.tb) > read_tol)
                        ++ctgs;
                }
                max_lines = std::max(max_lines, lines);
                max_ctgs = std::max(max_ctgs, ctgs);
                u.read_lines += lines;
                if (lines && detail)
                    pdet += strf("GAFAUDIT-PIECE\t%s\t%u\t%s[%lld,%lld)\t->[%lld,%lld)\t%llu GAF lines also read the target\n", lab.c_str(),
                                 unit, g_.name(p.node).c_str(), (long long)p.a, (long long)p.b, (long long)p.ta, (long long)p.tb,
                                 (unsigned long long)lines);
            }
        }
        if (u.viol_bp > 0) u.viol_units = 1;
        if (max_lines > 0) { u.read_units = 1; u.read_bp = u.bp; }
        else if (max_ctgs > 0) { u.ctg_read_units = 1; u.ctg_read_bp = u.bp; }
        if (detail) {
            const char* obs = u.viol_bp > 0 ? "harm" : (u.ok_bp > 0 ? "safe" : "unobserved");
            *detail += strf("GAFAUDIT\t%s\t%u\tpieces=%llu\tbp=%lld\tobs=%s\tviol_bp=%lld\tok_bp=%lld\tnone_bp=%lld\tmax_bad_haps=%u\t"
                            "read_lines=%llu\tread_ctgs=%llu\n",
                            lab.c_str(), unit, (unsigned long long)u.pieces, (long long)u.bp, obs, (long long)u.viol_bp, (long long)u.ok_bp,
                            (long long)u.unobserved_bp, u.max_bad_haps, (unsigned long long)max_lines, (unsigned long long)max_ctgs);
            *detail += pdet;
        }
        tot.add(u);
    }
    return tot;
}

} // namespace zip
