/*
  zip_align.cpp -- spec step 4: the aligner runner, PAF/CIGAR verification, blocks, the left-snap
  projection P, feasibility cuts, the chain DP (frames F/R), the island rule, rule U and the
  fragmentation gate.  See zip_align.hpp for coordinates, outcomes and the --inject-chains format.

  Runner (per distinct window = per distinct target, written once as w.fa):
    1. screen (pass 1):      minimap2 -x <preset> -N 50 -p 0.01 --secondary=yes -t 1 [-X] w.fa q.fa
       on the first --screen-sample queries in canonical order; if at least half pass, every query
       goes to pass 2, otherwise every query is screened and the failures are 'prefiltered'.
       A query passes if its records' nmatch (column 10) sums to >= N and one record has
       nmatch / block length (column 11) >= F (--prefilter N,F).
    2. pass 2:               the same with -c --eqx.
    There is no index file: each process indexes w.fa in memory (0.2 s for a 5 Mb window), because
    minimap2 -d does not check its writes and a short index silently maps nothing.
    Each phase splits its queries over up to -j processes (largest first onto the lightest batch);
    per-query PAF does not depend on batching (checked: minimap2 2.30, batched vs single, .mmi vs
    FASTA target, byte-identical) as long as the query keeps its name, because minimap2 seeds its
    tie-breaking hash with the name: queries are named q<k> by their index in the window's query
    list, so results do not depend on -j or -t.
    Failures: a non-zero exit or unparseable output -> the window's units are 'aligner-failed'
    (stderr tail logged); a signal (OOM kill included) -> the window is retried once, alone (it
    takes the whole -j budget); a second failure -> ZipError(EXIT_ALIGNER).  minimap2 reporting
    that it could not write or read its files (a full --tmpdir) -> ZipError(EXIT_IO); a process
    that cannot be started (process limits) -> ZipError(EXIT_ALIGNER): outcomes never depend on
    transient resources.  check_failure_rate(): more than 5% of windows failed ->
    ZipError(EXIT_ALIGNER).  No timeouts.  Once an abort is requested (another thread's fatal
    error, a stop signal) nothing new starts and the running processes get SIGTERM.
    Timing (stderr only): every process records when it asked for a -j slot, when it got one and
    when it ended, and minimap2's own "Real time" and "CPU" (its last stderr line).  A window's log
    line splits its wall time into minimap2 running and waiting for a slot (with none of the
    window's processes running), and gives minimap2's own figures; the calls made on a thread add
    their windows to its AlignTimingScope (main: one per site).

  Verification: every record is checked column by column against both sequences, and minimap2's
  nmatch (column 10) and block length (column 11) against a recount.  Column 11 leaves out every
  ambiguous base: =/X columns with an ambiguous base on either side, ambiguous query bases of I
  runs and ambiguous target bases of D runs (a deletion over a reference N gap).

  Chain (one per unit, on feasible sub-records; mirrors the design's Python prototype):
    elements: blocks cut to feasible runs, each >= --min-piece, with its own exact '=' count;
    identical sub-records (same coordinates and ops, from different records) are kept once;
    sorted by (qs, qe, ts, strand, te, record, block); the 500 best by '=' are kept.
    DP per frame, O(n^2): a predecessor a of b needs a.qe <= b.qs + G and
      frame F: ++ collinear (b.ts >= a.te - G); -- collinear (b.te <= a.ts + G) and above the floor
               (b.ts >= floor - G); a strand switch starts a new segment after the current
               segment's hull (b.ts >= seg.hi - G; floor := seg.hi)
      frame R: mirrored (-- collinear; ++ collinear and below the floor b.te <= floor + G; a switch
               needs b.te <= seg.lo + G; floor := seg.lo)
    Every DP end's chain gets the island rule; a frame's chain is the best one after the island
    rule, then the best DP score before it, then the lower chain target start, then the lower end
    index.  The higher-scoring frame wins, ties to F.
    U: C* = best chain.  For every element x outside C* whose bound (best prefix to x + relaxed
    suffix) reaches (1 - delta) S*, C_x = best chain forced through x; C_x ties if its part not
    shared with C* scores >= (1 - delta) x C*'s part not shared with C_x.  Among C* and its ties,
    keep the fewest out-of-place stretches E; differing strand strings -> 'ambiguous-strand';
    otherwise phase ties (one stretch of the window covered by this unit's sub-records with no gap
    > G) take the lowest target start, then the lowest query start; copy ties take the highest
    score, then the lowest target start, then the lowest query start.
    Fragmentation: internal stretches >= G on either axis > max(2, frag x aligned / 5000) refuses.
*/
#include "zip_align.hpp"

#include <dirent.h>

#include <array>
#include <chrono>
#include <map>
#include <set>
#include <unordered_map>

namespace zip {

namespace {

const size_t MAX_CHAIN_PARTS = 500;     // spec: keep the 500 highest-scoring sub-records
const int REPEAT_K = 16;                // k-mer size of the repeat-fraction diagnostic

// ================================================================ bases

// minimap2's seq_nt4_table: A/C/G/T(U) -> 0..3, everything else (N, IUPAC) -> 4; case-insensitive
struct Nt4Table {
    uint8_t t[256];
    Nt4Table() {
        for (int i = 0; i < 256; ++i) t[i] = 4;
        t[(uint8_t)'A'] = t[(uint8_t)'a'] = 0;
        t[(uint8_t)'C'] = t[(uint8_t)'c'] = 1;
        t[(uint8_t)'G'] = t[(uint8_t)'g'] = 2;
        t[(uint8_t)'T'] = t[(uint8_t)'t'] = t[(uint8_t)'U'] = t[(uint8_t)'u'] = 3;
    }
};
const Nt4Table NT4;
inline uint8_t nt4(char c) { return NT4.t[(uint8_t)c]; }
inline uint8_t nt4_comp(char c) {
    uint8_t x = nt4(c);
    return x < 4 ? (uint8_t)(3 - x) : (uint8_t)4;
}

inline void push_op(std::vector<CigarOp>& ops, char op, int64_t len) {
    if (len <= 0) return;
    if (!ops.empty() && ops.back().op == op) ops.back().len += (uint32_t)len;
    else {
        CigarOp c;
        c.len = (uint32_t)len;
        c.op = op;
        ops.push_back(c);
    }
}

// "12=3X5D..." (ops =, X, I, D; M too when allow_m); adjacent equal ops are merged
bool parse_cigar(const std::string& s, bool allow_m, std::vector<CigarOp>& ops) {
    ops.clear();
    uint64_t len = 0;
    bool have = false;
    for (char c : s) {
        if (c >= '0' && c <= '9') {
            len = len * 10 + (uint64_t)(c - '0');
            if (len > 0x7fffffffULL) return false;
            have = true;
            continue;
        }
        if (!have || len == 0) return false;
        if (!(c == '=' || c == 'X' || c == 'I' || c == 'D' || (allow_m && c == 'M'))) return false;
        push_op(ops, c, (int64_t)len);
        len = 0;
        have = false;
    }
    return !have && !ops.empty();
}

std::string cigar_str(const std::vector<CigarOp>& ops) {
    std::string s;
    for (const CigarOp& o : ops) {
        s += std::to_string(o.len);
        s.push_back(o.op);
    }
    return s;
}

// ================================================================ files and directories

// RAII temporary directory <parent>/rgfa-zip.XXXXXX, removed with everything in it (registered
// with the Runtime, so a stop signal removes it too)
class TempDir {
public:
    explicit TempDir(const std::string& parent) {
        std::string tmpl = parent + "/rgfa-zip.XXXXXX";
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        if (!mkdtemp(buf.data())) fail(EXIT_IO, strf("cannot create a temporary directory under %s: %s", parent.c_str(), strerror(errno)));
        path_ = buf.data();
        Runtime::get().add_path(path_, true);
    }
    ~TempDir() {
        if (!path_.empty()) {
            remove_tree(path_);
            Runtime::get().remove_path(path_);
        }
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
    const std::string& path() const { return path_; }
private:
    std::string path_;
};

void write_file(const std::string& path, const std::string& content) {
    int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
    if (fd < 0) fail(EXIT_IO, strf("cannot create %s: %s", path.c_str(), strerror(errno)));
    const char* p = content.data();
    size_t n = content.size();
    while (n > 0) {
        ssize_t w = ::write(fd, p, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            int e = errno;
            ::close(fd);
            fail(EXIT_IO, strf("write %s: %s", path.c_str(), strerror(e)));
        }
        p += w;
        n -= (size_t)w;
    }
    if (::close(fd) != 0) fail(EXIT_IO, strf("close %s: %s", path.c_str(), strerror(errno)));
}

bool read_file(const std::string& path, std::string& out) {
    out.clear();
    int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd < 0) return false;
    char buf[1 << 16];
    while (true) {
        ssize_t r = ::read(fd, buf, sizeof(buf));
        if (r < 0) {
            if (errno == EINTR) continue;
            ::close(fd);
            return false;
        }
        if (r == 0) break;
        out.append(buf, (size_t)r);
    }
    ::close(fd);
    return true;
}

std::vector<std::string> split_ws(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ' ' || c == '\t' || c == '\n') {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

std::string sanitize(const std::string& s) {
    std::string o;
    for (char c : s) o.push_back((isalnum((unsigned char)c) || c == '.' || c == '_' || c == '-') ? c : '_');
    return o;
}

// ================================================================ the process budget (-j)

// At most `total` minimap2 processes at once.  acquire_all() takes the whole budget for a retry
// that must run alone; while it waits, ordinary acquires are held back so it cannot starve.
// Waiting stops with Aborted once an abort is requested (checked every 100 ms).
class ProcessBudget {
public:
    explicit ProcessBudget(int total) : total_(std::max(1, total)), free_(std::max(1, total)) {}
    void acquire() {
        std::unique_lock<std::mutex> lk(m_);
        while (!(free_ > 0 && !exclusive_)) {
            if (aborting()) throw Aborted();
            cv_.wait_for(lk, std::chrono::milliseconds(100));
        }
        if (aborting()) throw Aborted();
        --free_;
    }
    void release() {
        { std::lock_guard<std::mutex> lk(m_); ++free_; }
        cv_.notify_all();
    }
    void acquire_all() {
        excl_.lock();
        std::unique_lock<std::mutex> lk(m_);
        exclusive_ = true;
        while (free_ != total_ || aborting()) {
            if (aborting()) {
                exclusive_ = false;
                lk.unlock();
                cv_.notify_all();
                excl_.unlock();
                throw Aborted();
            }
            cv_.wait_for(lk, std::chrono::milliseconds(100));
        }
        free_ = 0;
    }
    void release_all() {
        {
            std::lock_guard<std::mutex> lk(m_);
            free_ = total_;
            exclusive_ = false;
        }
        cv_.notify_all();
        excl_.unlock();
    }
    int total() const { return total_; }
private:
    std::mutex m_, excl_;
    std::condition_variable cv_;
    const int total_;
    int free_;
    bool exclusive_ = false;
};

// minimap2's own peak RSS ("[M::main] Real time: ...; Peak RSS: 0.537 GB", GiB) in KB; 0 if absent
long mm2_peak_rss_kb(const std::string& tail) {
    size_t p = tail.rfind("Peak RSS: ");
    if (p == std::string::npos) return 0;
    double gb = strtod(tail.c_str() + p + 10, nullptr);
    return gb > 0 ? (long)(gb * 1048576.0 + 0.5) : 0;
}

// minimap2's own wall and CPU seconds ("[M::main] Real time: 0.537 sec; CPU: 0.524 sec; ..."); 0 if absent
void mm2_times(const std::string& tail, double& real, double& cpu) {
    real = cpu = 0;
    size_t p = tail.rfind("Real time: ");
    if (p == std::string::npos) return;
    real = strtod(tail.c_str() + p + 11, nullptr);
    size_t c = tail.find("CPU: ", p);
    if (c != std::string::npos) cpu = strtod(tail.c_str() + c + 5, nullptr);
    if (!(real > 0)) real = 0;
    if (!(cpu > 0)) cpu = 0;
}

// the sorted union of intervals [a, b)
std::vector<std::pair<double, double>> iv_union(std::vector<std::pair<double, double>> v) {
    std::sort(v.begin(), v.end());
    std::vector<std::pair<double, double>> out;
    for (const auto& x : v) {
        if (x.second <= x.first) continue;
        if (!out.empty() && x.first <= out.back().second) out.back().second = std::max(out.back().second, x.second);
        else out.push_back(x);
    }
    return out;
}

AlignTimingScope*& timing_top() {
    thread_local AlignTimingScope* top = nullptr;
    return top;
}

// minimap2's error lines in a stderr tail ("[ERROR] ..." / "ERROR: ..."), else its last line
std::string mm2_error_lines(const std::string& tail) {
    std::string out, last;
    size_t p = 0;
    while (p < tail.size()) {
        size_t e = tail.find('\n', p);
        if (e == std::string::npos) e = tail.size();
        std::string line = tail.substr(p, e - p);
        p = e + 1;
        if (line.empty()) continue;
        last = line;
        if (line.find("ERROR") != std::string::npos) out += (out.empty() ? "" : " | ") + line;
    }
    return out.empty() ? last : out;
}

// Does minimap2's stderr report that it could not write its output or read its input -- an I/O
// problem of --tmpdir (full disk, quota, a vanished directory), not an alignment failure?
bool mm2_io_error(const std::string& tail) {
    for (const char* m : {"failed to write", "failed to read data", "failed to open file", "failed to map the query file"})
        if (tail.find(m) != std::string::npos) return true;
    for (int e : {ENOSPC, EDQUOT, EIO, EROFS})
        if (tail.find(strerror(e)) != std::string::npos) return true;
    return false;
}

// ================================================================ PAF

struct RawRec {
    uint32_t q = 0;                 // query index within the window's query list
    int64_t qlen = 0, qs = 0, qe = 0, tlen = 0, ts = 0, te = 0, nmatch = 0, blen = 0;
    char strand = '+';
    bool has_cigar = false;
    std::string cigar;
};

// Parse minimap2's PAF.  Queries are named q<k>, k = the query's index in its window's query list
// (n_queries of them); in_batch[k] says which of them this process was given.  False on
// unparseable output (the window is then aligner-failed).
bool parse_paf(const std::string& text, size_t n_queries, const std::vector<uint8_t>& in_batch, bool want_cigar, std::vector<RawRec>& out,
               std::string& why) {
    size_t p = 0;
    std::vector<std::string> f;
    uint64_t lineno = 0;
    while (p < text.size()) {
        size_t e = text.find('\n', p);
        if (e == std::string::npos) {
            why = "the last PAF line has no newline (truncated output)";
            return false;
        }
        std::string line = text.substr(p, e - p);
        p = e + 1;
        ++lineno;
        if (line.empty()) { why = strf("empty PAF line %llu", (unsigned long long)lineno); return false; }
        split_tabs(line, f);
        if (f.size() < 12) { why = strf("PAF line %llu has %zu columns", (unsigned long long)lineno, f.size()); return false; }
        RawRec r;
        int64_t qi = -1;
        if (f[0].size() < 2 || f[0][0] != 'q' || !parse_int64(f[0].substr(1), qi) || qi < 0 || (size_t)qi >= n_queries ||
            !in_batch[(size_t)qi]) {
            why = strf("PAF line %llu: unknown query '%s'", (unsigned long long)lineno, f[0].c_str());
            return false;
        }
        r.q = (uint32_t)qi;
        if (f[5] != "t") { why = strf("PAF line %llu: unknown target '%s'", (unsigned long long)lineno, f[5].c_str()); return false; }
        if (f[4] != "+" && f[4] != "-") { why = strf("PAF line %llu: strand '%s'", (unsigned long long)lineno, f[4].c_str()); return false; }
        r.strand = f[4][0];
        if (!parse_int64(f[1], r.qlen) || !parse_int64(f[2], r.qs) || !parse_int64(f[3], r.qe) || !parse_int64(f[6], r.tlen) ||
            !parse_int64(f[7], r.ts) || !parse_int64(f[8], r.te) || !parse_int64(f[9], r.nmatch) || !parse_int64(f[10], r.blen)) {
            why = strf("PAF line %llu: a numeric column does not parse", (unsigned long long)lineno);
            return false;
        }
        for (size_t k = 12; k < f.size(); ++k)
            if (f[k].compare(0, 5, "cg:Z:") == 0) { r.cigar = f[k].substr(5); r.has_cigar = true; }
        if (want_cigar && !r.has_cigar) { why = strf("PAF line %llu has no cg:Z: tag", (unsigned long long)lineno); return false; }
        out.push_back(std::move(r));
    }
    return true;
}

// Verify a record against both sequences and its own coordinates, and build the PafRecord
// (target shifted by tlo; '=' columns over ambiguous bases rewritten as 'X').  check_counts:
// minimap2's nmatch and block length must equal the recount.  False (with why) if inconsistent.
//
// minimap2's block length (PAF column 11; mm_update_extra in align.c) leaves out every ambiguous
// (non-ACGT) base: in =/X columns where either base is ambiguous, in I runs the ambiguous query
// bases, in D runs the ambiguous target bases.  So a D run over a reference N gap is shorter in
// column 11 than in the CIGAR.  Only this consistency check sees those bases; nmatch, the record's
// '=' count (an N is never '=': minimap2's N-vs-N '=' columns become 'X') and the identity
// '=' / ('=' + 'X' + bases in I/D runs shorter than G) are as the spec defines them, Ns included.
bool verify_record(const RawRec& r, const std::string& qseq, const std::string& tseq, int64_t tlo, bool allow_m, bool check_counts,
                   int64_t G, PafRecord& out, std::string& why) {
    const int64_t qlen = (int64_t)qseq.size(), tlen = (int64_t)tseq.size();
    if (r.qlen != qlen) { why = strf("query length %lld, PAF says %lld", (long long)qlen, (long long)r.qlen); return false; }
    if (r.tlen != tlen) { why = strf("target length %lld, PAF says %lld", (long long)tlen, (long long)r.tlen); return false; }
    if (!(0 <= r.qs && r.qs < r.qe && r.qe <= qlen)) { why = strf("query interval %lld-%lld", (long long)r.qs, (long long)r.qe); return false; }
    if (!(0 <= r.ts && r.ts < r.te && r.te <= tlen)) { why = strf("target interval %lld-%lld", (long long)r.ts, (long long)r.te); return false; }
    std::vector<CigarOp> ops;
    if (!r.has_cigar || !parse_cigar(r.cigar, allow_m, ops)) { why = "CIGAR does not parse as =/X/I/D"; return false; }
    int64_t qc = 0, tc = 0;
    for (const CigarOp& o : ops) {
        if (o.op != 'D') qc += o.len;
        if (o.op != 'I') tc += o.len;
    }
    if (qc != r.qe - r.qs || tc != r.te - r.ts) {
        why = strf("CIGAR consumes %lld query / %lld target bases, the record spans %lld / %lld", (long long)qc, (long long)tc,
                   (long long)(r.qe - r.qs), (long long)(r.te - r.ts));
        return false;
    }
    out = PafRecord();
    out.qs = r.qs;
    out.qe = r.qe;
    out.ts = tlo + r.ts;
    out.te = tlo + r.te;
    out.strand = r.strand;
    const bool fwd = r.strand == '+';
    int64_t qi = fwd ? r.qs : r.qe - 1;   // next query base ('-': read downward, complemented)
    int64_t ti = r.ts;
    int64_t unamb_cols = 0, small = 0;
    int64_t amb_ins = 0, amb_del = 0;     // ambiguous query bases in I runs, target bases in D runs
    for (const CigarOp& o : ops) {
        const int64_t L = o.len;
        if (o.op == '=' || o.op == 'X' || o.op == 'M') {
            for (int64_t k = 0; k < L; ++k) {
                uint8_t cq = fwd ? nt4(qseq[(size_t)qi]) : nt4_comp(qseq[(size_t)qi]);
                uint8_t ct = nt4(tseq[(size_t)ti]);
                bool eq = cq == ct;
                if (o.op == '=' && !eq) {
                    why = strf("'=' column at query %lld / target %lld: %c vs %c%s", (long long)qi, (long long)(tlo + ti), qseq[(size_t)qi],
                               tseq[(size_t)ti], fwd ? "" : " (query complemented)");
                    return false;
                }
                if (o.op == 'X' && eq) {
                    why = strf("'X' column at query %lld / target %lld: %c vs %c%s", (long long)qi, (long long)(tlo + ti), qseq[(size_t)qi],
                               tseq[(size_t)ti], fwd ? "" : " (query complemented)");
                    return false;
                }
                bool unamb = cq < 4 && ct < 4;
                if (unamb) ++unamb_cols;
                bool match = eq && unamb;
                push_op(out.ops, match ? '=' : 'X', 1);
                if (match) ++out.n_eq;
                else ++out.n_x;
                qi += fwd ? 1 : -1;
                ++ti;
            }
        } else if (o.op == 'I') {
            // the run's query bases: [qi, qi + L) ('+'), or qi, qi - 1, ..., qi - L + 1 ('-')
            for (int64_t k = 0; k < L; ++k)
                if (nt4(qseq[(size_t)(fwd ? qi + k : qi - k)]) > 3) ++amb_ins;
            push_op(out.ops, 'I', L);
            out.n_ins += L;
            if (L < G) small += L;
            qi += fwd ? L : -L;
        } else {
            for (int64_t k = 0; k < L; ++k)
                if (nt4(tseq[(size_t)(ti + k)]) > 3) ++amb_del;
            push_op(out.ops, 'D', L);
            out.n_del += L;
            if (L < G) small += L;
            ti += L;
        }
    }
    if (check_counts) {
        if (r.nmatch != out.n_eq) {
            why = strf("nmatch %lld, recount %lld", (long long)r.nmatch, (long long)out.n_eq);
            return false;
        }
        const int64_t blen = unamb_cols + (out.n_ins - amb_ins) + (out.n_del - amb_del);
        if (r.blen != blen) {
            why = strf("block length %lld, recount %lld (%lld ambiguous base(s) in I runs, %lld in D runs left out)", (long long)r.blen,
                       (long long)blen, (long long)amb_ins, (long long)amb_del);
            return false;
        }
    }
    int64_t den = out.n_eq + out.n_x + small;
    out.gci = den > 0 ? (double)out.n_eq / (double)den : 0.0;
    return true;
}

bool record_less(const PafRecord& a, const PafRecord& b) {
    if (a.qs != b.qs) return a.qs < b.qs;
    if (a.qe != b.qe) return a.qe < b.qe;
    if (a.ts != b.ts) return a.ts < b.ts;
    if (a.te != b.te) return a.te < b.te;
    if (a.strand != b.strand) return a.strand < b.strand;
    if (a.ops.size() != b.ops.size()) return a.ops.size() < b.ops.size();
    for (size_t i = 0; i < a.ops.size(); ++i) {
        if (a.ops[i].op != b.ops[i].op) return a.ops[i].op < b.ops[i].op;
        if (a.ops[i].len != b.ops[i].len) return a.ops[i].len < b.ops[i].len;
    }
    return false;
}
bool record_equal(const PafRecord& a, const PafRecord& b) { return !record_less(a, b) && !record_less(b, a); }

// sort records canonically and drop exact duplicates, so nothing depends on minimap2's output order
void canonical_records(std::vector<PafRecord>& recs) {
    std::sort(recs.begin(), recs.end(), record_less);
    recs.erase(std::unique(recs.begin(), recs.end(), record_equal), recs.end());
}

// ================================================================ parts

double part_identity(const ChainPart& p, int64_t G) {
    int64_t eq = 0, x = 0, small = 0;
    for (const CigarOp& o : p.ops) {
        if (o.op == '=') eq += o.len;
        else if (o.op == 'X') x += o.len;
        else if ((int64_t)o.len < G) small += o.len;
    }
    int64_t den = eq + x + small;
    return den > 0 ? (double)eq / (double)den : 0.0;
}

// The sub-part [qa, qb) of a block: target [P(qa), P(qb)) ('+') or [P(qb), P(qa)) ('-'), and the
// ops between them (a D run at the start point is inside, one at the end point is not; this is
// exactly P's left-snap convention).  n_eq is recounted.
ChainPart slice_part(const ChainPart& blk, int64_t qa, int64_t qb) {
    ChainPart p;
    p.record = blk.record;
    p.block = blk.block;
    p.strand = blk.strand;
    p.qs = qa;
    p.qe = qb;
    const bool fwd = blk.strand == '+';
    int64_t c0 = fwd ? qa - blk.qs : blk.qe - qb;
    int64_t c1 = fwd ? qb - blk.qs : blk.qe - qa;
    std::vector<int64_t> P = project_points(blk, {qa, qb});
    p.ts = fwd ? P[0] : P[1];
    p.te = fwd ? P[1] : P[0];
    int64_t q = 0, tsum = 0;
    for (const CigarOp& o : blk.ops) {
        const int64_t L = o.len;
        if (o.op == 'D') {
            if (c0 <= q && q < c1) { push_op(p.ops, 'D', L); tsum += L; }
            continue;
        }
        int64_t lo = std::max(q, c0), hi = std::min(q + L, c1);
        if (hi > lo) {
            push_op(p.ops, o.op, hi - lo);
            if (o.op != 'I') tsum += hi - lo;
            if (o.op == '=') p.n_eq += hi - lo;
        }
        q += L;
    }
    if (tsum != p.te - p.ts)
        fail(EXIT_INVARIANT, strf("internal: slice [%lld,%lld) of a block spans %lld target bases, P gives %lld", (long long)qa,
                                  (long long)qb, (long long)tsum, (long long)(p.te - p.ts)));
    return p;
}

// handles of the walk overlapping [qs, qe): (k, a0, a1) in walk coordinates
void overlapping_handles(const QueryLayout& ql, int64_t qs, int64_t qe, std::vector<std::array<int64_t, 3>>& out) {
    out.clear();
    if (qe <= qs || ql.walk.empty()) return;
    size_t k = (size_t)(std::upper_bound(ql.off.begin(), ql.off.end() - 1, qs) - ql.off.begin());
    k = k == 0 ? 0 : k - 1;
    for (; k < ql.walk.size() && ql.off[k] < qe; ++k) {
        int64_t a0 = std::max(ql.off[k], qs), a1 = std::min(ql.off[k + 1], qe);
        if (a1 > a0) out.push_back({{(int64_t)k, a0, a1}});
    }
}

// Cut a block to its maximal runs of consecutive pieces whose targets are feasible; keep runs of
// at least min_piece query bases.
void feasible_parts(const ChainPart& blk, const QueryLayout& ql, const FeasibleFn& feasible, int64_t min_piece,
                    std::vector<ChainPart>& out) {
    std::vector<std::array<int64_t, 3>> hs;
    overlapping_handles(ql, blk.qs, blk.qe, hs);
    if (hs.empty()) return;
    std::vector<int64_t> xs;
    xs.reserve(hs.size() * 2);
    for (auto& h : hs) { xs.push_back(h[1]); xs.push_back(h[2]); }
    std::vector<int64_t> P = project_points(blk, xs);
    const bool fwd = blk.strand == '+';
    int64_t run_a = -1, run_b = -1;
    auto close = [&]() {
        if (run_a >= 0 && run_b - run_a >= min_piece) out.push_back(slice_part(blk, run_a, run_b));
        run_a = run_b = -1;
    };
    for (size_t i = 0; i < hs.size(); ++i) {
        NodeId node = handle_node(ql.walk[(size_t)hs[i][0]]);
        int64_t ta = fwd ? P[2 * i] : P[2 * i + 1];
        int64_t tb = fwd ? P[2 * i + 1] : P[2 * i];
        bool ok = !feasible || feasible(node, ta, tb);
        if (ok) {
            if (run_a < 0) run_a = hs[i][1];
            run_b = hs[i][2];
        } else {
            close();
        }
    }
    close();
}

bool part_less(const ChainPart& a, const ChainPart& b) {
    if (a.qs != b.qs) return a.qs < b.qs;
    if (a.qe != b.qe) return a.qe < b.qe;
    if (a.ts != b.ts) return a.ts < b.ts;
    if (a.strand != b.strand) return a.strand < b.strand;
    if (a.te != b.te) return a.te < b.te;
    if (a.record != b.record) return a.record < b.record;
    return a.block < b.block;
}

std::string strand_string(const std::vector<ChainPart>& parts) {
    std::string s;
    for (const ChainPart& p : parts)
        if (s.empty() || s.back() != p.strand) s.push_back(p.strand);
    return s;
}

// internal unaligned stretches >= G between consecutive parts (query order), on either axis
int64_t internal_stretches(const std::vector<ChainPart>& parts, int64_t G) {
    std::vector<const ChainPart*> v;
    for (const ChainPart& p : parts) v.push_back(&p);
    std::sort(v.begin(), v.end(), [](const ChainPart* a, const ChainPart* b) {
        if (a->qs != b->qs) return a->qs < b->qs;
        if (a->qe != b->qe) return a->qe < b->qe;
        if (a->ts != b->ts) return a->ts < b->ts;
        return a->te < b->te;
    });
    int64_t n = 0;
    for (size_t i = 1; i < v.size(); ++i) {
        const ChainPart& a = *v[i - 1];
        const ChainPart& b = *v[i];
        if (b.qs - a.qe >= G || std::max(a.ts, b.ts) - std::min(a.te, b.te) >= G) ++n;
    }
    return n;
}

// ================================================================ the chainer

struct Chain {
    bool valid = false;
    std::vector<uint32_t> el;          // element indices (canonical order = query order)
    int frame = 0;                     // 0 = F, 1 = R
    int64_t score = 0, aq = 0, hlo = 0, qlo = 0, thi = 0;
    std::vector<std::pair<int64_t, int64_t>> islands;
    int64_t island_bp = 0;
};

class Chainer {
public:
    Chainer(const std::vector<ChainPart>& c, const Options& opt) : c_(c), n_(c.size()), G_(opt.G), B_(opt.b), ISL_(opt.island) {
        for (int fr = 0; fr < 2; ++fr) {
            alloc(base_[fr]);
            alloc(work_[fr]);
            fill(fr, -1, 0, base_[fr]);
        }
    }

    // the best chain (forced through element `force` when >= 0), island rule applied; invalid if
    // no frame keeps >= b aligned query
    Chain best(int force) {
        Chain best;
        for (int fr = 0; fr < 2; ++fr) {
            DP* d = &base_[fr];
            if (force >= 0) {
                work_[fr] = base_[fr];
                fill(fr, force, (size_t)force + 1, work_[fr]);
                d = &work_[fr];
            }
            // every DP end is a candidate chain; the island rule runs on each, and the frame's chain is
            // the one with the best score AFTER the island rule, then the best DP score before it (so
            // the DP's own chain, and its island-below-b report, win whenever the island rule does not
            // change the winner), then the lower chain target start, then the lower end index.  Taking
            // the best end before the island rule can pick a chain whose islands are then stripped,
            // below a chain the DP also holds (chr15:101.75 VNTR).
            Chain bf;
            int64_t bf_dp = 0;
            std::vector<uint32_t> el;
            for (size_t j = force >= 0 ? (size_t)force : 0; j < n_; ++j) {
                if (!d->ok[j]) continue;
                // the island rule can only lower a chain's score: skip ends that cannot beat the best so far
                if (bf.valid && d->sc[j] < bf.score) continue;
                if (bf.valid && d->sc[j] == bf.score && bf_dp > d->sc[j]) continue;
                el.clear();
                for (int k = (int)j; k >= 0; k = d->prv[(size_t)k]) el.push_back((uint32_t)k);
                std::reverse(el.begin(), el.end());
                if (force >= 0 && !std::binary_search(el.begin(), el.end(), (uint32_t)force)) continue;
                Chain ch;
                ch.frame = fr;
                island_rule(el, ch);
                if (ch.el.empty()) continue;
                bool better = !bf.valid || ch.score > bf.score ||
                              (ch.score == bf.score && (d->sc[j] > bf_dp || (d->sc[j] == bf_dp && ch.hlo < bf.hlo)));
                if (better) {
                    bf = std::move(ch);
                    bf_dp = d->sc[j];
                }
            }
            if (!bf.valid) continue;
            if (!best.valid || bf.score > best.score) best = std::move(bf);
        }
        if (best.valid && best.aq < B_) best = Chain();
        return best;
    }

    // upper bound on the score of any chain through element x (for rule U's pruning)
    std::vector<int64_t> bounds() const {
        std::vector<int64_t> suf(n_, 0), out(n_, 0);
        for (size_t i = n_; i-- > 0;) {
            int64_t best = 0;
            for (size_t j = i + 1; j < n_; ++j)
                if (c_[i].qe <= c_[j].qs + G_) best = std::max(best, suf[j]);
            suf[i] = c_[i].n_eq + best;
        }
        for (size_t x = 0; x < n_; ++x) out[x] = std::max(base_[0].sc[x], base_[1].sc[x]) + suf[x] - c_[x].n_eq;
        return out;
    }

private:
    struct DP {
        std::vector<int64_t> sc, flo, slo, shi, hlo;
        std::vector<int32_t> prv;
        std::vector<uint8_t> ok, has_flo;
    };
    void alloc(DP& d) const {
        d.sc.assign(n_, 0); d.flo.assign(n_, 0); d.slo.assign(n_, 0); d.shi.assign(n_, 0); d.hlo.assign(n_, 0);
        d.prv.assign(n_, -1); d.ok.assign(n_, 0); d.has_flo.assign(n_, 0);
    }

    // DP entries j >= from (entries < from are already filled); with force >= 0 only chains that
    // contain c[force] (a chain may start at j only if j <= force, and no transition skips force)
    void fill(int fr, int force, size_t from, DP& d) const {
        const int64_t TOL = G_;
        for (size_t j = from; j < n_; ++j) {
            const ChainPart& b = c_[j];
            bool start = force < 0 || (int64_t)j <= force;
            d.ok[j] = start ? 1 : 0;
            d.sc[j] = start ? b.n_eq : -1;
            d.prv[j] = -1;
            d.has_flo[j] = 0;
            d.flo[j] = 0;
            d.slo[j] = b.ts;
            d.shi[j] = b.te;
            d.hlo[j] = b.ts;
            size_t i0 = (force >= 0 && (int64_t)j > force) ? (size_t)force : 0;
            for (size_t i = i0; i < j; ++i) {
                if (!d.ok[i]) continue;
                const ChainPart& a = c_[i];
                if (a.qe > b.qs + TOL) continue;
                bool same = false;
                if (fr == 0) {
                    if (a.strand == '+' && b.strand == '+') {
                        if (b.ts < a.te - TOL) continue;
                        same = true;
                    } else if (a.strand == '-' && b.strand == '-') {
                        if (b.te > a.ts + TOL) continue;
                        if (d.has_flo[i] && b.ts < d.flo[i] - TOL) continue;
                        same = true;
                    } else if (b.ts < d.shi[i] - TOL) {
                        continue;
                    }
                } else {
                    if (a.strand == '-' && b.strand == '-') {
                        if (b.te > a.ts + TOL) continue;
                        same = true;
                    } else if (a.strand == '+' && b.strand == '+') {
                        if (b.ts < a.te - TOL) continue;
                        if (d.has_flo[i] && b.te > d.flo[i] + TOL) continue;
                        same = true;
                    } else if (b.te > d.slo[i] + TOL) {
                        continue;
                    }
                }
                int64_t s2 = d.sc[i] + b.n_eq;
                if (s2 > d.sc[j]) {
                    d.sc[j] = s2;
                    d.prv[j] = (int32_t)i;
                    d.ok[j] = 1;
                    d.hlo[j] = std::min(d.hlo[i], b.ts);
                    if (same) {
                        d.has_flo[j] = d.has_flo[i];
                        d.flo[j] = d.flo[i];
                        d.slo[j] = std::min(d.slo[i], b.ts);
                        d.shi[j] = std::max(d.shi[i], b.te);
                    } else {
                        d.has_flo[j] = 1;
                        d.flo[j] = fr == 0 ? d.shi[i] : d.slo[i];
                        d.slo[j] = b.ts;
                        d.shi[j] = b.te;
                    }
                }
            }
        }
    }

    // groups split where consecutive elements are >= --island apart on both axes; groups with
    // < b aligned query are removed
    void island_rule(const std::vector<uint32_t>& el, Chain& ch) const {
        std::vector<std::vector<uint32_t>> groups;
        for (size_t k = 0; k < el.size(); ++k) {
            if (k == 0) { groups.push_back({el[0]}); continue; }
            const ChainPart& a = c_[el[k - 1]];
            const ChainPart& b = c_[el[k]];
            int64_t tgap = std::max<int64_t>(0, std::max(a.ts, b.ts) - std::min(a.te, b.te));
            if (b.qs - a.qe >= ISL_ && tgap >= ISL_) groups.push_back({el[k]});
            else groups.back().push_back(el[k]);
        }
        ch.el.clear();
        ch.islands.clear();
        ch.island_bp = 0;
        for (const std::vector<uint32_t>& gr : groups) {
            int64_t aq = 0;
            for (uint32_t e : gr) aq += c_[e].qe - c_[e].qs;
            if (aq >= B_) ch.el.insert(ch.el.end(), gr.begin(), gr.end());
            else {
                ch.islands.push_back(std::make_pair(c_[gr.front()].qs, c_[gr.back()].qe));
                ch.island_bp += aq;
            }
        }
        ch.score = ch.aq = 0;
        ch.hlo = INT64_MAX;
        ch.thi = INT64_MIN;
        ch.qlo = ch.el.empty() ? 0 : c_[ch.el.front()].qs;
        for (uint32_t e : ch.el) {
            ch.score += c_[e].n_eq;
            ch.aq += c_[e].qe - c_[e].qs;
            ch.hlo = std::min(ch.hlo, c_[e].ts);
            ch.thi = std::max(ch.thi, c_[e].te);
        }
        ch.valid = !ch.el.empty();
    }

    const std::vector<ChainPart>& c_;
    size_t n_;
    int64_t G_, B_, ISL_;
    DP base_[2], work_[2];
};

std::vector<ChainPart> chain_parts(const std::vector<ChainPart>& c, const Chain& ch) {
    std::vector<ChainPart> v;
    for (uint32_t e : ch.el) v.push_back(c[e]);
    return v;
}

} // namespace

// ================================================================ P

int64_t project_left_snap(const ChainPart& p, int64_t x) {
    int64_t c = p.strand == '+' ? x - p.qs : p.qe - x;   // query bases consumed
    int64_t t = p.ts, q = 0;
    if (c <= 0) return t;
    for (const CigarOp& op : p.ops) {
        if (op.op == '=' || op.op == 'X' || op.op == 'M') {
            if (c <= q + (int64_t)op.len) return t + (c - q);
            q += op.len;
            t += op.len;
        } else if (op.op == 'I') {
            if (c <= q + (int64_t)op.len) return t;
            q += op.len;
        } else if (op.op == 'D') {
            t += op.len;
        }
    }
    return t;
}

std::vector<int64_t> project_points(const ChainPart& p, const std::vector<int64_t>& xs) {
    const size_t n = xs.size();
    std::vector<int64_t> out(n, p.ts);
    std::vector<std::pair<int64_t, size_t>> cs(n);
    for (size_t i = 0; i < n; ++i) cs[i] = std::make_pair(p.strand == '+' ? xs[i] - p.qs : p.qe - xs[i], i);
    std::sort(cs.begin(), cs.end());
    size_t i = 0;
    while (i < n && cs[i].first <= 0) out[cs[i++].second] = p.ts;
    int64_t t = p.ts, q = 0;
    for (const CigarOp& op : p.ops) {
        if (i >= n) break;
        const int64_t L = op.len;
        if (op.op == '=' || op.op == 'X' || op.op == 'M') {
            while (i < n && cs[i].first <= q + L) { out[cs[i].second] = t + (cs[i].first - q); ++i; }
            q += L;
            t += L;
        } else if (op.op == 'I') {
            while (i < n && cs[i].first <= q + L) { out[cs[i].second] = t; ++i; }
            q += L;
        } else if (op.op == 'D') {
            t += L;
        }
    }
    while (i < n) out[cs[i++].second] = t;
    return out;
}

// ================================================================ chaining API

QueryLayout QueryLayout::of(const Graph& g, const std::vector<Handle>& walk) {
    QueryLayout ql;
    ql.walk = walk;
    ql.off.resize(walk.size() + 1);
    int64_t o = 0;
    for (size_t k = 0; k < walk.size(); ++k) {
        ql.off[k] = o;
        o += g.len(handle_node(walk[k]));
    }
    ql.off[walk.size()] = o;
    return ql;
}

std::vector<ChainPart> record_blocks(const PafRecord& r, uint32_t record_index, int64_t G) {
    std::vector<ChainPart> out;
    const bool fwd = r.strand == '+';
    int64_t q = fwd ? r.qs : r.qe, t = r.ts, q0 = 0;
    bool open = false, aligned = false;
    ChainPart cur;
    uint32_t bi = 0;
    auto close = [&]() {
        if (!open) return;
        if (fwd) { cur.qs = q0; cur.qe = q; }
        else { cur.qs = q; cur.qe = q0; }
        cur.te = t;
        if (aligned) {
            cur.block = bi++;
            out.push_back(cur);
        }
        open = false;
    };
    for (const CigarOp& o : r.ops) {
        const int64_t L = o.len;
        if ((o.op == 'I' || o.op == 'D') && L >= G) {
            close();
            if (o.op == 'I') q += fwd ? L : -L;
            else t += L;
            continue;
        }
        if (!open) {
            cur = ChainPart();
            cur.record = record_index;
            cur.strand = r.strand;
            cur.ts = t;
            q0 = q;
            open = true;
            aligned = false;
        }
        push_op(cur.ops, o.op, L);
        if (o.op == '=' || o.op == 'X') {
            q += fwd ? L : -L;
            t += L;
            aligned = true;
            if (o.op == '=') cur.n_eq += L;
        } else if (o.op == 'I') {
            q += fwd ? L : -L;
        } else {
            t += L;
        }
    }
    close();
    return out;
}

int64_t chain_E(const std::vector<ChainPart>& parts, int64_t qlen, int64_t tlo, int64_t thi, int64_t G) {
    int64_t e = 0;
    std::vector<std::pair<int64_t, int64_t>> iv;
    for (int axis = 0; axis < 2; ++axis) {
        iv.clear();
        for (const ChainPart& p : parts) iv.push_back(axis == 0 ? std::make_pair(p.qs, p.qe) : std::make_pair(p.ts, p.te));
        std::sort(iv.begin(), iv.end());
        int64_t pos = axis == 0 ? 0 : tlo, end = axis == 0 ? qlen : thi;
        for (auto& x : iv) {
            if (x.first - pos > G) ++e;
            pos = std::max(pos, x.second);
        }
        if (end - pos > G) ++e;
    }
    std::string s = strand_string(parts);
    return e + (s.empty() ? 0 : (int64_t)s.size() - 1);
}

void make_pieces(const QueryLayout& ql, ChainResult& ch) {
    ch.pieces.clear();
    std::vector<std::array<int64_t, 3>> hs;
    std::vector<int64_t> xs;
    for (uint32_t pi = 0; pi < (uint32_t)ch.parts.size(); ++pi) {
        const ChainPart& p = ch.parts[pi];
        overlapping_handles(ql, p.qs, p.qe, hs);
        xs.clear();
        for (auto& h : hs) { xs.push_back(h[1]); xs.push_back(h[2]); }
        std::vector<int64_t> P = project_points(p, xs);
        const bool fwd = p.strand == '+';
        for (size_t i = 0; i < hs.size(); ++i) {
            size_t k = (size_t)hs[i][0];
            Handle h = ql.walk[k];
            int64_t len = ql.node_len(k), a0 = hs[i][1] - ql.off[k], a1 = hs[i][2] - ql.off[k];
            Piece pc;
            pc.part = pi;
            pc.walk_pos = (uint32_t)k;
            pc.node = handle_node(h);
            if (!handle_rev(h)) { pc.a = a0; pc.b = a1; }
            else { pc.a = len - a1; pc.b = len - a0; }
            pc.qa = hs[i][1];
            pc.qb = hs[i][2];
            pc.ta = fwd ? P[2 * i] : P[2 * i + 1];
            pc.tb = fwd ? P[2 * i + 1] : P[2 * i];
            pc.rel = (!handle_rev(h)) == fwd ? '+' : '-';
            ch.pieces.push_back(pc);
        }
    }
}

void chain_records(const QueryLayout& ql, int64_t tlo, int64_t thi, const FeasibleFn& feasible, const Options& opt, AlignResult& out) {
    out.chain = ChainResult();
    const int64_t G = opt.G;
    // records -> blocks -> feasible sub-records (each >= --min-piece)
    std::vector<ChainPart> parts;
    for (uint32_t ri = 0; ri < (uint32_t)out.records.size(); ++ri) {
        const PafRecord& r = out.records[ri];
        if (r.qe - r.qs < opt.min_piece || r.gci < opt.ident) continue;
        for (const ChainPart& blk : record_blocks(r, ri, G)) feasible_parts(blk, ql, feasible, opt.min_piece, parts);
    }
    if (parts.empty()) {
        out.outcome = "no-chain";
        return;
    }
    // identical sub-records from different records (same coordinates and ops) are one alignment, not
    // two placements: keep the one from the lowest (record, block), or rule U would see a spurious tie
    std::sort(parts.begin(), parts.end(), part_less);
    {
        std::vector<ChainPart> uniq;
        size_t run0 = 0;
        for (size_t i = 0; i < parts.size(); ++i) {
            const ChainPart& p = parts[i];
            if (!uniq.empty()) {
                const ChainPart& q = uniq[run0];
                if (!(p.qs == q.qs && p.qe == q.qe && p.ts == q.ts && p.te == q.te && p.strand == q.strand)) run0 = uniq.size();
            }
            bool dup = false;
            for (size_t k = run0; k < uniq.size() && !dup; ++k) {
                const ChainPart& q = uniq[k];
                if (q.qs == p.qs && q.qe == p.qe && q.ts == p.ts && q.te == p.te && q.strand == p.strand && q.ops.size() == p.ops.size()) {
                    dup = true;
                    for (size_t o = 0; o < q.ops.size() && dup; ++o)
                        if (q.ops[o].op != p.ops[o].op || q.ops[o].len != p.ops[o].len) dup = false;
                }
            }
            if (!dup) uniq.push_back(p);
        }
        parts.swap(uniq);
    }
    // the 500 highest-scoring sub-records
    out.chain.n_feasible = (uint32_t)parts.size();
    if (parts.size() > MAX_CHAIN_PARTS) {
        std::sort(parts.begin(), parts.end(), [](const ChainPart& a, const ChainPart& b) {
            if (a.n_eq != b.n_eq) return a.n_eq > b.n_eq;
            return part_less(a, b);
        });
        parts.resize(MAX_CHAIN_PARTS);
        out.chain.truncated = true;
    }
    std::sort(parts.begin(), parts.end(), part_less);
    out.chain.n_candidates = (uint32_t)parts.size();

    Chainer cz(parts, opt);
    Chain star = cz.best(-1);
    if (!star.valid) {
        out.outcome = "below-b";
        return;
    }
    // ---- rule U
    const double keep = 1.0 - opt.delta;
    std::vector<uint8_t> in_star(parts.size(), 0);
    for (uint32_t e : star.el) in_star[e] = 1;
    std::vector<int64_t> bound = cz.bounds();
    std::set<std::vector<uint32_t>> seen;
    std::vector<Chain> evaluated, tied;
    for (size_t x = 0; x < parts.size(); ++x) {
        if (in_star[x]) continue;
        if ((double)bound[x] + 1.0 < keep * (double)star.score) continue;   // cannot tie: tie => score >= (1-delta) S*
        Chain cx = cz.best((int)x);
        if (!cx.valid) continue;
        if (!seen.insert(cx.el).second) continue;
        bool subset = true;
        for (uint32_t e : cx.el)
            if (!in_star[e]) { subset = false; break; }
        if (subset) continue;
        std::vector<uint8_t> in_x(parts.size(), 0);
        for (uint32_t e : cx.el) in_x[e] = 1;
        int64_t own = 0, alt = 0;
        for (uint32_t e : star.el)
            if (!in_x[e]) own += parts[e].n_eq;
        for (uint32_t e : cx.el)
            if (!in_star[e]) alt += parts[e].n_eq;
        evaluated.push_back(cx);
        if ((double)alt >= keep * (double)own) tied.push_back(cx);
    }
    const int64_t qlen = ql.length();
    std::vector<Chain> cands;
    cands.push_back(star);
    for (Chain& c : tied) cands.push_back(c);
    std::vector<int64_t> E(cands.size());
    for (size_t i = 0; i < cands.size(); ++i) E[i] = chain_E(chain_parts(parts, cands[i]), qlen, tlo, thi, G);
    size_t pick = 0;
    std::string status = "unique";
    bool ambiguous = false;
    if (!tied.empty()) {
        int64_t m = *std::min_element(E.begin(), E.end());
        std::vector<size_t> best;
        std::set<std::string> strands;
        for (size_t i = 0; i < cands.size(); ++i)
            if (E[i] == m) {
                best.push_back(i);
                strands.insert(strand_string(chain_parts(parts, cands[i])));
            }
        if (strands.size() > 1) {
            ambiguous = true;
            status = "ambiguous-strand";
            pick = 0;
        } else {
            // phase ties: the survivors lie in one stretch of the window that this unit's
            // sub-records cover with no gap > G
            int64_t lo = INT64_MAX, hi = INT64_MIN;
            for (size_t i : best) { lo = std::min(lo, cands[i].hlo); hi = std::max(hi, cands[i].thi); }
            std::vector<std::pair<int64_t, int64_t>> iv;
            for (const ChainPart& p : parts)
                if (p.te > lo && p.ts < hi) iv.push_back(std::make_pair(p.ts, p.te));
            std::sort(iv.begin(), iv.end());
            int64_t pos = lo;
            bool phase = true;
            for (auto& x : iv) {
                if (x.first > pos + G) { phase = false; break; }
                pos = std::max(pos, x.second);
            }
            phase = phase && pos >= hi - G;
            pick = best[0];
            for (size_t k = 1; k < best.size(); ++k) {
                const Chain& a = cands[best[k]];
                const Chain& b = cands[pick];
                bool better;
                if (phase) better = a.hlo != b.hlo ? a.hlo < b.hlo : (a.qlo != b.qlo ? a.qlo < b.qlo : a.score > b.score);
                else better = a.score != b.score ? a.score > b.score : (a.hlo != b.hlo ? a.hlo < b.hlo : a.qlo < b.qlo);
                if (better) pick = best[k];
            }
            status = pick == 0 ? "tie-kept" : "tie-resolved";
        }
    }
    const Chain& P = cands[pick];
    if (!opt.dump_dir.empty()) {
        // role score E tied strands frame tlo thi qlo nparts (C* first, then every evaluated C_x)
        auto line = [&](const char* role, const Chain& c, bool tie) {
            std::vector<ChainPart> cp = chain_parts(parts, c);
            out.chain.u_trace += strf("%s\t%lld\t%lld\t%d\t%s\t%c\t%lld\t%lld\t%lld\t%zu%s\n", role, (long long)c.score,
                                      (long long)chain_E(cp, qlen, tlo, thi, G), tie ? 1 : 0, strand_string(cp).c_str(), c.frame == 0 ? 'F' : 'R',
                                      (long long)c.hlo, (long long)c.thi, (long long)c.qlo, cp.size(), c.el == P.el ? "\tpicked" : "");
        };
        line("star", star, false);
        for (const Chain& c : evaluated) {
            bool t = false;
            for (const Chain& x : tied) if (x.el == c.el) t = true;
            line("alt", c, t);
        }
    }
    // the best alternative: highest score among C* and every evaluated C_x other than the pick
    {
        const Chain* alt = nullptr;
        int64_t altE = -1;
        std::vector<const Chain*> pool;
        if (pick != 0) pool.push_back(&star);
        for (const Chain& c : evaluated)
            if (c.el != P.el) pool.push_back(&c);
        for (const Chain* c : pool) {
            int64_t e = chain_E(chain_parts(parts, *c), qlen, tlo, thi, G);
            if (!alt || c->score > alt->score || (c->score == alt->score && e < altE)) { alt = c; altE = e; }
        }
        if (alt) { out.chain.alt_score = alt->score; out.chain.alt_E = altE; }
    }
    ChainResult& ch = out.chain;
    ch.parts = chain_parts(parts, P);
    ch.frame = P.frame == 0 ? 'F' : 'R';
    ch.strands = strand_string(ch.parts);
    ch.score = P.score;
    ch.aligned_q = P.aq;
    ch.tie = status;
    ch.E = E[pick];
    ch.islands = P.islands;
    ch.island_bp = P.island_bp;
    ch.island_dropped = !P.islands.empty();
    ch.internal = internal_stretches(ch.parts, G);
    make_pieces(ql, ch);
    if (ambiguous) { out.outcome = "ambiguous-strand"; return; }
    if (opt.ties_refuse && status != "unique") { out.outcome = "tie-refused"; return; }
    int64_t limit = std::max<int64_t>(2, opt.frag * ch.aligned_q / 5000);
    if (ch.internal > limit) { out.outcome = "fragmented"; return; }
    out.outcome = "";
}

// ================================================================ the aligner

namespace {

// fraction of the window's non-N bases covered by a canonical 16-mer occurring at least twice
double repeat_fraction(const std::string& s) {
    const int k = REPEAT_K;
    if ((int64_t)s.size() < k) return 0.0;
    std::vector<std::pair<uint32_t, uint32_t>> km;
    km.reserve(s.size());
    uint32_t fw = 0, rv = 0;
    const uint32_t mask = 0xffffffffu;
    int valid = 0;
    int64_t nonN = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t c = nt4(s[i]);
        if (c > 3) { valid = 0; continue; }
        ++nonN;
        fw = ((fw << 2) | c) & mask;
        rv = (rv >> 2) | ((uint32_t)(3 - c) << (2 * (k - 1)));
        if (++valid >= k) km.push_back(std::make_pair(std::min(fw, rv), (uint32_t)(i + 1 - k)));
    }
    if (nonN == 0) return 0.0;
    std::sort(km.begin(), km.end());
    std::vector<int32_t> cover(s.size() + 1, 0);
    for (size_t i = 0; i < km.size();) {
        size_t j = i;
        while (j < km.size() && km[j].first == km[i].first) ++j;
        if (j - i >= 2)
            for (size_t x = i; x < j; ++x) { cover[km[x].second] += 1; cover[km[x].second + k] -= 1; }
        i = j;
    }
    int64_t cov = 0, run = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        run += cover[i];
        if (run > 0 && nt4(s[i]) < 4) ++cov;
    }
    return (double)cov / (double)nonN;
}

// one unit or job to align
struct Item {
    std::string key, alt_key;          // injection keys
    std::vector<Handle> walk;          // query walk
    AlignTarget target;
    FeasibleFn feasible;
    AlignResult* res = nullptr;
};

struct WinQuery {
    uint32_t item = 0;
    std::string seq;
};

struct WindowRun {
    bool failed = false, signaled = false;
    std::string err;
    std::vector<uint8_t> screened, prefiltered, aligned, bad_screen;
    std::vector<std::vector<RawRec>> recs;   // pass-2 records per query
    long peak_rss_kb = 0;
    AlignTiming timing;                      // its processes: running, waiting for a slot, minimap2's own times
};

// one minimap2 process on the steady clock: asked for a budget slot, got it, ended; minimap2's own times
struct ProcClock {
    double req = 0, acq = 0, end = 0;
    double real = 0, cpu = 0;
};

struct InjectLine {
    std::string type;
    int64_t qs = 0, qe = 0, ts = 0, te = 0;
    char strand = '+';
    std::string cigar;
    uint64_t lineno = 0;
};

} // namespace

struct Aligner::Impl {
    Options opt;
    std::string version;
    std::unique_ptr<TempDir> tmp;
    std::unique_ptr<ProcessBudget> budget;
    std::vector<std::string> extra;
    bool inject = false;
    std::map<std::string, std::vector<InjectLine>> injected;

    std::atomic<uint64_t> n_windows{0}, n_failed{0}, n_inconsistent{0}, n_processes{0}, n_signals{0}, n_retries{0};
    std::atomic<uint64_t> n_screened{0}, n_prefiltered{0}, n_pass2{0}, n_truncated{0}, n_repeat_skipped{0}, n_units{0}, n_confident{0};
    std::atomic<uint64_t> win_counter{0};
    std::atomic<long> peak_rss{0};
    std::atomic<int64_t> mm2_millis{0};
    std::atomic<int64_t> pass_millis[3] = {{0}, {0}, {0}};   // wall ms in minimap2 by pass: (unused), screen, pass 2
    // minimap2's own real and CPU ms (its last stderr line), and ms spent waiting for a budget slot,
    // summed over the processes
    std::atomic<int64_t> real_millis{0}, cpu_millis{0}, wait_millis{0};

    void note_rss(long kb) {
        long cur = peak_rss.load();
        while (kb > cur && !peak_rss.compare_exchange_weak(cur, kb)) {}
    }

    std::vector<std::string> mm2_args(int pass) const {
        // pass 0: index; 1: screen; 2: pass 2
        std::vector<std::string> a = {opt.minimap2};
        if (pass == 2) { a.push_back("-c"); a.push_back("--eqx"); }
        a.push_back("-x");
        a.push_back(opt.preset);
        if (pass != 0) {
            for (const char* s : {"-N", "50", "-p", "0.01", "--secondary=yes", "-t", "1"}) a.push_back(s);
        }
        a.insert(a.end(), extra.begin(), extra.end());
        return a;
    }

    // run one minimap2 process (holding a budget slot unless the caller holds the whole budget).
    // Its peak RSS is minimap2's own figure (wait4's includes rgfa-zip's memory; see ProcResult).
    // pc: when it asked for a slot, got it and ended, and minimap2's own real and CPU time.
    ProcResult run_mm2(const std::vector<std::string>& argv, const std::string& out_path, bool exclusive, int pass, ProcClock& pc) {
        ProcResult r;
        pc = ProcClock();
        pc.req = steady_seconds();
        if (exclusive) {
            pc.acq = pc.req;
            r = run_process(argv, out_path);
        } else {
            budget->acquire();
            pc.acq = steady_seconds();   // time the process, not the wait for a slot
            try {
                r = run_process(argv, out_path);
            } catch (...) {
                budget->release();
                throw;
            }
            budget->release();
        }
        pc.end = steady_seconds();
        ++n_processes;
        r.max_rss_kb = mm2_peak_rss_kb(r.err_tail);
        note_rss(r.max_rss_kb);
        mm2_times(r.err_tail, pc.real, pc.cpu);
        int64_t ms = (int64_t)((pc.end - pc.acq) * 1000.0);
        mm2_millis += ms;
        pass_millis[pass] += ms;
        wait_millis += (int64_t)((pc.acq - pc.req) * 1000.0);
        real_millis += (int64_t)(pc.real * 1000.0);
        cpu_millis += (int64_t)(pc.cpu * 1000.0);
        if (r.signaled) ++n_signals;
        return r;
    }

    // One phase over a subset of the window's queries, split over up to -j processes, each mapping
    // against the window's FASTA (minimap2 indexes it in memory: no index file can be written short).
    // Results go into wr (per query).  A signal is recorded in wr (the window is retried alone),
    // and so is any other non-zero exit or unparseable output (aligner-failed), except an I/O
    // failure of minimap2's own files in --tmpdir (ZipError EXIT_IO).  A process that cannot be
    // started or whose status cannot be read is a systemic failure (ZipError EXIT_ALIGNER).
    void run_phase(int pass, const std::string& dir, const std::string& tfa, const std::vector<WinQuery>& qs,
                   const std::vector<uint32_t>& idx, bool exclusive, WindowRun& wr, std::vector<std::vector<RawRec>>& recs) {
        if (idx.empty()) return;
        size_t k = std::min<size_t>((size_t)budget->total(), idx.size());
        // largest first onto the lightest batch (ties: lower batch); queries keep canonical order in a batch
        std::vector<uint32_t> order(idx);
        std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return qs[a].seq.size() > qs[b].seq.size(); });
        std::vector<std::vector<uint32_t>> batch(k);
        std::vector<int64_t> load(k, 0);
        for (uint32_t q : order) {
            size_t bi = 0;
            for (size_t b = 1; b < k; ++b)
                if (load[b] < load[bi]) bi = b;
            batch[bi].push_back(q);
            load[bi] += (int64_t)qs[q].seq.size();
        }
        for (auto& b : batch) std::sort(b.begin(), b.end());
        struct BatchOut {
            ProcResult pr;
            ProcClock pc;
            bool ran = false;            // run_mm2 returned (pc is set)
            bool parsed = false, unreadable = false;
            std::string why;
            std::vector<RawRec> recs;
            std::exception_ptr ex;
        };
        std::vector<BatchOut> outs(k);
        auto run_batch = [&](size_t b) {
            try {
                // A query's name is its index in the window's query list, never its position in
                // the batch: minimap2 seeds its tie-breaking hash with the query name (map.c,
                // mm_map_frag), so a batch-dependent name made records depend on -j.
                std::string fa, base = strf("%s/p%d.b%zu", dir.c_str(), pass, b);
                std::vector<uint8_t> in_batch(qs.size(), 0);
                for (size_t i = 0; i < batch[b].size(); ++i) {
                    fa += strf(">q%u\n", batch[b][i]);
                    fa += qs[batch[b][i]].seq;
                    fa.push_back('\n');
                    in_batch[batch[b][i]] = 1;
                }
                write_file(base + ".fa", fa);
                std::vector<std::string> argv = mm2_args(pass);
                argv.push_back(tfa);
                argv.push_back(base + ".fa");
                outs[b].pr = run_mm2(argv, base + ".paf", exclusive, pass, outs[b].pc);
                outs[b].ran = true;
                if (outs[b].pr.ok()) {
                    std::string text;
                    if (!read_file(base + ".paf", text)) {
                        outs[b].unreadable = true;
                        outs[b].why = strf("cannot read minimap2's output %s.paf: %s", base.c_str(), strerror(errno));
                    } else {
                        outs[b].parsed = parse_paf(text, qs.size(), in_batch, pass == 2, outs[b].recs, outs[b].why);
                    }
                }
                unlink((base + ".fa").c_str());
                unlink((base + ".paf").c_str());
            } catch (...) {
                outs[b].ex = std::current_exception();
            }
        };
        // batch 0 runs here; the others on threads, or here in turn when no thread can start
        size_t started = 1;
        std::vector<std::thread> th;
        th.reserve(k);
        for (size_t b = 1; b < k; ++b) {
            try {
                th.emplace_back(run_batch, b);
                started = b + 1;
            } catch (const std::exception& e) {
                ZLOG("warning: could not start an aligner thread (%s); running the window's remaining batches one by one", e.what());
                break;
            }
        }
        run_batch(0);
        for (size_t b = started; b < k; ++b) run_batch(b);
        for (std::thread& t : th) t.join();
        // the first real error wins over Aborted (which only says that someone else failed)
        bool aborted = false;
        for (size_t b = 0; b < k; ++b) {
            if (!outs[b].ex) continue;
            try {
                std::rethrow_exception(outs[b].ex);
            } catch (const Aborted&) {
                aborted = true;
            }
        }
        if (aborted || aborting()) throw Aborted();
        // systemic failures first: the aligner cannot be run, or its files cannot be written or read
        for (size_t b = 0; b < k; ++b) {
            const BatchOut& o = outs[b];
            if (!o.pr.spawned)
                fail(EXIT_ALIGNER, strf("cannot start minimap2 (%s): too many processes or too little memory? Lower -t/-j or raise "
                                        "the limits", o.pr.err_tail.c_str()));
            if (!o.pr.exited && !o.pr.signaled) fail(EXIT_ALIGNER, "cannot read the exit status of a minimap2 process");
            if (o.pr.exited && o.pr.exit_code != 0 && mm2_io_error(o.pr.err_tail))
                fail(EXIT_IO, strf("minimap2 could not write or read its files in --tmpdir %s (%s): %s", opt.tmpdir.c_str(),
                                   o.pr.describe().c_str(), mm2_error_lines(o.pr.err_tail).c_str()));
            if (o.unreadable) fail(EXIT_IO, o.why);
        }
        for (size_t b = 0; b < k; ++b) {
            BatchOut& o = outs[b];
            if (o.ran) {
                wr.timing.run.push_back(std::make_pair(o.pc.acq, o.pc.end));
                if (o.pc.acq > o.pc.req) wr.timing.wait.push_back(std::make_pair(o.pc.req, o.pc.acq));
                wr.timing.mm2_real += o.pc.real;
                wr.timing.mm2_cpu += o.pc.cpu;
                ++wr.timing.processes;
            }
            wr.peak_rss_kb = std::max(wr.peak_rss_kb, o.pr.max_rss_kb);
            if (o.pr.signaled) {
                wr.signaled = true;
                wr.err = strf("minimap2 %s", o.pr.describe().c_str());
            } else if (!o.pr.ok()) {
                wr.failed = true;
                if (wr.err.empty()) wr.err = strf("minimap2 %s: %s", o.pr.describe().c_str(), o.pr.err_tail.c_str());
            } else if (!o.parsed) {
                wr.failed = true;
                if (wr.err.empty()) wr.err = "unparseable minimap2 output: " + o.why;
            } else {
                for (RawRec& r : o.recs) recs[r.q].push_back(std::move(r));
            }
        }
    }

    // screen and pass 2 for one window.  There is no index file: every process indexes the window's
    // FASTA in memory (0.2 s for a 5 Mb window; output identical to a prebuilt .mmi), because
    // minimap2 -d does not check its writes -- a full --tmpdir left an empty or partial index and
    // the window's queries silently 'prefiltered'.
    WindowRun run_window(const std::string& tseq, const std::vector<WinQuery>& qs, bool exclusive) {
        WindowRun wr;
        const size_t n = qs.size();
        wr.screened.assign(n, 0);
        wr.prefiltered.assign(n, 0);
        wr.aligned.assign(n, 0);
        wr.bad_screen.assign(n, 0);
        wr.recs.assign(n, std::vector<RawRec>());
        std::string dir = strf("%s/w%llu", tmp->path().c_str(), (unsigned long long)++win_counter);
        if (mkdir(dir.c_str(), 0755) != 0) fail(EXIT_IO, strf("mkdir %s: %s", dir.c_str(), strerror(errno)));
        struct DirGuard {
            std::string d;
            ~DirGuard() { remove_tree(d); }
        } guard{dir};
        write_file(dir + "/t.fa", ">t\n" + tseq + "\n");
        const std::string tfa = dir + "/t.fa";
        // screen
        std::vector<std::vector<RawRec>> srec(n);
        auto passes = [&](uint32_t q) {
            int64_t sum = 0;
            double cov = 0;
            for (const RawRec& r : srec[q]) {
                if (r.qlen != (int64_t)qs[q].seq.size() || r.tlen != (int64_t)tseq.size() || !(0 <= r.qs && r.qs < r.qe && r.qe <= r.qlen) ||
                    !(0 <= r.ts && r.ts < r.te && r.te <= r.tlen) || r.nmatch < 0 || r.blen <= 0 || r.nmatch > r.blen)
                    wr.bad_screen[q] = 1;
                sum += r.nmatch;
                cov = std::max(cov, (double)r.nmatch / (double)std::max<int64_t>(1, r.blen));
            }
            return sum >= opt.prefilter_nmatch && cov >= opt.prefilter_cov;
        };
        std::vector<uint32_t> sample, rest, to_align;
        size_t ns = std::min<size_t>(n, (size_t)std::max(1, opt.screen_sample));
        for (uint32_t q = 0; q < n; ++q) (q < ns ? sample : rest).push_back(q);
        run_phase(1, dir, tfa, qs, sample, exclusive, wr, srec);
        if (wr.failed || wr.signaled) return wr;
        for (uint32_t q : sample) wr.screened[q] = 1;
        size_t npass = 0;
        for (uint32_t q : sample) npass += passes(q) ? 1 : 0;
        if (2 * npass >= sample.size()) {
            for (uint32_t q = 0; q < n; ++q) to_align.push_back(q);
        } else {
            run_phase(1, dir, tfa, qs, rest, exclusive, wr, srec);
            if (wr.failed || wr.signaled) return wr;
            for (uint32_t q : rest) wr.screened[q] = 1;
            for (uint32_t q = 0; q < n; ++q) {
                if (passes(q)) to_align.push_back(q);
                else wr.prefiltered[q] = 1;
            }
        }
        for (uint32_t q : to_align) wr.aligned[q] = 1;
        run_phase(2, dir, tfa, qs, to_align, exclusive, wr, wr.recs);
        return wr;
    }

    // ------------------------------------------------------------ verification and chaining of one item
    void finish_item(const Graph& g, Item& it, const std::string& qseq, const std::string& tseq, int64_t tlo, int64_t thi,
                     const std::vector<RawRec>& raw, bool check_counts, bool allow_m) {
        AlignResult& r = *it.res;
        r.records.clear();
        for (const RawRec& rr : raw) {
            PafRecord pr;
            std::string why;
            if (!verify_record(rr, qseq, tseq, tlo, allow_m, check_counts, opt.G, pr, why)) {
                ++n_inconsistent;
                ZLOG("paf-inconsistent: %s: %s", it.key.c_str(), why.c_str());
                r.outcome = "paf-inconsistent";
                r.records.clear();
                return;
            }
            r.records.push_back(std::move(pr));
        }
        canonical_records(r.records);
        QueryLayout ql = QueryLayout::of(g, it.walk);
        chain_records(ql, tlo, thi, it.feasible, opt, r);
        if (r.chain.truncated) ++n_truncated;
        if (r.outcome.empty()) ++n_confident;
        // after chaining, only the count (the report) and --dump need the records: drop their op
        // vectors (at chr1:2.65, 41,549 records and 99 M ops, 0.9 of rgfa-zip's 1.5 GB)
        r.n_records = (int64_t)r.records.size();
        if (opt.dump_dir.empty()) std::vector<PafRecord>().swap(r.records);
    }

    // ------------------------------------------------------------ injection
    void inject_item(const Graph& g, Item& it, const std::string& tseq, int64_t tlo, int64_t thi) {
        AlignResult& r = *it.res;
        r.attempted = true;
        auto f = injected.find(it.key);
        if (f == injected.end() && !it.alt_key.empty()) f = injected.find(it.alt_key);
        if (f == injected.end()) { r.outcome = "no-chain"; return; }
        const std::vector<InjectLine>& lines = f->second;
        std::string qseq = g.spell(it.walk);
        bool chain_mode = lines[0].type == "chain";
        r.records.clear();
        for (const InjectLine& L : lines) {
            if ((L.type == "chain") != chain_mode)
                fail(EXIT_INPUT, strf("--inject-chains line %llu: %s mixes 'record' and 'chain' lines", (unsigned long long)L.lineno, it.key.c_str()));
            RawRec rr;
            rr.qlen = (int64_t)qseq.size();
            rr.tlen = (int64_t)tseq.size();
            rr.qs = L.qs;
            rr.qe = L.qe;
            rr.ts = L.ts - tlo;
            rr.te = L.te - tlo;
            rr.strand = L.strand;
            rr.cigar = L.cigar;
            rr.has_cigar = true;
            PafRecord pr;
            std::string why;
            if (!verify_record(rr, qseq, tseq, tlo, true, false, opt.G, pr, why))
                fail(EXIT_INPUT, strf("--inject-chains line %llu (%s) does not fit the sequences: %s", (unsigned long long)L.lineno, it.key.c_str(), why.c_str()));
            r.records.push_back(std::move(pr));
        }
        canonical_records(r.records);
        r.n_records = (int64_t)r.records.size();
        QueryLayout ql = QueryLayout::of(g, it.walk);
        if (!chain_mode) {
            chain_records(ql, tlo, thi, it.feasible, opt, r);
            if (r.chain.truncated) ++n_truncated;   // as finish_item
            if (r.outcome.empty()) ++n_confident;
            if (opt.dump_dir.empty()) std::vector<PafRecord>().swap(r.records);   // as finish_item
            return;
        }
        ChainResult& ch = r.chain;
        ch = ChainResult();
        for (uint32_t ri = 0; ri < (uint32_t)r.records.size(); ++ri)
            for (const ChainPart& blk : record_blocks(r.records[ri], ri, opt.G)) ch.parts.push_back(slice_part(blk, blk.qs, blk.qe));
        std::sort(ch.parts.begin(), ch.parts.end(), part_less);
        if (ch.parts.empty()) { r.outcome = "no-chain"; return; }
        ch.frame = (ch.parts.size() > 1 && ch.parts.back().ts < ch.parts.front().ts) ? 'R' : 'F';
        ch.strands = strand_string(ch.parts);
        for (const ChainPart& p : ch.parts) { ch.score += p.n_eq; ch.aligned_q += p.qe - p.qs; }
        ch.tie = "injected";
        ch.internal = internal_stretches(ch.parts, opt.G);
        ch.E = chain_E(ch.parts, ql.length(), tlo, thi, opt.G);
        ch.n_candidates = (uint32_t)ch.parts.size();
        make_pieces(ql, ch);
        r.outcome = "";
        ++n_confident;
        if (opt.dump_dir.empty()) std::vector<PafRecord>().swap(r.records);
    }

    void read_inject(const std::string& path) {
        LineReader lr(path);
        std::string line;
        std::vector<std::string> f;
        while (lr.next(line)) {
            if (line.empty() || line[0] == '#') continue;
            split_tabs(line, f);
            uint64_t ln = lr.line_number();
            if (f.size() < 8) fail(EXIT_INPUT, strf("--inject-chains %s line %llu: expected 8 columns, got %zu", path.c_str(), (unsigned long long)ln, f.size()));
            InjectLine L;
            L.lineno = ln;
            L.type = f[1];
            if (L.type != "record" && L.type != "chain")
                fail(EXIT_INPUT, strf("--inject-chains line %llu: type must be 'record' or 'chain'", (unsigned long long)ln));
            if (!parse_int64(f[2], L.qs) || !parse_int64(f[3], L.qe) || !parse_int64(f[5], L.ts) || !parse_int64(f[6], L.te))
                fail(EXIT_INPUT, strf("--inject-chains line %llu: bad coordinates", (unsigned long long)ln));
            if (f[4] != "+" && f[4] != "-") fail(EXIT_INPUT, strf("--inject-chains line %llu: strand must be + or -", (unsigned long long)ln));
            L.strand = f[4][0];
            L.cigar = f[7];
            injected[f[0]].push_back(L);
        }
        lr.close();
        inject = true;
        ZLOG("--inject-chains: %zu key(s) from %s; minimap2 will not run", injected.size(), path.c_str());
    }

    // ------------------------------------------------------------ the engine: items -> windows -> results
    void run_items(const Graph& g, std::vector<Item>& items, const std::string& dump_label) {
        if (items.empty()) return;
        n_units += items.size();
        // group by target (one index per distinct window / target walk), in order of first item
        std::vector<std::vector<uint32_t>> wins;
        {
            std::map<std::string, size_t> wkey;
            for (uint32_t i = 0; i < (uint32_t)items.size(); ++i) {
                const AlignTarget& t = items[i].target;
                std::string k;
                if (t.ref) k = strf("r%d:%lld:%lld", t.sn, (long long)t.lo, (long long)t.hi);
                else {
                    k = "a";
                    for (Handle h : t.walk) k += strf(":%u", h);
                }
                auto it = wkey.find(k);
                if (it == wkey.end()) {
                    wkey[k] = wins.size();
                    wins.push_back({i});
                } else {
                    wins[it->second].push_back(i);
                }
            }
        }
        const bool dumping = !opt.dump_dir.empty();
        std::vector<std::string> dump_win(wins.size()), dump_rec(wins.size()), dump_chain(wins.size()), dump_u(wins.size());
        std::vector<std::pair<std::string, std::string>> dump_fa(wins.size());
        auto target_seq = [&](const AlignTarget& t, int64_t& tlo, int64_t& thi) {
            if (t.ref) {
                tlo = t.lo;
                thi = t.hi;
                return g.ref_seq(t.sn, t.lo, t.hi);
            }
            std::string s = g.spell(t.walk);
            tlo = 0;
            thi = (int64_t)s.size();
            return s;
        };
        // timing: every window's, added to the caller's AlignTimingScope (main's site, zip_alt's pass)
        AlignTiming* const sink = AlignTimingScope::current();
        AlignTiming calls;
        std::mutex calls_m;
        auto do_window = [&](size_t w) {
            const double t0 = steady_seconds();
            AlignTiming wt;
            const std::vector<uint32_t>& its = wins[w];
            const AlignTarget& T = items[its[0]].target;
            int64_t tlo = 0, thi = 0;
            std::string tseq = target_seq(T, tlo, thi);
            double rf = repeat_fraction(tseq);
            for (uint32_t i : its) items[i].res->repeat_frac = rf;
            if (inject) {
                for (uint32_t i : its) inject_item(g, items[i], tseq, tlo, thi);
            } else if (opt.max_repeat_frac >= 0 && rf > opt.max_repeat_frac) {
                for (uint32_t i : its) items[i].res->outcome = "repeat-skipped";
                n_repeat_skipped += its.size();
            } else {
                std::vector<WinQuery> qs(its.size());
                for (size_t k = 0; k < its.size(); ++k) {
                    qs[k].item = its[k];
                    qs[k].seq = g.spell(items[its[k]].walk);
                }
                ++n_windows;
                WindowRun wr = run_window(tseq, qs, false);
                wt.add(wr.timing);
                check_abort();     // a stop requested meanwhile also killed this window's processes: no retry
                if (wr.signaled) {
                    ZLOG("aligner: window %s (%zu queries) %s; retrying it alone", items[its[0]].key.c_str(), its.size(), wr.err.c_str());
                    ++n_retries;
                    const double r0 = steady_seconds();
                    budget->acquire_all();
                    const double r1 = steady_seconds();
                    wt.wait.push_back(std::make_pair(r0, r1));
                    wait_millis += (int64_t)((r1 - r0) * 1000.0);
                    try {
                        wr = run_window(tseq, qs, true);
                    } catch (...) {
                        budget->release_all();
                        throw;
                    }
                    budget->release_all();
                    wt.add(wr.timing);
                    check_abort();
                    if (wr.signaled || wr.failed)
                        fail(EXIT_ALIGNER, strf("minimap2 failed twice on window %s (%zu queries): %s", items[its[0]].key.c_str(), its.size(), wr.err.c_str()));
                }
                note_rss(wr.peak_rss_kb);
                if (wr.failed) {
                    ++n_failed;
                    ZLOG("aligner-failed: window %s (%zu queries): %s", items[its[0]].key.c_str(), its.size(), wr.err.c_str());
                    for (uint32_t i : its) {
                        AlignResult& r = *items[i].res;
                        r.attempted = true;
                        r.outcome = "aligner-failed";
                        r.aligner_err = wr.err;
                    }
                } else {
                    for (size_t k = 0; k < its.size(); ++k) {
                        Item& it = items[its[k]];
                        AlignResult& r = *it.res;
                        r.attempted = true;
                        r.screened = wr.screened[k] != 0;
                        if (r.screened) ++n_screened;
                        if (wr.bad_screen[k]) {
                            ++n_inconsistent;
                            ZLOG("paf-inconsistent: %s: a screen (pass 1) record disagrees with the query or window lengths", it.key.c_str());
                            r.outcome = "paf-inconsistent";
                            continue;
                        }
                        if (wr.prefiltered[k]) {
                            ++n_prefiltered;
                            r.outcome = "prefiltered";
                            continue;
                        }
                        ++n_pass2;
                        finish_item(g, it, qs[k].seq, tseq, tlo, thi, wr.recs[k], true, false);
                    }
                }
                // wall time, split into minimap2 running, waiting for a -j slot (none of the window's
                // processes running), and the rest (rgfa-zip's own work: queries, PAF, chaining);
                // minimap2's own real and CPU time are what the window costs
                const double secs = steady_seconds() - t0, run_s = wt.running(), wait_s = wt.waiting();
                const std::string line = strf("aligner: window %s: %zu queries, %.2f s (minimap2 running %.2f s, waiting for an aligner slot "
                                              "%.2f s); %llu minimap2 process(es): real %.2f s, CPU %.2f s, peak RSS %.3f GB (minimap2's own "
                                              "figures)",
                                              items[its[0]].key.c_str(), its.size(), secs, run_s, wait_s, (unsigned long long)wt.processes,
                                              wt.mm2_real, wt.mm2_cpu, wr.peak_rss_kb / 1048576.0);
                if (secs - wait_s > 60) ZLOG("%s", line.c_str());
                else ZDEBUG("%s", line.c_str());
            }
            {
                std::lock_guard<std::mutex> lk(calls_m);
                calls.add(wt);
            }
            if (dumping) {
                std::string& dw = dump_win[w];
                dw += strf("%s\tw%zu\t%s\t%lld\t%lld\t%zu\t%.4f\n", dump_label.c_str(), w, T.ref ? "ref" : "alt", (long long)tlo, (long long)thi,
                           its.size(), rf);
                dump_fa[w].first = ">" + dump_label + ".w" + std::to_string(w) + "\n" + tseq + "\n";
                for (uint32_t i : its) {
                    Item& it = items[i];
                    const AlignResult& r = *it.res;
                    dump_fa[w].second += ">" + it.key + "\n" + g.spell(it.walk) + "\n";
                    for (size_t k = 0; k < r.records.size(); ++k) {
                        const PafRecord& p = r.records[k];
                        dump_rec[w] += strf("%s\tw%zu\t%zu\t%lld\t%lld\t%c\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%.5f\t%s\n", it.key.c_str(), w, k,
                                            (long long)p.qs, (long long)p.qe, p.strand, (long long)p.ts, (long long)p.te, (long long)p.n_eq,
                                            (long long)p.n_x, (long long)p.n_ins, (long long)p.n_del, p.gci, cigar_str(p.ops).c_str());
                    }
                    const ChainResult& ch = r.chain;
                    std::string head = strf("%s\tw%zu\t%s\t%s\t%c\t%s\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%u\t%d", it.key.c_str(), w,
                                            r.outcome.empty() ? "confident" : r.outcome.c_str(), ch.tie.empty() ? "." : ch.tie.c_str(),
                                            ch.parts.empty() ? '.' : ch.frame, ch.strands.empty() ? "." : ch.strands.c_str(), (long long)ch.score,
                                            (long long)ch.aligned_q, (long long)ch.internal, (long long)ch.E, (long long)ch.alt_score,
                                            (long long)ch.alt_E, ch.n_candidates, ch.truncated ? 1 : 0);
                    if (!ch.u_trace.empty()) {
                        size_t p0 = 0;
                        while (p0 < ch.u_trace.size()) {
                            size_t e = ch.u_trace.find('\n', p0);
                            dump_u[w] += it.key + "\t" + ch.u_trace.substr(p0, e - p0 + 1);
                            p0 = e + 1;
                        }
                    }
                    if (ch.parts.empty()) dump_chain[w] += head + "\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
                    for (size_t k = 0; k < ch.parts.size(); ++k) {
                        const ChainPart& p = ch.parts[k];
                        dump_chain[w] += head + strf("\t%zu\t%u\t%u\t%lld\t%lld\t%lld\t%lld\t%c\t%lld\n", k, p.record, p.block, (long long)p.qs,
                                                     (long long)p.qe, (long long)p.ts, (long long)p.te, p.strand, (long long)p.n_eq);
                    }
                }
            }
        };
        // windows in parallel (up to -j at once; each window splits its own phases over -j processes)
        // A fatal error stops this site's other windows and (request_abort) every other site's work
        // and running process; Aborted only means that some other thread failed first.
        size_t nw = wins.size();
        int workers = (int)std::min<size_t>(nw, (size_t)std::max(1, opt.jobs));
        std::atomic<size_t> next(0);
        std::exception_ptr err;
        std::mutex em;
        std::atomic<bool> stop(false), saw_abort(false);
        auto worker = [&](int) {
            while (!stop.load()) {
                if (aborting()) {
                    saw_abort = true;
                    break;
                }
                size_t w = next.fetch_add(1);
                if (w >= nw) break;
                try {
                    do_window(w);
                } catch (const Aborted&) {
                    saw_abort = true;
                    stop = true;
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(em);
                        if (!err) err = std::current_exception();
                    }
                    stop = true;
                    Runtime::get().request_abort();
                }
            }
        };
        if (workers <= 1) worker(0);
        else run_threads(workers, worker);
        if (err) std::rethrow_exception(err);
        if (saw_abort) throw Aborted();
        if (sink) sink->add(calls);
        if (dumping) {
            std::string dir = opt.dump_dir + "/align";
            if (!make_dirs(dir)) fail(EXIT_IO, "cannot create " + dir);
            std::string base = dir + "/" + sanitize(dump_label);
            auto put = [&](const std::string& path, const std::string& header, const std::vector<std::string>& parts) {
                AtomicWriter w(path);
                w.write(header);
                for (const std::string& s : parts) w.write(s);
                w.commit();
            };
            put(base + ".windows.tsv", "label\twindow\ttarget\ttlo\tthi\tqueries\trepeat_frac\n", dump_win);
            put(base + ".records.tsv", "key\twindow\trecord\tqs\tqe\tstrand\tts\tte\tn_eq\tn_x\tn_ins\tn_del\tgci\tcigar\n", dump_rec);
            put(base + ".chains.tsv",
                "key\twindow\toutcome\ttie\tframe\tstrands\tscore\taligned_q\tinternal\tE\talt_score\talt_E\tcandidates\ttruncated\t"
                "part\trecord\tblock\tqs\tqe\tts\tte\tstrand\tn_eq\n",
                dump_chain);
            put(base + ".u.tsv", "key\trole\tscore\tE\ttied\tstrands\tframe\ttlo\tthi\tqlo\tparts\tpicked\n", dump_u);
            for (size_t w = 0; w < wins.size(); ++w) {
                put(strf("%s.w%zu.t.fa", base.c_str(), w), "", {dump_fa[w].first});
                put(strf("%s.w%zu.q.fa", base.c_str(), w), "", {dump_fa[w].second});
            }
        }
    }
};

// ================================================================ timing

void AlignTiming::add(const AlignTiming& o) {
    run.insert(run.end(), o.run.begin(), o.run.end());
    wait.insert(wait.end(), o.wait.begin(), o.wait.end());
    mm2_real += o.mm2_real;
    mm2_cpu += o.mm2_cpu;
    processes += o.processes;
}

double AlignTiming::running() const {
    double s = 0;
    for (const auto& x : iv_union(run)) s += x.second - x.first;
    return s;
}

double AlignTiming::waiting() const {
    const std::vector<std::pair<double, double>> w = iv_union(wait), r = iv_union(run);
    double s = 0;
    size_t j = 0;
    for (const auto& x : w) {
        while (j < r.size() && r[j].second <= x.first) ++j;
        double cur = x.first;
        for (size_t k = j; k < r.size() && r[k].first < x.second; ++k) {
            if (r[k].first > cur) s += r[k].first - cur;
            cur = std::max(cur, r[k].second);
        }
        if (x.second > cur) s += x.second - cur;
    }
    return s;
}

AlignTimingScope::AlignTimingScope() : parent_(timing_top()) { timing_top() = this; }

AlignTimingScope::~AlignTimingScope() {
    timing_top() = parent_;
    if (parent_) parent_->t_.add(t_);
}

AlignTiming* AlignTimingScope::current() {
    AlignTimingScope* s = timing_top();
    return s ? &s->t_ : nullptr;
}

double steady_seconds() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ================================================================ report rows

void fill_align_row(ReportRow& row, const AlignResult& r, int64_t G) {
    row.repeat_frac = r.repeat_frac;
    row.outcome = r.outcome.empty() ? "confident" : r.outcome;
    if (!r.attempted || r.outcome == "prefiltered" || r.outcome == "aligner-failed" || r.outcome == "repeat-skipped") return;
    row.records = r.n_records;
    const ChainResult& ch = r.chain;
    // the 500 sub-record cap decided which placements the chain could use: say so in the row
    const std::string trunc = ch.truncated ? strf("truncated:%u/%u-sub-records", ch.n_candidates, ch.n_feasible) : std::string();
    if (ch.parts.empty()) {
        if (!trunc.empty()) row.dropped = trunc;
        return;
    }
    std::string b;
    for (size_t k = 0; k < ch.parts.size(); ++k) {
        const ChainPart& p = ch.parts[k];
        if (k) b.push_back(';');
        b += strf("%lld-%lld:%lld-%lld:%c:%.4f", (long long)p.qs, (long long)p.qe, (long long)p.ts, (long long)p.te, p.strand, part_identity(p, G));
    }
    row.blocks = b;
    row.frame = std::string(1, ch.frame);
    row.label = ch.strands == "+" ? "FWD" : ch.strands == "-" ? "INV" : ch.strands;
    row.tie = ch.tie;
    row.alt_score = ch.alt_score;
    row.alt_E = ch.alt_E;
    row.internal = ch.internal;
    std::string d;
    for (auto& x : ch.islands) {
        if (!d.empty()) d.push_back(';');
        d += strf("island-below-b:q%lld-%lld", (long long)x.first, (long long)x.second);
    }
    if (!trunc.empty()) d += (d.empty() ? "" : ";") + trunc;
    if (!d.empty()) row.dropped = d;
}

// ================================================================ Aligner

Aligner::Aligner(const Options& opt) : impl_(new Impl) {
    Impl& I = *impl_;
    I.opt = opt;
    I.version = minimap2_version(opt.minimap2);
    I.extra = split_ws(opt.mm_extra);
    I.budget.reset(new ProcessBudget(opt.jobs));
    if (!opt.inject_chains.empty()) I.read_inject(opt.inject_chains);
    else I.tmp.reset(new TempDir(opt.tmpdir));
}

Aligner::~Aligner() {}

const std::string& Aligner::version() const { return impl_->version; }

void Aligner::align_reference(const Graph& g, const SiteData& sd, std::vector<AlignResult>& res, std::vector<ReportRow>& rows) {
    Impl& I = *impl_;
    const Site& s = *sd.site;
    const std::string label = s.label(g);
    if (res.size() != sd.units.size()) res.resize(sd.units.size());
    std::vector<Item> items;
    for (const Unit& u : sd.units) {
        if (!u.outcome.empty()) {
            res[u.id].outcome = u.outcome;
            continue;
        }
        Item it;
        it.key = label + ":ref:" + std::to_string(u.id);
        it.alt_key = label + ":ref:" + excursion_str(g, u.exc);
        it.walk = u.exc.alts;
        it.target.ref = true;
        it.target.sn = s.sn;
        it.target.lo = u.wlo;
        it.target.hi = u.whi;
        it.feasible = [&sd](NodeId n, int64_t a, int64_t b) { return sd.allowed_of(n).contains(a, b); };
        it.res = &res[u.id];
        items.push_back(std::move(it));
    }
    if (items.empty()) return;
    I.run_items(g, items, label);
    for (const Unit& u : sd.units)
        if (u.outcome.empty()) fill_align_row(rows[u.id], res[u.id], I.opt.G);
}

void Aligner::align_jobs(const Graph& g, const std::vector<AlignJob>& jobs, std::vector<AlignResult>& res) {
    Impl& I = *impl_;
    res.assign(jobs.size(), AlignResult());
    if (jobs.empty()) return;
    std::vector<Item> items(jobs.size());
    for (size_t i = 0; i < jobs.size(); ++i) {
        items[i].key = jobs[i].key;
        items[i].walk = jobs[i].query;
        items[i].target = jobs[i].target;
        items[i].feasible = jobs[i].feasible;
        items[i].res = &res[i];
    }
    I.run_items(g, items, jobs[0].key + ".jobs");
}

void Aligner::check_failure_rate() const {
    const Impl& I = *impl_;
    uint64_t w = I.n_windows.load(), f = I.n_failed.load();
    if (I.n_units.load() > 0) {
        ZLOG("aligner: %llu units, %llu windows (%llu failed, %llu retried), %llu minimap2 processes (%llu killed by a signal), "
             "%.1f s in minimap2 (screen %.1f, pass 2 %.1f), peak RSS %.2f GB (minimap2's own figures)",
             (unsigned long long)I.n_units.load(), (unsigned long long)w, (unsigned long long)f, (unsigned long long)I.n_retries.load(),
             (unsigned long long)I.n_processes.load(), (unsigned long long)I.n_signals.load(), I.mm2_millis.load() / 1000.0,
             I.pass_millis[1].load() / 1000.0, I.pass_millis[2].load() / 1000.0, I.peak_rss.load() / 1048576.0);
        ZLOG("aligner: minimap2's own figures, summed over its processes: real %.1f s, CPU %.1f s; %.1f s waiting for an aligner slot "
             "(-j %d), summed over the processes",
             I.real_millis.load() / 1000.0, I.cpu_millis.load() / 1000.0, I.wait_millis.load() / 1000.0, I.budget->total());
        ZLOG("aligner: %llu screened, %llu prefiltered, %llu to pass 2, %llu confident, %llu paf-inconsistent, %llu truncated at %zu "
             "sub-records, %llu repeat-skipped",
             (unsigned long long)I.n_screened.load(), (unsigned long long)I.n_prefiltered.load(), (unsigned long long)I.n_pass2.load(),
             (unsigned long long)I.n_confident.load(), (unsigned long long)I.n_inconsistent.load(), (unsigned long long)I.n_truncated.load(),
             MAX_CHAIN_PARTS, (unsigned long long)I.n_repeat_skipped.load());
    }
    if (w > 0 && f * 20 > w)
        fail(EXIT_ALIGNER, strf("minimap2 failed on %llu of %llu windows (more than 5%%): a broken aligner install?", (unsigned long long)f,
                                (unsigned long long)w));
}

uint64_t Aligner::windows() const { return impl_->n_windows.load(); }
uint64_t Aligner::windows_failed() const { return impl_->n_failed.load(); }
uint64_t Aligner::paf_inconsistent() const { return impl_->n_inconsistent.load(); }

std::string minimap2_version(const std::string& path) {
    if (path.empty()) fail(EXIT_INPUT, "no minimap2 given: -m <minimap2> is required (there is no PATH default)");
    ProcResult r = run_process({path, "--version"});
    if (!r.ok()) fail(EXIT_INPUT, strf("cannot run %s --version: %s %s", path.c_str(), r.describe().c_str(), r.err_tail.c_str()));
    std::string v = r.out;
    while (!v.empty() && (v.back() == '\n' || v.back() == '\r' || v.back() == ' ')) v.pop_back();
    if (v.empty()) fail(EXIT_INPUT, strf("%s --version printed nothing", path.c_str()));
    return v;
}

} // namespace zip
