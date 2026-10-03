/*
  zip_common.hpp -- shared types and utilities for rgfa-zip.

  Contents:
    - handle / side encoding (NodeId, Handle, Side)
    - exit codes and the ZipError exception every module throws to stop the run
    - logging to stderr
    - Options: every command-line option, with the spec's defaults
    - Runtime: stopping cleanly -- the abort flag, the registry of running children and of
      temporary paths (a fatal error or a stop signal stops every thread and child)
    - Semaphore (the -j aligner-process budget)
    - run_threads / parallel_for: threads that degrade to fewer (or none) when a thread cannot start
    - run_process: run a child without a shell, capture stdout/stderr tail, exit status, signal
    - LineReader: line-by-line input, gzip read through a `gzip -dc` child (no zlib dependency)
    - AtomicWriter: <path>.tmp, every write checked, fsync, rename (exit 5 on any I/O error)
    - small helpers: strf, revcomp, number parsing and formatting

  Header-only.  Every module includes it.
*/
#pragma once

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <exception>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <dirent.h>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <spawn.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace zip {

static const char* const ZIP_VERSION = "0.1";

// ---------------------------------------------------------------- handles and sides
//
// Nodes are dense indices into Graph::nodes (input nodes sorted by numeric id).
// A handle is a node in an orientation: h = 2*node + (reverse ? 1 : 0).
// A side is one end of a node: s = 2*node + (right ? 1 : 0).
// Walking handle (n,+) enters n.L and leaves n.R; walking (n,-) enters n.R and leaves n.L.
// A GFA link 'L u uo v vo' joins side u.R (uo '+') or u.L (uo '-') to side v.L (vo '+') or v.R (vo '-').

typedef uint32_t NodeId;
typedef uint32_t Handle;
typedef uint32_t Side;
static const uint32_t NONE = UINT32_MAX;

inline Handle make_handle(NodeId n, bool rev) { return (n << 1) | (rev ? 1u : 0u); }
inline NodeId handle_node(Handle h) { return h >> 1; }
inline bool handle_rev(Handle h) { return (h & 1u) != 0; }
inline Handle flip(Handle h) { return h ^ 1u; }
inline Side left_side(NodeId n) { return n << 1; }
inline Side right_side(NodeId n) { return (n << 1) | 1u; }
inline NodeId side_node(Side s) { return s >> 1; }
inline bool side_is_right(Side s) { return (s & 1u) != 0; }
// the side a walk leaves through when it walks h, and the side it enters through
inline Side exit_side(Handle h) { return handle_rev(h) ? left_side(handle_node(h)) : right_side(handle_node(h)); }
inline Side entry_side(Handle h) { return handle_rev(h) ? right_side(handle_node(h)) : left_side(handle_node(h)); }
// the handle that leaves through side s, and the handle that enters through side s
inline Handle handle_leaving(Side s) { return make_handle(side_node(s), !side_is_right(s)); }
inline Handle handle_entering(Side s) { return make_handle(side_node(s), side_is_right(s)); }

// ---------------------------------------------------------------- exit codes, errors

enum ExitCode {
    EXIT_OK = 0,
    EXIT_INPUT = 2,       // invalid input (graph, snarls, command line); nothing written
    EXIT_INVARIANT = 3,   // invariant failure (global placement assert, --strict drop)
    EXIT_ALIGNER = 4,     // systemic aligner failure
    EXIT_IO = 5           // output I/O error
};

// Every fatal condition is a ZipError carrying its exit code.  main() catches it, removes any
// temporary output and exits with the code, so nothing is written on a non-zero exit.
struct ZipError : public std::runtime_error {
    int code;
    ZipError(int c, const std::string& m) : std::runtime_error(m), code(c) {}
};

[[noreturn]] inline void fail(int code, const std::string& msg) { throw ZipError(code, msg); }

// ---------------------------------------------------------------- formatting

inline std::string vstrf(const char* fmt, va_list ap) {
    va_list ap2;
    va_copy(ap2, ap);
    int n = vsnprintf(nullptr, 0, fmt, ap2);
    va_end(ap2);
    if (n < 0) return std::string();
    std::string s((size_t)n + 1, '\0');
    vsnprintf(&s[0], s.size(), fmt, ap);
    s.resize((size_t)n);
    return s;
}

inline std::string strf(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
inline std::string strf(const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    std::string s = vstrf(fmt, ap);
    va_end(ap);
    return s;
}

// fixed-precision double, locale-independent enough for our purposes (C locale is never changed)
inline std::string fmt_fixed(double x, int prec) { return strf("%.*f", prec, x); }

// the shortest %g form that reads back as exactly x (0.95 -> "0.95", 0.95004 -> "0.95004"), so
// that different values never print alike
inline std::string fmt_exact(double x) {
    for (int p = 1; p <= 17; ++p) {
        std::string s = strf("%.*g", p, x);
        if (strtod(s.c_str(), nullptr) == x) return s;
    }
    return strf("%.17g", x);
}

// ---------------------------------------------------------------- logging (stderr only)

inline std::mutex& log_mutex() { static std::mutex m; return m; }
inline int& log_level() { static int v = 1; return v; }

inline void log_msg(int level, const char* fmt, ...) __attribute__((format(printf, 2, 3)));
inline void log_msg(int level, const char* fmt, ...) {
    if (level > log_level()) return;
    va_list ap;
    va_start(ap, fmt);
    std::string s = vstrf(fmt, ap);
    va_end(ap);
    std::lock_guard<std::mutex> lk(log_mutex());
    fprintf(stderr, "[rgfa-zip] %s\n", s.c_str());
    fflush(stderr);
}
#define ZLOG(...) ::zip::log_msg(1, __VA_ARGS__)
#define ZDEBUG(...) ::zip::log_msg(2, __VA_ARGS__)

// ---------------------------------------------------------------- sequences

inline char comp_base(char c) {
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
        case 'N': return 'N'; case 'n': return 'n';
        case 'R': return 'Y'; case 'Y': return 'R'; case 'K': return 'M'; case 'M': return 'K';
        case 'S': return 'S'; case 'W': return 'W'; case 'B': return 'V'; case 'V': return 'B';
        case 'D': return 'H'; case 'H': return 'D';
        case 'r': return 'y'; case 'y': return 'r'; case 'k': return 'm'; case 'm': return 'k';
        case 's': return 's'; case 'w': return 'w'; case 'b': return 'v'; case 'v': return 'b';
        case 'd': return 'h'; case 'h': return 'd';
        default: return c;
    }
}
inline std::string revcomp(const std::string& s) {
    std::string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) r[s.size() - 1 - i] = comp_base(s[i]);
    return r;
}

// ---------------------------------------------------------------- number parsing (strict)

inline bool parse_int64(const std::string& s, int64_t& out) {
    if (s.empty()) return false;
    std::string t;
    for (char c : s) if (c != ',') t.push_back(c);   // allow 1,000,000
    if (t.empty()) return false;
    char* end = nullptr;
    errno = 0;
    long long v = strtoll(t.c_str(), &end, 10);
    if (errno != 0 || !end || *end != '\0') return false;
    out = (int64_t)v;
    return true;
}
inline bool parse_double(const std::string& s, double& out) {
    if (s.empty()) return false;
    char* end = nullptr;
    errno = 0;
    double v = strtod(s.c_str(), &end);
    if (errno != 0 || !end || *end != '\0' || !std::isfinite(v)) return false;
    out = v;
    return true;
}
// "8G", "1500M", "all" (0), plain bytes
inline bool parse_mem(const std::string& s, int64_t& out) {
    if (s == "all") { out = 0; return true; }
    if (s.empty()) return false;
    int64_t mult = 1;
    std::string t = s;
    char last = t.back();
    if (last == 'K' || last == 'k') mult = 1024LL;
    else if (last == 'M' || last == 'm') mult = 1024LL * 1024;
    else if (last == 'G' || last == 'g') mult = 1024LL * 1024 * 1024;
    else if (last == 'T' || last == 't') mult = 1024LL * 1024 * 1024 * 1024;
    if (mult != 1) t.pop_back();
    double v = 0;
    if (!parse_double(t, v) || v < 0) return false;
    out = (int64_t)(v * (double)mult);
    return true;
}

// ---------------------------------------------------------------- options

enum class Stage { DETECT = 0, ALIGN = 1, PLAN = 2, FULL = 3 };
enum class WalksKind { CREATOR, WITNESS, GAF };

// --region SN:lo-hi: a site is processed if its SN equals SN or ends with "#SN", and its
// closed span [lo(site), hi(site)] meets [lo, hi].
struct Region {
    std::string sn;
    int64_t lo = 0, hi = 0;
};

// Every option of the spec's CLI table plus the debug options.  Defaults are the spec's.
struct Options {
    // inputs / outputs
    std::string gfa_path, snarls_path;
    std::string out_path;            // -o
    std::string report_path;         // -r
    std::string minimap2;            // -m (required; there is no PATH default)

    // parallelism and resources
    int threads = 1;                 // -t: sites in parallel
    int jobs = 0;                    // -j: concurrent aligner processes (0: = -t), lowered so j x 1.5 GB fits --mem
    int64_t mem = 0;                 // --mem bytes (0 = all physical memory)
    std::string tmpdir;              // --tmpdir (default $TMPDIR, else /tmp)

    // thresholds
    int64_t b = 5000;                // -b: min chain bp; also min unit query
    double ident = 0.95;             // -i: min gap-compressed identity
    int64_t G = 50;                  // -G: I/D runs >= G split records into blocks
    int64_t min_piece = 1000;        // --min-piece: min aligned query per record
    int64_t island = 1000;           // --island: gap that separates island groups
    double delta = 0.005;            // --delta: tie tolerance of rule U
    bool ties_refuse = false;        // --ties resolve|refuse
    int64_t frag = 1;                // --frag: max internal stretches per 5 kb (floor 2)
    int64_t prefilter_nmatch = 2000; // --prefilter nmatch,coverage
    double prefilter_cov = 0.2;
    int screen_sample = 8;           // --screen-sample
    int64_t max_pair = 5000000;      // --max-pair: per-pair query and window cap
    int64_t max_site_query = 50000000;   // --max-site-query: feasible query per site
    int64_t max_site_nodes = 200000;     // --max-site-nodes: interior size cap

    // aligner
    std::string preset = "asm20";    // -x
    std::string mm_extra;            // -X (split on whitespace)

    // passes and walks
    bool alt = true;                 // --no-alt turns alt-vs-alt off (v1)
    int alt_rounds = 3;              // --alt-rounds
    WalksKind walks = WalksKind::CREATOR;   // --walks creator|witness|gaf:FILE
    std::string gaf_path;
    double max_repeat_frac = -1;     // --max-repeat-frac (negative: off)

    // modes
    bool detect_only = false;        // --detect-only: plan and report, emit the input graph unchanged
    bool check = false;              // --check: V2 over every anchored path at sites with <= 5000 paths
    bool strict = false;             // --strict: exit 3 on any dropped (reverted) chain, or a failed --audit-gaf
    int64_t id_base = -1;            // --id-base (default: max input id + 1)
    std::string audit_gaf;           // --audit-gaf FILE: audit the zips against observed walks (the release gate)

    // debug
    std::vector<Region> regions;     // --region SN:lo-hi (repeatable)
    Stage stage = Stage::FULL;       // --stage detect|align|plan|full
    std::string inject_chains;       // --inject-chains FILE: bypass alignment with given chains
    std::string dump_dir;            // --dump DIR: units, records, chains as TSV/FASTA

    // derived
    std::string command_line;
    int64_t min_window() const { return (int64_t)std::ceil((double)b * ident - 1e-9); }
    // windows/queries per site of the reference pass need --max-pair on both sides
};

// ---------------------------------------------------------------- stopping cleanly
//
// The first fatal error in any thread, or a stop signal (SIGTERM, SIGINT, SIGHUP; main's signal
// thread), calls Runtime::request_abort().  Every worker loop and every aligner launch checks
// aborting() and stops by throwing Aborted -- not an error of its own: the first real error is the
// one reported -- and every running child (minimap2, gzip, bgzip) gets SIGTERM, so that nobody
// waits for in-flight work.  A child is registered from its spawn until it is reaped, and is
// unregistered while it is still a zombie (waitid WNOWAIT), so a pid in the registry can never
// belong to another process.  Temporary paths (the aligner's directory, the outputs' .tmp files)
// are registered too: a stop signal normally lets the run unwind (RAII removes them), and removes
// them itself only if that takes too long.

struct Aborted {};   // thrown by a worker that stops because another thread failed or a signal arrived

inline void remove_tree(const std::string& path) {
    DIR* d = opendir(path.c_str());
    if (d) {
        std::vector<std::string> names;
        while (struct dirent* e = readdir(d)) {
            std::string n = e->d_name;
            if (n != "." && n != "..") names.push_back(n);
        }
        closedir(d);
        for (const std::string& n : names) {
            std::string p = path + "/" + n;
            struct stat st;
            if (lstat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) remove_tree(p);
            else unlink(p.c_str());
        }
    }
    rmdir(path.c_str());
}

class Runtime {
public:
    // never destroyed: the signal thread may use it while main exits
    static Runtime& get() {
        static Runtime* r = new Runtime();
        return *r;
    }
    bool aborting() const { return abort_.load(); }
    // stop every worker and SIGTERM every running child (idempotent)
    void request_abort() {
        abort_.store(true);
        std::lock_guard<std::mutex> lk(m_);
        for (pid_t p : children_) kill(p, SIGTERM);
    }
    void add_child(pid_t p) {
        std::lock_guard<std::mutex> lk(m_);
        children_.push_back(p);
        if (abort_.load()) kill(p, SIGTERM);
    }
    void remove_child(pid_t p) {
        std::lock_guard<std::mutex> lk(m_);
        auto it = std::find(children_.begin(), children_.end(), p);
        if (it != children_.end()) children_.erase(it);
    }
    // temporary paths: removed by remove_paths() (a stop signal that cannot wait for the unwinding)
    void add_path(const std::string& p, bool dir) {
        std::lock_guard<std::mutex> lk(m_);
        paths_.push_back(std::make_pair(p, dir));
    }
    void remove_path(const std::string& p) {
        std::lock_guard<std::mutex> lk(m_);
        for (size_t i = paths_.size(); i-- > 0;)
            if (paths_[i].first == p) { paths_.erase(paths_.begin() + (long)i); break; }
    }
    void remove_paths() {
        std::vector<std::pair<std::string, bool>> v;
        {
            std::lock_guard<std::mutex> lk(m_);
            v = paths_;
        }
        for (const auto& p : v) {
            if (p.second) remove_tree(p.first);
            else unlink(p.first.c_str());
        }
    }
    int signal_received() const { return sig_.load(); }
    void set_signal(int s) { sig_.store(s); }
private:
    Runtime() {}
    std::atomic<bool> abort_{false};
    std::atomic<int> sig_{0};
    std::mutex m_;
    std::vector<pid_t> children_;
    std::vector<std::pair<std::string, bool>> paths_;
};

inline bool aborting() { return Runtime::get().aborting(); }
inline void check_abort() {
    if (aborting()) throw Aborted();
}

// ---------------------------------------------------------------- semaphore (aligner process budget)

class Semaphore {
public:
    explicit Semaphore(int n) : count_(n) {}
    void acquire() {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return count_ > 0; });
        --count_;
    }
    void release() {
        { std::lock_guard<std::mutex> lk(m_); ++count_; }
        cv_.notify_one();
    }
private:
    std::mutex m_;
    std::condition_variable cv_;
    int count_;
};

struct SemaphoreGuard {
    Semaphore& s;
    explicit SemaphoreGuard(Semaphore& s_) : s(s_) { s.acquire(); }
    ~SemaphoreGuard() { s.release(); }
    SemaphoreGuard(const SemaphoreGuard&) = delete;
    SemaphoreGuard& operator=(const SemaphoreGuard&) = delete;
};

// ---------------------------------------------------------------- threads

// Run worker(t) on n threads, t = 0..k-1 for the k threads that start.  A thread that cannot start
// (a process or thread limit) does not stop the run: it continues with the threads that started,
// and with none, worker(0) runs on the calling thread.  Callers distribute work dynamically and
// store results by index, so the outcome does not depend on how many threads ran.  worker must not
// throw.  Returns the number of threads used.
template <typename W>
int run_threads(int n, const W& worker) {
    std::vector<std::thread> pool;
    pool.reserve((size_t)std::max(0, n));
    for (int t = 0; t < n; ++t) {
        try {
            pool.emplace_back(worker, t);
        } catch (const std::exception& e) {
            log_msg(1, "warning: could not start thread %d of %d (%s); continuing with %d", t + 1, n, e.what(), std::max(1, t));
            break;
        }
    }
    if (pool.empty()) worker(0);
    for (std::thread& th : pool) th.join();
    return std::max<int>(1, (int)pool.size());
}

// Run fn(i, thread) for i in [0, n) on up to `threads` threads, taking indices in increasing order
// (dynamic scheduling).  Callers store results by index, so the outcome does not depend on the
// thread count.  The first exception thrown stops the loop (and requests an abort, so that other
// loops and the running children stop too) and is rethrown here; a loop stopped by an abort
// requested elsewhere throws Aborted.
template <typename F>
void parallel_for(size_t n, int threads, F&& fn) {
    if (threads <= 1 || n <= 1) {
        for (size_t i = 0; i < n; ++i) {
            check_abort();
            fn(i, 0);
        }
        return;
    }
    std::atomic<size_t> next(0);
    std::atomic<bool> stop(false), saw_abort(false);
    std::exception_ptr err;
    std::mutex em;
    int nt = (int)std::min<size_t>((size_t)threads, n);
    run_threads(nt, [&](int t) {
        while (!stop.load()) {
            if (aborting()) {
                saw_abort = true;
                break;
            }
            size_t i = next.fetch_add(1);
            if (i >= n) break;
            try {
                fn(i, t);
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
    });
    if (err) std::rethrow_exception(err);
    if (saw_abort) throw Aborted();
}

// ---------------------------------------------------------------- child processes

struct ProcResult {
    bool spawned = false;        // false: the program could not be started (spawn_error says why)
    int spawn_error = 0;
    bool exited = false;         // normal exit
    int exit_code = -1;
    bool signaled = false;       // killed by a signal (OOM kill shows up as SIGKILL)
    int signal = 0;
    // peak RSS from wait4.  A posix_spawn child shares rgfa-zip's memory until exec, and the kernel
    // carries that high-water mark over, so this is max(rgfa-zip's RSS, the child's): the aligner
    // replaces it with minimap2's own "Peak RSS" figure
    long max_rss_kb = 0;
    std::string out;             // captured stdout (only when stdout_path is empty)
    std::string err_tail;        // the last err_keep bytes of stderr
    bool ok() const { return spawned && exited && exit_code == 0; }
    std::string describe() const {
        if (!spawned) return "could not be started";
        if (signaled) return strf("killed by signal %d", signal);
        if (exited) return strf("exit status %d", exit_code);
        return "unknown status";
    }
};

// posix_spawn attributes for every child: an empty signal mask (rgfa-zip blocks the stop signals
// in its threads, and a mask survives exec) and default dispositions for SIGPIPE (which rgfa-zip
// ignores), SIGCHLD and the stop signals (an ignored disposition survives exec too)
struct SpawnAttr {
    posix_spawnattr_t a;
    SpawnAttr() {
        posix_spawnattr_init(&a);
        sigset_t none, def;
        sigemptyset(&none);
        sigemptyset(&def);
        for (int s : {SIGPIPE, SIGCHLD, SIGTERM, SIGINT, SIGHUP}) sigaddset(&def, s);
        posix_spawnattr_setsigmask(&a, &none);
        posix_spawnattr_setsigdefault(&a, &def);
        posix_spawnattr_setflags(&a, POSIX_SPAWN_SETSIGMASK | POSIX_SPAWN_SETSIGDEF);
    }
    ~SpawnAttr() { posix_spawnattr_destroy(&a); }
    SpawnAttr(const SpawnAttr&) = delete;
    SpawnAttr& operator=(const SpawnAttr&) = delete;
};

// posix_spawnp with SpawnAttr; a started child is registered with the Runtime (reap it with
// reap_child).  Returns 0 or the spawn error.
inline int spawn_child(pid_t& pid, const std::vector<std::string>& argv, const posix_spawn_file_actions_t* fa) {
    std::vector<char*> cargv;
    for (const std::string& s : argv) cargv.push_back(const_cast<char*>(s.c_str()));
    cargv.push_back(nullptr);
    SpawnAttr at;
    pid = -1;
    int rc = posix_spawnp(&pid, cargv[0], fa, &at.a, cargv.data(), environ);
    if (rc == 0) Runtime::get().add_child(pid);
    return rc;
}

// Wait for a registered child to exit, unregister it while it is still a zombie (so a registry
// pid is never reused), then reap it.  False if the child cannot be waited for (e.g. ECHILD: it
// was reaped elsewhere) -- its status is then unknown.
inline bool reap_child(pid_t pid, int& status, struct rusage* ru) {
    siginfo_t si;
    while (true) {
        memset(&si, 0, sizeof(si));
        if (waitid(P_PID, (id_t)pid, &si, WEXITED | WNOWAIT) == 0) break;
        if (errno == EINTR) continue;
        Runtime::get().remove_child(pid);
        return false;
    }
    Runtime::get().remove_child(pid);
    struct rusage tmp;
    memset(&tmp, 0, sizeof(tmp));
    while (wait4(pid, &status, 0, ru ? ru : &tmp) < 0) {
        if (errno != EINTR) return false;
    }
    return true;
}

// Run argv[0] (searched on PATH only if it has no '/') with stdin from /dev/null.  stdout goes to
// stdout_path when given (created/truncated), otherwise it is captured.  stderr is drained and its
// tail kept.  Never uses a shell.  Thread-safe (posix_spawn; all pipe fds are close-on-exec).
// Throws Aborted, without starting anything, once an abort was requested.
inline ProcResult run_process(const std::vector<std::string>& argv, const std::string& stdout_path = std::string(),
                              size_t err_keep = 8192) {
    ProcResult r;
    if (argv.empty()) return r;
    check_abort();
    int outp[2] = {-1, -1}, errp[2] = {-1, -1};
    if (stdout_path.empty() && pipe2(outp, O_CLOEXEC) != 0) fail(EXIT_IO, strf("pipe: %s", strerror(errno)));
    if (pipe2(errp, O_CLOEXEC) != 0) {
        if (outp[0] >= 0) { close(outp[0]); close(outp[1]); }
        fail(EXIT_IO, strf("pipe: %s", strerror(errno)));
    }
    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
    if (stdout_path.empty()) posix_spawn_file_actions_adddup2(&fa, outp[1], 1);
    else posix_spawn_file_actions_addopen(&fa, 1, stdout_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    posix_spawn_file_actions_adddup2(&fa, errp[1], 2);
    pid_t pid = -1;
    int rc = spawn_child(pid, argv, &fa);
    posix_spawn_file_actions_destroy(&fa);
    if (outp[1] >= 0) close(outp[1]);
    close(errp[1]);
    if (rc != 0) {
        if (outp[0] >= 0) close(outp[0]);
        close(errp[0]);
        r.spawn_error = rc;
        r.err_tail = strf("%s: %s", argv[0].c_str(), strerror(rc));
        return r;
    }
    r.spawned = true;
    // drain stdout (if captured) and stderr without deadlock
    std::string err;
    char buf[65536];
    struct pollfd pf[2];
    int nfd = 0;
    int out_i = -1, err_i = -1;
    if (outp[0] >= 0) { pf[nfd].fd = outp[0]; pf[nfd].events = POLLIN; out_i = nfd++; }
    pf[nfd].fd = errp[0]; pf[nfd].events = POLLIN; err_i = nfd++;
    int open_fds = nfd;
    while (open_fds > 0) {
        // once an abort is requested (the child got SIGTERM) stop draining: a grandchild of a
        // wrapper script may hold the pipes open long after the child itself is gone
        int pr = poll(pf, (nfds_t)nfd, 200);
        if (pr == 0) {
            if (aborting()) break;
            continue;
        }
        if (pr < 0) { if (errno == EINTR) continue; break; }
        for (int i = 0; i < nfd; ++i) {
            if (pf[i].fd < 0 || !(pf[i].revents & (POLLIN | POLLHUP | POLLERR))) continue;
            ssize_t n = read(pf[i].fd, buf, sizeof(buf));
            if (n < 0 && errno == EINTR) continue;
            if (n <= 0) { close(pf[i].fd); pf[i].fd = -1; --open_fds; continue; }
            if (i == out_i) r.out.append(buf, (size_t)n);
            else if (i == err_i) {
                err.append(buf, (size_t)n);
                if (err.size() > 4 * err_keep + 65536) err.erase(0, err.size() - err_keep);
            }
        }
    }
    for (int i = 0; i < nfd; ++i) if (pf[i].fd >= 0) close(pf[i].fd);
    if (err.size() > err_keep) err.erase(0, err.size() - err_keep);
    r.err_tail = err;
    int st = 0;
    struct rusage ru;
    memset(&ru, 0, sizeof(ru));
    if (!reap_child(pid, st, &ru)) {
        r.exited = false;
        return r;
    }
    r.max_rss_kb = ru.ru_maxrss;
    if (WIFEXITED(st)) { r.exited = true; r.exit_code = WEXITSTATUS(st); }
    else if (WIFSIGNALED(st)) { r.signaled = true; r.signal = WTERMSIG(st); }
    return r;
}

// ---------------------------------------------------------------- line input (plain or gzip)

// Reads a text file line by line.  A gzip file (detected by its magic bytes) is read through a
// `gzip -dc` child, so no zlib dependency is needed.  Any failure is a ZipError(EXIT_INPUT).
class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path), buf_(1 << 20) {
        int fd = open(path.c_str(), O_RDONLY | O_CLOEXEC);
        if (fd < 0) fail(EXIT_INPUT, strf("cannot open %s: %s", path.c_str(), strerror(errno)));
        unsigned char magic[2] = {0, 0};
        ssize_t n = pread(fd, magic, 2, 0);
        if (n == 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
            ::close(fd);
            int p[2];
            if (pipe2(p, O_CLOEXEC) != 0) fail(EXIT_INPUT, strf("pipe: %s", strerror(errno)));
            posix_spawn_file_actions_t fa;
            posix_spawn_file_actions_init(&fa);
            posix_spawn_file_actions_addopen(&fa, 0, "/dev/null", O_RDONLY, 0);
            posix_spawn_file_actions_adddup2(&fa, p[1], 1);
            int rc = spawn_child(pid_, {"gzip", "-dc", "--", path}, &fa);
            posix_spawn_file_actions_destroy(&fa);
            ::close(p[1]);
            if (rc != 0) {
                ::close(p[0]);
                pid_ = -1;
                fail(EXIT_INPUT, strf("cannot run gzip to read %s: %s", path.c_str(), strerror(rc)));
            }
            fd_ = p[0];
        } else {
            fd_ = fd;
        }
    }
    ~LineReader() {
        if (fd_ >= 0) ::close(fd_);
        if (pid_ > 0) {
            kill(pid_, SIGTERM);
            int st;
            reap_child(pid_, st, nullptr);
        }
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;

    // next line without its '\n' (and a trailing '\r'); false at end of input
    bool next(std::string& line) {
        line.clear();
        while (true) {
            if (pos_ < len_) {
                char* s = &buf_[pos_];
                char* e = (char*)memchr(s, '\n', len_ - pos_);
                if (e) {
                    line.append(s, (size_t)(e - s));
                    pos_ = (size_t)(e - &buf_[0]) + 1;
                    ++lineno_;
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    return true;
                }
                line.append(s, len_ - pos_);
                pos_ = len_;
            }
            if (eof_) {
                if (line.empty()) return false;
                ++lineno_;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            check_abort();   // once per buffer: a stop signal does not wait for a whole input file
            ssize_t n = read(fd_, &buf_[0], buf_.size());
            if (n < 0) {
                if (errno == EINTR) continue;
                fail(EXIT_INPUT, strf("read error on %s: %s", path_.c_str(), strerror(errno)));
            }
            if (n == 0) eof_ = true;
            pos_ = 0;
            len_ = (size_t)n;
        }
    }
    // close and check the decompressor's status (a status that cannot be read is a failure too)
    void close() {
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
        if (pid_ > 0) {
            int st = 0;
            bool waited = reap_child(pid_, st, nullptr);
            pid_ = -1;
            if (!waited) fail(EXIT_INPUT, strf("cannot read the exit status of gzip reading %s", path_.c_str()));
            if (!(WIFEXITED(st) && WEXITSTATUS(st) == 0))
                fail(EXIT_INPUT, strf("gzip failed while reading %s (corrupt or truncated?)", path_.c_str()));
        }
    }
    uint64_t line_number() const { return lineno_; }
    const std::string& path() const { return path_; }
private:
    std::string path_;
    int fd_ = -1;
    pid_t pid_ = -1;
    std::vector<char> buf_;
    size_t pos_ = 0, len_ = 0;
    bool eof_ = false;
    uint64_t lineno_ = 0;
};

// split a line on tabs (no copies of the line itself are kept by the caller)
inline void split_tabs(const std::string& line, std::vector<std::string>& out) {
    out.clear();
    size_t p = 0;
    while (true) {
        size_t q = line.find('\t', p);
        if (q == std::string::npos) { out.push_back(line.substr(p)); break; }
        out.push_back(line.substr(p, q - p));
        p = q + 1;
    }
}

// ---------------------------------------------------------------- atomic, checked output

// Writes <path>.tmp, checking every write; finish() flushes, fsyncs and closes; commit() renames
// onto <path> and fsyncs the directory; rollback() removes a committed <path> again (a later
// output failed to commit).  Destroying an uncommitted writer removes the temp file, which is also
// registered with the Runtime for a stop signal.  Any I/O error is a ZipError(EXIT_IO).
class AtomicWriter {
public:
    explicit AtomicWriter(const std::string& path) : path_(path), tmp_(path + ".tmp") {
        fd_ = open(tmp_.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
        if (fd_ < 0) fail(EXIT_IO, strf("cannot create %s: %s", tmp_.c_str(), strerror(errno)));
        Runtime::get().add_path(tmp_, false);
        buf_.reserve(cap_);
    }
    ~AtomicWriter() {
        if (fd_ >= 0) ::close(fd_);
        if (!committed_) unlink(tmp_.c_str());
        Runtime::get().remove_path(tmp_);
    }
    // remove the committed output (all-or-nothing commits of several outputs)
    void rollback() {
        if (committed_) unlink(path_.c_str());
    }
    bool committed() const { return committed_; }
    AtomicWriter(const AtomicWriter&) = delete;
    AtomicWriter& operator=(const AtomicWriter&) = delete;

    void write(const char* p, size_t n) {
        if (buf_.size() + n > cap_) flush_buf();
        if (n > cap_) { raw_write(p, n); return; }
        buf_.append(p, n);
    }
    void write(const std::string& s) { write(s.data(), s.size()); }
    void put(char c) { if (buf_.size() + 1 > cap_) flush_buf(); buf_.push_back(c); }

    void finish() {
        if (finished_) return;
        flush_buf();
        if (fsync(fd_) != 0) fail(EXIT_IO, strf("fsync %s: %s", tmp_.c_str(), strerror(errno)));
        if (::close(fd_) != 0) { fd_ = -1; fail(EXIT_IO, strf("close %s: %s", tmp_.c_str(), strerror(errno))); }
        fd_ = -1;
        finished_ = true;
    }
    void commit() {
        finish();
        if (rename(tmp_.c_str(), path_.c_str()) != 0)
            fail(EXIT_IO, strf("rename %s -> %s: %s", tmp_.c_str(), path_.c_str(), strerror(errno)));
        committed_ = true;
        std::string dir = ".";
        size_t sl = path_.rfind('/');
        if (sl != std::string::npos) dir = sl == 0 ? "/" : path_.substr(0, sl);
        int dfd = open(dir.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (dfd >= 0) {
            if (fsync(dfd) != 0 && errno != EINVAL && errno != EROFS)
                { ::close(dfd); fail(EXIT_IO, strf("fsync directory %s: %s", dir.c_str(), strerror(errno))); }
            ::close(dfd);
        }
    }
    const std::string& path() const { return path_; }
private:
    void raw_write(const char* p, size_t n) {
        while (n > 0) {
            ssize_t w = ::write(fd_, p, n);
            if (w < 0) {
                if (errno == EINTR) continue;
                fail(EXIT_IO, strf("write %s: %s", tmp_.c_str(), strerror(errno)));
            }
            p += w;
            n -= (size_t)w;
        }
    }
    void flush_buf() {
        if (!buf_.empty()) { raw_write(buf_.data(), buf_.size()); buf_.clear(); }
    }
    std::string path_, tmp_;
    int fd_ = -1;
    bool finished_ = false, committed_ = false;
    std::string buf_;
    size_t cap_ = 4u << 20;
};

// mkdir -p; false on failure
inline bool make_dirs(const std::string& path) {
    if (path.empty()) return false;
    std::string cur;
    size_t p = 0;
    while (p <= path.size()) {
        size_t q = path.find('/', p);
        if (q == std::string::npos) q = path.size();
        cur = path.substr(0, q);
        if (!cur.empty()) {
            struct stat st;
            if (stat(cur.c_str(), &st) != 0) {
                if (mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) return false;
            } else if (!S_ISDIR(st.st_mode)) {
                return false;
            }
        }
        p = q + 1;
    }
    return true;
}

inline int64_t physical_memory_bytes() {
    long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || psz <= 0) return 0;
    return (int64_t)pages * (int64_t)psz;
}

inline long self_peak_rss_kb() {
    struct rusage ru;
    memset(&ru, 0, sizeof(ru));
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss;
}

} // namespace zip
