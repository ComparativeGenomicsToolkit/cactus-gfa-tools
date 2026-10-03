/*
  zip_align.hpp -- spec step 4: aligner runner, PAF/CIGAR, blocks, projection, feasibility,
  chain DP, island rule, rule U and the fragmentation gate.

  Coordinates
    query   walk coordinates: offset in the unit's walk spelled in walk orientation (Unit::exc.alts),
            0 <= q <= Unit::query_bp.  Walk handle k covers [off_k, off_k + len(node)).
    target  reference pass: ABSOLUTE SN coordinates (the PAF's window-relative coordinates plus the
            window start), so they compare directly with Allowed and the window.
            alt-vs-alt (v2): offsets in the representative branch spelled in walk orientation
            (the virtual backbone; spec §5).
  A record/part on strand '-' aligns rc(query) to the target: its CIGAR is in target order and
  consumes the query from qe downward (PAF convention).

  Flow per site (align_reference): units with an empty outcome -> one minimap2 index per distinct
  window -> screen (pass 1, no CIGAR) -> pass 2 (-c --eqx) -> verify every record against both
  sequences (paf-inconsistent must stay 0) -> drop records < --min-piece or gap-compressed identity
  < -i -> blocks at I/D runs >= G -> project blocks onto walk nodes (left-snap P) -> cut records to
  feasible sub-records (targets inside Allowed(node)) -> chain DP (frames F/R) -> island rule ->
  U (ties, parsimony, strand, phase/copy ties) -> fragmentation gate -> AlignResult.

  Outcomes written by this module (AlignResult::outcome; the report row gets the same string, or
  "confident" when the chain passed every gate):
    ""                  a confident chain (the edit module's candidates)
    "prefiltered"       the screen (pass 1) found too little homology
    "aligner-failed"    minimap2 exited non-zero or wrote unparseable output for the window
    "paf-inconsistent"  a record disagrees with the sequences or its own coordinates (counter; must stay 0)
    "no-chain"          no feasible sub-record survived the filters
    "below-b"           no chain keeps >= b aligned query after the island rule
    "ambiguous-strand"  rule U: tied placements of equal parsimony differ in strand string
    "tie-refused"       --ties refuse and rule U found a tie
    "fragmented"        more internal unaligned stretches >= G than max(2, frag x aligned / 5 kb)
    "repeat-skipped"    --max-repeat-frac: the window is more repetitive than allowed (not aligned)
  Units the detect stage already decided keep their outcome: res[k].outcome = unit outcome, so that
  "" always means a confident chain.

  Normalisation: minimap2's --eqx writes '=' for any two bases with equal nt4 codes, so an N facing
  an N (or any IUPAC code facing another) is '='.  Verified records rewrite such columns as 'X':
  n_eq counts only identical A/C/G/T columns, which is minimap2's own nmatch (PAF column 10).
  minimap2's block length (column 11) leaves out every ambiguous base, those of I runs (query) and
  D runs (target) included; verification recounts it the same way.  Identity is unchanged by that:
  '=' / ('=' + 'X' + bases in I/D runs shorter than G), ambiguous bases counted as X or as gap.

  --dump DIR: per site, DIR/align/<site label>.{records,chains,windows}.tsv and the window/query
  FASTA (DIR/align/<site label>.w<k>.{t,q}.fa).  Main writes DIR/units.tsv and DIR/allowed.tsv.

  --inject-chains FILE (no minimap2 runs at all): TSV lines
      <key>  <record|chain>  <qs>  <qe>  <strand>  <ts>  <te>  <cigar>
  key = "<site label>:ref:<unit id>" (or "<site label>:ref:<excursion_str>") for reference units,
  AlignJob::key for jobs.  Query in walk coordinates, target in target coordinates (absolute for a
  reference window).  CIGAR ops =, X, I, D, or M (resolved to =/X against the sequences); every
  line is verified like a PAF record.  'record' lines go through the whole chain pipeline (filters,
  feasibility, DP, U, gates); 'chain' lines ARE the chain (split into blocks at I/D runs >= G, no
  filter, no feasibility cut, tie "injected").  A unit with no line gets "no-chain".
*/
#pragma once

#include "zip_report.hpp"

#include <functional>

namespace zip {

// ---------------------------------------------------------------- alignment data

struct CigarOp {
    uint32_t len = 0;
    char op = '=';                     // '=', 'X', 'I', 'D' (pass 2 runs with --eqx)
};

// A verified PAF record of one unit.
struct PafRecord {
    int64_t qs = 0, qe = 0;            // query interval, walk coordinates
    int64_t ts = 0, te = 0;            // target interval, target coordinates
    char strand = '+';
    int64_t n_eq = 0, n_x = 0;         // '=' and 'X' columns (n_eq: identical A/C/G/T only; see Normalisation)
    int64_t n_ins = 0, n_del = 0;      // inserted / deleted bases
    double gci = 0;                    // gap-compressed identity: '=' / ('=' + 'X' + bases in I/D runs < G)
    std::vector<CigarOp> ops;          // target order
};

// One chain element: a feasible sub-record, i.e. a stretch of one block (a record split at I/D
// runs >= G) whose pieces all land inside Allowed.  Blocks of different records are never mixed.
struct ChainPart {
    uint32_t record = 0;               // index into AlignResult::records
    uint32_t block = 0;                // block index within the record
    int64_t qs = 0, qe = 0;            // query interval (walk coordinates)
    int64_t ts = 0, te = 0;            // target interval (target coordinates)
    char strand = '+';
    int64_t n_eq = 0;                  // exact '=' count inside the part (the chain score)
    std::vector<CigarOp> ops;          // the part's own alignment, target order (all I/D runs < G)
};

// A piece: the part of one walk node that one chain part aligns, and its target under the left-snap
// projection P (unsnapped; snapping is the edit's step 2).
struct Piece {
    uint32_t part = 0;                 // index into ChainResult::parts
    uint32_t walk_pos = 0;             // index of the node's handle in Unit::exc.alts
    NodeId node = NONE;
    int64_t a = 0, b = 0;              // node-forward interval, 0 <= a < b <= len(node)
    int64_t qa = 0, qb = 0;            // the same interval in walk coordinates
    int64_t ta = 0, tb = 0;            // target [P(qa), P(qb)) for a '+' part, [P(qb), P(qa)) for '-'
    char rel = '+';                    // '+' if the walk orientation of the node equals the part's strand
};

struct ChainResult {
    std::vector<ChainPart> parts;      // in query order
    std::vector<Piece> pieces;         // every walk node overlapping every part, in query order
    char frame = 'F';                  // 'F': segment hulls increase along the query; 'R': decrease
    std::string strands;               // segment strand string, e.g. "+-+" (report label in walk orientation)
    int64_t score = 0;                 // sum of '=' over the parts
    int64_t aligned_q = 0;             // aligned query bp after the island rule
    std::string tie;                   // unique | tie-kept | tie-resolved | ambiguous-strand | injected
    int64_t alt_score = -1;            // best alternative chain's score, -1 if none
    int64_t alt_E = -1;                // and its E
    int64_t internal = 0;              // internal unaligned stretches >= G on either axis
    bool island_dropped = false;       // the island rule removed a group (island-below-b; that query stays alt)
    // ---- added by the align module
    int64_t E = -1;                    // out-of-place stretches of this chain (rule U's parsimony)
    std::vector<std::pair<int64_t, int64_t>> islands;   // query hulls of the groups the island rule removed
    int64_t island_bp = 0;             // their aligned query bp
    uint32_t n_candidates = 0;         // feasible sub-records offered to the DP (after the 500 cap)
    uint32_t n_feasible = 0;           // feasible sub-records before the cap
    bool truncated = false;            // more than 500 feasible sub-records: the lowest-scoring were dropped
                                       // (the report row's dropped column says "truncated:500/<n_feasible>-sub-records")
    std::string u_trace;               // --dump only: rule U's candidates (role, score, E, tied, strands, hull), one per line
};

// Result of aligning one unit (or one alt-vs-alt job).
struct AlignResult {
    // "" when a confident chain was found (chain is valid); otherwise the refusal (see the header)
    std::string outcome;
    bool attempted = false;            // the unit was handed to the aligner (or to --inject-chains)
    // verified records (all of them, before the --min-piece / -i filter).  The aligner keeps them
    // after chaining only with --dump (their per-column ops are most of rgfa-zip's memory at the
    // array tangles); n_records always holds their number (the report's records column)
    std::vector<PafRecord> records;
    int64_t n_records = 0;
    ChainResult chain;
    std::string aligner_err;           // stderr tail when aligner-failed
    // ---- added by the align module
    double repeat_frac = -1;           // repeat fraction of the target window (diagnostic), -1 if not computed
    bool screened = false;             // the query went through pass 1
};

// A generic job (alt-vs-alt, v2): a query walk against a target, with its own feasibility oracle.
struct AlignTarget {
    bool ref = true;                   // true: reference window; false: oriented alt walk (virtual backbone)
    int32_t sn = -1;                   // ref: SN and window [lo, hi), absolute
    int64_t lo = 0, hi = 0;
    std::vector<Handle> walk;          // alt: the target walk in walk orientation
};
struct AlignJob {
    std::string key;                   // stable name (dumps, --inject-chains): "<site>:<pass>:<id>"
    std::vector<Handle> query;         // query walk (alt handles, walk orientation)
    AlignTarget target;
    // may the piece of `node` land on target [ta, tb)?  (reference pass: Allowed(node) contains it)
    // An empty function allows everything.
    std::function<bool(NodeId node, int64_t ta, int64_t tb)> feasible;
};

// ---------------------------------------------------------------- timing (stderr only)

// Wall-clock accounting of aligner work under the -j budget.  The windows of all sites share at
// most -j minimap2 processes, so a window's or a site's wall time mixes running with waiting for a
// free process slot.  An AlignTiming records, for the windows it covers, when their minimap2
// processes ran and when one of them waited for a slot (intervals on the steady clock, seconds),
// and minimap2's own "Real time" and "CPU" (its last stderr line), summed over the processes.
struct AlignTiming {
    std::vector<std::pair<double, double>> run, wait;
    double mm2_real = 0, mm2_cpu = 0;
    uint64_t processes = 0;
    void add(const AlignTiming& o);
    double running() const;            // seconds during which at least one of the processes ran
    double waiting() const;            // seconds during which one waited for a slot and none ran
};

// While an AlignTimingScope lives, every Aligner call made on its thread adds the timing of its
// windows to it; a scope that closes adds its timing to the enclosing one.  main opens one per
// site, zip_alt one for the alt-vs-alt pass.
class AlignTimingScope {
public:
    AlignTimingScope();
    ~AlignTimingScope();
    AlignTimingScope(const AlignTimingScope&) = delete;
    AlignTimingScope& operator=(const AlignTimingScope&) = delete;
    const AlignTiming& timing() const { return t_; }
    static AlignTiming* current();     // the innermost open scope of this thread, or nullptr
private:
    AlignTiming t_;
    AlignTimingScope* parent_;
};

// the steady clock in seconds (AlignTiming's time base)
double steady_seconds();

// ---------------------------------------------------------------- the aligner

class Aligner {
public:
    // Checks the -m binary (version), creates an RAII temp directory under --tmpdir / $TMPDIR, and
    // sets the process budget to opt.jobs (already lowered for --mem by main).  Reads
    // --inject-chains (ZipError(EXIT_INPUT) on a malformed file).
    explicit Aligner(const Options& opt);
    ~Aligner();
    Aligner(const Aligner&) = delete;
    Aligner& operator=(const Aligner&) = delete;

    const std::string& version() const;

    // Reference pass of one site.  Every unit with an empty outcome is aligned against its window
    // and chained; res[k] and the align fields + outcome of rows[k] are filled, k = unit id
    // (res.size() == rows.size() == sd.units.size()).  Units with a non-empty outcome keep their
    // rows; their res[k].outcome is set to the unit's outcome.  With --inject-chains, chains come
    // from the file instead of minimap2.  Thread-safe: called concurrently by the -t site threads;
    // at most opt.jobs minimap2 processes run at any time.  A signal kills -> retry the window once,
    // alone; a second failure throws ZipError(EXIT_ALIGNER).
    void align_reference(const Graph& g, const SiteData& sd, std::vector<AlignResult>& res, std::vector<ReportRow>& rows);

    // Generic form for alt-vs-alt: align and chain every job under the same rules.  Jobs with the
    // same target share one index.  res is resized to jobs.size() and filled in job order.
    void align_jobs(const Graph& g, const std::vector<AlignJob>& jobs, std::vector<AlignResult>& res);

    // After every site: ZipError(EXIT_ALIGNER) when more than 5% of windows failed.  Also logs the
    // aligner summary (processes, screened, prefiltered, retries, truncations, peak RSS).
    void check_failure_rate() const;

    // stderr summary counters
    uint64_t windows() const;
    uint64_t windows_failed() const;
    uint64_t paf_inconsistent() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// The align fields of a report row from a result: records, blocks, frame, label, tie, alt_score,
// alt_E, internal, repeat_frac, island-below-b entries in dropped, and outcome ("confident" for a
// confident chain).  align_reference fills unit rows with it; main uses it for diagnostic rows.
void fill_align_row(ReportRow& row, const AlignResult& r, int64_t G);

// `<path> --version` trimmed (e.g. "2.30-r1287"); ZipError(EXIT_INPUT) if it cannot be run.
std::string minimap2_version(const std::string& path);

// Left-snap projection P(x): the target coordinate reached when the query bases of the part up to
// walk coordinate x have been consumed ('+': from qs upward; '-': from qe downward).  A D run that
// starts exactly at x is not consumed, so both sides of a junction land on one point.
// Requires qs <= x <= qe.  This is the one definition of P; the edit module uses it too.
int64_t project_left_snap(const ChainPart& p, int64_t x);

// P at many points at once (one pass over the ops); identical to project_left_snap at each point.
std::vector<int64_t> project_points(const ChainPart& p, const std::vector<int64_t>& xs);

// ---------------------------------------------------------------- chaining (exposed for tests and zip_alt)

// The query of one unit/job: its walk and the walk offset of every handle.
struct QueryLayout {
    std::vector<Handle> walk;          // walk orientation
    std::vector<int64_t> off;          // off[k] = walk offset of walk[k]; off[walk.size()] = query length
    static QueryLayout of(const Graph& g, const std::vector<Handle>& walk);
    int64_t length() const { return off.empty() ? 0 : off.back(); }
    int64_t node_len(size_t k) const { return off[k + 1] - off[k]; }
};

typedef std::function<bool(NodeId node, int64_t ta, int64_t tb)> FeasibleFn;

// Everything after verification, on one unit's verified records: drop records < --min-piece or
// gap-compressed identity < -i, split into blocks at I/D runs >= G, cut blocks to feasible
// sub-records (dropping those < --min-piece), keep the 500 best, chain DP (frames F/R), island
// rule, rule U, fragmentation gate; fills out.chain and out.outcome ("" when confident).
// out.records must already hold the verified records.  Target extent [tlo, thi) is the window
// (reference) or [0, len) of the target walk.  A null feasible function allows everything.
void chain_records(const QueryLayout& ql, int64_t tlo, int64_t thi, const FeasibleFn& feasible, const Options& opt,
                   AlignResult& out);

// The pieces of a chain: every walk node overlapping every part (P projection, node-forward
// intervals, rel).  chain_records calls it for the picked chain.
void make_pieces(const QueryLayout& ql, ChainResult& ch);

// Split a verified record into blocks at I/D runs >= G (blocks with no aligned column are dropped).
std::vector<ChainPart> record_blocks(const PafRecord& r, uint32_t record_index, int64_t G);

// E of a set of parts: unaligned query stretches > G, uncovered target stretches > G (within
// [tlo, thi)), plus strand switches of the strand string.
int64_t chain_E(const std::vector<ChainPart>& parts, int64_t qlen, int64_t tlo, int64_t thi, int64_t G);

} // namespace zip
