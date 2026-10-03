/*
  zip_site.hpp -- spec steps 1-3: snarls, sites, covering walks, units, kinds and Allowed.

  Step 1  read_snarls (real JSON), build_sites: top-level snarls with rank-0 boundaries on one SN;
          A = lower-SO boundary, window [lo, hi) = [end(A), start(B)); side-aware interior flood
          from A.R and B.L that never enters A or B; skips site:leak / site:ref-mismatch /
          site:too-big.
  Step 2  build_runs (global), CreatorWalks (creator excursion per run, witness fallback),
          WitnessWalks (fallback for every run), the WalkSource interface (GafWalks lives in
          zip_gaf), canonical frame, deduplication.
  Step 3  kinds F/I/BK/J, reference units, Allowed by monotone propagation to a fixpoint (done
          exactly, in O(H+E), over the SCC condensation of the site's alt-handle graph), the exact
          pre-filter (infeasible) and the per-site cap (site:capped).

  Everything after this step reads only SiteData: the site, its alt-handle graph, Allowed and the
  ordered unit list.
*/
#pragma once

#include "zip_graph.hpp"

#include <unordered_map>

namespace zip {

// ================================================================ snarls

struct Snarl {
    NodeId start = NONE, end = NONE;
    bool start_back = false, end_back = false;
    bool nested = false;             // the line has "parent"
    uint64_t line = 0;
};

// Read `vg snarls -n -P <ref> | vg view -Rj -` output (plain or gzip).  Every line is parsed as
// real JSON; only start/end "name" and "backward", and whether "parent" is present, are read.
// An unparseable line, a missing name, or a name not in the graph is ZipError(EXIT_INPUT).
std::vector<Snarl> read_snarls(const std::string& path, const Graph& g);

// ================================================================ sites

enum class SiteStatus : uint8_t { OK, NO_ALT, LEAK, REF_MISMATCH, TOO_BIG };
// "ok", "no-alt", "site:leak", "site:ref-mismatch", "site:too-big"
const char* site_status_name(SiteStatus s);

struct Site {
    uint32_t index = 0;              // position in the canonical site list (SN, lo, A id, B id)
    NodeId A = NONE, B = NONE;       // boundaries; A has the lower SO
    int32_t sn = -1;
    int64_t lo = 0, hi = 0;          // window [lo, hi) = [end(A), start(B)), absolute SN coordinates
    SiteStatus status = SiteStatus::OK;
    uint64_t snarl_line = 0;         // first snarls line naming this site
    std::vector<NodeId> interior;    // side-aware interior, sorted; excludes A and B (cleared for skipped sites)
    std::vector<NodeId> backbone;    // rank-0 interior nodes in SO order; they tile [lo, hi)
    std::vector<NodeId> alts;        // non-rank-0 interior nodes, sorted
    size_t n_interior = 0;           // interior size (kept for skipped sites)
    int64_t alt_bp = 0;              // alt bp of the interior (whole for site:too-big; a lower bound for a leak past the cap)

    bool contains(NodeId n) const { return std::binary_search(interior.begin(), interior.end(), n); }
    bool is_alt(NodeId n) const { return std::binary_search(alts.begin(), alts.end(), n); }
    // A, B or an interior rank-0 node: a reference end of an excursion
    bool is_anchor_node(const Graph& g, NodeId n) const { return n == A || n == B || (g.is_ref(n) && contains(n)); }
    std::string label(const Graph& g) const { return g.name(A) + ".." + g.name(B); }
};

struct SiteStats {
    uint64_t snarls = 0, nested = 0, nonref_boundary = 0, cross_sn = 0, same_node = 0, duplicate = 0;
    uint64_t sites = 0, ok = 0, no_alt = 0, leak = 0, ref_mismatch = 0, too_big = 0;
    uint64_t pure_insertion = 0;          // OK sites with lo == hi
    uint64_t orientation_disagrees = 0;   // boundary 'backward' flags disagree with the SO order (informational)
    // alt (non-rank-0) input nodes, and those in no site's interior (any status)
    uint64_t alt_nodes = 0, uncovered_nodes = 0;
    int64_t alt_bp = 0, uncovered_bp = 0;
};

// All sites (top-level snarls whose boundaries are rank-0 nodes on one SN) in canonical order, each
// flooded and classified, on opt.threads threads.  Every site is built (the leak rate is a
// property of the whole input, independent of --region).  Fails with EXIT_INPUT when more than 3
// sites and more than 1% of sites leak (one bad snarl makes up to three leaking sites: itself and
// the neighbours sharing its boundary nodes), or when two processed sites overlap.  Counts the alt
// nodes that lie in no site (main checks that the snarls cover the graph).
std::vector<Site> build_sites(const Graph& g, const std::vector<Snarl>& snarls, const Options& opt, SiteStats& st);

// ================================================================ the site's alt-handle graph

// A reference end of a link from/to an alt handle.
//   departure (reference handle before the alt handle): (r,+) leaves r.R at end(r); (r,-) leaves r.L at start(r)
//   arrival   (reference handle after the alt handle):  (r,+) enters r.L at start(r); (r,-) enters r.R at end(r)
// Coordinates are absolute SN coordinates (end(A) = lo, start(B) = hi).
struct Anchor {
    Handle ref = 0;              // global handle of A, B or a backbone node
    bool minus = false;          // sign '-'
    int64_t coord = 0;
};

// CSR over the site's alt handles.  Local node l is site.alts[l]; local handle lh = 2*l + rev.
struct SiteGraph {
    std::vector<NodeId> alt;                          // local -> global (== site.alts)
    std::unordered_map<NodeId, uint32_t> local;       // global alt node -> local index
    std::vector<uint32_t> succ_off, succ;             // alt -> alt steps (local handles)
    std::vector<uint32_t> pred_off, pred;
    std::vector<uint32_t> dep_off, arr_off;
    std::vector<Anchor> dep, arr;                     // departures into / arrivals from each local handle
    std::vector<uint32_t> scc;                        // SCC id of each local handle (ids in reverse topological order:
    uint32_t n_scc = 0;                               // an edge c1 -> c2 between SCCs has c2 < c1)

    size_t n_nodes() const { return alt.size(); }
    size_t n_handles() const { return alt.size() * 2; }
    Handle global(uint32_t lh) const { return make_handle(alt[lh >> 1], (lh & 1u) != 0); }
    uint32_t local_node(NodeId n) const {
        auto it = local.find(n);
        return it == local.end() ? NONE : it->second;
    }
    uint32_t local_handle(Handle h) const {
        uint32_t l = local_node(handle_node(h));
        return l == NONE ? NONE : 2 * l + (handle_rev(h) ? 1u : 0u);
    }
};

void build_site_graph(const Graph& g, const Site& s, SiteGraph& sg);

// ================================================================ Allowed

enum class AllowedStatus : uint8_t { NONE, OK, EMPTY, BLOCKED };
// "none", "ok", "empty", "blocked"
const char* allowed_status_name(AllowedStatus s);

// Allowed(n): the intersection of the windows of every anchored alt-only path through n, in
// either orientation.  BLOCKED: mixed signs (some path is a junction); EMPTY: insertion or loop
// (no target); NONE: no departure or no arrival reaches n.
struct AllowedIv {
    AllowedStatus status = AllowedStatus::NONE;
    bool minus = false;          // the paths are all '-' (canonical flip of a reverse-walked allele)
    int64_t lo = 0, hi = 0;      // absolute SN coordinates [lo, hi) (also set for EMPTY)
    bool observed = false;       // replaced by a complete WalkSource from observed walks
    // may a piece land on target [a, b)?
    bool contains(int64_t a, int64_t b) const { return status == AllowedStatus::OK && lo <= a && b <= hi; }
};

// Strict Allowed for every alt node of the site (indexed by local node).
void compute_allowed(const Site& s, const SiteGraph& sg, std::vector<AllowedIv>& out);

// ================================================================ runs (global)

struct Runs {
    std::vector<std::vector<NodeId>> runs;   // alt nodes of one SN, contiguous SO, ++ linked; canonical order (SR, SN, SO)
    std::vector<uint32_t> run_of;            // node -> run index (NONE for rank-0 nodes)
    std::vector<uint32_t> pos_in_run;        // node -> index within its run
};
Runs build_runs(const Graph& g);

// ================================================================ excursions

enum class Src : uint8_t { CREATOR = 0, GAF = 1, WITNESS = 2 };
const char* src_name(Src s);   // "creator", "gaf", "witness"

// who an excursion belongs to: a run (creator / witness), or an observed haplotype contig (GAF)
struct Owner {
    std::string contig;          // SN of the run, or the observed contig's name
    int32_t sr = -1;             // the run's rank (-1 if unknown)
    int64_t so = 0;              // SO of the run's first node, or the walk's start on the contig
    NodeId first = NONE;         // the run's first node (NONE for GAF owners)
    uint32_t run = NONE;         // run index (NONE for GAF owners)
    Src src = Src::CREATOR;
};
// (sr, contig, so, first, src)
bool owner_less(const Owner& a, const Owner& b);
bool owner_equal(const Owner& a, const Owner& b);

// A maximal run of alt handles between two reference handles of one site.
struct Excursion {
    Handle dep = 0, arr = 0;     // reference handles before and after (A, B or backbone)
    std::vector<Handle> alts;    // alt handles in walk orientation
    std::vector<Owner> owners;
    uint32_t weight = 0;         // haplotype support (GafWalks); 0 when unknown
    Src src = Src::CREATOR;      // best source among the owners (CREATOR < GAF < WITNESS)
};

// canonical frame: when dep and arr are both '-', reverse the run, flip every handle, and use
// (flip(arr), flip(dep)) as anchors.  Returns true if it flipped.
bool canonicalize(Excursion& e);
// lexicographic order on (dep, arr, alts): the canonical unit key
bool excursion_key_less(const Excursion& a, const Excursion& b);
bool excursion_key_equal(const Excursion& a, const Excursion& b);
// "s1+>s5+,s6->s2+" (dep > alts > arr)
std::string excursion_str(const Graph& g, const Excursion& e);

// ================================================================ walk sources

struct WalkStats {
    uint64_t runs = 0, creator_ok = 0, fallback = 0, witness_dag = 0, witness_bubble = 0, unanchored = 0;
    uint64_t fail_creator_links = 0, fail_own_ambiguous = 0, fail_creator_edge = 0, fail_cycle = 0, fail_other = 0;
    uint64_t excursions_raw = 0, excursions = 0, merged = 0;
    void add(const WalkStats& o);
};

struct SiteData;

// The covering walks of a site.  Everything downstream reads only excursion lists plus the Allowed
// oracle, so a new source needs no other change.
class WalkSource {
public:
    virtual ~WalkSource() {}
    virtual const char* name() const = 0;
    // Excursions of this site, in any frame and possibly repeated; prepare_site canonicalizes them,
    // merges identical ones (owners merged, weights summed) and orders them.  sd.site, sd.sg and
    // sd.allowed (strict) are ready when this is called.  Must be thread-safe (const).
    virtual std::vector<Excursion> excursions(const SiteData& sd, WalkStats& st) const = 0;
    // true: the source's walks are complete (observed haplotypes), so they also define Allowed for
    // the nodes they touch, and acceptance orders by support first
    virtual bool complete() const { return false; }
    // called after excursions() when complete(): narrow or replace Allowed of observed nodes
    // (allowed is indexed by local node of sd.sg)
    virtual void refine_allowed(const SiteData& sd, std::vector<AllowedIv>& allowed) const {
        (void)sd;
        (void)allowed;
    }
    // Acceptance rule (G), the whole-record rule: does a record of this source that walks node n
    // also read, anywhere in the same record, reference of the site's SN that overlaps [ta, tb) by
    // more than the GAF audit's read tolerance?  Zipping a piece of n onto [ta, tb) would then
    // collapse that haplotype's second copy.  Only a source of whole observed records (GafWalks)
    // can tell; Allowed cannot, as it sees one excursion at a time.  piece_bp is the length of the
    // piece of n.  Must be thread-safe (const).
    virtual bool reads_target(const Site& s, NodeId n, int64_t piece_bp, int64_t ta, int64_t tb) const {
        (void)s;
        (void)n;
        (void)piece_bp;
        (void)ta;
        (void)tb;
        return false;
    }
};

// One creator excursion per run, from its creator links (L lines whose SR equals the run's rank),
// extended outward to the first rank-0 handle: at a reused node, leave by the unique link whose SR
// equals the run's rank, else follow the reused node's own run to its creator link.  Ambiguity, a
// revisited handle (the extension loops), or a missing creator link sends the run to the witness
// fallback (src WITNESS).
class CreatorWalks : public WalkSource {
public:
    CreatorWalks(const Graph& g, const Runs& runs) : g_(g), runs_(runs) {}
    const char* name() const override { return "creator"; }
    std::vector<Excursion> excursions(const SiteData& sd, WalkStats& st) const override;
private:
    const Graph& g_;
    const Runs& runs_;
};

// The witness of every run: the anchored path through the run with the most bp not yet covered
// (two DP passes over the site's alt-handle graph in DFS topological order, back edges dropped;
// the shortest-bp anchored bubble when no such path exists; ties to lower SR, then name).
class WitnessWalks : public WalkSource {
public:
    WitnessWalks(const Graph& g, const Runs& runs) : g_(g), runs_(runs) {}
    const char* name() const override { return "witness"; }
    std::vector<Excursion> excursions(const SiteData& sd, WalkStats& st) const override;
private:
    const Graph& g_;
    const Runs& runs_;
};

// the runs with a node in this site, in canonical run order
std::vector<uint32_t> site_runs(const Runs& runs, const Site& s);

// ================================================================ units

enum class Kind : uint8_t { F, I, BK, J };
const char* kind_name(Kind k);   // "F", "I", "BK", "J"

// Kind of a canonical excursion and its window.  With dep = r1+ and arr = r2+: F if end(r1) <
// start(r2) (window [end(r1), start(r2))), I if equal, BK if greater (wlo > whi); J if the anchor
// signs differ (wlo = whi = 0).
Kind classify(const Graph& g, const Excursion& e, int64_t& wlo, int64_t& whi);

// One canonical, deduplicated excursion of a site.
struct Unit {
    uint32_t id = 0;             // index in SiteData::units (canonical order)
    Excursion exc;
    Kind kind = Kind::F;
    int64_t wlo = 0, whi = 0;    // window (absolute); F: wlo < whi; I: wlo == whi; BK: wlo > whi
    int64_t query_bp = 0;        // alt bp of the walk (the query spelled in walk orientation)
    int64_t feasible_bp = 0;     // bp of walk nodes whose Allowed is OK and overlaps the window (F only)
    bool ref_unit = false;       // F, query >= b, window >= ceil(b*i): a reference-pass unit
    // detect-stage outcome:
    //   ""              reference unit to align
    //   "pair-too-big"  reference unit with query or window > --max-pair
    //   "infeasible"    reference unit whose feasible bp < b (never aligned)
    //   "site:capped"   reference unit beyond the site's --max-site-query
    //   "no-window"     query >= b but no reference target: kind I, BK or J, or an F window < ceil(b*i)
    //   "small"         query < b (not reported)
    std::string outcome;
    int64_t window_bp() const { return kind == Kind::J ? 0 : whi - wlo; }
};

// Everything the later stages read about one site.
struct SiteData {
    const Site* site = nullptr;
    SiteGraph sg;
    std::vector<AllowedIv> allowed;      // per local node of sg
    std::vector<Unit> units;             // canonical order: (first owner SR, contig, SO, first node), then key
    bool complete_walks = false;         // the walk source is complete (GafWalks): acceptance orders by support first
    const WalkSource* walks = nullptr;   // the walk source (acceptance asks it rule (G)); outlives the SiteData
    WalkStats wstats;
    // Allowed of any node; a NONE interval for nodes that are not alt nodes of the site
    const AllowedIv& allowed_of(NodeId n) const;
};

// Spec steps 2-3 for one site: alt-handle graph, strict Allowed, excursions from `src`
// (canonical, deduplicated, ordered), kinds, windows, feasible bp, outcomes and the cap.
void prepare_site(const Graph& g, const Site& s, const WalkSource& src, const Options& opt, SiteData& sd);

} // namespace zip
