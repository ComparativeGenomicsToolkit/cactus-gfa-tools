/*
  zip_edit.hpp -- spec step 6: pieces, snapping, acceptance (A)-(D), cuts, splitter, pinch,
  transaction with validators V1-V7, --check, id allocation, and applying the plans.

  Contract
    - Every decision is made on the input graph; Graph::out() (CSR) is the input adjacency.
    - A site's plan is a pure function of (input site, accepted chain set), so sites are planned
      in parallel (one SitePlanner per site, on the site's thread) and the plans commute.
    - Until ids are allocated, a piece is identified by (original node, node-forward offset).
    - Ids are allocated only after every site is planned (apply_plans), in (SN, lo, original id,
      offset) order, from --id-base (default: Graph::max_id + 1).  SR is never edited.
    - Site boundaries A and B are never cut or re-pointed; every link with an end in the site keeps
      both ends on site pieces (interior, A.R, B.L).

  Target spaces
    A chain's target coordinates live in a target space: an oriented walk of nodes with
    coordinates.  The reference pass has one space per site: A+, the backbone, B+, in absolute SN
    coordinates (targets lie in [lo, hi]).  Alt-vs-alt (v2) adds one space per representative branch:
    u, the branch, v, in walk offsets (targets lie in [0, branch length]).  "L side of the piece
    starting at T" and "R side of the piece ending at T" are read in the space's walk orientation
    and converted to node sides (a '-' element swaps them), so the same pinch serves both passes.

  Testing the edit without alignment
    --inject-chains is read by zip_align (its 'chain' lines ARE the chain: no filter, no
    feasibility cut), so the edit module sees injected chains exactly like aligned ones.
    Debug hook (test Z29; an environment variable, not a CLI option):
        RGFA_ZIP_DEBUG_DROP_LINK="<site label>:<unit id>[,...]"
    leaves out the first translated link (canonical order) of that unit's chain whenever the chain
    is in the edit; the transaction must then revert that chain alone (reverted:V2).
*/
#pragma once

#include "zip_align.hpp"

namespace zip {

// ---------------------------------------------------------------- the plan of one site

// A piece of an input node: the piece of `node` that starts at node-forward offset `off`
// (off 0 for a node that is not cut).
struct PlanEnd {
    NodeId node = NONE;
    int64_t off = 0;
    bool right = false;                // R side (else L side)
};

// A link the plan adds.
struct PlanLink {
    PlanEnd a, b;
    int32_t sr = -1;
    uint32_t src = NONE;               // input link it comes from (overlap and tags are copied); NONE: split link
    uint8_t mode = 0;                  // 0: a re-targeted input link, written as the input wrote it
                                       // 1: a new piece-to-piece link, written a+ b+ (split; SR of the node)
                                       // 2: a translated link, written canonically (lower id first)
};

// An input node the plan cuts and/or deletes: its pieces in offset order.
struct PlanNode {
    NodeId node = NONE;
    std::vector<int64_t> off;          // piece start offsets (off[0] == 0)
    std::vector<int64_t> len;
    std::vector<uint8_t> keep;         // 0: the piece is zipped (deleted)
};

// A zipped piece, for reports and the GAF audit (zip_gaf's GafAuditPiece has the same fields).
struct ZippedPiece {
    uint32_t unit = NONE;              // unit id (reference pass) or alt candidate key (v2, with alt set)
    bool alt = false;                  // v2: target is an alt branch (ta/tb are branch offsets)
    NodeId node = NONE;
    int64_t a = 0, b = 0;              // node-forward
    int64_t ta = 0, tb = 0;            // target (absolute SN coordinates for the reference pass)
    char rel = '+';
};

// The edit of one site, decided on the input graph.
struct SitePlan {
    uint32_t site = NONE;              // Site::index
    int32_t sn = -1;                   // the site's SN and window start (id allocation order)
    int64_t lo = 0;
    uint32_t accepted = 0;             // chains committed
    uint32_t accepted_new = 0;         // ... of which zip bp of their own (the others only agree with an earlier chain's zip)
    uint32_t reverted = 0;             // chains dropped by the transaction (reverted:<check>)
    int64_t zipped_bp = 0;             // alt bp removed (zipped onto a target) = the sum of the zipped rows' kept_bp
    uint32_t g_chains = 0;             // rule (G) (observed walks): chains that lost a piece to it,
    uint32_t g_pieces = 0;             //   the pieces,
    int64_t g_bp = 0;                  //   and their bp
    int64_t target_bp = 0;             // target bp those pieces landed on
    uint32_t links_added = 0, links_removed = 0, cut_nodes = 0;
    std::vector<PlanNode> nodes;       // changed nodes, by id
    std::vector<uint32_t> removed_links;   // input links the plan replaces (sorted)
    std::vector<PlanLink> links;       // links the plan adds
    std::vector<ZippedPiece> zipped;   // every zipped piece (unit order, then walk order)
    bool empty() const { return nodes.empty() && links.empty() && removed_links.empty(); }
};

// ---------------------------------------------------------------- v2 hook

// One alt-vs-alt candidate (zip_alt): a member branch aligned against a representative branch.
// The chain's target coordinates are offsets in the representative branch spelled in walk
// orientation (the virtual backbone, spec section 5).  Privacy, disjointness and "not in series"
// are zip_alt's checks; acceptance applies (A), (C) and (D), with the representative branch as
// the target in place of (B).
struct AltCandidate {
    std::vector<Handle> query;         // member branch, walk orientation (alt handles)
    Handle u = 0, v = 0;               // bounds, walk orientation: the target walk is u, target..., v
    std::vector<Handle> target;        // representative branch between u and v, walk orientation
    const AlignResult* res = nullptr;  // outcome "" with a chain; target coordinates are branch offsets
    ReportRow* row = nullptr;          // kept_bp, dropped and outcome are filled in
    uint32_t key = 0;                  // stable tie-break (e.g. the job index); also ZippedPiece::unit
    int round = 1;
};

// ---------------------------------------------------------------- the planner

// Acceptance, transaction and validation for one site (spec section 6 steps 1-8).
class SitePlanner {
public:
    SitePlanner(const Graph& g, const SiteData& sd, const Options& opt);
    ~SitePlanner();
    SitePlanner(const SitePlanner&) = delete;
    SitePlanner& operator=(const SitePlanner&) = delete;

    // Reference pass.  Candidates are the units whose res[k].outcome is "" (a confident chain),
    // taken greedily in the spec's total order: kept bp, then sum '=', then the canonical unit key
    // (with a complete WalkSource, haplotype support first).  As in the prototype, kept bp is
    // weighted by node support: each piece counts (its bp) x (the excursion owners -- creator
    // runs or observed contigs -- through its node), so a chain that serves many haplotypes is
    // not blocked by a private one.  Pieces failing (A) one fate per
    // node interval, (B) Allowed after snapping, (C) disjoint targets, or (with observed walks)
    // (G) a GAF record through the node also reads the target are dropped; a chain is refused by
    // (D) or as trimmed-below-b.  Fills kept_bp (the bp its new pieces remove; 0 when it removes
    // nothing), dropped and outcome of rows[k] (indexed by unit id) for every candidate.
    void accept_reference(std::vector<AlignResult>& res, std::vector<ReportRow>& rows);

    // v2: alt-vs-alt candidates, after accept_reference, in the order given (zip_alt orders them).
    // Their rows stay owned by the caller and must outlive finish().  Candidates with an identical
    // target walk (u, representative branch, v) share one target space, so the members of one
    // group are compared by (A) and (C) like chains of one reference window, and the
    // representative stays their common target.  A candidate whose query, bound or target node is
    // used (node_used; a target node may be the target of its own space) gets "shares-nodes".
    void accept_alt(std::vector<AltCandidate>& cands);

    // A node zipped or used as a target so far -- query nodes of accepted chains, the reference
    // nodes a reference chain lands on, and the representative branch of an accepted alt chain:
    // alt-vs-alt never uses it again as a query, a target or a bound.
    bool node_used(NodeId n) const;

    // Build the edited site and run V1-V7; on failure re-plan incrementally in acceptance order,
    // dropping each chain whose addition fails (its row's outcome becomes "reverted:<check>").
    // Diagnostic rows go to extra_rows.  Returns the plan.
    SitePlan finish(std::vector<ReportRow>& rows, std::vector<ReportRow>& extra_rows);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// After every site is planned: allocate piece ids s<N> in (SN, lo, original id, offset) order and
// apply every plan to g (split nodes, translate links of deleted pieces, delete them).  The
// global placement assert (check_placement) runs after this, in main.
void apply_plans(Graph& g, std::vector<SitePlan>& plans, const Options& opt);

} // namespace zip
