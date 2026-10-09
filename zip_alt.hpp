/*
  zip_alt.hpp -- spec step 5: alt-vs-alt (v2).

  After the reference pass has decided (SitePlanner::accept_reference), alleles that did not zip to
  the reference are aligned, under the same rules, to a parallel branch that shares their bounding
  handles and is private to them.  It runs in the same site plan: one piece map, acceptance (A), (C)
  and (D), and the transaction (SitePlanner::accept_alt, then finish), and the same aligner
  (Aligner::align_jobs: alt targets in walk orientation, the same records, chain, island rule, rule U
  and fragmentation gate).  A node that either pass zipped or used as a target is never again a
  query, a target or a bound (SitePlanner::node_used).

  Definitions
    residual excursion  a unit of kind F, I or BK none of whose nodes (anchors included) the reference
                        pass zipped or used as a target; walk = dep, alts, arr in canonical
                        orientation (anchors '+')
    rank                the excursion's SR: the highest SR of its alt nodes, i.e. when the path came
                        into existence -- for a creator excursion its run's SR (reused nodes are
                        older); for an observed (GAF) excursion its insertion, not the ranks of the
                        contigs that walk it
    shared handles      of a pair R, M with rank(R) < rank(M): the handles both walks pass in the same
                        orientation, whose node occurs once in each walk -- the anchors when they are
                        identical, and any alt handles the two walks share
    branch              of M: a maximal run of M's handles between two consecutive shared handles u
                        and v (in the same order in R), none of them on R; its target is R's run
                        between u and v.  Identical anchors give u = dep, v = arr.
    group (u, v)        representative: the branch of the lowest-rank residual creator excursion that
                        runs through u then v (ties: the longest branch, then the owner's SN, SO and
                        first node, then the unit); fallback witnesses and GAF-only excursions are
                        never representatives.  Members: the branches at (u, v) of higher-rank
                        residual excursions, relative to the representative.  A unit needs member
                        branch >= b and target >= ceil(b*i) (smaller ones get no row), so a branch
                        shorter than ceil(b*i) is no round's representative.
    rounds              up to --alt-rounds: in round r every group takes the next distinct branch in
                        its order whose nodes and bounds are unused (the next-lowest-rank branch at
                        the same (u, v); excursions sharing one branch count once); members that did
                        not merge regroup under it.

  Eligibility of a member (before alignment; refusals are report rows)
    shares-nodes  some excursion of the site walks both a member-branch node and a representative-
                  branch node (the member is not node-disjoint from the representative branch), or a
                  node is already used
    not-private   a side-aware flood from the member branch that never crosses u or v reaches a
                  reference node, or touches u on a side other than u's outgoing side, or v on a
                  side other than v's incoming side
    in-series     an alt-only path joins the member branch and the representative branch, in either
                  orientation (u and v may be crossed)
    pair-too-big  member or target > --max-pair
    site:capped   beyond --max-site-query of member query aligned at this site (all rounds and the
                  near-parallel diagnostics together), in order of member bp (descending)

  Alignment: the member branch (query) against the representative branch (target), both spelled in
  walk orientation; every target base is a representative node (empty feasibility function).  The
  target coordinates are branch offsets: the virtual backbone, which SitePlanner converts to node
  offsets and sides by orientation.  Confident chains go to SitePlanner::accept_alt in a total order:
  haplotype support first (complete walks only), then aligned bp, then sum '=', then the candidate
  key (u, v, member branch).

  Near-parallel (diagnostic, never zipped): residual excursions of query >= b whose windows' starts
  and ends both lie within NEAR_BP (200) of each other's are clustered (single linkage); in a cluster
  with more than one anchor pair, the first excursion of each other anchor pair is reported against
  the cluster's representative (as above, by rank and length), unless they share a node.  The row
  (pass "diag") has outcome near-parallel:<overlap|disjoint>:<private|entangled> (windows overlap or
  touch; the member is private to its own anchors) and the alignment of the two alt runs.

  Report rows (pass "alt", round r): owners of the member branch's excursions, kind of the first,
  anchors = the bounds "u>v", window = "0-<representative branch bp>" (the virtual backbone),
  query_bp = member branch bp, the align columns, kept_bp / dropped / outcome from the planner, and
  representative = "<contig>:<first node> <branch handles>".

  --dump DIR adds DIR/alt/<site>.tsv: every candidate and near-parallel pair with both branches,
  keyed like the edit dump (the planner key: the site's unit count plus the alt row's index, so it
  never equals a reference unit id).  AlignJob keys (--dump file names, --inject-chains):
  "<site>:alt:<round>:<u>><member branch>><v>" (a branch of more than 8 handles abbreviated) and
  "<site>:near:<member unit id>".
*/
#pragma once

#include "zip_edit.hpp"

#include <deque>

namespace zip {

// Alt-vs-alt counters of one site (main sums and logs them).
struct AltStats {
    uint64_t residual = 0;            // residual excursions
    uint64_t pairs = 0;               // (representative candidate, member) pairs sharing >= 2 handles (scanned once per site)
    uint64_t groups = 0;              // (round, u, v) groups with at least one candidate row
    uint64_t candidates = 0;          // member branches with a row (sizes qualify), over all rounds
    uint64_t shares_nodes = 0, not_private = 0, in_series = 0, pair_too_big = 0, capped = 0;
    uint64_t aligned = 0, confident = 0;
    uint64_t zipped = 0, refused = 0, reverted = 0;   // final (after SitePlanner::finish)
    int64_t zipped_bp = 0;            // kept bp of the zipped alt chains
    uint64_t near_parallel = 0;       // near-parallel rows
    void add(const AltStats& o);
};

// Spec step 5 for one site, after SitePlanner::accept_reference and before SitePlanner::finish.
// Appends one row per candidate (and per near-parallel pair) to `rows`, which must not reallocate
// while the planner holds pointers to them (hence a deque).  Thread-safe across sites.
void alt_pass(const Graph& g, const SiteData& sd, Aligner* aligner, SitePlanner& planner, const Options& opt,
              std::deque<ReportRow>& rows, AltStats& st);

// After SitePlanner::finish: count the final outcomes of the alt rows into st.
void alt_tally(const std::deque<ReportRow>& rows, AltStats& st);

} // namespace zip
