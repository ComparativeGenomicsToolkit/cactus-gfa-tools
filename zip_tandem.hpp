/*
  zip_tandem.hpp -- --tandem: replace a long insertion that is a forward tandem copy of the reference
  beside it by one back link over the reference copy.

  The spec declines this ("no-window": copy collapse outside the window, and forward cycles).  It is
  an opt-in exception for copy-number variants: once minigraph stores an extra copy of an array
  unit as new sequence, it is a second route for every later haplotype with the same copy number
  (a module deletion edge plus the inserted copy spells the same thing as the reference), and they
  split between the routes.  A back link instead makes the extra copy a second pass over the
  reference copy, so every haplotype has one route for its copy number.

  Per site, after the site's plan (tandem_site, in the site threads):
    1. Candidates: units of kind I or F whose detect outcome is no-window or infeasible (never
       aligned to their own window), with query >= b and query - window >= b.  Anchors are '+'
       (the canonical frame); p1 = end(dep), p2 = start(arr).  A query that is mostly its own
       repeats (repeat_fraction > 0.5: a VNTR expansion, or several units at once) is not aligned
       ("tandem:self-repeat"): one copy of a long unit has little internal repetition, and these
       are the costly alignments.
    2. Each candidate's query is aligned (the Aligner, records kept) to the flank
       [p1 - X, p2 + X), X = 1.1 query + 1 kb, with no feasibility cut, and chained (frame F).
    3. Copy shapes, all '+', touching the insertion point (p1 or p2) within G:
         right    the chain covers the query and starts at p1 or p2:  U = ref[s, d), back link d -> p2
         left     the chain covers the query and ends at p1 or p2:    U = ref[c, e), back link p1 -> c
         rotated  the chain is the head (starting at p1 or p2) or     U = ref[s, d) + ref[c, e),
                  the tail (ending at p1 or p2); one record carries   back link d -> c
                  the rest, clipped where the chain starts or ends,
                  or the rest is shorter than --min-piece and matches
                  the reference beside the insertion point by edit
                  distance at identity >= -i
       The haplotype ref[..p1) + U + ref[p2..) is then spelled through the reference copy.  A back
       link joins the R side of the reference ending at its first position to the L side of the one
       starting at its second.  An F unit whose copy starts at p2 or ends at p1 lacks its window and
       also needs the skip link p1 -> p2.  The loop is >= b and <= 2 x query, and at most max(2G, 5%)
       of the query is unaligned.  Several copies of one unit in one insertion (multi-copy) are not
       handled.
    4. A copy is accepted with its whole alt component: the alt nodes of the site connected (over
       alt-alt links) to its walk.  Every unit through the component must be a verified copy, every
       node of it must be on one of their walks (a node no walk covers, e.g. a partial GAF record's
       or an unanchored run's, may carry novel sequence), and the site's plan must not touch it (cut
       or zip a node, or add or remove a link on one: an uncut alt-vs-alt representative keeps its
       member's links); otherwise the component is refused and stays as it is.  Differences between
       the copies and the reference copy (bubbles inside the component) are dropped: the graph keeps
       one copy, and cactus aligns each haplotype's bases to it.

  After apply_plans (apply_tandem, serial): delete each accepted component and its links, then add
  the links.  Back links go shortest loop first, and a loop with the period of one already added
  (within max(G, 1%)), overlapping it by at least half the period, whose phase difference reads the
  same one period later (the stretch between the two starts against the stretch one period on, at
  identity >= -i), is not added: copies of one array unit placed at different phases are one route.
  Loops at separate units of a long array are kept (each is a route for the same copy number).  A loop end that is not a reference node
  boundary snaps to one within G (unless that would close the loop), or splits that reference node
  (new pieces get the next free ids; the parent is deleted, as apply_plans does).  Link SR is the
  lowest rank in the component.  check_placement runs after.

  Report: the units of an accepted component get outcome "tandem" ("would-tandem" with
  --detect-only), frame "T", label right|left|rotated, the copy blocks in `blocks`, and the
  component's bp as kept_bp on its first unit.  A verified copy in a refused component gets
  "tandem-refused:mixed", "tandem-refused:unwalked" or "tandem-refused:edited".  A candidate that is
  not a copy keeps its outcome, with "tandem:<why>" in `dropped`.
*/
#pragma once

#include "zip_align.hpp"
#include "zip_edit.hpp"

#include <map>

namespace zip {

struct TandemStats {
    uint64_t candidates = 0;           // units offered (kind I/F, no-window or infeasible, query and net insertion >= b)
    uint64_t copies = 0;               // verified copies
    uint64_t right = 0, left = 0, rotated = 0;
    uint64_t components = 0;           // accepted components
    uint64_t refused_mixed = 0, refused_unwalked = 0, refused_edited = 0;   // refused components
    uint64_t units = 0;                // units in accepted components
    uint64_t nodes = 0;                // alt nodes removed
    int64_t bp = 0;                    // ... and their bp
    std::map<std::string, uint64_t> why;   // candidates that are not copies, by reason
    // apply_tandem
    uint64_t back_links = 0, skip_links = 0, existing = 0, equivalent = 0, snapped = 0, splits = 0;
    void add(const TandemStats& o);
};

// One accepted component: the alt nodes to delete and the links that replace them.
struct TandemEdit {
    int32_t sn = -1;
    std::vector<NodeId> remove;        // input alt nodes of the site, sorted
    std::vector<std::pair<int64_t, int64_t>> back;   // (from, to), from > to: R side of the reference ending at from -> L side of the one starting at to
    std::vector<std::pair<int64_t, int64_t>> skip;   // (p1, p2), p1 < p2: R side ending at p1 -> L side starting at p2
    int32_t sr = -1;                   // SR of the new links
};

// Steps 1-4 for one site (needs the site's plan and the aligner).  Rewrites the rows of its units
// and appends the site's edits.
void tandem_site(const Graph& g, const SiteData& sd, const SitePlan& plan, Aligner& aligner, const Options& opt,
                 std::vector<ReportRow>& unit_rows, std::vector<TandemEdit>& edits, TandemStats& st);

// The edit, after apply_plans: edits of every site, in site order.
void apply_tandem(Graph& g, const std::vector<const std::vector<TandemEdit>*>& edits, const Options& opt, TandemStats& st);

// Global identity of two short sequences: 1 - edit distance / max length (unit costs).  Exposed for
// tests; the caller keeps the product of the lengths small.
double global_identity(const std::string& a, const std::string& b);
// The same within a diagonal band of half-width w (0 when the lengths differ by more than w).
double banded_identity(const std::string& a, const std::string& b, int64_t w);

} // namespace zip
