/*
  zip_gaf.hpp -- GafWalks (v3): observed haplotype walks as the covering set, and the GAF audit
  (the spec's release gate).

  Source: a pre-zip graphmap GAF of the assemblies against THIS (unzipped) graph.  Paths may be in
  minigraph's stable coordinates (">SN:lo-hi<SN2:lo-hi", one step per run of consecutive nodes of
  one SN), in segment names (">s12<s34"), or a bare stable sequence name (a one-SN path whose
  alignment is columns 8-9).  SN names are matched exactly, then in canonical PanSN form (GAF
  "id=HG00320.1|JBHIKI010000013.1" = graph "HG00320#1#JBHIKI010000013.1", and back), then by their
  last '#' field.  Each record is projected onto handles of this graph exactly as the design
  prototype did: every step must start and end on node boundaries, and (beyond the prototype)
  every step-to-step move must be a link of the graph and the path length must equal column 7.
  Secondary records (tp:A:S) and unmapped ones are skipped; a record none of whose names is in the
  graph belongs to another graph and is skipped; any other record that cannot be projected is
  counted, its nodes are marked partial, and more than 1% of them (or no record touching the graph
  at all) is ZipError(EXIT_INPUT).  The GAF is read once, through `bgzip -@ t -dc` when the file is
  BGZF and bgzip is on PATH, else `gzip -dc`; only the 12 mandatory fields and the first tags of a
  line are kept, so multi-GB cg/ds tags cost only the decompression.

  Excursions: each projected walk is cut at rank-0 handles into (dep, alt handles, arr), put in
  the canonical frame (zip_site's canonicalize) and deduplicated across the whole GAF.  An
  excursion's owners are the contigs that walk it (Owner: contig = canonical name, sr = rank of the
  graph SN of that name or -1 if the contig inserted nothing, so = smallest query start of its
  records that walk the excursion, src GAF); its weight is the number of distinct haplotypes
  (sample#hap) among them.  Alt handles before the first or after the last rank-0 handle of a record
  (a record that ends inside an allele) are "partial": the haplotype touches those nodes but its
  excursion window is unknown.

  complete() is true:
    - units are the observed excursions of the site (both anchors on the site's reference nodes);
      an observed excursion that equals a creator excursion also carries that run's creator owner,
      and its source is then "creator" (the best source among its owners), otherwise "gaf";
      runs no observed walk confirms give no unit;
    - acceptance orders by support (weight) first (the edit module reads SiteData::complete_walks);
    - Allowed(n), for an alt node that observed excursions of the site cross, is the intersection
      of their windows: BLOCKED if one is a junction (J), EMPTY if one is an insertion or back-link
      (I, BK) or the intersection is empty, else OK (sign '-' only when every observed walk crosses
      n as '-').  A node that a partial walk touches, or that an unprojectable record may touch,
      keeps the strict all-path value (conservative), as does a node no observed walk touches.
      Observed walks are graph paths, so this never narrows strict, only widens or opens it.
    - WalkStats: runs = the site's runs; creator_ok = runs whose creator excursion is observed.

  Diagnostics (stderr, after a successful run): Allowed versus strict (same / widened / opened /
  kept strict and what that costs), the creator-versus-observed comparison of the design
  prototype (creator excursion equal to the creator contig's observed excursion, exactly and by
  anchors; observed excursions that are not creator excursions, by count and by contigs).
  --dump DIR adds gaf_excursions.tsv (every distinct observed excursion, its kind, window, weight
  and contigs) and gaf_contigs.tsv (projected records per contig).

  The audit (GafAudit) checks zipped pieces against the observed walks: 0 bp on a node that an
  observed excursion of kind I, BK or J crosses or whose observed window excludes the target, and
  0 GAF lines that walk the node and also read the target (overlap > read_tol).  It is the logic
  of the design prototype's GAF audit, on the same projection.

  Acceptance rule (G) (reads_target, v3): the audit's reading rule, applied when a piece is
  accepted.  Observed Allowed sees one excursion at a time, so a record that walks the node inside
  an ordinary excursion over the window and elsewhere reads the target (an inverted duplication
  on that haplotype) is invisible to it; (G) refuses such a piece.  It uses the index the audit
  uses: per record its merged rank-0 coverage, per alt node the records that walk it; a query is
  one binary search per record through the node.

  Thread-safety: the GAF is read in the constructor; every const method may be called from any
  number of site threads at once.
*/
#pragma once

#include "zip_site.hpp"

namespace zip {

// The projected GAF (defined in zip_gaf.cpp); shared by GafWalks and GafAudit.
class GafIndex;

// The reading rule's tolerances (the prototype audit's), shared by the audit and by rule (G): a record
// reads a piece's target when its reference coverage overlaps the target by more than
// GAF_READ_TOL bp; pieces shorter than GAF_MIN_READ_PIECE bp are not checked.
const int64_t GAF_READ_TOL = 200;
const int64_t GAF_MIN_READ_PIECE = 200;

// Read and project a GAF.  Malformed lines, a decompressor failure, or a GAF that does not match
// the graph (more than 1% of the records that touch it cannot be projected, or none touches it)
// are ZipError(EXIT_INPUT).  With opt.dump_dir set, gaf_excursions.tsv and gaf_contigs.tsv are
// written there as temporary files (gaf_dump_writers; main commits them with its outputs).
std::shared_ptr<const GafIndex> load_gaf(const Graph& g, const std::string& path, const Options& opt);

class GafWalks : public WalkSource {
public:
    // Reads the whole GAF (plain or gzip).  A malformed file is ZipError(EXIT_INPUT).
    GafWalks(const Graph& g, const Runs& runs, const std::string& path, const Options& opt);
    // Uses an index that is already loaded.
    GafWalks(const Graph& g, const Runs& runs, std::shared_ptr<const GafIndex> index, const Options& opt);
    ~GafWalks() override;     // logs the observed-walk summary (not when unwinding from an error)
    GafWalks(const GafWalks&) = delete;
    GafWalks& operator=(const GafWalks&) = delete;

    const char* name() const override { return "gaf"; }
    std::vector<Excursion> excursions(const SiteData& sd, WalkStats& st) const override;
    bool complete() const override;
    void refine_allowed(const SiteData& sd, std::vector<AllowedIv>& allowed) const override;
    // rule (G): a projected record that walks n also reads [ta, tb) of the site's SN by more than
    // GAF_READ_TOL (pieces of at least GAF_MIN_READ_PIECE bp); exactly the audit's reading rule
    bool reads_target(const Site& s, NodeId n, int64_t piece_bp, int64_t ta, int64_t tb) const override;

    std::shared_ptr<const GafIndex> index() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// The finished --dump files of a GafWalks' index (gaf_excursions.tsv, gaf_contigs.tsv): main
// commits them with the other outputs, all or nothing; uncommitted, they are removed.
std::vector<AtomicWriter*> gaf_dump_writers(const GafWalks& w);

// ---------------------------------------------------------------- the GAF audit (release gate)

struct ZippedPiece;   // zip_edit.hpp: SitePlan::zipped

// One zipped piece: the part [a, b) (node-forward) of alt node `node` that lands on target
// [ta, tb), absolute coordinates on the site's SN.  `unit` groups pieces into chains for the
// per-chain verdict lines.
struct GafAuditPiece {
    uint32_t unit = NONE;
    NodeId node = NONE;
    int64_t a = 0, b = 0;
    int64_t ta = 0, tb = 0;
};

struct GafAuditStats {
    uint64_t units = 0, pieces = 0;
    int64_t bp = 0;
    // window rule: an observed excursion of the site through the node is I/BK/J, or its window
    // (widened by window_tol on each side) does not contain the target
    uint64_t viol_units = 0, viol_pieces = 0;
    int64_t viol_bp = 0, ok_bp = 0, unobserved_bp = 0;
    uint32_t max_bad_haps = 0;
    // reading rule: a GAF line that walks the node also reads reference overlapping the target by
    // more than read_tol (pieces >= min_read_piece bp)
    uint64_t read_units = 0, read_lines = 0;
    int64_t read_bp = 0;
    // diagnostic only: another line of a contig that walks the node reads the target
    uint64_t ctg_read_units = 0;
    int64_t ctg_read_bp = 0;
    // alt-vs-alt pieces (targets are branch offsets): counted, not audited
    uint64_t alt_pieces = 0;
    int64_t alt_bp = 0;
    void add(const GafAuditStats& o);
    // the release gate: 0 bp violating the window rule and 0 lines reading the target
    bool clean() const { return viol_bp == 0 && read_lines == 0; }
    std::string summary() const;
};

class GafAudit {
public:
    // Loads its own index from `path`.
    GafAudit(const Graph& g, const std::string& path, const Options& opt);
    // Shares an index (e.g. GafWalks::index()).
    GafAudit(const Graph& g, std::shared_ptr<const GafIndex> index);
    ~GafAudit();
    GafAudit(const GafAudit&) = delete;
    GafAudit& operator=(const GafAudit&) = delete;

    // tolerances (the prototype audit used 200 for both); defaults: window 0, read 200, pieces >= 200 bp
    // (the read tolerances are rule (G)'s too)
    int64_t window_tol = 0, read_tol = GAF_READ_TOL, min_read_piece = GAF_MIN_READ_PIECE;

    // Audit the zipped pieces of one site.  `detail` (optional) receives one line per unit:
    //   GAFAUDIT <site> <unit> pieces=<n> bp=<bp> obs=<harm|safe|unobserved> viol_bp= ok_bp= none_bp=
    //   max_bad_haps= read_lines= read_ctgs=
    // followed by one "GAFAUDIT-PIECE ..." line per violating piece.
    GafAuditStats audit_site(const Site& s, const std::vector<GafAuditPiece>& pieces, std::string* detail) const;
    // The same for a site's plan (SitePlan::zipped of zip_edit): reference-pass pieces are audited,
    // alt-vs-alt pieces only counted.
    GafAuditStats audit_plan(const Site& s, const std::vector<ZippedPiece>& zipped, std::string* detail) const;

private:
    const Graph& g_;
    std::shared_ptr<const GafIndex> gi_;
};

} // namespace zip
