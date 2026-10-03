/*
  zip_report.hpp -- the report (TSV): one row per candidate unit, one row per skipped or capped
  site, written once after every decision, plus the --stage detect printout and --dump tables.

  Who fills which ReportRow fields:
    detect (zip_site / unit_row):  site, sn, lo, hi, pass, round, source, owners, kind, anchors,
                                   window, query_bp, feasible_bp, outcome
    align  (zip_align):            records, blocks, frame, label, tie, alt_score, alt_E, internal,
                                   repeat_frac, outcome (no-chain, below-b, ambiguous-strand, ...)
    edit   (zip_edit):             kept_bp, dropped, outcome (zipped, (D), trimmed-below-b, reverted:<check>)
    alt    (zip_alt, v2):          pass "alt", round, anchors = bounds "u>v", representative
  Unset fields print as ".".  Rows are written in site order, then in each site's row order, so the
  report is byte-identical across -t/-j.  Peak RSS and timings go to stderr, never here.
*/
#pragma once

#include "zip_site.hpp"

namespace zip {

struct ReportRow {
    // ---- identity (detect)
    std::string site;              // "s<A>..s<B>"
    std::string sn;                // reference SN of the site
    int64_t lo = -1, hi = -1;      // the site's window [lo, hi)
    std::string pass;              // "site" (skip/cap rows), "ref", "alt", "diag"
    int round = 0;                 // 0 for the reference pass; 1..--alt-rounds for alt-vs-alt
    std::string source;            // creator | witness | gaf
    std::string owners;            // owner runs "<contig>:<first node>", comma-separated
    std::string kind;              // F | I | BK | J
    std::string anchors;           // "dep>arr" handles, or alt-vs-alt bounds "u>v"
    std::string window;            // "lo-hi" (absolute)
    int64_t query_bp = -1;         // unit query bp; for site rows, the site's alt bp
    int64_t feasible_bp = -1;      // feasible query bp (Allowed OK inside the window)
    // ---- alignment (align module)
    int64_t records = -1;          // verified PAF records
    std::string blocks;            // chain blocks "qs-qe:ts-te:strand:identity", ';'-separated
    std::string frame;             // F | R
    std::string label;             // FWD | INV | the strand string (e.g. "+-+"), in walk orientation
    std::string tie;               // unique | tie-kept | tie-resolved | ambiguous-strand
    int64_t alt_score = -1;        // best alternative chain's score
    int64_t alt_E = -1;            // ... and its E (out-of-place stretches)
    int64_t internal = -1;         // internal unaligned stretches >= G
    double repeat_frac = -1;       // repeat fraction (diagnostic)
    // ---- acceptance (edit module)
    int64_t kept_bp = -1;          // alt bp this chain removes from the graph: its zipped pieces that no
                                   // earlier chain zipped (a node interval several compatible chains zip
                                   // counts once, for the first); 0 for a chain that zips nothing (refused,
                                   // reverted, or wholly compatible), so the column sums to the bp removed
    std::string dropped;           // dropped pieces with their rules, e.g. "s12[0-500):(C)", ';'-separated
                                   // (rules (A) (B) (C) (G) (empty) (short) (fold)); first the align
                                   // module's island-below-b groups and truncated:500/<n>-sub-records
    // ---- every stage
    std::string outcome;           // see the spec's "What it declines to do"; "candidate" at --stage detect
    std::string representative;    // alt-vs-alt: the representative branch
    // ---- not printed
    uint32_t unit = NONE;          // unit id within its site
};

// column names in output order
const std::vector<std::string>& report_columns();
// one TSV line (no newline)
std::string format_row(const ReportRow& r);

// owners as "<contig>:<first node>" (GAF owners: "<contig>"), comma-separated
std::string owners_str(const Graph& g, const std::vector<Owner>& owners);
// a row for a skipped site (site:leak, site:ref-mismatch, site:too-big); query_bp carries the alt bp
ReportRow site_skip_row(const Graph& g, const Site& s);
// the detect fields of a unit's row (outcome "" becomes "candidate")
ReportRow unit_row(const Graph& g, const SiteData& sd, const Unit& u);

// header comment lines: tool version, minimap2 version, and every option that affects the output
// (numbers exactly as given: two different settings never print alike)
std::vector<std::string> report_header(const Options& opt, const std::string& minimap2_version);
// write header + column line + rows (per site, in order)
void write_report(AtomicWriter& w, const std::vector<std::string>& header, const std::vector<const std::vector<ReportRow>*>& rows);

// --stage detect printout: a header, then per site one SITE line and one UNIT line per unit
std::string detect_header();
std::string detect_site_lines(const Graph& g, const Site& s, const SiteData* sd);

// --dump tables (one string per site; main concatenates them in site order)
std::string dump_units_header();
std::string dump_units(const Graph& g, const SiteData& sd);       // every unit, with its walk
std::string dump_allowed_header();
std::string dump_allowed(const Graph& g, const SiteData& sd);     // Allowed of every alt node

} // namespace zip
