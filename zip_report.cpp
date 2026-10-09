/*
  zip_report.cpp -- report rows, header and writer; --stage detect printout; --dump tables.
*/
#include "zip_report.hpp"

namespace zip {

namespace {

inline void put_str(std::string& out, const std::string& s) {
    if (s.empty()) { out.push_back('.'); return; }
    for (char c : s) out.push_back(c == '\t' || c == '\n' ? ' ' : c);
}
inline void put_int(std::string& out, int64_t v) {
    if (v < 0) out.push_back('.');
    else out += std::to_string(v);
}

} // namespace

const std::vector<std::string>& report_columns() {
    static const std::vector<std::string> cols = {
        "site", "sn", "lo", "hi", "pass", "round", "source", "owners", "kind", "anchors", "window",
        "query_bp", "feasible_bp", "records", "blocks", "frame", "label", "tie", "alt_score", "alt_E",
        "internal", "repeat_frac", "kept_bp", "dropped", "outcome", "representative"};
    return cols;
}

std::string format_row(const ReportRow& r) {
    std::string o;
    o.reserve(256);
    put_str(o, r.site); o.push_back('\t');
    put_str(o, r.sn); o.push_back('\t');
    put_int(o, r.lo); o.push_back('\t');
    put_int(o, r.hi); o.push_back('\t');
    put_str(o, r.pass); o.push_back('\t');
    o += std::to_string(r.round); o.push_back('\t');
    put_str(o, r.source); o.push_back('\t');
    put_str(o, r.owners); o.push_back('\t');
    put_str(o, r.kind); o.push_back('\t');
    put_str(o, r.anchors); o.push_back('\t');
    put_str(o, r.window); o.push_back('\t');
    put_int(o, r.query_bp); o.push_back('\t');
    put_int(o, r.feasible_bp); o.push_back('\t');
    put_int(o, r.records); o.push_back('\t');
    put_str(o, r.blocks); o.push_back('\t');
    put_str(o, r.frame); o.push_back('\t');
    put_str(o, r.label); o.push_back('\t');
    put_str(o, r.tie); o.push_back('\t');
    put_int(o, r.alt_score); o.push_back('\t');
    put_int(o, r.alt_E); o.push_back('\t');
    put_int(o, r.internal); o.push_back('\t');
    if (r.repeat_frac < 0) o.push_back('.'); else o += fmt_fixed(r.repeat_frac, 3);
    o.push_back('\t');
    put_int(o, r.kept_bp); o.push_back('\t');
    put_str(o, r.dropped); o.push_back('\t');
    put_str(o, r.outcome); o.push_back('\t');
    put_str(o, r.representative);
    return o;
}

std::string owners_str(const Graph& g, const std::vector<Owner>& owners) {
    std::string s;
    for (size_t i = 0; i < owners.size(); ++i) {
        if (i) s.push_back(',');
        s += owners[i].contig;
        if (owners[i].first != NONE) { s.push_back(':'); s += g.name(owners[i].first); }
    }
    return s;
}

ReportRow site_skip_row(const Graph& g, const Site& s) {
    ReportRow r;
    r.site = s.label(g);
    r.sn = g.sn_names[s.sn];
    r.lo = s.lo;
    r.hi = s.hi;
    r.pass = "site";
    r.window = strf("%lld-%lld", (long long)s.lo, (long long)s.hi);
    r.query_bp = s.alt_bp;
    r.outcome = site_status_name(s.status);
    return r;
}

ReportRow unit_row(const Graph& g, const SiteData& sd, const Unit& u) {
    const Site& s = *sd.site;
    ReportRow r;
    r.site = s.label(g);
    r.sn = g.sn_names[s.sn];
    r.lo = s.lo;
    r.hi = s.hi;
    r.pass = "ref";
    r.round = 0;
    r.source = src_name(u.exc.src);
    r.owners = owners_str(g, u.exc.owners);
    r.kind = kind_name(u.kind);
    r.anchors = g.handle_str(u.exc.dep) + ">" + g.handle_str(u.exc.arr);
    if (u.kind != Kind::J) r.window = strf("%lld-%lld", (long long)u.wlo, (long long)u.whi);
    r.query_bp = u.query_bp;
    if (u.kind == Kind::F) r.feasible_bp = u.feasible_bp;
    r.outcome = u.outcome.empty() ? "candidate" : u.outcome;
    r.unit = u.id;
    return r;
}

std::vector<std::string> report_header(const Options& opt, const std::string& minimap2_version) {
    std::vector<std::string> h;
    h.push_back(strf("#rgfa-zip\t%s", ZIP_VERSION));
    h.push_back("#minimap2\t" + (minimap2_version.empty() ? std::string("unknown") : minimap2_version));
    std::string walks = opt.walks == WalksKind::CREATOR ? "creator" : opt.walks == WalksKind::WITNESS ? "witness" : "gaf:" + opt.gaf_path;
    std::string o = strf("-b %lld -i %s -G %lld --min-piece %lld --island %lld --delta %s --ties %s --frag %lld --prefilter %lld,%s "
                         "--screen-sample %d --max-pair %lld --max-site-query %lld --max-site-nodes %lld -x %s",
                         (long long)opt.b, fmt_exact(opt.ident).c_str(), (long long)opt.G, (long long)opt.min_piece,
                         (long long)opt.island, fmt_exact(opt.delta).c_str(), opt.ties_refuse ? "refuse" : "resolve",
                         (long long)opt.frag, (long long)opt.prefilter_nmatch, fmt_exact(opt.prefilter_cov).c_str(),
                         opt.screen_sample, (long long)opt.max_pair, (long long)opt.max_site_query, (long long)opt.max_site_nodes,
                         opt.preset.c_str());
    if (!opt.mm_extra.empty()) o += " -X '" + opt.mm_extra + "'";
    o += opt.alt ? strf(" --alt-rounds %d", opt.alt_rounds) : std::string(" --no-alt");
    o += " --walks " + walks;
    if (opt.max_repeat_frac >= 0) o += " --max-repeat-frac " + fmt_exact(opt.max_repeat_frac);
    if (opt.detect_only) o += " --detect-only";
    if (opt.check) o += " --check";
    if (opt.strict) o += " --strict";
    if (opt.id_base >= 0) o += strf(" --id-base %lld", (long long)opt.id_base);
    for (const Region& rg : opt.regions) o += strf(" --region %s:%lld-%lld", rg.sn.c_str(), (long long)rg.lo, (long long)rg.hi);
    static const char* stages[] = {"detect", "align", "plan", "full"};
    if (opt.stage != Stage::FULL) o += std::string(" --stage ") + stages[(int)opt.stage];
    if (!opt.inject_chains.empty()) o += " --inject-chains " + opt.inject_chains;
    h.push_back("#options\t" + o);
    return h;
}

void write_report(AtomicWriter& w, const std::vector<std::string>& header, const std::vector<const std::vector<ReportRow>*>& rows) {
    for (const std::string& l : header) { w.write(l); w.put('\n'); }
    const std::vector<std::string>& cols = report_columns();
    std::string line = "#";
    for (size_t i = 0; i < cols.size(); ++i) {
        if (i) line.push_back('\t');
        line += cols[i];
    }
    line.push_back('\n');
    w.write(line);
    for (const std::vector<ReportRow>* v : rows) {
        if (!v) continue;
        for (const ReportRow& r : *v) {
            w.write(format_row(r));
            w.put('\n');
        }
    }
}

// ---------------------------------------------------------------- --stage detect

std::string detect_header() {
    return "#SITE\tsite\tsn\tlo\thi\tstatus\tinterior\talts\talt_bp\truns\tcreator_ok\tfallback\twitness\tunanchored\t"
           "excursions\tF\tI\tBK\tJ\tref_units\tto_align\talign_query_bp\talign_feasible_bp\tinfeasible\tpair_too_big\t"
           "capped\tno_window\tallowed_ok\tallowed_empty\tallowed_blocked\tallowed_none\n"
           "#UNIT\tsite\tunit\towners\tsource\tkind\tanchors\twindow\twindow_bp\tquery_bp\tfeasible_bp\tnodes\toutcome\n";
}

std::string detect_site_lines(const Graph& g, const Site& s, const SiteData* sd) {
    std::string o;
    std::string lab = s.label(g);
    uint64_t kinds[4] = {0, 0, 0, 0};
    uint64_t ref_units = 0, to_align = 0, infeasible = 0, too_big = 0, capped = 0, no_window = 0;
    int64_t aq = 0, af = 0;
    uint64_t al[4] = {0, 0, 0, 0};
    if (sd) {
        for (const Unit& u : sd->units) {
            kinds[(int)u.kind]++;
            if (u.ref_unit) ++ref_units;
            if (u.ref_unit && u.outcome.empty()) { ++to_align; aq += u.query_bp; af += u.feasible_bp; }
            if (u.outcome == "infeasible") ++infeasible;
            else if (u.outcome == "pair-too-big") ++too_big;
            else if (u.outcome == "site:capped") ++capped;
            else if (u.outcome == "no-window") ++no_window;
        }
        for (const AllowedIv& a : sd->allowed) al[(int)a.status]++;
    }
    const WalkStats ws = sd ? sd->wstats : WalkStats();
    o += strf("SITE\t%s\t%s\t%lld\t%lld\t%s\t%zu\t%zu\t%lld\t%llu\t%llu\t%llu\t%llu\t%llu\t%zu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t"
              "%lld\t%lld\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\t%llu\n",
              lab.c_str(), g.sn_names[s.sn].c_str(), (long long)s.lo, (long long)s.hi, site_status_name(s.status), s.n_interior,
              s.alts.size(), (long long)s.alt_bp, (unsigned long long)ws.runs, (unsigned long long)ws.creator_ok,
              (unsigned long long)ws.fallback, (unsigned long long)(ws.witness_dag + ws.witness_bubble),
              (unsigned long long)ws.unanchored, sd ? sd->units.size() : (size_t)0, (unsigned long long)kinds[0],
              (unsigned long long)kinds[1], (unsigned long long)kinds[2], (unsigned long long)kinds[3], (unsigned long long)ref_units,
              (unsigned long long)to_align, (long long)aq, (long long)af, (unsigned long long)infeasible, (unsigned long long)too_big,
              (unsigned long long)capped, (unsigned long long)no_window, (unsigned long long)al[(int)AllowedStatus::OK],
              (unsigned long long)al[(int)AllowedStatus::EMPTY], (unsigned long long)al[(int)AllowedStatus::BLOCKED],
              (unsigned long long)al[(int)AllowedStatus::NONE]);
    if (!sd) return o;
    for (const Unit& u : sd->units) {
        std::string win = u.kind == Kind::J ? std::string(".") : strf("%lld-%lld", (long long)u.wlo, (long long)u.whi);
        o += strf("UNIT\t%s\t%u\t%s\t%s\t%s\t%s>%s\t%s\t%lld\t%lld\t%lld\t%zu\t%s\n", lab.c_str(), u.id,
                  owners_str(g, u.exc.owners).c_str(), src_name(u.exc.src), kind_name(u.kind), g.handle_str(u.exc.dep).c_str(),
                  g.handle_str(u.exc.arr).c_str(), win.c_str(), (long long)(u.kind == Kind::J ? 0 : u.whi - u.wlo),
                  (long long)u.query_bp, (long long)(u.kind == Kind::F ? u.feasible_bp : -1), u.exc.alts.size(),
                  u.outcome.empty() ? "candidate" : u.outcome.c_str());
    }
    return o;
}

// ---------------------------------------------------------------- --dump

std::string dump_units_header() {
    return "site\tunit\tkind\tsource\tdep\tarr\twlo\twhi\tquery_bp\tfeasible_bp\toutcome\tweight\towners\twalk\towner_sr\n";
}

std::string dump_units(const Graph& g, const SiteData& sd) {
    std::string o;
    std::string lab = sd.site->label(g);
    for (const Unit& u : sd.units) {
        std::string walk;
        for (size_t i = 0; i < u.exc.alts.size(); ++i) {
            if (i) walk.push_back(',');
            walk += g.handle_str(u.exc.alts[i]);
        }
        o += strf("%s\t%u\t%s\t%s\t%s\t%s\t%lld\t%lld\t%lld\t%lld\t%s\t%u\t%s\t%s\t%d\n", lab.c_str(), u.id, kind_name(u.kind),
                  src_name(u.exc.src), g.handle_str(u.exc.dep).c_str(), g.handle_str(u.exc.arr).c_str(), (long long)u.wlo,
                  (long long)u.whi, (long long)u.query_bp, (long long)u.feasible_bp, u.outcome.empty() ? "candidate" : u.outcome.c_str(),
                  u.exc.weight, owners_str(g, u.exc.owners).c_str(), walk.c_str(), u.exc.owners.empty() ? -1 : (int)u.exc.owners[0].sr);
    }
    return o;
}

std::string dump_allowed_header() { return "site\tnode\tlen\tstatus\tsign\tlo\thi\n"; }

std::string dump_allowed(const Graph& g, const SiteData& sd) {
    std::string o;
    std::string lab = sd.site->label(g);
    for (size_t l = 0; l < sd.sg.alt.size(); ++l) {
        const AllowedIv& a = sd.allowed[l];
        bool iv = a.status == AllowedStatus::OK || a.status == AllowedStatus::EMPTY;
        o += strf("%s\t%s\t%lld\t%s\t%s\t%s\t%s\n", lab.c_str(), g.name(sd.sg.alt[l]).c_str(), (long long)g.len(sd.sg.alt[l]),
                  allowed_status_name(a.status), iv ? (a.minus ? "-" : "+") : ".", iv ? std::to_string(a.lo).c_str() : ".",
                  iv ? std::to_string(a.hi).c_str() : ".");
    }
    return o;
}

} // namespace zip
