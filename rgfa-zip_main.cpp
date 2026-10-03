/*
  rgfa-zip: zip minigraph alleles onto the reference they bypass, and onto parallel alleles.

  For each top-level vg snarl with reference boundaries, take one creator excursion per minigraph
  insertion run as the covering set of traversals, align each excursion's alt sequence (both
  strands) to the reference window it bypasses, and zip it only when the homology is confident and
  the zip is safe for every graph path (Allowed).  All decisions are made on the input graph and
  applied once; validators re-spell every excursion before a change commits; output is
  deterministic and written atomically.

  Versions by option: v1 = reference pass only (--no-alt); v2 = + alt-vs-alt (default);
  v3 = observed walks (--walks gaf:FILE), which also gives acceptance the whole-record rule (G).

  Logs: every per-window and per-site time separates minimap2 running from waiting for a -j slot,
  and gives minimap2's own real and CPU time (AlignTiming), so the log ranks real costs.

  This file: command line, orchestration, the per-site thread pool, outputs, the
  split-inversion-candidate diagnostic and the --audit-gaf wiring.
  Per site: prepare_site (walks, units, Allowed) -> Aligner::align_reference -> diagnostics ->
  SitePlanner::accept_reference -> [alt_pass, v2] -> SitePlanner::finish (transaction, V1-V7) ->
  GafAudit::audit_plan.  Then apply_plans (ids), check_placement, and the atomic writes.
    zip_graph   rGFA model, input checks, emit, global placement assert
    zip_site    snarls, sites, runs, walks, units, kinds, Allowed, pre-filter, caps
    zip_align   aligner, records, blocks, feasibility, chain, island rule, U, fragmentation
    zip_edit    acceptance, cuts, splitter, pinch, transaction, validators, ids
    zip_alt     alt-vs-alt (v2): residual excursions, branches, groups, eligibility, rounds
    zip_gaf     GafWalks, the GAF audit
    zip_report  report rows, detect printout, dumps
*/
#include "zip_common.hpp"
#include "zip_graph.hpp"
#include "zip_site.hpp"
#include "zip_report.hpp"
#include "zip_align.hpp"
#include "zip_edit.hpp"
#include "zip_alt.hpp"
#include "zip_gaf.hpp"

#include <chrono>
#include <deque>
#include <iostream>
#include <map>
#include <set>

using namespace zip;

namespace {

// ================================================================ command line

const char* const USAGE =
    "usage: rgfa-zip -m <minimap2> [options] <in.gfa> <snarls.json> -o out.gfa -r report.tsv\n"
    "\n"
    "Zip minigraph alleles onto the reference they bypass, and onto parallel alleles.\n"
    "snarls.json from:  vg snarls -n -P <ref> <in.gfa> | vg view -Rj -\n"
    "Inputs may be gzipped.\n"
    "\n"
    "Required:\n"
    "  -m PATH                  minimap2 binary; there is no PATH default (its version goes in the report)\n"
    "  -o FILE                  output rGFA (FILE.tmp, checked writes, fsync, rename)\n"
    "  -r FILE                  report TSV: one row per candidate, and per skipped or capped site\n"
    "Parallelism and resources:\n"
    "  -t INT                   sites in parallel [1]\n"
    "  -j INT                   concurrent aligner processes [= -t, lowered so j x 1.5 GB fits --mem]\n"
    "  --mem SIZE               memory for aligner processes, e.g. 16G [all]\n"
    "  --tmpdir DIR             temporary files [$TMPDIR]\n"
    "Thresholds:\n"
    "  -b INT                   min chain bp, and min unit query [5000]\n"
    "  -i FLOAT                 min gap-compressed identity [0.95]\n"
    "  -G INT                   I/D runs >= G split records into blocks [50]\n"
    "  --min-piece INT          min aligned query per record [1000]\n"
    "  --island INT             gap that separates island groups [1000]\n"
    "  --delta FLOAT            tie tolerance of the unique-placement rule [0.005]\n"
    "  --ties resolve|refuse    resolve array ties positionally, or refuse them [resolve]\n"
    "  --frag INT               max internal unaligned stretches per 5 kb, floor 2 [1]\n"
    "  --prefilter N,F          screen: records' nmatch sum >= N and one with nmatch/blocklen >= F [2000,0.2]\n"
    "  --screen-sample INT      queries screened first per window [8]\n"
    "  --max-pair INT           max query and max window per pair [5000000]\n"
    "  --max-site-query INT     max feasible query aligned per site [50000000]\n"
    "  --max-site-nodes INT     max interior nodes per site [200000]\n"
    "Aligner:\n"
    "  -x STR                   minimap2 preset [asm20]\n"
    "  -X STR                   extra minimap2 arguments, appended last [none]\n"
    "Passes and walks:\n"
    "  --no-alt                 reference pass only (v1); alt-vs-alt is on by default (v2)\n"
    "  --alt-rounds INT         alt-vs-alt rounds [3]\n"
    "  --walks SRC              covering walks: creator | witness | gaf:FILE (v3) [creator]; with gaf:FILE a\n"
    "                           piece is also refused (G) when a GAF record through its node reads its target\n"
    "  --max-repeat-frac FLOAT  skip windows more repetitive than FLOAT, with a report row [off]\n"
    "Modes:\n"
    "  --detect-only            decide and report everything; emit the input graph unchanged\n"
    "  --check                  validators also re-spell every anchored path at sites with <= 5000 paths\n"
    "  --strict                 exit 3 if the transaction drops any chain\n"
    "  --id-base INT            first id of new pieces [max input id + 1]\n"
    "  --audit-gaf FILE         audit the zips against the observed walks of a pre-zip GAF of this graph (the\n"
    "                           release gate: 0 bp contradicting a walk, 0 GAF lines reading the target);\n"
    "                           summary on stderr, per-chain verdicts in --dump DIR/gafaudit.tsv\n"
    "Debug:\n"
    "  --region SN:LO-HI        only sites overlapping this region (SN or its last #-field; repeatable)\n"
    "  --stage S                stop after detect | align | plan | full [full]\n"
    "  --inject-chains FILE     take chains from FILE instead of running minimap2\n"
    "  --dump DIR               write units, Allowed, records and chains under DIR\n"
    "  -v                       more log output (repeatable)\n"
    "  -h, --help               this help\n"
    "\n"
    "Exit codes: 0 ok; 2 invalid input; 3 invariant failure; 4 systemic aligner failure; 5 I/O.\n"
    "Nothing is written on a non-zero exit.\n";

[[noreturn]] void usage_error(const std::string& msg) { fail(EXIT_INPUT, msg + "\n(rgfa-zip -h lists the options)"); }

// set once every output has been renamed into place
std::atomic<bool> g_committed(false);

// a path made comparable: the realpath of the file, else of its directory plus the name
std::string canon_path(const std::string& p) {
    char buf[PATH_MAX];
    if (realpath(p.c_str(), buf)) return buf;
    size_t sl = p.rfind('/');
    std::string dir = sl == std::string::npos ? "." : (sl == 0 ? "/" : p.substr(0, sl));
    std::string name = sl == std::string::npos ? p : p.substr(sl + 1);
    if (realpath(dir.c_str(), buf)) return std::string(buf) + "/" + name;
    return p;
}

// The outputs (-o, -r, their .tmp files, the --dump files) may not be one another, an input, or a
// directory: an alias would let one write clobber another, and a directory would only fail at the
// rename, after other outputs were already in place.
void check_output_paths(const Options& opt) {
    std::vector<std::pair<std::string, std::string>> inputs = {{"the input graph", opt.gfa_path}, {"the snarls", opt.snarls_path}};
    if (!opt.gaf_path.empty()) inputs.push_back(std::make_pair("--walks gaf:", opt.gaf_path));
    if (!opt.audit_gaf.empty()) inputs.push_back(std::make_pair("--audit-gaf", opt.audit_gaf));
    if (!opt.inject_chains.empty()) inputs.push_back(std::make_pair("--inject-chains", opt.inject_chains));
    std::vector<std::pair<std::string, std::string>> outs;   // (what, path)
    if (!opt.out_path.empty()) {
        outs.push_back(std::make_pair("-o", opt.out_path));
        outs.push_back(std::make_pair("-o's temporary file", opt.out_path + ".tmp"));
    }
    if (!opt.report_path.empty()) {
        outs.push_back(std::make_pair("-r", opt.report_path));
        outs.push_back(std::make_pair("-r's temporary file", opt.report_path + ".tmp"));
    }
    std::string dump = opt.dump_dir.empty() ? std::string() : canon_path(opt.dump_dir);
    if (!dump.empty()) {
        for (const char* f : {"units.tsv", "allowed.tsv", "gafaudit.tsv", "gaf_excursions.tsv", "gaf_contigs.tsv"}) {
            outs.push_back(std::make_pair(strf("--dump's %s", f), opt.dump_dir + "/" + f));
            outs.push_back(std::make_pair(strf("--dump's %s.tmp", f), opt.dump_dir + "/" + f + ".tmp"));
        }
    }
    std::vector<std::string> canon;
    for (const auto& o : outs) canon.push_back(canon_path(o.second));
    for (size_t i = 0; i < outs.size(); ++i) {
        for (size_t k = 0; k < i; ++k)
            if (canon[i] == canon[k]) usage_error(strf("%s (%s) is also %s", outs[i].first.c_str(), outs[i].second.c_str(), outs[k].first.c_str()));
        for (const auto& in : inputs)
            if (canon[i] == canon_path(in.second)) usage_error(strf("%s (%s) is %s", outs[i].first.c_str(), outs[i].second.c_str(), in.first.c_str()));
        if (!dump.empty() && outs[i].first.compare(0, 6, "--dump") != 0)
            for (const char* sub : {"/align/", "/edit/", "/alt/"})
                if (canon[i].compare(0, dump.size() + strlen(sub), dump + sub) == 0)
                    usage_error(strf("%s (%s) lies in --dump's %s directory", outs[i].first.c_str(), outs[i].second.c_str(), sub));
        struct stat st;
        if (stat(outs[i].second.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
            usage_error(strf("%s (%s) is a directory", outs[i].first.c_str(), outs[i].second.c_str()));
    }
}

// rgfa-collapse options: errors that name the rgfa-zip equivalent
const std::map<std::string, std::string>& old_options() {
    static const std::string gone_alt = "is gone: queries are alt walks only, so there is no reference flank to discount";
    static const std::string tr = "is not part of rgfa-zip: tandem-repeat flattening stays in the frozen rgfa-collapse -T "
                                  "until it gets its own tool (spec decision 8)";
    static const std::map<std::string, std::string> m = {
        {"-a", "rgfa-collapse's -a " + gone_alt},
        {"--min-alt", "rgfa-collapse's --min-alt " + gone_alt},
        {"-C", "rgfa-collapse's -C is gone: nodes are cut at block ends, so only the aligned part is zipped"},
        {"--min-cover", "rgfa-collapse's --min-cover is gone: nodes are cut at block ends, so only the aligned part is zipped"},
        {"-M", "rgfa-collapse's -M is gone: one collinear chain per walk zips at most min(M,N) copies, and V7 checks removed bp <= target bp / i"},
        {"--max-removed", "rgfa-collapse's --max-removed is gone: one collinear chain per walk zips at most min(M,N) copies"},
        {"-A", "rgfa-collapse's -A N is now --alt-rounds N (alt-vs-alt is on by default; --no-alt turns it off)"},
        {"-D", "rgfa-collapse's -D is gone: strand is never a gate, forward and inverted homology are zipped by the same rules"},
        {"--duplications", "rgfa-collapse's --duplications is gone: strand is never a gate"},
        {"-c", "rgfa-collapse's -c is now --max-site-nodes N (interior nodes per site)"},
        {"--max-components", "rgfa-collapse's --max-components is now --max-site-nodes N (interior nodes per site)"},
        {"-L", "rgfa-collapse's -L is now --max-pair N (max query and window per pair)"},
        {"--max-traversal", "rgfa-collapse's --max-traversal is now --max-pair N"},
        {"-k", "rgfa-collapse's -k is gone: each window gets its own minimap2 index"},
        {"--chunk", "rgfa-collapse's --chunk is gone: each window gets its own minimap2 index"},
        {"--jobs", "use -j N (concurrent aligner processes)"},
        {"--threads", "rgfa-collapse's --threads (threads per minimap2) is gone: minimap2 runs with -t 1; -t N is sites in parallel"},
        {"-N", "rgfa-collapse's -N is gone: rgfa-zip runs minimap2 with -N 50 -p 0.01 --secondary=yes; override with -X"},
        {"--mm-secondary", "rgfa-collapse's --mm-secondary is gone: override minimap2's -N with -X '-N INT'"},
        {"-P", "rgfa-collapse's -P is gone: override minimap2's -p with -X '-p FLOAT'"},
        {"--mm-p", "rgfa-collapse's --mm-p is gone: override minimap2's -p with -X '-p FLOAT'"},
        {"--mm-preset", "use -x STR"},
        {"--mm-extra", "use -X STR"},
        {"--minimap2", "use -m PATH"},
        {"-z", "rgfa-collapse's -z (lastz) is gone: lastz is a test arbiter outside the tool"},
        {"--lastz", "rgfa-collapse's --lastz is gone: lastz is a test arbiter outside the tool"},
        {"-Z", "rgfa-collapse's -Z is gone with lastz"},
        {"--lastz-extra", "rgfa-collapse's --lastz-extra is gone with lastz"},
        {"-R", "rgfa-collapse's -R is gone: -r writes one row per candidate with its outcome"},
        {"--call-report", "rgfa-collapse's --call-report is gone: -r writes one row per candidate with its outcome"},
        {"--report", "use -r FILE"},
        {"-d", "rgfa-collapse's -d is now --detect-only"},
        {"--min-block", "use -b INT"},
        {"--min-ident", "use -i FLOAT"},
        {"-T", "rgfa-collapse's -T " + tr},
        {"--tr-flatten", "--tr-flatten " + tr},
        {"--tr-max-span", "--tr-max-span " + tr},
        {"--tr-max-allele", "--tr-max-allele " + tr},
        {"--tr-min-nodes", "--tr-min-nodes " + tr},
        {"--tr-min-frac", "--tr-min-frac " + tr},
        {"--tr-max-period", "--tr-max-period " + tr},
        {"--tr-report", "--tr-report " + tr},
    };
    return m;
}

enum OptId {
    O_M, O_O, O_R, O_T, O_J, O_B, O_I, O_G, O_MINPIECE, O_ISLAND, O_DELTA, O_TIES, O_FRAG, O_PREFILTER, O_SCREEN,
    O_MAXPAIR, O_MAXSITEQ, O_MAXSITEN, O_X, O_XX, O_NOALT, O_ALTROUNDS, O_WALKS, O_MAXREP, O_DETECTONLY, O_CHECK,
    O_STRICT, O_MEM, O_TMPDIR, O_IDBASE, O_AUDIT, O_REGION, O_STAGE, O_INJECT, O_DUMP, O_VERBOSE, O_HELP
};

struct OptDef {
    const char* name;
    OptId id;
    bool arg;
};

const OptDef OPTS[] = {
    {"-m", O_M, true}, {"-o", O_O, true}, {"-r", O_R, true}, {"-t", O_T, true}, {"-j", O_J, true},
    {"-b", O_B, true}, {"-i", O_I, true}, {"-G", O_G, true}, {"--min-piece", O_MINPIECE, true},
    {"--island", O_ISLAND, true}, {"--delta", O_DELTA, true}, {"--ties", O_TIES, true}, {"--frag", O_FRAG, true},
    {"--prefilter", O_PREFILTER, true}, {"--screen-sample", O_SCREEN, true}, {"--max-pair", O_MAXPAIR, true},
    {"--max-site-query", O_MAXSITEQ, true}, {"--max-site-nodes", O_MAXSITEN, true}, {"-x", O_X, true},
    {"-X", O_XX, true}, {"--no-alt", O_NOALT, false}, {"--alt-rounds", O_ALTROUNDS, true}, {"--walks", O_WALKS, true},
    {"--max-repeat-frac", O_MAXREP, true}, {"--detect-only", O_DETECTONLY, false}, {"--check", O_CHECK, false},
    {"--strict", O_STRICT, false}, {"--mem", O_MEM, true}, {"--tmpdir", O_TMPDIR, true}, {"--id-base", O_IDBASE, true},
    {"--audit-gaf", O_AUDIT, true},
    {"--region", O_REGION, true}, {"--stage", O_STAGE, true}, {"--inject-chains", O_INJECT, true}, {"--dump", O_DUMP, true},
    {"-v", O_VERBOSE, false}, {"-h", O_HELP, false}, {"--help", O_HELP, false},
};

int64_t need_int(const std::string& name, const std::string& v, int64_t lo, int64_t hi) {
    int64_t x = 0;
    if (!parse_int64(v, x) || x < lo || x > hi) usage_error(strf("%s: '%s' is not an integer in [%lld, %lld]", name.c_str(), v.c_str(), (long long)lo, (long long)hi));
    return x;
}
double need_double(const std::string& name, const std::string& v, double lo, double hi) {
    double x = 0;
    if (!parse_double(v, x) || x < lo || x > hi) usage_error(strf("%s: '%s' is not a number in [%g, %g]", name.c_str(), v.c_str(), lo, hi));
    return x;
}

Region parse_region(const std::string& v) {
    size_t c = v.rfind(':');
    if (c == std::string::npos || c == 0) usage_error("--region: expected SN:LO-HI, got '" + v + "'");
    Region r;
    r.sn = v.substr(0, c);
    std::string rng = v.substr(c + 1);
    size_t d = rng.find('-');
    if (d == std::string::npos) usage_error("--region: expected SN:LO-HI, got '" + v + "'");
    if (!parse_int64(rng.substr(0, d), r.lo) || !parse_int64(rng.substr(d + 1), r.hi) || r.lo < 0 || r.hi < r.lo)
        usage_error("--region: bad coordinates in '" + v + "'");
    return r;
}

// returns false when -h was given
bool parse_args(int argc, char** argv, Options& opt) {
    std::vector<std::string> pos;
    for (int i = 0; i < argc; ++i) {
        if (i) opt.command_line.push_back(' ');
        opt.command_line += argv[i];
    }
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--") {
            for (++i; i < argc; ++i) pos.push_back(argv[i]);
            break;
        }
        if (a.size() < 2 || a[0] != '-') { pos.push_back(a); continue; }
        std::string name = a, val;
        bool has_val = false;
        if (a.compare(0, 2, "--") == 0) {
            size_t eq = a.find('=');
            if (eq != std::string::npos) { name = a.substr(0, eq); val = a.substr(eq + 1); has_val = true; }
        } else if (a.size() > 2) {
            name = a.substr(0, 2);
            val = a.substr(2);
            has_val = true;
        }
        auto old = old_options().find(name);
        if (old != old_options().end()) usage_error(strf("option %s: %s", name.c_str(), old->second.c_str()));
        const OptDef* def = nullptr;
        for (const OptDef& d : OPTS)
            if (name == d.name) { def = &d; break; }
        if (!def) usage_error("unknown option " + name);
        if (def->arg) {
            if (!has_val) {
                if (i + 1 >= argc) usage_error(name + " needs a value");
                val = argv[++i];
            }
        } else if (has_val) {
            usage_error(name + " takes no value");
        }
        switch (def->id) {
            case O_M: opt.minimap2 = val; break;
            case O_O: opt.out_path = val; break;
            case O_R: opt.report_path = val; break;
            case O_T: opt.threads = (int)need_int(name, val, 1, 4096); break;
            case O_J: opt.jobs = (int)need_int(name, val, 1, 4096); break;
            case O_B: opt.b = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_I: opt.ident = need_double(name, val, 0.0, 1.0); if (opt.ident <= 0) usage_error("-i must be > 0"); break;
            case O_G: opt.G = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_MINPIECE: opt.min_piece = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_ISLAND: opt.island = need_int(name, val, 0, INT64_MAX / 4); break;
            case O_DELTA: opt.delta = need_double(name, val, 0.0, 0.999999); break;
            case O_TIES:
                if (val == "resolve") opt.ties_refuse = false;
                else if (val == "refuse") opt.ties_refuse = true;
                else usage_error("--ties: resolve or refuse");
                break;
            case O_FRAG: opt.frag = need_int(name, val, 0, INT64_MAX / 4); break;
            case O_PREFILTER: {
                size_t c = val.find(',');
                if (c == std::string::npos) usage_error("--prefilter: expected NMATCH,COVERAGE");
                opt.prefilter_nmatch = need_int(name, val.substr(0, c), 0, INT64_MAX / 4);
                opt.prefilter_cov = need_double(name, val.substr(c + 1), 0.0, 1.0);
                break;
            }
            case O_SCREEN: opt.screen_sample = (int)need_int(name, val, 1, 1000000); break;
            case O_MAXPAIR: opt.max_pair = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_MAXSITEQ: opt.max_site_query = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_MAXSITEN: opt.max_site_nodes = need_int(name, val, 1, INT64_MAX / 4); break;
            case O_X: opt.preset = val; break;
            case O_XX: opt.mm_extra = val; break;
            case O_NOALT: opt.alt = false; break;
            case O_ALTROUNDS: opt.alt_rounds = (int)need_int(name, val, 0, 1000); break;
            case O_WALKS:
                if (val == "creator") opt.walks = WalksKind::CREATOR;
                else if (val == "witness") opt.walks = WalksKind::WITNESS;
                else if (val.compare(0, 4, "gaf:") == 0 && val.size() > 4) { opt.walks = WalksKind::GAF; opt.gaf_path = val.substr(4); }
                else usage_error("--walks: creator, witness or gaf:FILE");
                break;
            case O_MAXREP: opt.max_repeat_frac = need_double(name, val, 0.0, 1.0); break;
            case O_DETECTONLY: opt.detect_only = true; break;
            case O_CHECK: opt.check = true; break;
            case O_STRICT: opt.strict = true; break;
            case O_MEM:
                if (!parse_mem(val, opt.mem)) usage_error("--mem: e.g. 16G, 1500M or all");
                break;
            case O_TMPDIR: opt.tmpdir = val; break;
            case O_IDBASE: opt.id_base = need_int(name, val, 0, INT64_MAX / 4); break;
            case O_AUDIT: opt.audit_gaf = val; break;
            case O_REGION: opt.regions.push_back(parse_region(val)); break;
            case O_STAGE:
                if (val == "detect") opt.stage = Stage::DETECT;
                else if (val == "align") opt.stage = Stage::ALIGN;
                else if (val == "plan") opt.stage = Stage::PLAN;
                else if (val == "full") opt.stage = Stage::FULL;
                else usage_error("--stage: detect, align, plan or full");
                break;
            case O_INJECT: opt.inject_chains = val; break;
            case O_DUMP: opt.dump_dir = val; break;
            case O_VERBOSE: ++log_level(); break;
            case O_HELP: return false;
        }
    }
    if (pos.size() != 2) usage_error(strf("expected <in.gfa> <snarls.json>, got %zu positional argument(s)", pos.size()));
    opt.gfa_path = pos[0];
    opt.snarls_path = pos[1];
    if (opt.minimap2.empty()) usage_error("-m <minimap2> is required (there is no PATH default)");
    if (opt.stage == Stage::FULL) {
        if (opt.out_path.empty()) usage_error("-o <out.gfa> is required");
        if (opt.report_path.empty()) usage_error("-r <report.tsv> is required");
    }
    check_output_paths(opt);
    // aligner processes: -j defaults to -t and is lowered so j x 1.5 GB fits --mem.  1.5 GB is
    // 1.5e9 bytes, as the spec and cactus (ZIP_ALIGNER_MEMORY) count it: with 1.5 GiB, cactus's
    // --mem of 1.5e9 x cores gave every job one aligner process too few.
    {
        int64_t mem = opt.mem > 0 ? opt.mem : physical_memory_bytes();
        const int64_t per = 1500000000LL;
        int j = opt.jobs > 0 ? opt.jobs : opt.threads;
        if (mem > 0) {
            int64_t jmax = std::max<int64_t>(1, mem / per);
            if (j > jmax) {
                ZLOG("lowering -j from %d to %lld so that j x 1.5 GB (1.5e9 bytes) fits --mem", j, (long long)jmax);
                j = (int)jmax;
            }
        }
        opt.jobs = j;
    }
    if (opt.tmpdir.empty()) {
        const char* t = getenv("TMPDIR");
        opt.tmpdir = (t && *t) ? t : "/tmp";
    }
    {
        struct stat st;
        if (stat(opt.tmpdir.c_str(), &st) != 0 || !S_ISDIR(st.st_mode) || access(opt.tmpdir.c_str(), W_OK) != 0)
            usage_error("--tmpdir " + opt.tmpdir + " is not a writable directory");
    }
    return true;
}

// ================================================================ per-site processing

bool in_regions(const Graph& g, const Site& s, const Options& opt) {
    if (opt.regions.empty()) return true;
    const std::string& sn = g.sn_names[s.sn];
    for (const Region& r : opt.regions) {
        bool match = sn == r.sn || (sn.size() > r.sn.size() && sn.compare(sn.size() - r.sn.size(), r.sn.size(), r.sn) == 0 &&
                                    sn[sn.size() - r.sn.size() - 1] == '#');
        if (match && s.hi >= r.lo && s.lo <= r.hi) return true;
    }
    return false;
}

struct SiteCounters {
    uint64_t units = 0, kinds[4] = {0, 0, 0, 0}, ref_units = 0, to_align = 0, infeasible = 0, pair_too_big = 0, capped = 0,
             no_window = 0, small = 0;
    int64_t ref_query = 0, align_query = 0, align_feasible = 0;
    uint64_t allowed[4] = {0, 0, 0, 0};
    void add(const SiteCounters& o) {
        units += o.units;
        for (int k = 0; k < 4; ++k) { kinds[k] += o.kinds[k]; allowed[k] += o.allowed[k]; }
        ref_units += o.ref_units; to_align += o.to_align; infeasible += o.infeasible; pair_too_big += o.pair_too_big;
        capped += o.capped; no_window += o.no_window; small += o.small;
        ref_query += o.ref_query; align_query += o.align_query; align_feasible += o.align_feasible;
    }
};

struct SiteResult {
    bool processed = false;
    std::vector<ReportRow> rows;        // final report rows of the site, in order
    SitePlan plan;
    SiteCounters cnt;
    WalkStats ws;
    std::string detect_text, dump_units, dump_allowed;
    GafAuditStats audit;
    std::string audit_detail;
    AltStats alt;
    uint64_t split_inv = 0;
    double seconds = 0;                 // wall time of the site
    AlignTiming timing;                 // its aligner work: minimap2 running, waiting for a -j slot, minimap2's own times
};

// ================================================================ diagnostics

// spec "Diagnostics only": split-inversion-candidate.  Two consecutive F excursions of one contig
// (adjacent in the contig's excursions of the site, ordered along the reference; a contig with a J
// excursion in the site is skipped, as it does not walk the reference forward between them) that
// both lack a chain in their own windows (no-chain, below-b or prefiltered), are each >= b, and are
// separated by <= 10 kb of reference.  Each half of a split inversion aligns to the other half's
// window, so both are reference units (windows >= ceil(b*i)); a long insertion with a small window
// never qualifies, which keeps the diagnostic cheap.  The merged query -- the first alt run, the
// reference the contig walks between (the first unit's arrival node through the second's departure
// node), the second alt run -- is aligned to the merged window [wlo1, whi2) with no feasibility
// cut; a confident chain with one '-' segment (consecutive '-' parts, collinear on '-') covering
// >= b of each alt run is reported.  Never zipped.
const int64_t SPLIT_INV_MAX_GAP = 10000;

bool lacks_chain(const std::string& outcome) { return outcome == "no-chain" || outcome == "below-b" || outcome == "prefiltered"; }

int64_t overlap_bp(int64_t a, int64_t b, int64_t c, int64_t d) { return std::max<int64_t>(0, std::min(b, d) - std::max(a, c)); }

void split_inversion_rows(const Graph& g, const SiteData& sd, const std::vector<AlignResult>& res, Aligner& aligner, const Options& opt,
                          std::vector<ReportRow>& out) {
    const Site& s = *sd.site;
    const std::string label = s.label(g);
    // each contig's excursions of the site, in reference order
    std::map<std::string, std::vector<uint32_t>> by_contig;
    std::set<std::string> walks_j;
    for (const Unit& u : sd.units) {
        std::set<std::string> cs;
        for (const Owner& o : u.exc.owners) cs.insert(o.contig);
        for (const std::string& c : cs) {
            if (u.kind == Kind::J) walks_j.insert(c);
            else by_contig[c].push_back(u.id);
        }
    }
    std::map<std::pair<uint32_t, uint32_t>, std::vector<std::string>> pairs;   // (u1, u2) -> contigs
    for (auto& kv : by_contig) {
        if (walks_j.count(kv.first)) continue;
        std::vector<uint32_t>& v = kv.second;
        std::sort(v.begin(), v.end(), [&](uint32_t x, uint32_t y) {
            const Unit& a = sd.units[x];
            const Unit& b = sd.units[y];
            int64_t a0 = std::min(a.wlo, a.whi), b0 = std::min(b.wlo, b.whi);
            if (a0 != b0) return a0 < b0;
            if (a.whi != b.whi) return a.whi < b.whi;
            return x < y;
        });
        for (size_t k = 0; k + 1 < v.size(); ++k) {
            const Unit& u1 = sd.units[v[k]];
            const Unit& u2 = sd.units[v[k + 1]];
            if (!u1.ref_unit || !u2.ref_unit) continue;      // F, query >= b, window >= ceil(b*i)
            if (!lacks_chain(res[u1.id].outcome) || !lacks_chain(res[u2.id].outcome)) continue;
            NodeId r0 = handle_node(u1.exc.arr), r1 = handle_node(u2.exc.dep);
            int64_t gap = u2.wlo - u1.whi;
            if (gap < 0 || gap > SPLIT_INV_MAX_GAP || g.start(r0) > g.start(r1)) continue;
            pairs[std::make_pair(u1.id, u2.id)].push_back(kv.first);
        }
    }
    if (pairs.empty()) return;
    std::vector<AlignJob> jobs;
    std::vector<std::pair<uint32_t, uint32_t>> which;
    std::vector<int64_t> mid_bp;
    const std::vector<NodeId>& ref = g.ref_nodes[(size_t)s.sn];
    for (const auto& kv : pairs) {
        const Unit& u1 = sd.units[kv.first.first];
        const Unit& u2 = sd.units[kv.first.second];
        AlignJob j;
        j.key = strf("%s:split:%u+%u", label.c_str(), u1.id, u2.id);
        j.query = u1.exc.alts;
        NodeId r0 = handle_node(u1.exc.arr), r1 = handle_node(u2.exc.dep);
        auto pos = [&](NodeId n) {
            return (size_t)(std::lower_bound(ref.begin(), ref.end(), n, [&](NodeId a, NodeId b) { return g.start(a) < g.start(b); }) - ref.begin());
        };
        int64_t mid = 0;
        for (size_t k = pos(r0); k < ref.size() && k <= pos(r1); ++k) {
            j.query.push_back(make_handle(ref[k], false));
            mid += g.len(ref[k]);
        }
        j.query.insert(j.query.end(), u2.exc.alts.begin(), u2.exc.alts.end());
        j.target.ref = true;
        j.target.sn = s.sn;
        j.target.lo = u1.wlo;
        j.target.hi = u2.whi;
        if (u1.query_bp + mid + u2.query_bp > opt.max_pair || u2.whi - u1.wlo > opt.max_pair) continue;
        jobs.push_back(std::move(j));
        which.push_back(kv.first);
        mid_bp.push_back(mid);
    }
    std::vector<AlignResult> jr;
    aligner.align_jobs(g, jobs, jr);
    for (size_t k = 0; k < jobs.size(); ++k) {
        const AlignResult& r = jr[k];
        if (!r.outcome.empty() || r.chain.parts.empty()) continue;
        const Unit& u1 = sd.units[which[k].first];
        const Unit& u2 = sd.units[which[k].second];
        const int64_t a1 = u1.query_bp, b2 = u1.query_bp + mid_bp[k], e2 = b2 + u2.query_bp;
        const std::vector<ChainPart>& parts = r.chain.parts;
        bool hit = false;
        for (size_t i = 0; i < parts.size() && !hit;) {
            if (parts[i].strand != '-') { ++i; continue; }
            int64_t c1 = 0, c2 = 0;
            size_t j = i;
            for (; j < parts.size() && parts[j].strand == '-' && (j == i || parts[j].te <= parts[j - 1].ts + opt.G); ++j) {
                c1 += overlap_bp(parts[j].qs, parts[j].qe, 0, a1);
                c2 += overlap_bp(parts[j].qs, parts[j].qe, b2, e2);
            }
            if (c1 >= opt.b && c2 >= opt.b) hit = true;
            i = j;
        }
        if (!hit) continue;
        ReportRow row;
        row.site = label;
        row.sn = g.sn_names[(size_t)s.sn];
        row.lo = s.lo;
        row.hi = s.hi;
        row.pass = "diag";
        row.round = 0;
        row.source = src_name(std::max(u1.exc.src, u2.exc.src));
        for (const std::string& c : pairs[which[k]]) row.owners += (row.owners.empty() ? "" : ",") + c;
        row.kind = "F";
        row.anchors = g.handle_str(u1.exc.dep) + ">" + g.handle_str(u2.exc.arr);
        row.window = strf("%lld-%lld", (long long)u1.wlo, (long long)u2.whi);
        row.query_bp = e2;
        fill_align_row(row, r, opt.G);
        row.outcome = "split-inversion-candidate";
        out.push_back(std::move(row));
    }
}

void process_site(const Graph& g, const Site& s, const WalkSource& src, Aligner* aligner, const GafAudit* audit, const Options& opt,
                  SiteResult& r) {
    const double t0 = steady_seconds();
    AlignTimingScope timing;            // every aligner call of this site (reference, diagnostics, alt-vs-alt)
    r.processed = true;
    SiteData sd;
    prepare_site(g, s, src, opt, sd);
    r.ws = sd.wstats;
    std::vector<ReportRow> unit_rows(sd.units.size()), extra_rows;
    for (const Unit& u : sd.units) {
        unit_rows[u.id] = unit_row(g, sd, u);
        SiteCounters& c = r.cnt;
        ++c.units;
        c.kinds[(int)u.kind]++;
        if (u.ref_unit) { ++c.ref_units; c.ref_query += u.query_bp; }
        if (u.ref_unit && u.outcome.empty()) { ++c.to_align; c.align_query += u.query_bp; c.align_feasible += u.feasible_bp; }
        if (u.outcome == "infeasible") ++c.infeasible;
        else if (u.outcome == "pair-too-big") ++c.pair_too_big;
        else if (u.outcome == "site:capped") ++c.capped;
        else if (u.outcome == "no-window") ++c.no_window;
        else if (u.outcome == "small") ++c.small;
    }
    for (const AllowedIv& a : sd.allowed) r.cnt.allowed[(int)a.status]++;
    if (opt.stage == Stage::DETECT) r.detect_text = detect_site_lines(g, s, &sd);
    if (!opt.dump_dir.empty()) {
        r.dump_units = dump_units(g, sd);
        r.dump_allowed = dump_allowed(g, sd);
    }
    std::vector<AlignResult> res;
    std::vector<ReportRow> diag_rows;
    std::deque<ReportRow> alt_rows;
    if (opt.stage >= Stage::ALIGN && aligner) {
        res.assign(sd.units.size(), AlignResult());
        aligner->align_reference(g, sd, res, unit_rows);
        split_inversion_rows(g, sd, res, *aligner, opt, diag_rows);
        r.split_inv = diag_rows.size();
    }
    if (opt.stage >= Stage::PLAN) {
        // one plan per site: the reference pass, then alt-vs-alt on the same planner, then the
        // transaction (V1-V7)
        SitePlanner planner(g, sd, opt);
        planner.accept_reference(res, unit_rows);
        // v2 (spec step 5, zip_alt): alt-vs-alt shares the planner's piece map, acceptance rules
        // and transaction, and the aligner; its rows sit in a deque because the planner keeps
        // pointers to them, and they are written after the site's unit rows.  --no-alt (v1)
        // skips it.
        if (opt.alt) alt_pass(g, sd, aligner, planner, opt, alt_rows, r.alt);
        r.plan = planner.finish(unit_rows, extra_rows);
        if (opt.alt) alt_tally(alt_rows, r.alt);
        if (audit) r.audit = audit->audit_plan(s, r.plan.zipped, &r.audit_detail);
        // --detect-only emits the input graph unchanged: what would have been zipped says so
        if (opt.detect_only) {
            for (ReportRow& row : unit_rows)
                if (row.outcome == "zipped") row.outcome = "would-zip";
            for (ReportRow& row : alt_rows)
                if (row.outcome == "zipped") row.outcome = "would-zip";
        }
    }
    for (ReportRow& row : unit_rows)
        if (row.outcome != "small") r.rows.push_back(std::move(row));
    for (ReportRow& row : alt_rows) r.rows.push_back(std::move(row));
    for (ReportRow& row : diag_rows) r.rows.push_back(std::move(row));
    for (ReportRow& row : extra_rows) r.rows.push_back(std::move(row));
    r.seconds = steady_seconds() - t0;
    r.timing = timing.timing();
}

// A site's wall time, split: minimap2 running (some process of the site running), waiting for an
// aligner slot (-j; none of its processes running), and rgfa-zip's own work; then minimap2's own
// real and CPU time, which is what the site's alignments cost whatever the contention.  At level 1
// only for sites that took more than 30 s outside the waits; with -v for every site that ran
// minimap2 or took more than 1 s.
void log_site_time(const Graph& g, const Site& s, const SiteResult& r) {
    const double run = r.timing.running(), wait = r.timing.waiting(), busy = r.seconds - wait;
    if (!(busy > 30 || (log_level() > 1 && (busy > 1 || r.timing.processes > 0)))) return;
    log_msg(busy > 30 ? 1 : 2,
            "site %s: %.1f s: minimap2 running %.1f s, waiting for an aligner slot %.1f s, rgfa-zip alone %.1f s; %llu minimap2 "
            "process(es): real %.1f s, CPU %.1f s (minimap2's own figures)",
            s.label(g).c_str(), r.seconds, run, wait, std::max(0.0, r.seconds - run - wait), (unsigned long long)r.timing.processes,
            r.timing.mm2_real, r.timing.mm2_cpu);
}

void write_stdout(const std::string& s) {
    size_t off = 0;
    while (off < s.size()) {
        size_t n = fwrite(s.data() + off, 1, s.size() - off, stdout);
        if (n == 0) fail(EXIT_IO, strf("write to stdout failed: %s", strerror(errno)));
        off += n;
    }
}

// ================================================================ run

// The snarls must cover the graph: every alt node lies in the interior of some top-level site
// (processed or skipped).  A truncated or partial snarls file leaves sites out without any other
// symptom, so more than SNARL_UNCOVERED_MAX of the alt bp outside every site is invalid input.  A
// valid decomposition leaves out only alt nodes that hang off the ends of the reference (tips at
// a chromosome end are in no snarl).
const double SNARL_UNCOVERED_MAX = 0.01;

void check_snarl_coverage(const Options& opt, const SiteStats& st) {
    double frac = st.alt_bp > 0 ? (double)st.uncovered_bp / (double)st.alt_bp : 0.0;
    ZLOG("snarls cover %llu of %llu alt nodes (%lld of %lld bp); %llu alt node(s), %lld bp (%.3f%%), lie in no site",
         (unsigned long long)(st.alt_nodes - st.uncovered_nodes), (unsigned long long)st.alt_nodes, (long long)(st.alt_bp - st.uncovered_bp),
         (long long)st.alt_bp, (unsigned long long)st.uncovered_nodes, (long long)st.uncovered_bp, 100.0 * frac);
    if (st.alt_nodes > 0 && st.snarls == 0)
        fail(EXIT_INPUT, strf("%s has no snarls, but the graph has %llu alt nodes: the snarls file is empty or truncated", opt.snarls_path.c_str(),
                              (unsigned long long)st.alt_nodes));
    if (frac > SNARL_UNCOVERED_MAX)
        fail(EXIT_INPUT, strf("%.2f%% of the graph's alt bp (%llu nodes, %lld bp) lie in no snarl of %s (at most %.0f%% may): the snarls file "
                              "is truncated, partial or of another graph",
                              100.0 * frac, (unsigned long long)st.uncovered_nodes, (long long)st.uncovered_bp, opt.snarls_path.c_str(),
                              100.0 * SNARL_UNCOVERED_MAX));
}

int run(Options& opt) {
    auto t0 = std::chrono::steady_clock::now();
    auto secs = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    ZLOG("rgfa-zip %s: %s", ZIP_VERSION, opt.command_line.c_str());
    std::string mm2_version = minimap2_version(opt.minimap2);
    ZLOG("minimap2 %s (%s); -t %d, -j %d, tmpdir %s", mm2_version.c_str(), opt.minimap2.c_str(), opt.threads, opt.jobs, opt.tmpdir.c_str());

    Graph g;
    read_rgfa(opt.gfa_path, g);
    std::vector<Snarl> snarls = read_snarls(opt.snarls_path, g);
    Runs runs = build_runs(g);
    SiteStats sst;
    std::vector<Site> sites = build_sites(g, snarls, opt, sst);
    snarls.clear();
    snarls.shrink_to_fit();
    ZLOG("%llu snarls: %llu nested, %llu with a non-rank-0 boundary, %llu across SNs, %llu with one boundary node, %llu duplicate sites",
         (unsigned long long)sst.snarls, (unsigned long long)sst.nested, (unsigned long long)sst.nonref_boundary,
         (unsigned long long)sst.cross_sn, (unsigned long long)sst.same_node, (unsigned long long)sst.duplicate);
    ZLOG("%llu sites: %llu with alt nodes (%llu pure insertion), %llu without, %llu site:leak, %llu site:ref-mismatch, %llu site:too-big"
         "%s; %zu runs; %.1f s",
         (unsigned long long)sst.sites, (unsigned long long)sst.ok, (unsigned long long)sst.pure_insertion,
         (unsigned long long)sst.no_alt, (unsigned long long)sst.leak, (unsigned long long)sst.ref_mismatch,
         (unsigned long long)sst.too_big,
         sst.orientation_disagrees ? strf(" (%llu with boundary flags against the SO order)", (unsigned long long)sst.orientation_disagrees).c_str() : "",
         runs.runs.size(), secs());
    check_snarl_coverage(opt, sst);

    // the walk source
    std::unique_ptr<WalkSource> src;
    if (opt.walks == WalksKind::CREATOR) src.reset(new CreatorWalks(g, runs));
    else if (opt.walks == WalksKind::WITNESS) src.reset(new WitnessWalks(g, runs));
    else src.reset(new GafWalks(g, runs, opt.gaf_path, opt));

    std::unique_ptr<Aligner> aligner;
    if (opt.stage >= Stage::ALIGN) aligner.reset(new Aligner(opt));
    if (!opt.dump_dir.empty() && !make_dirs(opt.dump_dir)) fail(EXIT_IO, "cannot create --dump directory " + opt.dump_dir);
    // the GAF audit (release gate): shares the walk source's index when it reads the same GAF
    std::unique_ptr<GafAudit> audit;
    if (!opt.audit_gaf.empty()) {
        if (opt.stage < Stage::PLAN) {
            ZLOG("note: --audit-gaf needs --stage plan or full; no audit");
        } else if (opt.walks == WalksKind::GAF && opt.audit_gaf == opt.gaf_path) {
            audit.reset(new GafAudit(g, static_cast<const GafWalks&>(*src).index()));
        } else {
            Options aopt = opt;
            aopt.dump_dir.clear();      // the walk source's GAF dumps are its own
            audit.reset(new GafAudit(g, opt.audit_gaf, aopt));
        }
    }

    // sites to process, and skip rows
    std::vector<uint32_t> proc;
    std::vector<uint8_t> selected(sites.size(), 0);
    for (const Site& s : sites) {
        if (!in_regions(g, s, opt)) continue;
        selected[s.index] = 1;
        if (s.status == SiteStatus::OK) proc.push_back(s.index);
    }
    std::vector<SiteResult> results(sites.size());
    std::atomic<size_t> done(0);
    const size_t report_every = std::max<size_t>(1, proc.size() / 10);
    parallel_for(proc.size(), opt.threads, [&](size_t i, int) {
        const Site& s = sites[proc[i]];
        process_site(g, s, *src, aligner.get(), audit.get(), opt, results[s.index]);
        log_site_time(g, s, results[s.index]);
        size_t d = ++done;
        if (d % report_every == 0 && proc.size() >= 20) ZLOG("%zu / %zu sites, %.1f s", d, proc.size(), secs());
    });
    if (aligner) aligner->check_failure_rate();

    // totals
    SiteCounters tot;
    WalkStats ws;
    GafAuditStats audit_tot;
    AltStats alt_tot;
    uint64_t accepted = 0, accepted_new = 0, reverted = 0, split_inv = 0, g_chains = 0, g_pieces = 0;
    int64_t zipped = 0, g_bp = 0;
    for (const SiteResult& r : results) {
        if (!r.processed) continue;
        tot.add(r.cnt);
        ws.add(r.ws);
        accepted += r.plan.accepted;
        accepted_new += r.plan.accepted_new;
        reverted += r.plan.reverted;
        zipped += r.plan.zipped_bp;
        g_chains += r.plan.g_chains;
        g_pieces += r.plan.g_pieces;
        g_bp += r.plan.g_bp;
        split_inv += r.split_inv;
        audit_tot.add(r.audit);
        alt_tot.add(r.alt);
    }
    uint64_t n_sel = 0;
    for (uint8_t x : selected) n_sel += x;
    ZLOG("processed %zu of %llu selected sites in %.1f s", proc.size(), (unsigned long long)n_sel, secs());
    ZLOG("runs %llu: creator %llu, witness fallback %llu (creator links %llu, own-link ambiguous %llu, creator edge %llu, cycle %llu, other %llu) "
         "-> %llu by DAG, %llu by bubble, %llu unanchored",
         (unsigned long long)ws.runs, (unsigned long long)ws.creator_ok, (unsigned long long)ws.fallback,
         (unsigned long long)ws.fail_creator_links, (unsigned long long)ws.fail_own_ambiguous, (unsigned long long)ws.fail_creator_edge,
         (unsigned long long)ws.fail_cycle, (unsigned long long)ws.fail_other, (unsigned long long)ws.witness_dag,
         (unsigned long long)ws.witness_bubble, (unsigned long long)ws.unanchored);
    ZLOG("excursions %llu distinct (%llu merged): F %llu, I %llu, BK %llu, J %llu",
         (unsigned long long)tot.units, (unsigned long long)ws.merged, (unsigned long long)tot.kinds[0],
         (unsigned long long)tot.kinds[1], (unsigned long long)tot.kinds[2], (unsigned long long)tot.kinds[3]);
    ZLOG("reference units %llu (query %lld bp): %llu to align (query %lld bp, feasible %lld bp), %llu infeasible, %llu pair-too-big, "
         "%llu site:capped; %llu no-window",
         (unsigned long long)tot.ref_units, (long long)tot.ref_query, (unsigned long long)tot.to_align, (long long)tot.align_query, (long long)tot.align_feasible,
         (unsigned long long)tot.infeasible, (unsigned long long)tot.pair_too_big, (unsigned long long)tot.capped,
         (unsigned long long)tot.no_window);
    ZLOG("Allowed over alt nodes: ok %llu, empty %llu, blocked %llu, none %llu", (unsigned long long)tot.allowed[(int)AllowedStatus::OK],
         (unsigned long long)tot.allowed[(int)AllowedStatus::EMPTY], (unsigned long long)tot.allowed[(int)AllowedStatus::BLOCKED],
         (unsigned long long)tot.allowed[(int)AllowedStatus::NONE]);
    if (aligner)
        ZLOG("aligner: %llu windows, %llu failed, %llu paf-inconsistent", (unsigned long long)aligner->windows(),
             (unsigned long long)aligner->windows_failed(), (unsigned long long)aligner->paf_inconsistent());
    if (aligner) ZLOG("diagnostics: %llu split-inversion-candidate(s)", (unsigned long long)split_inv);
    if (opt.alt && opt.stage >= Stage::PLAN)
        ZLOG("alt-vs-alt: %llu residual excursion(s), %llu pair(s), %llu group(s), %llu candidate(s): %llu shares-nodes, %llu not-private, "
             "%llu in-series, %llu pair-too-big, %llu site:capped; %llu aligned, %llu confident; %llu zipped (%lld bp), %llu refused, "
             "%llu reverted; %llu near-parallel",
             (unsigned long long)alt_tot.residual, (unsigned long long)alt_tot.pairs, (unsigned long long)alt_tot.groups,
             (unsigned long long)alt_tot.candidates, (unsigned long long)alt_tot.shares_nodes, (unsigned long long)alt_tot.not_private,
             (unsigned long long)alt_tot.in_series, (unsigned long long)alt_tot.pair_too_big, (unsigned long long)alt_tot.capped,
             (unsigned long long)alt_tot.aligned, (unsigned long long)alt_tot.confident, (unsigned long long)alt_tot.zipped,
             (long long)alt_tot.zipped_bp, (unsigned long long)alt_tot.refused, (unsigned long long)alt_tot.reverted,
             (unsigned long long)alt_tot.near_parallel);
    if (opt.walks == WalksKind::GAF && opt.stage >= Stage::PLAN)
        ZLOG("rule (G) (a GAF record through the piece's node also reads its target): %llu piece(s), %lld bp, refused in %llu chain(s)",
             (unsigned long long)g_pieces, (long long)g_bp, (unsigned long long)g_chains);
    // the bp are the sum of the report's kept_bp: a node interval that several compatible chains
    // zip is counted once, by the chain that zipped it first
    const std::string own = strf("%llu of the chain(s) zip bp of their own, %llu only agree with an earlier chain's zip",
                                 (unsigned long long)accepted_new, (unsigned long long)(accepted - accepted_new));
    if (opt.stage >= Stage::PLAN && opt.detect_only)
        ZLOG("would zip %llu chain(s), %lld bp; %llu reverted (--detect-only: the graph is written unchanged); %s", (unsigned long long)accepted,
             (long long)zipped, (unsigned long long)reverted, own.c_str());
    else if (opt.stage >= Stage::PLAN)
        ZLOG("accepted %llu chain(s), %lld bp zipped; %llu reverted; %s", (unsigned long long)accepted, (long long)zipped,
             (unsigned long long)reverted, own.c_str());
    if (audit) ZLOG("%s", audit_tot.summary().c_str());
    if (opt.strict && reverted > 0) fail(EXIT_INVARIANT, strf("--strict: the transaction dropped %llu chain(s)", (unsigned long long)reverted));
    if (opt.strict && audit && !audit_tot.clean()) fail(EXIT_INVARIANT, "--strict: the GAF audit failed");

    // the edit
    bool write_graph = !opt.out_path.empty() && (opt.stage == Stage::FULL || opt.detect_only);
    if (opt.stage == Stage::FULL && !opt.detect_only) {
        std::vector<SitePlan> plans;
        for (const SiteResult& r : results)
            if (r.processed) plans.push_back(r.plan);
        apply_plans(g, plans, opt);
    }
    if (write_graph) check_placement(g);
    if (!opt.out_path.empty() && !write_graph) ZLOG("note: --stage stops before the edit; no graph written to %s", opt.out_path.c_str());

    // report rows in site order
    std::vector<std::vector<ReportRow>> skip_rows(sites.size());
    std::vector<const std::vector<ReportRow>*> rows;
    for (const Site& s : sites) {
        if (!selected[s.index]) continue;
        if (s.status == SiteStatus::LEAK || s.status == SiteStatus::REF_MISMATCH || s.status == SiteStatus::TOO_BIG) {
            skip_rows[s.index].push_back(site_skip_row(g, s));
            rows.push_back(&skip_rows[s.index]);
        } else if (results[s.index].processed) {
            rows.push_back(&results[s.index].rows);
        }
    }

    // Outputs, all or nothing: every file is written, flushed and fsynced as <path>.tmp; the --stage
    // detect printout goes to stdout; only then is anything renamed into place, and if one rename
    // fails the outputs already renamed are removed again.
    std::unique_ptr<AtomicWriter> gw, rw;
    if (write_graph) {
        gw.reset(new AtomicWriter(opt.out_path));
        write_rgfa(*gw, g);
        gw->finish();
    }
    if (!opt.report_path.empty()) {
        rw.reset(new AtomicWriter(opt.report_path));
        write_report(*rw, report_header(opt, mm2_version), rows);
        rw->finish();
    }
    std::vector<std::unique_ptr<AtomicWriter>> dumps;
    if (!opt.dump_dir.empty()) {
        dumps.emplace_back(new AtomicWriter(opt.dump_dir + "/units.tsv"));
        dumps.back()->write(dump_units_header());
        for (const SiteResult& r : results) dumps.back()->write(r.dump_units);
        dumps.emplace_back(new AtomicWriter(opt.dump_dir + "/allowed.tsv"));
        dumps.back()->write(dump_allowed_header());
        for (const SiteResult& r : results) dumps.back()->write(r.dump_allowed);
        if (audit) {
            dumps.emplace_back(new AtomicWriter(opt.dump_dir + "/gafaudit.tsv"));
            dumps.back()->write("#" + audit_tot.summary() + "\n");
            for (const SiteResult& r : results) dumps.back()->write(r.audit_detail);
        }
        for (auto& d : dumps) d->finish();
    }
    if (opt.stage == Stage::DETECT) {
        std::string out = detect_header();
        write_stdout(out);
        for (const Site& s : sites) {
            if (!selected[s.index]) continue;
            if (results[s.index].processed) write_stdout(results[s.index].detect_text);
            else if (s.status != SiteStatus::NO_ALT) write_stdout(detect_site_lines(g, s, nullptr));
        }
        if (fflush(stdout) != 0) fail(EXIT_IO, strf("write to stdout failed: %s", strerror(errno)));
    }
    std::vector<AtomicWriter*> commit;
    if (gw) commit.push_back(gw.get());
    if (rw) commit.push_back(rw.get());
    for (auto& d : dumps) commit.push_back(d.get());
    if (opt.walks == WalksKind::GAF)
        for (AtomicWriter* w : gaf_dump_writers(static_cast<const GafWalks&>(*src))) commit.push_back(w);
    check_abort();
    try {
        for (AtomicWriter* w : commit) w->commit();
    } catch (...) {
        for (AtomicWriter* w : commit) w->rollback();
        throw;
    }
    g_committed = true;
    ZLOG("done in %.1f s, peak RSS %.2f GB", secs(), self_peak_rss_kb() / 1048576.0);
    return EXIT_OK;
}

// ================================================================ stop signals

// SIGTERM, SIGINT and SIGHUP are blocked in every thread (main blocks them before any thread
// starts; children get an empty mask) and taken by one thread with sigwait.  It records the signal
// and requests an abort: workers stop, running children get SIGTERM, and main unwinds as from an
// error -- RAII removes the temporary directory and the .tmp outputs -- and then re-raises the
// signal, so the exit status says what happened.  If main has not finished 10 s later (or a second
// stop signal arrives), the signal thread removes the registered temporary paths itself and
// re-raises.
[[noreturn]] void reraise(int sig) {
    signal(sig, SIG_DFL);
    sigset_t s;
    sigemptyset(&s);
    sigaddset(&s, sig);
    pthread_sigmask(SIG_UNBLOCK, &s, nullptr);
    raise(sig);
    _exit(128 + sig);
}

void start_signal_thread() {
    sigset_t set;
    sigemptyset(&set);
    for (int s : {SIGTERM, SIGINT, SIGHUP}) sigaddset(&set, s);
    pthread_sigmask(SIG_BLOCK, &set, nullptr);
    try {
        std::thread([set]() {
            int sig = 0;
            while (sigwait(&set, &sig) != 0) {}
            Runtime::get().set_signal(sig);
            Runtime::get().request_abort();
            struct timespec t = {0, 100000000};   // 100 ms
            for (int i = 0; i < 100; ++i)
                if (sigtimedwait(&set, nullptr, &t) > 0) break;   // a second stop signal: do not wait
            Runtime::get().remove_paths();
            reraise(sig);
        }).detach();
    } catch (const std::exception& e) {
        // no signal thread: stop signals keep their default action (nothing is cleaned up)
        pthread_sigmask(SIG_UNBLOCK, &set, nullptr);
        log_msg(1, "warning: could not start the signal thread (%s)", e.what());
    }
}

} // namespace

int main(int argc, char** argv) {
    // SIGPIPE: write errors are reported, not fatal.  SIGCHLD: an ignored disposition inherited from
    // the launcher makes the kernel reap children itself, and their exit status would be lost.
    signal(SIGPIPE, SIG_IGN);
    signal(SIGCHLD, SIG_DFL);
    start_signal_thread();
    Options opt;
    int code = EXIT_OK;
    std::string msg;
    try {
        if (!parse_args(argc, argv, opt)) {
            fputs(USAGE, stdout);
            return EXIT_OK;
        }
        code = run(opt);
    } catch (const ZipError& e) {
        msg = strf("[rgfa-zip] error: %s\n", e.what());
        code = e.code;
    } catch (const Aborted&) {
        msg = "[rgfa-zip] error: aborted\n";
        code = EXIT_INVARIANT;
    } catch (const std::bad_alloc&) {
        msg = "[rgfa-zip] error: out of memory\n";
        code = EXIT_INVARIANT;
    } catch (const std::exception& e) {
        msg = strf("[rgfa-zip] internal error: %s\n", e.what());
        code = EXIT_INVARIANT;
    }
    // a stop signal ends the run whatever happened while it unwound (errors caused by the stop,
    // e.g. a decompressor killed mid-file, are not the reason): RAII has cleaned up; re-raise it
    if (int sig = Runtime::get().signal_received()) {
        fprintf(stderr, "[rgfa-zip] stopped by signal %d (%s)%s\n", sig, strsignal(sig),
                g_committed.load() ? " after the outputs were written" : "; nothing written");
        fflush(stderr);
        reraise(sig);
    }
    if (!msg.empty()) fputs(msg.c_str(), stderr);
    return code;
}
