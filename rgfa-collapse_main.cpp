/*
  rgfa-collapse: find inverted (and optionally duplicated) alleles that minigraph stored as
  novel sequence in an rGFA, and rewrite them so the haplotype traverses the reference instead.

  minigraph sometimes cannot align a haplotype's inverted copy of a region back to the reference
  copy, and stores it as one or more alt nodes.  Downstream everything treats it as a large
  insertion: cactus-align cannot merge two separate backbone nodes, and vcfbub flattens the giant
  allele, so the inversion never reaches the VCF.  Collapsing it to a proper inversion edge fixes
  the representation at the root, and minigraph's mapper reuses such edges rather than
  re-creating the redundant node.

  METHOD.  Sites come from vg's snarl decomposition (a snarl is a pair of nodes with no way in or
  out except through them, so it is a self-contained unit of variation).  Within each snarl we
  build the reference traversal and one traversal per connected component of alt nodes, then
  align each alt traversal against the reference traversal.  A minus-strand block inside that
  alignment is an inversion.

  Aligning whole TRAVERSALS rather than individual alt nodes is what makes nested cases visible:
  an inversion occupying 36 kb of a 98 kb snarl scores 37% by any node- or snarl-level coverage
  test and is rejected, but shows plainly as an internal minus-strand block.

  Snarls are read from `vg snarls -n -P <ref> | vg view -Rj -`.

  Alignment is delegated to minimap2 (subprocess).
*/

#include <unistd.h>
#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <atomic>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>

#include "gfakluge.hpp"

using namespace std;

// ---------------------------------------------------------------- data model

struct Node {
    string name, seq, sn;
    int64_t len = 0, rank = -1, so = -1;
    vector<string> raw_tags;
    bool deleted = false;
};

struct Edge {
    string from, to;
    bool from_fwd = true, to_fwd = true;
    string overlap = "0M";
    vector<string> raw_tags;
    bool deleted = false;
    int64_t sr_rank = -1;          // set on synthesised inversion links
};

struct Snarl { string L, R; bool nested = false; };

// one alt allele of a snarl: the nodes, and where each sits in the traversal sequence
struct Trav {
    size_t snarl_i = 0;
    vector<string> nodes;
    vector<pair<int64_t,int64_t>> spans;   // [start,end) in traversal coords, parallel to nodes
    string seq;
};

struct Call {
    string sn;                     // reference sequence the snarl sits on
    int64_t ref_start = 0, ref_end = 0;
    bool inversion = false;
    vector<string> nodes;          // alt nodes the block touches, in traversal order
    string haps;
    int64_t block_bp = 0, node_bp = 0;
    double ident = 0;
    string L, R;                   // snarl boundaries
    size_t snarl_i = 0, trav_i = 0;
    // the alt side in traversal coordinates, so it can be cut at the block boundary: the whole
    // ordered component, each node's [start,end) in the traversal, and the block's extent
    vector<string> trav_nodes;
    vector<pair<int64_t,int64_t>> trav_spans;
    int64_t q_start = 0, q_end = 0;
    // The target is normally the reference, located by coordinate.  For a non-reference collapse
    // it is another traversal's nodes -- a "representative" -- so carry them and the block's
    // extent in that traversal's coordinates.  Empty target_nodes means the reference.
    vector<string> target_nodes;
    vector<pair<int64_t,int64_t>> target_spans;
    int64_t t_start = 0, t_end = 0;
};

// RAII temp directory
struct TmpDir {
    string path;
    TmpDir() {
        char t[] = "/tmp/rgfacol_XXXXXX";
        if (!mkdtemp(t)) { cerr << "[rgfa-collapse] cannot create temp dir\n"; exit(1); }
        path = t;
    }
    ~TmpDir() {
        // remove the files we know we made; leave the directory if anything else is in it
        for (const char* f : {"ref", "alt", "paf"}) {
            for (int i = 0; i < 4096; ++i) {
                string p = path + "/" + f + to_string(i) + ".fa";
                if (unlink(p.c_str()) != 0) break;
            }
        }
        rmdir(path.c_str());
    }
};

// ---------------------------------------------------------------- options

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <in.gfa> <snarls.json> > out.gfa\n"
         << "Collapse inverted / duplicated alleles that minigraph stored as novel sequence.\n\n"
         << "  snarls.json from:  vg snarls -n -P <ref> -t N in.gfa | vg view -Rj - > snarls.json\n\n"
         << "  -b, --min-block N      minimum aligned block to call [5000]\n"
         << "  -i, --min-ident F      minimum block identity [0.95]\n"
         << "  -a, --min-alt F        block must be >= F alt sequence, not reference flank [0.5]\n"
         << "  -C, --min-cover F      skip alleles whose block covers < F of their alt nodes.\n"
         << "                         Not needed now that alt nodes are split at the block\n"
         << "                         boundary -- only the inverted part is removed [0]\n"
         << "  -A, --alt-rounds N     after the reference pass, collapse traversals that match no\n"
         << "                         reference onto each other: N greedy rounds, each taking the\n"
         << "                         longest remaining traversal in a snarl as representative [0]\n"
         << "  -D, --duplications     also collapse forward duplicates (UNSAFE without a\n"
         << "                         copy-number gate: deletes real tandem expansions) [off]\n"
         << "  -c, --max-components N skip snarls with more than N alt components [64]\n"
         << "  -L, --max-traversal N  skip traversals longer than N bp [5000000]\n"
         << "  -k, --chunk N          snarls per minimap2 invocation [300]\n"
         << "  -j, --jobs N           concurrent minimap2 invocations [4]\n"
         << "  -t, --threads N        threads per minimap2 [2]\n"
         << "  -N, --mm-secondary N   minimap2 -N: alignments retained per query.  A site can\n"
         << "                         lose its self-hit when homologous sites share an index [50]\n"
         << "  -x, --mm-preset STR    minimap2 preset [asm20]\n"
         << "  -m, --minimap2 PATH    minimap2 binary [minimap2]\n"
         << "  -r, --report FILE      write a TSV of every call\n"
         << "  -d, --detect-only      report only, do not rewrite\n";
}

int main(int argc, char** argv) {
    int64_t min_block = 5000, max_trav = 5000000, chunk = 300;
    int max_comp = 64, n_jobs = 4, threads = 2, mm_N = 50, alt_rounds = 0;
    double min_ident = 0.95, min_alt = 0.5, min_cover = 0.0;
    bool do_dup = false, detect_only = false;
    string mm2 = "minimap2", preset = "asm20", report;

    int c;
    while (true) {
        static struct option lo[] = {
            {"min-block",required_argument,0,'b'},{"min-ident",required_argument,0,'i'},
            {"min-alt",required_argument,0,'a'},{"min-cover",required_argument,0,'C'},
            {"alt-rounds",required_argument,0,'A'},
            {"duplications",no_argument,0,'D'},
            {"max-components",required_argument,0,'c'},{"max-traversal",required_argument,0,'L'},
            {"chunk",required_argument,0,'k'},{"jobs",required_argument,0,'j'},
            {"threads",required_argument,0,'t'},{"mm-secondary",required_argument,0,'N'},
            {"mm-preset",required_argument,0,'x'},
            {"minimap2",required_argument,0,'m'},{"report",required_argument,0,'r'},
            {"detect-only",no_argument,0,'d'},{"help",no_argument,0,'h'},{0,0,0,0}};
        c = getopt_long(argc, argv, "b:i:a:C:A:Dc:L:k:j:t:N:x:m:r:dh", lo, 0);
        if (c == -1) break;
        switch (c) {
            case 'b': min_block = stol(optarg); break;
            case 'i': min_ident = stod(optarg); break;
            case 'a': min_alt = stod(optarg); break;
            case 'C': min_cover = stod(optarg); break;
            case 'A': alt_rounds = stoi(optarg); break;
            case 'D': do_dup = true; break;
            case 'c': max_comp = stoi(optarg); break;
            case 'L': max_trav = stol(optarg); break;
            case 'k': chunk = stol(optarg); break;
            case 'j': n_jobs = stoi(optarg); break;
            case 't': threads = stoi(optarg); break;
            case 'N': mm_N = stoi(optarg); break;
            case 'x': preset = optarg; break;
            case 'm': mm2 = optarg; break;
            case 'r': report = optarg; break;
            case 'd': detect_only = true; break;
            case 'h': default: help(argv); return c == 'h' ? 0 : 1;
        }
    }
    if (optind + 1 >= argc) { help(argv); return 1; }
    string gfa_path = argv[optind], snarl_path = argv[optind + 1];

    // ------------------------------------------------------------ read rGFA

    vector<Node> nodes;
    unordered_map<string,size_t> idx;
    vector<Edge> edges;
    function<void(const gfak::sequence_elem&)> vs = [&](const gfak::sequence_elem& s) {
        Node n; n.name = s.name; n.seq = s.sequence; n.len = s.sequence.size();
        for (const gfak::opt_elem& o : s.opt_fields) {
            if (o.key == "SR") n.rank = stol(o.val);
            else if (o.key == "SO") n.so = stol(o.val);
            else if (o.key == "SN") n.sn = o.val;
            n.raw_tags.push_back(o.key + ":" + o.type + ":" + o.val);
        }
        idx[n.name] = nodes.size(); nodes.push_back(move(n));
    };
    function<void(const gfak::edge_elem&)> ve = [&](const gfak::edge_elem& e) {
        Edge ed; ed.from = e.source_name; ed.to = e.sink_name;
        ed.from_fwd = e.source_orientation_forward; ed.to_fwd = e.sink_orientation_forward;
        for (auto& kv : e.tags) ed.raw_tags.push_back(kv.second.to_string());
        edges.push_back(move(ed));
    };
    gfak::GFAKluge kluge;
    kluge.for_each_sequence_line_in_file(gfa_path.c_str(), vs);
    kluge.for_each_edge_line_in_file(gfa_path.c_str(), ve);
    cerr << "[rgfa-collapse] " << nodes.size() << " nodes, " << edges.size() << " edges\n";

    auto NI = [&](const string& n)->const Node* {
        auto it = idx.find(n); return it == idx.end() ? nullptr : &nodes[it->second];
    };
    auto is_ref = [&](const string& n)->bool { auto p = NI(n); return p && p->rank == 0; };

    // plain directed adjacency, ignoring link orientation: an alt chain carrying a reverse link
    // is invisible to orientation-filtered adjacency, and those are the inversions we want
    unordered_map<string, vector<string>> succ, pred;
    unordered_map<string, vector<size_t>> node_edges;
    for (size_t i = 0; i < edges.size(); ++i) {
        succ[edges[i].from].push_back(edges[i].to);
        pred[edges[i].to].push_back(edges[i].from);
        node_edges[edges[i].from].push_back(i);
        if (edges[i].to != edges[i].from) node_edges[edges[i].to].push_back(i);
    }

    // reference nodes indexed per SN (SO is an offset within a reference sequence, so a single
    // SO-sorted list interleaves chromosomes on a whole-genome rGFA)
    map<string, vector<pair<int64_t,string>>> ref_by_sn;
    for (const Node& n : nodes)
        if (n.rank == 0 && n.so >= 0 && !n.sn.empty()) ref_by_sn[n.sn].push_back({n.so, n.name});
    for (auto& kv : ref_by_sn) sort(kv.second.begin(), kv.second.end());

    // ------------------------------------------------------------ read snarls

    // vg view -Rj output: protobuf JSON, keys sorted alphabetically.  With -n the snarl's own
    // boundaries carry "name" while the optional "parent" carries "node_id", so every line has
    // exactly two "name" fields -- end first, then start -- and "parent" marks a nested snarl.
    // Verified 7466/7466 lines on HPRC chr15.  Assert it rather than trust it.
    vector<Snarl> snarls;
    {
        ifstream sf(snarl_path);
        if (!sf) { cerr << "[rgfa-collapse] cannot open " << snarl_path << "\n"; return 1; }
        string line; int64_t lineno = 0;
        while (getline(sf, line)) {
            ++lineno;
            if (line.empty()) continue;
            vector<string> names;
            for (size_t p = line.find("\"name\":"); p != string::npos; p = line.find("\"name\":", p + 1)) {
                size_t a = line.find('"', p + 7);
                if (a == string::npos) break;
                size_t b = line.find('"', a + 1);
                if (b == string::npos) break;
                names.push_back(line.substr(a + 1, b - a - 1));
            }
            if (names.size() != 2) {
                cerr << "[rgfa-collapse] error: " << snarl_path << " line " << lineno
                     << " has " << names.size() << " \"name\" fields, expected 2.\n"
                     << "  Snarls must come from: vg snarls -n ... | vg view -Rj -\n";
                return 1;
            }
            Snarl s; s.R = names[0]; s.L = names[1];         // end first, start second
            s.nested = line.find("\"parent\":") != string::npos;
            snarls.push_back(move(s));
        }
    }
    size_t n_top = 0; for (auto& s : snarls) if (!s.nested) ++n_top;
    cerr << "[rgfa-collapse] " << snarls.size() << " snarls, " << n_top << " top-level\n";

    // ------------------------------------------------------------ build sites

    struct TravMeta { vector<string> nodes; int64_t pa = 0, pb = 0; };
    struct Job {
        size_t snarl_i = 0; string L, R, sn;
        int64_t lo = 0, hi = 0;
        vector<string> refp; int64_t reflen = 0;
        vector<TravMeta> travs;
    };
    vector<Job> jobs;
    atomic<int64_t> sk_bound(0), sk_diffseq(0), sk_nospan(0), sk_big(0),
                    sk_noalt(0), sk_noref(0), sk_manycomp(0), sk_longtrav(0);
    {
        size_t nthread = max(1u, thread::hardware_concurrency());
        vector<vector<Job>> per(nthread);
        atomic<size_t> next(0);
        vector<thread> pool;
        for (size_t t = 0; t < nthread; ++t) pool.emplace_back([&, t]() {
            while (true) {
                size_t si = next++;
                if (si >= snarls.size()) break;
                const Snarl& s = snarls[si];
                if (s.nested) continue;                       // nested sites are inside a top-level one
                const Node* LN = NI(s.L); const Node* RN = NI(s.R);
                if (!LN || !RN || LN->rank != 0 || RN->rank != 0) { ++sk_bound; continue; }
                if (LN->sn != RN->sn) { ++sk_diffseq; continue; }
                const Node* A = LN->so < RN->so ? LN : RN;
                const Node* B = LN->so < RN->so ? RN : LN;
                int64_t lo = A->so + A->len, hi = B->so;
                if (hi <= lo) { ++sk_nospan; continue; }

                // nodes strictly inside the snarl
                vector<string> inside; unordered_set<string> seen{s.L, s.R};
                vector<string> st;
                { auto it = succ.find(s.L);
                  if (it != succ.end()) for (auto& x : it->second)
                      if (x != s.R && seen.insert(x).second) st.push_back(x); }
                bool too_big = false;
                while (!st.empty()) {
                    if (inside.size() > 200000) { too_big = true; break; }
                    string x = st.back(); st.pop_back(); inside.push_back(x);
                    for (const auto* m : {&succ, &pred}) {
                        auto it = m->find(x);
                        if (it == m->end()) continue;
                        for (auto& y : it->second) if (seen.insert(y).second) st.push_back(y);
                    }
                }
                if (too_big) { ++sk_big; continue; }
                vector<string> alts;
                for (auto& n : inside) { auto p = NI(n); if (p && p->rank > 0) alts.push_back(n); }
                if (alts.empty()) { ++sk_noalt; continue; }

                Job j; j.snarl_i = si; j.L = A->name; j.R = B->name; j.sn = A->sn; j.lo = lo; j.hi = hi;
                auto& rv = ref_by_sn[A->sn];
                auto rit = lower_bound(rv.begin(), rv.end(), make_pair(lo, string()));
                for (; rit != rv.end() && rit->first < hi; ++rit) {
                    auto p = NI(rit->second);
                    if (p && rit->first >= lo && rit->first + p->len <= hi) {
                        j.refp.push_back(rit->second); j.reflen += p->len;
                    }
                }
                if (j.refp.empty()) { ++sk_noref; continue; }
                unordered_map<string,int64_t> roff;
                { int64_t acc = 0; for (auto& n : j.refp) { roff[n] = acc; acc += NI(n)->len; } }

                // connected components among the alt nodes
                unordered_set<string> aset(alts.begin(), alts.end()), done;
                vector<vector<string>> comps;
                for (auto& n : alts) {
                    if (!done.insert(n).second) continue;
                    vector<string> comp, s2{n};
                    while (!s2.empty()) {
                        string x = s2.back(); s2.pop_back(); comp.push_back(x);
                        for (const auto* m : {&succ, &pred}) {
                            auto it = m->find(x);
                            if (it == m->end()) continue;
                            for (auto& y : it->second)
                                if (aset.count(y) && done.insert(y).second) s2.push_back(y);
                        }
                    }
                    comps.push_back(move(comp));
                }
                if ((int)comps.size() > max_comp) { ++sk_manycomp; continue; }

                for (auto& comp : comps) {
                    // order the component along the graph, best effort
                    unordered_set<string> cs(comp.begin(), comp.end());
                    string start = comp[0];
                    for (auto& n : comp) {
                        bool head = true;
                        auto it = pred.find(n);
                        if (it != pred.end()) for (auto& p : it->second) if (cs.count(p)) { head = false; break; }
                        if (head) { start = n; break; }
                    }
                    vector<string> oc; unordered_set<string> vis;
                    string cur = start;
                    while (!cur.empty() && vis.insert(cur).second) {
                        oc.push_back(cur);
                        string nx;
                        auto it = succ.find(cur);
                        if (it != succ.end()) for (auto& y : it->second)
                            if (cs.count(y) && !vis.count(y)) { nx = y; break; }
                        cur = nx;
                    }
                    for (auto& n : comp) if (!vis.count(n)) oc.push_back(n);

                    // where does the component attach to the reference?
                    int64_t pa = 0, pb = j.reflen; bool haveA = false, haveB = false;
                    for (auto& n : oc) {
                        auto it = pred.find(n);
                        if (it != pred.end()) for (auto& p : it->second) {
                            auto r = roff.find(p);
                            if (r != roff.end()) { int64_t e = r->second + NI(p)->len;
                                pa = haveA ? max(pa, e) : e; haveA = true; }
                        }
                        auto it2 = succ.find(n);
                        if (it2 != succ.end()) for (auto& q : it2->second) {
                            auto r = roff.find(q);
                            if (r != roff.end()) { pb = haveB ? min(pb, r->second) : r->second; haveB = true; }
                        }
                    }
                    if (!haveA) pa = 0;
                    if (!haveB) pb = j.reflen;
                    if (pb < pa) { pa = 0; pb = j.reflen; }
                    int64_t body = 0; for (auto& n : oc) body += NI(n)->len;
                    if (pa + body + (j.reflen - pb) > max_trav) { ++sk_longtrav; continue; }
                    TravMeta tm; tm.nodes = move(oc); tm.pa = pa; tm.pb = pb;
                    j.travs.push_back(move(tm));
                }
                if (!j.travs.empty()) per[t].push_back(move(j));
            }
        });
        for (auto& th : pool) th.join();
        for (auto& v : per) for (auto& j : v) jobs.push_back(move(j));
        // Gathering from per-thread vectors is completion-ordered, so chunk composition -- and
        // therefore which alignments minimap2 reports -- would vary run to run.  Sort so the
        // output is reproducible.
        sort(jobs.begin(), jobs.end(), [](const Job& a, const Job& b) { return a.snarl_i < b.snarl_i; });
    }
    cerr << "[rgfa-collapse] " << jobs.size() << " site(s) to align; skipped:"
         << " boundary=" << sk_bound << " diff-seq=" << sk_diffseq << " no-span=" << sk_nospan
         << " too-big=" << sk_big << " no-alt=" << sk_noalt << " no-ref=" << sk_noref
         << " many-components=" << sk_manycomp << " long-traversal=" << sk_longtrav << "\n";

    // ------------------------------------------------------------ align, chunked and parallel

    TmpDir tmp;
    vector<Call> calls;
    mutex call_mx;
    {
        size_t nchunk = (jobs.size() + chunk - 1) / max<int64_t>(1, chunk);
        atomic<size_t> next(0); atomic<int64_t> done(0);
        vector<thread> pool;
        int J = max(1, min<int>(n_jobs, (int)max<size_t>(1, nchunk)));
        for (int t = 0; t < J; ++t) pool.emplace_back([&, t]() {
            while (true) {
                size_t ci = next++;
                if (ci >= nchunk) break;
                size_t b0 = ci * chunk, b1 = min(jobs.size(), b0 + (size_t)chunk);
                string rf = tmp.path + "/ref" + to_string(t) + ".fa";
                string qf = tmp.path + "/alt" + to_string(t) + ".fa";
                // A single global minimap2 index does not work: snarls are homologous to one
                // another, so with many targets a query's own site is often not among the
                // reported hits and a same-site filter then discards everything.  Each batch is
                // aligned against only its own references.
                {
                    ofstream tf(rf), qfs(qf);
                    for (size_t i = b0; i < b1; ++i) {
                        const Job& j = jobs[i];
                        string rs; rs.reserve(j.reflen);
                        for (auto& n : j.refp) rs += NI(n)->seq;
                        tf << ">r" << i << "\n" << rs << "\n";
                        for (size_t k = 0; k < j.travs.size(); ++k) {
                            const TravMeta& tm = j.travs[k];
                            string body; for (auto& n : tm.nodes) body += NI(n)->seq;
                            // reference prefix + the alt component + reference suffix
                            string ts = rs.substr(0, tm.pa) + body + rs.substr(tm.pb);
                            qfs << ">q" << i << "_" << k << "\n" << ts << "\n";
                        }
                    }
                }
                stringstream cmd;
                cmd << mm2 << " -cx " << preset << " -t " << threads << " -N 50 -p 0.01 "
                    << rf << " " << qf << " 2>/dev/null";
                FILE* pf = popen(cmd.str().c_str(), "r");
                if (!pf) { cerr << "[rgfa-collapse] error: cannot run minimap2\n"; exit(1); }
                char* line = nullptr; size_t cap = 0;
                vector<Call> local;
                while (getline(&line, &cap, pf) > 0) {
                    stringstream ss(line);
                    string q, st, tn; int64_t ql, qs, qe, tl, ts_, te, nm, al, mq;
                    if (!(ss >> q >> ql >> qs >> qe >> st >> tn >> tl >> ts_ >> te >> nm >> al >> mq)) continue;
                    if (q.empty() || q[0] != 'q') continue;
                    size_t us = q.find('_');
                    if (us == string::npos) continue;
                    size_t ji = stoul(q.substr(1, us - 1));
                    size_t ki = stoul(q.substr(us + 1));
                    if (tn != "r" + to_string(ji)) continue;          // must be its own site
                    if (ji >= jobs.size()) continue;
                    const Job& j = jobs[ji];
                    if (ki >= j.travs.size()) continue;
                    int64_t blk = qe - qs;
                    if (blk < min_block) continue;
                    double id = al ? (double)nm / al : 0.0;
                    if (id < min_ident) continue;
                    // the block must be mostly ALT sequence: a traversal aligns to its own
                    // reference by construction, so otherwise every traversal yields a large
                    // spurious forward "duplication" that merely clips an alt node at the edge
                    const TravMeta& tm = j.travs[ki];
                    int64_t off = tm.pa, altov = 0; vector<string> hit;
                    for (auto& n : tm.nodes) {
                        int64_t a = off, b = off + NI(n)->len; off = b;
                        int64_t ov = min(qe, b) - max(qs, a);
                        if (ov > 0) { altov += ov; hit.push_back(n); }
                    }
                    if (hit.empty() || (double)altov / blk < min_alt) continue;
                    Call cl;
                    { int64_t o2 = tm.pa;
                      for (auto& n : tm.nodes) { cl.trav_spans.push_back({o2, o2 + NI(n)->len}); o2 += NI(n)->len; } }
                    cl.trav_nodes = tm.nodes; cl.q_start = qs; cl.q_end = qe;
                    cl.sn = j.sn; cl.ref_start = j.lo + ts_; cl.ref_end = j.lo + te;
                    cl.inversion = (st == "-"); cl.nodes = hit; cl.block_bp = blk; cl.ident = id;
                    cl.L = j.L; cl.R = j.R; cl.snarl_i = ji; cl.trav_i = ki;
                    set<string> hs;
                    for (auto& n : hit) { cl.node_bp += NI(n)->len; if (!NI(n)->sn.empty()) hs.insert(NI(n)->sn); }
                    for (auto& h : hs) { if (!cl.haps.empty()) cl.haps += ","; cl.haps += h; }
                    local.push_back(move(cl));
                }
                free(line);
                int rc = pclose(pf);
                if (rc != 0) { cerr << "[rgfa-collapse] error: minimap2 exited " << rc << "\n"; exit(1); }
                { lock_guard<mutex> g(call_mx);
                  for (auto& c2 : local) calls.push_back(move(c2)); }
                int64_t d = ++done;
                if (d % 50 == 0 || d == (int64_t)nchunk)
                    cerr << "[rgfa-collapse]   " << d << "/" << nchunk << " chunk(s)\n";
            }
        });
        for (auto& th : pool) th.join();
    }

    // ------------------------------------------------------------ non-reference rounds
    //
    // A novel insertion carried by many haplotypes is stored as many near-identical alt nodes,
    // none matching the reference, so the pass above cannot see it (on HPRC chr15, 2.02 Mb of
    // non-reference sequence >=1 kb has no reference match at all; on chr21, 52%).  Compare those
    // traversals to EACH OTHER instead: greedily take the longest one in a snarl as the
    // representative, collapse whatever matches it, and repeat on what is left.  All-vs-all would
    // be quadratic; this is one batched pass per round over a shrinking set.
    //
    // Component BODIES are compared, not whole traversals -- traversals share reference flanks
    // and would align through those whether or not the alleles are redundant.
    if (alt_rounds > 0) {
        set<pair<size_t,size_t>> resolved;          // (site, traversal) already collapsed or kept
        for (auto& c : calls) resolved.insert({c.snarl_i, c.trav_i});
        auto body_len = [&](size_t i, size_t j) {
            int64_t L = 0; for (auto& n : jobs[i].travs[j].nodes) { auto q = NI(n); if (q) L += q->len; } return L;
        };
        for (int round = 1; round <= alt_rounds; ++round) {
            struct AltJob { size_t i, repr; vector<size_t> q; };
            vector<AltJob> aj;
            for (size_t i = 0; i < jobs.size(); ++i) {
                vector<size_t> open;
                for (size_t j = 0; j < jobs[i].travs.size(); ++j)
                    if (!resolved.count({i, j}) && body_len(i, j) >= min_block) open.push_back(j);
                if (open.size() < 2) continue;
                // longest wins; ties broken on the first node's name so rounds are reproducible
                size_t best = open[0];
                for (size_t j : open) {
                    int64_t a = body_len(i, j), b = body_len(i, best);
                    if (a > b || (a == b && jobs[i].travs[j].nodes[0] < jobs[i].travs[best].nodes[0]))
                        best = j;
                }
                AltJob a; a.i = i; a.repr = best;
                for (size_t j : open) if (j != best) a.q.push_back(j);
                aj.push_back(move(a));
            }
            if (aj.empty()) { cerr << "[rgfa-collapse] round " << round << ": nothing left to compare\n"; break; }

            size_t nchunk = (aj.size() + chunk - 1) / max<int64_t>(1, chunk);
            atomic<size_t> next(0);
            vector<thread> pool;
            vector<Call> round_calls;
            mutex rmx;
            int J = max(1, min<int>(n_jobs, (int)max<size_t>(1, nchunk)));
            for (int t = 0; t < J; ++t) pool.emplace_back([&, t]() {
                while (true) {
                    size_t ci = next++;
                    if (ci >= nchunk) break;
                    size_t b0 = ci * chunk, b1 = min(aj.size(), b0 + (size_t)chunk);
                    string rf = tmp.path + "/ref" + to_string(t) + ".fa";
                    string qf = tmp.path + "/alt" + to_string(t) + ".fa";
                    {
                        ofstream tf(rf), qfs(qf);
                        for (size_t x = b0; x < b1; ++x) {
                            const AltJob& a = aj[x];
                            string rs; for (auto& n : jobs[a.i].travs[a.repr].nodes) rs += NI(n)->seq;
                            tf << ">r" << x << "\n" << rs << "\n";
                            for (size_t qi : a.q) {
                                string qs2; for (auto& n : jobs[a.i].travs[qi].nodes) qs2 += NI(n)->seq;
                                qfs << ">q" << x << "_" << qi << "\n" << qs2 << "\n";
                            }
                        }
                    }
                    stringstream cmd;
                    cmd << mm2 << " -cx " << preset << " -t " << threads << " -N " << mm_N
                        << " -p 0.01 " << rf << " " << qf << " 2>/dev/null";
                    FILE* pf = popen(cmd.str().c_str(), "r");
                    if (!pf) { cerr << "[rgfa-collapse] error: cannot run minimap2\n"; exit(1); }
                    char* line = nullptr; size_t cap = 0;
                    vector<Call> local;
                    while (getline(&line, &cap, pf) > 0) {
                        stringstream ss(line);
                        string q, st, tn; int64_t ql, qs2, qe, tl, ts_, te, nm, al, mq;
                        if (!(ss >> q >> ql >> qs2 >> qe >> st >> tn >> tl >> ts_ >> te >> nm >> al >> mq)) continue;
                        if (q.empty() || q[0] != 'q') continue;
                        size_t us = q.find('_');
                        if (us == string::npos) continue;
                        size_t xi = stoul(q.substr(1, us - 1)), qi = stoul(q.substr(us + 1));
                        if (tn != "r" + to_string(xi) || xi >= aj.size()) continue;
                        const AltJob& a = aj[xi];
                        int64_t blk = qe - qs2;
                        if (blk < min_block) continue;
                        double id = al ? (double)nm / al : 0.0;
                        if (id < min_ident) continue;
                        const Job& jb = jobs[a.i];
                        Call cl;
                        { int64_t o2 = 0;
                          for (auto& n : jb.travs[qi].nodes) { cl.trav_spans.push_back({o2, o2 + NI(n)->len}); o2 += NI(n)->len; }
                          o2 = 0;
                          for (auto& n : jb.travs[a.repr].nodes) { cl.target_spans.push_back({o2, o2 + NI(n)->len}); o2 += NI(n)->len; } }
                        cl.trav_nodes = jb.travs[qi].nodes;
                        cl.target_nodes = jb.travs[a.repr].nodes;
                        cl.q_start = qs2; cl.q_end = qe; cl.t_start = ts_; cl.t_end = te;
                        cl.sn = jb.sn; cl.ref_start = jb.lo; cl.ref_end = jb.hi;
                        cl.inversion = (st == "-");
                        cl.L = jb.L; cl.R = jb.R; cl.snarl_i = a.i; cl.trav_i = qi;
                        cl.block_bp = blk; cl.ident = id;
                        int64_t altov = 0;
                        for (auto& pr : cl.trav_spans) {
                            int64_t ov = min(qe, pr.second) - max(qs2, pr.first);
                            if (ov > 0) altov += ov;
                        }
                        if ((double)altov / blk < min_alt) continue;
                        for (size_t k = 0; k < cl.trav_nodes.size(); ++k)
                            if (cl.trav_spans[k].first < qe && qs2 < cl.trav_spans[k].second)
                                cl.nodes.push_back(cl.trav_nodes[k]);
                        if (cl.nodes.empty()) continue;
                        set<string> hs;
                        for (auto& n : cl.nodes) { cl.node_bp += NI(n)->len; if (!NI(n)->sn.empty()) hs.insert(NI(n)->sn); }
                        for (auto& h : hs) { if (!cl.haps.empty()) cl.haps += ","; cl.haps += h; }
                        local.push_back(move(cl));
                    }
                    free(line);
                    int rc = pclose(pf);
                    if (rc != 0) { cerr << "[rgfa-collapse] error: minimap2 exited " << rc << "\n"; exit(1); }
                    { lock_guard<mutex> g(rmx); for (auto& c2 : local) round_calls.push_back(move(c2)); }
                }
            });
            for (auto& th : pool) th.join();
            // representatives are retired whether or not anything matched them, so a later round
            // cannot pick one again or treat it as a query
            for (auto& a : aj) resolved.insert({a.i, a.repr});
            int64_t added = 0, added_bp = 0;
            for (auto& c : round_calls) {
                if (!resolved.insert({c.snarl_i, c.trav_i}).second) continue;   // one call per traversal
                added_bp += c.node_bp; ++added;
                calls.push_back(move(c));
            }
            cerr << "[rgfa-collapse] round " << round << ": " << aj.size() << " site(s), "
                 << added << " traversal(s) collapsed onto a non-reference representative, "
                 << added_bp << " bp\n";
            if (added == 0) break;
        }
    }

    cerr << "[rgfa-collapse] " << calls.size() << " raw call(s)\n";

    // ------------------------------------------------------------ dedup

    // Nested snarls and per-haplotype alleles both produce several calls over the same reference
    // interval.  Merge overlapping intervals within a reference sequence so a locus is repaired
    // once, taking the union of the alt nodes each call named.
    struct Site {
        string sn, L, R; int64_t start = 0, end = 0; bool inversion = false;
        vector<string> nodes; set<string> haps; double ident = 0; int64_t block_bp = 0, node_bp = 0;
    };
    vector<Site> sites;
    {
        sort(calls.begin(), calls.end(), [](const Call& a, const Call& b) {
            if (a.sn != b.sn) return a.sn < b.sn;
            if (a.inversion != b.inversion) return a.inversion > b.inversion;
            return a.ref_start < b.ref_start;
        });
        for (auto& c : calls) {
            if (!sites.empty() && sites.back().sn == c.sn && sites.back().inversion == c.inversion
                && c.ref_start <= sites.back().end) {
                Site& s2 = sites.back();
                s2.end = max(s2.end, c.ref_end);
                for (auto& n : c.nodes) if (find(s2.nodes.begin(), s2.nodes.end(), n) == s2.nodes.end())
                    s2.nodes.push_back(n);
                if (!c.haps.empty()) { stringstream hs(c.haps); string h;
                    while (getline(hs, h, ',')) s2.haps.insert(h); }
                s2.ident = max(s2.ident, c.ident);
                s2.block_bp = max(s2.block_bp, c.block_bp);
            } else {
                Site s2; s2.sn = c.sn; s2.L = c.L; s2.R = c.R; s2.start = c.ref_start; s2.end = c.ref_end;
                s2.inversion = c.inversion; s2.nodes = c.nodes; s2.ident = c.ident; s2.block_bp = c.block_bp;
                if (!c.haps.empty()) { stringstream hs(c.haps); string h;
                    while (getline(hs, h, ',')) s2.haps.insert(h); }
                sites.push_back(move(s2));
            }
        }
        for (auto& s2 : sites) { s2.node_bp = 0; for (auto& n : s2.nodes) s2.node_bp += NI(n)->len; }
    }
    int64_t n_inv = 0, n_dup = 0, bp_inv = 0, bp_dup = 0;
    for (auto& s2 : sites) { if (s2.inversion) { ++n_inv; bp_inv += s2.node_bp; }
                             else { ++n_dup; bp_dup += s2.node_bp; } }
    cerr << "[rgfa-collapse] " << sites.size() << " site(s) after merge: "
         << n_inv << " inversion(s) / " << bp_inv << " bp, "
         << n_dup << " duplication(s) / " << bp_dup << " bp\n";

    if (!report.empty()) {
        ofstream rf(report);
        rf << "#seq\tstart\tend\ttype\thaplotypes\tnodes\tnode_bp\tblock_bp\tident\tsnarl_L\tsnarl_R\tnode_ids\n";
        for (auto& s2 : sites) {
            rf << s2.sn << "\t" << s2.start << "\t" << s2.end << "\t"
               << (s2.inversion ? "INV" : "DUP") << "\t";
            bool first = true; for (auto& h : s2.haps) { if (!first) rf << ","; rf << h; first = false; }
            rf << "\t" << s2.nodes.size() << "\t" << s2.node_bp << "\t" << s2.block_bp << "\t"
               << s2.ident << "\t" << s2.L << "\t" << s2.R << "\t";
            first = true; for (auto& n : s2.nodes) { if (!first) rf << ","; rf << n; first = false; }
            rf << "\n";
        }
        cerr << "[rgfa-collapse] wrote " << report << "\n";
    }
    if (detect_only) return 0;

    // ------------------------------------------------------------ repair

    // Split nodes so an inverted allele can be removed exactly.
    //
    // Two reasons.  On the REFERENCE side an inversion need not align to whole nodes -- on HPRC
    // chr15 one sits entirely inside a single 17,712 bp node, leaving no chain to reverse.  On
    // the ALT side the node holding the inversion often carries sequence beyond it: a yeast
    // chrXIV node is 23,061 bp for a 5,914 bp inverted block, and on chr15 seven of thirteen
    // sites would have had 243,585 bp of genuine novel sequence deleted with the inversion.
    //
    // After cutting both sides at the block boundaries the repair is exact:
    //
    //     before:  L -> [A_pre | A_inv | A_post] -> M   and   L -> ...R... -> M
    //     after:   L -> A_pre -> reverse(R) -> A_post -> M
    //
    // A_pre / A_post attach to the reversed reference chain; when the block covers the whole alt
    // component they are empty and the links come from the reference flanks instead.
    map<string, set<int64_t>> cuts;
    for (auto& c : calls) {
        if (!c.inversion && !do_dup) continue;
        auto& rv = ref_by_sn[c.sn];
        for (int64_t bp : {c.ref_start, c.ref_end}) {
            auto it = upper_bound(rv.begin(), rv.end(), make_pair(bp, string("\xff")));
            if (it == rv.begin()) continue;
            --it;
            const Node* p = NI(it->second);
            if (p && bp > p->so && bp < p->so + p->len) cuts[it->second].insert(bp - p->so);
        }
        for (size_t k = 0; k < c.trav_nodes.size(); ++k) {
            int64_t a0 = c.trav_spans[k].first, b0 = c.trav_spans[k].second;
            for (int64_t bp : {c.q_start, c.q_end})
                if (bp > a0 && bp < b0) cuts[c.trav_nodes[k]].insert(bp - a0);
        }
        // a non-reference target is cut the same way, so the wired chain is exactly the part
        // the query matched rather than the representative's whole component
        for (size_t k = 0; k < c.target_nodes.size(); ++k) {
            int64_t a0 = c.target_spans[k].first, b0 = c.target_spans[k].second;
            for (int64_t bp : {c.t_start, c.t_end})
                if (bp > a0 && bp < b0) cuts[c.target_nodes[k]].insert(bp - a0);
        }
    }
    // name -> the pieces it became, with each piece's offset range in the original
    map<string, vector<pair<int64_t,string>>> pieces;
    {
        int64_t n_split = 0, n_pieces = 0;
        for (auto& kv : cuts) {
            auto i2 = idx.find(kv.first);
            if (i2 == idx.end() || nodes[i2->second].deleted) continue;
            Node orig = nodes[i2->second];
            vector<int64_t> off(kv.second.begin(), kv.second.end());
            off.insert(off.begin(), 0); off.push_back(orig.len);
            vector<string> pl;
            for (size_t k = 0; k + 1 < off.size(); ++k) {
                Node pn;
                pn.name = orig.name + "." + to_string(k + 1);
                pn.seq = orig.seq.substr(off[k], off[k + 1] - off[k]);
                pn.len = pn.seq.size();
                pn.rank = orig.rank; pn.sn = orig.sn;
                pn.so = orig.so >= 0 ? orig.so + off[k] : -1;
                pn.raw_tags = {"LN:i:" + to_string(pn.len)};
                if (!pn.sn.empty()) pn.raw_tags.push_back("SN:Z:" + pn.sn);
                if (pn.so >= 0) pn.raw_tags.push_back("SO:i:" + to_string(pn.so));
                pn.raw_tags.push_back("SR:i:" + to_string(pn.rank));
                idx[pn.name] = nodes.size(); nodes.push_back(move(pn));
                pieces[orig.name].push_back({off[k], orig.name + "." + to_string(k + 1)});
                pl.push_back(orig.name + "." + to_string(k + 1));
                ++n_pieces;
            }
            nodes[i2->second].deleted = true;
            for (size_t ei : node_edges[kv.first]) {
                Edge& e = edges[ei];
                if (e.deleted) continue;
                if (e.to == kv.first && e.from == kv.first) { e.deleted = true; continue; }
                // re-register under the new endpoint, or deleting that piece later will not
                // mark this edge deleted and it is left pointing at a removed node
                if (e.to == kv.first)        { e.to = pl.front();  node_edges[e.to].push_back(ei); }
                else if (e.from == kv.first) { e.from = pl.back(); node_edges[e.from].push_back(ei); }
                for (auto& t : e.raw_tags) {
                    if (t.compare(0, 5, "L1:i:") == 0 && NI(e.from)) t = "L1:i:" + to_string(NI(e.from)->len);
                    else if (t.compare(0, 5, "L2:i:") == 0 && NI(e.to)) t = "L2:i:" + to_string(NI(e.to)->len);
                }
            }
            for (size_t k = 0; k + 1 < pl.size(); ++k) {
                Edge e; e.from = pl[k]; e.from_fwd = true; e.to = pl[k + 1]; e.to_fwd = true;
                e.sr_rank = orig.rank; e.raw_tags = {"SR:i:" + to_string(orig.rank)};
                edges.push_back(e);
                node_edges[e.from].push_back(edges.size() - 1);
                node_edges[e.to].push_back(edges.size() - 1);
            }
            ++n_split;
        }
        if (n_split) {
            cerr << "[rgfa-collapse] split " << n_split << " node(s) into " << n_pieces
                 << " piece(s) at block boundaries\n";
            ref_by_sn.clear();
            for (const Node& n : nodes)
                if (!n.deleted && n.rank == 0 && n.so >= 0 && !n.sn.empty())
                    ref_by_sn[n.sn].push_back({n.so, n.name});
            for (auto& kv : ref_by_sn) sort(kv.second.begin(), kv.second.end());
        }
    }
    // expand a node into its pieces, each tagged with its offset within the original
    auto expand = [&](const string& n) {
        vector<pair<int64_t,string>> out;
        auto it = pieces.find(n);
        if (it == pieces.end()) out.push_back({0, n});
        else out = it->second;
        return out;
    };

    vector<Edge> new_edges;
    int64_t did_inv = 0, did_dup = 0, bp_removed = 0, skipped_nochain = 0;
    for (auto& c : calls) {
        if (!c.inversion && !do_dup) continue;
        // alt pieces lying inside the block, in traversal order
        vector<string> del; string before, after;
        for (size_t k = 0; k < c.trav_nodes.size(); ++k) {
            int64_t base = c.trav_spans[k].first;
            for (auto& pr : expand(c.trav_nodes[k])) {
                auto pn = NI(pr.second);
                if (!pn) continue;
                int64_t a0 = base + pr.first, b0 = a0 + pn->len;
                if (a0 >= c.q_start && b0 <= c.q_end) del.push_back(pr.second);
                else if (b0 <= c.q_start) before = pr.second;          // last piece before
                else if (a0 >= c.q_end && after.empty()) after = pr.second;  // first piece after
            }
        }
        if (del.empty()) { ++skipped_nochain; continue; }
        {
            // The wiring is the same for both directions -- route the haplotype through the
            // reference the allele duplicates, forwards or reversed:
            //
            //   inversion:   X -> chain_last(-)   chain_first(-) -> Y
            //   duplication: X -> chain_first(+)  chain_last(+)  -> Y
            //
            // Links are needed in the forward case too once alt nodes are split: without them a
            // surviving A_pre is a dead end and A_post has nothing entering it.  (Deleting the
            // whole alt node needed no links, because the reference path L->R was already there.)
            vector<string> chain;
            if (c.target_nodes.empty()) {
                auto& rv = ref_by_sn[c.sn];
                auto it = lower_bound(rv.begin(), rv.end(), make_pair(c.ref_start, string()));
                for (auto k = it; k != rv.end() && k->first < c.ref_end; ++k) {
                    auto p = NI(k->second);
                    if (p && k->first >= c.ref_start && k->first + p->len <= c.ref_end) chain.push_back(k->second);
                }
            } else {
                // the representative's pieces covered by the matched block, in traversal order
                for (size_t k = 0; k < c.target_nodes.size(); ++k) {
                    int64_t base = c.target_spans[k].first;
                    for (auto& pr : expand(c.target_nodes[k])) {
                        auto pn = NI(pr.second);
                        if (!pn) continue;
                        int64_t a0 = base + pr.first, b0 = a0 + pn->len;
                        if (a0 >= c.t_start && b0 <= c.t_end) chain.push_back(pr.second);
                    }
                }
            }
            if (chain.empty()) { ++skipped_nochain; continue; }
            // entry / exit: the alt pieces either side of the removed run, else the snarl flanks
            string X = before, Y = after;
            if (X.empty()) { auto pl = expand(c.L); X = pl.back().second; }
            if (Y.empty()) { auto pl = expand(c.R); Y = pl.front().second; }
            if (!NI(X) || !NI(Y)) { ++skipped_nochain; continue; }
            int64_t rk = 0; for (auto& n : del) if (NI(n)) rk = max(rk, NI(n)->rank);
            const string& first_hop = c.inversion ? chain.back() : chain.front();
            const string& last_hop  = c.inversion ? chain.front() : chain.back();
            const bool orient = !c.inversion;              // '+' forward, '-' reversed
            Edge e1; e1.from = X; e1.from_fwd = true;   e1.to = first_hop; e1.to_fwd = orient; e1.sr_rank = rk;
            Edge e2; e2.from = last_hop; e2.from_fwd = orient; e2.to = Y;  e2.to_fwd = true;   e2.sr_rank = rk;
            new_edges.push_back(e1); new_edges.push_back(e2);
            if (c.inversion) ++did_inv; else ++did_dup;
        }
        for (auto& n : del) {
            auto i2 = idx.find(n); if (i2 == idx.end()) continue;
            Node& N = nodes[i2->second];
            if (N.deleted) continue;
            N.deleted = true; bp_removed += N.len;
            for (size_t ei : node_edges[n]) edges[ei].deleted = true;
        }
    }
    cerr << "[rgfa-collapse] repaired " << did_inv << " inversion allele(s)";
    if (do_dup) cerr << " and " << did_dup << " duplication(s)";
    cerr << "; " << bp_removed << " bp removed";
    if (skipped_nochain) cerr << "; " << skipped_nochain << " skipped (no chain or no alt piece inside the block)";
    cerr << "\n";

    // Collapsing can orphan alt nodes that only existed inside the redundant sequence -- nested
    // bubbles from other haplotypes that aligned to it.  They are redundant for the same reason,
    // and leaving them behind fragments the graph.
    {
        unordered_map<string, vector<string>> live;
        for (const Edge& e : edges) if (!e.deleted) { live[e.from].push_back(e.to); live[e.to].push_back(e.from); }
        for (const Edge& e : new_edges) { live[e.from].push_back(e.to); live[e.to].push_back(e.from); }
        unordered_set<string> seen; int64_t nc = 0, bc = 0;
        for (const Node& n0 : nodes) {
            if (n0.deleted || seen.count(n0.name)) continue;
            vector<string> comp, st{n0.name}; seen.insert(n0.name); bool anchored = false;
            while (!st.empty()) {
                string x = st.back(); st.pop_back(); comp.push_back(x);
                if (is_ref(x)) anchored = true;
                auto it = live.find(x);
                if (it == live.end()) continue;
                for (auto& y : it->second) {
                    auto yi = idx.find(y);
                    if (yi == idx.end() || nodes[yi->second].deleted) continue;
                    if (seen.insert(y).second) st.push_back(y);
                }
            }
            if (anchored) continue;
            for (auto& x : comp) {
                auto xi = idx.find(x); if (xi == idx.end()) continue;
                Node& N = nodes[xi->second];
                if (N.deleted) continue;
                N.deleted = true; ++nc; bc += N.len;
                for (size_t ei : node_edges[x]) edges[ei].deleted = true;
            }
        }
        if (nc) cerr << "[rgfa-collapse] cascaded " << nc << " orphaned node(s) / " << bc << " bp\n";
    }

    // ------------------------------------------------------------ emit

    { ifstream hf(gfa_path); string line;
      while (getline(hf, line)) { if (!line.empty() && line[0] == 'H') cout << line << "\n"; else break; } }
    for (auto& n : nodes) {
        if (n.deleted) continue;
        cout << "S\t" << n.name << "\t" << n.seq;
        for (auto& t : n.raw_tags) cout << "\t" << t;
        cout << "\n";
    }
    for (auto& e : edges) {
        if (e.deleted) continue;
        cout << "L\t" << e.from << "\t" << (e.from_fwd ? '+' : '-') << "\t"
             << e.to << "\t" << (e.to_fwd ? '+' : '-') << "\t" << e.overlap;
        for (auto& t : e.raw_tags) cout << "\t" << t;
        cout << "\n";
    }
    // A link emitted for one allele can reference a piece that a later allele removes (calls at
    // one locus share nodes).  Drop any whose endpoints did not survive.
    int64_t dropped = 0;
    for (auto& e : new_edges) {
        const Node* f = NI(e.from); const Node* t = NI(e.to);
        if (!f || !t || f->deleted || t->deleted) { e.deleted = true; ++dropped; }
    }
    if (dropped) cerr << "[rgfa-collapse] dropped " << dropped
                      << " synthesised link(s) whose endpoints were removed by another allele\n";
    for (auto& e : new_edges)
        if (!e.deleted)
        cout << "L\t" << e.from << "\t" << (e.from_fwd ? '+' : '-') << "\t"
             << e.to << "\t" << (e.to_fwd ? '+' : '-') << "\t0M\tSR:i:" << e.sr_rank
             << "\tL1:i:" << NI(e.from)->len << "\tL2:i:" << NI(e.to)->len << "\n";
    return 0;
}
