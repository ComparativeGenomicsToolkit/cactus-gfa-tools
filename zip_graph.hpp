/*
  zip_graph.hpp -- the rGFA model of rgfa-zip: nodes, links, CSR handle adjacency, input checks,
  emit, and the global placement assert.

  Model
    - Graph::nodes holds the input S lines sorted by numeric id (names are s<int>), so NodeId order
      is id order and nothing depends on the order of the input lines.  Sequences are kept.
    - Graph::links holds the input L lines in canonical order, each with its from/to sides as
      written, its overlap, its SR (or -1) and its raw tags (L1/L2 are recomputed at emit).
    - The CSR (out()) gives, for every handle, the handles one step away and the link used.  It is
      built from the input links and is never updated: every decision is made on the input graph.
    - Edits are applied once at the end by the edit module through add_node/add_link and the
      deleted flags; write_rgfa then emits the non-deleted nodes and links.

  Input checks (read_rgfa; any failure is ZipError(EXIT_INPUT) and nothing is written):
    every S line has SN, SO and SR; each SN has one rank; its rank-0 nodes tile it by SO; names are
    unique s<int>; every L end exists; links use '+'/'-' and overlap 0M or *; every node is
    reachable from rank 0 (the graph is placeable).  L lines without SR only warn.
*/
#pragma once

#include "zip_common.hpp"

#include <unordered_map>

namespace zip {

struct Node {
    int64_t id = 0;            // numeric id; the GFA name is "s<id>"
    std::string seq;
    std::string tags;          // raw tags after the sequence, tab-joined, emitted verbatim
    int32_t sn = -1;           // index into Graph::sn_names
    int32_t sr = -1;           // SR: rank (0 = reference)
    int64_t so = -1;           // SO: offset of the node's first base on its SN
    bool deleted = false;
    int64_t len() const { return (int64_t)seq.size(); }
};

struct Link {
    Side a = 0, b = 0;         // from-side and to-side as written ('L u uo v vo': a = u.R if uo is '+', else u.L;
                               // b = v.L if vo is '+', else v.R)
    int32_t sr = -1;           // SR tag; -1 when absent
    std::string overlap = "0M";
    std::string tags;          // raw tags after the overlap, tab-joined; L1/L2 are rewritten at emit
    bool deleted = false;
};

// A link equals its reverse, so its canonical key is the unordered pair of sides it joins.
struct LinkKey {
    Side lo = 0, hi = 0;
    bool operator<(const LinkKey& o) const { return lo != o.lo ? lo < o.lo : hi < o.hi; }
    bool operator==(const LinkKey& o) const { return lo == o.lo && hi == o.hi; }
    bool operator!=(const LinkKey& o) const { return !(*this == o); }
};
inline LinkKey link_key(Side a, Side b) { return a < b ? LinkKey{a, b} : LinkKey{b, a}; }
struct LinkKeyHash {
    size_t operator()(const LinkKey& k) const { return std::hash<uint64_t>()(((uint64_t)k.lo << 32) | k.hi); }
};

// One oriented step: from the handle whose list this is, over link `link`, to handle `to`.
struct Edge {
    Handle to;
    uint32_t link;
};

struct EdgeRange {
    const Edge* b;
    const Edge* e;
    const Edge* begin() const { return b; }
    const Edge* end() const { return e; }
    size_t size() const { return (size_t)(e - b); }
    bool empty() const { return b == e; }
};

class Graph {
public:
    std::vector<std::string> header_lines;      // H lines, verbatim, input order
    std::vector<Node> nodes;                    // input nodes sorted by id; pieces appended by the edit
    std::vector<Link> links;                    // input links in canonical order; new links appended
    std::vector<std::string> sn_names;          // SN table, sorted by name (sn index order = name order)
    std::vector<int32_t> sn_rank;               // rank of each SN
    std::vector<std::vector<NodeId>> ref_nodes; // per SN: its rank-0 nodes sorted by SO (empty for rank > 0)
    int64_t max_id = -1;                        // largest input id
    size_t n_input_nodes = 0, n_input_links = 0;
    std::vector<uint64_t> ref_hash;             // per SN: hash of the input rank-0 sequence in SO order (final assert)

    // ---- lookups
    NodeId find_id(int64_t id) const {
        auto it = id_index_.find(id);
        return it == id_index_.end() ? NONE : it->second;
    }
    // NONE when the name is not s<int> or not in the graph
    NodeId find_name(const std::string& name) const;
    int32_t find_sn(const std::string& sn) const;     // -1 if absent
    std::string name(NodeId n) const { return "s" + std::to_string(nodes[n].id); }
    std::string handle_str(Handle h) const { return name(handle_node(h)) + (handle_rev(h) ? "-" : "+"); }
    int64_t len(NodeId n) const { return nodes[n].len(); }
    bool is_ref(NodeId n) const { return nodes[n].sr == 0; }
    int64_t start(NodeId n) const { return nodes[n].so; }                    // SO
    int64_t end(NodeId n) const { return nodes[n].so + nodes[n].len(); }     // SO + LN

    // ---- adjacency over the INPUT links (both directions of every link)
    // out(h): the handles a walk can step to after h.  Predecessors of h are flip(e.to) for e in out(flip(h)).
    // Only input handles have adjacency; a handle of a node added by the edit gets an empty range.
    EdgeRange out(Handle h) const {
        if ((size_t)h + 1 >= csr_off_.size()) return EdgeRange{nullptr, nullptr};
        return EdgeRange{csr_.data() + csr_off_[h], csr_.data() + csr_off_[h + 1]};
    }
    int32_t link_sr(uint32_t link) const { return links[link].sr; }

    // ---- reference helpers (rank-0 nodes of one SN)
    // the rank-0 node of `sn` covering position pos, or NONE
    NodeId ref_node_at(int32_t sn, int64_t pos) const;
    // reference sequence [lo, hi) of `sn`
    std::string ref_seq(int32_t sn, int64_t lo, int64_t hi) const;

    // ---- sequences
    void append_handle_seq(std::string& out, Handle h) const;
    std::string spell(const std::vector<Handle>& walk) const;

    // ---- edit support (used once, after planning; the CSR is not updated)
    NodeId add_node(Node&& n);       // appended; its id is registered for find_id
    uint32_t add_link(Link&& l);     // appended

    // ---- construction (read_rgfa calls these)
    void build_index();              // id index
    void build_csr();                // CSR from the non-deleted links

private:
    std::unordered_map<int64_t, NodeId> id_index_;
    std::vector<uint32_t> csr_off_;
    std::vector<Edge> csr_;
};

// Parse "s<int>" (no sign, no leading zeros except "s0").  False if not that form.
bool parse_node_name(const std::string& name, int64_t& id);

// Read and check an rGFA (plain or gzip).  Any check failure is ZipError(EXIT_INPUT).
void read_rgfa(const std::string& path, Graph& g);

// Emit H lines, then non-deleted S lines by id, then non-deleted L lines by canonical key (on ids).
// L1/L2 tags present on a link are rewritten from the current endpoint lengths.
void write_rgfa(AtomicWriter& w, const Graph& g);

// The global placement assert on the (edited) graph: unique ids, no dangling or duplicate link,
// one rank per SN, rank-0 nodes of each SN still tile it, every node reachable from rank 0
// (undirected, as rgfa-split does).  Failure is ZipError(EXIT_INVARIANT).
void check_placement(const Graph& g);

// Tags for a new piece node / a new link (L1/L2 placeholders are filled in by write_rgfa).
std::string make_node_tags(int64_t len, const std::string& sn, int64_t so, int32_t sr);
std::string make_link_tags(int32_t sr);

} // namespace zip
