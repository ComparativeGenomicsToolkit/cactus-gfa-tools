#include <map>
#include <list>
#include <cassert>
#include "rgfa2contig.hpp"
#include "gfakluge.hpp"

/*

sort nodes by rank (ascending)

index all edges

for each node in sorted list
 if rank == 0
   node.contig = contig from gfa optional element
 else:
   scan all edges for rank i-1 adjacencies
   build consensus contig
*/
pair<unordered_map<int64_t, int64_t>, vector<string>> rgfa2contig(const string& gfa_path) {

    // map rank -> nodes
    map<int64_t, list<int64_t>> rank_to_nodes;

    // map node -> rank
    unordered_map<int64_t, int64_t> node_to_rank;

    // store edges (don't care about orientation, just connectivity)
    unordered_multimap<int64_t, int64_t> edges;

    // output contig ids
    vector<string> contigs;
    // used to compute above
    unordered_map<string, int64_t> contig_map;

    // output map of node id to contig id
    unordered_map<int64_t, int64_t> node_to_contig;

    // get the node ranks, and rank-0 contigs
    function<void(const gfak::sequence_elem&)> visit_seq = [&](const gfak::sequence_elem& gfa_seq) {
        int64_t gfa_id = node_id(gfa_seq.name);
        bool found_SN = false;
        bool found_SR = false;
        int64_t rank;
        string contig;
        // parse the SN (contig) and SR (rank) optional tags
        for (const gfak::opt_elem& oe : gfa_seq.opt_fields) {
            if (oe.key == "SN") {
                assert(found_SN == false);
                contig = oe.val;
                found_SN = true;
            } else if (oe.key == "SR") {
                assert(found_SR == false);
                rank = stol(oe.val);
                assert(rank >= 0);                               
                found_SR = true;
            }
        }
        assert(found_SN);
        assert(found_SR);                       
        // remember the rank
        rank_to_nodes[rank].push_back(gfa_id);
        node_to_rank[gfa_id] = rank;
        // remember the contig (if rank 0)
        if (rank == 0) {
            // get the contig id
            int64_t contig_id;
            if (contig_map.count(contig)) {
                contig_id = contig_map[contig];
            } else {
                contig_id = contig_map.size();
                contig_map[contig] = contig_id;
                contigs.push_back(contig);
            }
            node_to_contig[gfa_id] = contig_id;
        }
    };

    // get the edges
    function<void(const gfak::edge_elem&)> visit_edge = [&](const gfak::edge_elem& gfa_edge) {
        int64_t source_id = node_id(gfa_edge.source_name);
        int64_t sink_id = node_id(gfa_edge.sink_name);
        edges.insert(make_pair(source_id, sink_id));
        edges.insert(make_pair(sink_id, source_id));
    };

    // load the GFA into memory
    gfak::GFAKluge kluge;
    kluge.for_each_sequence_line_in_file(gfa_path.c_str(), visit_seq);
    kluge.for_each_edge_line_in_file(gfa_path.c_str(), visit_edge);

    // fill out the contigs by rank
    for (auto& rank_nodes : rank_to_nodes) {
        // rank 0 was added above
        if (rank_nodes.first > 0) {
            const int64_t& rank = rank_nodes.first;
            auto& nodes_at_rank = rank_nodes.second;
            // todo: clean up (this was a mistake in original design where i forgot that not every rank i node
            // connects to rank i-1)
            int64_t consecutive_pushes = 0;
            while (!nodes_at_rank.empty()) {
                int64_t node_id = nodes_at_rank.back();
                nodes_at_rank.pop_back();
                
                // contig_id -> count of times it's connected
                unordered_map<int64_t, int64_t> counts;
                auto edge_iterators = edges.equal_range(node_id);
                for (auto e = edge_iterators.first; e != edge_iterators.second; ++e) {
                    int64_t& other_id = e->second;
                    int64_t other_rank = node_to_rank[other_id];
                    if (other_rank < rank ||
                        (other_rank == rank && node_to_contig.count(other_id))) {
                        int64_t other_contig = node_to_contig[other_id];
                        ++counts[other_contig];
                    }
                }
                if (counts.size() == 0) {
                    // this node isn't connected to any nodes with rank -1, try it later
                    nodes_at_rank.push_front(node_id);
                    ++consecutive_pushes;
                    if (consecutive_pushes > nodes_at_rank.size()) {
                        cerr << "[error] Unable to assign contigs for the following nodes at rank " << rank << ":\n";
                        for (const auto& ni : nodes_at_rank) {
                            cerr << ni << endl;
                        }
                        exit(1);
                    }
                } else if (counts.size() > 1) {
                    cerr << "[error] Conflict found for node \"" << node_id << "\" with rank \"" << rank << ":\n";
                    for (auto& count_elem : counts) {
                        // we could use a heuristic to resolve. but do not expect this from minigraph output
                        cerr << "\tcontig=" << contigs[count_elem.first] << " count=" << count_elem.second << endl;
                    }
                    exit(1);                    
                } else {
                    assert(counts.size() == 1);
                    // set the contig to the unambiguous neighbour
                    node_to_contig[node_id] = counts.begin()->first;
                    consecutive_pushes = 0;
                }
            }
        }
    }

    return make_pair(node_to_contig, contigs);
}
