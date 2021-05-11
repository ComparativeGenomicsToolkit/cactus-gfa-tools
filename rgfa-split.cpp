#include <map>
#include <list>
#include <cassert>
#include "rgfa-split.hpp"
#include "gfakluge.hpp"
#include "pafcoverage.hpp"

//#define debug 1

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


pair<unordered_map<int64_t, int64_t>, vector<string>> load_contig_map(const string& contgs_path) {
    unordered_map<int64_t, int64_t> contig_map;
    vector<string> contigs;
    assert(false);
    return make_pair(contig_map, contigs);
}

unordered_map<string, int64_t> load_query_mask_stats(const string& bed_path) {
    ifstream bed_file(bed_path);
    string bed_line;
    unordered_map<string, int64_t> mask_stats;
    while (getline(bed_file, bed_line)) {
        vector<string> toks;
        split_delims(bed_line, "\t\n", toks);
        if (toks.size() > 2) {
            int64_t range_size = stol(toks[2]) - stol(toks[1]);
            mask_stats[toks[0]] += range_size;
        }
    }
    return mask_stats;
}

void set_other_contig(unordered_map<int64_t, int64_t>& contig_map,
                      vector<string>& contigs,
                      function<bool(const string&)> visit_contig,
                      const string& other_name) {

    // add the other contig
    int64_t other_idx = contigs.size();
    contigs.push_back(other_name);

    // change the mapping of all unselected contigs to point to it
    for (auto& query_ref : contig_map) {
        const string& ref_contig = contigs[query_ref.second];
        if (!visit_contig(ref_contig)) {
            query_ref.second = other_idx;
        }
    }    
}

void set_other_contig_paf(unordered_map<string, int64_t>& target_to_id,
                          vector<string>& contigs,
                          function<bool(const string&)> visit_contig,
                          const string& other_name) {

    // add the other contig
    int64_t other_idx = contigs.size();
    contigs.push_back(other_name);

    // change the mapping of all unselected contigs to point to it
    for (auto& target_id : target_to_id) {
        const string& ref_contig = target_id.first;
        if (!visit_contig(ref_contig)) {
            target_id.second = other_idx;
        }
    }    
}

/**
 * Use contigs identified above to split PAF
 */
void paf_split(const string& input_paf_path,
               function<int64_t(const string&)> name_to_refid,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_prefix,
               double min_query_coverage,
               double min_small_query_coverage,
               int64_t small_coverage_threshold,
               double min_query_uniqueness,
               int64_t min_query_chunk,
               bool allow_softclip,
               int64_t ambiguous_id,
               const string& reference_prefix,
               const unordered_map<string, int64_t>& mask_stats,
               int64_t max_gap_as_match,
               int64_t min_mapq,
               ostream& log_stream) {

    // first pass, figure out which contig aligns where
    ifstream input_paf_stream(input_paf_path);

    // map query_contig to [reference_contig -> coverage]
    unordered_map<string, map<int64_t, vector<CoverageInterval>>> coverage_map;
    
    // keep track of query lengths
    unordered_map<string, int64_t> query_lengths;

    string paf_line;
    while (getline(input_paf_stream, paf_line)) {
        vector<string> toks;
        split_delims(paf_line, "\t\n", toks);

        // parse the paf columns
        string& query_name = toks[0];
        int64_t query_length = stol(toks[1]);
        string& target_name = toks[5];
        int64_t matching_bases = stol(toks[9]);
        int64_t mapq = stol(toks[11]);

        // use the map to go from the target name (rgfa node id in this case) to t
        // the reference contig (ex chr20)
        int64_t reference_id = name_to_refid(target_name);
        
        // also count tiny indels between matches
        int64_t small_gap_bases = count_small_gap_bases(toks, max_gap_as_match);

        // zero out the coverage if mapq too small
        int64_t effective_coverage = mapq >= min_mapq ? matching_bases + small_gap_bases : 0;
            
        // add the coverage of this reference contig to this query contig
        // note: important to use matching_bases here instead of just the query interval
        //       to account for softclips which can have a big impact
        coverage_map[query_name][reference_id].emplace_back(stol(toks[2]), stol(toks[3]) - 1, make_pair(effective_coverage, 0));

        // store the query length (todo: we could save a few bytes by
        // sticking it in the coverage map somewhere)
        query_lengths[query_name] = query_length;
    }

    // merge the coverage intervals
    for (auto& query_coverage : coverage_map) {
        for (auto& ref_coverage : query_coverage.second) {
            auto intervals = ref_coverage.second;
            CoverageIntervalTree coverage_intervals(intervals);
            vector<CoverageInterval> merged_intervals;
            scan_coverage_intervals(coverage_intervals, max_gap_as_match, [&](int64_t s, int64_t e, pair<int64_t, int64_t> v) {
                    merged_intervals.emplace_back(s, e, v);
                });
            ref_coverage.second = std::move(merged_intervals);
        }
    }
    
    // use the coverage map to decide a unique mapping for each query
    // (if min_query_chunk is set, this is generalized to chunks of the query)
    unordered_map<string, CoverageIntervalTree> query_ref_map;

    if (min_query_chunk <= 0) {
        // this is the "old" way where the entire query contig is mapped to a single reference
        for (auto& query_coverage : coverage_map) {
            int64_t max_coverage = 0;
            int64_t max_id;
            int64_t next_coverage = 0;
            int64_t next_id;
        
            // find the highest coverage
            for (auto& ref_coverage : query_coverage.second) {
                int64_t total_ref_coverage = 0;
                for (auto& interval : ref_coverage.second) {
                    total_ref_coverage += interval.value.first;
                }
                if (total_ref_coverage > max_coverage) {
                    next_coverage = max_coverage;
                    next_id = max_id;
                    max_id = ref_coverage.first;
                    max_coverage = total_ref_coverage;
                } else if (total_ref_coverage > next_coverage) {
                    next_id = ref_coverage.first;
                    next_coverage = total_ref_coverage;
                }
            }
            // check if it's good enough
            int64_t query_length = query_lengths[query_coverage.first];
            if (mask_stats.count(query_coverage.first)) {
                // factor in the masking stats by subtracting from denominator
                int64_t masked_bases = mask_stats.at(query_coverage.first);
                assert(masked_bases <= query_length);
                // avoid degenerate cases by making sure at least half of query contig is unmasked
                // before applying correction
                if (masked_bases < query_length / 2) {
                    query_length -= masked_bases;
                }
            }
            double query_coverage_fraction = (double)max_coverage / query_length;
            double min_coverage = min_query_coverage;
            if (small_coverage_threshold > 0 && query_length < small_coverage_threshold) {
                min_coverage = min_small_query_coverage;
            }
            bool is_ref = !reference_prefix.empty() &&
                query_coverage.first.substr(0, reference_prefix.length()) == reference_prefix;
            if (!is_ref && (query_coverage_fraction < min_coverage || 
                            (next_coverage > 0 && max_coverage < (double)next_coverage * min_query_uniqueness))) {
                log_stream << "Query contig is ambiguous: ";
                max_id = ambiguous_id;
                assert(max_id >= 0 && max_id < contigs.size());
            } else {
                log_stream << "Assigned contig to " << contigs[max_id] << ": ";
            }
            log_stream << query_coverage.first 
                       << "  len=" << query_length << " cov=" << query_coverage_fraction << " (vs " << min_coverage << ") ";
            if (next_coverage > 0) {
                log_stream << "uf=" << ((double)max_coverage / next_coverage) << " (vs " << min_query_uniqueness << ")";
                log_stream << "\n Reference contig mappings:" << "\n";
                for (auto& ref_coverage : query_coverage.second) {
                    int64_t total_ref_coverage = 0;
                    for (auto& interval : ref_coverage.second) {
                        total_ref_coverage += interval.value.first;
                    }            
                    log_stream << "  " << contigs[ref_coverage.first] << ": " << total_ref_coverage << endl;
                }
            } else {
                log_stream << "uf= infinity (vs " << min_query_uniqueness << ")" << endl;
            }
            vector<CoverageInterval> intervals;
            intervals.emplace_back(0, query_lengths[query_coverage.first] - 1, make_pair(max_coverage, max_id));
            CoverageIntervalTree interval_tree(intervals);
            query_ref_map[query_coverage.first] = interval_tree;
        }
    } else {
        // just map query ranges exceeding chunk size to which ever reference contig they are aligned to
        for (auto& query_coverage : coverage_map) {
            vector<CoverageInterval> intervals;
            for (auto& ref_coverage : query_coverage.second) {
                for (auto& interval : ref_coverage.second) {
                    //if (interval.stop - interval.start + 1 >= min_query_chunk || query_lengths[query_coverage.first] < min_query_chunk) {
                        intervals.emplace_back(interval.start, interval.stop, make_pair(interval.value.first, ref_coverage.first));
                        //}
                }
            }
            CoverageIntervalTree interval_tree(intervals);
            // we can have overlapping query intervals that map to separate chromosomes.  rather than go through
            // the rabbit hole of chopping them up to enforce constency, we just drop all but the biggest.
            vector<CoverageInterval> non_overlapping_intervals;
            interval_tree.visit_all([&](const CoverageInterval& interval) {
                    vector<CoverageInterval> overlaps = interval_tree.findOverlapping(interval.start, interval.stop);
                    bool keep = true;
                    for (CoverageInterval& overlap : overlaps) {
                        if (overlap.stop - overlap.start > interval.stop - interval.start) {
                            log_stream << "Dropping PAF line as it overlaps larger query range that maps to different contig: "
                                       << query_coverage.first << "\t" << interval.start << "\t" << (interval.stop + 1) << "\t"
                                       << contigs[interval.value.second] << endl;                        
                            keep = false;
                            break;
                        }
                    }
                    if (keep) {
                        non_overlapping_intervals.push_back(interval);
                    } 
                });
            int64_t query_length = query_lengths[query_coverage.first];
            double min_coverage = min_query_coverage;
            int64_t masked_bases = 0;
            if (mask_stats.count(query_coverage.first)) {
                // factor in the masking stats by subtracting from denominator
                masked_bases = mask_stats.at(query_coverage.first);
                assert(masked_bases <= query_length);
                // avoid degenerate cases by making sure at least half of query contig is unmasked
                // before applying correction
                if (masked_bases >= query_length / 2) {
                    masked_bases = 0;
                }
            }

            if (small_coverage_threshold > 0 && query_lengths[query_coverage.first] < small_coverage_threshold) {
                min_coverage = min_small_query_coverage;
            }
            smooth_query_intervals(query_coverage.first, query_length, masked_bases,
                                   non_overlapping_intervals, min_coverage, min_query_uniqueness, min_query_chunk,
                                   contigs, allow_softclip, log_stream);
            
            interval_tree = CoverageIntervalTree(non_overlapping_intervals);
            
            query_ref_map[query_coverage.first] = interval_tree;
        }
        // add in complement intervals as ambiguous
        for (auto& query_ref : query_ref_map) {
            int64_t query_len = query_lengths[query_ref.first];
            // todo: make better
            vector<bool> covered(query_len, false);
            vector<CoverageInterval> intervals;
            query_ref.second.visit_all([&](const CoverageInterval& interval) {
                    for (int64_t p = interval.start; p <= interval.stop; ++p) {
                        covered[p] = true;
                    }
                    intervals.push_back(interval);
                });
            size_t intervals_size = intervals.size();
            int64_t start = -1;
            for (int64_t i = 0; i < covered.size(); ++i) {
                if (!covered[i] && start == -1) {
                    start = i;
                } else if ((covered[i] || i == covered.size() -1) && start >=0) {
                    int64_t stop = covered[i] ? i - 1: i;
                    intervals.emplace_back(start, stop, make_pair(0, ambiguous_id));
                    start = -1;
                }
            }
            if (intervals.size() > intervals_size) {
                query_ref.second = CoverageIntervalTree(intervals);
            }
#ifdef debug
            cerr << "intervals for " << query_ref.first << endl;
            query_ref.second.visit_all([&](const CoverageInterval& interval) {
                    cerr << " " << interval << endl;
                });
            cerr << "." << endl;
#endif
        }        
    }
    
    coverage_map.clear();
    query_lengths.clear();
        
    // second pass, do the splitting
    input_paf_stream.clear();
    input_paf_stream.seekg(0, ios::beg);
    
    // note: we're assuming a small number of reference contigs (ie 23), so we can afford to open file
    // for each. 
    unordered_map<int64_t, ofstream*> out_files;

    // load up the query contigs for downstream fasta splitting
    unordered_map<int64_t, unordered_set<string> > query_map;

    // keep track of every unique taret
    unordered_set<string> target_set;

    // keep track of pafs written (pad out with empty files to help downstream scripts)
    vector<bool> pafs_written(contigs.size(), false);

    while (getline(input_paf_stream, paf_line)) {
        vector<string> toks;
        split_delims(paf_line, "\t\n", toks);

        // parse the paf columns
        string& query_name = toks[0];
        int64_t query_start = stol(toks[2]);
        int64_t query_end = stol(toks[3]);
        string& target_name = toks[5];
        target_set.insert(target_name);

        // use the map to go from the target name (rgfa node id in this case) to
        // the reference contig (ex chr20)
        int64_t target_reference_id = name_to_refid(target_name);

        assert(query_ref_map.count(query_name));
        CoverageIntervalTree& intervals = query_ref_map.at(query_name);
        vector<CoverageInterval> overlaps = intervals.findOverlapping(query_start, query_end - 1);

        if (overlaps.size() > 1) {
            // the only way for this to happen is if the paf line corresponds to a query overlap that gets
            // filtered out.  all we can do is erase it, as even putting it in the ambiguous pile can cause
            // trouble later on
            continue;
        }
        assert(overlaps.size() == 1);
        int64_t reference_id = overlaps[0].value.second;
        const string& reference_contig = contigs[reference_id];

        // do both the query and reference sequences fall in the same chromosome, and we wnat to visit that
        // chromosome?  if so, we write the paf line, otherwise it's effectively filtered out
        if ((reference_id == target_reference_id && visit_contig(reference_contig)) ||
            (ambiguous_id >= 0 && reference_contig == contigs[ambiguous_id])) {
            ofstream*& out_paf_stream = out_files[reference_id];
            if (out_paf_stream == nullptr) {
                string out_paf_path = output_prefix + reference_contig + ".paf";
                out_paf_stream = new ofstream(out_paf_path);
                pafs_written[reference_id] = true;
                assert(out_files.size() < 100);
                if (!(*out_paf_stream)) {
                    cerr << "error: unable to open output paf file: " << out_paf_path << endl;
                    exit(1);
                }
            }
            apply_paf_query_offsets(toks, overlaps[0].start, overlaps[0].stop);
            for (size_t i = 0; i < toks.size(); ++i) {
                if (i > 0) {
                    *out_paf_stream << "\t";
                }
                *out_paf_stream << toks[i];
            }
            *out_paf_stream << "\n";
            // remember this query contig for future fasta splitting
            query_map[reference_id].insert(query_name);
        } 
        
    }

    // write some empty pafs
    for (size_t i = 0; i < pafs_written.size(); ++i) {
        if (!pafs_written[i]) {
            string out_paf_path = output_prefix + contigs[i] + ".paf";
            ofstream empty_paf(out_paf_path);
        }
    }

    // clean up the files
    for (auto& ref_stream : out_files) {
        delete ref_stream.second;
    }
    out_files.clear();

    // write the query_contigs
    for (auto& ref_queries : query_map) {
        const string& reference_contig = contigs[ref_queries.first];
        string out_contigs_path = output_prefix + reference_contig + ".fa_contigs";
        ofstream out_contigs_stream(out_contigs_path);
        if (!out_contigs_stream) {
            cerr << "error: unable to open output contigs path: " << out_contigs_path << endl;
            exit(1);
        }
        for (const string& query_name : ref_queries.second) {
            out_contigs_stream << query_name << "\n";
        }
        out_contigs_stream.close();
    }

    // write the target contigs
    // start by sorting by reference contig
    vector<string> mg_contigs;
    mg_contigs.reserve(target_set.size());
    for (const auto& target_name : target_set) {
        mg_contigs.push_back(target_name);
    }
    std::sort(mg_contigs.begin(), mg_contigs.end(), [&](const string& a, const string& b) {
            return contigs[name_to_refid(a)] < contigs[name_to_refid(b)];
        });
    int64_t prev_ref_contig_id = -1;
    ofstream out_contigs_stream;
    for (const auto& target_name : mg_contigs) {
        int64_t reference_contig_id = name_to_refid(target_name);
        const string& reference_contig = contigs[name_to_refid(target_name)];
        if (visit_contig(reference_contig) ||
            (ambiguous_id >= 0 && reference_contig == contigs[ambiguous_id])) {  
            if (reference_contig_id != prev_ref_contig_id) {
                string out_contigs_path = output_prefix + reference_contig + ".fa_contigs";
                if (out_contigs_stream.is_open()) {
                    out_contigs_stream.close();
                }
                out_contigs_stream.open(out_contigs_path, std::ios_base::app);
                if (!out_contigs_stream) {
                    cerr << "error: unable to open output contigs path: " << out_contigs_path << endl;
                    exit(1);
                }
                prev_ref_contig_id = reference_contig_id;
            }
            out_contigs_stream << target_name << endl;
        }
    }
}

void gfa_split(const string& rgfa_path,
               const unordered_map<int64_t, int64_t>& contig_map,
               const vector<string>& contigs,
               function<bool(const string&)> visit_contig,
               const string& output_prefix) {

    ifstream input_gfa_stream(rgfa_path);
    assert(input_gfa_stream);

    // note: we're assuming a small number of reference contigs (ie 23), so we can afford to open file
    // for each. 
    unordered_map<int64_t, ofstream*> out_files;

    string gfa_line;
    while (getline(input_gfa_stream, gfa_line)) {
        vector<string> toks;
        split_delims(gfa_line, "\t\n", toks);

        const string* ref_contig = nullptr;
        int64_t reference_id;
        if (toks[0] == "S") {
            int64_t seq_id = node_id(toks[1]);
            assert(contig_map.count(seq_id));
            reference_id = contig_map.at(seq_id);
            ref_contig = &contigs[reference_id];
        } else if (toks[0] == "L") {
            int64_t seq_id = node_id(toks[1]);
            assert(contig_map.count(seq_id));            
            reference_id = contig_map.at(seq_id);
            int64_t sink_seq_id = node_id(toks[3]);
            assert(contig_map.count(sink_seq_id));            
            int64_t sink_reference_id = contig_map.at(sink_seq_id);
            assert(sink_reference_id == reference_id);
            ref_contig = &contigs[reference_id];
        }
        if (ref_contig != nullptr && visit_contig(*ref_contig)) {
            ofstream*& out_gfa_stream = out_files[reference_id];
            if (out_gfa_stream == nullptr) {
                string out_gfa_path = output_prefix + *ref_contig + ".gfa";
                out_gfa_stream = new ofstream(out_gfa_path);
                assert(out_files.size() < 100);
                if (!(*out_gfa_stream)) {
                    cerr << "error: unable to open output gfa file: " << out_gfa_path << endl;
                    exit(1);
                }
            }
            *out_gfa_stream << gfa_line << "\n";
        }
    }

    // clean up the files
    for (auto& ref_stream : out_files) {
        delete ref_stream.second;
    }
    out_files.clear();
}

int64_t count_small_gap_bases(const vector<string>& toks, int64_t max_gap_as_match) {

    bool after_match = false;
    int64_t running_ins = 0;
    int64_t running_del = 0;
    int64_t total_gap = 0;
    for (int i = 12; i < toks.size(); ++i) {
        if (toks[i].substr(0, 5) == "cg:Z:") {
            for_each_cg(toks[i], [&](const string& val, const string& cat) {
                    int64_t len = stol(val);
                    if (cat == "M") {
                        if (after_match && running_ins < max_gap_as_match && running_del < max_gap_as_match) {
                            total_gap += running_ins;
                        }
                        running_ins = 0;
                        running_del = 0;
                        after_match = true;
                    } else if (cat == "I") {
                        running_ins += len;
                    } else {
                        assert(cat == "D");
                        running_del += len;
                    }
                });
        }
    }

    return total_gap;
}

void scan_coverage_intervals(CoverageIntervalTree& interval_tree, int64_t padding, function<void(int64_t, int64_t, pair<int64_t, int64_t>)> fn) {
    unordered_set<const CoverageInterval*> visited;
    // go through every interval, and all its overlaps once
    interval_tree.visit_all([&](const CoverageInterval& interval) {
            if (!visited.count(&interval)) {
                // collect a set of all overlapping intervals (taking into account padding)
                // and mark them all as visited
                visited.insert(&interval);
                vector<const CoverageInterval*> overlaps = {&interval};
                int64_t idx_to_search = 0;
                // loop here to collect all transitive overlaps
                while (idx_to_search < overlaps.size()) {
                    const CoverageInterval* to_search = overlaps[idx_to_search++];
                    interval_tree.visit_overlapping(to_search->start - padding, to_search->stop + padding, [&](const CoverageInterval& overlapping_interval) {
                            if (!visited.count(&overlapping_interval)) {
                                overlaps.push_back(&overlapping_interval);
                                visited.insert(&overlapping_interval);
                            }
                        });
                }
                // merge up the set, using weighted average to determine a rough approximation of the merged coverage
                size_t total_coverage_numerator = 0;
                size_t total_coverage_denominator = 0;
                int64_t start = numeric_limits<int64_t>::max();
                int64_t end = -1;
                for (auto& overlap_interval : overlaps) {
                    start = std::min(start, overlap_interval->start);
                    end = std::max(end, overlap_interval->stop);
                    total_coverage_numerator += overlap_interval->value.first;
                    total_coverage_denominator += overlap_interval->stop - overlap_interval->start + 1;
                }
                double coverage_density = (double)total_coverage_numerator / (double)total_coverage_denominator;
                fn(start, end, make_pair((end - start + 1) * coverage_density, 0));
            }
        });
}


void smooth_query_intervals(const string& query_name, int64_t query_length, int64_t masked_bases,
                            vector<CoverageInterval>& intervals, double min_coverage, double min_uniqueness,
                            int64_t min_chunk, const vector<string>& ref_contigs, bool allow_softclip,
                            ostream& log_stream) {

    if (intervals.empty()) {
        return;
    }
    // total up coverage by reference contig
    map<int64_t, int64_t> coverage_by_contig;
    for (const CoverageInterval& interval : intervals) {
        coverage_by_contig[interval.value.second] += interval.value.first;
    }

    // find the top two
    pair<int64_t, int64_t> top(-1, -1);
    pair<int64_t, int64_t> next(-1, -1);
    for (const auto& cov : coverage_by_contig) {
        if (cov.second > top.second) {
            next = top;
            top = cov;
        } else if (cov.second > next.second) {
            next = cov;
        }
    }

    // look for regions to possibly clip out, by finding runs of consecutive
    // mappings to non-top contigs
    vector<vector<size_t>> clip_candidates;
    if (min_chunk > 0) {
        int64_t ref = -1;
        for (size_t i = 0; i < intervals.size(); ++i) {
            const CoverageInterval& interval = intervals[i];
            if (interval.value.second != top.first) {
                if (interval.value.first != ref || clip_candidates.empty()) {
                    clip_candidates.push_back({});
                }
                clip_candidates.back().push_back(i);
            }
            ref = interval.value.second;
        }
    }
    vector<CoverageInterval> clip_intervals;
    unordered_set<size_t> clip_set;
    int64_t total_clip_length = 0;
    for (size_t i = 0; i < clip_candidates.size(); ++i) {
        int64_t min_pos = query_length;
        int64_t max_pos = -1;
        int64_t max_interval_length = 0;
        int64_t total_coverage = 0;
        for (size_t j = 0; j < clip_candidates[i].size(); ++j) {
            const CoverageInterval& interval = intervals[clip_candidates[i][j]];
            max_interval_length = std::max(max_interval_length, interval.stop - interval.start +1);
            min_pos = std::min(min_pos, interval.start);
            max_pos = std::max(max_pos, interval.stop);
            total_coverage += interval.value.first;
        }
        // a candidate clip at beginning or end will have to own that sequence
        if (clip_candidates[i][0] == 0) {
            min_pos = 0;
        }
        if (clip_candidates[i].back() == intervals.size() - 1) {
            max_pos = query_length - 1;
        }
        if (max_interval_length > min_chunk &&
            (double)total_coverage / (max_pos - min_pos + 1) >= min_coverage) {
            int64_t ref_contig = intervals[clip_candidates[i][0]].value.second;
            assert(ref_contig != top.first);
            clip_intervals.push_back(CoverageInterval(min_pos, max_pos, make_pair(total_coverage, ref_contig)));
            clip_set.insert(clip_candidates[i][0]);
            total_clip_length += max_pos - min_pos + 1;
        }
    }

    // smooth out the intervals.
    vector<CoverageInterval> smooth_intervals;    
    size_t next_clip_idx = 0;
    size_t prev_top = intervals.size();
    for (size_t i = 0; i < intervals.size(); ++i) {
        if (clip_set.count(i)) {
            smooth_intervals.push_back(clip_intervals.at(next_clip_idx++));            
        } else if (intervals[i].value.second == top.first) {
            if (!smooth_intervals.empty() && smooth_intervals.back().value.second == top.first &&
                (i - 1 == prev_top && intervals[i].start - intervals[prev_top].stop < min_chunk)) {
                smooth_intervals.back().stop = intervals[i].stop;
                smooth_intervals.back().value.first += intervals[i].value.first;
            } else {
                smooth_intervals.push_back(intervals[i]);
            }
            prev_top = i;
        }
    }

    // filter out small bits
    vector<CoverageInterval> filtered_intervals;
    int64_t min_len_filter = std::min(min_chunk, int64_t(query_length * min_coverage));
    for (size_t i = 0; i < smooth_intervals.size(); ++i) {
        if (smooth_intervals[i].value.second != top.first || 
            smooth_intervals[i].stop - smooth_intervals[i].start > min_len_filter) {
            filtered_intervals.push_back(smooth_intervals[i]);
        } else {
            log_stream << "Unable to smooth small fragment: " << query_name << " " << smooth_intervals[i].start << "-" << smooth_intervals[i].stop << " -> "
                       << ref_contigs[smooth_intervals[i].value.second] << endl;
            top.second -= smooth_intervals[i].value.first;
        }
    }
    smooth_intervals = filtered_intervals;

    // clamp intervals to ends of query.  if allow_softclip is set, then we only do
    // so if they are within "min_chunk" (and count as softclip otherwise)
    int64_t softclip = 0;
    if (!smooth_intervals.empty()) {
        if (allow_softclip) {
            if (smooth_intervals[0].start <= min_chunk) {
                smooth_intervals[0].start = 0;
            } else {
                softclip += smooth_intervals[0].start;
            }
            if (smooth_intervals.back().stop > query_length - min_chunk) {
                smooth_intervals.back().stop = query_length - 1;
            } else {
                softclip += query_length - smooth_intervals.back().stop -1;
            }
        } else {
            smooth_intervals[0].start = 0;
            smooth_intervals.back().stop = query_length - 1;
        }
    }
    
    for (size_t i = 0; i < smooth_intervals.size(); ++i) {
        if (i > 0 && smooth_intervals[i].value.second == top.first && smooth_intervals[i].start != smooth_intervals[i-1].stop + 1) {
            smooth_intervals[i].start = smooth_intervals[i-1].stop + 1;
        }
        if (i < smooth_intervals.size() - 1 && smooth_intervals[i].value.second == top.first && smooth_intervals[i].stop != smooth_intervals[i+1].start - 1) {
            smooth_intervals[i].stop = smooth_intervals[i+1].start - 1;
        }
    }

    // we are now ready to make our decision based on coverage (note that the clip candidates have already passed,
    // but using them is dependent on the rest of the contig)

    // TODO : proper logging
#ifdef debug
    log_stream << "adjusted coverage = " << top.second << " / max(" << query_length << " - " << masked_bases << " - " << total_clip_length << " - " << softclip << ", " << top.second << ")" << endl;
#endif
    double adjusted_coverage = 0;
    if (top.second > 0) {
        adjusted_coverage = (double)top.second / std::max(query_length - std::max(masked_bases, softclip) - total_clip_length, top.second);
    }
    if (adjusted_coverage > min_coverage) {
        log_stream << "Assigning contig " << query_name << " with adjusted covarege " << adjusted_coverage << " vs " << min_coverage << " " << query_name << " to:\n";
        for (const auto& interval : smooth_intervals) {
            log_stream << interval.start << "-" << interval.stop << " -> " << ref_contigs[interval.value.second] << "(" << interval.value.first << ")" << endl;
        }
        swap(intervals, smooth_intervals);
    } else {
        log_stream << "Leaving " << query_name << " as ambigious with adjusted covarege " << adjusted_coverage << " vs " << min_coverage << " " << endl;
#ifdef debug
        log_stream << "input" << endl;
        for (const auto& interval : intervals) {
            log_stream << interval.start << "-" << interval.stop << " -> " << ref_contigs[interval.value.second] << "(" << interval.value.first << ")" << endl;
        }        
        log_stream << "smooth" << endl;
        for (const auto& interval : smooth_intervals) {
            log_stream << interval.start << "-" << interval.stop << " -> " << ref_contigs[interval.value.second] << "(" << interval.value.first << ")" << endl;
        }
#endif
        intervals.clear();
    }
}

void apply_paf_query_offsets(vector<string>& paf_toks, int64_t query_fragment_start, int64_t query_fragment_end) {

    int64_t query_length = stol(paf_toks[1]);

    if (query_fragment_end - query_fragment_start + 1 == query_length) {
        assert(query_fragment_start == 0);
        // nothing to do, as the fragment spans the entire query sequence
        return;
    }

    int64_t query_start = stol(paf_toks[2]);
    // note, paf coordinates are 0-based end exclusive, but internally we're always using
    // 0-based inclusive.  
    int64_t query_end = stol(paf_toks[3]);

    tuple<string, int64_t, int64_t> parsed_query_name = parse_faidx_subpath(paf_toks[0]);
    string& query_name = get<0>(parsed_query_name);

    // apply adjustments to convert back to coordinates on the original contig
    if (get<1>(parsed_query_name) > 0) {
        query_start += get<1>(parsed_query_name);
        query_end += get<1>(parsed_query_name);
    }

    assert(query_fragment_start <= query_start && query_fragment_end >= query_end - 1);

    int64_t delta = query_fragment_start;
    int64_t adjusted_query_start = query_start - delta;
    int64_t adjusted_query_end = query_end - delta;


    paf_toks[0] = make_faidx_subpath(query_name, query_fragment_start, query_fragment_end);
    paf_toks[1] = to_string(query_fragment_end - query_fragment_start + 1);
    paf_toks[2] = to_string(adjusted_query_start);
    paf_toks[3] = to_string(adjusted_query_end);
}

// todo: harmonize with vg::Paths::parse_subpath_name
tuple<string, int64_t, int64_t>  parse_faidx_subpath(const string& name) {
    size_t tag_offset = name.rfind(":");
    if (tag_offset == string::npos) {
        return make_tuple(name, 0, -1);
    } else {
        string offset_str = name.substr(tag_offset + 1, name.length() - tag_offset - 2);
        size_t sep_offset = offset_str.find("-");
        assert(sep_offset != string::npos && sep_offset > 0);
        int64_t start_val = stol(offset_str.substr(0, sep_offset)) - 1;
        int64_t end_val = stol(offset_str.substr(sep_offset + 1)) - 1;
        return make_tuple(name.substr(0, tag_offset), start_val, end_val);
    }
}

// todo: harmonize with vg::Paths::make_subpath_name
string make_faidx_subpath(const string& name, int64_t start, int64_t end) {
    stringstream ss;
    ss << name << ":" << (start + 1) << "-" << (end + 1);
    return ss.str();    
}
