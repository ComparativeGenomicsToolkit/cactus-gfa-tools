#include "paf2stable.hpp"
#include "pafcoverage.hpp"

#define debug

void update_stable_mapping_info(const vector<string>& paf_toks,
                                unordered_map<string, int64_t>& query_name_to_id,
                                vector<pair<string, int64_t>>& query_id_to_info,
                                unordered_map<string, vector<StableInterval>>& target_to_intervals) {

    // update the query info
    const string& query_name = paf_toks[0];
    int64_t query_length = stol(paf_toks[1]);
    int64_t query_id;
    if (query_name_to_id.count(query_name)) {
        query_id = query_name_to_id[query_name];
    } else {
        query_id = query_id_to_info.size();
        query_name_to_id[query_name] = query_id;
        query_id_to_info.push_back(make_pair(query_name, query_length));
    }

    // update the mapping interval info for the target
    const string& target_name = paf_toks[5];    
    vector<StableInterval>& target_intervals = target_to_intervals[target_name];

    bool is_reverse = paf_toks[4] == "-";
    // todo: reverse strand!!! (ugh)
    assert(!is_reverse);
   
    int64_t query_pos = stol(paf_toks[2]);
    int64_t target_pos = stol(paf_toks[7]);

    vector<pair<string, string>> cigars;
    for (int i = 12; i < paf_toks.size(); ++i) {
        if (paf_toks[i].substr(0, 5) == "cg:Z:") {
            for_each_cg(paf_toks[i], [&](const string& val, const string& cat) {
                    cigars.push_back(make_pair(val, cat));
                });
        }
    }

    for (const auto& vc : cigars) {
        const string& val = vc.first;
        const string& cat = vc.second;
        int64_t len = stol(val);
        if (cat == "M") {
#ifdef debug
            cerr << "adding interval for " << target_name << ": " << target_pos << "-" << (target_pos + len - 1)
                 << " ==> " << query_name << ": " << query_pos << ", rev=" << is_reverse << endl;
#endif             
            StableInterval interval(target_pos, target_pos + len - 1,
                                    make_tuple(query_id, query_pos, is_reverse));
            target_intervals.push_back(interval);
        }
        if (cat == "M") {
            query_pos += len;
            target_pos += len;
        } else if (cat == "I") {
            query_pos += len;
        } else if (cat == "D") {
            target_pos += len;
        } else {
            assert(false);
        }
    }
}

unordered_map<string, StableIntervalTree> create_interval_trees(unordered_map<string, vector<StableInterval>>& target_to_intervals) {
    unordered_map<string, StableIntervalTree> target_to_interval_trees;
    for (auto& kv : target_to_intervals) {
        vector<StableInterval>& intervals = kv.second;
        StableIntervalTree interval_tree(intervals);
        // intervals that are not contained in another interval
        vector<StableInterval> intervals_no_nested;
        for (StableInterval& interval : intervals) {
            vector<StableInterval> overlapping = interval_tree.findOverlapping(interval.start, interval.stop);
            bool is_nested = false;
            bool found_self = false;
            for (size_t i = 0; i < overlapping.size() && !is_nested; ++i) {
                if (!found_self && interval.start == overlapping[i].start && interval.stop == overlapping[i].stop &&
                    interval.value == overlapping[i].value) {
                    found_self = true;
                } else {
                    is_nested = interval.start >= overlapping[i].start && interval.stop <= overlapping[i].stop;
                }
            }
            if (!is_nested) {
                intervals_no_nested.push_back(interval);
            }
        }
        target_to_interval_trees[kv.first] = StableIntervalTree(intervals_no_nested);
        intervals.clear();

#ifdef debug
        cerr << "Interval Tree (" << kv.first << "):";
        for (auto interval : intervals_no_nested) {
            cerr << endl << "   " << interval;
        }
        cerr << endl;
#endif

    }
    target_to_intervals.clear();
    return target_to_interval_trees;
}


void paf_to_stable(const vector<string>& paf_toks,
                   const vector<pair<string, int64_t>>& query_id_to_info,
                   const unordered_map<string, StableIntervalTree> target_to_interval_tree) {

    const string& target_name = paf_toks[5];
    int64_t target_start = stol(paf_toks[7]);
    int64_t target_end = stol(paf_toks[8]);
    bool is_reverse = paf_toks[4] == "-";
    // todo: reverse strand!!! (ugh)
    assert(!is_reverse);

    // find the mapping of our target sequence to query sequence(s)
    const StableIntervalTree& interval_tree = target_to_interval_tree.at(target_name);

    int64_t query_pos = stol(paf_toks[2]);
    int64_t target_pos = stol(paf_toks[7]);

    vector<pair<string, string>> cigars;
    for (int i = 12; i < paf_toks.size(); ++i) {
        if (paf_toks[i].substr(0, 5) == "cg:Z:") {
            for_each_cg(paf_toks[i], [&](const string& val, const string& cat) {
                    cigars.push_back(make_pair(val, cat));
                });
        }
    }

    for (const auto& vc : cigars) {
        const string& val = vc.first;
        const string& cat = vc.second;
        int64_t len = stol(val);
        if (cat == "M") {

            vector<StableInterval> overlapping_intervals = interval_tree.findOverlapping(target_pos, target_pos + len - 1);

            // these are the non-overlapping intervals span the target (but don't necessarily cover it completely)
            vector<pair<int64_t, int64_t>> target_intervals;
            int64_t cur_target_offset = target_start;
                        
            // assumption: the overlapping intervals are returned in sorted order
            for (StableInterval& overlapping_interval : overlapping_intervals) {
                int64_t clip_start = max(overlapping_interval.start, cur_target_offset);
                int64_t clip_stop = min(overlapping_interval.stop, target_end - 1);
                if (clip_stop >= clip_start) {
                    // write the new paf line for the subinterval
                    make_paf_line_for_interval(paf_toks, query_id_to_info, query_pos, target_pos, len,
                                               overlapping_interval, clip_start, clip_stop);
                    cur_target_offset = clip_stop;
                }
            }

        }
        if (cat == "M") {
            query_pos += len;
            target_pos += len;
        } else if (cat == "I") {
            query_pos += len;
        } else if (cat == "D") {
            target_pos += len;
        } else {
            assert(false);
        }
    }
}

void make_paf_line_for_interval(const vector<string>& paf_toks,
                                const vector<pair<string, int64_t>>& query_id_to_info,
                                int64_t match_start_query,
                                int64_t match_start_target,
                                int64_t match_length,
                                const StableInterval& overlapping_interval,
                                int64_t target_start, 
                                int64_t target_stop) {

#ifdef debug
    cerr << "make sub interval msq=" << match_start_query << " mst=" << match_start_target << " len=" << match_length
         << " ols=" << overlapping_interval.start << " olp=" << overlapping_interval.stop << " ts=" << target_start
         << " tp=" << target_stop << endl;
#endif

    // offset with respect to query_pos and target_pos (which indicate the beginning of the Match block)
    int64_t delta = target_start - match_start_target;
    assert(delta >= 0);

    int64_t output_block_length = target_stop - target_start + 1;
    
    int64_t output_query_start = match_start_query + delta; 
    int64_t output_query_end = output_query_start + output_block_length;

    // pull out the mapping information from the interval
    // this is where the target interval ends up in the stable sequence
    const pair<string, int64_t>& mapped_interval_info = query_id_to_info.at(get<0>(overlapping_interval.value));
    const int64_t& mapped_interval_start = get<1>(overlapping_interval.value);
    const bool& mapped_interval_reversed = get<2>(overlapping_interval.value);

    int64_t output_target_start = mapped_interval_start + delta;
    int64_t output_target_end = output_target_start + output_block_length;
    
    cout << paf_toks[0] << "\t"
         << paf_toks[1] << "\t"
         << output_query_start << "\t"
         << output_query_end << "\t"
         << "+" << "\t"
         << mapped_interval_info.first << "\t"
         << mapped_interval_info.second << "\t"
         << output_target_start << "\t"
         << output_target_end << "\t"
         << output_block_length << "\t"
         << output_block_length << "\t"
         << paf_toks[11] << "\t"
         << "cg:Z:" << output_block_length << "M"
         << "\n";
}

