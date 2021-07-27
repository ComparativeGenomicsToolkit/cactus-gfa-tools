#include "paf2stable.hpp"
#include "pafcoverage.hpp"

#define debug

void update_stable_mapping_info(const vector<string>& paf_toks,
                                unordered_map<string, int64_t>& query_name_to_id,
                                vector<pair<string, int64_t>>& query_id_to_info,
                                unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals) {

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
    pair<int64_t, vector<StableInterval>>& target_info = target_to_intervals[target_name];
    target_info.first = stol(paf_toks[6]);
    vector<StableInterval>& target_intervals = target_info.second;

    bool is_reverse = paf_toks[4] == "-";
   
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

    if (is_reverse) {
        std::reverse(cigars.begin(), cigars.end());
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

unordered_map<string, StableIntervalTree> create_interval_trees(unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals) {
    unordered_map<string, StableIntervalTree> target_to_interval_trees;
    for (auto& kv : target_to_intervals) {
        vector<StableInterval>& intervals = kv.second.second;
        int64_t target_size = kv.second.first;

        // take a quick pass just to remove duplicates, as there are often many
        std::sort(intervals.begin(), intervals.end(), StableIntervalTree::IntervalStartCmp());
        vector<StableInterval> unique_intervals;
        for (size_t i = 0; i < intervals.size(); ++i) {
            if (i == 0 || intervals[i].start != intervals[i-1].start || intervals[i].stop != intervals[i-1].stop) {
                unique_intervals.push_back(intervals.at(i));
            }
        }
        cerr << "Unique filter reduces from " << intervals.size() << " to " << unique_intervals.size() << endl;
        intervals = std::move(unique_intervals);
        
        // cut at all interval ends (a cut point cuts to the *right* of its position)
        set<int64_t> cut_points;
        for (StableInterval& interval : intervals) {
            if (interval.start > 0) {
                cut_points.insert(interval.start - 1);
            }
            if (interval.stop < target_size - 1) {
                cut_points.insert(interval.stop);
            }
        }
        vector<StableInterval> clipped_intervals;
        for (StableInterval& interval : intervals) {
            clip_interval(interval, target_size, cut_points, clipped_intervals);
        }
        cerr << "Clipped filter reduces from " << intervals.size() << " to " << clipped_intervals.size() << endl;
        intervals = std::move(clipped_intervals);        

        // take a quick pass just to remove duplicates after clipping
        std::sort(intervals.begin(), intervals.end(), StableIntervalTree::IntervalStartCmp());
        unique_intervals.clear();
        for (size_t i = 0; i < intervals.size(); ++i) {
            if (i == 0 || intervals[i].start != intervals[i-1].start || intervals[i].stop != intervals[i].stop) {
                unique_intervals.push_back(intervals.at(i));
            }
        }
        cerr << "Unique filter2 reduces from " << intervals.size() << " to " << unique_intervals.size() << endl;
        for (auto xx : unique_intervals) {
            cerr << xx << endl;
        }
        intervals = std::move(unique_intervals);

#ifdef debug
        cerr << "Interval Tree (" << kv.first << "):";
        for (auto interval : intervals) {
            cerr << endl << "   " << interval;
        }
        cerr << endl;
#endif
        target_to_interval_trees[kv.first] = StableIntervalTree(intervals);
        intervals.clear();

    }
    target_to_intervals.clear();
    return target_to_interval_trees;
}

void clip_interval(const StableInterval& interval, int64_t target_size,
                   const set<int64_t>& cut_points, vector<StableInterval>& clipped_intervals) {

    assert(interval.stop >= interval.start);
    if (interval.stop == interval.start) {
        // can't cut a size-1 interval
        clipped_intervals.push_back(interval);
        return;
    }

    // fish relevant cut points out of the set (todo: use iterators throughotu)
    auto i = cut_points.lower_bound(interval.start);
    auto j = cut_points.upper_bound(interval.stop - 1);
    vector<int64_t> cut_positions;
    for (auto k = i; k != j; ++k) {
        cut_positions.push_back(*k);
    }
    cerr << "qyuery cuts on " << interval << " gives " << cut_positions.size() << " results" << endl;
    // make sure last point is a cut point no matter what so
    // we can clip everything in a loop
    if (cut_positions.empty() || cut_positions.back() != interval.stop) {
        cut_positions.push_back(interval.stop);
    }
    // pull out the cut intervals
    vector<pair<int64_t, int64_t>> new_intervals;
    int64_t cur = interval.start;
    for (int64_t cp : cut_positions) {
        assert(cur <= cp);
        new_intervals.push_back(make_pair(cur, cp));
        cur = cp+1;
    }
    if (new_intervals.size() == 1 && new_intervals[0].first == interval.start && new_intervals[0].second == interval.stop) {
        // nothing touched, just return the interval
        clipped_intervals.push_back(interval);
        return;
    }

    // make the clipped intervals
    // todo: just do this in above pass
    bool is_reverse = get<2>(interval.value);
    for (auto& ni : new_intervals) {
        int64_t stable_offset;
        if (is_reverse) {
            stable_offset = get<1>(interval.value) + target_size - 1 - ni.second;
        } else {
            stable_offset = get<1>(interval.value) + ni.first;
        }
        clipped_intervals.emplace_back(ni.first, ni.second, make_tuple(get<0>(interval.value),
                                                                       stable_offset,
                                                                       is_reverse));
        assert(clipped_intervals.back().stop >= clipped_intervals.back().start);
    }
}


void paf_to_stable(const vector<string>& paf_toks,
                   const vector<pair<string, int64_t>>& query_id_to_info,
                   const unordered_map<string, StableIntervalTree> target_to_interval_tree) {

    const string& target_name = paf_toks[5];
    int64_t target_size = stol(paf_toks[6]);
    int64_t target_start = stol(paf_toks[7]);
    int64_t target_end = stol(paf_toks[8]);
    bool is_reverse = paf_toks[4] == "-";

    // find the mapping of our target sequence to query sequence(s)
    const StableIntervalTree& interval_tree = target_to_interval_tree.at(target_name);

    int64_t query_pos = stol(paf_toks[2]);
    int64_t target_pos = target_start;

    vector<pair<string, string>> cigars;
    for (int i = 12; i < paf_toks.size(); ++i) {
        if (paf_toks[i].substr(0, 5) == "cg:Z:") {
            for_each_cg(paf_toks[i], [&](const string& val, const string& cat) {
                    cigars.push_back(make_pair(val, cat));
                });
        }
    }

    if (is_reverse) {
        std::reverse(cigars.begin(), cigars.end());
    }
    
    for (const auto& vc : cigars) {
        const string& val = vc.first;
        const string& cat = vc.second;
        int64_t len = stol(val);
        if (cat == "M") {

            vector<StableInterval> overlapping_intervals = interval_tree.findOverlapping(target_pos, target_pos + len - 1);
            // we mostly do this to make it easier to debug
            std::sort(overlapping_intervals.begin(), overlapping_intervals.end(), StableIntervalTree::IntervalStartCmp());

            // these intervals must, by definition, exactly cover the whole match block (because we built and clipped on
            // every match block in the set
            assert(!overlapping_intervals.empty());
            assert(overlapping_intervals[0].start == target_pos);
            assert(overlapping_intervals.back().stop == target_pos + len - 1);

            // these are the non-overlapping intervals span the target (but don't necessarily cover it completely)
            vector<pair<int64_t, int64_t>> target_intervals;
            int64_t cur_target_offset = target_start;
                        
            // assumption: the overlapping intervals are returned in sorted order
            for (size_t i = 0; i < overlapping_intervals.size(); ++i) {
                StableInterval& overlapping_interval = overlapping_intervals[i];
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

    bool is_reverse = mapped_interval_reversed != (paf_toks[4] == "-");
    
    cout << paf_toks[0] << "\t"
         << paf_toks[1] << "\t"
         << output_query_start << "\t"
         << output_query_end << "\t"
         << (is_reverse ? "-" : "+") << "\t"
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
