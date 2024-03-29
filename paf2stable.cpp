#include "paf2stable.hpp"
#include "pafcoverage.hpp"

//#define debug

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
   
    int64_t query_start = stol(paf_toks[2]);
    int64_t target_start = stol(paf_toks[7]);
    int64_t target_end = stol(paf_toks[8]);
    
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

    int64_t target_offset = 0;
    int64_t query_offset = 0;
    int64_t query_pos;
    int64_t target_pos;
    for (const auto& vc : cigars) {
        const string& val = vc.first;
        const string& cat = vc.second;
        int64_t len = stol(val);
        if (cat == "M") {
            query_pos = query_start + query_offset;
            // if we're in reverse coordinates, we need to measure from the end
            if (is_reverse) {
                target_pos = target_end - len - target_offset;
            } else {
                target_pos = target_start + target_offset;
            }
#ifdef debug
            cerr << "adding interval for " << target_name << ": " << target_pos << "-" << (target_pos + len - 1)
                 << " ==> " << query_name << ": " << query_pos << ", rev=" << is_reverse << endl;
#endif
            StableInterval interval(target_pos, target_pos + len - 1,
                                    make_tuple(query_id, query_pos, is_reverse));
            target_intervals.push_back(interval);
            
            query_offset += len;
            target_offset += len;
        } else if (cat == "I") {
            query_offset += len;
        } else if (cat == "D") {
            target_offset += len;
        } else {
            assert(false);
        }
    }
}

void create_interval_trees(unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals) {
    for (auto& kv : target_to_intervals) {
        vector<StableInterval>& intervals = kv.second.second;
        int64_t target_size = kv.second.first;

        // sort on start coordinate, which is good enough for our soon-to-be non-overlapping intervals
        // we break ties using query name to avoid fragmentation (ie always make same choice in dupe sets
        // when possible)
        function<bool(const StableInterval&, const StableInterval&)> stable_interval_less =
            [](const StableInterval& i1, const StableInterval& i2) {
            if (i1.start == i2.start) {
                return get<0>(i1.value) < get<0>(i2.value);
            } else {
                return i1.start < i2.start;
            }
        };
        
        // take a quick pass just to remove duplicates, as there are often many
        std::sort(intervals.begin(), intervals.end(), stable_interval_less);
        vector<StableInterval> unique_intervals;
        for (size_t i = 0; i < intervals.size(); ++i) {
            if (i == 0 || intervals[i].start != intervals[i-1].start || intervals[i].stop != intervals[i-1].stop) {
                unique_intervals.push_back(intervals.at(i));
            }
        }
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
            clip_interval(interval, cut_points, clipped_intervals);
        }
        intervals = std::move(clipped_intervals);        

        // take a quick pass just to remove duplicates after clipping
        std::sort(intervals.begin(), intervals.end(), stable_interval_less);
        unique_intervals.clear();
        for (size_t i = 0; i < intervals.size(); ++i) {
            if (i == 0 || intervals[i].start != intervals[i-1].start || intervals[i].stop != intervals[i].stop) {
                unique_intervals.push_back(intervals.at(i));
            }
        }        
        intervals = std::move(unique_intervals);

#ifdef debug
        cerr << "Interval Tree (" << kv.first << "):";
        for (auto interval : intervals) {
            cerr << endl << "   " << interval;
        }
        cerr << endl;
#endif
    }
}

void clip_interval(const StableInterval& interval,
                   const set<int64_t>& cut_points, vector<StableInterval>& clipped_intervals) {

    assert(interval.stop >= interval.start);
    if (interval.stop == interval.start) {
        // can't cut a size-1 interval
        clipped_intervals.push_back(interval);
        return;
    }

    // fish relevant cut points out of the set (todo: use iterators throughout)
    auto i = cut_points.lower_bound(interval.start);
    auto j = cut_points.upper_bound(interval.stop - 1);
    vector<int64_t> cut_positions;
    for (auto k = i; k != j; ++k) {
        cut_positions.push_back(*k);
    }
#ifdef debug
    cerr << "query cuts on " << interval << " gives " << cut_positions.size() << " results" << endl;
#endif
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

    int64_t interval_size = interval.stop - interval.start + 1;

    // make the clipped intervals
    // todo: just do this in above pass
    bool is_reverse = get<2>(interval.value);
    for (auto& ni : new_intervals) {
        int64_t stable_offset;
        if (is_reverse) {
            stable_offset = get<1>(interval.value) + interval_size - 1 - (ni.first - interval.start) - (ni.second - ni.first);
        } else {
            stable_offset = get<1>(interval.value) + (ni.first - interval.start);
        }
        clipped_intervals.emplace_back(ni.first, ni.second, make_tuple(get<0>(interval.value),
                                                                       stable_offset,
                                                                       is_reverse));
        assert(clipped_intervals.back().stop >= clipped_intervals.back().start);
#ifdef debug
        cerr << "Adding clipped " << clipped_intervals.back() << endl;
#endif
    }
}


size_t paf_to_stable(const vector<string>& paf_toks,
                     const vector<pair<string, int64_t>>& query_id_to_info,
                     const unordered_map<string, pair<int64_t, vector<StableInterval>>>& target_to_intervals) {

    int64_t query_start = stol(paf_toks[2]);
    const string& target_name = paf_toks[5];
    int64_t target_size = stol(paf_toks[6]);
    int64_t target_start = stol(paf_toks[7]);
    int64_t target_end = stol(paf_toks[8]);
    bool is_reverse = paf_toks[4] == "-";
    size_t lines_written = 0;

    // find the mapping of our target sequence to query sequence(s)
    const vector<StableInterval>& intervals = target_to_intervals.at(target_name).second;

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
    int64_t target_offset = 0;
    int64_t query_offset = 0;
    int64_t query_pos;
    int64_t target_pos;
    vector<string> out_toks_buffer;
    for (const auto& vc : cigars) {
        const string& val = vc.first;
        const string& cat = vc.second;
        int64_t len = stol(val);
        if (cat == "M") {

            query_pos = query_start + query_offset;
            // if we're in reverse coordinates, we need to measure from the end
            if (is_reverse) {
                target_pos = target_end - len - target_offset;
            } else {
                target_pos = target_start + target_offset;
            }

            // pull all overlapping intervals out of our sorted list
            StableInterval query_interval(target_pos, target_pos, make_tuple(0, 0, false));
            auto lb = std::lower_bound(intervals.begin(), intervals.end(), query_interval, StableIntervalTree::IntervalStartCmp());
            query_interval = StableInterval(target_pos + len - 1, target_pos + len - 1, make_tuple(0, 0, false));
            auto ub = std::upper_bound(intervals.begin(), intervals.end(), query_interval, StableIntervalTree::IntervalStartCmp());
            vector<StableInterval> overlapping_intervals(lb, ub);
            
            // these intervals must, by definition, exactly cover the whole match block (because we built and clipped on
            // every match block in the set
            assert(!overlapping_intervals.empty());
            assert(overlapping_intervals[0].start == target_pos);
            assert(overlapping_intervals.back().stop == target_pos + len - 1);

            if (is_reverse) {
                std::reverse(overlapping_intervals.begin(), overlapping_intervals.end());
            }

            int64_t total_block_length = 0;
            for (size_t i = 0; i < overlapping_intervals.size(); ++i) {
                StableInterval& overlapping_interval = overlapping_intervals[i];
                if (i > 0 && !is_reverse) {
                    // expect exact coverage (see above)
                    assert(overlapping_interval.start == overlapping_intervals[i-1].stop + 1);
                }
                lines_written += make_paf_line_for_interval(paf_toks, query_id_to_info, overlapping_interval, query_pos + total_block_length, out_toks_buffer);
                total_block_length += overlapping_interval.stop - overlapping_interval.start + 1;
            }
            assert(total_block_length == len);
            query_offset += len;
            target_offset += len;
        } else if (cat == "I") {
            query_offset += len;
        } else if (cat == "D") {
            target_offset += len;
        } else {
            assert(false);
        }
    }
    // flush the buffer
    if (!out_toks_buffer.empty()) {
        for (size_t i = 0; i < out_toks_buffer.size(); ++i) {
            cout << out_toks_buffer[i] << (i < out_toks_buffer.size() - 1 ? "\t" : "\n");
        }
        ++lines_written;
    }            
    return lines_written;
}

size_t make_paf_line_for_interval(const vector<string>& paf_toks,
                                  const vector<pair<string, int64_t>>& query_id_to_info,
                                  const StableInterval& overlapping_interval,
                                  int64_t query_pos,
                                  vector<string>& out_toks_buffer) {

    // pull out the mapping information from the interval
    // this is where the target interval ends up in the stable sequence
    const pair<string, int64_t>& mapped_interval_info = query_id_to_info.at(get<0>(overlapping_interval.value));
    const int64_t& mapped_interval_start = get<1>(overlapping_interval.value);
    const bool& mapped_interval_reversed = get<2>(overlapping_interval.value);

    int64_t block_length = overlapping_interval.stop - overlapping_interval.start + 1;

    bool is_reverse = mapped_interval_reversed != (paf_toks[4] == "-");

    vector<string> paf_line =  {
        paf_toks[0],
        paf_toks[1],
        to_string(query_pos),
        to_string(query_pos + block_length),
        (is_reverse ? "-" : "+"),
        mapped_interval_info.first,
        to_string(mapped_interval_info.second),
        to_string(mapped_interval_start),
        to_string(mapped_interval_start + block_length),
        to_string(block_length),
        to_string(block_length),
        paf_toks[11],
        "cg:Z:" + to_string(block_length) + "M"
    };

    bool merged = false;

    if (!out_toks_buffer.empty() && paf_line[4] == out_toks_buffer[4] && paf_line[5] == out_toks_buffer[5]) {
        // queries should never change
        assert(paf_line[0] == out_toks_buffer[0]);
        assert(paf_line[11] == out_toks_buffer[11]);
        // todo: a lot of avoidable string <==> integer conversion going on, but keeps refactor simple for now
        int64_t prev_query_start = stol(out_toks_buffer[2]);
        int64_t prev_query_end = stol(out_toks_buffer[3]);
        int64_t prev_target_start = stol(out_toks_buffer[7]);
        int64_t prev_target_end = stol(out_toks_buffer[8]);
        int64_t cur_query_start = stol(paf_line[2]);
        int64_t cur_query_end = stol(paf_line[3]);
        int64_t cur_target_start = stol(paf_line[7]);
        int64_t cur_target_end = stol(paf_line[8]);

        if (paf_line[4] == "+" && prev_query_end <= cur_query_start && prev_target_end <= cur_target_start) {
            // we merge paf_line into out_toks_buffer on the forward strand
            int64_t query_delta = cur_query_start - prev_query_end;
            int64_t target_delta = cur_target_start - prev_target_end;
            out_toks_buffer[3] = paf_line[3];
            out_toks_buffer[8] = paf_line[8];
            out_toks_buffer[9] = to_string(cur_query_end - prev_query_start); // not exact, but not used
            out_toks_buffer[10] = to_string(stol(out_toks_buffer[10]) + (cur_query_end - cur_query_start));
            if (query_delta > 0) {
                out_toks_buffer[12] += to_string(query_delta) + "I";
            }
            if (target_delta > 0) {
                out_toks_buffer[12] += to_string(target_delta) + "D";
            }
            out_toks_buffer[12] += to_string(cur_query_end - cur_query_start) + "M";
            merged = true;
            
        } else if (paf_line[4] == "-" && prev_query_end <= cur_query_start && cur_target_end <= prev_target_start) {
            // we merge paf_line into out_toks_buffer on the reverse strand            
            int64_t query_delta = cur_query_start - prev_query_end;
            int64_t target_delta = prev_target_start - cur_target_end;
            out_toks_buffer[3] = paf_line[3];
            out_toks_buffer[7] = paf_line[7];
            out_toks_buffer[9] = to_string(cur_query_end - prev_query_start); // not exact, but not used
            out_toks_buffer[10] = to_string(stol(out_toks_buffer[10]) + (cur_query_end - cur_query_start));
            if (query_delta > 0) {
                paf_line[12] += to_string(query_delta) + "I";
            }
            if (target_delta > 0) {
                paf_line[12] += to_string(target_delta) + "D";
            }            
            out_toks_buffer[12] = paf_line[12] + out_toks_buffer[12].substr(5);
            merged = true;            
        }
    }

    size_t ret = 0;
    if (!merged) {
        if (!out_toks_buffer.empty()) {
            for (size_t i = 0; i < out_toks_buffer.size(); ++i) {
                cout << out_toks_buffer[i] << (i < out_toks_buffer.size() - 1 ? "\t" : "\n");
            }
            ++ret;
        }
        out_toks_buffer = paf_line;
    }
    return ret;
}

