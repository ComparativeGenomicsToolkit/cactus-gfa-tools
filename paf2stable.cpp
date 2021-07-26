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
    
    for (int i = 12; i < paf_toks.size(); ++i) {
        if (paf_toks[i].substr(0, 5) == "cg:Z:") {
            for_each_cg(paf_toks[i], [&](const string& val, const string& cat) {
                    int64_t len = stol(val);
                    if (cat == "M") {
#ifdef debug
                        cerr << "adding interval for " << target_name << ": " << target_pos << "-" << (target_pos + len)
                             << " ==> " << query_name << ": " << query_pos << ", rev=" << is_reverse << endl;
#endif             
                        StableInterval interval(target_pos, target_pos + len,
                                                make_tuple(query_id, query_pos, is_reverse));
                        target_intervals.push_back(interval);
                    }
                    if (cat == "M" || cat == "I") {
                        query_pos += len;
                    } else if (cat == "D") {
                        target_pos += len;
                    } else {
                        assert(false);
                    }
                });
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
                   const unordered_map<string, int64_t>& query_name_to_id,
                   const vector<pair<string, int64_t>>& query_id_to_info,
                   const unordered_map<string, StableIntervalTree> target_to_interval_tree) {
}
