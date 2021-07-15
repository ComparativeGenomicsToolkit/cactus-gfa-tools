/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */

#include "mzgaf2paf.hpp"
#include "mzgaf.hpp"
#include "pafcoverage.hpp"
#include "gfakluge.hpp"
#include <limits>
#include <fstream>

//#define debug

using namespace gafkluge;
using namespace std;

struct MatchBlock {
    int64_t query_start;
    int64_t query_end;
    int64_t target_start;
    int64_t target_end; // could just store length instead of ends, but use to sanity check
};

size_t mzgaf2paf(const MzGafRecord& gaf_record,
                 const GafRecord& parent_record,
                 ostream& paf_stream,
                 int64_t min_gap,
                 int64_t min_match_length,
                 MZMap& mz_map,
                 double universal_filter,
                 QueryCoverage& query_coverage,
                 int64_t min_overlap_len,                 
                 const string& target_prefix,
                 unordered_map<string, tuple<string, int64_t, int64_t>>& stable_lookup) {

    // paf coordinates are always on forward strand. but the mz output coordinates for the target
    // can apparently be on the reverse strand, so we flip them as needed
    int64_t paf_target_start = gaf_record.target_start;
    int64_t paf_target_end = gaf_record.target_end;
    if (gaf_record.is_reverse) {
        paf_target_start = gaf_record.target_length - gaf_record.target_end;
        paf_target_end = gaf_record.target_length - gaf_record.target_start;
    }
  
    assert(gaf_record.target_mz_offsets.size() == gaf_record.query_mz_offsets.size());

    // if we have a mz map, determine the acceptable count for each minimzer to consider "universal"
    MZCount* mz_counts = nullptr;
    if (universal_filter > 0) {
        auto it = mz_map.find(gaf_record.target_name);
        assert(it != mz_map.end());
        mz_counts = &it->second;
    }

    TwoBitVec* cov_vec = min_overlap_len > 0 ? &query_coverage[parent_record.query_name] : nullptr;
    if (cov_vec && cov_vec->size() == 0) {
        assert(parent_record.block_length < min_overlap_len);
        cov_vec = nullptr;
    }
    
    // turn the offsets vectors into match blocks, applying gap and inconsistency filters
    
    // positions as we move along (relative to gaf_record.query/target_starts)
    int64_t query_pos = 0;
    int64_t target_pos = 0;

    // cigar matches
    vector<MatchBlock> matches;
    matches.reserve(gaf_record.num_minimizers);

    for (size_t i = 0; i < gaf_record.num_minimizers; ++i) {

        MatchBlock match = {query_pos, query_pos + gaf_record.kmer_size,
                            target_pos, target_pos + gaf_record.kmer_size};

        // optional universal check
        bool universal = true;
        if (mz_counts != nullptr) {
            // index on forward strand
            int64_t mz_idx = !gaf_record.is_reverse ? gaf_record.target_start + target_pos :
                gaf_record.target_length - gaf_record.target_start - target_pos - gaf_record.kmer_size;
            assert(universal == 1 || (mz_counts->at(mz_idx).first > 0 && mz_counts->at(mz_idx).second > 0));
            // proportion of mapped sequences that have this exact minimizer
            // todo: this is sample-unaware.  so if two sequences in a given sample have this
            // minimizer, it's okay as far as the filter's concerned.
            // update: this is sample-aware if the samples are different files and u==1
            // by way of lowering the numerator in the combine_mz funciton below
            float mz_frac = (float)mz_counts->at(mz_idx).first / (float)mz_counts->at(mz_idx).second;
            // mz_frac can be > 1 due to 0-offsets in the minigraph mz lists.  
            universal = mz_frac >= universal_filter && mz_frac <= 1.;
        }

        // optional mz_overlap check: filter out minimizer if in query region covered more than once,
        // or if the query region is covered once and we're too small
        if (cov_vec != nullptr) {
            // todo: this is all pretty slow (borrowing the per-base code of universal filter)
            // may have to switch to interval queries at scale
            for (int64_t i = match.query_start; i < match.query_end && universal; ++i) {
                size_t coverage = cov_vec->get(gaf_record.query_start + i);
                if (coverage > 1 || (coverage == 1 && parent_record.block_length < min_overlap_len)) {
                    // we hijack the unversal flag, as the filters have the same effect
                    universal = false;
                }
            }
        }

        if (matches.empty()) {
            if (universal) {
                matches.push_back(match);
            }
        } else {
            
            // compute the overlap with the previous minimizer
            int64_t query_delta = match.query_start - matches.back().query_end;
            int64_t target_delta = match.target_start - matches.back().target_end;

            if (query_delta == target_delta && query_delta <= 0) {
                // extend adjacent minimizers (todo: do we want to always do this or control with option?)
                if (universal) {
                    matches.back().query_end = match.query_end;
                    matches.back().target_end = match.target_end;
                }
            } else if (query_delta < 0 || target_delta < 0) {
                // drop inconsistent minimizers
                // todo: do we want to accept universal minimizers even if inconsistent with non-universal minimizers
                //       which we'd drop?
                matches.pop_back();
            } else if (query_delta >= min_gap && target_delta >= min_gap) {
                // add if passes gap filter
                if (universal) {
                    // apply length filter
                    if (min_match_length > 0 && !matches.empty() && matches.back().query_end - matches.back().query_start < min_match_length) {
                        matches.pop_back();
                    }
                    matches.push_back(match);
                }
            }
        }

        // advance our position
        if (i < gaf_record.num_minimizers - 1) {
            query_pos += gaf_record.query_mz_offsets[i];
            target_pos += gaf_record.target_mz_offsets[i];
        }
    }
    
    // apply length filter
    if (min_match_length > 0 && !matches.empty() && matches.back().query_end - matches.back().query_start < min_match_length) {
        matches.pop_back();
    }

    // turn our matches into cigar (todo: go straight to string!)
    vector<string> cigar;
    cigar.reserve(matches.size() * 3 + 2);

    // need for output
    int64_t total_matches = 0;
    // used for sanity checks
    int64_t total_deletions = 0;
    int64_t total_insertions = 0;

    // don't allow cigars to start on indel, as sometimes they can really inflate block
    // lengths (big softclips can come when hardmasking centromeres going into minigraph)
    int64_t leading_insertions = 0;
    int64_t leading_deletions = 0;
    if (!matches.empty() && matches[0].query_start > 0) {
        total_insertions += matches[0].query_start;
        leading_insertions = matches[0].query_start;
    }
    if (!matches.empty() && matches[0].target_start > 0) {
        total_deletions += matches[0].target_start;
        leading_deletions = matches[0].target_start; 
    }
    
    for (size_t i = 0; i < matches.size(); ++i) {
        // match
        int64_t match_size = matches[i].query_end - matches[i].query_start;
        assert(match_size == matches[i].target_end - matches[i].target_start);
        cigar.push_back(std::to_string(match_size) + "M");
        total_matches += match_size;
        if (i < matches.size() - 1) {
            // insertion before next match
            int64_t insertion_size = matches[i+1].query_start - matches[i].query_end;
            assert(insertion_size >= min_gap);
            cigar.push_back(std::to_string(insertion_size) + "I");
            total_insertions += insertion_size;
            // deletion before next match
            int64_t deletion_size = matches[i+1].target_start - matches[i].target_end;
            assert(deletion_size >= min_gap);
            cigar.push_back(std::to_string(deletion_size) + "D");
            total_deletions += deletion_size;
        }
    }

    // don't allow cigars to end on indel, as sometimes they can really inflate block
    // lengths (big softclips can come when hardmasking centromeres going into minigraph)    
    int64_t query_length = gaf_record.query_end - gaf_record.query_start;
    int64_t leftover_insertions = query_length - (total_insertions + total_matches);
    if (leftover_insertions) {
        assert(matches.empty() || matches.back().query_end + leftover_insertions == query_length);
    }
    int64_t target_length = gaf_record.target_end - gaf_record.target_start;
    int64_t leftover_deletions = target_length - (total_deletions + total_matches);
    if (leftover_deletions) {
        assert(matches.empty() || matches.back().target_end + leftover_deletions == target_length);
    }
    assert(leftover_insertions >= 0 && leftover_deletions >= 0);
    if (gaf_record.is_reverse) {
        swap(leading_deletions, leftover_deletions);
    }

    // optional stable coordinate conversion
    string stable_target_name = gaf_record.target_name;
    int64_t stable_target_length = gaf_record.target_length;
    int64_t stable_offset = 0;    
    if (!stable_lookup.empty()) {
        tuple<string, int64_t, int64_t>& stable_info = stable_lookup.at(gaf_record.target_name);
        stable_target_name = get<0>(stable_info);
        stable_target_length = get<2>(stable_info);
        stable_offset =  get<1>(stable_info);
    }

    if (!matches.empty()) {
        // output the paf columns
        paf_stream << parent_record.query_name << "\t"
                   << parent_record.query_length << "\t"
                   << (gaf_record.query_start + leading_insertions) << "\t"
                   << (gaf_record.query_end - leftover_insertions) << "\t"
                   << (gaf_record.is_reverse ? "-" : "+") << "\t"
                   << target_prefix << stable_target_name << "\t"
                   << stable_target_length << "\t"
                   << (paf_target_start + leading_deletions + stable_offset) << "\t"
                   << (paf_target_end - leftover_deletions + stable_offset) << "\t"
                   << total_matches << "\t"
                   << (gaf_record.target_end - gaf_record.target_start - leftover_deletions - leading_deletions) << "\t" // fudged
                   << parent_record.mapq << "\t" << "cg:Z:";

        // and the cigar
        if (gaf_record.is_reverse) {
            for (vector<string>::reverse_iterator ci = cigar.rbegin(); ci != cigar.rend(); ++ci) {
                paf_stream << *ci;
            }
        } else {
            for (vector<string>::iterator ci = cigar.begin(); ci != cigar.end(); ++ci) {
                paf_stream << *ci;
            }
        }

        paf_stream << "\n";
    }

    return total_matches;
}


// update the counts for one mapping of query to target
void update_mz_map(const gafkluge::MzGafRecord& gaf_record,
                   const gafkluge::GafRecord& parent_record,
                   MZMap& mz_map,
                   int64_t min_mapq,
                   int64_t min_block_len,
                   int64_t min_node_len,
                   bool node_based_universal) {

    // find our target in the map
    MZCount& mz_counts = mz_map[gaf_record.target_name];

    // init the count array
    if (mz_counts.empty()) {
        mz_counts.resize(gaf_record.target_length, make_pair(0, 0));
    }

    // increment the mapping counts (regardless of thresholds)
    int64_t paf_target_start = gaf_record.target_start;
    int64_t paf_target_end = gaf_record.target_end;
    if (gaf_record.is_reverse) {
        paf_target_start = gaf_record.target_length - gaf_record.target_end;
        paf_target_end = gaf_record.target_length - gaf_record.target_start;
    }
    int64_t range_start = paf_target_start;
    int64_t range_end = paf_target_end;
    if (node_based_universal) {
        // todo: if we want to do things node based, we can use a much simpler structure that stores
        // one coutn per node.  keeping range just to preserve flexibility for now
        range_start = 0;
        range_end = gaf_record.target_length;
    }
    for (int64_t i = range_start; i < range_end; ++i) {
        ++mz_counts[i].second;        
    }

    // increment each minimzer (if the record passes the thresholds)
    if (gaf_record.num_minimizers > 0 &&
        parent_record.mapq >= min_mapq &&
        (parent_record.query_length <= min_block_len || parent_record.block_length >= min_block_len) &&
        gaf_record.target_length >= min_node_len) {
        
        int64_t target_pos = gaf_record.target_start;

        for (size_t i = 0; i < gaf_record.num_minimizers; ++i) {

            int64_t mz_idx = !gaf_record.is_reverse ? target_pos :
                gaf_record.target_length - target_pos - gaf_record.kmer_size;

            // note, mz_counts.first can be higher than .second because 0-offsets are apparently a thing in the
            // minigraph output
            ++mz_counts[mz_idx].first;

#ifdef debug
            cerr << (gaf_record.is_reverse ? "r" : "f") << "-increment i= " << i << " mzi=" << mz_idx << " to="
                 << mz_counts[mz_idx].first << "," << mz_counts[mz_idx].second << endl;
#endif
            // advance our position
            if (i < gaf_record.num_minimizers - 1) {
                target_pos += gaf_record.target_mz_offsets[i];
            }
        }
    }
}

void combine_mz_maps(MZMap& map1, MZMap& map2, bool reset_multiple_counts_to_0) {
    for (auto& kv : map1) {

        MZCount& mz_counts1 = kv.second;
        MZCount& mz_counts2 = map2[kv.first];

        // init the count array
        if (mz_counts2.empty() && !mz_counts1.empty()) {
            mz_counts2.resize(mz_counts1.size());
        }

        // add the counts
        for (int64_t i = 0; i < mz_counts1.size(); ++i) {
            mz_counts2[i].first += mz_counts1[i].first;
            mz_counts2[i].second += mz_counts1[i].second;
            if (reset_multiple_counts_to_0 && (mz_counts1[i].first > 1 || mz_counts1[i].second > 1)) {
                // this minimizer was covered or was present more than once in the current GAF file
                // we set its hits to 0 (if reset_multple flag is set) 
                mz_counts2[i].first = 0;
            }
        }
        // remove from map1 to keep memory down
        mz_counts1.clear();   
    }
}

void update_query_coverage(const gafkluge::GafRecord& parent_record,
                           QueryCoverage& query_coverage) {
    
    TwoBitVec& cov_vec = query_coverage[parent_record.query_name];
    if (cov_vec.size() == 0) {
        cov_vec.resize(parent_record.query_length);
    }
    bool check = false;
    for (size_t i = parent_record.query_start; i < parent_record.query_end; ++i) {
        cov_vec.increment(i);
    }
}

unordered_map<string, tuple<string, int64_t, int64_t>> build_stable_lookup(const string& fa_table_path,
                                                                           const string& rgfa_path) {
    cerr << "building" << endl;
    // load in the fasta header table: original contig name -> (cut contig name, event, length)
    unordered_map<string, tuple<string, string, size_t>> fa_header_table;
    if (!fa_table_path.empty()) {
        string buffer;
        ifstream fa_lengths_stream(fa_table_path);
        if (!fa_lengths_stream) {
            cerr << "[mzgaf2paf] error: unable to open input lengths: " << fa_table_path << endl;
            exit(1);
        }
        while (getline(fa_lengths_stream, buffer)) {
            vector<string> toks;
            split_delims(buffer, "\t\n", toks);
            if (toks.size() == 4) {
                assert(fa_header_table.count(toks[0]) == 0);
                fa_header_table[toks[0]] = make_tuple(toks[1], toks[2], stol(toks[3]));
            } else if (!toks.empty()) {
                cerr << "[mzgaf2paf] error: Unable to parse fasta header table line: " << buffer << endl;
                exit(1);
            }
        }
        assert(!fa_header_table.empty());
    }

    // load in the rgfa information
    unordered_map<string, tuple<string, int64_t, int64_t>> lookup_table;
    function<void(const gfak::sequence_elem&)> visit_seq = [&](const gfak::sequence_elem& gfa_seq) {
        bool found_SN = false;
        bool found_SR = false;
        bool found_SO = false;
        int64_t rank;
        int64_t offset;
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
            }  else if (oe.key == "SO") {
                assert(found_SO == false);
                offset = stol(oe.val);
                assert(offset >= 0);                               
                found_SO = true;
            }
        }
        assert(found_SN);
        assert(found_SR);
        assert(found_SO);

        auto& header_info = fa_header_table.at(contig);
        string cactus_name = "id=" + get<1>(header_info) + "|" + get<0>(header_info);        
        
        lookup_table[gfa_seq.name] = make_tuple(cactus_name, offset, get<2>(header_info));
    };

    // load the GFA
    gfak::GFAKluge kluge;
    kluge.for_each_sequence_line_in_file(rgfa_path.c_str(), visit_seq);

    assert(!lookup_table.empty());
    return lookup_table;
}
