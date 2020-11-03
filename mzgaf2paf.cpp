/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */

#include "mzgaf2paf.hpp"
#include "mzgaf.hpp"
#include <limits>

//#define debug

using namespace gafkluge;
using namespace std;

struct MatchBlock {
    int64_t query_start;
    int64_t query_end;
    int64_t target_start;
    int64_t target_end; // could just store length instead of ends, but use to sanity check
};

void mzgaf2paf(const MzGafRecord& gaf_record,
               const GafRecord& parent_record,
               ostream& paf_stream,
               size_t min_gap,
               MZMap& mz_map,
               double universal_filter,
               const string& target_prefix) {


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
            assert(mz_counts->at(mz_idx).first > 0 && mz_counts->at(mz_idx).second > 0);
            // proportion of mapped sequences that have this exact minimizer
            // todo: this is sample-unaware.  so if two sequences in a given sample have this
            // minimizer, it's okay as far as the filter's concerned. 
            float mz_frac = (float)mz_counts->at(mz_idx).first / (float)mz_counts->at(mz_idx).second;
            // mz_frac can be > 1 due to 0-offsets in the minigraph mz lists.  
            universal = mz_frac >= universal_filter && mz_frac <= 1.;
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

    // turn our matches into cigar (todo: go straight to string!)
    vector<string> cigar;
    cigar.reserve(matches.size() * 3 + 2);

    // need for output
    int64_t total_matches = 0;
    // used for sanity checks
    int64_t total_deletions = 0;
    int64_t total_insertions = 0;

    // todo: kind of silly to start cigar with an indel -- should just clip
    if (!matches.empty() && matches[0].query_start > 0) {
        cigar.push_back(std::to_string(matches[0].query_start) + "I");
        total_insertions += matches[0].query_start;
    }
    if (!matches.empty() && matches[0].target_start > 0) {
        cigar.push_back(std::to_string(matches[0].target_start) + "D");
        total_deletions += matches[0].target_start;
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

    // todo: kind of silly to end cigar with an indel -- should just clip
    int64_t query_length = gaf_record.query_end - gaf_record.query_start;
    int64_t leftover_insertions = query_length - (total_insertions + total_matches);
    if (leftover_insertions) {
        assert(matches.empty() || matches.back().query_end + leftover_insertions == query_length);
        cigar.push_back(std::to_string(leftover_insertions) + "I");
    }
    int64_t target_length = gaf_record.target_end - gaf_record.target_start;
    int64_t leftover_deletions = target_length - (total_deletions + total_matches);
    if (leftover_deletions) {
        assert(matches.empty() || matches.back().target_end + leftover_deletions == target_length);
        cigar.push_back(std::to_string(leftover_deletions) + "D");
    }

    if (!matches.empty()) {
        // output the paf columns
        paf_stream << parent_record.query_name << "\t"
                   << parent_record.query_length << "\t"
                   << gaf_record.query_start << "\t"
                   << gaf_record.query_end << "\t"
                   << (gaf_record.is_reverse ? "-" : "+") << "\t"
                   << target_prefix << gaf_record.target_name << "\t"
                   << gaf_record.target_length << "\t"
                   << paf_target_start << "\t"
                   << paf_target_end << "\t"
                   << total_matches << "\t"
                   << (gaf_record.target_end - gaf_record.target_start) << "\t" // fudged
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
}


// update the counts for one mapping of query to target
void update_mz_map(const gafkluge::MzGafRecord& gaf_record,
                   const gafkluge::GafRecord& parent_record,
                   MZMap& mz_map) {

    // find our target in the map
    MZCount& mz_counts = mz_map[gaf_record.target_name];

    // init the count array
    if (mz_counts.empty()) {
        mz_counts.resize(gaf_record.target_length, make_pair(0, 0));
    }

    // increment the mapping counts
    int64_t paf_target_start = gaf_record.target_start;
    int64_t paf_target_end = gaf_record.target_end;
    if (gaf_record.is_reverse) {
        paf_target_start = gaf_record.target_length - gaf_record.target_end;
        paf_target_end = gaf_record.target_length - gaf_record.target_start;
    }
    for (int64_t i = paf_target_start; i < paf_target_end; ++i) {
        ++mz_counts[i].second;        
    }
    
    // increment each minimzer
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

