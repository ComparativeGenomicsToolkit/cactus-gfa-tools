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

void mzgaf2paf(const MzGafRecord& gaf_record, ostream& paf_stream, const string& target_prefix) {

    paf_stream << gaf_record.query_name << "\t"
               << gaf_record.query_length << "\t"
               << gaf_record.query_start << "\t"
               << gaf_record.query_end << "\t"
               << (gaf_record.is_reverse ? "-" : "+") << "\t"
               << target_prefix << gaf_record.target_name << "\t"
               << gaf_record.target_length << "\t"
               << gaf_record.target_start << "\t"
               << gaf_record.target_end << "\t";

    // do the cigar string
    assert(gaf_record.target_mz_offsets.size() == gaf_record.query_mz_offsets.size());

    stringstream cigar_stream;
    cigar_stream << "cg:Z:";

    // positions as we move along (relative to gaf_record.query/target_starts)
    int64_t query_pos = 0;
    int64_t target_pos = 0;
    // positions of our current match block (relative as above
    int64_t query_start = 0;
    int64_t query_end = gaf_record.kmer_size;
    int64_t prev_query_end = numeric_limits<int64_t>::max();
    int64_t target_start = 0;
    int64_t target_end = gaf_record.kmer_size;
    int64_t prev_target_end = numeric_limits<int64_t>::max();
    // need for output
    int64_t total_matches = 0;
    // used for sanity checks
    int64_t total_deletions = 0;
    int64_t total_insertions = 0;

    // every minimzer will either be a new cigar match record, or will extend the current one
    for (size_t i = 0; i < gaf_record.num_minimizers; ++i) {

#ifdef debug
        cerr << "[" << i << "]: query_pos = " << query_pos << " target_pos = " << target_pos << " (query_start = " << query_start
             << " query_end = " << query_end << ") (target_start = " << target_start << " target_end = " << target_end << ")" << endl;
#endif

        // compute the overlap with the previous minimizer
        int64_t query_delta = query_pos - query_end;
        int64_t target_delta = target_pos - target_end;

        if (query_delta == target_delta && query_delta <= 0) {
            // if the deltas are the same and both negative, we can just extend the current block
            query_end += gaf_record.kmer_size + query_delta;
            target_end += gaf_record.kmer_size + target_delta;
        }  else {
            int64_t min_delta = min((int64_t)0, min(query_delta, target_delta));
            if (min_delta < 0) {
                // there is an inconsistent overlap, which implies conflicting alignment.  we cut the blocks
                // to leave the conflicting bits unaligned
                query_end += min_delta;
                target_end += min_delta;
#ifdef debug
                cerr << "query_delta=" << query_delta << " target_delta=" << target_delta << " min_delta=" << min_delta << endl;
#endif
            }
                
            // we are going to make a new block, let's output the previous hit into the cigar
            assert(query_end - query_start == target_end - target_start);
            int64_t match_size = query_end - query_start;
            if (match_size > 0) {
                // catch up on query
                if (query_start > prev_query_end) {
                    cigar_stream << (query_start - prev_query_end) << "I";
                    total_insertions += (query_start - prev_query_end);
                }
                // catch up on target
                if (target_start > prev_target_end) {
                    cigar_stream << (target_start - prev_target_end) << "D";
                    total_deletions += (target_start - prev_target_end);
                }
                cigar_stream << match_size << "M";
                total_matches += match_size;
                // store end of block to compute indels needed before next
                prev_query_end = query_end;
                prev_target_end = target_end;
            }

#ifdef debug
            cerr << "  print previous block as " << (query_end - query_start) << "M" << endl;
#endif

            // start new block (cutting the overlap off the front with min_delta)
            query_start = query_pos - min_delta;
            query_end = query_pos + gaf_record.kmer_size;
            target_start = target_pos - min_delta;
            target_end = target_pos + gaf_record.kmer_size;
        }

        // advance our position
        if (i < gaf_record.num_minimizers - 1) {
            query_pos += gaf_record.query_mz_offsets[i];
            target_pos += gaf_record.target_mz_offsets[i];
        }
    }

    assert(gaf_record.query_start + query_end == gaf_record.query_end);
    assert(gaf_record.target_start + target_end == gaf_record.target_end);

    // output the last block    
    assert(query_end - query_start == target_end - target_start);
    int64_t match_size = query_end - query_start;
    if (match_size > 0) {
        // catch up on query
        if (query_start > prev_query_end) {
            cigar_stream << (query_start - prev_query_end) << "I";
            total_insertions += (query_start - prev_query_end);
        }
        // catch up on target
        if (target_start > prev_target_end) {
            cigar_stream << (target_start - prev_target_end) << "D";
            total_deletions += (target_start - prev_target_end);
        }
        cigar_stream << match_size << "M";
        total_matches += match_size;
    }

    assert(total_insertions + total_matches == gaf_record.query_end - gaf_record.query_start);
    assert(total_deletions + total_matches == gaf_record.target_end - gaf_record.target_start);

    // do the last 3 columns the cigar
    paf_stream << total_matches << "\t"
               << (gaf_record.target_end - gaf_record.target_start) << "\t" // fudged
               << 255 << "\t"
               << cigar_stream.str()
               << "\n";
}

