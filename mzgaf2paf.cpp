/**
 * mzgaf2paf.hpp: Make base level pairwise alignemnts from minigraph --write-mz output with the object of using them as
 *                anchors for other graph methods
 */

#include "mzgaf2paf.hpp"
#include "mzgaf.hpp"

#define debug

using namespace gafkluge;
using namespace std;

void mzgaf2paf(const MzGafRecord& gaf_record, ostream& paf_stream, const string& target_prefix) {

    paf_stream << gaf_record.query_name << "\t"
               << gaf_record.query_length << "\t"
               << gaf_record.query_start << "\t"
               << gaf_record.query_end << "\t"
               << (gaf_record.is_reverse ? "+" : "-") << "\t"
               << target_prefix << gaf_record.target_name << "\t"
               << gaf_record.target_length << "\t"
               << gaf_record.target_start << "\t"
               << gaf_record.target_end << "\t";

    // do the cigar string
    assert(gaf_record.target_mz_offsets.size() == gaf_record.query_mz_offsets.size());

    stringstream cigar_stream;
    cigar_stream << "cg:Z:";

    // queue up our first minimizer block, which starts at 0 by definition (all relative to gaf_record coordinates)
    int64_t query_pos = 0;
    int64_t query_start = 0;
    int64_t query_end = gaf_record.kmer_size;
    int64_t target_pos = 0;
    int64_t target_start = 0;
    int64_t target_end = gaf_record.kmer_size;

    int64_t total_matches = 0;
    
    for (size_t i = 0; i < gaf_record.query_mz_offsets.size(); ++i) {
        // the offset from the gaf, is relative to the first position of the previous block
        query_pos += gaf_record.query_mz_offsets[i];
        target_pos += gaf_record.target_mz_offsets[i];

#ifdef debug
        cerr << "[" << i << "]: query_pos = " << query_pos << " target_pos = " << target_pos << " (query_start = " << query_start
             << " query_end = " << query_end << ") (target_start = " << target_start << " target_end = " << target_end << ")" << endl;
#endif

        // compute the overlap with the previous minimizer
        int64_t query_delta = query_pos - query_end;
        int64_t target_delta = target_pos - target_end;

        if (query_delta <= 0 || target_delta <= 0) {
            // there is an overlap
            int64_t min_delta = std::max(query_delta, target_delta);
            // and we can extend in both query and target
            if (min_delta <= 0) {
                query_end += gaf_record.kmer_size + min_delta;
                target_end += gaf_record.kmer_size + min_delta;
#ifdef debug
                cerr << "  overlap found, extending by delta " << (gaf_record.kmer_size + min_delta - 1) << " gives query_end = " << query_end << " target_end = " << target_end << endl;
#endif
            } else {
                // i'm not sure if this can happen?
                cerr << "Warning [mzgaf2paf] : inconsistent overlap at position " << i << " of line with query " << gaf_record.query_name << ":" << gaf_record.query_start << endl;
                // figure this out if it does
                exit(1);
            }
        } else {
            // there is no overlap, let's output the previous hit into the cigar
            assert(query_end - query_start == target_end - target_start);
            int64_t match_size = query_end - query_start;
            cigar_stream << match_size << "M";
            total_matches += match_size;

#ifdef debug
            cerr << "  print previous block as " << (query_end - query_start) << "M" << endl;
#endif

            // output the deltas as indels
            if (query_delta > 0) {
                cigar_stream << "I" << query_delta;
            }
            if (target_delta > 0) {
                cigar_stream << "D" << target_delta;
            }

            // start new block
            query_start = query_pos;
            query_end = query_start + gaf_record.kmer_size;
            target_start = target_pos;
            target_end = target_start + gaf_record.kmer_size;
#ifdef debug
            cerr << "  starting new block query_start = " << query_start << " query_end = " << query_end << endl;
#endif
        }
    }

    #ifdef debug
        cerr << "[final]: (query_start = " << query_start
             << " query_end = " << query_end << ") (target_start = " << target_start << " target_end = " << target_end << ")" << endl;
#endif


    // output the last block
    assert(gaf_record.query_start + query_end == gaf_record.query_end);
    assert(gaf_record.target_start + target_end == gaf_record.target_end);

    assert(query_end - query_start == target_end - target_start);
    int64_t match_size = query_end - query_start;
    cigar_stream << match_size << "M";
    total_matches += match_size;

    // do the last 3 columns the cigar
    paf_stream << total_matches << "\t"
               << (gaf_record.target_end - gaf_record.target_start) << "\t" // fudged
               << 255 << "\t"
               << cigar_stream.str()
               << "\n";
}

