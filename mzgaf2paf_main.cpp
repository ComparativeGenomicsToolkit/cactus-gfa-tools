#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <fstream>

#include "mzgaf2paf.hpp"

using namespace std;
using namespace gafkluge;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <gaf> [gaf2] [gaf3] [...] > output.paf" << endl
       << "Convert minigraph --write-mz output(s) to PAF" << endl
       << endl
       << "options: " << endl
       << "    -p, --target-prefix PREFIX          Prepend all target (graph) contig names with this prefix" << endl
       << "    -b, --min-block-length N            Ignore records with block length (GAF col 11) (only applies if query length > N)< N [0]" << endl      
       << "    -q, --min-mapq N                    Ignore records with MAPQ (GAF col 12) < N [0]" << endl
       << "    -g, --min-gap N                     Filter so that reported minimizer matches have >=N bases between them [0]" << endl
       << "    -m, --min-match-len N               Only write matches (formed by overlapping/adjacent mz chains) with length < N" << endl
       << "    -u, --universal-mz FLOAT            Filter minimizers that appear in fewer than this fraction of alignments to target [0]" << endl
       << "    -n, --node-based-universal          Universal computed on entire node instead of mapped region" << endl
       << "    -s, --min-node-length N             Ignore minimizers on GAF nodes of length < N [0]" << endl
       << "    -i, --strict-unversal               Count mapq and block length filters against universal (instead of ignoring)" << endl
       << "    -o, --min-overlap-length N          If >= query regions with size >= N overlap, ignore the query region.  If 1 query region with size >= N overlaps any regions of size <= N, ignore the smaller ones only. (0 = disable) [0]" << endl;
}    

int main(int argc, char** argv) {

    string target_prefix;
    int64_t min_block_len = 0;
    int64_t min_mapq = 0;
    int64_t min_gap = 0;
    int64_t min_match_length = 0;
    double universal_filter = 0.;
    bool node_based_universal = false;
    // toggle on file-based filtering, which will catch and filter out cases where the same minimizer
    // is touched more than once in a file.
    bool file_based_filter = false;
    int64_t min_node_len = 0;
    bool strict_universal = false;
    int64_t min_overlap_len = 0;
       
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"target-prefix", required_argument, 0, 'p'},
            {"min-block-length", required_argument, 0, 'b'},
            {"min-mapq", required_argument, 0, 'q'},
            {"min-gap", required_argument, 0, 'g'},
            {"min-match-len", required_argument, 0, 'm'},
            {"universal-mz", required_argument, 0, 'u'},
            {"node-based-universal", no_argument, 0, 'n'},
            {"min-node-length", required_argument, 0, 's'},
            {"strict-unversal", no_argument, 0, 'i'},
            {"min-overlap-length", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hp:b:q:g:m:u:ns:io:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            target_prefix = optarg;
            break;
        case 'b':
            min_block_len = std::stol(optarg);
            break;
        case 'q':
            min_mapq = std::stol(optarg);
            break;
        case 'g':
            min_gap = std::stol(optarg);
            break;
        case 'm':
            min_match_length = std::stol(optarg);
            break;
        case 'u':
            universal_filter = std::stof(optarg);
            break;
        case 'n':
            node_based_universal = true;
            break;
        case 's':
            min_node_len = std::stol(optarg);
            break;
        case 'i':
            strict_universal = true;
            break;
        case 'o':
            min_overlap_len = std::stol(optarg);
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }

    // Parse the positional argument
    if (optind >= argc) {
        cerr << "[mzgaf2paf] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    vector<string> in_paths;
    int stdin_count = 0;
    while (optind < argc) {
        in_paths.push_back(argv[optind++]);
        if (in_paths.back() == "-") {
            ++stdin_count;
        }
    }

    if ((universal_filter > 0 || min_overlap_len > 0) && stdin_count > 0) {
        cerr << "[mzgaf2paf] error: nonzero values for -u and -o require 2 passes, so input cannot be streamed in with -" << endl;
        return 1;
    }

    if (stdin_count > 1) {
        cerr << "[mzgaf2paf] error: only one input can be piped with -" << endl;
        return 1;
    }

    // keep global counts of minimizers (used only for the universal filter)
    MZMap mz_map;

    // per query region counts (todo: can save some memory with an interval structure)
    QueryCoverage query_coverage;

    size_t total_match_length = 0;
    size_t total_target_block_length = 0;
    size_t total_records = 0;

    for (const string& in_path : in_paths) {

        ifstream in_file;
        istream* in_stream;
        if (in_path == "-") {
            in_stream = &cin;
        } else {
            in_file.open(in_path);
            if (!in_file) {
                cerr << "[mzgaf2paf] error: unable to open input: " << in_path << endl;
                return 1;
            }
            in_stream = &in_file;
        }

        // optional first pass counts (very inefficiently atm) how many queries each minimizer appears in
        if (universal_filter > 0 || min_overlap_len > 0) {
            // per file counts (todo: we can save some memory by using fewer bits here in some cases)
            MZMap file_mz_map;
            scan_mzgaf(*in_stream,
                       [&](MzGafRecord& gaf_record, GafRecord& parent_record) {
                           // todo: buffer and parallelize?
                           if (strict_universal ||
                               (gaf_record.num_minimizers > 0 &&
                                parent_record.mapq >= min_mapq &&
                                (parent_record.query_length <= min_block_len || parent_record.block_length >= min_block_len) &&
                                gaf_record.target_length >= min_node_len)) {                        
                               update_mz_map(gaf_record, parent_record, file_mz_map, min_mapq, min_block_len, min_node_len, node_based_universal);
                           }
                       },
                       [&](GafRecord& parent_record) {
                           if (min_overlap_len > 0 && parent_record.block_length >= min_overlap_len) {
                               update_query_coverage(parent_record, query_coverage);
                           }
                       });
            
            // go back to the beginning by resetting the stream
            in_stream->clear();
            in_stream->seekg(0, ios::beg);

            // combine the file-map into the overall map.  if universal_filter is 1, then
            // minimizers counted more than once in the file get zapped to 0
            combine_mz_maps(file_mz_map, mz_map, file_based_filter);
        }
    
        scan_mzgaf(*in_stream, [&](MzGafRecord& gaf_record, GafRecord& parent_record) {
                // todo: buffer and parallelize?
                if (gaf_record.num_minimizers > 0 &&
                    parent_record.mapq >= min_mapq &&
                    (parent_record.query_length <= min_block_len || parent_record.block_length >= min_block_len) &&
                    gaf_record.target_length >= min_node_len) {

                    total_match_length += mzgaf2paf(gaf_record,
                                                    parent_record,
                                                    cout,
                                                    min_gap,
                                                    min_match_length,
                                                    mz_map,
                                                    universal_filter,
                                                    query_coverage,
                                                    min_overlap_len,
                                                    target_prefix);
                    
                    total_target_block_length += gaf_record.target_end - gaf_record.target_start;
                    ++total_records;
                }
            });

    }
    
    cerr << "Converted " << total_records << " recs with " << total_match_length << " bp of cigar Matches over " << total_target_block_length
         << " bp of alignments to target (" << ((double)(total_match_length) / total_target_block_length) <<")" << endl;
    
    return 0;
}
