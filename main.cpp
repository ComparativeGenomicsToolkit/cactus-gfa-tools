#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <fstream>

#include "mzgaf2paf.hpp"

using namespace std;
using namespace gafkluge;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <gaf> > output.paf" << endl
       << "Convert minigraph --write-mz output to PAF" << endl
       << endl
       << "options: " << endl
       << "    -p, --target-prefix PREFIX          Prepend all target (graph) contig names with this prefix" << endl
       << "    -b, --min-block-length N            Ignore records with block length (GAF col 11) < N [100000]" << endl      
       << "    -q, --min-mapq N                    Ignore records with MAPQ (GAF col 12) < N [5]" << endl
       << "    -g, --min-gap N                     Filter so that reported minimizer matches have >=N bases between them [0]" << endl
       << "    -m, --min-match-len N               Only write matches (formed by overlapping/adjacent mz chains) with length < N" << endl
       << "    -u, --universal-mz FLOAT            Filter minimizers that appear in fewer than this fraction of alignments to target [0]" << endl;
}    

int main(int argc, char** argv) {

    string target_prefix;
    int64_t min_block_len = 100000;
    int64_t min_mapq = 5;
    int64_t min_gap = 0;
    int64_t min_match_length = 0;
    double universal_filter = 0.;
    
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
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hp:b:q:g:m:u:",
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

    string in_path = argv[optind++];    

    if (universal_filter > 0 && in_path == "-") {
        cerr << "[mzgaf2paf] error: -u option requires 2 passes, so input cannot be streamed in with -" << endl;
        return 1;
    }

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

    ostream& out_stream = cout;

    // optional first pass counts (very inefficiently atm) how many queries each minimizer appears in
    MZMap mz_map;
    if (universal_filter > 0) {
        scan_mzgaf(*in_stream, [&](MzGafRecord& gaf_record, GafRecord& parent_record) {
                // todo: buffer and parallelize?
                if (gaf_record.num_minimizers > 0 &&
                    parent_record.mapq >= min_mapq &&
                    parent_record.block_length >= min_block_len) {

                    update_mz_map(gaf_record, parent_record, mz_map);
                }
            });

        // go back to the beginning by resetting the stream
        in_stream->clear();
        in_stream->seekg(0, ios::beg);
    }

    size_t total_match_length = 0;
    size_t total_target_block_length = 0;
    
    scan_mzgaf(*in_stream, [&](MzGafRecord& gaf_record, GafRecord& parent_record) {
            // todo: buffer and parallelize?
            if (gaf_record.num_minimizers > 0 &&
                parent_record.mapq >= min_mapq &&
                parent_record.block_length >= min_block_len) {

                total_match_length += mzgaf2paf(gaf_record, parent_record, out_stream, min_gap, min_match_length, mz_map, universal_filter, target_prefix);
                total_target_block_length += gaf_record.target_end - gaf_record.target_start;
            }
        });


    cerr << "Converted " << total_match_length << " bp of cigar Matches over " << total_target_block_length
         << " bp of alignments to target (" << ((double)(total_match_length) / total_target_block_length) <<")" << endl;
    
    return 0;
}
