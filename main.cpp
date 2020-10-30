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
       << "    -g, --min-gap N                     Filter so that reported minimizer matches have >=N bases between them [500]" << endl;
}    

int main(int argc, char** argv) {

    string target_prefix;
    int64_t min_block_len = 100000;
    int64_t min_mapq = 5;
    int64_t min_gap = 500;
    
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
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hp:b:q:g:",
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

    scan_mzgaf(*in_stream, [&](MzGafRecord& gaf_record, GafRecord& parent_record) {
            // todo: buffer and parallelize?
            if (gaf_record.num_minimizers > 0 &&
                parent_record.mapq >= min_mapq &&
                parent_record.block_length >= min_block_len) {
                
                mzgaf2paf(gaf_record, parent_record, out_stream, min_gap, target_prefix);
            }
        });

    return 0;
}
