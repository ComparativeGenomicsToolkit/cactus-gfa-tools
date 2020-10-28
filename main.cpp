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
       << "    -p, --target-prefix PREFIX          Prepend all target (graph) contig names with this prefix" << endl;
}    

int main(int argc, char** argv) {

    string target_prefix;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"threads", required_argument, 0, 't'},
            {"help", no_argument, 0, 'h'},
            {"target-prefix", required_argument, 0, 'p'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hp:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            target_prefix = optarg;
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

    scan_mzgaf(*in_stream, [&](MzGafRecord& gaf_record) {
            // todo: buffer and parallelize?
            if (gaf_record.num_minimizers > 0) {
                mzgaf2paf(gaf_record, out_stream, target_prefix);
            }
        });

    return 0;
}
