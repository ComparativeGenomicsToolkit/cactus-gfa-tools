#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>

#include "pafcoverage.hpp"

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <paf> [paf2] [paf3] [...]" << endl
       << "Print some PAF coverages statistics for query sequences" << endl
       << endl
       << "options: " << endl
       << "    -p, --query-prefix PREFIX           Only look at query sequences with given prefix" << endl
       << "    -g, --print-gaps                    Print gaps in coverage in BED format" << endl;
}    

int main(int argc, char** argv) {

    string query_prefix;
    bool print_gaps = false;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"query-prefix", required_argument, 0, 'q'},
            {"print-gaps", no_argument, 0, 'g'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hp:g",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            query_prefix = optarg;
            break;
        case 'g':
            print_gaps = true;
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
        cerr << "[pafcoverage] error: too few arguments" << endl;
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

    if (stdin_count > 1) {
        cerr << "mzgaf2paf] error: only one input can be piped with -" << endl;
        return 1;
    }

    // coverage stats go here
    CoverageMap coverage_map;

    for (const string& in_path : in_paths) {

        ifstream in_file;
        istream* in_stream;
        if (in_path == "-") {
            in_stream = &cin;
        } else {
            in_file.open(in_path);
            if (!in_file) {
                cerr << "[pafcoverage] error: unable to open input: " << in_path << endl;
                return 1;
            }
            in_stream = &in_file;
        }

        string buffer;
        while (getline(*in_stream, buffer)) {
            if (buffer.substr(0, query_prefix.length()) == query_prefix) {
                update_coverage_map(buffer, coverage_map);
            }
        }
    }

    // print the bed
    if (print_gaps) {
        print_coverage_gaps_as_bed(coverage_map, cout);
    } else {
        print_coverage_summary(coverage_map, cout);
    }
        

    return 0;
}
