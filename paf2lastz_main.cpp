#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "paf2lastz.hpp"

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <paf> [paf2] [paf3] [...] > output.cigar" << endl
       << "Convert PAF(s) with cg cigars to LASTZ cigars" << endl
       << endl
       << "options: " << endl
       << "    -q, --mapq-score          Take score from MAPQ field (PAF column 12) instead of AS tag" << endl
       << "    -s, --secondary-file      Separate out secondaries (tp tag == S) and write them to given path" << endl;
}    

int main(int argc, char** argv) {

    bool mapq_score = false;
    string secondary_path;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"mapq-score", no_argument, 0, 'q'},
            {"secondary-file", required_argument, 0, 's'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hqs:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'q':
            mapq_score = true;
            break;
        case 's':
            secondary_path = optarg;
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
        cerr << "[paf2lastz] error: too few arguments" << endl;
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
        cerr << "[paf2lastz] error: only one input can be piped with -" << endl;
        return 1;
    }

    ofstream secondary_file;
    if (!secondary_path.empty()) {
        secondary_file.open(secondary_path);
        if (!secondary_file) {
            cerr << "[paf2lastz] error: could not open secondary-file: " << secondary_path << endl;
            return 1;
        }
    }

    for (const string& in_path : in_paths) {

        ifstream in_file;
        istream* in_stream;
        if (in_path == "-") {
            in_stream = &cin;
        } else {
            in_file.open(in_path);
            if (!in_file) {
                cerr << "[paf2lastz] error: unable to open input: " << in_path << endl;
                return 1;
            }
            in_stream = &in_file;
        }

        string buffer;
        while (getline(*in_stream, buffer)) {
            pair<string, bool> lastz_record = paf2lastz(buffer, mapq_score);
            if (lastz_record.second == true && !secondary_path.empty()) {
                secondary_file << lastz_record.first << "\n";
            } else {
                cout << lastz_record.first << "\n";
            }
        }
    }
    
    return 0;
}
