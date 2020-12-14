#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>

#include "rgfa2contig.hpp"

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <rgfa>" << endl
       << "Partition rGFA nodes into reference contigs" << endl
       << endl;
}    

int main(int argc, char** argv) {

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

        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
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

    if (optind != argc - 1) {
        cerr << "[rfa2contig] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    // get stream to gfa (not supporting compression)
    string in_path = argv[optind];

    ifstream in_file;
    istream* in_stream;
    if (in_path == "-") {
        in_stream = &cin;
        // thanks tinygfa!
        cerr << "[rgfa2contig] error: - (stdin) not supported" << endl;
        return 1;
    } else {
        in_file.open(in_path);
        if (!in_file) {
            cerr << "[pafcoverage] error: unable to open input: " << in_path << endl;
            return 1;
        }
        in_stream = &in_file;
    }

    // get the parittion of GFA nodes -> reference contig
    pair<unordered_map<int64_t, int64_t>, vector<string>> partition = rgfa2contig(in_file);

    // print the partition as a tsv
    for (auto& node_contig : partition.first) {
        cout << "S" << node_contig.first << "\t" << partition.second[node_contig.second] << "\n";
    }

    return 0;
}
