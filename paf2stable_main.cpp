
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include "paf2stable.hpp"
#include "pafcoverage.hpp"

//#define debug

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <paf>" << endl
       << "Replace every target sequence with a query sequence (preserving all transitive mappings between queries)" << endl
       << endl;
}    

int main(int argc, char** argv) {

    int64_t min_length = 1;
    int64_t padding = 100;
    bool validate = false;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
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

    // Parse the positional arguments
    if (optind >= argc ) {
        cerr << "[paf2stable] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    string in_paf_path = argv[optind++];

    if (optind < argc - 1) {
        cerr << "[paf2stable] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    ifstream paf_file(in_paf_path);
    if (!paf_file) {
        cerr << "[paf2stable] error: Unable to open input PAF file, \"" << in_paf_path << "\"" << endl;
        return 1;
    }

    // first pass: build the mapping table from target to query intervals
    unordered_map<string, int64_t> query_name_to_id;
    vector<pair<string, int64_t>> query_id_to_info;
    unordered_map<string, pair<int64_t, vector<StableInterval>>> target_to_intervals;

    string buffer;
    while (getline(paf_file, buffer)) {
        // split into array of tokens
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);        

        if (toks.size() < 12) {
            throw runtime_error("too few tokens in PAF line: " + buffer);
        }

        update_stable_mapping_info(toks, query_name_to_id, query_id_to_info, target_to_intervals);
    }

    // make the interval tree from the intervals
    unordered_map<string, StableIntervalTree> target_to_interval_tree = create_interval_trees(target_to_intervals);

    // second pass: output the paf with all targets replaced by queries
    paf_file.close();
    paf_file.open(in_paf_path);
    assert(paf_file);
    while (getline(paf_file, buffer)) {
        // split into array of tokens
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);        

        paf_to_stable(toks, query_id_to_info, target_to_interval_tree);
    }
    
    return 0;
}
