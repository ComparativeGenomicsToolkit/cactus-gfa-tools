
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>

#include "mzgaf2paf.hpp"
#include "pafcoverage.hpp"

//#define debug

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <rgfa> <header_table> <paf>" << endl
       << "Convert PAF from minigraph to stable coordinates using mapping in rgfa" << endl
       << endl
       << "options: " << endl
       << "    -o, --only-query         Only convert targets that appear as a query " << endl
       << endl;
}    

int main(int argc, char** argv) {

    bool only_query = false;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"only-query", no_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "ho",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'o':
            only_query = true;
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

    // Parse the positional arguments
    if (optind >= argc + 2) {
        cerr << "[mask] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    string in_rgfa_path = argv[optind++];
    string in_table_path = argv[optind++];
    string in_paf_path = argv[optind++];

    if (optind < argc - 1) {
        cerr << "[pafmask] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    // Buid the lookup table
    unordered_map<string, tuple<string, int64_t, int64_t>> stable_lookup = build_stable_lookup(in_table_path, in_rgfa_path);

    // open the paf
    ifstream in_paf_file(in_paf_path);
    if (!in_paf_file) {
        cerr << "[paf2stable] error: unable to open paf: " << in_paf_path << endl;
        return 1;
    }

    // load up the query contigs
    unordered_set<string> query_set;
    string buffer;
    while (getline(in_paf_file, buffer)) {
        // split into array of tokens
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);

        if (toks.size() < 12) {
            throw runtime_error("too few tokens in PAF line: " + buffer);
        }

        string query_name = toks[0];
        query_set.insert(query_name);
    }
        
    in_paf_file.close();
    in_paf_file.open(in_paf_path);

    // convert the PAF
    while (getline(in_paf_file, buffer)) {
        // split into array of tokens
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);
        string target_name = toks[5];
        int64_t target_start = stol(toks[7]);
        int64_t target_end = stol(toks[8]);
        // cut off cactus_prefix
        if (target_name.compare(0, 3, "id=") == 0) {
            size_t bar_pos = target_name.find("|");
            assert(bar_pos != string::npos && bar_pos < target_name.length() - 1);
            target_name = target_name.substr(bar_pos + 1);
        }
        
        tuple<string, int64_t, int64_t>& stable_info = stable_lookup.at(target_name);
        toks[5] = get<0>(stable_info);
        toks[6] = std::to_string(get<2>(stable_info));
        toks[7] = std::to_string(target_start + get<1>(stable_info));
        toks[8] = std::to_string(target_end + get<1>(stable_info));
        
        if (only_query || query_set.count(toks[5])) {
            cout << toks[0];
            for (size_t i = 1; i < toks.size(); ++i) {
                cout << "\t" << toks[i];
            }
        } else {
            cout << buffer;
        }
        cout << endl;
    }
    
    return 0;
}
