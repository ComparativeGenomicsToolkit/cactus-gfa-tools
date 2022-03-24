#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <algorithm>
#include <sys/stat.h>
#include "gfakluge.hpp"
#include "pafcoverage.hpp"

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <GFA>" << endl
       << "Create PAF from rGFA tags, representing the alignment of reference contig intervals to GFA nodes.  Input cannot be stdin." << endl
       << "options: " << endl
       << "    -r, --max-rank N                    Process nodes with rank <= N [0]" << endl
       << "    -q, --query-lengths FILE            Tab-separated file listing query contig lengths" << endl
       << "    -T, --target-prefix STRING          Prefix all target (reference) contig names with STRING" << endl
       << "    -P, --query-prefix STRING           Prefix all query contig names with STRING" << endl
       << "    -i, --ignore-queries FILE           Ignore all sequences found in paf FILE" << endl
       << endl;
}    

int main(int argc, char** argv) {

    string rgfa_path;
    int64_t max_rank = 0;
    string query_lengths_path;
    string query_prefix;
    string target_prefix;
    string ignore_paf_path;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"max-rank", required_argument, 0, 'r'},
            {"query-lengths", required_argument, 0, 'q'},
            {"target-prefix", required_argument, 0, 'T'},
            {"query-prefix", required_argument, 0, 'P'},
            {"ignore-queries", required_argument, 0, 'i'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hr:q:T:P:i:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            max_rank = stol(optarg);
            break;
        case 'q':
            query_lengths_path = optarg;
            break;
        case 'T':
            target_prefix = optarg;
            break;
        case 'P':
            query_prefix = optarg;
            break;
        case 'i':
            ignore_paf_path = optarg;
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

    if (optind != argc - 1) {
        cerr << "[rgfa2paf] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    rgfa_path = argv[optind++];

    if (max_rank > 0 && query_lengths_path.empty()) {
        cerr << "[rgfa2paf] error: Query lengths (-q) must be provided when max rank > 0" << endl;
        return 1;
    }
        
    // check input
    {
        ifstream gfa_stream(rgfa_path);
        if (!gfa_stream) {
            cerr << "[rgfa2paf] error: Unable to read file: " << rgfa_path << endl;
            return 1;
        }
    }

    unordered_map<string, int64_t> query_lengths;
    
    // read query lengths from file
    if (!query_lengths_path.empty()) {
        ifstream query_lengths_stream(query_lengths_path);
        if (!query_lengths_stream) {
            cerr << "[rgfa2paf] error: Unable to read query lengths file: " << query_lengths_path << endl;
            return 1;
        }
        string buffer;
        while (getline(query_lengths_stream, buffer)) {
            vector<string> toks;
            split_delims(buffer, "\t\n", toks);
            if (toks.size() > 1) {
                string contig = toks[0];
                if (contig.compare(0, query_prefix.length(), query_prefix) != 0) {
                    contig = query_prefix + contig;
                }                
                query_lengths[contig] = stol(toks[1]);
            }
        }
    }

    // read the queries we want to ignore
    unordered_set<string> ignore_set;
    if (!ignore_paf_path.empty()) {
        ifstream ignore_paf_stream(ignore_paf_path);
        if (!ignore_paf_stream) {
            cerr << "[rgfa2paf] error: Unable to read query ignore PAF file: " << ignore_paf_path << endl;
            return 1;
        }
        string buffer;
        while (getline(ignore_paf_stream, buffer)) {
            vector<string> toks;
            split_delims(buffer, "\t\n", toks);
            if (toks.size() > 1) {
                string contig = toks[0];
                if (contig.compare(0, query_prefix.length(), query_prefix) != 0) {
                    contig = query_prefix + contig;
                }                                
                ignore_set.insert(contig);
            }
        }
    }

    // gfa pass 1: get total query lengths
    function<void(const gfak::sequence_elem&)> visit_seq_1 = [&](const gfak::sequence_elem& gfa_seq) {
        bool found_SN = false;
        bool found_SR = false;
        bool found_SO = false;
        int64_t rank;
        string contig;
        int64_t offset;
        // parse the SN (contig), SR (rank) and SO (offset) optional tags
        for (const gfak::opt_elem& oe : gfa_seq.opt_fields) {
            if (oe.key == "SN") {
                assert(found_SN == false);
                contig = oe.val;
                found_SN = true;
            } else if (oe.key == "SR") {
                assert(found_SR == false);
                rank = stol(oe.val);
                assert(rank >= 0);                               
                found_SR = true;
            } else if (oe.key == "SO") {
                assert(found_SO == false);
                offset = stol(oe.val);
                assert(offset >= 0);
                found_SO = true;
            }
        }
        assert(found_SN);
        assert(found_SR);
        assert(found_SO);
        if (contig.compare(0, query_prefix.length(), query_prefix) != 0) {
            contig = query_prefix + contig;
        }

        if (rank <= max_rank) {
            query_lengths[contig] += gfa_seq.sequence.length();
        }
    };

    // pass 2: convert rank 0 node to paf line
    function<void(const gfak::sequence_elem&)> visit_seq_2 = [&](const gfak::sequence_elem& gfa_seq) {
        bool found_SN = false;
        bool found_SR = false;
        bool found_SO = false;
        int64_t rank;
        string contig;
        int64_t offset;
        // parse the SN (contig), SR (rank) and SO (offset) optional tags
        for (const gfak::opt_elem& oe : gfa_seq.opt_fields) {
            if (oe.key == "SN") {
                assert(found_SN == false);
                contig = oe.val;
                found_SN = true;
            } else if (oe.key == "SR") {
                assert(found_SR == false);
                rank = stol(oe.val);
                assert(rank >= 0);                               
                found_SR = true;
            } else if (oe.key == "SO") {
                assert(found_SO == false);
                offset = stol(oe.val);
                assert(offset >= 0);
                found_SO = true;
            }
        }
        assert(found_SN);
        assert(found_SR);
        assert(found_SO);
        if (contig.compare(0, query_prefix.length(), query_prefix) != 0) {
            contig = query_prefix + contig;
        }

        string target = gfa_seq.name;
        if (target.compare(0, target_prefix.length(), target_prefix) != 0) {
            target = target_prefix + target;
        }

        if (rank <= max_rank && !ignore_set.count(contig)) {
            // emit the paf line
            cout << contig << "\t"
                 << query_lengths[contig] << "\t"
                 << offset << "\t"
                 << offset + gfa_seq.sequence.length() << "\t"
                 << "+" << "\t"
                 << target << "\t"
                 << gfa_seq.sequence.length() << "\t"
                 << "0" << "\t"
                 << gfa_seq.sequence.length()  << "\t"
                 << gfa_seq.sequence.length() << "\t"
                 << gfa_seq.sequence.length() << "\t"
                 << "60" << "\t"
                 << "cg:Z:" << gfa_seq.sequence.length() << "M"
                 << "\n";                
        }
    };
    
    gfak::GFAKluge kluge;
    if (query_lengths_path.empty()) {
        kluge.for_each_sequence_line_in_file(rgfa_path.c_str(), visit_seq_1);
    }
    kluge.for_each_sequence_line_in_file(rgfa_path.c_str(), visit_seq_2);

    return 0;
}
