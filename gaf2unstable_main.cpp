#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <set>
#include "paf2stable.hpp"
#include "pafcoverage.hpp"
#include "gfakluge.hpp"
#include "gafkluge.hpp"

//#define debug

using namespace std;
using namespace gafkluge;

struct MGSeq {
   string name;
   int64_t offset;
   int64_t length;
};

static inline bool operator<(const MGSeq& s1, const MGSeq& s2) {
    return s1.offset < s2.offset;
}
static inline ostream& operator<<(ostream& os, const MGSeq& s1) {
    os << "mg{name=" << s1.name << ",offset=" << s1.offset << ",len=" << s1.length << "}";
    return os;
}

// maps a stable contig name to a set of minigraph sequences that are sorted by their
// offset in the stable sequence.
static unordered_map<string, set<MGSeq>> get_unstable_mapping(const string& rgfa_path) {
    unordered_map<string,set< MGSeq>> unstable_mapping;
    if (!ifstream(rgfa_path)) {
        cerr << "[gaf2unstable] error: Could not open " << rgfa_path << endl;
        exit(1);
    }
        
    gfak::GFAKluge kluge;    
    kluge.for_each_sequence_line_in_file(rgfa_path.c_str(), [&](const gfak::sequence_elem& gfa_seq) {
            bool found_SN = false;
            bool found_SO = false;
            MGSeq mg_seq;
            mg_seq.length = gfa_seq.sequence.length();
            mg_seq.name = gfa_seq.name;
            string contig;
            for (const gfak::opt_elem& oe : gfa_seq.opt_fields) {
                if (oe.key == "SN") {
                    assert(found_SN == false);
                    contig = oe.val;
                    found_SN = true;
                } else if (oe.key == "SO") {
                    assert(found_SO == false);
                    mg_seq.offset = stol(oe.val);
                    assert(mg_seq.offset >= 0);
                    found_SO = true;
                }
            }
            assert(found_SN);
            assert(found_SO);

            unstable_mapping[contig].insert(mg_seq);

        });
    return unstable_mapping;
}

static vector<pair<MGSeq, pair<int64_t, int64_t>>>  get_unstable_interval(const unordered_map<string, set<MGSeq>>& lookup,
                                                                          const string& stable_contig, int64_t start, int64_t end) {
    assert(lookup.count(stable_contig));
    auto& contig_idx = lookup.at(stable_contig);
    // find the node that overlaps start
    MGSeq query;
    query.offset = start;
    auto interval_start = contig_idx.upper_bound(query);
    assert(interval_start != contig_idx.begin());
    --interval_start;
    // find the node that overlaps end
    query.offset = end;
    auto interval_end = contig_idx.lower_bound(query);
    assert(interval_end != contig_idx.begin());

    vector<pair<MGSeq, pair<int64_t, int64_t>>> unstable_intervals;

    int64_t ui_len = 0;
    for (auto i = interval_start; i != interval_end; ++i) {
        unstable_intervals.push_back(make_pair(*i, make_pair(0, i->length)));
        ui_len += i->length;
    }

    // clip the endpoints
    if (unstable_intervals[0].first.offset > start) {
        unstable_intervals[0].second.first = unstable_intervals[0].first.offset - start;
        ui_len -= unstable_intervals[0].second.first;
    }
    
    if (ui_len > end - start) {
        unstable_intervals.back().second.second -= ui_len - (end -start);
        ui_len -= ui_len - (end -start);
        assert(unstable_intervals.back().second.second > 0);
    }
    if (end != numeric_limits<int64_t>::max()) {
        assert(ui_len == end - start);
    }
    return unstable_intervals;
}

static void gaf2unstable(const unordered_map<string, set<MGSeq>>& lookup, GafRecord& gaf_record) {
    vector<GafStep> unstable_path;
    for (auto& step : gaf_record.path) {        
        auto unstable_interval = get_unstable_interval(lookup, step.name, step.start,
                                                       step.is_interval ? step.end : numeric_limits<int64_t>::max());
        if (step.is_reverse) {
            std::reverse(unstable_interval.begin(), unstable_interval.end());            
        }

        for (const auto& unstable_frag : unstable_interval) {
            GafStep unstable_step;
            unstable_step.name = unstable_frag.first.name;
            unstable_step.is_reverse = step.is_reverse;
            unstable_step.is_stable = false;
            unstable_step.start = unstable_frag.second.first;
            unstable_step.end = unstable_frag.second.second;
            unstable_step.is_interval = unstable_step.start != 0 || unstable_step.end != unstable_frag.first.length;
            unstable_path.push_back(unstable_step);
        }

#ifdef debug
        cerr << "Step " << step << " maps to ";
        for (const auto& s : unstable_interval) {
            cerr << s.first << " " << s.second.first << "," << s.second.second << "; ";
        }
        cerr << endl;        
#endif        
    }
    gaf_record.path = unstable_path;
}

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <gaf> " << endl
       << "Replace stable sequences in path steps, ex >chr1:500-1000, with the unstable graph node names, ex >s1:1-100>s2:100-600" << endl
       << endl
       << "options: " << endl
       << "    -g, --rGFA FILE           (uncompressed) minigraph rGFA, required to look up unstable mappings" << endl
       << "    -o, --out-lengths FILE    Output lengths of all minigraph sequences in given file (can be passed to gaf2paf)" << endl;
}    

int main(int argc, char** argv) {

    string rgfa_path;
    string node_lengths_path;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"rgfa", required_argument, 0, 'g'},
            {"out-lengths", required_argument, 0, '0'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hg:o:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'h':
        case 'g':
            rgfa_path = optarg;
            break;
        case 'o':
            node_lengths_path = optarg;
            break;
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
        cerr << "[gaf2unstable] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    string in_gaf_path = argv[optind++];

    if (optind < argc - 1) {
        cerr << "[gaf2unstable] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    if (rgfa_path.empty()) {
        cerr << "[gaf2unstable] error: -g option required" << endl;
        return 1;
    }

    // open the gaf file
    ifstream in_gaf_file;
    istream* in_gaf_stream;
    if (in_gaf_path == "-") {
        in_gaf_stream = &cin;
    } else {
        in_gaf_file.open(in_gaf_path);
        if (!in_gaf_file) {
            cerr << "[gaf2unstable] error: unable to open input: " << in_gaf_path << endl;
            return 1;
        }
        in_gaf_stream = &in_gaf_file;
    }

    // load the gfa
    auto lookup = get_unstable_mapping(rgfa_path);

    if (!node_lengths_path.empty()) {
        ofstream node_lengths_file(node_lengths_path);
        if (!node_lengths_file) {
            cerr << "[gaf2unstable] error: unable to open output: " << node_lengths_path << endl;
            return 1;
        }
        for (const auto& cs : lookup) {
            for (const auto& s : cs.second) {
                node_lengths_file << s.name << "\t" << s.length << "\n";
            }
        }
    }

    // process the gaf
    GafRecord gaf_record;
    string line_buffer;
    while (getline(*in_gaf_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            // skip -S stuff
            continue;
        }
        parse_gaf_record(line_buffer, gaf_record);
        gaf2unstable(lookup, gaf_record);
        cout << gaf_record << "\n";
    }

    return 0;
}
