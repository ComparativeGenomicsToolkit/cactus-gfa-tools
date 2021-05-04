#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>

#include "pafcoverage.hpp"
#include "rgfa-split.hpp"

#define debug

using namespace std;

// make a map of sequence name -> interval tree from a bed file
// intervals get merged up if they overlap or are within padding of each other
static unordered_map<string, CoverageIntervalTree> load_bed(istream& bed_stream, int64_t padding);

// cut interval_b out of interval_a, any peices of a that remain get added to out_fragments
// the number of such pieces is 0 (a in b), 1 (b overlaps one end of a) or 2 (b in a)
static void interval_subtract(CoverageInterval& interval_a, CoverageInterval& interval_b, vector<CoverageInterval>& out_fragments);

// subtract all masked intervals from a paf line and output what's left
static void mask_paf_line(const string& paf_line, int64_t min_length, const unordered_map<string, CoverageIntervalTree>& ref_to_intervals);

// output paf line(s) corresponding to given sub-interval of the original paf line
static void clip_paf(const vector<string>& toks, const string& query_name, int64_t query_length, int64_t query_start, int64_t query_end,
                     CoverageInterval& interval, int64_t min_length);


void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <paf> <bed>" << endl
       << "Cut masked regions out of a paf file" << endl
       << endl
       << "options: " << endl
       << "    -m, --min-length N           Remove any remaining intervals less than N bp" << endl
       << "    -p, --padding N              Merge up bed intervals close than this [100]" << endl;
}    

int main(int argc, char** argv) {

    int64_t min_length = 1;
    int64_t padding = 100;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"min-length", required_argument, 0, 'm'},
            {"padding", required_argument, 0, 'p'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hm:p:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'm':
            min_length = stol(optarg);
            break;
        case 'p':
            padding = stol(optarg);
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
    if (optind >= argc + 1) {
        cerr << "[mask] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    string in_paf_path = argv[optind++];
    string in_bed_path = argv[optind++];

    if (optind < argc - 1) {
        cerr << "[pafmask] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    // open the paf
    istream* in_paf;
    ifstream in_paf_file;
    if (in_paf_path == "-") {
        in_paf = &cin;
    } else {
        in_paf_file.open(in_paf_path);
        if (!in_paf_file) {
            cerr << "[pafmask] error: unable to open paf: " << in_paf_path << endl;
            return 1;
        }
        in_paf = &in_paf_file;
    }

    // load the bed
    ifstream in_bed_file(in_bed_path);
    if (!in_bed_file) {
        cerr << "[pafmask] error: unable to open bed: " << in_bed_path << endl;
        return 1;
    }
    unordered_map<string, CoverageIntervalTree> ref_to_intervals = load_bed(in_bed_file, padding);

    string buffer;
    while (getline(*in_paf, buffer)) {
        mask_paf_line(buffer, min_length, ref_to_intervals);
    }

    return 0;
}

unordered_map<string, CoverageIntervalTree> load_bed(istream& bed_stream, int64_t padding) {
    // load bed
    unordered_map<string, vector<CoverageInterval>> intervals;
    string buffer;
    while (getline(bed_stream, buffer)) {
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);
        if (toks.size() >= 3) {
            string& name = toks[0];
            int64_t start = stol(toks[1]);
            int64_t end = stol(toks[2]);
            intervals[name].emplace_back(start, end, make_pair(0, 0));
        }
    }

    // make the trees merging up anything that overlaps within padding
    unordered_map<string, CoverageIntervalTree> trees;

    for (auto& ref_intervals : intervals) {
        CoverageIntervalTree interval_tree(ref_intervals.second);
        vector<CoverageInterval> merged_intervals;
        scan_coverage_intervals(interval_tree, padding, [&](int64_t start, int64_t stop, pair<int64_t, int64_t> coverage) {
                merged_intervals.emplace_back(start, stop, coverage);
            });
        trees[ref_intervals.first] = merged_intervals;
    }

    return trees;
}

void mask_paf_line(const string& paf_line, int64_t min_length, const unordered_map<string, CoverageIntervalTree>& ref_to_intervals) {
    // split into array of tokens
    vector<string> toks;
    split_delims(paf_line, "\t\n", toks);

    if (toks.size() < 12) {
        throw runtime_error("too few tokens in PAF line: " + paf_line);
    }

    string& query_name = toks[0];
    int64_t query_length = stol(toks[1]);
    int64_t query_start = stol(toks[2]);
    int64_t query_end = stol(toks[3]) - 1;

    vector<CoverageInterval> overlapping_intervals;
    unordered_map<string, CoverageIntervalTree>::const_iterator map_it = ref_to_intervals.find(query_name);
    if (map_it != ref_to_intervals.end()) {
        overlapping_intervals = map_it->second.findOverlapping(query_start, query_end);
    }
#ifdef debug
    cerr << "Found " << overlapping_intervals.size() << " overlaps for paf line " << query_name << " " << query_start << "-" << query_end << ":" << endl;
    for (auto& overlap : overlapping_intervals) {
        cerr << overlap << endl;
    }
#endif

    if (overlapping_intervals.empty()) {
        // nothing to mask
        cout << paf_line << "\n";
    }

    vector<CoverageInterval> remaining_intervals;
    remaining_intervals.emplace_back(query_start, query_end, pair<int64_t, int64_t>(0, 0));

    // todo: we can do this more efficiently (but in practice, there's only ever going to be a small number of overlaps)
    for (CoverageInterval& overlap : overlapping_intervals) {
        vector<CoverageInterval> cut_intervals;
        for (CoverageInterval& paf_interval : remaining_intervals) {
#ifdef debug
            cerr << "subtracting " << paf_interval << " minus " << overlap << endl;
#endif
            interval_subtract(paf_interval, overlap, cut_intervals);
        }
        remaining_intervals = cut_intervals;
    }

    std::sort(remaining_intervals.begin(), remaining_intervals.end(), [](const CoverageInterval& a, const CoverageInterval& b) {
            return a.start < b.start;
        });

    for (CoverageInterval& remaining_interval : remaining_intervals) {
        clip_paf(toks, query_name, query_length, query_start, query_end, remaining_interval, min_length);        
    }
}

void interval_subtract(CoverageInterval& interval_a, CoverageInterval& interval_b, vector<CoverageInterval>& out_fragments) {

    if (interval_b.start <= interval_a.start && interval_b.stop >= interval_a.stop) {
        // there's nothing left!
#ifdef debug
        cerr << "subtraction deltes everything" << endl;
#endif
        return;
    }

    if (interval_b.start >= interval_a.start && interval_b.start < interval_a.stop) {
        // b overlaps to the right
        out_fragments.emplace_back(interval_a.start, interval_b.start - 1, make_pair(0, 0));
#ifdef debug
        cerr << "subtraction adds left fragment " << out_fragments.back() << endl;
#endif        
    }

    if (interval_b.stop >= interval_a.start && interval_b.stop < interval_a.stop) {
        // b overlaps to the left
        out_fragments.emplace_back(interval_b.stop + 1, interval_a.stop, make_pair(0, 0));
#ifdef debug
        cerr << "subtraction adds right fragment " << out_fragments.back() << endl;
#endif
    }

}

static void clip_paf(const vector<string>& toks, const string& query_name, int64_t query_length, int64_t query_start, int64_t query_end,
                     CoverageInterval& interval, int64_t min_length) {
    cerr << "Clipping " << interval << " out of " << query_name << " " << query_length << " " << query_start << " " << query_end << endl;
}
