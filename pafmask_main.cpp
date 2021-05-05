#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include "pafcoverage.hpp"
#include "rgfa-split.hpp"

//#define debug

using namespace std;

// make a map of sequence name -> interval tree from a bed file
// intervals get merged up if they overlap or are within padding of each other
static unordered_map<string, CoverageIntervalTree> load_bed(istream& bed_stream, int64_t padding);

// cut interval_b out of interval_a, any peices of a that remain get added to out_fragments
// the number of such pieces is 0 (a in b), 1 (b overlaps one end of a) or 2 (b in a)
static void interval_subtract(CoverageInterval& interval_a, CoverageInterval& interval_b,
                              vector<CoverageInterval>& out_fragments);

// subtract all masked intervals from a paf line and output what's left
static void mask_paf_line(const string& paf_line, int64_t min_length, const unordered_map<string,
                          CoverageIntervalTree>& ref_to_intervals, bool validate);

// output paf line(s) corresponding to given sub-interval of the original paf line
static void clip_paf(const vector<string>& toks, const string& query_name, int64_t query_length, int64_t query_start, int64_t query_end,
                     CoverageInterval& interval, int64_t min_length, bool validate);

// make sure every homology in the fragment_paf corresponds to a homology in toks
static void validate_paf(const vector<string>& toks, const string& fragment_paf);

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <paf> <bed>" << endl
       << "Cut masked regions out of a paf file" << endl
       << endl
       << "options: " << endl
       << "    -m, --min-length N           Remove any remaining intervals less than N bp" << endl
       << "    -p, --padding N              Merge up bed intervals close than this [100]" << endl
       << "    -v, --validate               Validate every cigar to make sure it's consistent with input" << endl;
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
            {"min-length", required_argument, 0, 'm'},
            {"padding", required_argument, 0, 'p'},
            {"validate", no_argument, 0, 'v'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hm:p:v",
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
        case 'v':
            validate = true;
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
        mask_paf_line(buffer, min_length, ref_to_intervals, validate);
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

void mask_paf_line(const string& paf_line, int64_t min_length, const unordered_map<string,
                   CoverageIntervalTree>& ref_to_intervals, bool validate) {
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
#ifndef debug
        // no need to proceed, but it is a sanity check to clip nothing out and get something valid
        // so leave it in when debugging.
        return;
#endif
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
        clip_paf(toks, query_name, query_length, query_start, query_end, remaining_interval, min_length, validate);        
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

    if (interval_b.start > interval_a.start && interval_b.start < interval_a.stop) {
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
                     CoverageInterval& interval, int64_t min_length, bool validate) {
#ifdef debug
    cerr << "Clipping " << interval << " out of " << query_name << " " << query_length << " " << query_start << " " << query_end << endl;
#endif
    if (interval.stop - interval.start + 1 < min_length) {
        return;
    }

    int64_t target_start = stol(toks[7]);
    int64_t target_end = stol(toks[8]);

    //int64_t start_delta = toks[4] == "+" ? interval.start - query_start : query_end - interval.stop;
    int64_t start_delta = interval.start - query_start;
    int64_t new_length = interval.stop - interval.start + 1; 

    // do the cigar
    int64_t query_offset = 0; // where we are in the cigar (in query coordinates)
    int64_t query_len = 0; // how much cigar we've written (in query coordinates)
    int64_t target_offset = 0; // where we are in the cigar (in target coordinates)
    int64_t target_len = 0; // how much cigar we've written (in target coordinates)
    int64_t target_start_offset = -1;
        
    vector<string> new_cigar_toks;
    int64_t new_match_len = 0;
    int64_t new_total_len = 0;
    bool in_range = false;

    vector<pair<int64_t, char>> cigar_toks;    
    for (int i = 12; i < toks.size(); ++i) {
        if (toks[i].substr(0, 5) == "cg:Z:") {
            // todo: quadratic alert: we are scanning the full cigar here
            for_each_cg(toks[i], [&](const string& val, const string& cat) {
                    assert(cat == "M" || cat == "I" || cat =="D");
                    cigar_toks.push_back(make_pair(stol(val), cat[0]));
                });
            break;
        }
    }

    // cigars are backwards if reverse strand
    if (toks[4] == "-") {
        std::reverse(cigar_toks.begin(), cigar_toks.end());
    }

    for (pair<int64_t, char>& cigar_tok : cigar_toks) {
        int64_t len = cigar_tok.first;
        char cat = cigar_tok.second;
        if (cat == 'M' || cat == 'I') {
            in_range = query_offset + len > start_delta && query_len < new_length;
                    
            int64_t left_clip = 0;
            if (in_range && query_offset + len > start_delta && query_offset < start_delta) {
                // we need to start this cigar on a cut
                left_clip = start_delta - query_offset;                            
            }

            int64_t right_clip = 0;
            if (in_range && query_len + len - left_clip > new_length) {
                // we need to end this cigar on a cut                            
                right_clip = query_len + len - left_clip - new_length;
            }

            if (in_range) {
                // emit the adjusted cigar
                int64_t adj_len = len - left_clip - right_clip;
                new_cigar_toks.push_back(to_string(adj_len) + string(1, cat));
                // add bases the the query length
                query_len += adj_len;
                // add the match bases
                if (cat == 'M') {
                    new_match_len += adj_len;
                }
                // add the total bases
                new_total_len += adj_len;
                // advance the query
                query_offset += len;
                // also need to adjust the target for a match
                if (cat == 'M') {
                    target_len += adj_len;
                }
                if (target_start_offset == -1) {
                    target_start_offset = target_offset + (cat == 'M' ? left_clip : 0);
                }
            }
            // advance offsets
            if (cat == 'M') {
                target_offset += len;
            }
            query_offset += len;

            if (in_range) {
                in_range = query_len < new_length;
            }
                        
        } else if (cat == 'D') {
            if (in_range) {
                new_cigar_toks.push_back(to_string(len) + "D");
                target_len += len;
            }
            target_offset += len;
        } else {
            assert(false);
        }
    }

    // cigars are backwards if reverse strand
    if (toks[4] == "-") {
        std::reverse(new_cigar_toks.begin(), new_cigar_toks.end());
    }
    
    assert(target_start_offset >= 0);
    if (toks[4] == "+") {
#ifdef debug      
        cerr << "target_start = " << target_start << " + " << target_start_offset << endl;
#endif
        target_start += target_start_offset;
#ifdef debug
        cerr << "target_end = " << target_start << " + " << target_len << endl;
#endif
        target_end = target_start + target_len;
    } else {
#ifdef debug
        cerr << "target_end = " << target_end << " - " << target_start_offset << endl;
#endif
        target_end = target_end - target_start_offset;
#ifdef debug
        cerr << "target_start = " << target_end << " - 1 - " << target_len << endl;
#endif
        target_start = target_end - target_len;
    }
    stringstream out_stream;
#ifdef debug
    cerr << "Orig\n";
    for (const string& tok : toks) {
        cerr << tok << "\t";
    }
    cerr << endl;
#endif
    out_stream << query_name << "\t" << query_length << "\t" << interval.start << "\t" << (interval.stop + 1) << "\t"
               << toks[4] << "\t" << toks[5] << "\t" << toks[6] << "\t" << target_start << "\t" << target_end
               << "\t" << new_match_len << "\t" << new_total_len << "\t" << toks[11] << "\t" << "cg:Z:";
    for (const auto& new_cigar_tok : new_cigar_toks) {
        out_stream << new_cigar_tok;
    }
    out_stream << "\n";

    cout << out_stream.str();

    if (validate) {
        validate_paf(toks, out_stream.str());
    }

}

// make sure every homology in the fragment_paf corresponds to a homology in toks
// this doesn't check for homolodies in toks that *should* be in fragment_paf, but it
// should be sufficient to catch glaring bugs (ie with reverse strand)
void validate_paf(const vector<string>& toks, const string& fragment_paf) {

    function<unordered_map<int64_t, int64_t>(const vector<string>&)> extract_homologies = [](const vector<string>& paf_toks) {
        unordered_map<int64_t, int64_t> homos;
        int64_t query_pos = stol(paf_toks[2]);
        int64_t target_start = stol(paf_toks[7]);
        int64_t target_end = stol(paf_toks[8]) - 1;
        int64_t target_offset = 0;
        
        for (int i = 12; i < paf_toks.size(); ++i) {
            if (paf_toks[i].substr(0, 5) == "cg:Z:") {
                for_each_cg(paf_toks[i], [&](const string& val, const string& cat) {
                        int64_t len = stol(val);
                        if (cat == "I") {
                            query_pos += len;
                        } else if (cat == "D") {
                            target_offset += len;
                        } else if (cat == "M") {
                            for (int64_t j = 0; j < len; ++j) {
                                if (paf_toks[4] == "+") {
                                    homos[query_pos + j] = target_start + target_offset + j;
                                } else {
                                    assert(paf_toks[4] == "-");
                                    homos[query_pos + j] = target_end - (target_offset + j); 
                                }
                            }
                            query_pos += len;
                            target_offset += len;
                        } else {
                            assert(false);
                        }
                    });
            }
        }
        return homos;
    };

    vector<string> frag_toks;
    split_delims(fragment_paf, "\t\n", frag_toks);
    assert(frag_toks.size() >= 12);

    unordered_map<int64_t, int64_t> homologies = extract_homologies(toks);
    unordered_map<int64_t, int64_t> frag_homologies = extract_homologies(frag_toks);

    int64_t frag_query_length = stol(frag_toks[1]);
    int64_t frag_query_start = stol(frag_toks[2]);
    int64_t frag_query_end = stol(frag_toks[3]) - 1;

    int64_t frag_target_length = stol(frag_toks[6]);
    int64_t frag_target_start = stol(frag_toks[7]);
    int64_t frag_target_end = stol(frag_toks[8]) -1;

    bool good = true;
    for (int64_t q = frag_query_start; q < frag_query_end; ++q) {
        int64_t frag_tgt = frag_homologies.count(q) ? frag_homologies[q] : -1;
        int64_t orig_tgt = homologies.count(q) ? homologies[q] : -1;
#ifdef debug
        cerr << "query pos " << q << " -> frag: " << frag_tgt << " orig: " << orig_tgt << endl;
#endif
        assert(frag_tgt == orig_tgt);

        if (frag_tgt != -1) {
            assert(frag_tgt >= 0);
            assert(frag_tgt >= frag_target_start);
            assert(frag_tgt <= frag_target_end);
            assert(frag_tgt < frag_target_length);
            assert(q < frag_query_length);
        }
    }
}
