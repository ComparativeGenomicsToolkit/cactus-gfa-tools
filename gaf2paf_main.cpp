#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <unordered_map>
#include <list>

#include "gafkluge.hpp"
#include "paf.hpp"

#define debug

using namespace std;
using namespace gafkluge;

static unordered_map<string, int64_t> get_len_map(const string& lengths_path) {
    unordered_map<string, int64_t> len_map;
    ifstream lengths_file(lengths_path);
    if (!lengths_file) {
        cerr << "[gaf2paf] error: unable to open " << lengths_path << endl;
        exit(1);
    }
    string line_buffer;
    vector<string> toks;
    while (getline(lengths_file, line_buffer)) {
        toks.clear();
        split_delims(line_buffer, "\t", toks);
        if (toks.size() > 1) {
            len_map[toks[0]] = stol(toks[1]);
        }
    }
#ifdef debug
    cerr << "length map " << endl;
    for (const auto& xx : len_map) {
        cerr << " " << xx.first << " ==> " << xx.second << endl;
    }
#endif
    return len_map;
}

typedef pair<char, int64_t> Cig;
typedef list<Cig> Cigar;

static inline bool consumes_query(const Cig& c) {
    return c.first == 'M' || c.first == 'I' || c.first == 'S' || c.first == '=' || c.first == 'X';
}

static inline bool consumes_target(const Cig& c) {
    return  c.first == 'M' || c.first == 'D' || c.first == 'N' || c.first == '=' || c.first == 'X';
}

// cut the cigar at pos removing cut_len. cut_len is inserted as a new element after pos
void cigar_cut(Cigar& cigar, Cigar::iterator pos, int64_t cut_len) {
    assert (cut_len > 0);
    int64_t remainder = pos->second - cut_len;
    assert(remainder > 0);
    auto pos2 = pos;
    ++pos2;
    Cig new_item = make_pair(pos->first, cut_len);
    cigar.insert(pos2, new_item);
    pos->second = remainder;
};

// get the next "target_len" bases worth of cigar, starting at pos, clipping at the end if necessary
pair<Cigar::iterator, Cigar::iterator> cigar_next_by_target(Cigar& cigar, Cigar::iterator pos, int64_t target_len) {
    int64_t cur_len = 0;
    auto pos2 = pos;
    for (pos2 = pos; pos2 != cigar.end() && cur_len < target_len; ++pos2) {
        if (consumes_target(*pos2)) {
            cur_len += pos2->second;
        }
    }
    if (cur_len != target_len) {
        assert(cur_len > target_len);
        int64_t cut_len = cur_len - target_len;
        --pos2;
        cur_len -= pos2->second;
        cigar_cut(cigar, pos2, cut_len);
        cur_len += pos2->second;
        ++pos2;
    }
    assert(cur_len == target_len);
    return make_pair(pos, pos2);
}

/* convert a GAF line to a PAF line */
static void gaf2paf(const GafRecord& gaf_record, const unordered_map<string, int64_t>& len_map, ostream& os) {
    // load up the cg_cigar
    Cigar cigar;
    for_each_cg(gaf_record, [&](const char& c, const int64_t& s) {
            cigar.push_back(make_pair(c, s));
        });

    // make a template output paf record
    PafLine paf_record;
    paf_record.query_name = gaf_record.query_name;
    paf_record.query_len = gaf_record.query_length;
    paf_record.strand = gaf_record.strand;

    Cigar::iterator cigar_pos = cigar.begin();
    
    int64_t query_base_count = 0; // keep track of bases in query
    int64_t target_base_count = 0; // and target
 
    // for every GAF step
    for (int64_t step_idx = 0; step_idx < gaf_record.path.size(); ++ step_idx) {
        auto step = gaf_record.path[step_idx];
        
        assert(step.is_stable);
        // if the step is just a chromosome, we shimmy it into an interval so we treat consistently
        if (!step.is_interval) {
            step.start = gaf_record.path_start;
            step.end = gaf_record.path_end;
            assert(gaf_record.path.size() == 1);
        } else {
            // if the step is the first step, we tack on the path offset
            if (step_idx == 0) {
                step.start += gaf_record.path_start;
            }
            // ditto last
            if (step_idx == gaf_record.path.size() - 1) {
                //cerr << "tbc " << target_base_count << " step " << step.start <<"," << step.end << " pathend " << gaf_record.path_end << endl;
                assert(target_base_count + (step.end - step.start) >= gaf_record.path_end - gaf_record.path_start);
                step.end = step.start + (gaf_record.path_end - gaf_record.path_start - target_base_count);
            }
        }

        // gobble up the step's worth of target bases from the cigar
        // we use target, because that's the only measure we have -- it's embedded in the step
        auto cig_range = cigar_next_by_target(cigar, cigar_pos, step.end - step.start);

        // turn the cigar into a string
        stringstream cig_string;        
        int64_t cig_query_bases = 0;
        int64_t cig_target_bases = 0;
        paf_record.num_matching = 0;
        paf_record.num_bases  = 0;
        // todo: may need to reverse it!!
        for (auto i = cig_range.first; i != cig_range.second; ++i) {
            if (consumes_query(*i)) {
                cig_query_bases += i->second;
            }
            if (consumes_target(*i)) {
                cig_target_bases += i->second;
            }
            if (i->first == 'M' || i->first == '=') {
                paf_record.num_matching += i->second;
            }
            paf_record.num_bases += i->second;
            cig_string << i->second << i->first;
        }

        // make a new paf record
        paf_record.query_start = gaf_record.query_start + query_base_count;
        paf_record.query_end = paf_record.query_start + cig_query_bases;
        paf_record.target_start = step.start;
        paf_record.target_end = step.end;
        assert(step.end - step.start == cig_target_bases);
        
        paf_record.strand = gaf_record.strand; //todo: must combine with step reverse tag

        // go fishing for the target length
        if (!len_map.count(step.name)) {
            cerr << "[gaf2paf] error: unable to find " << step.name << " in lengths map" << endl;
            exit(1);
        }
        paf_record.target_name = step.name;
        paf_record.target_len = len_map.at(step.name);
        
        paf_record.mapq = gaf_record.mapq;

        // output the record
        os << paf_record;
        
        // todo: are there other optional tags we want to preserve? most would need to be recomputed to be
        // valid on alignment subregion
        if (gaf_record.opt_fields.count("tp")) {
            const auto& tp = gaf_record.opt_fields.at("tp");
            os << "\ttp:" << tp.first << ":" << tp.second;
        }

        // output the cigar last
        os << "\tcg:Z:" << cig_string.str();
        
        os << "\n";
        
        // advance our counters
        query_base_count += cig_query_bases;
        target_base_count += cig_target_bases;
        cigar_pos = cig_range.second;
    }
}

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <gaf> [gaf2] [gaf3] [...] > output.paf" << endl
         << "Convert minigraph GAF to PAF" << endl
         << endl
         << "options: " << endl
         << "    -l, --lengths FILE      TSV with contig length as first two columns (.fai will do)." << endl;
}    

int main(int argc, char** argv) {

    string rgfa_path;
    string lengths_path;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"lengths", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:l:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'l':
            lengths_path = optarg;
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
        cerr << "[gaf2paf] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    if (lengths_path.empty()) {
        cerr << "gaf2paf] error: -l must be specified to produce valid PAF" << endl;
        return 1;        
    }    

    auto len_map = get_len_map(lengths_path);

    vector<string> in_paths;
    int stdin_count = 0;
    while (optind < argc) {
        in_paths.push_back(argv[optind++]);
        if (in_paths.back() == "-") {
            ++stdin_count;
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
                cerr << "[mzgaf2paf] error: unable to open input: " << in_path << endl;
                return 1;
            }
            in_stream = &in_file;
        }

        GafRecord gaf_record;
        string line_buffer;
        while (getline(*in_stream, line_buffer)) {
            parse_gaf_record(line_buffer, gaf_record);
            gaf2paf(gaf_record, len_map, cout);
        }
    }
        
    return 0;
}
