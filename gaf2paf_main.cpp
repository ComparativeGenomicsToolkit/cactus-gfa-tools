/*
  Convert GAF from minigraph -c to PAF (in stable coordinates)
 */

#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <list>
#include <cassert>
#include <cmath>

#include "gafkluge.hpp"
#include "paf.hpp"
#include "pafcoverage.hpp"

//#define debug

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

/* convert a GAF line to PAF lines, one per path step (gaf2paf_split in paf.hpp).  A record that
   gaffilter -x cut carries the query intervals it keeps in a kq:Z: tag: each line is then cut down
   to them, and its tags still describe the whole parent record (gm/gl/gi), which is what the
   block-length filter downstream has to see */
static void gaf2paf(const GafRecord& gaf_record, const unordered_map<string, int64_t>& len_map, ostream& os) {
    vector<gaf2paf_split::KeptInterval> kept;
    bool cut = false;
    if (gaf_record.opt_fields.count("kq")) {
        if (!gaf2paf_split::parse_kept_intervals(gaf_record.opt_fields.at("kq").second, kept)) {
            cerr << "[gaf2paf] error: malformed kq:Z: tag on a record of " << gaf_record.query_name << ": "
                 << gaf_record.opt_fields.at("kq").second << endl;
            exit(1);
        }
        cut = true;
    }
    gaf2paf_split::for_each_paf_line(gaf_record, len_map, [&](const PafLine& paf_record, const string& cigar) {
            if (!cut) {
                gaf2paf_split::print_paf_line(os, gaf_record, paf_record, cigar);
            } else {
                gaf2paf_split::cut_paf_line(paf_record, cigar, kept, [&](const PafLine& piece, const string& piece_cigar) {
                        gaf2paf_split::print_paf_line(os, gaf_record, piece, piece_cigar);
                    });
            }
        });
}

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <gaf> [gaf2] [gaf3] [...] > output.paf" << endl
         << "Convert minigraph GAF to PAF" << endl
         << endl
         << "options: " << endl
         << "    -l, --lengths FILE      TSV with contig length as first two columns (.fai will do)." << endl
         << endl
         << "A record with a kq:Z:s1-e1,s2-e2,... tag (written by gaffilter -x on a record it cut) keeps only those" << endl
         << "query intervals: each of its lines is cut down to them, keeping the parent record's gm/gl/gi." << endl;
}    

int main(int argc, char** argv) {
    // the result goes to standard output, and a failed write there is
    // otherwise reported to nobody
    check_stdout_at_exit();


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
        cerr << "[gaf2paf] error: -l must be specified to produce valid PAF" << endl;
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
                cerr << "[gaf2paf] error: unable to open input: " << in_path << endl;
                return 1;
            }
            in_stream = &in_file;
        }

        GafRecord gaf_record;
        string line_buffer;
        while (getline(*in_stream, line_buffer)) {
            if (line_buffer[0] == '*') {
                // skip -S stuff
                continue;
            }
            parse_gaf_record(line_buffer, gaf_record);
            if (!gaf_record.opt_fields.count("cg")) {
                cerr << "[gaf2paf] error: cg cigar not found. This tool only works on output of minigraph -c" << endl;
                return 1;
            }
            if (gaf_record.strand == '-') {
                gaf2paf_split::flip_gaf(gaf_record, len_map);
            }
            gaf2paf(gaf_record, len_map, cout);
        }
    }
        
    return 0;
}
