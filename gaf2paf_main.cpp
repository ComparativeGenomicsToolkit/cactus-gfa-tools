#include <unistd.h>
#include <getopt.h>
#include <fstream>

#include "gafkluge.hpp"
#include "paf.hpp"

#define debug

using namespace std;
using namespace gafkluge;

static void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <gaf> [gaf2] [gaf3] [...] > output.paf" << endl
       << "Convert minigraph GAF to PAF" << endl
       << endl
       << "options: " << endl;
}    

static inline bool consumes_query(char c) {
    return c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X';
}

static inline bool consumes_target(char c) {
    return  c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X';
}

/* convert a GAF line to a PAF line */
static void gaf2paf(const GafRecord& gaf_record, ostream& os) {
    // load up the cg_cigar
    vector<pair<char, int64_t>> cigar;
    for_each_cg(gaf_record, [&](const char& c, const int64_t& s) {
            cigar.push_back(make_pair(c, s));
        });

    // make a template output paf record
    PafLine paf_record;
    paf_record.query_name = gaf_record.query_name;
    paf_record.query_len = gaf_record.query_length;
    paf_record.strand = gaf_record.strand;

    // apparently the cigar doesn't cut evenly at step ends (todo: verify with heng?)
    // so we keep the remainder here
    pair<char, int64_t> remainder = make_pair('M', 0);
    auto next_remainder = remainder;

    int64_t cigar_offset = 0; // relative to cigar array
    int64_t query_base_count = 0; // keep track of bases in query
    int64_t target_base_count = 0; // and target
 
    // for every GAF step
    for (int64_t step_idx = 0; step_idx < gaf_record.path.size(); ++ step_idx) {
        auto step = gaf_record.path[step_idx];
#ifdef debug
        cerr << " -- gaf step step " << step << endl;
        cerr << " -- target base count " << target_base_count << endl;
#endif
        assert(step.is_stable);
        // if the step is just a chromosome, we shimmy it into an interval so we treat consistently
        if (!step.is_interval) {
            step.start = gaf_record.path_start;
            step.end = gaf_record.path_end;
            assert(gaf_record.path.size() == 1);
        }

        // counts for just this step
        int64_t step_target_offset = step.start;
        if (step_idx == 0 && step.is_interval) {
            // tack on the path_start
            step_target_offset += gaf_record.path_start;
        }
        int64_t step_target_base_count = 0;
        int64_t step_query_offset = query_base_count;
        
        int64_t cigar_idx;
        int64_t cigar_matches = 0;
        int64_t cigar_length = 0;
        if (remainder.second) {
            --cigar_idx; // there's some bases left over in the remainder, go back one
        }
        // scan through the cigar until we exceed the target bases in the step
        for (cigar_idx = cigar_offset; cigar_idx < cigar.size() && step_target_offset + step_target_base_count < step.end; ++cigar_idx) {
            pair<char, int64_t> cigar_step = cigar[cigar_idx];
#ifdef debug
//            cerr << " --- cigar step " << cigar_idx << "/" << cigar.size() << " " << cigar_step.first
//                 << " " << cigar_step.second <<  " cur target len = " << cur_target_len << " vs " << step_target_len << endl;
#endif
            // subtract the remainder from the first step
            int64_t corrected_step_length = (cigar_idx == cigar_offset) ? cigar_step.second - remainder.second : cigar_step.second;
            if (consumes_query(cigar_step.first)) {
                query_base_count += corrected_step_length;
            }
            if (consumes_target(cigar_step.first)) {
                step_target_base_count += corrected_step_length;
                target_base_count += corrected_step_length;
            }
            if (cigar_step.first == 'M' || cigar_step.first == '=') {
                cigar_matches += corrected_step_length;
            }
            cigar_length += corrected_step_length;
        }

        // cigar entries can apparently span steps, so we cut by keeping trakc of the "remainder"
        // which is clipped off this paf line and reserved for the next paf line
        auto next_remainder = cigar.at(cigar_idx - 1);
        next_remainder.second = step_target_offset + step_target_base_count - step.end;
        if (next_remainder.second) {
            assert(consumes_target(remainder.first));
            cerr << "remainder " << next_remainder.second << endl;
        }
        if (next_remainder.second) {
            // we need to roll back that next remainder
            if (consumes_query(next_remainder.first)) {
                query_base_count -= next_remainder.second;                
            }
            if (consumes_target(next_remainder.first)) {
                step_target_base_count -= next_remainder.second;
                target_base_count -= next_remainder.second;
            }
            if (next_remainder.first == 'M' || next_remainder.first == '=') {
                cigar_matches -= next_remainder.second;
            }
            cigar_length -= next_remainder.second;
        }

        paf_record.query_start = gaf_record.query_start + step_query_offset;
        paf_record.query_end = paf_record.query_start + query_base_count;
        paf_record.strand = gaf_record.strand; //todo: must combine with step reverse tag
                
        paf_record.target_name = step.name;
        paf_record.target_start = step_target_offset;; //todo: may need to involve strand?
        paf_record.target_end = step_target_offset + step_target_base_count;
        paf_record.target_len = paf_record.target_end; // this is a wrong value, but we can't know the right one

        paf_record.num_matching = cigar_matches;
        paf_record.num_bases  = cigar_length;
        paf_record.mapq = gaf_record.mapq;

        // the main fields
        os << paf_record << "\tcg:Z:";
        // the cigar, prepending the remainder if necessary
        if (remainder.second) {
            os << remainder.second << remainder.first;
        }
        for (int64_t i = cigar_offset; i < cigar_idx; ++i) {
            os << (cigar[i].second - (next_remainder.first ? i == cigar_idx-1 : 0)) << cigar[i].first;
        }
        // todo: are there other optional tags we want to preserve? most would need to be recomputed to be
        // valid on alignment subregion
        if (gaf_record.opt_fields.count("tp")) {
            const auto& tp = gaf_record.opt_fields.at("tp");
            os << "\ttp:" << tp.first << ":" << tp.second;
        }
        remainder = next_remainder;
        cigar_offset = cigar_idx;
        os << "\n";
    }
    cerr << "tbc " << target_base_count << " == " << gaf_record.path_end << "-" << gaf_record.path_start << "="
         << (gaf_record.path_end-gaf_record.path_start) << endl;
    //assert(target_base_count == gaf_record.path_end - gaf_record.path_start);
}

int main(int argc, char** argv) {

    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:",
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
        cerr << "[gaf2paf] error: too few arguments" << endl;
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
            gaf2paf(gaf_record, cout);
        }

    }
        
    return 0;
}
