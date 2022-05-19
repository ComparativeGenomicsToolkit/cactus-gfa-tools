/*
  Filter GAF records according to query overlap
 */

#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <list>
#include <cassert>

#include "gafkluge.hpp"
#include "paf.hpp"
#include "IntervalTree.h"

//#define debug

using namespace std;
using namespace gafkluge;

typedef IntervalTree<int64_t, const GafRecord*> GafIntervalTree;
typedef GafIntervalTree::interval GafInterval;


static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <gaf> > output.gaf" << endl
         << "Filter GAF record if its query interval overlaps another query interval and\n"
         << "  1) the record is secondary and the overlapping record is primary or\n"
         << "  1) the record's MAPQ is lower than {ratio, see -r} times the overlapping record's MAPQ or\n "
         << "  3) the record's query interval is less than {ratio, see -r} times larger than the overlapping record's query interval" << endl
         << endl
         << "options: " << endl
         << "    -r, --ratio N      If two query blocks overlap, and one is Nx bigger than the other, the bigger one is kept (otherwise both deleted) [2]" << endl;
}    

int main(int argc, char** argv) {

    double ratio = 2.;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"ratio", required_argument, 0, 'r'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:r:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            ratio = stof(optarg);
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
        cerr << "[gaffilter] error: too few arguments" << endl;
        help(argv);
        return 1;
    }
    
    string gaf_path = argv[optind++];
    
    // open the gaf file
    ifstream in_file;
    istream* in_stream;
    if (gaf_path == "-") {
        in_stream = &cin;
    } else {
        in_file.open(gaf_path);
        if (!in_file) {
            cerr << "[gaffilter] error: unable to open input: " << gaf_path << endl;
            return 1;
        }
        in_stream = &in_file;
    }

    // just load the gaf into memory
    vector<GafRecord> gaf_records;
    string line_buffer;
    while (getline(*in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            // skip -S stuff
            continue;
        }
        GafRecord gaf_record;
        parse_gaf_record(line_buffer, gaf_record);
        gaf_records.push_back(gaf_record);
    }    

    // make an interval tree for each query sequence
    unordered_map<string, vector<GafInterval>> gaf_intervals;
    for (const auto& gaf_record : gaf_records) {
        GafInterval gaf_interval(gaf_record.query_start, gaf_record.query_end - 1, &gaf_record);
        gaf_intervals[gaf_record.query_name].push_back(gaf_interval);
    }
    unordered_map<string, GafIntervalTree*> gaf_trees;
    for (const auto& qi : gaf_intervals) {
        gaf_trees[qi.first] = new GafIntervalTree(qi.second);
    }
    gaf_intervals.clear();


    int64_t filter_count = 0;
    int64_t filter_len_count = 0;

    // simple algorithm:
    // for each record, scan it's overlaps and flag it if it finds anything
    // that overlaps that isn't ratio X smaller.
    // this is a really inefficient in worst-case (where everything overlaps) but that's not at all what we expect
    for (int64_t i = 0; i < gaf_records.size(); ++i) {
        bool is_filtered = false;
        int64_t query_len = gaf_records[i].query_end - gaf_records[i].query_start;
        bool is_primary = !gaf_records[i].opt_fields.count("tp") || gaf_records[i].opt_fields.at("tp").second == "P";
        vector<GafInterval> overlapping = gaf_trees[gaf_records[i].query_name]->findOverlapping(
            gaf_records[i].query_start, gaf_records[i].query_end);
        for (const auto& ogi : overlapping) {
            if (ogi.value != &gaf_records[i]) {
                bool overlap_primary = !ogi.value->opt_fields.count("tp") || ogi.value->opt_fields.at("tp").second == "P";
                int64_t overlap_len = ogi.value->query_end - ogi.value->query_start;
                double orat = (double)query_len / ((double)overlap_len + 0.000001);
                double mrat = (double)gaf_records[i].mapq / ((double)ogi.value->mapq + 0.000001);
                // primary / secondary comparison takes priority
                if (overlap_primary && !is_primary) {
#ifdef debug
                    cerr << "filtering record " << i << "\n " << gaf_records[i] << "\ndue to primary overlap with\n "
                         << *ogi.value << endl;
#endif                    
                    is_filtered = true;
                } else {
                    // next is mapq
                    if (1. / mrat >= ratio) {
#ifdef debug
                        cerr << "filtering record " << i << "\n " << gaf_records[i] << "\n because its mapq is dominiated by " << (1./mrat) << " by\n "
                         << *ogi.value << endl;
#endif                                            
                        is_filtered = true;
                        // then interval length
                    } else if (mrat <= ratio && 1. / orat >= ratio) {
#ifdef debug
                        cerr << "filtering record " << i << "\n " << gaf_records[i] << "\n because its query interval is dominiated by " << (1./orat) << " by\n "
                         << *ogi.value << endl;
#endif                    
                        is_filtered = true;
                    }
                }
                if (is_filtered) {
                    ++filter_count;
                    filter_len_count += query_len;
                    break;
                }
            }
        }
        if (!is_filtered) {
            cout << gaf_records[i] << "\n";
        }
    }

    for (auto qt : gaf_trees) {
        delete qt.second;
    }
    
    cerr << "[gaffilter]: filtered " << filter_count << " / " << gaf_records.size() << ". total query intervals filtered: " << filter_len_count << endl;
    return 0;
}
