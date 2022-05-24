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

// TODO: we have two ways of filtering below, dominates and dominates_mzgaf()
// need to figure out which one is better (all signs point to the old one, which,
// if I remember, is Gospel From Heng.  If that's the case we can just forget the new one
// if not here, then at least hide the options in cactus config.

// test if one record "dominates" another, using primary/secondary, mapq, block length in that order
static bool dominates(const GafRecord& gaf1, const GafRecord& gaf2, double ratio) {
    bool primary1 = !gaf1.opt_fields.count("tp") || gaf1.opt_fields.at("tp").second == "P";
    bool primary2 = !gaf2.opt_fields.count("tp") || gaf2.opt_fields.at("tp").second == "P";
    // empty interval can't dominate
    if (gaf1.query_start >= gaf1.query_end) {
        return false;
    } else if (gaf2.query_start >= gaf2.query_end) {
        return true;
    }
    if (primary1 && !primary2) {
        return true;
    } else if (primary2 && !primary1) {
        return false;
    }
    if ((double)gaf1.mapq / ((double)gaf2.mapq + 0.000001) >= ratio) {
        return true;
    } else if ((double)gaf2.mapq / ((double)gaf1.mapq + 0.000001) >= ratio) {
        return false;
    }
    if ((double)gaf1.block_length / ((double)gaf2.block_length + 0.000001) >= ratio) {
        return true;
    } else if ((double)gaf2.block_length / ((double)gaf1.block_length + 0.000001) >= ratio) {
        return false;
    }    
    return false;
}

// test if one record "dominates" another, using same logic as mzgaf2paf 
static bool dominates_mzgaf2paf(const GafRecord& gaf1, const GafRecord& gaf2, int64_t query_overlap_threshold) {
    return (gaf1.block_length >= query_overlap_threshold && gaf2.block_length < query_overlap_threshold ||
            gaf1.block_length < query_overlap_threshold && gaf2.block_length < query_overlap_threshold);
}

// measure the overlap
static int64_t overlap_size(const GafRecord& gaf1, const GafRecord& gaf2) {
    int64_t ostart = std::max(gaf1.query_start, gaf2.query_start);
    int64_t oend = std::min(gaf1.query_end, gaf2.query_end);
    assert(oend >= ostart);
    return oend - ostart;
}

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <gaf> > output.gaf" << endl
         << "Filter GAF record if its query interval overlaps another query interval and\n"
         << "  1) the record is secondary and the overlapping record is primary or\n"
         << "  2) the record's MAPQ is lower than {ratio, see -r} times the overlapping record's MAPQ or\n"
         << "  3) the record's block length is less than {ratio, see -r} times larger than the overlapping record's block length (and its MAPQ isn't higher)" << endl
         << "  Also: the -o option can be used to mimic mzgaf2paf's query overlap filter" << endl
         << endl
         << "options: " << endl
         << "    -r, --ratio N                   If two query blocks overlap, and one is Nx bigger than the other, the bigger one is kept (otherwise both deleted) [0]" << endl
         << "    -m, --min-overlap N             Ignore overlaps that consitute <N% of the length [0]" << endl
         << "    -o, --min-overlap-length N      If >= 2 query regions with size >= N overlap, ignore the query region.  If 1 query region with size >= N overlaps any regions of size <= N, ignore the smaller ones only. Works separate to -r/-m but can be used in conjunction with them to combine the two filters (0 = disable) [0]" << endl
         << "    -q, --min-mapq N                Don't let an interval with MAPQ < N cause something to be filtered out" << endl
         << "    -b, --min-block-length N           Don't let an interval with block length < N cause something to be filtered out" << endl
         << "    -p, --paf                       Input is PAF, not GAF" << endl;
}    

int main(int argc, char** argv) {

    double ratio = 0.;
    double min_overlap_pct = 0.;
    int64_t min_overlap_len = 0;
    int64_t min_block_len = 0;
    int64_t min_mapq = 0;    
    
    int c;
    bool is_paf = false;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"ratio", required_argument, 0, 'r'},
            {"min-overlap", required_argument, 0, 'm'},
            {"min-overlap-length", required_argument, 0, 'o'},
            {"min-block-length", required_argument, 0, 'b'},
            {"min-mapq", required_argument, 0, 'q'},            
            {"paf", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h:r:m:po:b:q:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            ratio = stof(optarg);
            break;
        case 'm':
            min_overlap_pct = stof(optarg);
            break;
        case 'o':
            min_overlap_len = std::stol(optarg);
            break;            
        case 'p':
            is_paf = true;
            break;
        case 'b':
            min_block_len = std::stol(optarg);
            break;
        case 'q':
            min_mapq = std::stol(optarg);
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

    if ((ratio == 0) && (min_overlap_len == 0)) {
        cerr << "[gaffilter] error: at least one of -r or -o must be used to specify filter" << endl;
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

    // shimmy in paf support post hoc (at the cost of storing a dummy gaf record list in memory!)
    vector<PafLine> paf_records;
    function<string(const GafRecord&)> print_record = [&](const GafRecord& gaf_record) {
        stringstream ss;
        if (is_paf) {
            // hack alert: hijack path_length with offset in paf_records
            ss << paf_records[gaf_record.path_length];
        } else {
            ss << gaf_record;
        }
        return ss.str();
    };

    // just load the gaf into memory
    vector<GafRecord> gaf_records;
    string line_buffer;
    while (getline(*in_stream, line_buffer)) {
        if (line_buffer[0] == '*') {
            // skip -S stuff
            continue;
        }        
        GafRecord gaf_record;
        if (is_paf) {
            PafLine paf_record = parse_paf_line(line_buffer);
            paf_records.push_back(paf_record);
            // just copy what we (might) need
            gaf_record.query_name = paf_record.query_name;
            gaf_record.query_length = paf_record.query_len;
            gaf_record.query_start = paf_record.query_start;
            gaf_record.query_end = paf_record.query_end;
            gaf_record.strand = paf_record.strand;
            gaf_record.mapq = paf_record.mapq;
            if (paf_record.opt_fields.count("gl")) {
                gaf_record.block_length = stol(paf_record.opt_fields.at("gl").second);
            } else {
                gaf_record.block_length = paf_record.num_bases;
            }
            if (paf_record.opt_fields.count("tp")) {
                gaf_record.opt_fields["tp"] = paf_record.opt_fields.at("tp");
            }
            if (paf_record.opt_fields.count("rc")) {
                gaf_record.opt_fields["rc"] = paf_record.opt_fields.at("rc");
            }
            // hack alert: hijack path_length with offset in paf_records
            gaf_record.path_length = paf_records.size() - 1;
        } else {
            parse_gaf_record(line_buffer, gaf_record);
        }
        gaf_records.push_back(gaf_record);
    }
    cerr << "[gaffilter]: Loaded " << gaf_records.size() << (is_paf ? " PAF" : " GAF") << " records" << endl;

    // make an interval tree for each query sequence
    unordered_map<string, vector<GafInterval>> gaf_intervals;
    for (const auto& gaf_record : gaf_records) {
        int64_t end_point = gaf_record.query_end;
        if (end_point > gaf_record.query_start) {
            // interval tree expects closed coordinates.  but it also expects the end point >= start point
            --end_point;
        }
        GafInterval gaf_interval(gaf_record.query_start, gaf_record.query_end - 1, &gaf_record);
        gaf_intervals[gaf_record.query_name].push_back(gaf_interval);
    }
    unordered_map<string, GafIntervalTree*> gaf_trees;
    for (const auto& qi : gaf_intervals) {
        gaf_trees[qi.first] = new GafIntervalTree(qi.second);
    }
    gaf_intervals.clear();
    cerr << "[gaffilter]: Constructed interval trees" << endl;


    int64_t filter_count = 0;
    int64_t filter_len_count = 0;
        
    // simple algorithm:
    // for each record, scan its overlaps and flag it if it finds anything
    // that overlaps that isn't ratio X smaller.
    // this is a really inefficient in worst-case (where everything overlaps) but that's not at all what we expect
    for (int64_t i = 0; i < gaf_records.size(); ++i) {
        int64_t end_point = gaf_records[i].query_end;
        if (end_point > gaf_records[i].query_start) {
            // interval tree expects closed coordinates.  but it also expects the end point >= start point
            --end_point;
        }
        vector<GafInterval> overlapping;
        string ref_contig;
        if (gaf_records[i].opt_fields.count("rc")) {
            ref_contig = gaf_records[i].opt_fields.at("rc").second;
        }
        gaf_trees[gaf_records[i].query_name]->visit_overlapping(gaf_records[i].query_start, end_point, [&](const GafInterval& interval) {
                // filter self alignments 
                if (interval.value != &gaf_records[i] &&
                    //and mapq/block length failing alignments
                    interval.value->mapq >= min_mapq && interval.value->block_length >= min_block_len) {
                    string overlap_contig;
                    if (interval.value->opt_fields.count("rc")) {
                        overlap_contig = interval.value->opt_fields.at("rc").second;
                    }
                    // also ignore things that map to different contigs
                    if (ref_contig == overlap_contig || ref_contig.empty() || overlap_contig.empty()) {
                        int64_t overlap_bases = overlap_size(gaf_records[i], *interval.value);
                        // filter overlaps that are too small to matter (via min_overlap_pct)
                        if (gaf_records[i].block_length == 0 ||
                            (double)overlap_bases / (double)gaf_records[i].block_length >= min_overlap_pct) {
                            overlapping.push_back(interval);
                        }
                    }
                }
            });
        bool is_dominant = true;
        for (const auto& ogi : overlapping) {
            if (ratio) {
                is_dominant = dominates(gaf_records[i], *ogi.value, ratio);
            }
            if (is_dominant && min_overlap_len) {
                is_dominant = dominates_mzgaf2paf(gaf_records[i], *ogi.value, min_overlap_len);
            }
            if (!is_dominant) {
                break;
            }
        }
        if (is_dominant) {
            cout << print_record(gaf_records[i]) << "\n";
        } else {
            ++filter_count;
            if (is_paf) {
                filter_len_count += paf_records[i].num_bases;
            } else {
                filter_len_count += gaf_records[i].block_length;
            }
#ifdef debug
            cerr << "\nfiltering record " << i << " (" << &gaf_records[i] << ") because it doesn't dominate its "
                 << overlapping.size() << " overlaps\n  " << print_record(gaf_records[i]) << endl;
            int64_t ocount = 0;
            for (const auto& ogi : overlapping) {
                cerr << "overlap " << ocount++ << " (" << ogi.value << "):\n  " << print_record(*ogi.value) << endl;
            }
#endif
        }
    }

    for (auto qt : gaf_trees) {
        delete qt.second;
    }
    
    cerr << "[gaffilter]: filtered " << filter_count << " / " << gaf_records.size() << ". total block lengths filtered: " << filter_len_count << endl;
    return 0;
}
