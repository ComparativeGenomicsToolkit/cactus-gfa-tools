#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <sys/stat.h>
#include "rgfa-split.hpp"
#include "pafcoverage.hpp"

using namespace std;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options]" << endl
       << "Partition rGFA nodes into reference contigs.  Input must be uncompressed GFA (not stdin)" << endl
       << "input options: " << endl
       << "    -g, --rgfa FILE                         rGFA to use as baseline for contig splitting (if not defined, minmap2 output assumed)" << endl
       << "    -m, --input-contig-map FILE             Use tsv map (computed with -M) instead of rGFA" << endl
       << "    -p, --paf FILE                          PAF file to split" << endl
       << "    -B, --bed FILE                          BED file.  Used to subtract out softmasked regions when computing coverage (multiple allowed)" << endl
       << "output options: " << endl 
       << "    -b, --output-prefix PREFIX              All output files will be of the form <PREFIX><contig>.paf/.fa_contigs" << endl
       << "    -M, --output-contig-map FILE            Output rgfa node -> contig map to this file" << endl
       << "    -G, --split-gfa                         Split the input GFA too and output <PREFIX><config>.gfa files" << endl
       << "contig selection options: " << endl
       << "    -q, --contig-prefix PREFIX              Only process contigs beginning with PREFIX" << endl
       << "    -c, --contig-name NAME                  Only process NAME (multiple allowed)" << endl
       << "    -C, --contig-file FILE                  Path to list of contigs to process" << endl
       << "    -o, --other-name NAME                   Lump all contigs not selected by above options into single reference with name NAME" << endl
       << "contig assignment ambiguity handling options: " << endl
       << "    -n, --min-query-coverage FLOAT          At least this fraction of input contig must align to reference contig for it to be assigned (can repeat)" << endl
       << "    -T, --small-coverage-threshold N        Used to toggle between the coverage thresholds (-n). Should have one-fewer value than -n" << endl
       << "    -Q, --min-query-uniqueness FLOAT        The ratio of the number of query bases aligned to the chosen ref contig vs the next best ref contig must exceed this threshold to not be considered ambigious" << endl
       << "    -u, --min-query-chunk N                 I a query interval of >= N bp aligns to a reference with sufficient coverage, cut it out.  Disabled when 0. [0]" << endl
       << "    -s, --allow-softclip                    Allow softclipping with -u" << endl
       << "    -P, --max-gap N                         Count cigar gaps of length <= N towards coverage" << endl
       << "    -a, --ambiguous-name NAME               All query contigs that do not meet min coverage (-n) assigned to single reference with name NAME" << endl
       << "    -A, --min-mapq N                        Don't count PAF lines with MAPQ<N towards coverage" << endl
       << "    -r, --reference-prefix PREFIX           Don't apply ambiguity filters to query contigs with this prefix" << endl
       << "    -L, --log FILE                          Keep track of filtered and assigned contigs in given file [stderr]" << endl;
}    

int main(int argc, char** argv) {

    // input
    string rgfa_path;
    string input_contig_map_path;
    string input_paf_path;
    string bed_path;

    // output
    string output_prefix;
    string output_contig_map_path;
    bool split_gfa = false;

    // contig selection
    string contig_prefix;
    unordered_set<string> contig_names;
    string contig_names_path;
    size_t selection_options = 0;
    string other_name;

    // ambiguity handling
    vector<double> min_query_coverages;
    vector<int64_t> small_coverage_thresholds;
    double min_query_uniqueness = 0;
    int64_t min_query_chunk = 0;
    bool allow_softclip = false;
    int64_t max_gap = 0;
    string ambiguous_name;
    string reference_prefix;
    int64_t min_mapq = 0;
    string log_path;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"rgfa", required_argument, 0, 'g'},
            {"input-contig-map", required_argument, 0, 'm'},
            {"paf", required_argument, 0, 'p'},
            {"bed", required_argument, 0, 'B'},
            {"output-prefix", required_argument, 0, 'b'},
            {"output-contig-map", required_argument, 0, 'M'},
            {"split-gfa", no_argument, 0, 'G'},
            {"contig-prefix", required_argument, 0, 'q'},
            {"contig-name", required_argument, 0, 'c'},
            {"contig-file", required_argument, 0, 'C'},
            {"other-name", required_argument, 0, 'o'},            
            {"min-query-coverage", required_argument, 0, 'n'},
            {"small-coverage-threshold", required_argument, 0, 'T'},
            {"min-query-uniqueness", required_argument, 0, 'Q'},
            {"min-query-chunk", required_argument, 0, 'u'},
            {"allow-softlicp", no_argument, 0, 's'},            
            {"max-gap", required_argument, 0, 'P'},
            {"ambiguous-name", required_argument, 0, 'a'},
            {"min-mapq", required_argument, 0, 'A'},
            {"reference-prefix", required_argument, 0, 'r'},
            {"log", required_argument, 0, 'L'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hg:m:p:B:b:M:Gq:c:C:o:n:T:Q:u:sP:a:A:r:L:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            rgfa_path = optarg;
            break;
        case 'm':
            input_contig_map_path = optarg;
            break;
        case 'p':
            input_paf_path = optarg;
            break;
        case 'B':
            bed_path = optarg;
            break;
        case 'b':
            output_prefix = optarg;
            break;
        case 'M':
            output_contig_map_path = optarg;
            break;
        case 'G':
            split_gfa = true;
            break;
        case 'q':
            contig_prefix = optarg;
            break;
        case 'c':
            contig_names.insert(optarg);
            break;
        case 'C':
            contig_names_path = optarg;
            break;
        case 'o':
            other_name = optarg;
            break;            
        case 'n':
            min_query_coverages.push_back(stof(optarg));
            break;
        case 'T':
            small_coverage_thresholds.push_back(stol(optarg));
            break;
        case 'Q':
            min_query_uniqueness = stof(optarg);
            break;
        case 'u':
            min_query_chunk = stol(optarg);
            break;
        case 's':
            allow_softclip = true;
            break;
        case 'P':
            max_gap = stol(optarg);
            break;
        case 'a':
            ambiguous_name = optarg;
            break;
        case 'A':
            min_mapq = stol(optarg);
            break;
        case 'r':
            reference_prefix = optarg;
            break;
        case 'L':
            log_path = optarg;
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

    if (optind != argc - 0) {
        cerr << "[rgfa-split] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    if (selection_options > 1) {
        cerr << "[rgfa-split] error: -q, -c and -C are mutually exclusive" << endl;
        return 1;
    }    

    if (rgfa_path == "-") {
        cerr << "[rgfa-split] error: - (stdin) not supported for rgfa" << endl;
        return 1;
    };

    if (!output_prefix.empty() && output_prefix.back() == '/') {
        // note: not checking for error here (it'll show up downstream if directory doesn't exist)
        mkdir(output_prefix.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if ((!min_query_coverages.empty()  || min_query_uniqueness > 1) && ambiguous_name.empty()) {
        cerr << "[rgfa-split] error: ambiguous name must be set with -a when using -n or -Q" << endl;
        return 1;
    }

    if ((!min_query_coverages.empty() || !small_coverage_thresholds.empty()) &&
        small_coverage_thresholds.size() != min_query_coverages.size() - 1) {
        cerr << "[rgfa-split] error: If there are K min coverages specified with -n, there must be K-1 thresholds"
            " specified with -T" << endl;
        return 1;
    }
    for (size_t i = 1; i < small_coverage_thresholds.size(); ++i) {
        if (small_coverage_thresholds[i] <= small_coverage_thresholds[i-1]) {
            cerr << "[rgfa-split] error: Thresholds specified by -T must be in strictly increasing order" << endl;
            return 1;
        }
    }
    
    ofstream log_fstream;
    if (!log_path.empty()) {
        log_fstream.open(log_path);
        if (!log_fstream) {
            cerr << "[rgfa-split] error: Cannot open log file " << log_path << endl;
            return 1;
        }
    }
    ostream& log_stream = !log_path.empty() ? log_fstream : cerr;

    function<void(const string&)> check_ifile = [&](const string& path) {
        ifstream in_stream(path);
        if (!in_stream) {
            cerr << "[rgfa-split] error: unable to open input file \"" << path << "\"" << endl;
            exit(1);
        }
    };
    
    // get the partition of GFA nodes -> reference contig
    pair<unordered_map<int64_t, int64_t>, vector<string>> partition;
    unordered_map<string, int64_t> target_to_id;
    if (!rgfa_path.empty()) {
        // compute from minigraph output
        check_ifile(rgfa_path);
        partition = rgfa2contig(rgfa_path);
    } else if (!input_contig_map_path.empty()) {
        // load table
        check_ifile(input_contig_map_path);
        partition = load_contig_map(input_contig_map_path);
    } else {
        // load a paf file int the existing structs (bit of a hack)
        check_ifile(input_paf_path);
        ifstream paf_file(input_paf_path);
        string paf_line;
        while (getline(paf_file, paf_line)) {
            vector<string> toks;
            split_delims(paf_line, "\t\n", toks);
            if (toks.size() > 5 && !target_to_id.count(toks[5])) {
                target_to_id[toks[5]] = partition.second.size();
                partition.second.push_back(toks[5]);
            }
        }   
    }
    
    // output the contig map if path given
    if (!output_contig_map_path.empty()) {
        ofstream output_contig_map_file(output_contig_map_path);
        if (!output_contig_map_file) {
            cerr << "[rgfa-split] error: unable to open output contig map file \"" << output_contig_map_path << "\"" << endl;
            return 1;
        }
        for (auto& node_contig : partition.first) {
            output_contig_map_file << "S" << node_contig.first << "\t" << partition.second[node_contig.second] << "\n";
        }
    }

    // load the contig names file
    if (!contig_names_path.empty()) {
        ifstream contig_names_file(contig_names_path);
        if (!contig_names_file) {
            cerr << "[rgfa-split] error: unable to open contig names file \"" << contig_names_path << "\"" << endl;
            return 1;
        }
        string buffer;
        while (getline(contig_names_file, buffer)) {
            vector<string> toks;
            split_delims(buffer, "\t\n", toks);
            if (!toks.empty() && !toks[0].empty()) {
                contig_names.insert(toks[0]);
            }
        }
    }

    // load up contig selection options 
    function<bool(const string&)> visit_contig = [&](const string& name) -> bool {
        if (!contig_names.empty()) {
            return contig_names.count(name);
        } else if (!contig_prefix.empty()) {
            return name.substr(0, contig_prefix.length()) == contig_prefix;
        } else {
            return true;
        }
    };

    // make a coverage threshold map
    map<int64_t, double> cov_threshold_map;
    if (min_query_coverages.empty()) {
        // no threshold
        cov_threshold_map[numeric_limits<int64_t>::max()] = 0.;
    } else {
        for (size_t i = 0; i < min_query_coverages.size() - 1; ++i) {
            cov_threshold_map[small_coverage_thresholds[i]] = min_query_coverages[i];
        }
        cov_threshold_map[numeric_limits<int64_t>::max()] = min_query_coverages.back();
    }

    // lump unselected contigs into the "other" category if desired
    int64_t other_id = -1;
    if (!other_name.empty()) {
        if (target_to_id.empty()) {
            set_other_contig(partition.first, partition.second, visit_contig, other_name);
        } else {
            set_other_contig_paf(target_to_id, partition.second, visit_contig, other_name);
        }
        visit_contig = [&](const string&) -> bool {
            return true;
        };
        other_id = partition.second.size();
    }
    
    // set the ambigious name
    int64_t ambiguous_id = -1;
    if (!ambiguous_name.empty()) {
        assert(std::find(partition.second.begin(), partition.second.end(), ambiguous_name) == partition.second.end());
        ambiguous_id = partition.second.size();
        partition.second.push_back(ambiguous_name);
    }

    // make sure the bed file exists
    unordered_map<string, int64_t> query_mask_stats;
    if (!bed_path.empty()) {
        check_ifile(bed_path);
        query_mask_stats = load_query_mask_stats(bed_path);
    }

    // split the paf into one paf per reference contig.  also output a list of query contig names
    // alongside than can be used to split the fasta with samtools faidx
    if (!input_paf_path.empty()) {
        check_ifile(input_paf_path);
        // toggle how we want to handle target names, depending if the PAF comes from minigraph or minimap2
        function<int64_t(const string&)> name_to_refid;
        if (!rgfa_path.empty()) {
            name_to_refid = [&](const string& target_name) {
                // use the map to go from the target name (rgfa node id in this case) to t
                // the reference contig (ex chr20)
                int64_t target_id = node_id(target_name);
                assert(partition.first.count(target_id));
                int64_t reference_id = partition.first.at(target_id);
                return reference_id;
            };
        } else {
            name_to_refid = [&](const string& target_name) {
                // just use a trivial map expeting the contig names from minimap2
                assert(target_to_id.count(target_name));
                return target_to_id[target_name];
            };
        }
        paf_split(input_paf_path, name_to_refid, partition.second, visit_contig, output_prefix,
                  cov_threshold_map, min_query_uniqueness,
                  min_query_chunk, allow_softclip, other_id,
                  ambiguous_id, reference_prefix, query_mask_stats, max_gap, min_mapq, log_stream);
    }

    // split the gfa
    if (!rgfa_path.empty() && split_gfa) {
        gfa_split(rgfa_path, partition.first, partition.second, visit_contig, output_prefix);
    }

    return 0;
}
