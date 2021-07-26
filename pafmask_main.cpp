#include <unistd.h>
#include <getopt.h>

#include "pafmask.hpp"

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
    size_t masked_bases = 0;
    while (getline(*in_paf, buffer)) {
        masked_bases += mask_paf_line(buffer, min_length, ref_to_intervals, validate);
    }

    cerr << "[pafmask]: clipped out: " << masked_bases << " bp" << endl;

    return 0;
}

