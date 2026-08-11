/**
 * pafcoverage.hpp: Get some stats from cigar pafs in order to see what kind of anchors they make for Cactus
 *                 (not done within mzgaf2paf so it can be used on paf's from other sources too)
 */


#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <ostream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <functional>
#include "paf.hpp"

using namespace std;

/** Report a failed write to standard output when the program exits.
 *
 * These tools stream their result to standard output and leave the caller to
 * redirect it.  A C++ stream reports a failed write -- a full disk, a quota, a
 * read-only mount -- by setting a flag that nothing is obliged to read, and the
 * flush that happens at exit discards any error, so without this a truncated
 * paf, gaf or bed and a complete one are indistinguishable and the tool still
 * exits successfully.  A short paf is still a valid paf.
 *
 * Call once, at the top of main.  Inline, so that the tools which do not link
 * pafcoverage.o can use it too.
 */
inline void check_stdout_at_exit() {
    static bool installed = false;
    if (installed) {
        return;
    }
    installed = true;
    // the flush has to happen inside the handler, because atexit runs before
    // the runtime flushes the streams, and _Exit is used because calling exit()
    // from an atexit handler is undefined
    atexit([]() {
        cout.flush();
        if (cout.fail() || cout.bad()) {
            cerr << "error: failed to write to standard output, so its contents are incomplete. "
                 << "Check the free space, the quota and the permissions on the file system "
                 << "holding it." << endl;
            _Exit(EXIT_FAILURE);
        }
    });
}

/** Flush and close a file we have written, returning false if anything failed.
 *
 * Closing is the last point at which buffered data reaches the operating
 * system, so a write that fails there is reported nowhere else;  leaving it to
 * the destructor discards the failure.
 */
inline bool close_output_file(ofstream& out) {
    if (not out.is_open()) {
        return true;
    }
    out.flush();
    out.close();
    return not (out.fail() or out.bad());
}

// map a sequence name to all its covered bases (not caring about depth, so just using bools)
typedef unordered_map<string, vector<bool>> CoverageMap;

/** update bases covered in the query sequence of a paf line
 */
void update_coverage_map(const string& paf_line, CoverageMap& coverage_map);

/** print some stats
 */
void print_coverage_summary(const CoverageMap& coverage_map, ostream& out);

/** print bed of coverage gaps
 */
void print_coverage_gaps_as_bed(const CoverageMap& coverage_map, ostream& out, int64_t min_gap_length);

