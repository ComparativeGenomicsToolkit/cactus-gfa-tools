/*
 * paf-one-to-one: make a PAF one-to-one in target space, haplotype by haplotype.
 *
 * gaffilter resolves two alignments that claim the same stretch of a *query*.  This is the other
 * direction: two stretches of one haplotype that are aligned to the same stretch of a *target*.
 * In a graph built from such a PAF the haplotype visits that target twice, which is right for a
 * real duplication, where the two copies are about as similar to the target as each other, and
 * wrong for a paralog that has been placed on the wrong copy, which stands out by its mismatches.
 *
 * The haplotype of a line is its query name up to the first '|'.  Wherever two or more lines of one
 * haplotype cover the same position of the same target, the lines whose local mismatch rate there
 * is within --margin of the best one are kept, and the stretch is cut out of the others, exactly,
 * by CIGAR.  The rest of every line is kept.
 *
 * Local mismatch rate of a line over a contested stretch [p, q) = X / (= + X), taken from
 *   1. the line's own CIGAR over target [p - pad, q + pad) (in whole --bin columns), if that holds
 *      at least --min-aligned aligned bases; otherwise
 *   2. the lines of its parent record (same query, gm:i: and gl:i:, as gaf2paf writes them) that
 *      start in the same --window of the query, if they do; otherwise
 *   3. the whole parent record.
 *
 * The input is read twice and so must be a file.  Lines come out in input order.
 */

#include <unistd.h>
#include <getopt.h>
#include <fcntl.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

#include "pafcoverage.hpp"

using namespace std;

static const char* PROG = "paf-one-to-one";

struct Params {
    double margin = 0.01;
    int64_t pad = 25000;
    int64_t min_aligned = 5000;
    int64_t window = 50000;
    int64_t bin = 1000;
    int64_t min_mapq = 0;
    int64_t min_block = 0;
    double min_ident = 0;
    vector<string> refs;
    string target_prefix;
    string log_path;
};

// one PAF line, as much of it as the decision needs
struct Row {
    int64_t offset;     // of the line in the file
    int64_t ts, te;     // target interval
    int64_t qs;         // query start
    uint32_t length;    // of the line, without its newline
    uint32_t hap;
    uint32_t target;
    uint32_t rec;       // parent record
    bool takes_part;
};

typedef pair<char, int64_t> Op;

static inline int64_t floor_div(int64_t a, int64_t b) {
    int64_t d = a / b;
    return (a % b != 0 && ((a < 0) != (b < 0))) ? d - 1 : d;
}

[[noreturn]] static void die(const string& msg) {
    cerr << "[" << PROG << "] error: " << msg << endl;
    exit(1);
}

// the tab-separated fields of a line, as (start, length) into it
static void split_fields(const string& line, vector<pair<size_t, size_t>>& fields) {
    fields.clear();
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == string::npos) {
            fields.push_back(make_pair(start, line.length() - start));
            break;
        }
        fields.push_back(make_pair(start, tab - start));
        start = tab + 1;
    }
}

static int64_t parse_int(const string& line, const pair<size_t, size_t>& f, int64_t line_no, const char* what) {
    if (f.second == 0) {
        die("line " + to_string(line_no) + ": empty " + what);
    }
    int64_t v = 0;
    size_t i = f.first, end = f.first + f.second;
    bool neg = false;
    if (line[i] == '-') {
        neg = true;
        ++i;
    }
    for (; i < end; ++i) {
        if (line[i] < '0' || line[i] > '9') {
            die("line " + to_string(line_no) + ": " + what + " is not a number: " + line.substr(f.first, f.second));
        }
        v = v * 10 + (line[i] - '0');
    }
    return neg ? -v : v;
}

// the field holding the given tag (e.g. "cg:Z:"), or -1
static int find_tag(const string& line, const vector<pair<size_t, size_t>>& fields, const char* tag) {
    size_t n = strlen(tag);
    for (size_t i = 12; i < fields.size(); ++i) {
        if (fields[i].second >= n && line.compare(fields[i].first, n, tag) == 0) {
            return (int)i;
        }
    }
    return -1;
}

// the cigar in [start, start + length) of line, in target order
static void parse_cigar(const string& line, size_t start, size_t length, vector<Op>& ops, int64_t line_no) {
    ops.clear();
    int64_t n = 0;
    bool digits = false;
    for (size_t i = start; i < start + length; ++i) {
        char c = line[i];
        if (c >= '0' && c <= '9') {
            n = n * 10 + (c - '0');
            digits = true;
        } else {
            if (!digits || (c != '=' && c != 'X' && c != 'M' && c != 'I' && c != 'D')) {
                die("line " + to_string(line_no) + ": cannot read the cigar at \"" + line.substr(i, 10) +
                    "\" (expected =, X, M, I and D operations)");
            }
            if (n > 0) {
                ops.push_back(make_pair(c, n));
            }
            n = 0;
            digits = false;
        }
    }
    if (digits) {
        die("line " + to_string(line_no) + ": the cigar ends in a number");
    }
}

static inline bool on_target(char c) { return c != 'I'; }
static inline bool on_query(char c) { return c != 'D'; }
static inline bool is_match(char c) { return c == '=' || c == 'X' || c == 'M'; }

// a line's mismatches and aligned bases in --bin columns of its target, as prefix sums
struct Prefix {
    int64_t b0;
    vector<int64_t> cx, ca;
};

// the part of a line's alignment over target [a, b)
struct Piece {
    vector<Op> ops;
    int64_t a = 0, b = 0;   // target interval, after any stripping
    int64_t q0 = 0;         // query offset of its start, in the cigar's direction
    bool any = false;
};

/* the sub-alignment of a line (target [ts, te), cigar ops) over target [a, b).  An insertion is the
   piece's if it lies strictly inside [a, b).  With strip, the indels left dangling at either end go,
   as they do when gaf2paf cuts a line, and a and b move in over any deletion that went */
static Piece cut_piece(const vector<Op>& ops, int64_t ts, int64_t a, int64_t b, bool strip) {
    Piece piece;
    piece.a = a;
    piece.b = b;
    int64_t tpos = ts, qoff = 0;
    bool have_q0 = false;
    for (const Op& op : ops) {
        if (op.first == 'I') {
            if (a < tpos && tpos < b) {
                if (!have_q0) {
                    piece.q0 = qoff;
                    have_q0 = true;
                }
                piece.ops.push_back(op);
            }
            qoff += op.second;
            continue;
        }
        int64_t s = tpos, e = tpos + op.second;
        int64_t lo = max(s, a), hi = min(e, b);
        if (lo < hi) {
            if (!have_q0) {
                piece.q0 = qoff + (op.first != 'D' ? lo - s : 0);
                have_q0 = true;
            }
            piece.ops.push_back(make_pair(op.first, hi - lo));
        }
        if (op.first != 'D') {
            qoff += op.second;
        }
        tpos = e;
        if (tpos >= b) {
            break;
        }
    }
    if (strip) {
        size_t first = 0;
        while (first < piece.ops.size() && !is_match(piece.ops[first].first)) {
            if (piece.ops[first].first == 'I') {
                piece.q0 += piece.ops[first].second;
            } else {
                piece.a += piece.ops[first].second;
            }
            ++first;
        }
        size_t last = piece.ops.size();
        while (last > first && !is_match(piece.ops[last - 1].first)) {
            if (piece.ops[last - 1].first == 'D') {
                piece.b -= piece.ops[last - 1].second;
            }
            --last;
        }
        piece.ops = vector<Op>(piece.ops.begin() + first, piece.ops.begin() + last);
    }
    // merge runs of one operation
    vector<Op> merged;
    for (const Op& op : piece.ops) {
        if (!merged.empty() && merged.back().first == op.first) {
            merged.back().second += op.second;
        } else {
            merged.push_back(op);
        }
    }
    piece.ops.swap(merged);
    for (const Op& op : piece.ops) {
        piece.any = piece.any || is_match(op.first);
    }
    return piece;
}

static int64_t query_length(const vector<Op>& ops) {
    int64_t n = 0;
    for (const Op& op : ops) {
        if (on_query(op.first)) {
            n += op.second;
        }
    }
    return n;
}

struct Trim {
    uint32_t row;       // the line cut
    uint32_t keeper;    // the best line there
    int64_t p, q;       // target interval
    double mm, best;
    int src;            // which of the three rates mm is
};

static void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <paf> > output.paf" << endl
         << "Make a PAF one-to-one in target space, one haplotype at a time.  A haplotype is a query name up to" << endl
         << "its first '|'.  Wherever two or more lines of one haplotype cover the same position of the same target," << endl
         << "the lines whose local mismatch rate (X / (= + X)) there is less than --margin above the best are kept, and" << endl
         << "the stretch is cut out of the others by CIGAR.  The rest of each line is kept.  The local rate is the" << endl
         << "line's own, over the stretch and --pad either side; where that is under --min-aligned bases, it is" << endl
         << "that of its parent record (same query, gm:i: and gl:i:) in the same --window of the query, or else" << endl
         << "of the whole record.  The input is read twice, so it must be a file." << endl
         << endl
         << "options: " << endl
         << "    -g, --margin F               Keep every line whose mismatch rate is less than this above the best [0.01]" << endl
         << "    -r, --ref NAME               Never cut the lines of queries NAME|* or id=NAME|* (may be repeated)" << endl
         << "    -t, --target-prefix PREFIX   Only lines whose target name starts with PREFIX take part" << endl
         << "    -q, --min-mapq N             Only lines with MAPQ >= N take part [0]" << endl
         << "    -b, --min-block-length N     Only lines whose record's block length (gl:i:) is >= N take part," << endl
         << "                                 or whose query is no longer than N [0]" << endl
         << "    -i, --min-identity F         Only lines with identity >= F take part (the lower of the line's" << endl
         << "                                 and its record's, gi:f:) [0]" << endl
         << "    -p, --pad N                  Target bases either side of a stretch that its rate is taken over [25000]" << endl
         << "    -a, --min-aligned N          Fewest aligned bases a rate is taken from [5000]" << endl
         << "    -w, --window N               Query window of the record-level rate [50000]" << endl
         << "    -B, --bin N                  Target column width the line's own rate is counted in [1000]" << endl
         << "    -l, --log FILE               Write one line per stretch cut to FILE" << endl
         << endl
         << "Lines that do not take part (-t, -q, -b, -i) are written as they are: they neither cut other lines" << endl
         << "nor are cut.  -q, -b and -i are those of the line filter that follows in cactus, so that a line it" << endl
         << "is about to drop does not cost another line its alignment." << endl;
}

int main(int argc, char** argv) {
    check_stdout_at_exit();

    Params P;
    int c;
    optind = 1;
    while (true) {
        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"margin", required_argument, 0, 'g'},
            {"ref", required_argument, 0, 'r'},
            {"target-prefix", required_argument, 0, 't'},
            {"min-mapq", required_argument, 0, 'q'},
            {"min-block-length", required_argument, 0, 'b'},
            {"min-identity", required_argument, 0, 'i'},
            {"pad", required_argument, 0, 'p'},
            {"min-aligned", required_argument, 0, 'a'},
            {"window", required_argument, 0, 'w'},
            {"bin", required_argument, 0, 'B'},
            {"log", required_argument, 0, 'l'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "hg:r:t:q:b:i:p:a:w:B:l:", long_options, &option_index);
        if (c == -1) {
            break;
        }
        try {
            switch (c) {
            case 'g': P.margin = stod(optarg); break;
            case 'r': P.refs.push_back(optarg); break;
            case 't': P.target_prefix = optarg; break;
            case 'q': P.min_mapq = stol(optarg); break;
            case 'b': P.min_block = stol(optarg); break;
            case 'i': P.min_ident = stod(optarg); break;
            case 'p': P.pad = stol(optarg); break;
            case 'a': P.min_aligned = stol(optarg); break;
            case 'w': P.window = stol(optarg); break;
            case 'B': P.bin = stol(optarg); break;
            case 'l': P.log_path = optarg; break;
            case 'h':
            case '?':
                help(argv);
                exit(1);
            default:
                abort();
            }
        } catch (const exception&) {
            die(string("cannot read the value of -") + (char)c + ": " + optarg);
        }
    }
    if (optind != argc - 1) {
        help(argv);
        return 1;
    }
    if (P.margin < 0 || P.pad < 0 || P.min_aligned < 1 || P.window < 1 || P.bin < 1) {
        die("--margin and --pad cannot be negative, and --min-aligned, --window and --bin must be at least 1");
    }
    string paf_path = argv[optind];
    if (paf_path == "-") {
        die("the input is read twice and so cannot be standard input");
    }

    // ---- pass 1: every line's coordinates, and each record's mismatch counts ----
    ifstream paf(paf_path);
    if (!paf) {
        die("unable to open paf: " + paf_path);
    }
    vector<Row> rows;
    unordered_map<string, uint32_t> hap_ids, target_ids, rec_ids;
    vector<string> hap_names;
    vector<int64_t> rec_x, rec_al;
    unordered_map<uint64_t, pair<int64_t, int64_t>> win_sums;   // (record, query window) -> X, aligned
    string line, key;
    vector<pair<size_t, size_t>> fields;
    int64_t offset = 0, line_no = 0, lines_out_of_scope = 0;
    vector<int64_t> hap_m_columns;   // M columns of each haplotype's lines that take part
    while (getline(paf, line)) {
        ++line_no;
        int64_t line_offset = offset;
        offset += line.length() + 1;
        Row row;
        row.offset = line_offset;
        row.length = line.length();
        row.takes_part = false;
        row.hap = row.target = row.rec = 0;
        row.ts = row.te = row.qs = 0;
        if (line.empty() || line[0] == '#') {
            rows.push_back(row);
            continue;
        }
        split_fields(line, fields);
        if (fields.size() < 12) {
            die("line " + to_string(line_no) + " has " + to_string(fields.size()) + " columns, and a PAF line has at least 12");
        }
        if (line.length() > 0xFFFFFFFFull) {
            die("line " + to_string(line_no) + " is longer than 4 GiB");
        }
        int64_t query_len = parse_int(line, fields[1], line_no, "query length");
        row.qs = parse_int(line, fields[2], line_no, "query start");
        row.ts = parse_int(line, fields[7], line_no, "target start");
        row.te = parse_int(line, fields[8], line_no, "target end");
        int64_t matches = parse_int(line, fields[9], line_no, "matches");
        int64_t block = parse_int(line, fields[10], line_no, "block length");
        int64_t mapq = parse_int(line, fields[11], line_no, "MAPQ");

        int cg = find_tag(line, fields, "cg:Z:");
        int gm = find_tag(line, fields, "gm:i:");
        int gl = find_tag(line, fields, "gl:i:");
        int gi = find_tag(line, fields, "gi:f:");

        // the haplotype: the query name up to its first '|'
        size_t bar = line.find('|', fields[0].first);
        size_t hap_len = (bar == string::npos || bar >= fields[0].first + fields[0].second) ? fields[0].second : bar - fields[0].first;
        key.assign(line, fields[0].first, hap_len);
        auto hap_it = hap_ids.find(key);
        if (hap_it == hap_ids.end()) {
            hap_it = hap_ids.insert(make_pair(key, (uint32_t)hap_names.size())).first;
            hap_names.push_back(key);
            hap_m_columns.push_back(0);
        }
        row.hap = hap_it->second;
        key.assign(line, fields[5].first, fields[5].second);
        row.target = target_ids.insert(make_pair(key, (uint32_t)target_ids.size())).first->second;

        // the parent record: gaf2paf gives every line of a GAF record that record's matches and
        // block length, and nothing else names it
        key.assign(line, fields[0].first, fields[0].second);
        key += '\t';
        if (gm >= 0) {
            key.append(line, fields[gm].first + 5, fields[gm].second - 5);
        }
        key += '\t';
        if (gl >= 0) {
            key.append(line, fields[gl].first + 5, fields[gl].second - 5);
        }
        auto rec_it = rec_ids.find(key);
        if (rec_it == rec_ids.end()) {
            rec_it = rec_ids.insert(make_pair(key, (uint32_t)rec_x.size())).first;
            rec_x.push_back(0);
            rec_al.push_back(0);
        }
        row.rec = rec_it->second;

        // does it take part
        bool in_scope = cg >= 0 && row.te > row.ts;
        if (in_scope && !P.target_prefix.empty() &&
            line.compare(fields[5].first, min(fields[5].second, P.target_prefix.length()), P.target_prefix) != 0) {
            in_scope = false;
        }
        if (in_scope) {
            // the line filter cactus applies afterwards (cactus_graphmap.filter_paf)
            double ident = (double)matches / ((double)block + 0.00000001);
            if (gi >= 0) {
                ident = min(ident, atof(line.c_str() + fields[gi].first + 5));
            }
            int64_t record_block = gl >= 0 ? parse_int(line, make_pair(fields[gl].first + 5, fields[gl].second - 5), line_no, "gl:i:") : -1;
            in_scope = mapq >= P.min_mapq && ident >= P.min_ident &&
                (record_block < 0 || query_len <= P.min_block || record_block >= P.min_block);
        }
        row.takes_part = in_scope;
        if (!in_scope) {
            ++lines_out_of_scope;
        }

        // the record's counts are of all of its lines, whether or not they take part
        if (cg >= 0) {
            int64_t x = 0, al = 0, n = 0;
            for (size_t i = fields[cg].first + 5; i < fields[cg].first + fields[cg].second; ++i) {
                char ch = line[i];
                if (ch >= '0' && ch <= '9') {
                    n = n * 10 + (ch - '0');
                } else {
                    if (ch == '=' || ch == 'X' || ch == 'M') {
                        al += n;
                        if (ch == 'X') {
                            x += n;
                        } else if (ch == 'M' && in_scope) {
                            hap_m_columns[row.hap] += n;
                        }
                    }
                    n = 0;
                }
            }
            rec_x[row.rec] += x;
            rec_al[row.rec] += al;
            auto& win = win_sums[((uint64_t)row.rec << 32) | (uint64_t)(row.qs / P.window)];
            win.first += x;
            win.second += al;
        }
        if (rows.size() >= 0xFFFFFFFFull) {
            die("more than 4 billion lines");
        }
        rows.push_back(row);
    }
    paf.close();

    vector<bool> hap_is_ref(hap_names.size(), false);
    // M columns that could have been judged: not those of a reference, whose lines are never cut
    // (rgfa2paf writes the reference's own as one M each)
    int64_t m_columns = 0;
    for (size_t h = 0; h < hap_names.size(); ++h) {
        for (const string& ref : P.refs) {
            if (hap_names[h] == ref || hap_names[h] == "id=" + ref) {
                hap_is_ref[h] = true;
            }
        }
        if (!hap_is_ref[h]) {
            m_columns += hap_m_columns[h];
        }
    }

    // ---- the contested stretches ----
    int fd = open(paf_path.c_str(), O_RDONLY);
    if (fd < 0) {
        die("unable to open paf: " + paf_path);
    }
    auto read_line = [&](uint32_t i, string& out) {
        out.resize(rows[i].length);
        size_t done = 0;
        while (done < out.size()) {
            ssize_t got = pread(fd, &out[done], out.size() - done, rows[i].offset + done);
            if (got <= 0) {
                die("unable to read " + paf_path + " again: it must be a file that does not change");
            }
            done += got;
        }
    };
    auto load_ops = [&](uint32_t i, string& buffer, vector<pair<size_t, size_t>>& flds, vector<Op>& ops) {
        read_line(i, buffer);
        split_fields(buffer, flds);
        int cg = find_tag(buffer, flds, "cg:Z:");
        assert(cg >= 0);
        parse_cigar(buffer, flds[cg].first + 5, flds[cg].second - 5, ops, (int64_t)i + 1);
    };

    unordered_map<uint32_t, Prefix> prefixes;
    string buffer;
    vector<Op> ops;
    auto prefix_of = [&](uint32_t i) -> const Prefix& {
        auto it = prefixes.find(i);
        if (it != prefixes.end()) {
            return it->second;
        }
        load_ops(i, buffer, fields, ops);
        const Row& row = rows[i];
        Prefix& prefix = prefixes[i];
        prefix.b0 = row.ts / P.bin;
        size_t n_bins = (row.te - 1) / P.bin - prefix.b0 + 1;
        vector<int64_t> x(n_bins, 0), al(n_bins, 0);
        int64_t tpos = row.ts;
        for (const Op& op : ops) {
            if (op.first == 'I') {
                continue;
            }
            if (op.first == 'D') {
                tpos += op.second;
                continue;
            }
            int64_t e = tpos + op.second;
            if (e > row.te) {
                die("line " + to_string((int64_t)i + 1) + ": the cigar runs past the target end");
            }
            while (tpos < e) {
                int64_t k = tpos / P.bin;
                int64_t s = min(e, (k + 1) * P.bin);
                al[k - prefix.b0] += s - tpos;
                if (op.first == 'X') {
                    x[k - prefix.b0] += s - tpos;
                }
                tpos = s;
            }
        }
        prefix.cx.assign(n_bins + 1, 0);
        prefix.ca.assign(n_bins + 1, 0);
        for (size_t k = 0; k < n_bins; ++k) {
            prefix.cx[k + 1] = prefix.cx[k] + x[k];
            prefix.ca[k + 1] = prefix.ca[k] + al[k];
        }
        return prefix;
    };
    // the local mismatch rate of line i over target [p, q), and which of the three it is
    auto local_rate = [&](uint32_t i, int64_t p, int64_t q) -> pair<double, int> {
        const Prefix& prefix = prefix_of(i);
        int64_t n = (int64_t)prefix.ca.size() - 1;
        int64_t lo = max((int64_t)0, floor_div(p - P.pad, P.bin) - prefix.b0);
        int64_t hi = min(n, floor_div(q + P.pad - 1, P.bin) - prefix.b0 + 1);
        if (hi > lo && prefix.ca[hi] - prefix.ca[lo] >= P.min_aligned) {
            return make_pair((double)(prefix.cx[hi] - prefix.cx[lo]) / (double)(prefix.ca[hi] - prefix.ca[lo]), 1);
        }
        const Row& row = rows[i];
        auto win = win_sums.find(((uint64_t)row.rec << 32) | (uint64_t)(row.qs / P.window));
        if (win != win_sums.end() && win->second.second >= P.min_aligned) {
            return make_pair((double)win->second.first / (double)win->second.second, 2);
        }
        return make_pair(rec_al[row.rec] ? (double)rec_x[row.rec] / (double)rec_al[row.rec] : 0.0, 3);
    };

    // the lines that take part, grouped by haplotype and target
    vector<uint32_t> order;
    for (size_t i = 0; i < rows.size(); ++i) {
        if (rows[i].takes_part && !hap_is_ref[rows[i].hap]) {
            order.push_back(i);
        }
    }
    sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        const Row& ra = rows[a];
        const Row& rb = rows[b];
        if (ra.hap != rb.hap) return ra.hap < rb.hap;
        if (ra.target != rb.target) return ra.target < rb.target;
        return a < b;
    });

    vector<Trim> trims;
    int64_t contested_bp = 0, contested_bp_all_kept = 0, groups = 0;
    int64_t src_bp[4] = {0, 0, 0, 0};
    struct Event {
        int64_t pos;
        int start;      // ends sort before starts at one position
        uint32_t row;
    };
    vector<Event> events;
    vector<uint32_t> active;
    vector<pair<double, int>> rates;
    for (size_t g0 = 0; g0 < order.size();) {
        size_t g1 = g0 + 1;
        while (g1 < order.size() && rows[order[g1]].hap == rows[order[g0]].hap && rows[order[g1]].target == rows[order[g0]].target) {
            ++g1;
        }
        if (g1 - g0 >= 2) {
            ++groups;
            events.clear();
            for (size_t k = g0; k < g1; ++k) {
                events.push_back({rows[order[k]].ts, 1, order[k]});
                events.push_back({rows[order[k]].te, 0, order[k]});
            }
            sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
                if (a.pos != b.pos) return a.pos < b.pos;
                if (a.start != b.start) return a.start < b.start;
                return a.row < b.row;
            });
            active.clear();
            int64_t prev = 0;
            for (const Event& ev : events) {
                if (active.size() >= 2 && ev.pos > prev) {
                    int64_t p = prev, q = ev.pos;
                    rates.clear();
                    double best = 0;
                    uint32_t keeper = 0;
                    for (size_t k = 0; k < active.size(); ++k) {
                        rates.push_back(local_rate(active[k], p, q));
                        if (k == 0 || rates[k].first < best || (rates[k].first == best && active[k] < keeper)) {
                            best = rates[k].first;
                            keeper = active[k];
                        }
                    }
                    contested_bp += q - p;
                    bool all_kept = true;
                    for (size_t k = 0; k < active.size(); ++k) {
                        if (!(rates[k].first < best + P.margin)) {
                            all_kept = false;
                            trims.push_back({active[k], keeper, p, q, rates[k].first, best, rates[k].second});
                            src_bp[rates[k].second] += q - p;
                        }
                    }
                    if (all_kept) {
                        contested_bp_all_kept += q - p;
                    }
                }
                if (ev.start) {
                    active.push_back(ev.row);
                } else {
                    active.erase(find(active.begin(), active.end(), ev.row));
                }
                prev = ev.pos;
            }
            // a line is in one group only, so its columns are not needed again
            prefixes.clear();
        }
        g0 = g1;
    }

    // each cut line's stretches, in target order (a line's stretches never overlap: the sweep
    // visits every position of a target once)
    sort(trims.begin(), trims.end(), [](const Trim& a, const Trim& b) {
        if (a.row != b.row) return a.row < b.row;
        return a.p < b.p;
    });

    // ---- pass 2: write the lines, cutting the ones that lost a stretch ----
    ofstream log;
    if (!P.log_path.empty()) {
        log.open(P.log_path);
        if (!log) {
            die("unable to open log: " + P.log_path);
        }
        log << "haplotype\tquery\tq_start\tq_end\tstrand\ttarget\tt_start\tt_end\tmismatch_rate\trate_source\tbest_rate\tkeeper_query\tkeeper_q_start\n";
    }
    paf.open(paf_path);
    if (!paf) {
        die("unable to open paf: " + paf_path);
    }
    int64_t rows_cut = 0, rows_removed = 0, pieces_out = 0, pieces_dropped = 0, bp_cut = 0, lines_out = 0;
    vector<int64_t> hap_bp(hap_names.size(), 0);
    size_t t0 = 0;
    string keeper_line, out_line;
    vector<pair<size_t, size_t>> keeper_fields;
    vector<pair<int64_t, int64_t>> cut, keep;
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!getline(paf, line) || line.length() != rows[i].length) {
            die(paf_path + " changed while it was being read");
        }
        size_t t1 = t0;
        while (t1 < trims.size() && trims[t1].row == i) {
            ++t1;
        }
        if (t1 == t0) {
            cout << line << "\n";
            ++lines_out;
            continue;
        }
        const Row& row = rows[i];
        split_fields(line, fields);
        int cg = find_tag(line, fields, "cg:Z:");
        parse_cigar(line, fields[cg].first + 5, fields[cg].second - 5, ops, (int64_t)i + 1);
        int64_t qs = row.qs, qe = parse_int(line, fields[3], (int64_t)i + 1, "query end");
        bool reverse = line[fields[4].first] == '-';
        {
            int64_t tlen = 0;
            for (const Op& op : ops) {
                tlen += on_target(op.first) ? op.second : 0;
            }
            if (tlen != row.te - row.ts || query_length(ops) != qe - qs) {
                die("line " + to_string((int64_t)i + 1) + ": the cigar does not span the line's query and target intervals");
            }
        }
        ++rows_cut;
        // merge the stretches, and log them
        cut.clear();
        for (size_t t = t0; t < t1; ++t) {
            const Trim& trim = trims[t];
            if (!cut.empty() && trim.p <= cut.back().second) {
                cut.back().second = max(cut.back().second, trim.q);
            } else {
                cut.push_back(make_pair(trim.p, trim.q));
            }
            if (log.is_open()) {
                Piece piece = cut_piece(ops, row.ts, trim.p, trim.q, false);
                int64_t ql = query_length(piece.ops);
                int64_t qa = reverse ? qe - piece.q0 - ql : qs + piece.q0;
                read_line(trim.keeper, keeper_line);
                split_fields(keeper_line, keeper_fields);
                char mm[32], best[32];
                snprintf(mm, sizeof(mm), "%.4f", trim.mm);
                snprintf(best, sizeof(best), "%.4f", trim.best);
                log << hap_names[row.hap] << "\t" << line.substr(fields[0].first, fields[0].second) << "\t" << qa << "\t" << qa + ql << "\t"
                    << (reverse ? '-' : '+') << "\t" << line.substr(fields[5].first, fields[5].second) << "\t" << trim.p << "\t" << trim.q << "\t"
                    << mm << "\t" << trim.src << "\t" << best << "\t"
                    << keeper_line.substr(keeper_fields[0].first, keeper_fields[0].second) << "\t"
                    << keeper_line.substr(keeper_fields[2].first, keeper_fields[2].second) << "\n";
            }
        }
        for (const auto& iv : cut) {
            bp_cut += iv.second - iv.first;
            hap_bp[row.hap] += iv.second - iv.first;
        }
        // what is left
        keep.clear();
        int64_t last = row.ts;
        for (const auto& iv : cut) {
            if (iv.first > last) {
                keep.push_back(make_pair(last, iv.first));
            }
            last = iv.second;
        }
        if (row.te > last) {
            keep.push_back(make_pair(last, row.te));
        }
        int64_t written = 0;
        for (const auto& iv : keep) {
            Piece piece = cut_piece(ops, row.ts, iv.first, iv.second, true);
            if (!piece.any) {
                ++pieces_dropped;
                continue;
            }
            int64_t ql = query_length(piece.ops), tl = 0, matches = 0, block = 0;
            for (const Op& op : piece.ops) {
                tl += on_target(op.first) ? op.second : 0;
                matches += (op.first == '=' || op.first == 'M') ? op.second : 0;
                block += op.second;
            }
            assert(tl == piece.b - piece.a);
            int64_t qa = reverse ? qe - piece.q0 - ql : qs + piece.q0;
            out_line.clear();
            for (size_t f = 0; f < fields.size(); ++f) {
                if (f > 0) {
                    out_line += '\t';
                }
                switch (f) {
                case 2: out_line += to_string(qa); break;
                case 3: out_line += to_string(qa + ql); break;
                case 7: out_line += to_string(piece.a); break;
                case 8: out_line += to_string(piece.b); break;
                case 9: out_line += to_string(matches); break;
                case 10: out_line += to_string(block); break;
                default:
                    if ((int)f == cg) {
                        out_line += "cg:Z:";
                        for (const Op& op : piece.ops) {
                            out_line += to_string(op.second);
                            out_line += op.first;
                        }
                    } else {
                        out_line.append(line, fields[f].first, fields[f].second);
                    }
                }
            }
            cout << out_line << "\n";
            ++lines_out;
            ++pieces_out;
            ++written;
        }
        if (written == 0) {
            ++rows_removed;
        }
        t0 = t1;
    }
    close(fd);
    if (log.is_open()) {
        log.close();
        if (log.fail()) {
            die("failed to write log: " + P.log_path);
        }
    }

    int64_t haps_cut = 0;
    for (int64_t bp : hap_bp) {
        haps_cut += bp > 0;
    }
    cerr << "[" << PROG << "]: " << rows.size() << " lines in (" << lines_out_of_scope << " not taking part), " << lines_out << " out; "
         << groups << " haplotype x target groups of two or more lines" << endl
         << "[" << PROG << "]: " << contested_bp << " bp covered more than once by a haplotype, " << contested_bp_all_kept
         << " bp of it kept whole (every line within " << P.margin << " of the best)" << endl
         << "[" << PROG << "]: cut " << bp_cut << " bp out of " << rows_cut << " lines of " << haps_cut << " haplotypes ("
         << rows_removed << " lines removed entirely, " << pieces_out << " pieces written, " << pieces_dropped
         << " pieces with no aligned base left dropped)" << endl
         << "[" << PROG << "]: cut stretches by rate source: line's own " << src_bp[1] << " bp, record window " << src_bp[2]
         << " bp, whole record " << src_bp[3] << " bp" << endl;
    if (m_columns > 0) {
        cerr << "[" << PROG << "] warning: " << m_columns << " aligned columns of the lines taking part are M operations, which do not "
             << "say whether they match: they count as matches, so lines without =/X operations cannot be judged" << endl;
    }
    return 0;
}
