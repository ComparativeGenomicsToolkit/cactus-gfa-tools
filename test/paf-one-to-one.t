#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH

plan tests 16

# one line per call: query q_start q_end strand t_start t_end mapq cigar [record-block-length]
# (query length 1000000, target T of length 1000, matches and block length counted from the cigar)
paf_line() {
    python3 -c '
import sys, re
q, qs, qe, strand, ts, te, mapq, cg = sys.argv[1:9]
gl = sys.argv[9] if len(sys.argv) > 9 else "500000"
ops = [(int(n), o) for n, o in re.findall(r"(\d+)([=XID])", cg)]
assert sum(n for n, o in ops if o != "I") == int(te) - int(ts), "target span"
assert sum(n for n, o in ops if o != "D") == int(qe) - int(qs), "query span"
print("\t".join([q, "1000000", qs, qe, strand, "T", "1000", ts, te, str(sum(n for n, o in ops if o == "=")),
                 str(sum(n for n, o in ops)), mapq, "tp:A:P", "gm:i:" + gl, "gl:i:" + gl, "gi:f:1", "cg:Z:" + cg]))' "$@"
}
# every test counts mismatches base by base, over the contested stretch alone
O="-a 1 -p 0 -B 1"

# a clean line and a second copy of the same haplotype, 10% off over the 50 bp they share
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.paf
paf_line "id=S1|c2" 500 650 + 150 300 60 20=5X125= >> o2o.paf
paf-one-to-one o2o.paf $O -l o2o.log > o2o.out 2> o2o.err
is $? 0 "runs"
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.expect
paf_line "id=S1|c2" 550 650 + 200 300 60 100= >> o2o.expect
diff o2o.out o2o.expect
is $? 0 "the worse copy loses the stretch the two share, and keeps the rest"
is "$(tail -n +2 o2o.log | cut -f 1-8,10)" "$(printf 'id=S1\tid=S1|c2\t500\t550\t+\tT\t150\t200\t1')" "the cut is logged with its query and target intervals"

# the same on the reverse strand: the cigar runs along the target, the query backwards
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.paf
paf_line "id=S1|c2" 500 650 - 150 300 60 20=5X125= >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.expect
paf_line "id=S1|c2" 500 600 - 200 300 60 100= >> o2o.expect
diff o2o.out o2o.expect
is $? 0 "a reverse-strand line is cut at the right end of its query"

# a deletion and an insertion left dangling by the cut go, and the coordinates follow
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.paf
paf_line "id=S1|c2" 500 640 + 150 300 60 20=5X25=10D90= >> o2o.paf
paf_line "id=S2|c1" 0 100 + 100 200 60 100= >> o2o.paf
paf_line "id=S2|c2" 500 657 + 150 300 60 20=5X25=7I100= >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.expect
paf_line "id=S1|c2" 550 640 + 210 300 60 90= >> o2o.expect
paf_line "id=S2|c1" 0 100 + 100 200 60 100= >> o2o.expect
paf_line "id=S2|c2" 557 657 + 200 300 60 100= >> o2o.expect
diff o2o.out o2o.expect
is $? 0 "indels dangling at a cut are stripped"

# a line cut in the middle comes out as two
paf_line "id=S1|c1" 0 50 + 200 250 60 50= > o2o.paf
paf_line "id=S1|c2" 500 700 + 100 300 60 100=10X90= >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
paf_line "id=S1|c1" 0 50 + 200 250 60 50= > o2o.expect
paf_line "id=S1|c2" 500 600 + 100 200 60 100= >> o2o.expect
paf_line "id=S1|c2" 650 700 + 250 300 60 50= >> o2o.expect
diff o2o.out o2o.expect
is $? 0 "a line cut in its middle comes out as two"

# a real duplication: both copies about as close to the target as each other (0.5% apart; a line
# is cut when its rate is the margin or more above the best)
paf_line "id=S1|c1" 0 200 + 100 300 60 200= > o2o.paf
paf_line "id=S1|c2" 500 700 + 100 300 60 100=1X99= >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
diff o2o.out o2o.paf
is $? 0 "two copies within the margin are both kept"
paf-one-to-one o2o.paf $O -g 0.005 > o2o.out 2> o2o.err
is "$(cat o2o.out | wc -l)" 1 "and not when the margin is the difference between them"

# other haplotypes, and other targets, are none of a line's business
paf_line "id=S1|c1" 0 100 + 100 200 60 100= > o2o.paf
paf_line "id=S2|c1" 500 600 + 100 200 60 50=10X40= >> o2o.paf
paf_line "id=S1|c2" 500 600 + 100 200 60 50=10X40= | awk 'BEGIN{FS=OFS="\t"} {$6 = "U"; print}' >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
diff o2o.out o2o.paf
is $? 0 "lines of different haplotypes or on different targets do not compete"

# the reference is never cut
paf_line "id=REF|c1" 0 100 + 100 200 60 100= > o2o.paf
paf_line "id=REF|c2" 500 600 + 100 200 60 50=10X40= >> o2o.paf
paf-one-to-one o2o.paf $O -r REF > o2o.out 2> o2o.err
diff o2o.out o2o.paf
is $? 0 "-r exempts a haplotype"

# a line the line filter is about to drop does not take part
paf_line "id=S1|c1" 0 100 + 100 200 3 100= > o2o.paf
paf_line "id=S1|c2" 500 600 + 100 200 60 50=10X40= >> o2o.paf
paf-one-to-one o2o.paf $O > o2o.out 2> o2o.err
is "$(cat o2o.out | wc -l)" 1 "a low-MAPQ line cuts another by default"
paf-one-to-one o2o.paf $O -q 5 > o2o.out 2> o2o.err
diff o2o.out o2o.paf
is $? 0 "but not with -q above its MAPQ"
paf_line "id=S1|c1" 0 100 + 100 200 60 100= 2000 > o2o.paf
paf_line "id=S1|c2" 500 600 + 100 200 60 50=10X40= >> o2o.paf
paf-one-to-one o2o.paf $O -b 250000 > o2o.out 2> o2o.err
diff o2o.out o2o.paf
is $? 0 "nor a line of a record shorter than -b"

# M says nothing about mismatches: said once, and not about a reference's lines
paf_line "id=S1|c1" 0 100 + 100 200 60 100= | sed 's/100=$/100M/' > o2o.paf
paf_line "id=S1|c2" 500 600 + 100 200 60 50=10X40= >> o2o.paf
paf-one-to-one o2o.paf $O > /dev/null 2> o2o.err
is "$(grep -c 'warning: 100 aligned columns' o2o.err)" 1 "M operations are warned about"
sed -i 's/id=S1/id=REF/' o2o.paf
paf-one-to-one o2o.paf $O -r REF > /dev/null 2> o2o.err
is "$(grep -c warning o2o.err)" 0 "but not on the reference's lines"

cat o2o.paf | paf-one-to-one - $O > /dev/null 2> o2o.err
is $? 1 "standard input is refused"

rm -f o2o.paf o2o.out o2o.err o2o.log o2o.expect
