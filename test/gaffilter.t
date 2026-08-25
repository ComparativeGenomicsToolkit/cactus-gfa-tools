#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=..:$PATH
PATH=../bin:$PATH

plan tests 20

# all fixtures are written here rather than checked in: they are a few lines each, and having the
# expected geometry visible next to the assertion is the point of this file.  no external tools.
WORK=gaffilter.tmp
rm -rf $WORK && mkdir -p $WORK

# a PAF line.  cols: qname qlen qstart qend strand tname tlen tstart tend matches blocklen mapq tags
paf_line() { # qname qstart qend tname mapq rc [extra]
    printf '%s\t1000000\t%s\t%s\t+\t%s\t2000000\t0\t%s\t%s\t%s\t%s\trc:Z:%s\tcg:Z:%s=\n' \
        "$1" "$2" "$3" "$4" "$(($3-$2))" "$(($3-$2))" "$(($3-$2))" "$5" "$6" "$(($3-$2))"
}
# a GAF line. cols: qname qlen qstart qend strand path plen pstart pend matches blocklen mapq tags
gaf_line() { # qname qstart qend path mapq rc
    printf '%s\t1000000\t%s\t%s\t+\t>%s\t2000000\t0\t%s\t%s\t%s\t%s\trc:Z:%s\n' \
        "$1" "$2" "$3" "$4" "$(($3-$2))" "$(($3-$2))" "$(($3-$2))" "$5" "$6"
}
count() { wc -l < "$1" | tr -d ' '; }

########################################################################
# the rc:Z: exemption, and what -x changes about it
########################################################################

# two records on the same query, decisively different block lengths, on DIFFERENT ref contigs
paf_line q1     0 100000 t1 10 chrA  > $WORK/cross.paf
paf_line q1 50000  52000 t2 10 chrB >> $WORK/cross.paf

gaffilter $WORK/cross.paf -p -r 5 -m 0 -q 5 > $WORK/cross.x0 2>/dev/null
is "$(count $WORK/cross.x0)" 2 "without -x, records on different ref contigs are never compared"

gaffilter $WORK/cross.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/cross.x60 2>/dev/null
is "$(count $WORK/cross.x60)" 1 "with -x, an off-contig record can remove an overlapping one"
is "$(awk '{print $3}' $WORK/cross.x60)" 0 "the record that survives is the dominant one"

# -x is a ceiling on the LOSER's mapq, and the bound is strict
paf_line q2     0 100000 t1 60 chrA  > $WORK/bound.paf
paf_line q2 50000  52000 t2 59 chrB >> $WORK/bound.paf
gaffilter $WORK/bound.paf -p -r 5 -m 0 -q 5 -x 59 > $WORK/bound.59 2>/dev/null
is "$(count $WORK/bound.59)" 2 "a loser at exactly -x N is not removable (strict <)"
gaffilter $WORK/bound.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/bound.60 2>/dev/null
is "$(count $WORK/bound.60)" 1 "a loser one below -x N is removable"

########################################################################
# a cross-contig tie must never delete either side, at ANY ratio.
# dominates() is antisymmetric only for ratio > 1, so this is the case that
# regressed: at -r <= 1 the mapq arm returned true in both directions and the
# two records deleted each other.
########################################################################

paf_line q3 0 1000 t1 30 chrA  > $WORK/tie.paf
paf_line q3 0 1000 t2 30 chrB >> $WORK/tie.paf

gaffilter $WORK/tie.paf -p -r 5 -m 0 -x 60 > $WORK/tie.r5 2>/dev/null
is "$(count $WORK/tie.r5)" 2 "cross-contig tie deletes neither side at -r 5"
gaffilter $WORK/tie.paf -p -r 1 -m 0 -x 60 > $WORK/tie.r1 2>/dev/null
is "$(count $WORK/tie.r1)" 2 "cross-contig tie deletes neither side at -r 1"
gaffilter $WORK/tie.paf -p -r 0.9 -m 0 -x 60 > $WORK/tie.r09 2>/dev/null
is "$(count $WORK/tie.r09)" 2 "cross-contig tie deletes neither side at -r 0.9"

########################################################################
# a missing MAPQ parses to -1, which is below every threshold.  unknown
# confidence is not low confidence, so it must not be removable by -x.
########################################################################

gaf_line q4     0 100000 s1 '*' chrA  > $WORK/mq.gaf
gaf_line q4 50000  52000 s2 60  chrB >> $WORK/mq.gaf
gaffilter $WORK/mq.gaf -r 5 -m 0 -x 1 > $WORK/mq.out 2>/dev/null
is "$(count $WORK/mq.out)" 2 "a record with no MAPQ is not removable by -x"

########################################################################
# -P: this tool has no notion of a reference, but its caller exempts one from
# every other filter, so without -P the overlap pass quietly undoes that.
########################################################################

# the reference's own on-target alignment at mapq 0, against an off-target one at 60
paf_line 'id=REF|chr1' 100000 200000 t1  0 chr1  > $WORK/ref.paf
paf_line 'id=REF|chr1'  90000 900000 t2 60 chr5 >> $WORK/ref.paf

gaffilter $WORK/ref.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/ref.noP 2>/dev/null
is "$(count $WORK/ref.noP)" 1 "without -P the reference loses its on-target alignment"
gaffilter $WORK/ref.paf -p -r 5 -m 0 -q 5 -x 60 -P 'id=REF|' > $WORK/ref.P 2>/dev/null
is "$(count $WORK/ref.P)" 2 "with -P the reference keeps both"

# -P must not spare anything else.  note records only ever compete within one query name (the
# interval trees are keyed by it), so protecting a prefix protects whole queries and nothing more
paf_line 'id=REF|chr1' 100000 200000 t1  0 chr1  > $WORK/mix.paf
paf_line 'id=REF|chr1'  90000 900000 t2 60 chr5 >> $WORK/mix.paf
paf_line 'other|ctg'        0 100000 t1 60 chrA >> $WORK/mix.paf
paf_line 'other|ctg'    50000  52000 t2 10 chrB >> $WORK/mix.paf
gaffilter $WORK/mix.paf -p -r 5 -m 0 -q 5 -x 60 -P 'id=REF|' > $WORK/mix.out 2>/dev/null
is "$(count $WORK/mix.out)" 3 "-P spares the protected query without sparing the rest"
is "$(grep -c '^other|ctg' $WORK/mix.out)" 1 "the unprotected query is still filtered normally"

########################################################################
# -m keeps the deletion proportionate to the overlap.  gaffilter's unit of
# removal is the whole line, so at -m 0 a sliver of contested sequence deletes
# an alignment that is almost entirely uncontested.
########################################################################

# the 100kb record is the one at risk here: it loses on MAPQ to a 200bp record they share only
# 100bp with.  -m is measured against the record being judged, so 100/100000 spares it while
# 100/200 still lets the short record be judged on its own (mostly contested) length.
paf_line q5     0 100000 t1  5 chrA  > $WORK/sliver.paf
paf_line q5 99900 100100 t2 60 chrB >> $WORK/sliver.paf

gaffilter $WORK/sliver.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/sliver.m0 2>/dev/null
is "$(count $WORK/sliver.m0)" 1 "at -m 0 a 100bp overlap deletes a 100kb record whole"
gaffilter $WORK/sliver.paf -p -r 5 -m 0.25 -q 5 -x 60 > $WORK/sliver.m25 2>/dev/null
is "$(count $WORK/sliver.m25)" 2 "at -m 0.25 the mostly-uncontested record survives"
is "$(awk '$3==0' $WORK/sliver.m25 | wc -l | tr -d ' ')" 1 "and it is the long record that -m saved"

########################################################################
# same-contig behaviour is untouched by any of the above
########################################################################

paf_line q6     0 100000 t1 10 chrA  > $WORK/same.paf
paf_line q6 50000  52000 t2 10 chrA >> $WORK/same.paf
gaffilter $WORK/same.paf -p -r 5 -m 0 -q 5 > $WORK/same.x0 2>/dev/null
gaffilter $WORK/same.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/same.x60 2>/dev/null
is "$(diff -q $WORK/same.x0 $WORK/same.x60 > /dev/null; echo $?)" 0 "-x does not change same-contig results"

########################################################################
# the filter may trim, but it must not perforate.  each cross-contig removal
# is individually justified, but one with alignment still standing on both
# sides deletes the middle of a query and invents a breakpoint.  -m cannot see
# this: it judges one pair at a time, and a hole is a property of the run.
########################################################################

# q7: a 40kb record loses to an off-contig one that outlives it, with kept
# alignment before and after -- removing it would leave a 30kb hole
paf_line q7      0 100000 t1 60 chrA  > $WORK/hole.paf
paf_line q7 100000 140000 t2 10 chrA >> $WORK/hole.paf
paf_line q7 130000 260000 t3 60 chrB >> $WORK/hole.paf
paf_line q7 260000 360000 t4 60 chrA >> $WORK/hole.paf
gaffilter $WORK/hole.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/hole.out 2>/dev/null
is "$(count $WORK/hole.out)" 4 "a cross-contig removal that would perforate is not made"
is "$(awk -F'\t' '$3==100000' $WORK/hole.out | wc -l | tr -d ' ')" 1 "the rescued record is the one in the middle"

# q8: the same loss at the START of the query, where nothing survives before it,
# is an end trim and must still be removed
paf_line q8      0  40000 t1 10 chrA  > $WORK/trim2.paf
paf_line q8   5000  40000 t2 60 chrB >> $WORK/trim2.paf
paf_line q8  40000 140000 t3 60 chrA >> $WORK/trim2.paf
gaffilter $WORK/trim2.paf -p -r 5 -m 0 -q 5 -x 60 > $WORK/trim2.out 2>/dev/null
is "$(count $WORK/trim2.out)" 2 "the same loss at a query end is still removed"

rm -rf $WORK
