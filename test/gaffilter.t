#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../:$PATH

plan tests 39

# A three-record query: two long alignments that overlap at a seam, plus a bystander.  The seam is
# what the filter is for; the flanks are what -t stops it from taking as well.
mkdir -p trim_tmp
cat > trim_tmp/in.gaf <<'EOF'
q	3000	0	1200	+	>n1:0-1200	1200	0	1200	1200	1200	60	cg:Z:1200=
q	3000	1000	2200	+	>n2:0-1200	1200	0	1200	1200	1200	60	cg:Z:1200=
q	3000	2500	3000	+	>n3:0-500	500	0	500	500	500	60	cg:Z:500=
EOF

# without -t, a tie deletes BOTH overlapping records and the bystander survives alone
gaffilter trim_tmp/in.gaf -r 5 -m 0 2>/dev/null > trim_tmp/notrim.gaf
is $(wc -l < trim_tmp/notrim.gaf) 1 "without -t a tie deletes both records"
is $(awk '$3==2500' trim_tmp/notrim.gaf | wc -l) 1 "the record with no overlap is kept"

# with -t, both survive, each having given up only the contested 200 bp
gaffilter trim_tmp/in.gaf -r 5 -m 0 -t -e 0 2>/dev/null > trim_tmp/trim.gaf
is $(wc -l < trim_tmp/trim.gaf) 3 "with -t no record is deleted"
is $(awk '$1=="q" && $3==0 && $4==1000' trim_tmp/trim.gaf | wc -l) 1 "the first record is trimmed to the seam"
is $(awk '$1=="q" && $3==1200 && $4==2200' trim_tmp/trim.gaf | wc -l) 1 "the second record is trimmed to the seam"
is $(grep -c 'cg:Z:1000=' trim_tmp/trim.gaf) 2 "the cigars are cut to match"

# the cut must follow the coordinates: block length and matches are recomputed, not carried over
is $(awk '$3==0 && $4==1000 {print $11}' trim_tmp/trim.gaf) 1000 "block length is recomputed"
is $(awk '$3==0 && $4==1000 {print $10}' trim_tmp/trim.gaf) 1000 "matches are recomputed"

# A record that dominates its overlap is not trimmed at all, so it keeps its span and no hole
# opens.  This is why the default never has a hole to close: a hole needs EVERY spanning record to
# have given the span up, and none of them gave it up to a record it dominates.
cat > trim_tmp/dom.gaf <<'EOF'
q	20000	0	14000	+	>n1:0-14000	14000	0	14000	14000	14000	60	cg:Z:14000=
q	20000	600	2000	+	>n2:0-1400	1400	0	1400	1400	1400	60	cg:Z:1400=
EOF
gaffilter trim_tmp/dom.gaf -r 5 -m 0 -t -e 0 -g 100 2>trim_tmp/dom.err > trim_tmp/dom.out
is $(awk '$3==0 && $4==14000' trim_tmp/dom.out | wc -l) 1 "a dominating record is not trimmed, so no hole opens"
is $(grep -c 'hole(s)' trim_tmp/dom.err) 0 "and nothing is reported as a hole"

# an ambiguous span is NOT closed by guessing.  A hole clips sequence out; a wrong placement puts
# a wrong alignment in, so an ambiguous span is left as a hole rather than guessed at.
cat > trim_tmp/hole.gaf <<'EOF'
q	3000	0	1400	+	>n1:0-1400	1400	0	1400	1400	1400	60	cg:Z:1400=
q	3000	600	2000	+	>n2:0-1400	1400	0	1400	1400	1400	60	cg:Z:1400=
EOF
gaffilter trim_tmp/hole.gaf -r 5 -m 0 -t -e 0 -g 100 2>trim_tmp/hole.err > trim_tmp/hole.out
is $(grep -c 'rescued 0 contested spans' trim_tmp/hole.err) 1 "an ambiguous hole is not closed by guessing"
is $(grep -c 'ambiguous' trim_tmp/hole.err) 1 "and the hole it leaves is reported, not silent"

# -R opts in to closing it anyway, on the filter's own ordering rather than raw block length
gaffilter trim_tmp/hole.gaf -r 5 -m 0 -t -e 0 -g 100 --close-holes 2>trim_tmp/weak.err > /dev/null
is $(grep -c 'rescued 1 contested spans' trim_tmp/weak.err) 1 "--close-holes closes an ambiguous hole"

# a secondary must never take a span from a primary, whatever the block lengths say
cat > trim_tmp/sec.gaf <<'EOF'
q	4000	0	3000	+	>n1:0-3000	3000	0	3000	3000	3000	60	cg:Z:3000=	tp:A:P
q	4000	0	3200	+	>n2:0-3200	3200	0	3200	3200	3200	60	cg:Z:3200=	tp:A:S
q	4000	1000	2000	+	>n3:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=	tp:A:P
EOF
gaffilter trim_tmp/sec.gaf -r 5 -m 0 -t -e 0 -g 100 -R 2>/dev/null > trim_tmp/sec.out
is $(awk '$3==1000 && $4==2000 && /tp:A:S/' trim_tmp/sec.out | wc -l) 0 "a secondary does not take the span from a primary"

# a hole spanned by no single record cannot be closed by any one claimant; it must still be counted
cat > trim_tmp/two.gaf <<'EOF'
q	4000	0	1500	+	>n1:0-1500	1500	0	1500	1500	1500	60	cg:Z:1500=
q	4000	500	1500	+	>n2:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=
q	4000	1500	3000	+	>n3:0-1500	1500	0	1500	1500	1500	60	cg:Z:1500=
q	4000	1500	2500	+	>n4:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=
EOF
gaffilter trim_tmp/two.gaf -r 5 -m 0 -t -e 0 -g 100 2>trim_tmp/two.err > /dev/null
is $(grep -c 'spanned by no single record' trim_tmp/two.err) 1 "a hole no single record spans is reported"

# a hole shorter than -g is left alone, since cactus-graphmap-join will not split a path on it
gaffilter trim_tmp/hole.gaf -r 5 -m 0 -t -e 0 -g 10000 2>trim_tmp/nohole.err > /dev/null
is $(grep -c 'rescued 0 contested spans' trim_tmp/nohole.err) 1 "a hole shorter than -g is not rescued"
is $(grep -c 'hole(s)' trim_tmp/nohole.err) 0 "and a sub-threshold hole is not reported as one"

# -t must change nothing when there is nothing to trim.  (Comparing the default against itself
# would pass no matter what -t did, which is what this test used to do.)
cat > trim_tmp/disjoint.gaf <<'EOF'
q	9000	0	1000	+	>n1:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=
q	9000	2000	3000	+	>n2:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=
r	9000	0	1000	+	>n3:0-1000	1000	0	1000	1000	1000	60	cg:Z:1000=
EOF
gaffilter trim_tmp/disjoint.gaf -r 5 -m 0 -q 5 -b 250000 -i 0.5 2>/dev/null > trim_tmp/a.gaf
gaffilter trim_tmp/disjoint.gaf -r 5 -m 0 -q 5 -b 250000 -i 0.5 -t -e 0 2>/dev/null > trim_tmp/b.gaf
is $(cmp -s trim_tmp/a.gaf trim_tmp/b.gaf && echo same) "same" "-t is byte-identical when no record loses an overlap"
is $(wc -l < trim_tmp/b.gaf) 3 "and nothing is dropped"
gaffilter trim_tmp/in.gaf -r 5 -m 0 -t -e 0 -p 2>trim_tmp/paf.err > /dev/null || true
is $(grep -c 'cannot be used with -p' trim_tmp/paf.err) 1 "-t is refused with -p"

# a record with an empty query interval yields no fragment to trim, but it is still a survivor.
# deciding what to emit from the fragments rather than the verdict silently dropped it.
cat > trim_tmp/empty.gaf <<'EOF'
q	1000	500	500	+	>n1:0-10	10	0	10	10	10	60	cg:Z:10=
q	1000	0	400	+	>n2:0-400	400	0	400	400	400	60	cg:Z:400=
EOF
gaffilter trim_tmp/empty.gaf -r 5 -m 0 -t -e 0 2>/dev/null > trim_tmp/empty.out
is $(wc -l < trim_tmp/empty.out) 2 "a kept record with an empty query interval is not dropped by -t"

# The lengths file -t needs is written by the upstream process of the same pipe
# (gaf2unstable -o), so it does not exist when gaffilter starts.  Reading it at startup raced and
# every mapping job died; it has to be read once the input is exhausted.  These paths name bare
# nodes across TWO steps, so the length cannot be inferred from the path length column and the
# file is genuinely required.
cat > trim_tmp/bare.gaf <<'EOF'
q	3000	0	1200	+	>n1>n2	2400	0	1200	1200	1200	60	cg:Z:1200=
q	3000	1000	2200	+	>n3>n4	2400	0	1200	1200	1200	60	cg:Z:1200=
EOF
printf 'n1	1200
n2	1200
n3	1200
n4	1200
' > trim_tmp/lengths.src

rm -f trim_tmp/late.tsv
( sleep 1; cp trim_tmp/lengths.src trim_tmp/late.tsv; cat trim_tmp/bare.gaf )   | gaffilter - -r 5 -m 0 -t -e 0 -l trim_tmp/late.tsv 2>/dev/null > trim_tmp/late.out
is $(awk '$3==0 && $4==1000' trim_tmp/late.out | wc -l) 1 "a lengths file that only appears after the input is still read"

# the cut drops the path steps it left behind, so the record stays canonical for gaf2paf
is "$(awk '$3==0 && $4==1000 {print $6 "/" $7}' trim_tmp/late.out)" ">n1/1200" "the trimmed path drops the step the cut left behind"

# without -l a bare multi-step path cannot be rebased; that must degrade to the old behaviour
gaffilter trim_tmp/bare.gaf -r 5 -m 0 -t -e 0 2>trim_tmp/nolen.err > trim_tmp/nolen.out
is $(grep -c 'could not be cut' trim_tmp/nolen.err) 1 "without -l an uncuttable record is deleted whole, with a warning"

# -l is documented across these tools as accepting a .fai, which has five columns.  Reading it
# with >> took columns 3 and 4 of line 1 as the next pair and silently produced a garbage map.
printf 'n1\t1200\t10\t60\t61\nn2\t1200\t20\t60\t61\nn3\t1200\t30\t60\t61\nn4\t1200\t40\t60\t61\n' > trim_tmp/lengths.fai
gaffilter trim_tmp/bare.gaf -r 5 -m 0 -t -e 0 -l trim_tmp/lengths.fai 2>trim_tmp/fai.err > trim_tmp/fai.out
is $(grep -c 'Loaded 4 node lengths' trim_tmp/fai.err) 1 "a .fai-shaped lengths file is read, not mis-parsed"
is $(awk '$3==0 && $4==1000' trim_tmp/fai.out | wc -l) 1 "and the trim still happens with it"

# with only SOME lengths known, a multi-step path's missing length must not be invented: doing so
# satisfies the total-length check by construction and can drop the wrong steps
printf 'a0\t100\n' > trim_tmp/part.tsv
cat > trim_tmp/inf.gaf <<'EOF'
q	5000	0	200	+	>x>a0	300	0	200	200	200	1	cg:Z:200=
q	5000	0	150	+	>blk:0-150	150	0	150	150	150	60	cg:Z:150=
EOF
gaffilter trim_tmp/inf.gaf -r 5 -m 0 -t -e 0 -g -1 -l trim_tmp/part.tsv 2>/dev/null > trim_tmp/inf.out
is $(grep -c '>x' trim_tmp/inf.out) 0 "a multi-step path with an unknown step length is not cut on a guess"

# -e widens each contested span on both sides: the bases butting up against an overlap are the
# least trustworthy part of the alignment.  Clipped to the record's own span, so only the interior
# border actually moves, and the dominating record -- which was never contested -- does not move.
cat > trim_tmp/edge.gaf <<'EOF'
q	200000	0	40000	+	>n1:0-40000	40000	0	40000	40000	40000	60	cg:Z:40000=
q	200000	39000	46000	+	>n2:0-7000	7000	0	7000	7000	7000	60	cg:Z:7000=
EOF
gaffilter trim_tmp/edge.gaf -r 5 -m 0 -t -e 0 2>/dev/null > trim_tmp/e0.out
gaffilter trim_tmp/edge.gaf -r 5 -m 0 -t -e 5000 2>/dev/null > trim_tmp/e5.out
is $(awk '$3==40000 && $4==46000' trim_tmp/e0.out | wc -l) 1 "-e 0 trims exactly the contested span"
is $(awk '$3==45000 && $4==46000' trim_tmp/e5.out | wc -l) 1 "-e 5000 trims 5000 further past the border"
is $(awk '$3==0 && $4==40000' trim_tmp/e5.out | wc -l) 1 "and the record that was never contested does not move"

# in a tie both sides pull back, so the hole widens by 2*-e rather than -e
cat > trim_tmp/edgetie.gaf <<'EOF'
q	200000	33000	40000	+	>n1:0-7000	7000	0	7000	7000	7000	60	cg:Z:7000=
q	200000	39000	46000	+	>n2:0-7000	7000	0	7000	7000	7000	60	cg:Z:7000=
EOF
gaffilter trim_tmp/edgetie.gaf -r 5 -m 0 -t -e 1000 2>/dev/null > trim_tmp/et.out
is $(awk '$3==33000 && $4==38000' trim_tmp/et.out | wc -l) 1 "a tie pulls the left record back by -e too"
is $(awk '$3==41000 && $4==46000' trim_tmp/et.out | wc -l) 1 "and the right record forward by -e"

# -Q: a record that loses an overlap AND is poorly placed in its own right keeps the old
# treatment, deleted whole.  The overlap was always a signal about the record, not just about the
# overlapping part, and the flanks of a record that is wrong along its length are not worth having.
cat > trim_tmp/lowq.gaf <<'EOF'
q	200000	0	40000	+	>n1:0-40000	40000	0	40000	40000	40000	60	cg:Z:40000=
q	200000	39000	46000	+	>n2:0-7000	7000	0	7000	7000	7000	9	cg:Z:7000=
EOF
gaffilter trim_tmp/lowq.gaf -r 5 -m 0 -t -e 0 -Q 0 2>/dev/null > trim_tmp/q0.out
gaffilter trim_tmp/lowq.gaf -r 5 -m 0 -t -e 0 -Q 20 2>trim_tmp/q20.err > trim_tmp/q20.out
is $(awk '$3==40000 && $4==46000' trim_tmp/q0.out | wc -l) 1 "-Q 0 trims a low-mapq loser like any other"
is $(wc -l < trim_tmp/q20.out) 1 "-Q 20 deletes it whole instead"
is $(awk '$3==0 && $4==40000' trim_tmp/q20.out | wc -l) 1 "and the record that won is untouched"
is $(grep -c 'their own mapq is under 20' trim_tmp/q20.err) 1 "and says so"

# the bar applies to the record's OWN mapq, not the winner's: a mapq-60 loser is still trimmed
cat > trim_tmp/hiq.gaf <<'EOF'
q	200000	0	40000	+	>n1:0-40000	40000	0	40000	40000	40000	60	cg:Z:40000=
q	200000	39000	46000	+	>n2:0-7000	7000	0	7000	7000	7000	60	cg:Z:7000=
EOF
gaffilter trim_tmp/hiq.gaf -r 5 -m 0 -t -e 0 -Q 20 2>/dev/null > trim_tmp/hiq.out
is $(awk '$3==40000 && $4==46000' trim_tmp/hiq.out | wc -l) 1 "a loser that is itself mapq 60 is still trimmed under -Q 20"

# getopt_long resolves unique prefixes, so a new long option can silently break an old
# abbreviation.  --r was a unique prefix of --ratio until a --rescue-weak was added next to it,
# which turned a working invocation into a fatal parse error.  The option is now --close-holes.
cat > trim_tmp/one.gaf <<'EOF'
q1	100	0	50	+	>s1:0-50	50	0	50	50	50	60	cg:Z:50=
EOF
gaffilter trim_tmp/one.gaf --r 5 > trim_tmp/abbrev.out 2>/dev/null
is $? 0 "--r is still an unambiguous abbreviation of --ratio"
is $(wc -l < trim_tmp/abbrev.out) 1 "and still filters"

rm -rf trim_tmp
