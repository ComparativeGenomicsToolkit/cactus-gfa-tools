#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../:$PATH

plan tests 51

mkdir -p stock_tmp

# the stock rule: a tie deletes BOTH overlapping records, and a bystander survives alone
cat > stock_tmp/in.gaf <<'EOF'
q	3000	0	1200	+	>n1:0-1200	1200	0	1200	1200	1200	60	cg:Z:1200=
q	3000	1000	2200	+	>n2:0-1200	1200	0	1200	1200	1200	60	cg:Z:1200=
q	3000	2500	3000	+	>n3:0-500	500	0	500	500	500	60	cg:Z:500=
EOF
gaffilter stock_tmp/in.gaf -r 5 -m 0 2>/dev/null > stock_tmp/out.gaf
is $(wc -l < stock_tmp/out.gaf) 1 "a tie deletes both records"
is $(awk '$3==2500' stock_tmp/out.gaf | wc -l) 1 "the record with no overlap is kept"

# the trim mode (-t and its options) was removed in favour of -x: each of its options fails, and
# says what to use instead
for opt in -t --trim "-e 5000" "--trim-edge 5000" "-Q 20" "-g 100" "-l x.tsv" -R; do
    gaffilter stock_tmp/in.gaf -r 5 -m 0 $opt > stock_tmp/removed.out 2> stock_tmp/removed.err
    is "$? $(grep -c 'has been removed. Use -x/--exact' stock_tmp/removed.err) $(wc -l < stock_tmp/removed.out)" "1 1 0" "the removed $opt fails with a pointer to -x"
done

# getopt_long resolves unique prefixes, so a new long option can silently break an old
# abbreviation.  --r was a unique prefix of --ratio until a --rescue-weak was added next to it,
# which turned a working invocation into a fatal parse error.
cat > stock_tmp/one.gaf <<'EOF'
q1	100	0	50	+	>s1:0-50	50	0	50	50	50	60	cg:Z:50=
EOF
gaffilter stock_tmp/one.gaf --r 5 > stock_tmp/abbrev.out 2>/dev/null
is $? 0 "--r is still an unambiguous abbreviation of --ratio"
is $(wc -l < stock_tmp/abbrev.out) 1 "and still filters"

rm -rf stock_tmp

# ---- -x/--exact: per-segment resolution.  Six 50 kb reference nodes on one chromosome; every
# record aligns along them, so placements, the backbone and the one-to-one test are all checkable
# by hand.
mkdir -p exact_tmp
for i in 1 2 3 4 5 6; do
    printf "r$i\t50000\tid=REF|chrA\t%d\t0\n" $(( (i - 1) * 50000 ))
done > exact_tmp/nodes.tsv
X="-x -r 5 -m 0.25 -q 5 -b 0 -i 0.5 --exact-nodes exact_tmp/nodes.tsv"

# tier: D (40 kb) ties C over 20 kb, half its block, so the stock rule deletes it whole.  -x keeps
# the 10 kb nobody else claims.  D also beats X on MAPQ, but X is not demoted, so D yields to it.
cat > exact_tmp/tier.gaf <<'EOF'
id=S.1|c1	210000	0	100000	+	>r1>r2	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	210000	80000	120000	+	>r2>r3	100000	30000	70000	40000	40000	60	cg:Z:40000=	rc:Z:chrA
id=S.1|c1	210000	110000	210000	+	>r3>r4>r5	150000	10000	110000	100000	100000	10	cg:Z:100000=	rc:Z:chrA
EOF
gaffilter exact_tmp/tier.gaf -r 5 -m 0.25 -q 5 -b 0 -i 0.5 2>/dev/null > exact_tmp/tier.stock
gaffilter exact_tmp/tier.gaf $X --exact-plan exact_tmp/tier.plan --exact-summary exact_tmp/tier.summary 2>exact_tmp/tier.err > exact_tmp/tier.out
is $(wc -l < exact_tmp/tier.stock) 2 "-x tier: the stock filter deletes the demoted record whole"
is $(wc -l < exact_tmp/tier.out) 3 "-x keeps it"
is $(grep -c kq:Z exact_tmp/tier.out) 1 "only the demoted record is cut"
is $(awk '$3==80000' exact_tmp/tier.out | grep -o 'kq:Z:[^[:space:]]*') "kq:Z:100000-110000" "it keeps only the span no non-demoted record claims"
is $(awk '$3==80000 {print $10 "/" $11}' exact_tmp/tier.out) "40000/40000" "a cut record is printed whole (its block length is the parent's)"
is $(awk '$3==80000 {print $7}' exact_tmp/tier.plan) demoted "the plan marks it demoted"
is "$(cat exact_tmp/tier.summary)" "$(grep '^\[gaffilter\]: -x' exact_tmp/tier.err | sed 's/^\[gaffilter\]: //')" "--exact-summary holds the summary printed on stderr"

# gaf2paf applies kq:Z: per line, and every line keeps the parent's gl/gm
printf 'r1\t50000\nr2\t50000\nr3\t50000\nr4\t50000\nr5\t50000\nr6\t50000\n' > exact_tmp/lens.tsv
gaf2paf exact_tmp/tier.out -l exact_tmp/lens.tsv > exact_tmp/tier.paf
is "$(awk '$6=="r3" && $3==100000' exact_tmp/tier.paf | cut -f 3,4,8,9)" "$(printf '100000\t110000\t0\t10000')" "gaf2paf cuts the line to the kept span"
is $(awk '$3==100000 && $4==110000' exact_tmp/tier.paf | grep -c 'gl:i:40000') 1 "and the cut line keeps the parent's block length"
is $(awk '$3 >= 80000 && $4 <= 100000 && /gl:i:40000/' exact_tmp/tier.paf | wc -l) 0 "and nothing outside the kept span is printed"

# floor: the same remainder at 90% identity is dropped
sed 's/cg:Z:40000=/cg:Z:20000=1000X19000=/; s/40000\t40000\t60/39000\t40000\t60/' exact_tmp/tier.gaf > exact_tmp/floor.gaf
gaffilter exact_tmp/floor.gaf $X --exact-log exact_tmp/floor.log 2>/dev/null > exact_tmp/floor.out
is $(wc -l < exact_tmp/floor.out) 2 "-x drops a remainder under the identity floor"
is $(grep -c '^floor' exact_tmp/floor.log) 1 "and logs it"

# R4 (one-to-one): a remainder placed on reference its contig already covers is dropped
cat > exact_tmp/r4.gaf <<'EOF'
id=S.1|c1	200000	0	100000	+	>r1>r2	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	200000	80000	160000	+	>r1>r2	100000	0	80000	80000	80000	60	cg:Z:80000=	rc:Z:chrA
EOF
gaffilter exact_tmp/r4.gaf $X --exact-log exact_tmp/r4.log 2>/dev/null > exact_tmp/r4.out
is $(wc -l < exact_tmp/r4.out) 1 "-x drops a remainder that would place reference twice"
is $(grep -c '^collide' exact_tmp/r4.log) 1 "as a one-to-one collision"
# ...and reference covered by another contig of the same haplotype counts too
cat > exact_tmp/r4b.gaf <<'EOF'
id=S.1|c0	100000	0	100000	+	>r1>r2	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	200000	0	100000	+	>r3>r4	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	200000	80000	160000	+	>r1>r2	100000	0	80000	80000	80000	60	cg:Z:80000=	rc:Z:chrA
EOF
gaffilter exact_tmp/r4b.gaf $X --exact-log exact_tmp/r4b.log 2>/dev/null > exact_tmp/r4b.out
is $(grep -c '^collide' exact_tmp/r4b.log) 1 "-x counts reference the haplotype's other contigs cover"

# isolation: a remainder inverted against the backbone at the contig end would make a one-sided
# junction.  It keeps its sequence, but backs off the junction by the 31 kb gap
cat > exact_tmp/iso.gaf <<'EOF'
id=S.1|c1	200000	0	100000	+	>r1>r2	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	200000	80000	160000	-	>r4>r5	100000	20000	100000	80000	80000	60	cg:Z:80000=	rc:Z:chrA
EOF
gaffilter exact_tmp/iso.gaf $X --exact-log exact_tmp/iso.log --junctions exact_tmp/iso.j 2>/dev/null > exact_tmp/iso.out
is $(awk '$3==80000' exact_tmp/iso.out | grep -o 'kq:Z:[^[:space:]]*') "kq:Z:131000-160000" "-x isolates a one-sided excursion by the gap"
is $(grep -c 'isolate:R2:end-joined' exact_tmp/iso.log) 1 "and logs why"
gaffilter exact_tmp/iso.gaf $X --gap 21000 2>/dev/null | awk '$3==80000' | grep -o 'kq:Z:[^[:space:]]*' > exact_tmp/iso21.kq
is $(cat exact_tmp/iso21.kq) "kq:Z:121000-160000" "--gap sets how far"

# the guard: with --exact-ratio 2, w beats l on MAPQ (3x) where -r 5 would tie.  l is the more
# similar copy over the whole contested span, so the guard cancels w's win there
gen_w() { printf '40000='; for i in $(seq 1200); do printf '49=1X'; done; }
printf "id=S.1|c1\t300000\t0\t100000\t+\t>r1>r2\t100000\t0\t100000\t98800\t100000\t60\tcg:Z:%s\trc:Z:chrA\n" "$(gen_w)" > exact_tmp/guard.gaf
printf "id=S.1|c1\t300000\t40000\t300000\t+\t>r1>r2>r3>r4>r5>r6\t300000\t40000\t300000\t260000\t260000\t20\tcg:Z:260000=\trc:Z:chrA\n" >> exact_tmp/guard.gaf
gaffilter exact_tmp/guard.gaf $X --exact-ratio 2 2>/dev/null > exact_tmp/noguard.out
gaffilter exact_tmp/guard.gaf $X --exact-ratio 2 --guard --exact-log exact_tmp/guard.log 2>/dev/null > exact_tmp/guard.out
is $(awk '$3==0' exact_tmp/noguard.out | grep -c kq:Z) 0 "without the guard the lower ratio gives w the contested span"
is $(awk '$3==0' exact_tmp/guard.out | grep -o 'kq:Z:[^[:space:]]*') "kq:Z:0-40000" "the guard cancels it where l is the more similar"
is $(awk '$3==40000' exact_tmp/guard.out | grep -o 'kq:Z:[^[:space:]]*') "kq:Z:100000-300000" "and l does not get it either: it lost on MAPQ"
is $(grep -c '^guard-applied' exact_tmp/guard.log) 1 "the veto is logged as applied"
gaffilter exact_tmp/guard.gaf $X --guard 2>/dev/null > exact_tmp/guard5.out
gaffilter exact_tmp/guard.gaf $X 2>/dev/null > exact_tmp/plain5.out
is $(cmp -s exact_tmp/guard5.out exact_tmp/plain5.out && echo same) same "the guard does nothing when --exact-ratio is -r"

# the result does not depend on the input order
tac exact_tmp/tier.gaf > exact_tmp/tier.rev.gaf
gaffilter exact_tmp/tier.rev.gaf $X 2>/dev/null | sort > exact_tmp/tier.rev.out
is $(sort exact_tmp/tier.out | cmp -s - exact_tmp/tier.rev.out && echo same) same "-x does not depend on the input order"

# -x is GAF-only
gaffilter exact_tmp/tier.gaf $X -p 2>exact_tmp/xp.err > /dev/null || true
is $(grep -c 'cannot be used with -p or -o' exact_tmp/xp.err) 1 "-x is refused with -p"

# isolation never cuts stock-anchored sequence from a backbone record at a junction that is the
# stock chain's own.  The backbone b is new only at its start, an off-reference node n0 that the
# line filter drops (0% identity), far from its junction to the inverted record v at the contig
# end.  v is anchored whole by the stock chain, so nothing at the junction is new and the junction
# stays as the stock filter has it (this was the yeast SK1 chrI case)
printf 'n0\t30000\tid=S.0|x\t0\t1\n' > exact_tmp/side.nodes.tsv
cat exact_tmp/nodes.tsv >> exact_tmp/side.nodes.tsv
XS="-x -r 5 -m 0.25 -q 5 -b 0 -i 0.5 --exact-nodes exact_tmp/side.nodes.tsv"
cat > exact_tmp/side.gaf <<'EOF'
id=S.1|c1	180000	0	130000	+	>n0>r1>r2	130000	0	130000	100000	130000	60	cg:Z:30000X100000=	rc:Z:chrA
id=S.1|c1	180000	130000	180000	-	>r5	50000	0	50000	50000	50000	60	cg:Z:50000=	rc:Z:chrA
EOF
gaffilter exact_tmp/side.gaf -r 5 -m 0.25 -q 5 -b 0 -i 0.5 2>/dev/null > exact_tmp/side.stock
gaffilter exact_tmp/side.gaf $XS --exact-log exact_tmp/side.log 2>/dev/null > exact_tmp/side.out
is $(grep -c '^a0-junction' exact_tmp/side.log) 1 "-x leaves a junction the stock chain made (a0-junction)"
is $(grep -c '^isolate' exact_tmp/side.log) 0 "and isolates nothing"
is $(cmp -s exact_tmp/side.out exact_tmp/side.stock && echo same) same "so it keeps exactly what the stock filter keeps"
# ...but a backbone that is new AT the junction is still trimmed there
sed 's/>n0>r1>r2/>r1>r2>n0/; s/cg:Z:30000X100000=/cg:Z:100000=30000X/' exact_tmp/side.gaf > exact_tmp/side2.gaf
gaffilter exact_tmp/side2.gaf $XS --exact-log exact_tmp/side2.log 2>/dev/null > exact_tmp/side2.out
is $(grep -c 'isolate:R2:end-joined' exact_tmp/side2.log) 1 "-x isolates a backbone that is new at the junction"
is $(awk '$3==0' exact_tmp/side2.out | grep -o 'kq:Z:[^[:space:]]*') "kq:Z:0-99000" "by trimming its junction end back by the gap"

# eligibility (-i) reads the identity as the line filter downstream does: gaf2paf's gi:f:, rounded
# to 3 places.  l is 49,960/100,000 = 0.4996, which rounds to 0.5: it competes (and loses its
# remainder to the identity floor) rather than passing through uncontested, its lines kept anyway
cat > exact_tmp/round.gaf <<'EOF'
id=S.1|c1	200000	0	100000	+	>r1>r2	100000	0	100000	100000	100000	60	cg:Z:100000=	rc:Z:chrA
id=S.1|c1	200000	50000	150000	+	>r3>r4	100000	0	100000	49960	100000	10	cg:Z:49960=50040X	rc:Z:chrA
EOF
gaffilter exact_tmp/round.gaf $X --exact-plan exact_tmp/round.plan 2>/dev/null > exact_tmp/round.out
is $(awk '$3==50000' exact_tmp/round.out | wc -l) 0 "-x: a record whose gi:f: rounds up to -i competes"
is $(awk '$3==50000 {print $7}' exact_tmp/round.plan) demoted "and is demoted, not passed through as ineligible"
sed 's/49960	100000	10	cg:Z:49960=50040X/49900	100000	10	cg:Z:49900=50100X/' exact_tmp/round.gaf > exact_tmp/round2.gaf
gaffilter exact_tmp/round2.gaf $X --exact-plan exact_tmp/round2.plan 2>/dev/null > exact_tmp/round2.out
is $(awk '$3==50000 {print $7}' exact_tmp/round2.plan) ineligible "one whose gi:f: is under -i is still ineligible"

# a reused GAF resolved against a since-extended graph (cactus --inGAF) can carry path steps its
# alignment never enters; cactus trims them only after gaffilter, so -x has to cope with them
printf 'id=S.1|c1\t200000\t0\t100000\t+\t>r1>r2>r3>r4>r5\t250000\t60000\t160000\t100000\t100000\t60\tcg:Z:100000=\trc:Z:chrA\n' > exact_tmp/regran.gaf
gaffilter exact_tmp/regran.gaf $X > exact_tmp/regran.out 2> exact_tmp/regran.err
is $? 0 "-x takes a record whose path offsets reach past its end steps"
is "$(cat exact_tmp/regran.out)" "$(cat exact_tmp/regran.gaf)" "and passes it on unchanged"

# a duplicated record is a warning, not an error: the copies tie, and both go, as in the stock filter
cat exact_tmp/tier.gaf exact_tmp/tier.gaf > exact_tmp/dup.gaf
gaffilter exact_tmp/dup.gaf $X > exact_tmp/dup.out 2> exact_tmp/dup.err
is $? 0 "-x takes duplicated records"
is $(wc -l < exact_tmp/dup.out) $(gaffilter exact_tmp/dup.gaf -r 5 -m 0.25 -q 5 -b 0 -i 0.5 2>/dev/null | wc -l) "and keeps what the stock filter keeps of them"

# option values that cannot work are refused rather than hanging or placing everything
gaffilter exact_tmp/guard.gaf $X --exact-ratio 2 --guard --guard-chunk 0 > /dev/null 2> exact_tmp/bad.err
is $? 1 "-x refuses --guard-chunk 0"
gaffilter exact_tmp/tier.gaf $X --exact-ratio 0 > /dev/null 2> exact_tmp/bad.err
is $? 1 "-x refuses --exact-ratio 0"

rm -rf exact_tmp
