#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH
# the tools this repository builds, as gaffilter.t does (../bin only exists in a cactus build)
PATH=../:$PATH

plan tests 19

gzip -dc  hpp-20-2M/CHM13.fa.gz > CHM13.fa
gzip -dc  hpp-20-2M/hg38.fa.gz > hg38.fa
gzip -dc  hpp-20-2M/hg38-rev.fa.gz > hg38-rev.fa
gzip -dc pafmask/chr20.bed.gz > chr20.bed
gzip -dc hpp-20-2M/CHM13.fa.gz hpp-20-2M/hg38.fa.gz hpp-20-2M/hg38-rev.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > all.fa
samtools faidx all.fa

# validate the validator
minimap2 hpp-20-2M/CHM13.fa.gz hpp-20-2M/hg38.fa.gz -c -A 30 -B 10 > CHM13_hg38.paf
python ./verify_matches.py CHM13_hg38.paf CHM13.fa hg38.fa --min-identity 0.75
is $? 0 "validator accepts minimap paf with identiy 0.75"

rm -f CHM13_hg38.paf

# make the graph
minigraph -cxggs -l10k hpp-20-2M/CHM13.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > hpp-20-2M.gfa
gfatools gfa2fa hpp-20-2M.gfa > hpp-20-2M.gfa.fa
samtools faidx hpp-20-2M.gfa.fa

# align CHM back to it
minigraph -cxasm -t $(nproc) -K4g --inv=no  hpp-20-2M.gfa hpp-20-2M/CHM13.fa.gz > CHM13.gaf
gaf2paf CHM13.gaf -l all.fa.fai > CHM13.paf
is $? 0 "gaf2paf doesn't crash on simple forward alignment"
python ./verify_matches.py CHM13.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for very simple forward alignment"
gaf2unstable CHM13.gaf -g hpp-20-2M.gfa | gaf2paf - -l hpp-20-2M.gfa.fa.fai > CHM13.u.paf
is $? 0 "gaf2paf doesn't crash on simple unstable forward alignment"
python ./verify_matches.py CHM13.u.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for very simple unstable forward alignment"

rm -f  CHM13.paf CHM13.u.paf

# now try a simple reverse case with stable coordinates
minigraph -cxggs -l10k hg38-rev.fa hpp-20-2M/CHM13.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > hg38-rev.gfa
minigraph -cxasm -t $(nproc) -K4g --inv=no  hg38-rev.gfa hg38.fa > hg38-rev.gaf
gaf2paf hg38-rev.gaf -l all.fa.fai > hg38-rev.paf
python ./verify_matches.py hg38-rev.paf CHM13.fa all.fa
is $? 0 "paf checks out for simple reverse strand"
gaf2unstable  hg38-rev.gaf -g hg38-rev.gfa -o hg38-rev-gfa.lengths > hg38-rev.u.gaf
gaf2paf hg38-rev.u.gaf -l  hg38-rev-gfa.lengths > hg38-rev.u.paf
is $? 0 "gaf2paf doesn't crash on simple reverse alignment"
gfatools gfa2fa  hg38-rev.gfa >  hg38-rev.gfa.fa
python ./verify_matches.py hg38-rev.u.paf  all.fa hg38-rev.gfa.fa  
is $? 0 "paf checks out for simple unstable reverse strand"

rm -f hg38-rev.gfa hg38-rev.gaf  \hg38-rev.paf hg38-rev.u.gaf hg38-rev.u.paf hg38-rev.gfa.fa

# align a new sequence (hg38) to it
minigraph -cxasm -t $(nproc) -K4g --inv=no  hpp-20-2M.gfa hpp-20-2M/hg38.fa.gz > hg38.gaf
gaf2paf hg38.gaf -l all.fa.fai > hg38.paf
is $? 0 "gaf2paf doesn't crash on hg38 alignment"
python ./verify_matches.py hg38.paf hg38.fa all.fa
is $? 0 "paf checks out for hg38 alignment"
gaf2unstable hg38.gaf -g hpp-20-2M.gfa | gaf2paf - -l hpp-20-2M.gfa.fa.fai > hg38.u.paf
is $? 0 "gaf2paf doesn't crash on unstable hg38 alignment"
python ./verify_matches.py hg38.u.paf hpp-20-2M.gfa.fa all.fa
is $? 0 "paf checks out for unstable hg38 alignment"

rm -f  hg38.paf hg38.u.paf

rm -f hpp-20-2M.gfa hpp-20-2M.gfa.fa CHM13.fa hg38.fa hg38-rev.fa CHM13.gaf hg38.gaf chr20.bed all.fa all.fa.fai

# kq:Z: (gaffilter -x): a record keeps only the listed query intervals, cut line by line.  One
# record over three 100 bp nodes: n1 all '=', n2 with a 5 bp deletion at query 130, n3 ending in
# 50 mismatches
mkdir -p kq_tmp
printf 'n1\t100\nn2\t100\nn3\t100\n' > kq_tmp/lens.tsv
rec="q	1000	0	295	+	>n1>n2>n3	300	0	300	245	300	60	cg:Z:130=5D115=50X	tp:A:P"
printf '%s\n' "$rec" > kq_tmp/whole.gaf
gaf2paf kq_tmp/whole.gaf -l kq_tmp/lens.tsv > kq_tmp/whole.paf
printf '%s\tkq:Z:0-100,120-140,250-295\n' "$rec" > kq_tmp/cut.gaf
gaf2paf kq_tmp/cut.gaf -l kq_tmp/lens.tsv > kq_tmp/cut.paf
is "$(head -1 kq_tmp/cut.paf)" "$(head -1 kq_tmp/whole.paf)" "kq: a line inside a kept interval is printed verbatim"
is "$(awk '$6=="n2"' kq_tmp/cut.paf | cut -f 3,4,8,9,10,11)" "$(printf '120\t140\t20\t45\t20\t25')" "kq: a cut line gets its own coordinates and counts"
is $(awk '$6=="n2"' kq_tmp/cut.paf | grep -o 'cg:Z:[^[:space:]]*') "cg:Z:10=5D10=" "kq: a deletion strictly inside the kept interval stays"
is $(awk '$6=="n2"' kq_tmp/cut.paf | grep -c 'gm:i:245	gl:i:300') 1 "kq: a cut line keeps the parent's gm/gl"
is $(awk '$6=="n3"' kq_tmp/cut.paf | wc -l) 0 "kq: a piece with no '=' column left makes no line"
printf '%s\tkq:Z:100-130\n' "$rec" > kq_tmp/edge.gaf
is "$(gaf2paf kq_tmp/edge.gaf -l kq_tmp/lens.tsv | grep -o 'cg:Z:[^[:space:]]*')" "cg:Z:30=" "kq: a deletion at a cut end is dropped"
printf '%s\tkq:Z:100-oops\n' "$rec" > kq_tmp/bad.gaf
gaf2paf kq_tmp/bad.gaf -l kq_tmp/lens.tsv > /dev/null 2> kq_tmp/bad.err
is $? 1 "kq: a malformed tag is an error"
rm -rf kq_tmp
