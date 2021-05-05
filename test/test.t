#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 22

gzip -dc  hpp-20-2M/CHM13.fa.gz > CHM13.fa
gzip -dc  hpp-20-2M/hg38.fa.gz > hg38.fa
gzip -dc  hpp-20-2M/hg38-rev.fa.gz > hg38-rev.fa
gzip -dc pafmask/chr20.bed.gz > chr20.bed

# validate the validator
minimap2 hpp-20-2M/CHM13.fa.gz hpp-20-2M/hg38.fa.gz -c -A 30 -B 10 > CHM13_hg38.paf
python ./verify_matches.py CHM13_hg38.paf CHM13.fa hg38.fa --min-identity 0.75
is $? 0 "validator accepts minimap paf with identiy 0.75"

rm -f CHM13_hg38.paf

# make the graph
minigraph -xggs -l10k hpp-20-2M/CHM13.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > hpp-20-2M.gfa
gfatools gfa2fa hpp-20-2M.gfa > hpp-20-2M.gfa.fa

# align CHM back to it
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/CHM13.fa.gz > CHM13.gaf
mzgaf2paf CHM13.gaf > CHM13.paf
is $? 0 "mzgaf2paf doesn't crash on simple forward alignment"
python ./verify_matches.py CHM13.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for very simple forward alignment"

# clip it with pafmask
pafmask CHM13.paf chr20.bed -v > CHM13.mask.paf
python ./verify_matches.py CHM13.mask.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "masked paf checks out for very simple forward alignment"

rm -f  CHM13.paf CHM13.mask.paf

# extract the PAF with rgfa2paf
rgfa2paf hpp-20-2M.gfa > hpp-20-2M.paf
python ./verify_matches.py hpp-20-2M.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out when extracted from reference contig using rgfa2paf"
samtools faidx CHM13.fa
rgfa2paf hpp-20-2M.gfa -q CHM13.fa.fai > hpp-20-2M.paf.q
python ./verify_matches.py hpp-20-2M.paf.q CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out when extracted from reference contig using rgfa2paf and fai for query lengths"

rm -f  hpp-20-2M.paf CHM13.fa.fai hpp-20-2M.paf.q

# align a new sequence (hg38) to it
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/hg38.fa.gz > hg38.gaf
mzgaf2paf hg38.gaf > hg38.paf
is $? 0 "mzgaf2paf doesn't crash on hg38 alignment"
python ./verify_matches.py hg38.paf hg38.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for hg38 alignment"

# clip it with pafmask
pafmask hg38.paf chr20.bed -v > hg38.mask.paf
python ./verify_matches.py hg38.mask.paf hg38.fa hpp-20-2M.gfa.fa
is $? 0 "masked paf checks out for hg38 alignment"

rm -f  hg38.paf hg38.mask.paf

# repeat without gap filter
mzgaf2paf CHM13.gaf -g 0 > CHM13.paf
is $? 0 "mzgaf2paf doesn't crash on simple forward alignment with -g 0"
python ./verify_matches.py CHM13.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for very simple forward alignment with -g 0"

# repeat without gap filter
mzgaf2paf hg38.gaf -g 0 > hg38.paf
is $? 0 "mzgaf2paf doesn't crash on hg38 alignment with -g 0"
python ./verify_matches.py hg38.paf hg38.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for hg38 alignment with -g 0"

# try both at once
cat CHM13.gaf hg38.gaf > CHM13_hg38.gaf
mzgaf2paf CHM13_hg38.gaf -g 0 > CHM13_hg38.paf
cat CHM13.paf hg38.paf > CHM13_hg38_cat.paf
diff CHM13_hg38.paf CHM13_hg38_cat.paf
is $? 0 "same output when catting input as when catting output"

mzgaf2paf CHM13.gaf hg38.gaf -g 0 > CHM12_hg38_1.paf
diff CHM12_hg38_1.paf CHM13_hg38.paf
is $? 0 "same output when passing multiple input files doing separately and catting"

rm -f CHM13_hg38.paf CHM13_hg38_cat.paf CHM12_hg38_1.paf
rm -f hg38.paf CHM13.paf

# test the universal filter doesn't apply when it shouldn't
mzgaf2paf hg38.gaf -u 0 > hg38_u0.paf
is $? 0 "mzgaf2paf doesn't crash on hg38 alignment with -u 0"
mzgaf2paf hg38.gaf -u 1 > hg38_u1.paf
is $? 0 "mzgaf2paf doesn't crash on hg38 alignment with -u 1"
diff hg38_u0.paf hg38_u1.paf
is $? 0 "universal filter has no effect on single sample"

rm -f hg38_u0.paf hg38_u1.paf

# test that the reverse strand is being treated correctly
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/hg38-rev.fa.gz > hg38-rev.gaf
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/hg38.fa.gz > hg38.gaf
cat hg38-rev.gaf hg38.gaf > hg38-rf.gaf
mzgaf2paf hg38-rf.gaf -u 1 | grep -v reverse> hg38-rf.paf
mzgaf2paf hg38.gaf -u 0 > hg38.paf
#diff hg38-rf.paf hg38.paf
# this works on smaller tests, but it seems that minigraph minimizers are not perfectly symmetrical
# all the time, which means that this test doesn't work
#is $? 0 "universal filter doesn't obviously mess up on reverse strand"

rm -f hg38-rev.gaf hg38-rf.paf hg38.paf

# test that the universal filter runs without crashing
mzgaf2paf CHM13_hg38.gaf -u 0 > CHM13_hg38_u0.paf
is $? 0 "mzgaf2paf doesn't crash on CHM13_hg38 alignment with -u 0"
mzgaf2paf CHM13_hg38.gaf -u 1 > CHM13_hg38_u1.paf
is $? 0 "mzgaf2paf doesn't crash on CHM13_hg38 alignment with -u 1"
diff CHM13_hg38_u0.paf CHM13_hg38_u1.paf > /dev/null
isnt $? 0 "universal filter has effect on two samples"
cat CHM13.fa hg38.fa > CHM13_hg38.fa
python ./verify_matches.py CHM13_hg38_u1.paf CHM13_hg38.fa hpp-20-2M.gfa.fa
is $? 0 "universal filter produces valid paf"

rm -f CHM13_hg38.gaf CHM13_hg38_u0.paf CHM13_hg38_u1.paf CHM13_hg38.fa

rm -f hpp-20-2M.gfa hpp-20-2M.gfa.fa CHM13.fa hg38.fa hg38-rev.fa CHM13.gaf hg38.gaf chr20.bed


