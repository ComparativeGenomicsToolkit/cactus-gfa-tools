#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 2

# make the graph
minigraph -xggs -l10k hpp-20-2M/CHM13.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > hpp-20-2M.gfa
gfatools gfa2fa hpp-20-2M.gfa > hpp-20-2M.gfa.fa

# align CHM back to it
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/CHM13.fa.gz > CHM13.gaf
mzgaf2paf CHM13.gaf > CHM13.paf
gzip -dc  hpp-20-2M/CHM13.fa.gz > CHM13.fa
python ./verify_matches.py CHM13.paf CHM13.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for very simple forward alignment"

rm -f  CHM13.gaf CHM13.paf CHM13.fa

# align a new sequence (hg38) to it
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/hg38.fa.gz > hg38.gaf
mzgaf2paf hg38.gaf > hg38.paf
gzip -dc  hpp-20-2M/hg38.fa.gz > hg38.fa
python ./verify_matches.py hg38.paf hg38.fa hpp-20-2M.gfa.fa
is $? 0 "paf checks out for hg38 alignment"

rm -f  hg38.gaf hg38.paf hg38.fa

rm -f hpp-20-2M.gfa hpp-20-2M.gfa.fa
