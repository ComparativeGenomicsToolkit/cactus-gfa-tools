#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 3

gzip -dc  hpp-20-2M/CHM13.fa.gz > CHM13.fa
gzip -dc  hpp-20-2M/hg38.fa.gz > hg38.fa
gzip -dc  hpp-20-2M/hg38-rev.fa.gz > hg38-rev.fa
gzip -dc pafmask/chr20.bed.gz > chr20.bed

gzip -dc pafmask/chr20.paf.gz | pafmask - chr20.bed -v > /dev/null
is $? 0 "pafmask self validates on non-trivial input"

minimap2 hpp-20-2M/CHM13.fa.gz hpp-20-2M/hg38.fa.gz -c -A 30 -B 10 > CHM13_hg38.paf
pafmask CHM13_hg38.paf chr20.bed > CHM13_hg38.mask.paf
python ./verify_matches.py CHM13_hg38.mask.paf CHM13.fa hg38.fa
is $? 0 "validator accepts pafmask output"

rm -f CHM13_hg38.paf CHM13_hg38.mask.paf

minimap2 hpp-20-2M/CHM13.fa.gz hpp-20-2M/hg38-rev.fa.gz -c -A 30 -B 10 > CHM13_hg38-rev.paf
pafmask CHM13_hg38-rev.paf chr20.bed > CHM13_hg38-rev.mask.paf
python ./verify_matches.py CHM13_hg38-rev.mask.paf CHM13.fa hg38-rev.fa
is $? 0 "validator accepts pafmask output with reverse hg38"

rm -f CHM13_hg38-rev.paf CHM13_hg38-rev.mask.paf

rm -f CHM13.fa hg38.fa hg38-rev.fa  chr20.bed
