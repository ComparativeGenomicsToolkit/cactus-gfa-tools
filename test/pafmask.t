#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 3

gzip -dc pafmask/chr20.bed.gz > chr20.bed

gzip -dc pafmask/chr20.paf.gz | pafmask - chr20.bed -v > /dev/null
is $? 0 "pafmask self validates on non-trivial input"

rm -f  chr20.bed

minimap2 pafmask/pm1.fa pafmask/pm2.fa -c > pm.paf
python ./verify_matches.py pm.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "validator accepts minimap2 output because there are indels and no snps"

printf "B\t0\t60\n" > pm1.bed
pafmask pm.paf pm1.bed > pm1.mask.paf
python ./verify_matches.py pm1.mask.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "reverse strand clip 1 produces valid paf"

printf "B\t60\t100\n" > pm2.bed
pafmask pm.paf pm2.bed > pm2.mask.paf
python ./verify_matches.py pm2.mask.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "reverse strand clip 2 produces valid paf"

minimap2 pafmask/pm2.fa pafmask/pm1.fa -c > mp.paf
python ./verify_matches.py mp.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "validator accepts minimap2 output because there are indels and no snps"

printf "A\t0\t60\n" > mp1.bed
pafmask mp.paf mp1.bed > mp1.mask.paf
python ./verify_matches.py mp1.mask.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "reverse strand clip 1 produces valid paf"

printf "A\t60\t100\n" > mp2.bed
pafmask mp.paf mp2.bed > mp2.mask.paf
python ./verify_matches.py mp2.mask.paf pafmask/pm1.fa pafmask/pm2.fa
is $? 0 "reverse strand clip 2 produces valid paf"

