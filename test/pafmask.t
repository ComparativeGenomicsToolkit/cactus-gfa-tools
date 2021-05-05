#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../deps/hal:$PATH

plan tests 1

gzip -dc pafmask/chr20.bed.gz > chr20.bed

gzip -dc pafmask/chr20.paf.gz | pafmask - chr20.bed -v > /dev/null
is $? 0 "pafmask self validates on non-trivial input"


rm -f  chr20.bed
