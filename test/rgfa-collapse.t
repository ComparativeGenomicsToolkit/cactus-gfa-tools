#!/usr/bin/env bash

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../:$PATH

plan tests 20

# rgfa-collapse reads its sites from `vg snarls -n | vg view -Rj`, but only the two boundary names
# of each line (end first, then start), so the snarls are written by hand here and the test needs
# no vg.  Every case is its own reference contig (chrA ... chrF) with its own snarl, in one graph.
T=collapse_tmp
rm -rf $T && mkdir -p $T

python3 - $T <<'PYEOF'
import random, sys
out = sys.argv[1]
random.seed(7)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
U = 'CAGGTCCGCCTCACAGACCCACTTCTGCAGCCCCCC'          # a 36 bp VNTR unit
S, L = [], []
def seg(name, seq, sn, so, rank):
    S.append('S\t%s\t%s\tLN:i:%d\tSN:Z:%s\tSO:i:%d\tSR:i:%d' % (name, seq, len(seq), sn, so, rank))
def link(a, ao, b, bo): L.append('L\t%s\t%s\t%s\t%s\t0M' % (a, ao, b, bo))
snarls = []
def ref3(c, mid):
    """reference contig c as r1 (300 bp flank), r2 (mid, may be empty), r3 (300 bp flank)"""
    f1, f3 = rnd(300), rnd(300)
    seg(c + 'r1', f1, 'GRCh38#0#' + c, 0, 0)
    if mid:
        seg(c + 'r2', mid, 'GRCh38#0#' + c, 300, 0)
        link(c + 'r1', '+', c + 'r2', '+'); link(c + 'r2', '+', c + 'r3', '+')
    else:
        link(c + 'r1', '+', c + 'r3', '+')
    seg(c + 'r3', f3, 'GRCh38#0#' + c, 300 + len(mid), 0)
    snarls.append((c + 'r3', c + 'r1'))
def tangle(c, pieces, hap):
    """a chain of alt nodes from r1 to r3 with one branch, i.e. a minigraph STR/VNTR mess"""
    names = []
    off = 0
    for i, p in enumerate(pieces):
        n = '%sa%d' % (c, i + 1); names.append(n)
        seg(n, p, hap, off, 1); off += len(p)
    link(c + 'r1', '+', names[0], '+')
    for a, b in zip(names, names[1:]): link(a, '+', b, '+')
    link(names[2], '+', names[5], '+')                    # a branch: two paths through the site
    link(names[-1], '+', c + 'r3', '+')
    return names

# A: a VNTR tangle -- reference holds two copies, the haplotype nine in ten nodes -> flattened
ref3('chrA', U * 2); tangle('chrA', [U] * 9, 'HG1#1#ctgA')
# B: the same topology, all unique sequence -> not a tandem repeat, left alone
ref3('chrB', rnd(72)); tangle('chrB', [rnd(36) for _ in range(9)], 'HG1#1#ctgB')
# C: a simple two-allele VNTR bubble (2 nodes inside) -> under --tr-min-nodes, left alone
ref3('chrC', U * 2)
seg('chrCa1', U * 4, 'HG1#1#ctgC', 0, 1); link('chrCr1', '+', 'chrCa1', '+'); link('chrCa1', '+', 'chrCr3', '+')
# D: a pure insertion (boundaries adjacent on the reference) that is a VNTR tangle -> flattened
ref3('chrD', ''); tangle('chrD', [U] * 10, 'HG1#1#ctgD')
# E: a VNTR tangle whose longest allele is 20 copies (720 bp) -> flattened unless the cap is lower
ref3('chrE', U * 2); tangle('chrE', [U * 2] * 10, 'HG1#1#ctgE')
# F: an inversion allele: the haplotype carries the reference's middle 6 kb reverse-complemented,
# stored as novel sequence.  The alignment pass should find it and rewire it.
X = rnd(6000)
ref3('chrF', X)
seg('chrFa1', rc(X), 'HG1#1#ctgF', 0, 1); link('chrFr1', '+', 'chrFa1', '+'); link('chrFa1', '+', 'chrFr3', '+')

with open(out + '/g.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n' + '\n'.join(S) + '\n' + '\n'.join(L) + '\n')
with open(out + '/snarls.json', 'w') as f:
    for end, start in snarls:
        f.write('{"end": {"name": "%s"}, "start": {"name": "%s"}}\n' % (end, start))
PYEOF

segs() { grep -c "^S	$2" $1; }                         # $1 graph, $2 name prefix
has_link() { awk -F'\t' -v a=$2 -v b=$3 '$1=="L" && (($2==a && $4==b) || ($2==b && $4==a)) {f=1} END {print f+0}' $1; }

# ---- -T
rgfa-collapse -T --tr-report $T/tr.tsv -b 1000 $T/g.gfa $T/snarls.json > $T/t.gfa 2> $T/t.err
is $? 0 "rgfa-collapse -T runs"
is "$(grep -o 'flattened [0-9]* tandem-repeat site' $T/t.err)" "flattened 3 tandem-repeat site" "-T flattens the three tandem-repeat tangles (A, D, E)"
is "$(segs $T/t.gfa chrAa)/$(segs $T/t.gfa chrAr)" "0/3" "A: every alt node of the VNTR tangle is gone, the reference is not"
is $(segs $T/t.gfa chrAr) 3 "A: the reference path is intact"
is "$(has_link $T/t.gfa chrAr1 chrAr2)$(has_link $T/t.gfa chrAr2 chrAr3)" "11" "A: and its links"
is $(segs $T/t.gfa chrBa) 9 "B: a unique-sequence tangle is left alone"
is $(segs $T/t.gfa chrCa) 1 "C: a two-node VNTR bubble is under --tr-min-nodes and left alone"
is "$(segs $T/t.gfa chrDa)/$(segs $T/t.gfa chrDr)" "0/2" "D: a pure-insertion VNTR tangle is flattened"
is $(has_link $T/t.gfa chrDr1 chrDr3) 1 "D: leaving the reference's own link across the site"
is "$(segs $T/t.gfa chrEa)/$(segs $T/t.gfa chrEr)" "0/3" "E: a 720 bp VNTR allele is under the default --tr-max-allele and flattened"
is "$(grep -vc '^#' $T/tr.tsv 2> /dev/null)" 3 "--tr-report has one row per flattened site"
is "$(awk -F'\t' '$1 ~ /chrA$/ {print $9}' $T/tr.tsv)" "36" "--tr-report finds the 36 bp unit"

# ---- gates
rgfa-collapse -T --tr-max-allele 500 -b 1000 $T/g.gfa $T/snarls.json > $T/cap.gfa 2> /dev/null
is $(segs $T/cap.gfa chrEa) 10 "E: with --tr-max-allele 500 the 720 bp allele is left alone"
is "$(segs $T/cap.gfa chrAa)/$(segs $T/cap.gfa chrAr)" "0/3" "A: while the 324 bp one is still flattened"
rgfa-collapse -b 1000 $T/g.gfa $T/snarls.json > $T/noT.gfa 2> /dev/null
is "$(segs $T/noT.gfa chrAa)$(segs $T/noT.gfa chrDa)$(segs $T/noT.gfa chrEa)" "91010" "without -T no tangle is touched"
rgfa-collapse -T -d --tr-report $T/d.tsv -b 1000 $T/g.gfa $T/snarls.json > $T/d.gfa 2> $T/d.err
is "$(wc -c < $T/d.gfa | tr -d ' ')/$(grep -vc '^#' $T/d.tsv)" "0/3" "-d writes no graph but still reports the three sites"

# ---- the inversion pass (needs minimap2)
is "$(segs $T/t.gfa chrFa)/$(segs $T/t.gfa chrFr)" "0/3" "F: the inverted allele stored as novel sequence is removed"
is "$(awk -F'\t' '$1=="L" && (($2=="chrFr1" && $4=="chrFr2" && $5=="-") || ($2=="chrFr2" && $4=="chrFr1" && $3=="+")) {f=1} END {print f+0}' $T/t.gfa)" 1 "F: replaced by an inversion link into the reference copy"
is "$(segs $T/noT.gfa chrFa)/$(segs $T/noT.gfa chrFr)" "0/3" "F: the inversion pass does not depend on -T"
is "$(grep -c 'not reachable\|refusing to emit' $T/t.err)/$(grep -c '^S' $T/t.gfa)" "0/27" "the placement check accepts the result: the 17 reference nodes and the 10 alt nodes of B and C"

rm -rf $T
