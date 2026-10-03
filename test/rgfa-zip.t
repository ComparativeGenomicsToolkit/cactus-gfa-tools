#!/usr/bin/env bash
# rgfa-zip: the spec's synthetic cases (rgfa-zip-spec.md, Tests), end to end with real minimap2.
#
# Z1-Z43; Z24 (with variants for copies in series, a non-private member, two members of one group
# and rounds), Z24b, Z24c, Z40 and Z41 are the alt-vs-alt (v2) cases.  The fixtures are
# generated below (seeded); snarls are hand-written top-level snarls, so the test needs no vg.  Alt
# sequences carry no substitution within 60 bp of a segment end, and the bases at segment junctions
# are chosen so that no alignment can extend across a junction by a matching base, so minimap2 ends
# its records at the junctions.  (One exception: in A rc(B) C, minimap2 2.30 chains A and C and can
# stop A's record ~100 bp short; that stretch then stays alt as a small bubble.  Z5-Z9 therefore
# check exact cut points on the same chains injected with --inject-chains, as do Z27, Z38 and Z39.)
# Checks are exact strings: S and L lines, report columns, cut points, and every haplotype walk
# spelled before and after the edit (its image must spell the sequence the zips predict).  Z19 also
# runs with GafWalks.  Z36 is the SMN shape scaled down (copies 150 kb apart); Z43's satellite is
# synthetic and checks the screen (the real chr20:31.07 HSat3 window is too big for a fixture).
# The last two blocks are regression tests of the review fixes: overlaps at strand switches trimmed
# from the right end (GC-1; Z5o, Z8o, Z5f, Z5r), tips, fold-backs and leaks next to a site (GC-2,
# RD-2; Z1t, Z1h, Z23z), duplicate links (GC-3; Z1d), and the robustness findings RD-1 to RD-16
# (the Z28 fake minimap2 grows modes for a short index, a full disk, a slow failure, a hang that a
# SIGTERM must end, and a report path that turns into a directory).  The final block covers the
# second round: N gaps inside I and D runs, which minimap2 leaves out of PAF column 11 (Zn1, Zn2;
# they were paf-inconsistent), the report's kept_bp counting a compatibly zipped interval once
# (Z13), rule (G) refusing a piece that a GAF record also reads elsewhere (Zg), log times that
# separate waiting for the -j budget from running and quote minimap2's own times (Z31 with a fake
# minimap2; Z24 for the alt-vs-alt pass and the planning time), and the 500 sub-record cap shown in
# the report (Ztr).
#
# minimap2 comes from MINIMAP2, else PATH; the expected strings were made with minimap2 2.30, the
# version cactus pins.  Outputs go to $RGFA_ZIP_TEST_DIR/rgfa-zip.t.tmp (default ./rgfa-zip.t.tmp),
# removed when every test passes.

BASH_TAP_ROOT=./bash-tap
. ${BASH_TAP_ROOT}/bash-tap-bootstrap

PATH=../bin:$PATH
PATH=../:$PATH

MM2=${MINIMAP2:-$(command -v minimap2)}
if [ -z "$MM2" ] || ! "$MM2" --version > /dev/null 2>&1; then
    plan skip_all "rgfa-zip's test needs minimap2 (on PATH or in MINIMAP2)"
    exit 0
fi
MM2=$(cd "$(dirname "$MM2")" && pwd)/$(basename "$MM2")

plan tests 145

T=${RGFA_ZIP_TEST_DIR:-.}/rgfa-zip.t.tmp
rm -rf $T && mkdir -p $T/tmp
export T=$(cd $T && pwd)

python3 - $T <<'PYEOF'
#!/usr/bin/env python3
# rgfa-zip.t fixtures.  One directory per case under T: in.gfa, snarls.json (hand-written top-level
# snarls, so the test needs no vg), walks.tsv (label, input walk, its sequence, and the sequence its
# image must spell after the edit), optional inject.tsv and expected.gfa.  Every case is seeded.
# Alt sequences are mutated copies with no substitution within 60 bp of a segment end, and the bases
# at segment junctions are chosen so that no alignment can extend across a junction by a matching
# base (junctions_ok; search() takes the first seed that passes).
import os, sys, json, random

T = sys.argv[1]
REF = 'REF#0#chr1'
COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
def comp(c): return c.translate(COMP)

class Case:
    def __init__(self, name, seed):
        self.name, self.seed, self.r = name, seed, random.Random(seed)
        self.nodes, self.links, self.snarls, self.walks, self.inject = [], [], [], [], []
        self.files = {}
    def rnd(self, n): return ''.join(self.r.choice('ACGT') for _ in range(n))
    def mut(self, s, rate, margin=60):
        s = list(s)
        for i in range(margin, len(s) - margin):
            if self.r.random() < rate: s[i] = self.r.choice([c for c in 'ACGT' if c != s[i]])
        return ''.join(s)
    def mutn(self, s, n, margin=60):
        s = list(s)
        for p in self.r.sample(range(margin, len(s) - margin), n): s[p] = self.r.choice([c for c in 'ACGT' if c != s[p]])
        return ''.join(s)
    def node(self, seq, sn, so, sr):
        self.nodes.append((len(self.nodes) + 1, seq, sn, so, sr))
        return 's%d' % len(self.nodes)
    def ref(self, seqs, sn=REF, snarl=True):
        ids, so = [], 0
        for s in seqs:
            ids.append(self.node(s, sn, so, 0)); so += len(s)
        for a, b in zip(ids, ids[1:]): self.link(a, '+', b, '+', 0)
        if snarl: self.snarls.append((ids[0], ids[-1]))
        return ids
    def alt(self, seq, contig='HG1#1#c1', sr=1, so=0): return self.node(seq, contig, so, sr)
    def chain_run(self, seqs, contig, sr, dep, arr):
        """an alt run (contiguous SO, ++ links) from dep+ to arr+ with creator links"""
        ids, so = [], 0
        for s in seqs:
            ids.append(self.alt(s, contig, sr, so)); so += len(s)
        for a, b in zip(ids, ids[1:]): self.link(a, '+', b, '+', sr)
        self.link(dep, '+', ids[0], '+', sr); self.link(ids[-1], '+', arr, '+', sr)
        return ids
    def link(self, a, ao, b, bo, sr=None): self.links.append((a, ao, b, bo, sr))
    def seq(self, n): return self.nodes[int(n[1:]) - 1][1]
    def hseq(self, h): return self.seq(h[:-1]) if h[-1] == '+' else rc(self.seq(h[:-1]))
    def spell(self, walk): return ''.join(self.hseq(h) for h in walk)
    def walk(self, label, walk, after=None):
        self.walks.append((label, walk, self.spell(walk), after if after is not None else self.spell(walk)))
    def chain(self, key, parts):
        for qs, qe, st, ts, te, cg in parts: self.inject.append('%s\tchain\t%d\t%d\t%s\t%d\t%d\t%s' % (key, qs, qe, st, ts, te, cg))
    def gfa(self, no_sr=False, drop_sn=None, only_sn=None):
        out = ['H\tVN:Z:1.0']
        keep = set()
        for i, s, sn, so, sr in self.nodes:
            if only_sn and not any(sn.endswith(x) for x in only_sn): continue
            keep.add('s%d' % i)
            tags = 'LN:i:%d\tSN:Z:%s\tSO:i:%d\tSR:i:%d' % (len(s), sn, so, sr)
            if drop_sn == 's%d' % i: tags = 'LN:i:%d\tSO:i:%d\tSR:i:%d' % (len(s), so, sr)
            out.append('S\ts%d\t%s\t%s' % (i, s, tags))
        for a, ao, b, bo, sr in self.links:
            if a not in keep or b not in keep: continue
            out.append('L\t%s\t%s\t%s\t%s\t0M' % (a, ao, b, bo) + ('' if no_sr or sr is None else '\tSR:i:%d' % sr))
        return out
    def write(self, sub='', shuffle=False, **kw):
        d = os.path.join(T, self.name + sub)
        os.makedirs(d, exist_ok=True)
        lines = self.gfa(**kw)
        kept = set(l.split('\t')[1] for l in lines if l.startswith('S\t'))
        sn = [json.dumps({"start": {"name": a}, "end": {"name": b}}) for a, b in self.snarls if a in kept and b in kept]
        if shuffle:
            rr = random.Random(11)
            body = lines[1:]; rr.shuffle(body); lines = lines[:1] + body
            random.Random(13).shuffle(sn)
        open(os.path.join(d, 'in.gfa'), 'w').write('\n'.join(lines) + '\n')
        open(os.path.join(d, 'snarls.json'), 'w').write('\n'.join(sn) + '\n')
        with open(os.path.join(d, 'walks.tsv'), 'w') as f:
            for label, walk, before, after in self.walks:
                f.write('%s\t%s\t%s\t%s\n' % (label, ','.join(walk), before, after))
        if self.inject:
            open(os.path.join(d, 'inject.tsv'), 'w').write('\n'.join(self.inject) + '\n')
        for name, text in self.files.items():
            open(os.path.join(d, name), 'w').write(text)
        return d

def junctions_ok(Q, W, segs):
    """segs: (qs, qe, ts, te, strand) in query order, targets relative to W.  False if some record's
    alignment could extend across a segment end by a matching base."""
    for qs, qe, ts, te, s in segs:
        for qi, ti_p, ti_m in ((qe, te, ts - 1), (qs - 1, ts - 1, te)):
            if not (0 <= qi < len(Q)): continue
            ti, qb = (ti_p, Q[qi]) if s == '+' else (ti_m, comp(Q[qi]))
            if 0 <= ti < len(W) and qb == W[ti]: return False
    return True

def search(build, seed, tries=200):
    """the first seed from `seed` whose case passes its junction check"""
    for k in range(tries):
        c = build(seed + k)
        if c is not None: return c
    raise RuntimeError('no seed with clean junctions')

def gfa_expected(nodes, links):
    """rgfa-zip's output text for these nodes and links (S by id, L by canonical side key)"""
    out = ['H\tVN:Z:1.0']
    for i, seq, sn, so, sr in sorted(nodes):
        out.append('S\ts%d\t%s\tLN:i:%d\tSN:Z:%s\tSO:i:%d\tSR:i:%d' % (i, seq, len(seq), sn, so, sr))
    def key(l):
        a, ao, b, bo, sr = l
        sa, sb = (a, 1 if ao == '+' else 0), (b, 0 if bo == '+' else 1)
        return (min(sa, sb), max(sa, sb))
    for a, ao, b, bo, sr in sorted(links, key=key):
        out.append('L\ts%d\t%s\ts%d\t%s\t0M\tSR:i:%d' % (a, ao, b, bo, sr))
    return '\n'.join(out) + '\n'

cases = []
def case(f):
    cases.append(f)
    return f

# ---------------------------------------------------------------- Z1-Z3: in-place inversion, reverse form
def z1(name, form, alt_is_r2):
    c = Case(name, 101)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000)])
    R = c.seq(r2)
    a1 = c.mut(rc(R), 0.01)                       # Z1's allele: ~ rc(r2)
    if form == 'fwd':
        a = c.alt(a1); c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1); w = [r1 + '+', a + '+', r3 + '+']
    else:
        a = c.alt(rc(a1) if alt_is_r2 else c.mut(rc(R), 0.01))
        c.link(r3, '-', a, '+', 1); c.link(a, '+', r1, '-', 1); w = [r1 + '+', a + '-', r3 + '+']
    inv = form == 'fwd' or alt_is_r2
    c.walk('HG1', w, c.seq(r1) + (rc(R) if inv else R) + c.seq(r3))
    c.walk('REF', [r1 + '+', r2 + '+', r3 + '+'])
    links = [(1, '+', 2, '+', 0), (2, '+', 3, '+', 0)] + ([(1, '+', 2, '-', 1), (2, '-', 3, '+', 1)] if inv else [])
    c.files['expected.gfa'] = gfa_expected([x for x in c.nodes if x[0] != 4], links)
    return c

@case
def Z1(): return z1('Z1', 'fwd', False)
@case
def Z2(): return z1('Z2', 'rev', True)
@case
def Z3(): return z1('Z3', 'rev', False)

@case
def Z4():
    c = Case('Z4', 404)
    r = c.ref([c.rnd(2000), c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(2000)])
    a = c.alt(c.mut(c.seq(r[2]), 0.02))
    c.link(r[1], '+', a, '+', 1); c.link(a, '+', r[3], '+', 1)
    c.walk('HG1', [r[0] + '+', r[1] + '+', a + '+', r[3] + '+', r[4] + '+'], ''.join(c.seq(x) for x in r))
    return c

# ---------------------------------------------------------------- Z5-Z9: segment orders
def segcase(name, seed, lens, order, mutate=0.005):
    """reference r1 2k, r2 = the concatenated segments, r3 2k; one alt = the segments in `order`"""
    def build(s):
        c = Case(name, s)
        segs = {k: c.rnd(L) for k, L in zip('ABC', lens)}
        R = ''.join(segs[k] for k in 'ABC'[:len(lens)])
        off, o = {}, 0
        for k in 'ABC'[:len(lens)]: off[k] = (o, o + len(segs[k])); o += len(segs[k])
        r1, r2, r3 = c.ref([c.rnd(2000), R, c.rnd(2000)])
        alt, after, parts, q = '', '', [], 0
        for tok in order.split():
            k, inv = tok[-1], tok.startswith('r')
            s = rc(segs[k]) if inv else segs[k]
            parts.append((q, q + len(s), off[k][0], off[k][1], '-' if inv else '+'))
            alt += c.mut(s, mutate); after += s; q += len(s)
        if not junctions_ok(alt, R, parts): return None
        a = c.alt(alt)
        c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
        c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + after + c.seq(r3))
        c.walk('REF', [r1 + '+', r2 + '+', r3 + '+'])
        # the same chain, injected (exact segment ends: minimap2 may leave ~100 bp unaligned at a junction)
        c.chain('s1..s3:ref:s1+>s4+>s3+', [(qs, qe, st, 2000 + ts, 2000 + te, '%dM' % (qe - qs)) for qs, qe, ts, te, st in parts])
        return c
    return search(build, seed)

@case
def Z5(): return segcase('Z5', 500, (20000, 15000, 25000), 'A rB C')
@case
def Z6(): return segcase('Z6', 600, (20000, 3000, 25000), 'A rB C')
@case
def Z7(): return segcase('Z7', 700, (5000, 30000, 5000), 'A rB C')
@case
def Z8(): return segcase('Z8', 800, (20000, 15000, 25000), 'rC B rA')

@case
def Z9():
    # rc(T1).rc(T2) in place of T1.T2: only the larger segment (rc(T2)) zips; rc(T1) stays alt
    def build(s):
        c = Case('Z9', s)
        T1, T2 = c.rnd(10000), c.rnd(20000)
        r1, r2, r3 = c.ref([c.rnd(2000), T1 + T2, c.rnd(2000)])
        alt = c.mut(rc(T1), 0.004) + c.mut(rc(T2), 0.004)
        if not junctions_ok(alt, T1 + T2, [(0, 10000, 0, 10000, '-'), (10000, 30000, 10000, 30000, '-')]): return None
        a = c.alt(alt)
        c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
        c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + alt[:10000] + rc(T2) + c.seq(r3))
        # Z9b: the two-segment chain injected: refused by (D)
        c.chain('s1..s3:ref:s1+>s4+>s3+', [(0, 10000, '-', 2000, 12000, '10000M'), (10000, 30000, '-', 12000, 32000, '20000M')])
        return c
    return search(build, 900)

# ---------------------------------------------------------------- Z10-Z12: indels inside a copy
def indel_case(name, seed, kind):
    def build(s):
        c = Case(name, s)
        R = c.rnd(20000)
        r1, r2, r3 = c.ref([c.rnd(2000), R, c.rnd(2000)])
        ins = c.rnd(300)
        if kind == 'inv-ins':          # Z10: an inversion split into two '-' records by a 300 bp insertion
            alt = c.mut(rc(R)[:10000], 0.005) + ins + c.mut(rc(R)[10000:], 0.005)
            segs = [(0, 10000, 10000, 20000, '-'), (10300, 20300, 0, 10000, '-')]
            after = rc(R)[:10000] + ins + rc(R)[10000:]
        elif kind == 'fwd-ins':        # Z11: a 300 bp insertion inside a 20 kb in-place copy
            alt = c.mut(R[:10000], 0.005) + ins + c.mut(R[10000:], 0.005)
            segs = [(0, 10000, 0, 10000, '+'), (10300, 20300, 10000, 20000, '+')]
            after = R[:10000] + ins + R[10000:]
        else:                          # Z12: a 300 bp deletion inside a copy
            alt = c.mut(R[:10000], 0.005) + c.mut(R[10300:], 0.005)
            segs = [(0, 10000, 0, 10000, '+'), (10000, 19700, 10300, 20000, '+')]
            after = R[:10000] + R[10300:]
        if not junctions_ok(alt, R, segs): return None
        a = c.alt(alt)
        c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
        c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + after + c.seq(r3))
        c.walk('REF', [r1 + '+', r2 + '+', r3 + '+'])
        c.ins = ins
        return c
    return search(build, seed)

@case
def Z10(): return indel_case('Z10', 1010, 'inv-ins')
@case
def Z11(): return indel_case('Z11', 1110, 'fwd-ins')
@case
def Z12(): return indel_case('Z12', 1210, 'del')

@case
def Z13():
    # IGK shape: HG1's inverted run a1 a2 a3 (3 bp of junk first); HG2 enters it at a2 through a 4 bp node e
    def build(s):
        c = Case('Z13', s)
        R = c.rnd(36010)
        r1, r2, r3 = c.ref([c.rnd(2000), R, c.rnd(2000)])
        RR = rc(R)
        j3, e4 = c.rnd(3), c.rnd(4)
        q1 = j3 + c.mut(RR[:12000], 0.003) + RR[12000:12010] + c.mut(RR[12010:], 0.003)
        q2 = e4 + RR[12000:12010] + c.mut(RR[12010:], 0.003)
        if not junctions_ok(q1, R, [(3, 36013, 0, 36010, '-')]) or not junctions_ok(q2, R, [(4, 24014, 0, 24010, '-')]): return None
        a1 = c.alt(q1[:12003], 'HG1#1#c1', 1, 0)
        a2 = c.alt(q1[12003:12013], 'HG1#1#c1', 1, 12003)
        a3 = c.alt(q1[12013:], 'HG1#1#c1', 1, 12013)
        e = c.alt(e4, 'HG2#1#c2', 2, 0)
        c.link(r1, '+', a1, '+', 1); c.link(a1, '+', a2, '+', 1); c.link(a2, '+', a3, '+', 1); c.link(a3, '+', r3, '+', 1)
        c.link(r1, '+', e, '+', 2); c.link(e, '+', a2, '+', 2)
        c.walk('HG1', [r1 + '+', a1 + '+', a2 + '+', a3 + '+', r3 + '+'], c.seq(r1) + j3 + RR + c.seq(r3))
        c.walk('HG2', [r1 + '+', e + '+', a2 + '+', a3 + '+', r3 + '+'], c.seq(r1) + e4 + RR[12000:] + c.seq(r3))
        return c
    return search(build, 1300)

# ---------------------------------------------------------------- Z14-Z15: no reference target
@case
def Z14():
    # tandem copy beside its source (I), after a 200 bp / 2 kb spacer (I), and back-link tandems (BK):
    # five contigs of one graph, one site each
    c = Case('Z14', 1400)
    for k, (spacer, bk) in enumerate(((0, False), (200, False), (2000, False), (200, True), (2000, True))):
        sn = 'REF#0#chr%d' % (k + 1)
        ids = c.ref([c.rnd(2000), c.rnd(5000)] + ([c.rnd(spacer)] if spacer else []) + [c.rnd(2000)], sn=sn)
        t = c.alt(c.mut(c.seq(ids[1]), 0.01), 'HG1#1#c%d' % (k + 1), 1, 0)
        if bk: c.link(ids[-2], '+', t, '+', 1); c.link(t, '+', ids[1], '+', 1)
        else: c.link(ids[-2], '+', t, '+', 1); c.link(t, '+', ids[-1], '+', 1)
    return c

@case
def Z15():
    # 16p11.2 shape: s1939+ s6141- s1940+ (a tandem copy of the r1|r2 junction, walked '-')
    c = Case('Z15', 1500)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(13000), c.rnd(2000)])
    x = c.alt(c.mut(rc(c.seq(r2)[:12878] + c.seq(r1)[-117:]), 0.0085))      # walked '-': r2's start, then r1's end
    c.link(r1, '+', x, '-', 1); c.link(x, '-', r2, '+', 1)
    c.walk('HG1', [r1 + '+', x + '-', r2 + '+', r3 + '+'])
    return c

# ---------------------------------------------------------------- Z16-Z17: copy number
@case
def Z16():
    c = Case('Z16', 1600)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000)])
    U = c.seq(r2)
    cp = [c.mut(U, 0.01) for _ in range(3)]
    a = c.alt(''.join(cp))
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + U + cp[1] + cp[2] + c.seq(r3))
    c.rest = cp[1] + cp[2]
    return c

def array_case(name, seed, unit, n_ref, n_alt, ndiff=None):
    c = Case(name, seed)
    U = c.rnd(unit)
    cp = lambda: c.mutn(U, ndiff, margin=0) if ndiff else c.mut(U, 0.002, margin=0)
    W = [cp() for _ in range(n_ref)]
    A = [cp() for _ in range(n_alt)]
    r1, r2, r3 = c.ref([c.rnd(2000), ''.join(W), c.rnd(2000)])
    a = c.alt(''.join(A))
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    m = min(n_ref, n_alt)
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + ''.join(W[:m]) + ''.join(A[m:]) + c.seq(r3))
    return c

@case
def Z17a(): return array_case('Z17a', 1701, 6000, 3, 2)
@case
def Z17b(): return array_case('Z17b', 1702, 6000, 3, 5)
@case
def Z17c(): return array_case('Z17c', 1703, 2500, 13, 13, ndiff=1)

@case
def Z18():
    # dispersed copy whose homolog lies outside its window (in r4)
    c = Case('Z18', 1800)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(8000), c.rnd(2000)])
    a = c.alt(c.mut(c.seq(r[3])[1000:7000], 0.01))
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[2], '+', 1)
    c.walk('HG1', [r[0] + '+', a + '+', r[2] + '+'])
    return c

# ---------------------------------------------------------------- Z19: pooled shape, creator and GafWalks
def gaf_line(c, qname, walk):
    steps, plen = [], 0
    for h in walk:
        i, s, sn, so, sr = c.nodes[int(h[1:-1]) - 1]
        steps.append('%s%s:%d-%d' % ('>' if h[-1] == '+' else '<', sn, so, so + len(s)))
        plen += len(s)
    return '\t'.join(map(str, [qname, plen, 0, plen, '+', ''.join(steps), plen, 0, plen, plen, plen, 60, 'tp:A:P', 'cg:Z:%d=' % plen]))

@case
def Z19():
    # H1 inverts r2 (a); H2 inserts c at 8000 (r2 c r3); H3 walks r1 e c a r3.  The graph also holds the
    # pooled path r2 c a r3 (an insertion through a), which no haplotype walks: strict Allowed(a) is
    # empty, so a is infeasible; observed walks (a GAF without the pooled path) open it.
    c = Case('Z19', 1900)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000)])
    a = c.alt(c.mut(rc(c.seq(r2)), 0.01), 'HG1#1#c1', 1)
    cc = c.alt(c.rnd(100), 'HG2#1#c2', 2)
    e = c.alt(c.rnd(100), 'HG3#1#c3', 3)
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    c.link(r2, '+', cc, '+', 2); c.link(cc, '+', r3, '+', 2)
    c.link(r1, '+', e, '+', 3); c.link(e, '+', cc, '+', 3); c.link(cc, '+', a, '+', 3)
    W = {'HG1#1#c1': [r1 + '+', a + '+', r3 + '+'], 'HG2#1#c2': [r2 + '+', cc + '+', r3 + '+'],
         'HG3#1#c3': [r1 + '+', e + '+', cc + '+', a + '+', r3 + '+']}
    c.files['walks.gaf'] = ''.join(gaf_line(c, q, w) + '\n' for q, w in W.items())
    c.files['pooled.gaf'] = c.files['walks.gaf'] + gaf_line(c, 'HG4#1#c4', [r2 + '+', cc + '+', a + '+', r3 + '+']) + '\n'
    R = c.seq(r2)
    c.walk('HG1', W['HG1#1#c1'], c.seq(r1) + rc(R) + c.seq(r3))
    c.walk('HG2', W['HG2#1#c2'])
    c.walk('HG3', W['HG3#1#c3'], c.seq(r1) + c.seq(e) + c.seq(cc) + rc(R) + c.seq(r3))
    return c

# ---------------------------------------------------------------- Z20-Z22: shared nodes, loops, overlaps
@case
def Z20():
    # 16p13.11 shape: an inverted copy s shared by two walks; HG2 adds a forward copy f after it
    c = Case('Z20', 2000)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(10000), c.rnd(2000)])
    R = c.seq(r2)
    s = c.alt(c.mut(rc(R), 0.01), 'HG1#1#c1', 1)
    f = c.alt(c.mut(R, 0.002), 'HG2#1#c2', 2)
    c.link(r1, '+', s, '+', 1); c.link(s, '+', r3, '+', 1)
    c.link(s, '+', f, '+', 2); c.link(f, '+', r3, '+', 2)
    c.walk('HG1', [r1 + '+', s + '+', r3 + '+'], c.seq(r1) + rc(R) + c.seq(r3))
    c.walk('HG2', [r1 + '+', s + '+', f + '+', r3 + '+'], c.seq(r1) + rc(R) + c.seq(f) + c.seq(r3))
    return c

@case
def Z21():
    # GRCh38 16p13.11 shape: x2 is also walked inside an inversion loop r3- x2+ r3+
    c = Case('Z21', 2100)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(12000), c.rnd(2000)])
    R = c.seq(r2)
    x1, x2 = c.chain_run([c.mut(R[:6000], 0.004), c.mut(R[6000:], 0.004)], 'HG1#1#c1', 1, r1, r3)
    c.link(r3, '-', x2, '+', 2)
    c.walk('HG1', [r1 + '+', x1 + '+', x2 + '+', r3 + '+'], c.seq(r1) + R[:6000] + c.seq(x2) + c.seq(r3))
    c.walk('loop', [r3 + '-', x2 + '+', r3 + '+'])
    return c

@case
def Z21b():
    # the whole allele inside the loop: Allowed blocked, infeasible, never aligned
    c = Case('Z21b', 2101)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(12000), c.rnd(2000)])
    x = c.alt(c.mut(c.seq(r2), 0.004))
    c.link(r1, '+', x, '+', 1); c.link(x, '+', r3, '+', 1); c.link(r3, '-', x, '+', 2)
    c.walk('HG1', [r1 + '+', x + '+', r3 + '+'])
    return c

@case
def Z22():
    # overlapping inversions of two haplotypes (T11)
    c = Case('Z22', 2200)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(6000), c.rnd(6000), c.rnd(2000)])
    S2, S3, S4 = c.seq(r[1]), c.seq(r[2]), c.seq(r[3])
    a = c.alt(c.mut(rc(S2 + S3), 0.004), 'HG1#1#c1', 1)
    b = c.alt(c.mut(rc(S3 + S4), 0.006), 'HG2#1#c2', 2)
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[3], '+', 1)
    c.link(r[1], '+', b, '+', 2); c.link(b, '+', r[4], '+', 2)
    c.walk('HG1', [r[0] + '+', a + '+', r[3] + '+', r[4] + '+'], c.seq(r[0]) + rc(S2 + S3) + S4 + c.seq(r[4]))
    c.walk('HG2', [r[0] + '+', r[1] + '+', b + '+', r[4] + '+'], c.seq(r[0]) + S2 + rc(S3 + S4) + c.seq(r[4]))
    return c

@case
def Z23():
    # leak geometry: the second hand-written snarl is not a snarl of this graph (b.R links back to r3.L)
    c = Case('Z23', 2300)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(2000), c.rnd(2000)], snarl=False)
    a = c.alt(c.rnd(6000), 'HG1#1#c1', 1)
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[2], '+', 1)
    b = c.alt(c.rnd(3000), 'HG2#1#c2', 2)
    c.link(r[2], '+', b, '+', 2); c.link(b, '+', r[4], '+', 2); c.link(b, '+', r[2], '+', 3)
    c.snarls = [(r[0], r[2]), (r[2], r[4])]
    return c

# ---------------------------------------------------------------- Z25-Z27
@case
def Z25():
    # split inversion with a 600 bp forward centre: the contig is r1 rc(B2) M rc(B1) r5 over reference
    # r1 B1 M B2 r5; each half has no homology in its own window
    def build(s):
        c = Case('Z25', s)
        B1, M, B2 = c.rnd(10000), c.rnd(600), c.rnd(8000)
        r = c.ref([c.rnd(2000), B1, M, B2, c.rnd(2000)])
        h1, h2 = c.mut(rc(B2), 0.003), c.mut(rc(B1), 0.003)
        if not junctions_ok(h1 + M + h2, B1 + M + B2, [(0, 8000, 10600, 18600, '-'), (8600, 18600, 0, 10000, '-')]): return None
        a1 = c.alt(h1, 'HG1#1#c1', 1, 0)
        a2 = c.alt(h2, 'HG1#1#c1', 1, 8600)
        c.link(r[0], '+', a1, '+', 1); c.link(a1, '+', r[2], '+', 1)
        c.link(r[2], '+', a2, '+', 1); c.link(a2, '+', r[4], '+', 1)
        c.walk('HG1', [r[0] + '+', a1 + '+', r[2] + '+', a2 + '+', r[4] + '+'])
        return c
    return search(build, 2500)

@case
def Z26():
    # an excursion inside a native inversion, walked r4- x+ r2-
    c = Case('Z26', 2600)
    r = c.ref([c.rnd(2000), c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(2000)])
    S3 = c.seq(r[2])
    x = c.alt(c.mut(rc(S3), 0.005))
    c.link(r[3], '-', x, '+', 1); c.link(x, '+', r[1], '-', 1)
    c.walk('HG1', [r[3] + '-', x + '+', r[1] + '-'], rc(c.seq(r[3])) + rc(S3) + rc(c.seq(r[1])))
    return c

def z27(name, form):
    # projection across a 10 bp D and a 10 bp I at piece junctions (injected chain)
    c = Case(name, 2700)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000)])
    R = c.seq(r2)
    ins = c.rnd(10)
    x1s, x2s, x3s = R[:2000], R[2010:4000], ins + R[4000:]
    cg = '2000M10D1990M10I2000M'
    if form == 'fwd':
        x = c.chain_run([x1s, x2s, x3s], 'HG1#1#c1', 1, r1, r3)
        walk = [r1 + '+'] + [n + '+' for n in x] + [r3 + '+']
        c.chain('s1..s3:ref:s1+>s4+,s5+,s6+>s3+', [(0, 6000, '+', 2000, 8000, cg)])
    elif form == 'rev':
        # the same allele written by a reverse-strand contig: run y1 y2 y3 = rc(x3) rc(x2) rc(x1)
        y = [c.alt(rc(x3s), 'HG1#1#c1', 1, 0), c.alt(rc(x2s), 'HG1#1#c1', 1, 2010), c.alt(rc(x1s), 'HG1#1#c1', 1, 4000)]
        c.link(r3, '-', y[0], '+', 1); c.link(y[0], '+', y[1], '+', 1); c.link(y[1], '+', y[2], '+', 1); c.link(y[2], '+', r1, '-', 1)
        walk = [r1 + '+', y[2] + '-', y[1] + '-', y[0] + '-', r3 + '+']
        c.chain('s1..s3:ref:s1+>s6-,s5-,s4->s3+', [(0, 6000, '+', 2000, 8000, cg)])
    else:
        # walked forward but inverted: a '-' chain
        z = c.chain_run([rc(x3s), rc(x2s), rc(x1s)], 'HG1#1#c1', 1, r1, r3)
        walk = [r1 + '+'] + [n + '+' for n in z] + [r3 + '+']
        c.chain('s1..s3:ref:s1+>s4+,s5+,s6+>s3+', [(0, 6000, '-', 2000, 8000, cg)])
    c.walk('HG1', walk, c.seq(r1) + (rc(R) if form == 'minus' else R) + c.seq(r3))
    return c

@case
def Z27a(): return z27('Z27a', 'fwd')
@case
def Z27b(): return z27('Z27b', 'rev')
@case
def Z27c(): return z27('Z27c', 'minus')

@case
def Z27e():
    # a whole 5 bp node inside a 5 bp insertion has no aligned target base: it stays alt
    c = Case('Z27e', 2728)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000)])
    R = c.seq(r2)
    ins = c.rnd(5)
    x = c.chain_run([R[:3000], ins, R[3000:]], 'HG1#1#c1', 1, r1, r3)
    c.chain('s1..s3:ref:s1+>s4+,s5+,s6+>s3+', [(0, 6005, '+', 2000, 8000, '3000M5I3000M')])
    c.walk('HG1', [r1 + '+'] + [n + '+' for n in x] + [r3 + '+'], c.seq(r1) + R[:3000] + ins + R[3000:] + c.seq(r3))
    return c

# ---------------------------------------------------------------- Z28: aligner failures (25 one-window sites)
@case
def Z28():
    c = Case('Z28', 2800)
    mark = c.rnd(40)
    prev = c.node(c.rnd(2000), REF, 0, 0); so = 2000
    for k in range(25):
        w = c.rnd(10000)
        if k == 7: w = w[:5000] + mark + w[5040:]
        wn = c.node(w, REF, so, 0); so += len(w)
        nx = c.node(c.rnd(2000), REF, so, 0); so += 2000
        c.link(prev, '+', wn, '+', 0); c.link(wn, '+', nx, '+', 0)
        a = c.alt(c.mut(w, 0.01), 'HG%03d#1#c%d' % (k + 1, k + 1), k + 1)
        c.link(prev, '+', a, '+', k + 1); c.link(a, '+', nx, '+', k + 1)
        c.snarls.append((prev, nx))
        prev = nx
    c.files['mark'] = mark + '\n'
    return c

@case
def Z29():
    # two chains; a debug flag drops one translated link of the second
    c = Case('Z29', 2900)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(6000), c.rnd(2000)])
    S2, S3 = c.seq(r[1]), c.seq(r[2])
    a = c.alt(c.mut(rc(S2), 0.004), 'HG1#1#c1', 1)
    b = c.alt(c.mut(rc(S3), 0.004), 'HG2#1#c2', 2)
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[2], '+', 1)
    c.link(r[1], '+', b, '+', 2); c.link(b, '+', r[3], '+', 2)
    c.walk('HG1', [r[0] + '+', a + '+', r[2] + '+'], c.seq(r[0]) + rc(S2) + S3)
    c.walk('HG2', [r[1] + '+', b + '+', r[3] + '+'])
    return c

# ---------------------------------------------------------------- Z31: determinism (three contigs)
@case
def Z31():
    c = Case('Z31', 3100)
    def inv(sn, contig, sr, lens, order):
        segs = [c.rnd(L) for L in lens]
        r = c.ref([c.rnd(2000), ''.join(segs), c.rnd(2000)], sn=sn)
        alt = ''.join(c.mut(rc(segs[i]) if inv_ else segs[i], 0.004) for i, inv_ in order)
        a = c.alt(alt, contig, sr)
        c.link(r[0], '+', a, '+', sr); c.link(a, '+', r[2], '+', sr)
    inv('REF#0#chr1', 'HG1#1#c1', 1, (8000, 6000, 6000), [(0, False), (1, True), (2, False)])
    inv('REF#0#chr2', 'HG1#1#c2', 1, (20000,), [(0, True)])
    inv('REF#0#chr3', 'HG2#1#c3', 2, (7000, 9000), [(1, True), (0, True)])
    return c

# ---------------------------------------------------------------- Z34-Z37: rule U and the gates
def flip(name, k):
    # U at 10 kb and rc(U') at 60 kb of a 100 kb window, U' = U with 12 substitutions; the alt carries
    # U' 's base at k of the 12 sites: '+' onto U and '-' onto rc(U') score within 2 bases
    c = Case(name, 3400)
    U = c.rnd(6000)
    pos = sorted(c.r.sample(range(100, 5900), 12))
    Up = list(U)
    for p in pos: Up[p] = c.r.choice([x for x in 'ACGT' if x != U[p]])
    Up = ''.join(Up)
    w = c.rnd(100000)
    w = w[:10000] + U + w[16000:60000] + rc(Up) + w[66000:]
    V = list(U)
    for p in pos[:k]: V[p] = Up[p]
    r1, r2, r3 = c.ref([c.rnd(2000), w, c.rnd(2000)])
    a = c.alt(c.rnd(1000) + ''.join(V) + c.rnd(1000))
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'])
    return c

@case
def Z34a(): return flip('Z34a', 5)
@case
def Z34b(): return flip('Z34b', 7)

def hitch(name, alone):
    # window 120 kb holds B = [20k, 40k) and a 4 kb element P = [90k, 94k); the alt is rc(B') (20 kb),
    # 10 kb novel, P' (1.2% divergent), 5 kb novel
    def build(s):
        c = Case(name, s)
        w = c.rnd(120000)
        B, P = w[20000:40000], w[90000:94000]
        head = c.rnd(20000) if alone else c.mut(rc(B), 0.003)
        alt = head + c.rnd(10000) + c.mut(P, 0.012, margin=0) + c.rnd(5000)
        if not alone and not junctions_ok(alt, w, [(0, 20000, 20000, 40000, '-')]): return None
        r1, r2, r3 = c.ref([c.rnd(2000), w, c.rnd(2000)])
        a = c.alt(alt)
        c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
        c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + (alt if alone else rc(B) + alt[20000:]) + c.seq(r3))
        return c
    return search(build, 3500)

@case
def Z35(): return hitch('Z35', False)
@case
def Z35b(): return hitch('Z35b', True)

@case
def Z36():
    # SMN shape, scaled: two distinct 20 kb copies 150 kb apart in one window; the alt is 0.24% closer to
    # the right one (80 and 32 substitutions): a copy tie, resolved by score
    c = Case('Z36', 3600)
    U = c.rnd(20000)
    L, Rt = c.mutn(U, 80), c.mutn(U, 32)
    w = c.rnd(30000) + L + c.rnd(150000) + Rt + c.rnd(30000)
    r1, r2, r3 = c.ref([c.rnd(2000), w, c.rnd(2000)])
    a = c.alt(U)
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + Rt + c.seq(r3))
    return c

@case
def Z37():
    # VNTR expansion aligned out of phase: 13 insertions of 80 bp in 14.5 kb
    c = Case('Z37', 3700)
    W = c.rnd(14560)
    alt = ''
    for k in range(14):
        alt += c.mut(W[k * 1040:(k + 1) * 1040], 0.003, margin=0)
        if k < 13: alt += c.rnd(80)
    r1, r2, r3 = c.ref([c.rnd(2000), W, c.rnd(2000)])
    a = c.alt(alt)
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'])
    return c

# ---------------------------------------------------------------- Z38-Z39 (injected chains)
@case
def Z38():
    # two chains share node s, projected 33 bp apart; HG2's private neighbour b is re-anchored
    c = Case('Z38', 3800)
    r1, r2, r3 = c.ref([c.rnd(2000), c.rnd(20000), c.rnd(2000)])
    R = c.seq(r2)
    a, s = c.chain_run([c.mut(R[:10000], 0.002), c.mut(R[10000:], 0.002)], 'HG1#1#c1', 1, r1, r3)
    b = c.alt(c.mut(R[100:10033], 0.004), 'HG2#1#c2', 2)
    c.link(r1, '+', b, '+', 2); c.link(b, '+', s, '+', 2)
    c.chain('s1..s3:ref:s1+>s4+,s5+>s3+', [(0, 20000, '+', 2000, 22000, '20000M')])
    c.chain('s1..s3:ref:s1+>s6+,s5+>s3+', [(0, 19933, '+', 2100, 22000, '9933M33I9967M')])
    c.walk('HG1', [r1 + '+', a + '+', s + '+', r3 + '+'], c.seq(r1) + R + c.seq(r3))
    c.walk('HG2', [r1 + '+', b + '+', s + '+', r3 + '+'], c.seq(r1) + R[100:] + c.seq(r3))
    return c

@case
def Z39():
    # a 3 bp first piece; a rank-0 boundary 30 bp inside the block start: no inward snap
    c = Case('Z39', 3900)
    r1, r2a, r2b, r3 = c.ref([c.rnd(2000), c.rnd(130), c.rnd(5870), c.rnd(2000)])
    R = c.seq(r2a) + c.seq(r2b)
    x, y = c.chain_run([R[100:103], c.mut(R[103:], 0.003)], 'HG1#1#c1', 1, r1, r3)
    c.chain('s1..s4:ref:s1+>s5+,s6+>s4+', [(0, 5900, '+', 2100, 8000, '5900M')])
    c.walk('HG1', [r1 + '+', x + '+', y + '+', r3 + '+'], c.seq(r1) + R[100:] + c.seq(r3))
    return c

def z39ctl(name, start):
    # positive controls: a block start 3 bp from the window start / from a node boundary snaps to it
    c = Case(name, 3940)
    r1, r2a, r2b, r3 = c.ref([c.rnd(2000), c.rnd(130), c.rnd(5870), c.rnd(2000)])
    R = c.seq(r2a) + c.seq(r2b)
    off = start - 2000
    x = c.alt(c.mut(R[off:], 0.003))
    c.link(r1, '+', x, '+', 1); c.link(x, '+', r3, '+', 1)
    L = 6000 - off
    c.chain('s1..s4:ref:s1+>s5+>s4+', [(0, L, '+', start, 8000, '%dM' % L)])
    c.walk('HG1', [r1 + '+', x + '+', r3 + '+'], c.seq(r1) + (R if start == 2003 else R[130:]) + c.seq(r3))
    return c

@case
def Z39b(): return z39ctl('Z39b', 2003)
@case
def Z39c(): return z39ctl('Z39c', 2127)

@case
def Z42():
    # identical copies X, one outside Allowed(a) (another haplotype enters a from s3)
    c = Case('Z42', 4200)
    X = c.rnd(10000)
    r = c.ref([c.rnd(2000), c.rnd(5000) + X + c.rnd(5000), c.rnd(1000), c.rnd(5000) + X + c.rnd(5000), c.rnd(2000)])
    a = c.alt(c.mut(X, 0.003))
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[4], '+', 1); c.link(r[2], '+', a, '+', 2)
    c.walk('HG1', [r[0] + '+', a + '+', r[4] + '+'], c.seq(r[0]) + X + c.seq(r[4]))
    return c

@case
def Z43():
    # a satellite window (50 kb of N, then a 30 kb tandem array of a 60 bp unit) with 20 queries of
    # 100 kb, each a tandem array of its own variant of the unit: all fail the screen
    c = Case('Z43', 4300)
    base = 'CATTC' * 12
    def array(n, unit):
        s = []
        while len(s) * 60 < n: s.append(c.mut(unit, 0.03, margin=0))
        return ''.join(s)[:n]
    r1, r2, r3 = c.ref([c.rnd(2000), 'N' * 50000 + array(30000, c.mut(base, 0.25, margin=0)), c.rnd(2000)])
    for k in range(20):
        c.chain_run([array(100000, c.mut(base, 0.25, margin=0))], 'HG%03d#1#c%d' % (k + 1, k + 1), k + 1, r1, r3)
    return c

# ---------------------------------------------------------------- Z24, Z40, Z41: alt-vs-alt (v2)
def z24(name, rform, mform, rate, extra=None):
    # parallel insertions between adjacent reference nodes r1 r2 (kind I): R by rank 1, M by rank 3.
    # form 'rev': written by a reverse-strand contig (L r2 - X +, L X + r1 -; walked r1+ X- r2+).
    # extra 'series': a contig walks r1 R M r2 (copies in series); 'nonprivate': M also links to r3
    c = Case(name, 2400)
    refs = c.ref([c.rnd(2000), c.rnd(2000)] + ([c.rnd(2000), c.rnd(2000)] if extra == 'nonprivate' else []))
    r1, r2 = refs[0], refs[1]
    Rw = c.rnd(6000)
    Mw = c.mut(Rw, rate)
    R = c.alt(Rw if rform == 'fwd' else rc(Rw), 'HG1#1#c1', 1)
    M = c.alt(Mw if mform == 'fwd' else rc(Mw), 'HG3#1#c3', 3)
    for node, form, sr in ((R, rform, 1), (M, mform, 3)):
        if form == 'fwd': c.link(r1, '+', node, '+', sr); c.link(node, '+', r2, '+', sr)
        else: c.link(r2, '-', node, '+', sr); c.link(node, '+', r1, '-', sr)
    if extra == 'series': c.link(R, '+', M, '+', 4)
    if extra == 'nonprivate': c.link(M, '+', refs[2], '+', 5)
    tail = [x + '+' for x in refs[1:]]
    merged = extra is None
    c.walk('HG1', [r1 + '+', R + ('+' if rform == 'fwd' else '-')] + tail)
    c.walk('HG3', [r1 + '+', M + ('+' if mform == 'fwd' else '-')] + tail,
           c.seq(r1) + (Rw if merged else Mw) + ''.join(c.seq(x) for x in refs[1:]))
    return c

@case
def Z24(): return z24('Z24', 'fwd', 'fwd', 0.02)
@case
def Z24s(): return z24('Z24s', 'fwd', 'fwd', 0.02, 'series')
@case
def Z24p(): return z24('Z24p', 'fwd', 'fwd', 0.02, 'nonprivate')
@case
def Z24b(): return z24('Z24b', 'rev', 'rev', 0.02)
@case
def Z24c(): return z24('Z24c', 'rev', 'fwd', 0.005)

def z40(name, inv=False, nonprivate=False):
    # nested bubble: rank 1 inserts a b d (one run), rank 5 inserts x between a and d (reusing them);
    # x ~ b (or rc(b)).  The two walks share a and d: the branches x and b are bounded by alt nodes.
    c = Case(name, 4000)
    r1, r2 = c.ref([c.rnd(2000), c.rnd(2000)])
    A, B, D = c.rnd(2000), c.rnd(6000), c.rnd(2000)
    a, b, d = c.chain_run([A, B, D], 'HG1#1#c1', 1, r1, r2)
    x = c.alt(c.mut(rc(B) if inv else B, 0.01), 'HG5#1#c5', 5)
    c.link(a, '+', x, '+', 5); c.link(x, '+', d, '+', 5)
    if nonprivate: c.link(x, '+', r2, '+', 6)
    c.walk('HG1', [r1 + '+', a + '+', b + '+', d + '+', r2 + '+'])
    c.walk('HG5', [r1 + '+', a + '+', x + '+', d + '+', r2 + '+'],
           None if nonprivate else c.seq(r1) + A + (rc(B) if inv else B) + D + c.seq(r2))
    return c

@case
def Z40(): return z40('Z40')
@case
def Z40b(): return z40('Z40b', inv=True)
@case
def Z40p(): return z40('Z40p', nonprivate=True)

@case
def Z41():
    # HP shape (16q22): R is an insertion at p = end(r1); M, a homologous copy, bypasses r3: its window
    # is [p+90, p+114).  Different anchors: never grouped, reported as near-parallel (disjoint windows)
    c = Case('Z41', 4100)
    r1, r2, r3, r4 = c.ref([c.rnd(2000), c.rnd(90), c.rnd(24), c.rnd(2000)])
    S = c.rnd(6000)
    R = c.alt(S, 'HG1#1#c1', 1)
    M = c.alt(c.mut(S, 0.01), 'HG2#1#c2', 2)
    c.link(r1, '+', R, '+', 1); c.link(R, '+', r2, '+', 1)
    c.link(r2, '+', M, '+', 2); c.link(M, '+', r4, '+', 2)
    c.walk('HG1', [r1 + '+', R + '+', r2 + '+', r3 + '+', r4 + '+'])
    c.walk('HG2', [r1 + '+', r2 + '+', M + '+', r4 + '+'])
    return c

@case
def Z24m():
    # one group, two members: ranks 2 and 3 copy rank 1; both merge onto the rank-1 representative
    c = Case('Z24m', 2410)
    r1, r2 = c.ref([c.rnd(2000), c.rnd(2000)])
    S = c.rnd(6000)
    seqs = [S, c.mut(S, 0.02), c.mut(S, 0.02)]
    xs = [c.alt(seqs[k], 'HG%d#1#c%d' % (k + 1, k + 1), k + 1) for k in range(3)]
    for k, x in enumerate(xs):
        c.link(r1, '+', x, '+', k + 1); c.link(x, '+', r2, '+', k + 1)
        c.walk('HG%d' % (k + 1), [r1 + '+', x + '+', r2 + '+'], c.seq(r1) + S + c.seq(r2))
    return c

@case
def Z24r():
    # rounds: rank 1's allele is unrelated, rank 3 copies rank 2: round 1 aligns ranks 2 and 3 to rank 1
    # (no chain); round 2 regroups rank 3 under the next-lowest-SR branch (rank 2) and merges it
    c = Case('Z24r', 2420)
    r1, r2 = c.ref([c.rnd(2000), c.rnd(2000)])
    U, S = c.rnd(6000), c.rnd(6000)
    seqs = [U, S, c.mut(S, 0.02)]
    xs = [c.alt(seqs[k], 'HG%d#1#c%d' % (k + 1, k + 1), k + 1) for k in range(3)]
    for k, x in enumerate(xs):
        c.link(r1, '+', x, '+', k + 1); c.link(x, '+', r2, '+', k + 1)
        c.walk('HG%d' % (k + 1), [r1 + '+', x + '+', r2 + '+'], c.seq(r1) + (S if k == 2 else seqs[k]) + c.seq(r2))
    return c

@case
def Z24t():
    # the lowest-SR allele is 3 kb, shorter than ceil(b*i): it can be the target of no unit, so it is
    # no round's representative, and rank 3 merges onto rank 2 in round 1
    c = Case('Z24t', 2430)
    r1, r2 = c.ref([c.rnd(2000), c.rnd(2000)])
    S = c.rnd(6000)
    seqs = [c.rnd(3000), S, c.mut(S, 0.02)]
    xs = [c.alt(seqs[k], 'HG%d#1#c%d' % (k + 1, k + 1), k + 1) for k in range(3)]
    for k, x in enumerate(xs):
        c.link(r1, '+', x, '+', k + 1); c.link(x, '+', r2, '+', k + 1)
        c.walk('HG%d' % (k + 1), [r1 + '+', x + '+', r2 + '+'], c.seq(r1) + (S if k == 2 else seqs[k]) + c.seq(r2))
    return c

# ---------------------------------------------------------------- regressions from code review
def overlap_case(name, seed, order, parts, middle):
    """segcase's 20/15/25 kb shape with an injected chain whose parts overlap at a strand switch;
    middle(A, B, C, alt) is what the image spells between r1 and r3"""
    c = segcase(name, seed, (20000, 15000, 25000), order)
    R = c.seq('s2')
    A, B, C = R[:20000], R[20000:35000], R[35000:]
    c.inject = []
    c.chain('s1..s3:ref:s1+>s4+>s3+', parts)
    c.walks = []
    c.walk('HG1', ['s1+', 's4+', 's3+'], c.seq('s1') + middle(A, B, C, c.seq('s4')) + c.seq('s3'))
    return c

@case
def Z5o():
    # GC-1, frame F '+' then '-': A's record runs 5 bp past the A|rc(B) junction (into rc(B) on the query,
    # into B on the target), so the '-' part's target overlaps A's at its BOTTOM, which is its query END:
    # it is trimmed there (it used to vanish whole).  A's 5 extra bases spell B[:5]; 5 bp stay alt.
    return overlap_case('Z5o', 500, 'A rB C',
                        [(0, 20005, '+', 2000, 22005, '20005M'), (20000, 35000, '-', 22000, 37000, '15000M'),
                         (35000, 60000, '+', 37000, 62000, '25000M')],
                        lambda A, B, C, alt: A + B[:5] + rc(B[:14995]) + C)

@case
def Z8o():
    # GC-1, frame R '-' then '+': rc(C)'s record runs 3 bp into B; B's '+' part is below rc(C)'s target,
    # overlapping it at its TOP, which is its query END
    return overlap_case('Z8o', 800, 'rC B rA',
                        [(0, 25003, '-', 36997, 62000, '25003M'), (25000, 40000, '+', 22000, 37000, '15000M'),
                         (40000, 60000, '-', 2000, 22000, '20000M')],
                        lambda A, B, C, alt: rc(C) + rc(B[-3:]) + B[3:] + rc(A))

@case
def Z5f():
    # a part whose target lies inside an earlier part's is a fold: reported as dropped "(fold)", its query alt
    return overlap_case('Z5f', 500, 'A rB C',
                        [(0, 20000, '+', 2000, 22000, '20000M'), (20000, 35000, '-', 7000, 22000, '15000M'),
                         (35000, 60000, '+', 37000, 62000, '25000M')],
                        lambda A, B, C, alt: A + alt[20000:35000] + C)

@case
def Z5r():
    # GC-1 with real minimap2: A rc(B) over A B where B = mh core rc(mh), a 12 bp microhomology at the
    # breakpoint: A's '+' record runs 12 bp into rc(B) and overlaps the '-' record's target at its query end
    c = Case('Z5r', 5100)
    A, mh = c.rnd(20000), c.rnd(12)
    B = mh + c.rnd(15000 - 24) + rc(mh)
    r1, r2, r3 = c.ref([c.rnd(2000), A + B, c.rnd(2000)])
    a = c.alt(c.mut(A, 0.003) + c.mut(rc(B), 0.003, margin=100))
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    # the image: A and mh (A's record), rc(B[12:]) (the '-' part, its top snapped 12 bp to the window
    # end: absorbed), then the last 12 bp of the allele, which stay alt
    c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + A + mh + rc(B[12:]) + rc(mh) + c.seq(r3))
    return c

@case
def Z23z():
    # GC-2: Z23 with a zippable first allele (an inverted copy of r2).  The left site's B.L side (r3.L)
    # reaches b, which links to r3.R (its B.R): it leaks too, and both sites are skipped
    c = Case('Z23z', 2300)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(2000), c.rnd(2000)], snarl=False)
    a = c.alt(rc(c.seq(r[1])), 'HG1#1#c1', 1)
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[2], '+', 1)
    b = c.alt(c.rnd(3000), 'HG2#1#c2', 2)
    c.link(r[2], '+', b, '+', 2); c.link(b, '+', r[4], '+', 2); c.link(b, '+', r[2], '+', 3)
    c.snarls = [(r[0], r[2]), (r[2], r[4])]
    return c

def z1x(name, extra):
    # Z1 plus a 300 bp node z on B.L only, as vg puts it inside the snarl s1..s3: a tip (z.R -> s3.L)
    # or a fold-back walked s3- z+ s3+ (GC-2 / RD-2: the unedited site model used to fail, exit 3);
    # or the creator link a -> s3 written a second time, reversed and without SR (GC-3)
    c = z1(name, 'fwd', False)
    if extra == 'dup':
        c.link('s3', '-', 's4', '-')
        return c
    del c.files['expected.gfa']
    z = c.alt(c.rnd(300), 'HG2#1#c2', 2)
    if extra == 'hairpin': c.link('s3', '-', z, '+', 2)
    c.link(z, '+', 's3', '+', 2)
    return c

@case
def Z1t(): return z1x('Z1t', 'tip')
@case
def Z1h(): return z1x('Z1h', 'hairpin')
@case
def Z1d(): return z1x('Z1d', 'dup')

# ---------------------------------------------------------------- regressions from the genome-wide sweep
def ngap_case(name, seed, kind):
    """paf-inconsistent: minimap2 2.30 leaves ambiguous bases out of PAF column 11 inside I and D runs,
    not only in =/X columns.  'del': the allele lacks 106 bp of the window holding 10 Ns (a D run >= G)
    and 30 bp holding 5 Ns (a D run < G, inside a block), and carries 20 Ns facing 20 Ns of the window
    (minimap2 writes them '=': they are never '=' here).  'ins': an inverted allele ('-' strand) with a
    70 bp insertion holding 10 Ns (an I run >= G) and a 20 bp one holding 4 (< G, inside a block).
    The bases next to every deleted or inserted stretch differ from its ends, so no gap can shift."""
    def build(s):
        c = Case(name, s)
        if kind == 'del':
            X = c.rnd(10000)
            X = X[:5000] + 'N' * 20 + X[5020:]
            D1, U = c.rnd(53) + 'N' * 10 + c.rnd(43), c.rnd(4000)
            V, T = c.rnd(12) + 'N' * 5 + c.rnd(13), c.rnd(5864)
            if X[-1] == D1[-1] or U[0] == D1[0] or U[-1] == V[-1] or T[0] == V[0]: return None
            W, alt = X + D1 + U + V + T, X + U + T
        else:
            X, W2a, W2b = c.rnd(10000), c.rnd(5000), c.rnd(5000)
            ins, ins2 = c.rnd(30) + 'N' * 10 + c.rnd(30), c.rnd(8) + 'N' * 4 + c.rnd(8)
            if X[-1] == ins[-1] or W2a[0] == ins[0] or W2a[-1] == ins2[-1] or W2b[0] == ins2[0]: return None
            W, alt = X + W2a + W2b, rc(X + ins + W2a + ins2 + W2b)
        r1, r2, r3 = c.ref([c.rnd(2000), W, c.rnd(2000)])
        a = c.alt(alt)
        c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
        # the image: the reference under the blocks; a long gap is a deletion edge (del) or stays alt
        # (ins); a short one is absorbed
        mid = X + U + V + T if kind == 'del' else rc(X + ins + W2a + W2b)
        c.walk('HG1', [r1 + '+', a + '+', r3 + '+'], c.seq(r1) + mid + c.seq(r3))
        return c
    return search(build, seed)

@case
def Zn1(): return ngap_case('Zn1', 5500, 'del')
@case
def Zn2(): return ngap_case('Zn2', 5600, 'ins')

@case
def Zg():
    # rule (G), v3: HG1's record walks a (~ rc(r2)) inside an ordinary excursion over r2, folds back at
    # r4 through y (r4+ y+ r4-) and reads r2 again: zipping a would collapse its haplotype's second copy.
    # Observed Allowed sees one excursion at a time; only the whole record shows it.
    c = Case('Zg', 5700)
    r = c.ref([c.rnd(2000), c.rnd(6000), c.rnd(2000), c.rnd(2000), c.rnd(2000)], snarl=False)
    a = c.alt(c.mut(rc(c.seq(r[1])), 0.01), 'HG1#1#c1', 1)
    c.link(r[0], '+', a, '+', 1); c.link(a, '+', r[2], '+', 1)
    y = c.alt(c.rnd(300), 'HG2#1#c2', 2)
    c.link(r[3], '+', y, '+', 2); c.link(y, '+', r[3], '-', 2)
    c.snarls = [(r[0], r[2]), (r[2], r[4])]
    straight = [r[0] + '+', a + '+', r[2] + '+', r[3] + '+', r[4] + '+']
    loop = [r[0] + '+', a + '+', r[2] + '+', r[3] + '+', y + '+', r[3] + '-', r[2] + '-', r[1] + '-', r[0] + '-']
    c.files['loop.gaf'] = gaf_line(c, 'HG1#1#c1', loop) + '\n'
    c.files['straight.gaf'] = gaf_line(c, 'HG1#1#c1', straight) + '\n'
    c.walk('HG1', straight, c.seq(r[0]) + rc(c.seq(r[1])) + ''.join(c.seq(x) for x in r[2:]))
    return c

@case
def Ztr():
    # 520 copies of a 100 bp unit and one alt copy, aligned to every copy by 520 injected records: more
    # feasible sub-records than the 500 the chain keeps (run with -b 100 --min-piece 100)
    c = Case('Ztr', 5800)
    U = c.rnd(100)
    r1, r2, r3 = c.ref([c.rnd(2000), U * 520, c.rnd(2000)])
    a = c.alt(U)
    c.link(r1, '+', a, '+', 1); c.link(a, '+', r3, '+', 1)
    for k in range(520): c.inject.append('s1..s3:ref:0\trecord\t0\t100\t+\t%d\t%d\t100M' % (2000 + 100 * k, 2100 + 100 * k))
    return c

def main():
    os.makedirs(T, exist_ok=True)
    for f in cases:
        c = f()
        c.write()
        if c.name == 'Z1':
            c.write('_nosr', no_sr=True)                  # Z30: L lines without SR
            c.write('_nosn', drop_sn='s4')                # Z30: an S line without SN
        if c.name == 'Z40b':
            c.write('_shuf', shuffle=True)                # alt-vs-alt determinism
        if c.name == 'Z31':
            c.write('_shuf', shuffle=True)
            for k in (1, 2, 3):
                c2 = Case('Z31', 3100); c2.nodes, c2.links, c2.snarls = c.nodes, c.links, c.snarls
                c2.write('_chr%d' % k, only_sn=['#chr%d' % k, '#c%d' % k])

main()
PYEOF

cat > $T/q.py <<'PYEOF'
#!/usr/bin/env python3
# rgfa-zip.t queries on case D (the directory $T/D: in.gfa, walks.tsv) and one run R of it
# (D/R.gfa, D/R.tsv, D/R.dump/).  Each prints one line, compared to an exact string by the test.
#   S D R                 output segment names, by id
#   L D R                 output links 'a+b-' as written, sorted
#   newL D R              output links that join different input sides than any input link, as written
#   extra D R             output links off the reference path (all but ++ links between abutting rank-0
#                         segments of one SN), as written, sorted
#   row D R CONTIG COL..  columns of the report row(s) of CONTIG's unit (first owner), space-joined
#   col D R COL [PASS]    one column of every report row (of one pass), space-joined
#   walks D R             'ok' if every walk spells its sequence in the input and its predicted
#                         sequence in the output; else the labels that fail
#   cuts D R NODE         offsets at which input node NODE was cut
#   same D R              'same' if the output has the input's S and L lines (as sets), else 'differ'
#   fscc D R              non-trivial SCCs of the forward-strand graph (rank 0 '+', alt '+')
#   back D R              rank-0 back edges
#   alt D R               alt segments of the output as name:SN:SO:LN
#   removed D R           alt bp of the input minus alt bp of the output
#   pieces D R UNIT [ST]  the zipped pieces of UNIT in the --dump (state ST, default new): ta-tb
#   units D R COL         one column of every unit of the --dump
#   canon GFA MAXID       the graph with new pieces (id > MAXID) renamed by (SN, SO, LN), one line per
#                         record, sorted (for per-SN against whole-genome comparisons)
#   ident D LABEL MIN     identity of walk LABEL's input sequence and the sequence its image spells
#                         after the edit (same length: substitutions only), then '>=MIN' or '<MIN'
import sys, os, collections, glob

COMP = str.maketrans('ACGTN', 'TGCAN')
def rc(s): return s.translate(COMP)[::-1]
FL = {'+': '-', '-': '+'}

class G:
    def __init__(self, path):
        self.S, self.tags, self.L, self.Slines = collections.OrderedDict(), {}, [], {}
        for l in open(path):
            t = l.rstrip('\n').split('\t')
            if t[0] == 'S':
                self.S[t[1]] = t[2]
                self.Slines[t[1]] = '\t'.join(t)
                self.tags[t[1]] = dict((x[:2], x[5:]) for x in t[3:])
            elif t[0] == 'L':
                self.L.append((t[1], t[2], t[3], t[4], '\t'.join(t)))
        self.adj = collections.defaultdict(list)
        for a, ao, b, bo, _ in self.L:
            self.adj[(a, ao)].append((b, bo))
            self.adj[(b, FL[bo])].append((a, FL[ao]))
    def hseq(self, h): return self.S[h[0]] if h[1] == '+' else rc(self.S[h[0]])
    def spells(self, seq):
        """a walk of this graph that spells exactly seq, from a node start to a node end"""
        stack = []
        for n, s in self.S.items():
            for o in '+-':
                hs = s if o == '+' else rc(s)
                if len(hs) <= len(seq) and seq.startswith(hs): stack.append(((n, o), len(hs)))
        seen = set()
        while stack:
            h, pos = stack.pop()
            if (h, pos) in seen: continue
            seen.add((h, pos))
            if pos == len(seq): return True
            for h2 in self.adj.get(h, []):
                hs = self.hseq(h2)
                if pos + len(hs) <= len(seq) and seq.startswith(hs, pos): stack.append((h2, pos + len(hs)))
        return False
    def walk_spell(self, walk):
        """the sequence of an input walk 's1+,s4-', or None if a step has no link"""
        hs = [(w[:-1], w[-1]) for w in walk.split(',')]
        for x, y in zip(hs, hs[1:]):
            if y not in self.adj.get(x, []): return None
        return ''.join(self.hseq(h) for h in hs)
    def keys(self):
        """links as sets of sides (a link equals its reverse)"""
        out = set()
        for a, ao, b, bo, _ in self.L:
            sa = (a, 'R' if ao == '+' else 'L'); sb = (b, 'L' if bo == '+' else 'R')
            out.add(frozenset([sa, sb]))
        return out

def report(path):
    cols, rows = None, []
    for l in open(path):
        l = l.rstrip('\n')
        if l.startswith('#site'): cols = l[1:].split('\t'); continue
        if l.startswith('#'): continue
        rows.append(dict(zip(cols, l.split('\t'))))
    return rows

def main():
    a = sys.argv[1:]
    cmd = a[0]
    if cmd == 'canon':
        g, maxid = G(a[1]), int(a[2])
        name = {}
        for n in g.S:
            t = g.tags[n]
            name[n] = n if int(n[1:]) <= maxid else 'new:%s:%s:%s' % (t['SN'], t['SO'], t['LN'])
        lines = ['S\t%s\t%s' % (name[n], '\t'.join(g.Slines[n].split('\t')[2:])) for n in g.S]
        for x, xo, y, yo, l in g.L:
            r1, r2 = (name[x], xo, name[y], yo), (name[y], FL[yo], name[x], FL[xo])
            lines.append('L\t%s\t%s\t%s\t%s\t%s' % (min(r1, r2) + ('\t'.join(l.split('\t')[5:]),)))
        print('\n'.join(sorted(lines)))
        return
    D, R = os.path.join(os.environ.get('T', '.'), a[1]), a[2]
    out_path = os.path.join(D, R + '.gfa')
    if cmd == 'S':
        g = G(out_path)
        print(' '.join(sorted(g.S, key=lambda n: int(n[1:]))))
    elif cmd == 'L':
        g = G(out_path)
        print(' '.join(sorted('%s%s%s%s' % (x, xo, y, yo) for x, xo, y, yo, _ in g.L)))
    elif cmd == 'newL':
        g, i = G(out_path), G(os.path.join(D, 'in.gfa'))
        ik = i.keys()
        res = []
        for x, xo, y, yo, _ in g.L:
            sa = (x, 'R' if xo == '+' else 'L'); sb = (y, 'L' if yo == '+' else 'R')
            if frozenset([sa, sb]) not in ik: res.append('%s%s%s%s' % (x, xo, y, yo))
        print(' '.join(sorted(res)))
    elif cmd == 'extra':
        g = G(out_path)
        res = []
        for x, xo, y, yo, _ in g.L:
            tx, ty = g.tags[x], g.tags[y]
            if xo == yo and tx['SR'] == '0' and ty['SR'] == '0' and tx['SN'] == ty['SN']:
                a_, b_ = (x, y) if xo == '+' else (y, x)
                if int(g.tags[a_]['SO']) + len(g.S[a_]) == int(g.tags[b_]['SO']): continue
            res.append('%s%s%s%s' % (x, xo, y, yo))
        print(' '.join(sorted(res)))
    elif cmd == 'row':
        contig, want = a[3], a[4:]
        rows = [r for r in report(os.path.join(D, R + '.tsv')) if r['pass'] == 'ref' and r['owners'].split(',')[0].split(':')[0] == contig]
        print(' | '.join(' '.join(r[c] for c in want) for r in rows))
    elif cmd == 'col':
        rows = report(os.path.join(D, R + '.tsv'))
        if len(a) > 4: rows = [r for r in rows if r['pass'] == a[4]]
        print(' '.join(r[a[3]] for r in rows))
    elif cmd == 'walks':
        g, i = G(out_path), G(os.path.join(D, 'in.gfa'))
        bad = []
        for l in open(os.path.join(D, 'walks.tsv')):
            label, walk, before, after = l.rstrip('\n').split('\t')
            if i.walk_spell(walk) != before: bad.append(label + ':before')
            elif not g.spells(after): bad.append(label + ':after')
        print(' '.join(bad) if bad else 'ok')
    elif cmd == 'cuts':
        g, i = G(out_path), G(os.path.join(D, 'in.gfa'))
        t = i.tags[a[3]]
        so, L = int(t['SO']), len(i.S[a[3]])
        pcs = sorted(int(x['SO']) - so for n, x in g.tags.items() if x['SN'] == t['SN'] and so <= int(x['SO']) < so + L)
        print(' '.join(str(p) for p in pcs if p != 0))
    elif cmd == 'same':
        g, i = G(out_path), G(os.path.join(D, 'in.gfa'))
        print('same' if set(g.Slines.values()) == set(i.Slines.values()) and g.keys() == i.keys() and len(g.L) == len(i.L) else 'differ')
    elif cmd == 'fscc':
        g = G(out_path)
        adj = collections.defaultdict(list)
        for x, xo, y, yo, _ in g.L:
            if xo == '+' and yo == '+': adj[x].append(y)
            if xo == '-' and yo == '-': adj[y].append(x)
        idx, low, st, on, cnt, n = {}, {}, [], set(), [0], [0]
        sys.setrecursionlimit(1000000)
        def sc(v):
            idx[v] = low[v] = cnt[0]; cnt[0] += 1; st.append(v); on.add(v)
            for w in adj[v]:
                if w not in idx: sc(w); low[v] = min(low[v], low[w])
                elif w in on: low[v] = min(low[v], idx[w])
            if low[v] == idx[v]:
                comp = []
                while True:
                    w = st.pop(); on.discard(w); comp.append(w)
                    if w == v: break
                if len(comp) > 1 or v in adj[v]: n[0] += 1
        for v in list(g.S):
            if v not in idx: sc(v)
        print(n[0])
    elif cmd == 'back':
        g = G(out_path)
        n = 0
        for x, xo, y, yo, _ in g.L:
            tx, ty = g.tags[x], g.tags[y]
            if tx['SR'] != '0' or ty['SR'] != '0' or tx['SN'] != ty['SN']: continue
            rx, ry = xo == '+', yo == '-'
            if rx == ry: continue
            Rn, Ln = (x, y) if rx else (y, x)
            if int(g.tags[Ln]['SO']) < int(g.tags[Rn]['SO']) + len(g.S[Rn]): n += 1
        print(n)
    elif cmd == 'alt':
        g = G(out_path)
        print(' '.join('%s:%s:%s:%d' % (n, g.tags[n]['SN'], g.tags[n]['SO'], len(g.S[n])) for n in g.S if g.tags[n]['SR'] != '0'))
    elif cmd == 'removed':
        g, i = G(out_path), G(os.path.join(D, 'in.gfa'))
        altbp = lambda x: sum(len(x.S[n]) for n in x.S if x.tags[n]['SR'] != '0')
        print(altbp(i) - altbp(g))
    elif cmd == 'pieces':
        st = a[4] if len(a) > 4 else 'new'
        res = []
        for f in sorted(glob.glob(os.path.join(D, R + '.dump', 'edit', '*.tsv'))):
            for l in open(f):
                t = l.rstrip('\n').split('\t')
                if t[0] == 'piece' and t[1] == a[3] and t[10] == st: res.append('%s:%s-%s' % (t[2], t[7], t[8]))
        print(' '.join(res))
    elif cmd == 'units':
        rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(D, R + '.dump', 'units.tsv'))]
        c = rows[0].index(a[3])
        print(' '.join(r[c] for r in rows[1:]))
    elif cmd == 'ident':
        mn = float(a[3])
        for l in open(os.path.join(D, 'walks.tsv')):
            label, walk, before, after = l.rstrip('\n').split('\t')
            if label != R: continue
            same = sum(1 for x, y in zip(before, after) if x == y)
            idn = same / float(max(len(before), len(after)))
            print('%.4f %s%s' % (idn, '>=' if idn >= mn else '<', a[3]))
    else:
        sys.exit('unknown query ' + cmd)

main()
PYEOF

q() { python3 $T/q.py "$@"; }
# zip CASE RUN [options]: rgfa-zip on CASE's graph; outputs CASE/RUN.{gfa,tsv,err,dump}, exit code in CASE/RUN.exit
zip() {
    local c=$1 r=$2
    shift 2
    rm -rf $T/$c/$r.dump $T/$c/$r.gfa $T/$c/$r.tsv
    rgfa-zip -m $MM2 --tmpdir $T/tmp --dump $T/$c/$r.dump -o $T/$c/$r.gfa -r $T/$c/$r.tsv "$@" $T/$c/in.gfa $T/$c/snarls.json 2> $T/$c/$r.err
    echo $? > $T/$c/$r.exit
}
ex() { cat $T/$1/$2.exit; }
inj() { echo --inject-chains $T/$1/inject.tsv; }
none() { [ -e $T/$1/$2.gfa ] || [ -e $T/$1/$2.tsv ] && echo written || echo nothing; }

# ---- Z1-Z4: an allele that is an inverted or forward copy of its window
zip Z1 run
is "$(ex Z1 run)/$(cmp -s $T/Z1/run.gfa $T/Z1/expected.gfa && echo same)" "0/same" "Z1: exact S/L lines: a deleted, L s1 + s2 - and L s2 - s3 + added, 0 splits"
is "$(q row Z1 run HG1#1#c1 outcome label kept_bp)" "zipped INV 6000" "Z1: report row zipped INV 6000"
is "$(q walks Z1 run)" "ok" "Z1: the haplotype walk spells r1 rc(r2) r3, the reference walk r1 r2 r3"
zip Z2 run
is "$(ex Z2 run)/$(cmp -s $T/Z2/run.gfa $T/Z1/run.gfa && echo same)" "0/same" "Z2: minigraph's reverse form (L r3 - a +, L a + r1 -) gives output byte-identical to Z1"
is "$(q walks Z2 run)" "ok" "Z2: walks spell"
zip Z3 run
is "$(ex Z3 run)/$(cmp -s $T/Z3/run.gfa $T/Z3/expected.gfa && echo same)" "0/same" "Z3: exact S/L lines: the reverse-written forward copy is deleted, no new link"
is "$(q row Z3 run HG1#1#c1 outcome label kept_bp)/$(q walks Z3 run)" "zipped FWD 6000/ok" "Z3: zipped FWD 6000; walks spell"
zip Z4 run
is "$(q row Z4 run HG1#1#c1 outcome label kept_bp)/$(q S Z4 run)/$(q extra Z4 run)" "zipped FWD 6000/s1 s2 s3 s4 s5/" "Z4: a 2% divergent copy attached to interior nodes is zipped; no link off the reference"
is "$(q walks Z4 run)" "ok" "Z4: substitutions absorbed: the walk spells the reference"

# ---- Z5-Z9: segment orders (real alignments, then the same chains injected for exact cut points)
for c in Z5 Z6 Z7 Z8 Z9; do zip $c run; zip $c inj $(inj $c); done
is "$(q row Z5 run HG1#1#c1 outcome frame label)/$(q walks Z5 run)/$(q back Z5 run)/$(q fscc Z5 run)" "zipped F +-+/ok/0/0" "Z5: A rc(B) C (20/15/25 kb): frame F '+-+', walks spell, no back edge, no forward cycle"
is "$(q row Z5 inj HG1#1#c1 outcome kept_bp)/$(q cuts Z5 inj s2)/$(q extra Z5 inj)/$(q S Z5 inj)" "zipped 60000/20000 35000/s5+s6- s6-s7+/s1 s3 s5 s6 s7" "Z5 (injected): r2 cut at 20k and 35k, exactly 2 new links, new pieces s5-s7 in offset order"
is "$(q walks Z5 inj)" "ok" "Z5 (injected): walks spell"
zip Z5 idb $(inj Z5) --id-base 1000
is "$(q S Z5 idb)" "s1 s3 s1000 s1001 s1002" "Z5: --id-base 1000 names the new pieces s1000-s1002"
zip Z5 idbad $(inj Z5) --id-base 4
is "$(ex Z5 idbad)/$(none Z5 idbad)" "2/nothing" "Z5: --id-base not above the largest input id: exit 2, nothing written"
is "$(q row Z6 run HG1#1#c1 outcome frame label kept_bp)/$(q walks Z6 run)" "zipped F +-+ 48000/ok" "Z6: |B| = 3 kb is not an island (no hairpin rule): '+-+' zipped"
is "$(q cuts Z6 inj s2)/$(q extra Z6 inj)" "20000 23000/s5+s6- s6-s7+" "Z6 (injected): r2 cut at 20k and 23k, 2 new links"
is "$(q row Z7 run HG1#1#c1 outcome frame label)/$(q walks Z7 run)" "zipped F +-+/ok" "Z7: |B| > |A| + |C|: '+-+' kept"
is "$(q cuts Z7 inj s2)/$(q extra Z7 inj)" "5000 35000/s5+s6- s6-s7+" "Z7 (injected): r2 cut at 5k and 35k, 2 new links"
is "$(q row Z8 run HG1#1#c1 outcome frame label)/$(q walks Z8 run)/$(q fscc Z8 run)" "zipped R -+-/ok/0" "Z8: rc(C) B rc(A): frame R '-+-'; 0 forward SCCs"
is "$(q cuts Z8 inj s2)/$(q extra Z8 inj)/$(q fscc Z8 inj)" "20000 35000/s1+s7- s3-s5+ s5+s6- s6-s7+/0" "Z8 (injected): r2 cut at 20k and 35k, 4 inversion links, 0 forward SCCs"
is "$(q row Z9 run HG1#1#c1 outcome label kept_bp)/$(q alt Z9 run)/$(q fscc Z9 run)" "zipped INV 20000/s7:HG1#1#c1:0:10000/0" "Z9: rc(T1) rc(T2): only the larger segment zips, rc(T1) stays alt; 0 forward SCCs"
is "$(q walks Z9 run)" "ok" "Z9: walks spell"
is "$(q row Z9 inj HG1#1#c1 outcome)/$(q same Z9 inj)" "(D)/same" "Z9: both segments injected need a back edge: refused by (D), output = input"

# ---- Z10-Z13: indels inside a copy; IGK
zip Z10 run
is "$(q row Z10 run HG1#1#c1 outcome label kept_bp)/$(q cuts Z10 run s2)/$(q alt Z10 run)" "zipped INV 20000/10000/s7:HG1#1#c1:10000:300" "Z10: an inversion split by a 300 bp insertion: one '-' segment; the 300 bp stays alt at its projected point"
is "$(q walks Z10 run)" "ok" "Z10: walks spell"
zip Z11 run
is "$(q row Z11 run HG1#1#c1 outcome label kept_bp)/$(q cuts Z11 run s2)/$(q extra Z11 run)" "zipped FWD 20000/10000/s5+s7+ s7+s6+" "Z11: a 300 bp insertion in a copy is kept as an alt bubble"
is "$(q walks Z11 run)" "ok" "Z11: walks spell (the insertion included)"
zip Z12 run
is "$(q row Z12 run HG1#1#c1 outcome kept_bp)/$(q cuts Z12 run s2)/$(q extra Z12 run)/$(q alt Z12 run)" "zipped 19700/10000 10300/s5+s7+/" "Z12: a 300 bp deletion becomes a deletion edge of exactly 300 bp"
is "$(q walks Z12 run)" "ok" "Z12: walks spell"
zip Z13 run
is "$(q row Z13 run HG1#1#c1 outcome label kept_bp) | $(q row Z13 run HG2#1#c2 outcome label kept_bp)" "zipped INV 36010 | zipped INV 0" "Z13: IGK shape: both chains accepted; the second removes nothing of its own (kept_bp 0: its pieces were zipped by the first)"
is "$(q pieces Z13 run 1 new)/$(q pieces Z13 run 1 compatible)" "/s5:26000-26010 s6:2000-26000" "Z13: the second haplotype's pieces are all compatible and add nothing"
is "$(q alt Z13 run)/$(q extra Z13 run)" "s7:HG2#1#c2:0:4 s10:HG1#1#c1:0:3/s1+s10+ s1+s7+ s3-s8+ s7+s8- s9+s10-" "Z13: the 4 bp entry re-attaches at the projected point (s7+ s8-), the 3 bp junk stays"
is "$(q fscc Z13 run)/$(q back Z13 run)/$(q walks Z13 run)" "0/0/ok" "Z13: no loop, no back edge; both walks spell"

# ---- Z14-Z15: no reference target
zip Z14 run
is "$(q units Z14 run kind)/$(q col Z14 run outcome)/$(q same Z14 run)" "I I I BK BK/no-window no-window no-window no-window no-window/same" "Z14: tandem copies (I) and back-link tandems (BK, 200 bp and 2 kb spacers) are untouched; 0 cuts"
zip Z15 run
is "$(q units Z15 run kind)/$(q same Z15 run)" "I/same" "Z15: 16p11.2 shape (r1+ x- r2+) is kind I: untouched, no link added"

# ---- Z16-Z18: copy number, dispersed copies
zip Z16 run
is "$(q row Z16 run HG1#1#c1 outcome kept_bp)/$(q cuts Z16 run s2)/$(q alt Z16 run)/$(q extra Z16 run)" "zipped 6000//s5:HG1#1#c1:6000:12000/s2+s5+ s5+s3+" "Z16: 3 copies over 1: one zipped left-anchored, the other two stay one insertion at the window end"
is "$(q walks Z16 run)" "ok" "Z16: walks spell (copy number 3 kept)"
for c in Z17a Z17b Z17c; do zip $c run; done
is "$(q row Z17a run HG1#1#c1 outcome kept_bp)/$(q cuts Z17a run s2)/$(q extra Z17a run)/$(q walks Z17a run)" "zipped 12000/12000/s5+s3+/ok" "Z17: 2 units over 3: units 1-2 zipped plus a deletion edge over unit 3"
is "$(q row Z17b run HG1#1#c1 outcome kept_bp)/$(q alt Z17b run)/$(q extra Z17b run)/$(q walks Z17b run)" "zipped 18000/s5:HG1#1#c1:18000:12000/s2+s5+ s5+s3+/ok" "Z17: 5 units over 3: units 1-3 zipped, 4-5 kept as one insertion at the window end"
is "$(q row Z17c run HG1#1#c1 outcome kept_bp)/$(q S Z17c run)/$(q walks Z17c run)" "zipped 32500/s1 s2 s3/ok" "Z17: a 13-copy array with 1-base differences zips in place"
zip Z18 run
is "$(q row Z18 run HG1#1#c1 outcome)/$(q same Z18 run)" "prefiltered/same" "Z18: a dispersed copy whose homolog lies outside its window is untouched"

# ---- Z19: the pooled shape, creator walks and GafWalks
zip Z19 creator
is "$(q units Z19 creator outcome)/$(q same Z19 creator)" "infeasible small infeasible/same" "Z19: pooled shape: strict Allowed(a) is empty, a's unit is infeasible; 0 cuts"
zip Z19 gaf --walks gaf:$T/Z19/walks.gaf
is "$(q row Z19 gaf HG1#1#c1 outcome label)/$(q newL Z19 gaf)/$(q S Z19 gaf)" "zipped INV/s1+s2- s2+s5- s2-s3+/s1 s2 s3 s5 s6" "Z19: GafWalks without the pooled path: a zipped as an inversion (its links translated onto r2-)"
is "$(q walks Z19 gaf)" "ok" "Z19: GafWalks: walks spell"
zip Z19 pooled --walks gaf:$T/Z19/pooled.gaf
is "$(ex Z19 pooled)/$(q same Z19 pooled)" "0/same" "Z19: GafWalks with a haplotype on the pooled path: infeasible again"

# ---- Z20-Z23
zip Z20 run
is "$(q row Z20 run HG1#1#c1 outcome kept_bp) | $(q row Z20 run HG2#1#c2 outcome dropped)" "zipped 10000 | trimmed-below-b s5[0-10000):(C)" "Z20: 16p13.11 shape: exactly one zipped; the forward copy is refused by (C)"
is "$(q alt Z20 run)/$(q walks Z20 run)" "s5:HG2#1#c2:0:10000/ok" "Z20: the forward copy stays alt; walks spell"
zip Z21 run
is "$(q row Z21 run HG1#1#c1 outcome kept_bp feasible_bp)/$(q alt Z21 run)/$(q walks Z21 run)" "zipped 6000 6000/s5:HG1#1#c1:6000:6000/ok" "Z21: the node also walked inside an inversion loop is refused by (B) (blocked before chaining)"
zip Z21b run
is "$(q row Z21b run HG1#1#c1 outcome)/$(q same Z21b run)" "infeasible/same" "Z21: a whole allele inside the loop is infeasible and untouched"
zip Z22 run
is "$(q row Z22 run HG1#1#c1 outcome label) | $(q row Z22 run HG2#1#c2 outcome label)" "zipped INV | zipped INV" "Z22: overlapping inversions of two haplotypes: both zipped"
is "$(q extra Z22 run)/$(q fscc Z22 run)/$(q walks Z22 run)" "s1+s3- s2+s4- s2-s4+ s3-s5+/0/ok" "Z22: 4 inversion links (handle cycles allowed), 0 forward SCCs; walks spell"
zip Z23 run
is "$(ex Z23 run)/$(q col Z23 run outcome site)/$(q col Z23 run query_bp site)/$(q same Z23 run)" "0/site:leak site:leak/9000 3000/same" "Z23: leak geometry: the site is skipped with a site:leak row carrying its alt bp, and so is its neighbour, whose B.L side reaches the leaking node; output = input"

# ---- Z25-Z27
zip Z25 run
is "$(q col Z25 run outcome)/$(q same Z25 run)" "prefiltered prefiltered split-inversion-candidate/same" "Z25: split inversion with a 600 bp forward centre: not zipped; split-inversion-candidate"
is "$(q col Z25 run label diag)/$(q col Z25 run query_bp diag)" "INV/18600" "Z25: the merged query (both halves and the centre) has one '-' segment"
zip Z26 run
is "$(q units Z26 run dep) $(q units Z26 run walk) $(q units Z26 run arr) $(q units Z26 run kind)" "s2+ s6- s4+ F" "Z26: an excursion walked inside a native inversion (r4- x+ r2-) is flipped to the canonical frame"
is "$(q row Z26 run HG1#1#c1 outcome label)/$(q extra Z26 run)/$(q S Z26 run)/$(q walks Z26 run)" "zipped FWD//s1 s2 s3 s4 s5/ok" "Z26: zipped, no new link; spellable"
for c in Z27a Z27b Z27c Z27e; do zip $c run $(inj $c); done
is "$(q pieces Z27a run 0)/$(q S Z27a run)/$(q extra Z27a run)" "s4:2000-4000 s5:4000-6000 s6:6000-8000/s1 s2 s3/" "Z27: projection across a D and an I at piece junctions: pieces abut exactly; no cut, no link"
is "$(q pieces Z27b run 0)/$(q S Z27b run)/$(q extra Z27b run)" "s6:2000-4000 s5:4000-6000 s4:6000-8000/s1 s2 s3/" "Z27: the same allele written by a reverse-strand contig"
is "$(q pieces Z27c run 0)/$(q cuts Z27c run s2)/$(q extra Z27c run)" "s4:6000-8000 s5:4000-6000 s6:2000-4000//s1+s2- s2-s3+" "Z27: walked forward but inverted: no cut, 2 inversion links only"
is "$(q walks Z27a run)/$(q walks Z27b run)/$(q walks Z27c run)/$(q back Z27a run)/$(q back Z27c run)" "ok/ok/ok/0/0" "Z27: walks spell; no 1-bp back edge"
is "$(q row Z27e run HG1#1#c1 outcome dropped)/$(q alt Z27e run)/$(q cuts Z27e run s2)" "zipped s5[0-5):(empty)/s5:HG1#1#c1:3000:5/3000" "Z27: a 5 bp node inside a 5 bp insertion has no aligned target base: it stays alt"

# ---- Z28: aligner failures, with a fake minimap2 that misbehaves on the window holding a marker
cat > $T/fake_mm2.sh <<'SHEOF'
#!/bin/bash
if [ "$1" = "--version" ]; then exec "$FAKE_REAL" --version; fi
tfa=""
for a in "$@"; do case "$a" in */t.fa) tfa=$a;; */t.mmi) tfa=${a%.mmi}.fa;; esac; done
hit=0
if [ -n "$tfa" ] && grep -q "$FAKE_MARK" "$tfa"; then hit=1; fi
idx=0; for a in "$@"; do [ "$a" = "-d" ] && idx=1; done
case "$FAKE_MODE" in
  fail-one) if [ $hit = 1 ]; then echo "fake: failing on purpose" >&2; exit 1; fi ;;
  fail-all) echo "fake: failing on purpose" >&2; exit 1 ;;
  kill-once) if [ $hit = 1 ] && [ $idx = 0 ] && [ ! -e "$FAKE_STATE/killed" ]; then touch "$FAKE_STATE/killed"; kill -9 $$; fi ;;
  kill-always) if [ $hit = 1 ] && [ $idx = 0 ]; then kill -9 $$; fi ;;
  # RD-1: an index (-d) whose write silently failed: minimap2 exits 0 and leaves it empty
  empty-index) if [ $idx = 1 ]; then touch "$FAKE_STATE/indexed"; "$FAKE_REAL" "$@"; st=$?
                 for a in "$@"; do case "$a" in *.mmi) : > "$a";; esac; done; exit $st; fi ;;
  # RD-4: the marked window's PAF goes to a full disk (minimap2: "failed to write the results", exit 1)
  enospc) if [ $hit = 1 ] && [ $idx = 0 ]; then exec "$FAKE_REAL" "$@" > /dev/full; fi ;;
  # RD-7: the marked window is killed every time; every other alignment takes a while
  kill-slow) if [ $hit = 1 ] && [ $idx = 0 ]; then kill -9 $$; fi; sleep 4 ;;
  # RD-8: alignments that never finish (a stop signal has to end them)
  hang) if [ $idx = 0 ]; then exec -a "rgfa-zip-test-hang-$FAKE_TAG" sleep 61; fi ;;
  # RD-10: during the run, something makes the report path a non-empty directory
  mkdir-report) mkdir -p "$FAKE_REPORT/x" ;;
  # logs: every alignment holds its slot 0.5 s longer and reports known minimap2 times as its last line
  timing) if [ $idx = 0 ]; then "$FAKE_REAL" "$@"; st=$?; sleep 0.5
            echo "[M::main] Real time: 1.250 sec; CPU: 0.750 sec; Peak RSS: 0.010 GB" >&2; exit $st; fi ;;
esac
exec "$FAKE_REAL" "$@"
SHEOF
chmod +x $T/fake_mm2.sh
export FAKE_REAL=$MM2 FAKE_MARK=$(cat $T/Z28/mark) FAKE_STATE=$T/Z28
zip Z28 clean
is "$(ex Z28 clean)/$(q col Z28 clean outcome | tr ' ' '\n' | sort | uniq -c | tr -s ' ' | tr '\n' ';')" "0/ 25 zipped;" "Z28: the clean run zips all 25 one-window sites"
FAKE_MODE=fail-one zip Z28 failone -m $T/fake_mm2.sh
is "$(ex Z28 failone)/$(q col Z28 failone outcome | tr ' ' '\n' | sort | uniq -c | tr -s ' ' | tr '\n' ';')/$(grep -c 'failing on purpose' $T/Z28/failone.err)" "0/ 1 aligner-failed; 24 zipped;/1" "Z28: exit 1 on one window (4%): that unit aligner-failed, exit 0"
FAKE_MODE=fail-all zip Z28 failall -m $T/fake_mm2.sh
is "$(ex Z28 failall)/$(none Z28 failall)" "4/nothing" "Z28: exit 1 on every window: exit 4, nothing written"
rm -f $T/Z28/killed
FAKE_MODE=kill-once zip Z28 killonce -m $T/fake_mm2.sh
is "$(ex Z28 killonce)/$(cmp -s $T/Z28/killonce.gfa $T/Z28/clean.gfa && cmp -s $T/Z28/killonce.tsv $T/Z28/clean.tsv && echo same)/$(grep -c 'retrying it alone' $T/Z28/killonce.err)" "0/same/1" "Z28: SIGKILL once: the window is retried alone; output identical to the clean run"
FAKE_MODE=kill-always zip Z28 killalways -m $T/fake_mm2.sh
is "$(ex Z28 killalways)/$(none Z28 killalways)" "4/nothing" "Z28: SIGKILL always: exit 4, nothing written"
is "$(ls $T/tmp | wc -l)" "0" "Z28: no temporary files left behind, even after exit 4"
# review fixes (robustness): resource and I/O failures never become quiet outcomes
FAKE_MODE=empty-index zip Z28 emptyidx -m $T/fake_mm2.sh
is "$(ex Z28 emptyidx)/$(cmp -s $T/Z28/emptyidx.gfa $T/Z28/clean.gfa && cmp -s $T/Z28/emptyidx.tsv $T/Z28/clean.tsv && echo same)/$([ -e $T/Z28/indexed ] && echo indexed || echo no-index)" "0/same/no-index" "RD-1: no index file is written (minimap2 -d does not check its writes; an empty index made every query 'prefiltered'): output identical to the clean run"
FAKE_MODE=enospc zip Z28 enospc -m $T/fake_mm2.sh
is "$(ex Z28 enospc)/$(none Z28 enospc)/$(grep -c 'could not write or read its files in --tmpdir' $T/Z28/enospc.err)/$(ls $T/tmp | wc -l)" "5/nothing/1/0" "RD-4: minimap2 cannot write its PAF (full disk): exit 5, nothing written, not an aligner-failed row"
win0=$(awk '$1 == "S" && $2 == "s2" {print substr($3, 5001, 40)}' $T/Z28/in.gfa)
t0=$(date +%s%N)
FAKE_MODE=kill-slow FAKE_MARK=$win0 zip Z28 killslow -m $T/fake_mm2.sh -t 3 -j 3
t1=$(date +%s%N)
is "$(ex Z28 killslow)/$(none Z28 killslow)/$(ls $T/tmp | wc -l)" "4/nothing/0" "RD-7: the first window is killed twice while the others align slowly: exit 4, nothing written"
ok $(( (t1 - t0) / 1000000 < 12000 )) "RD-7: the failure stops the other sites and their aligners at once ($(( (t1 - t0) / 1000000 )) ms; 16 s when the sites in flight ran on)"
FAKE_MODE=mkdir-report FAKE_REPORT=$T/Z28/mkrep.tsv zip Z28 mkrep -m $T/fake_mm2.sh
is "$(ex Z28 mkrep)/$([ -e $T/Z28/mkrep.gfa ] && echo graph || echo no-graph)/$(grep -c 'rename' $T/Z28/mkrep.err)" "5/no-graph/1" "RD-10: the report cannot be renamed into place after the graph was: the graph is removed again, exit 5"
rm -rf $T/Z28/mkrep.tsv
FAKE_MODE=hang FAKE_TAG=$$ rgfa-zip -m $T/fake_mm2.sh -t 2 -j 2 --tmpdir $T/tmp -o $T/Z28/hang.gfa -r $T/Z28/hang.tsv $T/Z28/in.gfa $T/Z28/snarls.json 2> $T/Z28/hang.err &
hpid=$!
sleep 2
kill -TERM $hpid
wait $hpid
echo $? > $T/Z28/hang.exit
sleep 0.3
is "$(ex Z28 hang)/$(none Z28 hang)/$(ls $T/tmp | wc -l)/$(pgrep -f "rgfa-zip-test-hang-$$" | wc -l)/$(grep -c 'stopped by signal 15' $T/Z28/hang.err)" "143/nothing/0/0/1" "RD-8: SIGTERM: the aligners are stopped, the temporary directory removed, nothing written; the exit status is the signal's"

# ---- Z29: a debug flag drops one translated link of the second of two chains
zip Z29 clean
is "$(q col Z29 clean outcome)" "zipped zipped" "Z29: without the flag both chains zip"
RGFA_ZIP_DEBUG_DROP_LINK=s1..s4:1 zip Z29 drop
is "$(q row Z29 drop HG1#1#c1 outcome) | $(q row Z29 drop HG2#1#c2 outcome)" "zipped | reverted:V2" "Z29: the chain with the dropped link is reverted:V2; the other commits"
is "$(q S Z29 drop)/$(q walks Z29 drop)" "s1 s2 s3 s4 s6/ok" "Z29: a deleted, b untouched; walks spell"
RGFA_ZIP_DEBUG_DROP_LINK=s1..s4:1 zip Z29 strict --strict
is "$(ex Z29 strict)/$(none Z29 strict)" "3/nothing" "Z29: --strict: exit 3, nothing written"

# ---- Z30: L lines without SR (witness fallback, the same edits as Z1); an S line without SN
zip Z1_nosr run
is "$(ex Z1_nosr run)/$(q S Z1_nosr run)/$(q L Z1_nosr run)/$(q col Z1_nosr run source)" "0/s1 s2 s3/$(q L Z1 run)/witness" "Z30: L lines without SR: the witness fallback makes the same edit as Z1"
zip Z1_nosn run
is "$(ex Z1_nosn run)/$(none Z1_nosn run)" "2/nothing" "Z30: an S line without SN: exit 2, nothing written"

# ---- Z31: determinism (three contigs, three sites)
zip Z31 t1 -t 1 -j 1
zip Z31 t3 -t 3 -j 3
zip Z31_shuf t2 -t 2 -j 2
is "$(q col Z31 t1 outcome)" "zipped zipped zipped" "Z31: the three chains zip"
is "$(cmp -s $T/Z31/t1.gfa $T/Z31/t3.gfa && cmp -s $T/Z31/t1.tsv $T/Z31/t3.tsv && echo same)" "same" "Z31: -t 1 -j 1 and -t 3 -j 3: byte-identical graph and report"
is "$(cmp -s $T/Z31/t1.gfa $T/Z31_shuf/t2.gfa && cmp -s $T/Z31/t1.tsv $T/Z31_shuf/t2.tsv && echo same)" "same" "Z31: shuffled S/L/snarl lines at -t 2 -j 2: byte-identical graph and report"
maxid() { awk '$1=="S" {n = substr($2, 2) + 0; if (n > m) m = n} END {print m}' $1; }
for k in 1 2 3; do zip Z31_chr$k run; done
is "$(for k in 1 2 3; do q canon $T/Z31_chr$k/run.gfa $(maxid $T/Z31_chr$k/in.gfa); done | LC_ALL=C sort | md5sum)" "$(q canon $T/Z31/t1.gfa $(maxid $T/Z31/in.gfa) | LC_ALL=C sort | md5sum)" "Z31: per-SN runs equal the whole-input run after canonical renaming of new pieces"
is "$(for k in 1 2 3; do grep -v '^#' $T/Z31_chr$k/run.tsv; done | md5sum)" "$(grep -v '^#' $T/Z31/t1.tsv | md5sum)" "Z31: per-SN report rows equal the whole-input report rows"

# ---- Z32: a second run on the output zips nothing
for c in Z1 Z5 Z10 Z13 Z16 Z22; do
    mkdir -p $T/Z32_$c && cp $T/$c/run.gfa $T/Z32_$c/in.gfa && cp $T/$c/snarls.json $T/Z32_$c/
    zip Z32_$c run
done
is "$(for c in Z1 Z5 Z10 Z13 Z16 Z22; do echo -n "$(ex Z32_$c run):$(q col Z32_$c run outcome | tr ' ' '\n' | grep -c '^zipped$'):$(q same Z32_$c run) "; done)" "0:0:same 0:0:same 0:0:same 0:0:same 0:0:same 0:0:same " "Z32: run on its own output (Z1, Z5, Z10, Z13, Z16, Z22), rgfa-zip zips nothing"

# ---- Z33: old option strings are errors naming the equivalent
o33() { rgfa-zip -m $MM2 "$@" -o $T/x.gfa -r $T/x.tsv $T/Z1/in.gfa $T/Z1/snarls.json 2>&1 > /dev/null; echo "exit $?"; }
is "$(o33 -D | tail -1)/$(o33 -D | grep -c 'strand is never a gate')" "exit 2/1" "Z33: -D is an error: strand is never a gate"
is "$(o33 -A 2 | grep -c -- '--alt-rounds')/$(o33 -c 5 | grep -c -- '--max-site-nodes')/$(o33 -M 2 | grep -c 'collinear chain')/$(o33 --threads 4 | grep -c 'sites in parallel')/$(o33 -T | grep -c 'rgfa-collapse -T')" "1/1/1/1/1" "Z33: -A, -c, -M, --threads and -T name their rgfa-zip equivalents"

# ---- Z34-Z37: rule U and the gates
for c in Z34a Z34b Z35 Z35b Z36 Z37; do zip $c run; done
is "$(q row Z34a run HG1#1#c1 outcome) | $(q row Z34b run HG1#1#c1 outcome)/$(q same Z34a run)/$(q same Z34b run)" "ambiguous-strand | ambiguous-strand/same/same" "Z34: flip5 and flip7: '+' onto U and '-' onto rc(U') within 2 bases: ambiguous-strand"
is "$(q row Z35 run HG1#1#c1 outcome label kept_bp)/$(q row Z35 run HG1#1#c1 dropped | cut -d: -f1)/$(q alt Z35 run)" "zipped INV 20000/island-below-b/s8:HG1#1#c1:20000:19000" "Z35: hitch: only the 20 kb inversion zips; the 4 kb paralog island stays alt"
is "$(q walks Z35 run)/$(q row Z35b run HG1#1#c1 outcome)" "ok/below-b" "Z35: walks spell; control (the paralog alone) is not confident"
is "$(q row Z36 run HG1#1#c1 outcome kept_bp tie)/$(q row Z36 run HG1#1#c1 blocks | cut -d: -f1-3)/$(q walks Z36 run)" "zipped 20000 tie-kept/0-20000:202000-222000:+/ok" "Z36: two copies 150 kb apart, the alt 0.24% closer to the right one: zipped onto the right copy"
is "$(q row Z37 run HG1#1#c1 outcome internal)/$(q same Z37 run)" "fragmented 13/same" "Z37: a VNTR expansion aligned out of phase (13 stretches in 14.5 kb) is fragmented; no cut"

# ---- Z38-Z39 (injected chains)
for c in Z38 Z39 Z39b Z39c; do zip $c run $(inj $c); done
is "$(q col Z38 run outcome)/$(q pieces Z38 run 1 new)/$(q pieces Z38 run 1 compatible)" "zipped zipped/s6:2100-12000/s5:12033-22000" "Z38: two chains share s5, projected 33 bp apart (compatible): the neighbour is re-anchored to the accepted image, both kept"
is "$(q back Z38 run)/$(q fscc Z38 run)/$(q walks Z38 run)" "0/0/ok" "Z38: no back edge, no forward cycle; walks spell"
is "$(q pieces Z39 run 0)/$(q walks Z39 run)" "s5:2100-2103 s6:2103-8000/ok" "Z39: a 3 bp first piece with a rank-0 boundary 30 bp inside: no inward snap"
is "$(q pieces Z39b run 0)/$(q pieces Z39c run 0)/$(q cuts Z39b run s2)$(q cuts Z39b run s3)$(q cuts Z39c run s2)$(q cuts Z39c run s3)" "s5:2000-8000/s5:2130-8000/" "Z39: block starts 3 bp from the window start and from a node boundary snap to them; no cut"

# ---- Z42-Z43
zip Z42 run
is "$(q row Z42 run HG1#1#c1 outcome tie)/$(q row Z42 run HG1#1#c1 blocks | cut -d: -f1-3)/$(q walks Z42 run)" "zipped unique/0-10000:28000-38000:+/ok" "Z42: identical copies, one infeasible: zipped onto the feasible one, unique"
t0=$(date +%s%N)
zip Z43 run
t1=$(date +%s%N)
is "$(q col Z43 run outcome ref | tr ' ' '\n' | sort | uniq -c | tr -s ' ')/$(q col Z43 run outcome alt | tr ' ' '\n' | sort | uniq -c | tr -s ' ')/$(q col Z43 run round alt | tr ' ' '\n' | uniq -c | tr -s ' ' | tr '\n' ';')/$(grep -c ' 0 to pass 2,' $T/Z43/run.err)" " 20 prefiltered/ 54 prefiltered/ 19 1; 18 2; 17 3;/1" "Z43: a satellite window with 20 queries of 100 kb: all prefiltered by the screen, none reaches pass 2; so are the alt-vs-alt members (19, 18 and 17 over three rounds)"
ok $(( (t1 - t0) / 1000000 < 10000 )) "Z43: under 10 s ($(( (t1 - t0) / 1000000 )) ms)"

# ---- Z24, Z40, Z41: alt-vs-alt (v2, the default)
zip Z24 run
is "$(q col Z24 run outcome alt)/$(q col Z24 run round alt)/$(q col Z24 run label alt)/$(q col Z24 run kept_bp alt)/$(q col Z24 run anchors alt)/$(q col Z24 run representative alt)" "zipped/1/FWD/6000/s1+>s2+/HG1#1#c1:s3 s3+" "Z24: parallel insertions with SR 1 and SR 3 at 2% divergence: SR 3 merged onto SR 1 in round 1"
is "$(q S Z24 run)/$(q L Z24 run)/$(q newL Z24 run)/$(q walks Z24 run)" "s1 s2 s3/s1+s2+ s1+s3+ s3+s2+//ok" "Z24: the member is deleted, no new link; its walk spells the representative"
zip Z24 v1 --no-alt
is "$(q same Z24 v1)/$(q col Z24 v1 pass)" "same/ref ref" "Z24: --no-alt (v1): no alt-vs-alt row, output = input"
zip Z24s run
is "$(q col Z24s run outcome alt)/$(q same Z24s run)" "in-series/same" "Z24: copies in series (a walk r1 R M r2): in-series, not merged"
zip Z24p run
is "$(q col Z24p run outcome alt)/$(q same Z24p run)" "not-private/same" "Z24: a member that also links to another reference node: not-private, not merged"
zip Z24m run
is "$(q col Z24m run outcome alt)/$(q col Z24m run round alt)/$(q S Z24m run)/$(q walks Z24m run)" "zipped zipped/1 1/s1 s2 s3/ok" "Z24: both members of one group merge onto the representative"
zip Z24r run
is "$(q col Z24r run outcome alt)/$(q col Z24r run round alt)/$(q col Z24r run representative alt | cut -d' ' -f5-)/$(q S Z24r run)/$(q walks Z24r run)" "prefiltered prefiltered zipped/1 1 2/HG2#1#c2:s4 s4+/s1 s2 s3 s4/ok" "Z24: rounds: no chain onto the unrelated SR 1 branch; round 2 regroups SR 3 under SR 2 and merges it"
zip Z24r r1 --alt-rounds 1
is "$(q col Z24r r1 outcome alt)/$(q same Z24r r1)" "prefiltered prefiltered/same" "Z24: --alt-rounds 1: no second round"
zip Z24t run
is "$(q col Z24t run outcome alt)/$(q col Z24t run round alt)/$(q col Z24t run representative alt)/$(q S Z24t run)/$(q walks Z24t run)" "zipped/1/HG2#1#c2:s4 s4+/s1 s2 s3 s4/ok" "Z24: a lowest-SR branch shorter than ceil(b*i) is no representative: SR 3 merges onto SR 2 in round 1"
zip Z24b run
is "$(q col Z24b run outcome alt)/$(q col Z24b run label alt)/$(q S Z24b run)/$(q L Z24b run)/$(q newL Z24b run)/$(q walks Z24b run)" "zipped/FWD/s1 s2 s3/s1+s2+ s2-s3+ s3+s1-//ok" "Z24b: both alleles written by reverse-strand contigs: merged, no new link"
zip Z24c run
is "$(q col Z24c run outcome alt)/$(q col Z24c run representative alt)/$(q newL Z24c run)/$(q walks Z24c run)/$(q ident Z24c HG3 0.99)" "zipped/HG1#1#c1:s3 s3-//ok/0.9973 >=0.99" "Z24c: reverse-written representative, forward member: walk-side mapping, the member's walk spells at >= 99%, no inversion edge"
zip Z40 run
is "$(q col Z40 run outcome alt)/$(q col Z40 run anchors alt)/$(q col Z40 run label alt)/$(q col Z40 run representative alt)" "zipped/s3+>s5+/FWD/HG1#1#c1:s3 s4+" "Z40: nested bubble (a b d by rank 1, a c d by rank 5): the branch c, bounded by the alt nodes a and d, merges onto b"
is "$(q S Z40 run)/$(q newL Z40 run)/$(q walks Z40 run)" "s1 s2 s3 s4 s5//ok" "Z40: c pinched onto b: c deleted, no new link; walks spell"
zip Z40b run
is "$(q col Z40b run outcome alt)/$(q col Z40b run label alt)/$(q newL Z40b run)/$(q fscc Z40b run)/$(q walks Z40b run)" "zipped/INV/s3+s4- s4-s5+/0/ok" "Z40: c ~ rc(b): pinched inverted (two inversion links), 0 forward SCCs; walks spell"
zip Z40b t3 -t 3 -j 3
zip Z40b_shuf t2 -t 2 -j 2
is "$(cmp -s $T/Z40b/run.gfa $T/Z40b/t3.gfa && cmp -s $T/Z40b/run.tsv $T/Z40b/t3.tsv && cmp -s $T/Z40b/run.gfa $T/Z40b_shuf/t2.gfa && cmp -s $T/Z40b/run.tsv $T/Z40b_shuf/t2.tsv && echo same)" "same" "Z40: -t 3 -j 3 and shuffled S/L/snarl lines at -t 2 -j 2: byte-identical graph and report"
zip Z40p run
is "$(q col Z40p run outcome alt)/$(q same Z40p run)" "not-private/same" "Z40: with an extra link from c to a reference node: not-private, not merged"
zip Z41 run
is "$(q col Z41 run outcome diag)/$(q col Z41 run label diag)/$(q col Z41 run anchors diag)/$(q col Z41 run window diag)/$(q col Z41 run pass | tr ' ' '\n' | grep -c alt)/$(q same Z41 run)" "near-parallel:disjoint:private/FWD/s2+>s4+/2090-2114/0/same" "Z41: HP shape (insertion at p, member window [p+90, p+114)): not merged; reported near-parallel, disjoint"

# ---- regressions from code review: graph correctness
for c in Z5o Z8o Z5f; do zip $c inj $(inj $c); done
is "$(q row Z5o inj HG1#1#c1 outcome label kept_bp dropped)/$(q pieces Z5o inj 0)/$(q alt Z5o inj)" "zipped +-+ 59995 ./s4:2000-22005 s4:22005-36995 s4:37000-62000/s9:HG1#1#c1:34995:5" "GC-1: A's record runs 5 bp into rc(B): the '-' part is trimmed at its query end (its target's bottom), not lost"
is "$(q walks Z5o inj)/$(q fscc Z5o inj)/$(q back Z5o inj)" "ok/0/0" "GC-1: walks spell; no forward cycle, no back edge"
is "$(q row Z8o inj HG1#1#c1 outcome label kept_bp dropped)/$(q pieces Z8o inj 0)/$(q alt Z8o inj)/$(q walks Z8o inj)" "zipped -+- 59997 ./s4:36997-62000 s4:22003-36997 s4:2000-22000/s9:HG1#1#c1:39997:3/ok" "GC-1: frame R: rc(C) runs 3 bp into B; B's '+' part is trimmed at its query end (its target's top)"
is "$(q row Z5f inj HG1#1#c1 outcome kept_bp dropped)/$(q walks Z5f inj)" "zipped 45000 s4[20000-35000):(fold)/ok" "GC-1: a part whose target lies inside an earlier part's is reported dropped as a fold; its query stays alt"
zip Z5r run
is "$(q row Z5r run HG1#1#c1 outcome label kept_bp)/$(q alt Z5r run)/$(q walks Z5r run)" "zipped +- 34988/s7:HG1#1#c1:34988:12/ok" "GC-1: real minimap2, A rc(B) with a 12 bp breakpoint microhomology: the inversion is zipped (it used to stay alt), 12 bp stay alt"
zip Z23z run
is "$(ex Z23z run)/$(q col Z23z run outcome site)/$(q same Z23z run)" "0/site:leak site:leak/same" "GC-2: Z23 with a zippable allele beside the leaking site: both leaking sites are skipped (it was exit 3)"
for c in Z1t Z1h; do zip $c run; done
is "$(ex Z1t run)/$(q row Z1t run HG1#1#c1 outcome label kept_bp)/$(q extra Z1t run)/$(q walks Z1t run)" "0/zipped INV 6000/s1+s2- s2-s3+ s5+s3+/ok" "GC-2: a tip on B.L only (inside vg's snarl) is interior: the site zips (it was exit 3)"
is "$(ex Z1h run)/$(q row Z1h run HG1#1#c1 outcome label kept_bp)/$(q extra Z1h run)/$(q walks Z1h run)" "0/zipped INV 6000/s1+s2- s2-s3+ s3-s5+ s5+s3+/ok" "RD-2: a fold-back on B.L (s3- z+ s3+) is interior: the site zips (it was exit 3)"
zip Z1d run
is "$(ex Z1d run)/$(q col Z1d run source)/$(cmp -s $T/Z1d/run.gfa $T/Z1/expected.gfa && echo same)/$(grep -c 'duplicate link' $T/Z1d/run.err)" "0/creator/same/1" "GC-3: a creator link also written reversed without SR: the copy with SR is kept, the run keeps its creator excursion; same edit as Z1"

# ---- regressions from code review: robustness, resources, report and command line
zip Z21 cap --stage detect --max-site-query 10000
is "$(q row Z21 cap HG1#1#c1 outcome query_bp feasible_bp)" "site:capped 12000 6000" "RD-5: --max-site-query caps the query of the units aligned (12 kb here), not their feasible bp (6 kb)"
rgfa-zip -m $MM2 --tmpdir $T/tmp -o $T/Z35/nodump.gfa -r $T/Z35/nodump.tsv $T/Z35/in.gfa $T/Z35/snarls.json 2> /dev/null
is "$(cmp -s $T/Z35/nodump.tsv $T/Z35/run.tsv && cmp -s $T/Z35/nodump.gfa $T/Z35/run.gfa && echo same)/$(q col Z35 run records)" "same/2" "RD-6: without --dump the records are dropped after chaining: graph and report (records column) as with --dump"
o9() { rgfa-zip -m $MM2 --stage detect "$@" -r $T/x9.tsv $T/Z1/in.gfa $T/Z1/snarls.json 2>&1 > /dev/null | grep -c 'lowering -j'; }
is "$(o9 -t 2 -j 2 --mem 3000000000)/$(o9 -t 2 -j 2 --mem 2999999999)" "0/1" "RD-9: --mem counts 1.5 GB per aligner as 1.5e9 bytes, as cactus does"
zip Z1 detonly --detect-only
is "$(ex Z1 detonly)/$(q row Z1 detonly HG1#1#c1 outcome)/$(q same Z1 detonly)/$(grep -c 'would zip 1 chain' $T/Z1/detonly.err)" "0/would-zip/same/1" "RD-12: --detect-only reports would-zip, not zipped, and says so in the log; the graph is unchanged"
zip Z13 big --stage detect --max-site-nodes 2
is "$(q col Z13 big outcome site)/$(q col Z13 big query_bp site)" "site:too-big/36017" "RD-13: a site:too-big row carries the site's whole alt bp, not what the flood saw before the cap"
python3 -c "import signal, os, sys; signal.signal(signal.SIGCHLD, signal.SIG_IGN); os.execvp(sys.argv[1], sys.argv[1:])" rgfa-zip -m $MM2 --tmpdir $T/tmp -o $T/Z1/sigchld.gfa -r $T/Z1/sigchld.tsv $T/Z1/in.gfa $T/Z1/snarls.json 2> $T/Z1/sigchld.err
echo $? > $T/Z1/sigchld.exit
is "$(ex Z1 sigchld)/$(cmp -s $T/Z1/sigchld.gfa $T/Z1/run.gfa && echo same)" "0/same" "RD-14: launched with SIGCHLD ignored: the children's exit status is still read (it was exit 2)"
mkdir -p $T/Z28_trunc $T/Z1_nosnarls
cp $T/Z28/in.gfa $T/Z28_trunc/ && head -12 $T/Z28/snarls.json > $T/Z28_trunc/snarls.json
cp $T/Z1/in.gfa $T/Z1_nosnarls/ && : > $T/Z1_nosnarls/snarls.json
zip Z28_trunc run
zip Z1_nosnarls run
is "$(ex Z28_trunc run)/$(none Z28_trunc run)/$(grep -c 'lie in no snarl' $T/Z28_trunc/run.err)/$(ex Z1_nosnarls run)/$(none Z1_nosnarls run)" "2/nothing/1/2/nothing" "RD-15: a truncated (12 of 25 snarls) or empty snarls file leaves alt nodes in no site: exit 2"
zip Z1 hdr --stage detect -i 0.95004 --delta 0.00004 --prefilter 2000,0.2004
is "$(grep '^#options' $T/Z1/hdr.tsv | grep -o -- '-i [^ ]*\|--delta [^ ]*\|--prefilter [^ ]*' | tr '\n' ' ')" "-i 0.95004 --delta 4e-05 --prefilter 2000,0.2004 " "RD-16: the report header prints option values exactly"
a10() { rgfa-zip -m $MM2 "$@" $T/Z1/in.gfa $T/Z1/snarls.json > /dev/null 2>&1; echo $?; }
mkdir -p $T/Z1/al_dump $T/Z1/al_dir
is "$(a10 -o $T/Z1/al.tsv.tmp -r $T/Z1/al.tsv)/$(a10 -o $T/Z1/al_dump/units.tsv -r $T/Z1/al2.tsv --dump $T/Z1/al_dump)/$(a10 -o $T/Z1/al3.gfa -r $T/Z1/al_dir)/$(a10 -o $T/Z1/in.gfa -r $T/Z1/al4.tsv)/$(ls $T/Z1/al.tsv* $T/Z1/al2.tsv* $T/Z1/al3.gfa* $T/Z1/al4.tsv* 2>/dev/null | wc -l)$(ls -A $T/Z1/al_dump | wc -l)" "2/2/2/2/00" "RD-10: outputs that alias one another (-o is -r's .tmp; -o is a --dump file), a directory, or an input: exit 2 before anything runs"

# ---- regressions from the genome-wide sweep
# paf-inconsistent: ambiguous bases in I and D runs are not in PAF column 11 (minimap2 2.30)
zip Zn1 run
is "$(ex Zn1 run)/$(grep -c 'paf-inconsistent:' $T/Zn1/run.err)/$(q row Zn1 run HG1#1#c1 outcome label records kept_bp blocks)" "0/0/zipped FWD 1 19864 0-10000:2000-12000:+:0.9980;10000-19864:12106-22000:+:0.9970" "paf-inconsistent: a 106 bp D run over 10 window Ns (and a 30 bp one over 5) is consistent; the 20 N-facing-N columns are X, not '=' (0.9980); the short D counts in the identity (0.9970)"
is "$(q cuts Zn1 run s2)/$(q walks Zn1 run)" "10000 10106/ok" "paf-inconsistent: the long D is a deletion edge over the window's Ns, the short one is absorbed; walks spell"
zip Zn2 run
is "$(ex Zn2 run)/$(grep -c 'paf-inconsistent:' $T/Zn2/run.err)/$(q row Zn2 run HG1#1#c1 outcome label kept_bp blocks)" "0/0/zipped INV 20020 0-10020:12000-22000:-:0.9980;10090-20090:2000-12000:-:1.0000" "paf-inconsistent: '-' strand, a 70 bp I run holding 10 query Ns (and a 20 bp one holding 4, inside a block: 0.9980) is consistent"
is "$(q alt Zn2 run)/$(q walks Zn2 run)" "s7:HG1#1#c1:10020:70/ok" "paf-inconsistent: the 70 bp insertion stays alt at its projected point; walks spell"
# kept_bp: a node interval that several compatible chains zip counts once (Z13: the second chain's pieces are compatible)
is "$(q col Z13 run kept_bp | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i != ".") s += $i; print s}')/$(q removed Z13 run)/$(grep -o 'accepted [0-9]* chain(s), [0-9]* bp zipped; 0 reverted; [0-9]* of the chain(s) zip bp of their own, [0-9]* only agree' $T/Z13/run.err)" "36010/36010/accepted 2 chain(s), 36010 bp zipped; 0 reverted; 1 of the chain(s) zip bp of their own, 1 only agree" "kept_bp: the report's column sums to the bp zipped and to the alt bp the graph lost; the log counts the compatible chain apart"
# rule (G), the whole-record rule (v3): a GAF record through the node also reads the target
zip Zg loop --walks gaf:$T/Zg/loop.gaf --audit-gaf $T/Zg/loop.gaf
is "$(ex Zg loop)/$(q row Zg loop HG1#1#c1 outcome kept_bp dropped)/$(q same Zg loop)/$(grep -o 'rule (G) .*' $T/Zg/loop.err | cut -d: -f2)/$(grep -o 'GAF audit: .*' $T/Zg/loop.err | grep -o 'PASS\|FAIL')" "0/trimmed-below-b 0 s6[0-6000):(G)/same/ 1 piece(s), 6000 bp, refused in 1 chain(s)/PASS" "(G): HG1's record walks a in an excursion over r2, then folds back at r4 and reads r2: a's piece is refused, output = input, the GAF audit passes"
zip Zg straight --walks gaf:$T/Zg/straight.gaf --audit-gaf $T/Zg/straight.gaf
zip Zg creator --audit-gaf $T/Zg/loop.gaf
is "$(q row Zg straight HG1#1#c1 outcome label kept_bp)/$(q walks Zg straight)/$(grep -o 'GAF audit: .*' $T/Zg/straight.err | grep -o 'PASS\|FAIL') | $(q row Zg creator HG1#1#c1 outcome kept_bp)/$(grep -o '[0-9]* GAF lines in [0-9]* chains' $T/Zg/creator.err)/$(grep -o 'GAF audit: .*' $T/Zg/creator.err | grep -o 'PASS\|FAIL')" "zipped INV 6000/ok/PASS | zipped 6000/1 GAF lines in 1 chains/FAIL" "(G): the same record without the fold-back zips a; creator walks (no rule (G)) zip it too, and the audit's reading rule fails on exactly the record (G) refuses"
# logs: waiting for the -j budget apart from running, and minimap2's own times (three sites, -t 3 -j 1:
# every alignment holds the only slot 0.5 s longer, and reports real 1.25 s, CPU 0.75 s)
FAKE_MODE=timing zip Z31 timing -m $T/fake_mm2.sh -t 3 -j 1 -v
is "$(ex Z31 timing)/$(cmp -s $T/Z31/t1.gfa $T/Z31/timing.gfa && cmp -s $T/Z31/t1.tsv $T/Z31/timing.tsv && echo same)/$(grep -o "minimap2's own figures, summed over its processes: real [0-9.]* s, CPU [0-9.]* s" $T/Z31/timing.err)" "0/same/minimap2's own figures, summed over its processes: real 7.5 s, CPU 4.5 s" "logs: the run's minimap2 real and CPU time are minimap2's own figures (6 processes x 1.25 s / 0.75 s); output as at -t 1 -j 1"
python3 - $T/Z31/timing.err > $T/Z31/timing.check <<'PYEOF'
import re, sys
win = re.compile(r'aligner: window \S+: (\d+) queries, ([0-9.]+) s \(minimap2 running ([0-9.]+) s, waiting for an aligner slot ([0-9.]+) s\); '
                 r'(\d+) minimap2 process\(es\): real ([0-9.]+) s, CPU ([0-9.]+) s')
site = re.compile(r'site \S+: ([0-9.]+) s: minimap2 running ([0-9.]+) s, waiting for an aligner slot ([0-9.]+) s, rgfa-zip alone ([0-9.]+) s; '
                  r'(\d+) minimap2 process\(es\): real ([0-9.]+) s, CPU ([0-9.]+) s')
W = [tuple(float(x) for x in m.groups()) for m in map(win.search, open(sys.argv[1])) if m]
S = [tuple(float(x) for x in m.groups()) for m in map(site.search, open(sys.argv[1])) if m]
own = all(w[4] == 2 and w[5] == 2.5 and w[6] == 1.5 for w in W) and all(s[4] == 2 and s[5] == 2.5 and s[6] == 1.5 for s in S)
split = all(w[2] >= 1.0 and w[1] - w[3] - w[2] < 1.0 for w in W) and all(abs(s[0] - s[1] - s[2] - s[3]) < 0.25 for s in S)
print('windows %d sites %d own %s waited %s split %s' % (len(W), len(S), own, max([w[3] for w in W] + [0]) >= 0.4, split))
PYEOF
is "$(cat $T/Z31/timing.check)" "windows 3 sites 3 own True waited True split True" "logs (-v): every window and site line gives minimap2's own times (2 processes: 2.5 s / 1.5 s); a window that waited for the only slot says so, and its time outside the waits is its two alignments, not the queue"
# Z24 aligns only in its alt-vs-alt pass (two alignments of >= 0.5 s): that pass's line gives its minimap2
# time, and edit planning is the planner's own time (it used to include the alt pass and its waits)
FAKE_MODE=timing zip Z24 timing -m $T/fake_mm2.sh -v
is "$(ex Z24 timing)/$(cmp -s $T/Z24/run.gfa $T/Z24/timing.gfa && echo same)/$(grep -o 'alt-vs-alt [0-9.]* s, of which [0-9.]* s waiting for an aligner slot .*minimap2 running [0-9.]* s, CPU [0-9.]* s' $T/Z24/timing.err | grep -o 'CPU [0-9.]* s')/$(grep -o 'edit planning [0-9.]* s' $T/Z24/timing.err | awk '{print ($3 < 0.5) ? "planner only" : "includes the alt pass"}')" "0/same/CPU 1.5 s/planner only" "logs: the alt-vs-alt pass reports its minimap2 time; edit planning no longer counts the alt pass's alignments"
# the 500 sub-record cap is in the report row (it was a stderr counter only)
zip Ztr run $(inj Ztr) -b 100 --min-piece 100
is "$(ex Ztr run)/$(q row Ztr run HG1#1#c1 dropped)/$(grep -c ' 1 truncated at 500 sub-records' $T/Ztr/run.err)" "0/truncated:500/520-sub-records/1" "truncation: a unit with 520 feasible sub-records says in its row that the chain was chosen among the best 500"

# keep the outputs of a failed run (or with KEEP set) for inspection
failed=0
for i in $(seq 1 $_bt_current_test); do [ "${_bt_test_ok[$i]}" = 1 ] || failed=1; done
if [ -z "$KEEP" ] && [ $failed = 0 ] && [ "$_bt_current_test" = "$(expected_tests)" ]; then rm -rf $T; fi
