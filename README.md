# cactus-gfa-tools

Command-line utilitites required for the [Cactus Pangenome Pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pangenome pipeline.

## Build

```
git clone https://github.com/ComparativeGenomicsToolkit/cactus-gfa-tools.git
cd cactus-gfa-tools && make
```

### Run tests

Requires:
* [minigraph](https://github.com/lh3/minigraph)
* [gfatools](https://github.com/lh3/gfatools)
* [minimap2](https://github.com/lh3/minimap2)

```
export PATH=$(pwd):$PATH
make test
```

## Tools Included

* [paf2lastz](#paf2lastz)
* [mzgaf2paf](#mzgaf2paf)
* [rgfa-split](#rgfa-split)
* pafcoverage (for debugging only)

### paf2lastz

Convert PAF with cg cigars to LASTZ cigars

This is a re-implentation of the following [paftools](https://github.com/lh3/minimap2/blob/master/misc/paftools.js) command:
```
paftools view -f lastz-cigar
```
that
* is standalone (ie does not require Javascript)
* incorporates @Robin-Rounthwaite's [fix](https://github.com/Robin-Rounthwaite/reference-based-cactus-aligner/blob/master/src/paf_to_lastz.py#L49-L71)
* provides option (`-q`) to use the PAF MAPQ as the score

usage:
```
paf2lastz a.paf > a.cigar
```
where `a.paf` was created by, for example, `minimap2 -c`

### mzgaf2paf

Convert [minigraph](https://github.com/lh3/minigraph) output from GAF to PAF, where PAF records represent pairwise aligments between the target and nodes in the graph.  The output alignments have cigar strings, and are based on the minimizer offsets obtained when using `minigraph -S --write-mz`.

The use case is that a set of such alignments can be used as anchors by [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) to form its initial graph. 


#### Example

```
cd test
### Build a graph
minigraph -xggs -l10k hpp-20-2M/CHM13.fa.gz hpp-20-2M/HG003.fa.gz hpp-20-2M/HG004.fa.gz > hpp-20-2M.gfa
### Align a contig
minigraph -xasm -t $(nproc) -K4g --inv=no -S --write-mz hpp-20-2M.gfa hpp-20-2M/hg38.fa.gz > hg38.gaf
### Convert to PAF
mzgaf2paf hg38.gaf > hg38.paf
```

This GAF from the input
```
hg38.chr20  2833756  60006 2352627  +  <CHM13.CHM13_Super-Scaffold_117:2060700-2381813>HG003.HG003_h1tg000030l:316552-316571<CHM13.CHM13_Super-Scaffold_117:2060427-2060561<CHM13.CHM13_Super-Scaffold_117:2060086-2060231>HG003.HG003_h1tg000030l:316802-316827<CHM13.CHM13_Super-Scaffold_117:1930843-2059970>HG003.HG003_h1tg000030l:445928-445984<CHM13.CHM13_Super-Scaffold_117:1518313-1930656>HG004.HG004_h2tg000013l:858527-858576<CHM13.CHM13_Super-Scaffold_117:1517794-1517951>HG003.HG003_h1tg000030l:858750-858779<CHM13.CHM13_Super-Scaffold_117:1517556-1517656>HG003.HG003_h1tg000030l:858852-858997<CHM13.CHM13_Super-Scaffold_117:1385660-1516344>HG003.HG003_h1tg000030l:989703-989806<CHM13.CHM13_Super-Scaffold_117:1075838-1385660<CHM13.CHM13_Super-Scaffold_117:1031686-1075523<CHM13.CHM13_Super-Scaffold_117:588846-1030548<CHM13.CHM13_Super-Scaffold_117:488710-588690>HG004.HG004_h2tg000013l:1888505-1892245<CHM13.CHM13_Super-Scaffold_117:0-475698  2369008  100466   2368992  2228396  2303074  60 tp:A:P   cm:i:395124 s1:i:2187209   s2:i:1092   dv:f:0.0021
*  <s43  97 12 0  6  92 709771   709857   19 4,2,10,7,2,9,10,8,4,6,5 4,2,10,7,2,9,10,8,4,6,5                                                                                                                                
```
Will be converted to the following PAF
```
hg38.chr20  2833756  709771   709857   -  s43   97 5  91 86 86 255   cg:Z:86M                       
```

A variety of filters are available.  Use `mzgaf2paf -h` to list them.

### rgfa-split

Use [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) to assign query contigs from a PAF to reference contigs from the rGFA.  It performs the following steps.

1. Using a BFS from the rank-0 nodes, assign each other node to the stable sequence (reference contig) of its nearest rank-0 node
2. Use the PAF to determine the alignment coverage for each query contig against each reference contig
3. Assign each query contig to the reference contig to which it most aligns

Filters are provided to only assign contigs if they pass specificity and uniqueness filters (and label as "ambiguous" if they don't).  See `rgfa-split -h` for a full list of options. 
