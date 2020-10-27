# mzgaf2paf

Converts [minigraph](https://github.com/lh3/minigraph) output from GAF to PAF, where PAF records represent pairwise aligments between the target and nodes in the graph.  The output alignments have cigar strings, and are based on the minimizer offsets obtained when using `minigraph -S --write-mz`.

The use case is that a set of such alignments can be used as anchors by [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) to form its initial graph. 

## Build

```
git clone https://github.com/glennhickey/mzgaf2paf.git
cd mzgaf2paf && make
```

## Run tests

Requires:
* [minigraph](https://github.com/lh3/minigraph)
* [gfatools](https://github.com/lh3/gfatools)

```
export PATH=$(pwd):$PATH
cd test
prove -v test.t
```

## Example

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
hg38.chr20  2833756  709771   709857   -  s43   97 6  92 86 86 255   cg:Z:86M
```

