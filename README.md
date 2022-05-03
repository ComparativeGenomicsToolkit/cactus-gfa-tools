# cactus-gfa-tools

Command-line utilitites required for the [Cactus Pangenome Pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).

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

* [gaf2paf](#gaf2paf)
* [gaf2unstable](#gaf2unstable)
* [paf2lastz](#paf2lastz)
* [mzgaf2paf](#mzgaf2paf)
* [rgfa-split](#rgfa-split)
* [rgfa2paf](#rgfa2paf)
* pafcoverage (for debugging only)
* [pafmask](#pafmask)

### gaf2paf

Converts [minigraph](https://github.com/lh3/minigraph) output from GAF to PAF. Requires `minigraph` be run with `-c` to produce cigars. The output PAF is in stable sequence space.

**Important** The lengths of the target sequences are required to produce a valid PAF. This information is not in the input GAF, therefore the `-l` option must be used to specifiy a table of sequence sizes.  This file can be produced by concatenating the results of `samtools faidx` for each fasta file. Example:

```
minigraph -cxggs 1.fa 2.fa 3.fa > graph.gfa
minigraph -cxasm graph.gfa 1.fa  > 1.gaf
for i in 1 2 3 ; do samtools faidx $i.fa; done
cat 1.fa.fai 2.fa.fai 3.fa.fai > lengths.tsv
gaf2paf 1.gaf -l lengths.tsv > 1.paf
```

### gaf2unstable

By default, [gaf2paf](#gaf2paf) writes the output PAF in stable coordinates (same as the GAF). This means the target sequences are contigs like `chr1` and not minigraph nodes (like `s1`). The Minigraph-Cactus pipeline requires paf in unstable coordinates for best performance (for the time being). `gaf2unstable` can be used to switch over to unstable coordinates as follows. Note that its `-o` option is used to make the lengths table, so using `samtools faidx` as above is no longer required. But instead of the lengths, it needs to read rGFA file from `minigraph`. Example:

```
minigraph -cxggs 1.fa 2.fa 3.fa > graph.gfa
minigraph -cxasm graph.gfa 1.fa  > 1.gaf
gaf2unstable 1.gaf -g graph.gfa -o node-lengths.tsv > 1u.gaf
gaf2paf 1u.gaf -l node-lengths.tsv > 1u.paf

```

### paf2lastz

*(Cactus now works with PAF natively, so this tool is no longer needed)*

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

*(Now that minigraph can produce cigars with `-c`, this tool is no longer needed. To convert GAF-with-cigar output, see [gaf2paf](#gaf2paf) and # [gaf2unstable](#gaf2unstable) above)*

Convert [minigraph](https://github.com/lh3/minigraph) output from GAF to PAF, where PAF records represent pairwise aligments between the query and nodes in the graph.  The output alignments have cigar strings, and are based on the minimizer offsets obtained when using `minigraph -S --write-mz`.

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

### paf2stable

Converts a PAF where the targets are in node-space (ie the output of mzgaf2paf) into a stable PAF in sequence space (comparable to the output of gaf2paf)

### rgfa-split

Use [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) to assign query contigs from a PAF to reference contigs from the rGFA.  It performs the following steps.

1. Using a BFS from the rank-0 nodes, assign each other node to the stable sequence (reference contig) of its nearest rank-0 node
2. Use the PAF to determine the alignment coverage for each query contig against each reference contig
3. Assign each query contig to the reference contig to which it most aligns

Filters are provided to only assign contigs if they pass specificity and uniqueness filters (and label as "ambiguous" if they don't).  See `rgfa-split -h` for a full list of options. 

### rgfa2paf

Generate a PAF from the rank-0 [rGFA tags](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md) in the given rGFA file. Each PAF line will represent an exact alignment between the contig range from the tags to the given node.  Note that the query lengths are inferred from the rGFA and will be wrong unless the file is complete.


### pafmask

Clip out query intervals supplied in a BED file from a given PAF file, updating cg cigars appropriately.  Intervals can by glommed together with `-p` and small fragments filtered out with `-m`.  This can be used, for example, to prevent homologies within centromeric regions from making it into the graph.  In some cases, this works better than masking out the centromeres before generating the paf with minigraph or minimap2: by aligning through the centromeres, the regions nearby end up with better alignments.   
