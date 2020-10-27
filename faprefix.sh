# add a prefix to a fasta file
#
# example to make consistent paf and fasta where all graph sequences (ie contig names) get the "mg_ancors" prefix
# 
# gfatools gfa2fa graph.gfa | faprefix.sh mg_anchors > graph.gfa.fa
# mzgaf2paf aln.gaf -p mg_anchors > aln.paf
#
#!/bin/bash
cat - | sed -e "s/^>\(.*\)/>${1}\1/g"
