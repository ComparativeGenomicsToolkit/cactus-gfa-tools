#!/usr/bin/env python3

import os
import sys
import argparse
import string
import re
from Bio import SeqIO

        
def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("paf", type=str,
                        help="paf whose cigar strings we want to check")
    parser.add_argument("fasta1", type=str,
                        help="fasta1")
    parser.add_argument("fasta2", type=str,
                        help="fasta2")
    parser.add_argument("--gaf", action="store_true",
                        help="expect gaf instead of paf")
    args = args[1:]

    return parser.parse_args(args)

# make sure that our cigar matches are exact matches (as we'd expect from the minimizer input)
def check_cigar(paf_line, fa_dict):
    toks = paf_line.rstrip().split("\t")
    cigar = toks[-1]
    assert cigar[:4] == "cg:Z"

    query_name = toks[0]
    assert query_name in fa_dict
    query_seq = fa_dict[query_name]
    assert len(query_seq) == int(toks[1])
    
    target_name = toks[5]
    assert target_name in fa_dict
    target_seq = fa_dict[target_name]
    assert len(target_seq) == int(toks[6])

    assert toks[4] in ('-', '+')
    if toks[4] == '-':
        target_seq = target_seq.reverse_complement()

    query_pos = int(toks[2])
    target_pos = int(toks[7])

    cigar_toks = re.findall('([0-9]+)(M|D|I)', cigar[4:])

    for cig_len, cig_type in cigar_toks:
        if cig_type == "M":
            query_e = query_pos + int(cig_len)
            query_frag = query_seq[query_pos:query_e]
            target_e = target_pos + int(cig_len)
            target_frag = target_seq[target_pos:target_e]
            if query_frag.upper() != target_frag.upper():
                sys.stderr.write("Validation Error\n\t{}\n".format(paf_line))
                sys.stderr.write("\tCigar : {}{} :\n\tquery[{}:{}] = \"{}\"\n\ttarget[{}:{}] = \"{}\"\n".format(
                    cig_len, cig_type, query_pos, query_e, query_frag, target_pos, target_e, target_frag)) 
            assert query_frag.upper() == target_frag.upper()
        if cig_type != "I":
            target_pos += int(cig_len)
        if cig_type != "D":
            query_pos += int(cig_len)

    assert query_pos == int(toks[3])
    assert target_pos == int(toks[8])

def check_mz_offsets(gaf_toks, query_name, fa_dict):

    # parse the line
    target_name = gaf_toks[1]
    assert target_name[0] in ('<', '>')
    target_reversed = target_name[0] == '<'
    target_name = target_name[1:]
    num_minis = int(gaf_toks[3])

    if num_minis == 0:
        return

    target_start = int(gaf_toks[5])
    target_end = int(gaf_toks[6])
    query_start = int(gaf_toks[7])
    query_end = int(gaf_toks[8])
    kmer_size = int(gaf_toks[9])
    target_offsets = [int(x) for x in gaf_toks[10].split(",")]
    query_offsets = [int(x) for x in gaf_toks[11].split(",")]
    assert len(target_offsets) == len(query_offsets)
    
    assert query_name in fa_dict
    query_seq = fa_dict[query_name]
    
    assert target_name in fa_dict
    target_seq = fa_dict[target_name]
    if target_reversed:
        target_seq = target_seq.reverse_complement()

    query_pos = query_start
    target_pos = target_start
    for i in range(num_minis):
        query_frag = query_seq[query_pos:query_pos + kmer_size]
        target_frag = target_seq[target_pos:target_pos + kmer_size]
        is_match = query_frag.upper() == target_frag.upper()
        sys.stderr.write("[{}] qpos={} tpos={} {} {} {} {}\n".format(i, query_pos, target_pos, query_frag, "==" if is_match else "!=", target_frag, " *** Mismatch *** " if not is_match else ""))
        assert is_match
        if i < num_minis - 1:
            query_pos += query_offsets[i]
            target_pos += target_offsets[i]
    
def main(args):

    options = parse_args(args)

    # load and index the fasta sequences by name
    fa_dict = {}
    def load_fa(fa_path):
        with open(fa_path, 'r') as fa_file:
            for seq_record in SeqIO.parse(fa_file, 'fasta'):
                fa_dict[seq_record.id] = seq_record.seq

    load_fa(options.fasta1)
    load_fa(options.fasta2)

    with open(options.paf, "r") as aln_file:
        if options.gaf:
            # check the GAF output from minigraph to make sure the minimizers match up between query and target
            for line in aln_file:
                toks = line.rstrip().split()
                if toks[0] != '*':
                    query_name = toks[0]
                else:
                    check_mz_offsets(toks, query_name, fa_dict)
        else:
            # check the PAF output from the coverted GAM to make sure the matches in the cigar strings are exact
            for line in aln_file:
                check_cigar(line, fa_dict)

    print("OK!")
                                    

if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
