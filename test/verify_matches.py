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
    parser.add_argument("--min-identity", type=float, default=1.0,
                        help="minimum identity for matches (len > 100) for cigar comparison (for validation on minimap2 output only)")
    args = args[1:]

    return parser.parse_args(args)

# compute pct identity
def pct_identity(s1, s2, ignore_n=False):
    assert len(s1) == len(s2)
    same = 0
    for a,b in zip(s1, s2):
        if a == b or (ignore_n and (a == 'N' or b == 'N')):
            same += 1
    return float(same) / float(len(s1))
    
# make sure that our cigar matches are exact matches (as we'd expect from the minimizer input)
def check_cigar(paf_line, fa_dict, min_identity):
    toks = paf_line.rstrip().split("\t")
    cigar = toks[-1]
    assert cigar[:4] == "cg:Z"

    query_start = int(toks[2])
    query_end = int(toks[3])
    target_start = int(toks[7])
    target_end = int(toks[8])
    
    query_name = toks[0]
    if query_name not in fa_dict:
        raise RuntimeError("Query name {} not found in fasta".format(query_name))
    query_seq = fa_dict[query_name][query_start:query_end]
    assert len(query_seq) == query_end-query_start
    assert len(fa_dict[query_name]) == int(toks[1])
    
    target_name = toks[5]
    assert target_name in fa_dict
    target_seq = fa_dict[target_name][target_start:target_end]
    assert len(target_seq) == target_end - target_start
    assert len(fa_dict[target_name]) == int(toks[6])

    assert toks[4] in ('-', '+')
    if toks[4] == '-':
        target_seq = target_seq.reverse_complement()

    query_pos = 0
    target_pos = 0

    cigar_toks = re.findall('([0-9]+)(=|X|M|D|I)', cigar[4:])
    if toks[4] == '-':
        cigar_toks = reversed(cigar_toks)

    for cig_len, cig_type in cigar_toks:
        if cig_type in ["M", "="]:
            query_e = query_pos + int(cig_len)
            query_frag = query_seq[query_pos:query_e]
            target_e = target_pos + int(cig_len)
            target_frag = target_seq[target_pos:target_e]
            iden = pct_identity(query_frag.upper(), target_frag.upper(), ignore_n = min_identity < 1)
            if (min_identity == 1 and iden < 1) or (len(query_frag) > 100 and iden < min_identity):
                sys.stderr.write("Validation Error iden={} < min={}\n\t{}\n".format(iden, min_identity, paf_line))
                sys.stderr.write("\tCigar : {}{} :\n\tquery[{}:{}] = \"{}\"\n\ttarget[{}:{}] = \"{}\"\n".format(
                    cig_len, cig_type, query_pos, query_e, query_frag, target_pos, target_e, target_frag))
                sys.exit(1)
        if cig_type != "I":
            target_pos += int(cig_len)
        if cig_type != "D":
            query_pos += int(cig_len)

    assert query_pos == query_end - query_start
    assert target_pos == target_end - target_start

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

    line_count = 0
    with open(options.paf, "r") as aln_file:
        if options.gaf:
            # check the GAF output from minigraph to make sure the minimizers match up between query and target
            for line in aln_file:
                line_count += 1
                toks = line.rstrip().split()
                if toks[0] != '*':
                    query_name = toks[0]
                else:
                    check_mz_offsets(toks, query_name, fa_dict)
        else:
            # check the PAF output from the coverted GAM to make sure the matches in the cigar strings are exact
            for line in aln_file:
                line_count += 1
                check_cigar(line, fa_dict, options.min_identity)

    if line_count > 0:
        print("OK!")
    else:
        raise RuntimeError("Empty Input")
        

if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
