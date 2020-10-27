#!/usr/bin/env python3

import os
import sys
import argparse
import string
import re

# fasta reading and reversing from https://github.com/ComparativeGenomicsToolkit/sonLib/blob/master/src/sonLib/bioio.py
def _getFileHandle(fileHandleOrFile, mode="r"):
    if isinstance(fileHandleOrFile, "".__class__):
        return open(fileHandleOrFile, mode)
    else:
        return fileHandleOrFile

def fastaRead(fileHandleOrFile):
    """iteratively yields a sequence for each '>' it encounters, ignores '#' lines
    """
    fileHandle = _getFileHandle(fileHandleOrFile)
    line = fileHandle.readline()
    valid_chars = {x for x in string.ascii_letters + "-"}
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = []
            while line != '' and line[0] != '>':
                line = re.sub("[\\s]+", "", line)
                if len(line) > 0 and line[0] != '#':
                    seq.append(line)
                line = fileHandle.readline()
            seq = "".join(seq)
            try:
                assert all(x in valid_chars for x in seq)
            except AssertionError:
                bad_chars = {x for x in seq if x not in valid_chars}
                raise RuntimeError("Invalid FASTA character(s) see in fasta sequence: {}".format(bad_chars))
            yield name, seq
        else:
            line = fileHandle.readline()
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()

def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    dNA = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c' }
    def fn(i):
        if i in dNA:
            return dNA[i]
        return i
    return "".join([ fn(i) for i in seq ])

        
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
    args = args[1:]

    return parser.parse_args(args)

# make sure that our cigar matches are exact matches (as we'd expect from the minimizer input)
def check_cigar(paf_line, fa_dict):
    toks = paf_line.rstrip().split("\t")
    cigar = toks[-1]
    print (toks)
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
        target_seq = reverseComplement(target_seq)

    query_pos = int(toks[2])
    target_pos = int(toks[7])

    cigar_toks = re.findall('([0-9]+)(M|D|I)', cigar[4:])

    for cig_len, cig_type in cigar_toks:
        if cig_type == "M":
            query_e = query_pos + int(cig_len)
            query_frag = query_seq[query_pos:query_e]
            target_e = target_pos + int(cig_len)
            target_frag = target_seq[target_pos:target_e]
            if query_frag != target_frag:
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
    
def main(args):

    options = parse_args(args)

    # load and index the fasta sequences by name
    fa_dict = {}
    def load_fa(in_path):
        for seqname, seq in fastaRead(in_path):
            assert seqname not in fa_dict
            fa_dict[seqname.split(" ")[0]] = seq
    load_fa(options.fasta1)
    load_fa(options.fasta2)

    with open(options.paf, "r") as paf_file:
        for paf_line in paf_file:

            check_cigar(paf_line, fa_dict)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
    
