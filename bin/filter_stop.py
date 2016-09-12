#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Filter out stop codons

Usage:
     filter_stop.py <input.fa >output.fa

'''
from __future__ import print_function

import os
import sys
import argparse
from Bio import SeqIO


def main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser()
    p.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    p.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    a = p.parse_args()

    sample = []
    n = skip = 0
    for r in SeqIO.parse(a.infile, 'fasta'):
        n += 1
        p = r.seq.translate(cds=False, to_stop=True)
        if len(p)*3 != len(r.seq):
            skip += 1
            continue
        sample.append(r)
        
    SeqIO.write(sample, a.outfile, 'fasta')
    print("Skipped {}/{} ({} %)".format(skip, n, 100.0*(float(skip)/n)), file=sys.stderr)



if __name__ == "__main__":
   main(sys.argv[1:])
   


