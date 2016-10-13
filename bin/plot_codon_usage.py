#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Plot the codon usage of various fasta files

This reads an aggregated results.csv file and outputs a single
interactive bar chart (html format, using plotly JS library) that
shows the codon usage patterns in each of the sampled simulated
lineages.

Each row in the csv file represents one sampled lineage.  The dN/dS
ratios is extracted from the row along with the path to the FASTA
sample.  The fasta file is parsed and the distribution of codons is
calculated.  Each sample contributes one set of vertical bars to the
chart.  Samples with dN/dS ratio >900 (indicating failure) are grouped
apart from those with more "reasonable" dN/dS ratios.

You can see an example of the output from this script at 
https://cswarth.github.io/santa-dnds/codon_usage.html

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np

import plotly
from plotly import tools
import plotly.graph_objs as go

# cribbed from an otherise useless piece of code in biopython,
# https://github.com/biopython/biopython/blob/589bc075ace24cb29bcd69c9789963df7febeb4d/Bio/SeqUtils/CodonUsage.py#L13

CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

SynonymousCodons = {
    'CYS': ['TGT', 'TGC'],
    'ASP': ['GAT', 'GAC'],
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
    'GLN': ['CAA', 'CAG'],
    'MET': ['ATG'],
    'ASN': ['AAC', 'AAT'],
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
    'LYS': ['AAG', 'AAA'],
    'STOP': ['TAG', 'TGA', 'TAA'],
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
    'PHE': ['TTT', 'TTC'],
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
    'ILE': ['ATC', 'ATA', 'ATT'],
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
    'HIS': ['CAT', 'CAC'],
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
    'TRP': ['TGG'],
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
    'GLU': ['GAG', 'GAA'],
    'TYR': ['TAT', 'TAC']
}

# order dict keys so synonymous codons for Arginine appear first.
key = SynonymousCodons['ARG'] + [k for k in CodonsDict.keys() if k not in SynonymousCodons['ARG']]

# keys = [ 'CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA',
#          'CTT', 'ATG', 'AAG', 'AAA', 'ATC', 'AAC', 'ATA', 'CCT', 'ACT', 'AGC', 'ACA', 'CAT', 'AAT', 'ATT', 'CTG',
#          'CTA', 'CTC', 'CAC', 'ACG', 'CAA', 'AGT', 'CAG', 'CCG', 'CCC', 'TAT', 'GGT', 'TGT', 'CCA', 'TCT', 'GAT',
#          'TTT', 'TGC', 'GGG', 'TAG', 'GGA', 'TAA', 'GGC', 'TAC', 'TTC', 'TCG', 'TTA', 'TTG', 'TCC', 'GAA', 'TCA',
#          'GCA', 'GTA', 'GCC', 'GTC', 'GCG', 'GTG', 'GAG', 'GTT', 'GCT', 'ACC', 'TGA', 'GAC', 'TGG' ]

def count_codons(handle):
    # make the codon dictionary local
    codon_count = CodonsDict.copy()

    # iterate over sequence and count all the codons in the FastaFile.
    for cur_record in SeqIO.parse(handle, "fasta"):
        # make sure the sequence is lower case
        if str(cur_record.seq).islower():
            dna_sequence = str(cur_record.seq).upper()
        else:
            dna_sequence = str(cur_record.seq)
        for i in range(0, len(dna_sequence), 3):
            codon = dna_sequence[i:i + 3]
            if codon in codon_count:
                codon_count[codon] += 1
            else:
                raise TypeError("illegal codon %s in gene: %s" % (codon, cur_record.id))
    return codon_count

def count_codons_fasta(fname):
    with open(fname, 'r') as fh:
        return count_codons(fh)


def main(argv=sys.argv[1:]):
    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-o', '--output', default=None,
                            help='name of output graph file')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                            help='verbose output')
    a = parser.parse_args()


    
    def mkBarPlot(s):
        path = os.path.join(os.path.dirname(s['path']), 'sample_dedup.fa')
        d = count_codons_fasta(path)
        trace = go.Bar(x=keys, y=[d[k] for k in keys],
                           name=str(s['omega']),
                           marker=dict(
                               color="aqua" if s['omega'] < 999 else 'red',
                               line=dict(
                                   color="aqua" if s['omega'] < 999 else 'red',
                                   width=1.5),
                               ),
                      opacity=0.6, hoverinfo='none')
        return trace

    # Read the result file
    # group by dn/ds ratio  (those with riduculous values vs those without)
    # plot all the normal ones first followed by all the bad ones.
    # normal should be in on color, bad in another.
    # Name each trace after the file it came from plus the dnds value
    df = pd.read_csv("build/results.csv")

    df = df.sort_values(['omega'], ascending=[True])
    df = df.dropna(axis=0, how='any')
    traces = df.apply(mkBarPlot, axis=1)

    layout = go.Layout(title='Codon usage according to dN/dS ratio',
                        barmode='group',
                        bargroupgap=0.01)
    fig = go.Figure(data=list(traces), layout=layout)

    url = plotly.offline.plot(fig, show_link=False, auto_open=False)

    print(url)



if __name__ == "__main__":
   main(sys.argv[1:])
   


