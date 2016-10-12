#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''

Parse output of codeml to find dN/dS ratio for whole tree.

$ bin/parseresults.py -o build/results.csv $(find build -name results.txt)

'''
from __future__ import print_function

import os
import re
import sys
import collections
import argparse
import pandas as pd
from Bio import SeqIO

# A named tuple holding dN/dS ratio for a single pair of sequences
OmegaPairRecord = collections.namedtuple('OmegaPairRecord', 'path id1 id2 lnl dNdS dN dS')


# generator function to position filehandle for reading pairwise dN/dS ratios
# takes a filehandle to a results file from `codeml`.
# Yields the handle positioned at the start of the pairwise comparison section.
# after return the filehandle is exhausted.
def start_group(fhandle):
    # look for this phrase in the results file.
    start = re.compile('pairwise comparison, codon frequencies')
    for line in fhandle:
        if start.match(line):
            yield fhandle

            # only one pairwise section per file,
            # so return once we have found one.
            break
    return 


def next_record(fhandle):
    # pre-compile patterns used to extract portions of the record.
    # a pairwise record is spread across five lines in the results file.
    
    #     <id1>                                           <id2>
    # 2 (patient1_5000_228|1M|XXX|XXX|2011_11_10) ... 1 (patient1_5000_183|1M|XXX|XXX|2011_11_10)
    idpat = re.compile("\d+\s+\((?P<id1>\S+)\)\s+\.\.\.\s+\d+\s+\((?P<id2>\S+)\)$")
    
    #       <lnl>
    # lnL = -689.849468
    lnlpat = re.compile("lnL\s+=\s+(?P<lnl>-\d+.\d+)$")
    
    #                                            <dnds>       <dN>         <dS>
    # t= 0.1910  S=    96.3  N=   371.7  dN/dS=  0.4847  dN = 0.0523  dS = 0.1078
    dndspat = re.compile("t=\s+\d+\.\d+\s+S=\s+\d+\.\d+\s+N=\s+\d+.\d+\s+dN/dS=\s+(?P<dnds>\d+\.\d+)\s+dN\s+=\s+(?P<dN>\d+\.\d+)\s+dS\s+=\s+(?P<dS>\d+\.\d+)$")
    
    for line in fhandle:
        ids = idpat.match(line)
        if ids is not None:
            # read the next 4 lines as part of one record
            lnl = lnlpat.match(fhandle.next())
            fhandle.next()
            fhandle.next()
            dnds = dndspat.match(fhandle.next())
            yield Record( fhandle.name, ids.group('id1'), ids.group('id2'),
                          float(lnl.group('lnl')),
                          float(dnds.group('dnds')), float(dnds.group('dN')), float(dnds.group('dS')))


# parse single Omega value from results file that applies to whole tree.
# take a filehandle to a results.txt file, returns a float omega value, or None.
# filehandle is exhausted after return.
def extract_omega(fhandle):
    # regex pattern to match "omega (dN/dS) =  0.47112"
    omegapat = re.compile("omega\s+\(dN/dS\)\s+=\s+(?P<omega>\d+\.\d+)")
    data = fhandle.read()
    omega = omegapat.search(data)
    if omega is not None:
        return float(omega.group("omega"))
    return None

# A named tuple holding dN/dS ratio for a whole tree.
OmegaRecord = collections.namedtuple('OmegaRecord', 'path mut model fitness generation replicate omega count')


# generator function that yields one OmegaTreeRecord for every result file parsed.
# pass it a list (or generator) of result path names and get back a stream of OmegaTreeRecord tuples
def parse_results(files, verbose=False):
    pat = re.compile(r'(?P<path>[^/]*/[^_]*_(?P<mut>[^/]*)/(?P<model>[^/]*)/[^_]*_(?P<fitness>[^/]*)/[^_]*_(?P<replicate>[^/]*)/[^_]*_(?P<generation>[^/]*)/.*)')
    for path in files:
        if verbose:
            print("parsing {}".format(path))
        with open(path, 'r') as fh:
            omega = extract_omega(fh)
        with open(os.path.join(os.path.dirname(path), 'sample_dedup.fa'), 'r') as fh:
            records = list(SeqIO.parse(fh, "fasta"))
            count = len(records)

        # parse the path to extract pertinant information.
        #
        # some path components look like "var_value", e.g. "gen_20" or "fit_0.100",
        # where var is a string name and value is an integer or float.
        # strip off the first part and retain the numeric part.
        m = pat.match(path)
        if m:
            yield OmegaRecord(omega=omega, count=count, **m.groupdict())
        
def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''
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
    parser.add_argument('results', nargs='+', help='result files', type=existing_file)
    a = parser.parse_args()

    records = [r for r in parse_results(a.results, a.verbose)]
    df = pd.DataFrame.from_records(records, columns=records[0]._fields)
    if a.output is not None:
        df.to_csv(a.output, index=False)
    if a.verbose or a.output is None:
        print(df)
    
if __name__ == "__main__":
   main(sys.argv[1:])

