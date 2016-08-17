#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''


'''
from __future__ import print_function

import os
import re
import sys
import collections
import argparse
import pandas as pd

Record = collections.namedtuple('Record', 'id1 id2 lnl dNdS dN dS')
Orecord = collections.namedtuple('Orecord', 'fitness rep omega')


def start_group(fhandle):
    start = re.compile('pairwise comparison, codon frequencies')
    for line in fhandle:
        if start.match(line):
            yield fhandle
    return 


def next_record(fhandle):
    # pre-compile patterns we will use to match portions of the record.
    
    # 2 (patient1_5000_228|1M|XXX|XXX|2011_11_10) ... 1 (patient1_5000_183|1M|XXX|XXX|2011_11_10)
    idpat = re.compile("\d+\s+\((?P<id1>\S+)\)\s+\.\.\.\s+\d+\s+\((?P<id2>\S+)\)$")
    # lnL = -689.849468
    lnlpat = re.compile("lnL\s+=\s+(?P<lnl>-\d+.\d+)$")
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
            yield Record( ids.group('id1'), ids.group('id2'),
                          float(lnl.group('lnl')),
                          float(dnds.group('dnds')), float(dnds.group('dN')), float(dnds.group('dS')))


def extract_omega(fhandle):
    # regex pattern to match "omega (dN/dS) =  0.47112"
    omegapat = re.compile("omega\s+\(dN/dS\)\s+=\s+(?P<omega>\d+\.\d+)")
    data = fhandle.read()
    omega = omegapat.search(data)
    if omega is not None:
        return float(omega.group("omega"))
    return None

ORecord = collections.namedtuple('ORecord', 'fitness rep omega')

def parse_results(files):
    for path in files:
        # parse the path to extract pertinant information
        tokens = path.split(os.path.sep)
        tokens = [re.sub(r"[^_]*_", '', t) for t in tokens]
        with open(path, 'r') as fh:
            omega = extract_omega(fh)
            yield ORecord( float(tokens[2]), int(tokens[3]), omega)
        
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

    parser.add_argument('-o', '--output', default='graph.jpg',
                            help='name of output graph file')
    parser.add_argument('results', nargs='+', help='result files', type=existing_file)
    a = parser.parse_args()

    records = [r for r in parse_results(a.results)]
    df = pd.DataFrame.from_records(records, columns=ORecord._fields)
    print(df)
    
if __name__ == "__main__":
   main(sys.argv[1:])

