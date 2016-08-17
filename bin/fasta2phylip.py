#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import sys
import dendropy
import argparse

def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('input', help='fasta file')
    parser.add_argument('output', help='phylip output filename')

    a = parser.parse_args(argv)

    d = dendropy.DataSet.get(
        path=a.input,
        schema="fasta",
        label=None,
        taxon_namespace=None,
        ignore_unrecognized_keyword_arguments=False,
        data_type="dna",
        )
    d.write(
        path=a.output,
        schema="phylip",
        strict=False,
        spaces_to_underscores=False,
        force_unique_taxon_labels=False,
        suppress_missing_taxa=False,
        ignore_unrecognized_keyword_arguments=False,
        )

if __name__ == "__main__":
   main(sys.argv[1:])
   


