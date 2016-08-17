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

    tree = dendropy.Tree.get(
        path=a.input,
        schema="nexus",
        label=None,
        taxon_namespace=None,
        collection_offset=None,
        tree_offset=None,
        rooting="force-unrooted",
        edge_length_type=float,
        suppress_edge_lengths=False,
        extract_comment_metadata=True,
        store_tree_weights=False,
        finish_node_fn=None,
        case_sensitive_taxon_labels=False,
        preserve_underscores=False,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=False,
        terminating_semicolon_required=True,
        ignore_unrecognized_keyword_arguments=False
        )
    tree.write(
        path=a.output,
        schema='newick',
        suppress_leaf_taxon_labels=False,
        suppress_leaf_node_labels=True,
        suppress_internal_taxon_labels=False,
        suppress_internal_node_labels=False,
        suppress_rooting=True,
        suppress_edge_lengths=False,
        unquoted_underscores=False,
        preserve_spaces=False,
        store_tree_weights=False,
        taxon_token_map=None,
        suppress_annotations=True,
        annotations_as_nhx=False,
        suppress_item_comments=True,
        node_label_element_separator=' ',
        node_label_compose_fn=None,
        edge_label_compose_fn=None,
        real_value_format_specifier='',
        ignore_unrecognized_keyword_arguments=False
        )

if __name__ == "__main__":
   main(sys.argv[1:])
   


