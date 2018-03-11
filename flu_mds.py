#! /usr/bin/env python3

import sys
from argparse import ArgumentParser
from dendropy import TreeList
from metrics import get_unrooted_vector

parser = ArgumentParser('Make MDS for tropical and USA flu trees')
parser.add_argument('-t', action='store_true',
                    help='Recompute vectors from trees. Reads vector files otherwise')
args = parser.parse_args()

if args.t:
    print('Processing tropical trees', file=sys.stderr)
    tropical_trees = TreeList.get_from_path('data/flu_tropical.nwk',
                                            schema='newick')
    tropical_vectors = [get_unrooted_vector(x) for x in tropical_trees]
    with open('data/tropical_vectors', mode='w') as tv:
        for vector in tropical_vectors:
            print(vector, file=tv)
    print('Done\nProcessing USA trees', file=sys.stderr)
    usa_trees = TreeList.get_from_path('data/flu_usa.nwk')
    usa_vectors = [get_unrooted_vector(x) for x in usa_trees]
    with open('data/usa_vectors', mode='w') as uv:
        for vector in usa_vectors:
            print(vector, file=uv)
    print('Done', file=sys.stderr)
else:
    print('Do nothing')

