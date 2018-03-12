#! /usr/bin/env python3.6

import sys
from argparse import ArgumentParser
from dendropy import TreeList
from metrics import get_unrooted_vector, get_rooted_vector
from time import time

parser = ArgumentParser('Make MDS for tropical and USA flu trees')
parser.add_argument('-t', action='store_true',
                    help='Recompute vectors from trees. Reads vector files otherwise')
args = parser.parse_args()

if args.t:
    print('Processing tropical trees', file=sys.stderr)
    tropical_trees = TreeList.get_from_path('data/flu_tropical.nwk',
                                            schema='newick')
    print('Loaded trees', file=sys.stderr)
    with open('data/tropical_vectors', mode='w') as tv:
        for tree in tropical_trees:
            start = time()
            vector = get_unrooted_vector(tree)
            print(f'Processed a tree in {time()-start} seconds', file=sys.stderr)
            print(vector, file=tv)
    print('Done\nProcessing USA trees', file=sys.stderr)
    usa_trees = TreeList.get_from_path('data/flu_usa.nwk')
    print('Loaded trees', file=sys.stderr)
    with open('data/usa_vectors', mode='w') as uv:
        for tree in usa_trees:
            start = time()
            vector = get_unrooted_vector(tree)
            print(f'Processed a tree in {time()-start} seconds',
                  file=sys.stderr)
            print(vector, file=uv)
    print('Done', file=sys.stderr)
else:
    print('Do nothing')

