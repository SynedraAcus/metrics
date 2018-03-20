#! /usr/bin/env python3

from argparse import ArgumentParser
from dendropy import TreeList
from metrics import get_rooted_vector, get_unrooted_vector
from sys import stderr
from time import time


def write_tree(tree, function, filename):
    """
    Get a vector for a given tree and write it into a file.
    
    A separate function for multiprocessing later.
    :param tree:
    :return:
    """
    start = time()
    vector = function(tree)
    with open(filename, mode='w') as outfile:
        for x in vector:
            print(str(x), file=outfile)
    print('Processed a tree in {}'.format(str(time()-start)), file=stderr)
            
            
parser = ArgumentParser('Return CP- or CPM-vectors for a set of trees\n'+
                        'The vectors are written to a separate file each,\n'+
                        'named {tree_file}.tree{tree_number}.vector')
parser.add_argument('-t', type=str, help='Tree file in Newick format')
parser.add_argument('-u', action='store_true',
                    help='Produce unrooted (CPM) labelling')
parser.add_argument('-d', type=str, help='Output directory')
args = parser.parse_args()

file_mask = args.t.split('.')[0]+'_tree{}.vector'
trees = TreeList.get_from_path(args.t, schema='newick')
print('Loaded {} trees'.format(len(trees)), file=stderr)
counter = 0
f = args.u and get_unrooted_vector or get_rooted_vector
for tree in trees:
    counter += 1
    write_tree(tree, f, file_mask.format(str(counter)))
