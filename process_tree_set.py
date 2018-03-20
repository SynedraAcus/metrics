#! /usr/bin/env python3

from argparse import ArgumentParser
from dendropy import TreeList
from metrics import get_rooted_vector, get_unrooted_vector
from multiprocessing import Pool
from os import getpid, cpu_count
from sys import stderr
from time import time


def write_tree(tree, func, filename):
    """
    Get a vector for a given tree and write it into a file.
    
    A separate function for multiprocessing later.
    :param tree:
    :return:
    """
    # Unpacking an argument tuple. Which is a tuple because of Pool.map()
    start = time()
    vector = func(tree)
    with open(filename, mode='w') as outfile:
        for x in vector:
            print(str(x), file=outfile)
    print('Processed vector {} in {} seconds by {}'.format(filename,
                                                           str(time()-start),
                                                           getpid()),
          file=stderr)
            
            
parser = ArgumentParser('Return CP- or CPM-vectors for a set of trees\n'+
                        'The vectors are written to a separate file each,\n'+
                        'named {tree_file}.tree{tree_number}.vector')
parser.add_argument('-t', type=str, help='Tree file in Newick format')
parser.add_argument('-u', action='store_true',
                    help='Produce unrooted (CPM) labelling')
parser.add_argument('--processes', type=int, default=0,
                    help='Number of processes. Defaults to processor number')
args = parser.parse_args()

start = time()
process_count = args.processes if args.processes else cpu_count()
print('Using {} processes'.format(process_count), file=stderr)
file_mask = args.t.split('.')[0]+'_tree{}.vector'
trees = TreeList.get_from_path(args.t, schema='newick')
print('Loaded {} trees'.format(len(trees)), file=stderr)
counter = 0
f = args.u and get_unrooted_vector or get_rooted_vector
func_args = [(trees[i], f, file_mask.format(str(i))) for i in range(len(trees))]
p = Pool(2)
_ = p.starmap(write_tree, func_args, chunksize=1)
print('Processed {} trees in {} seconds using {} processes'.format(
                                                                str(len(trees)),
                                                                time()-start,
                                                                process_count),
      file=stderr)
