#! /usr/bin/env python3.6

from argparse import ArgumentParser
from dendropy import TreeList
from metrics import annotate_rooted_tree, annotate_unrooted_tree
from multiprocessing import Pool
from os import getpid, cpu_count
from sys import stderr
from time import time


def write_tree(tree, func, filename, hashing):
    """
    Get a vector for a given tree and write it into a file.
    :param tree:
    :return:
    """
    # Unpacking an argument tuple. Which is a tuple because of Pool.map()
    start = time()
    func(tree, hashing=hashing)
    with open(filename, mode='w') as outfile:
        for node in tree.preorder_node_iter():
            if func == annotate_unrooted_tree:
                if not isinstance(node.annotations['CPM-labels'].value, int):
                    # The int thing is for skipping the root value
                    for key in node.annotations['CPM-labels'].value:
                        print(str(node.annotations['CPM-labels'].value[key]),
                              file=outfile)
            else:
                print(str(node.annotations['CP-label'].value),
                      file=outfile)
    print('Processed vector {} in {} seconds by {}'.format(filename,
                                                           str(time()-start),
                                                           getpid()),
          file=stderr)
            
            
parser = ArgumentParser('Return CP- or CPM-vectors for a set of trees\n'+
                        'The vectors are written to a separate file each,\n'+
                        'named {tree_file}.tree_{tree_number}.vector')
parser.add_argument('-t', type=str, help='Tree file in Newick format')
parser.add_argument('-u', action='store_true',
                    help='Produce unrooted (CPM) labelling')
parser.add_argument('--hash', action='store_true',
                    help='Produce hashed labelling')
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
f = args.u and annotate_unrooted_tree or annotate_rooted_tree
func_args = [(trees[i], f, file_mask.format(str(i)), args.hash) for i in range(len(trees))]
p = Pool(process_count)
_ = p.starmap(write_tree, func_args, chunksize=1)
print('Processed {} trees in {} seconds using {} processes'.format(
                                                                str(len(trees)),
                                                                time()-start,
                                                                process_count),
      file=stderr)
