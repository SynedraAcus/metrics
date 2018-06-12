#! /usr/bin/env python3

from argparse import ArgumentParser
from glob import glob
from collections import OrderedDict
from gmpy2 import mpz
from sklearn import manifold
from metrics import euclidean
import os
import numpy as np


def dict_from_file(file):
    """
    Return {value:count} dict for file lines
    :param file:
    :return:
    """
    r = {}
    for line in open(file):
        v = mpz(line)
        if v in r:
            r[v] += 1
        else:
            r[v] = 1
    return r
            

parser = ArgumentParser("""Run MDS of vector sets.
Expects each vector set to be stored as a collection of *.vector files in a
separate directory""")
parser.add_argument('-d', type=str, nargs='*', help=""""
                    Vector directories.
                    If absent, data will be loaded from --data_filename.
                    """)
parser.add_argument('-o', type=str, help='Image output filename',
                    default='mds.svg')
parser.add_argument('-z', action='store_true',
                    help='Process zero values')
parser.add_argument('--no_draw', action='store_true',
                    help="""
                    Do not draw MDS results. Store the embedding coordinates and
                    color data in text files instead, to be drawn at other
                    machine.
                    """)
parser.add_argument('--data_filename', type=str, default='distances',
                    help='Filename base for data dump')
args = parser.parse_args()

lengths = None
coords = None
if args.d:
    # If data were supplied
    lengths = OrderedDict()
    files = []
    for dir in args.d:
        if not os.path.exists(dir):
            raise ValueError('Nonetexistent directory {}'.format(dir))
        dir_files = glob(dir+'/*.vector')
        lengths[dir] = len(dir_files)
        files += dir_files
        # for file in files:
        #     vectors.append(dict_from_file(file))
    diss = np.ndarray(shape=(sum(lengths.values()), sum(lengths.values())),
                      dtype=np.float32)
    for index, file in enumerate(files):
        diss[index, index] = 0
        vector = dict_from_file(file)
        for index2 in range(index+1, len(files)):
            vector2 = dict_from_file(files[index2])
            diss[index][index2] = np.float32(euclidean(vector, vector2,
                                                       process_zeroes=args.z))
            diss[index2][index] = diss[index, index2]
    mds = manifold.MDS(dissimilarity='precomputed')
    coords = mds.fit(diss).embedding_
    if args.no_draw:
        # Data exist, but need to be dumped, not stored
        with open(args.data_filename+'.lengths', mode='w') as lenfile:
            for x in lengths:
                print(x+'\t'+str(lengths[x]), file=lenfile)
        np.savetxt(args.data_filename+'.coords', coords)
        print('Data written to {}'.format(args.data_filename))
        
if not args.no_draw:
    # If something is to be drawn, `lengths` and `coords` should be set.
    if not lengths and not coords:
        lengths = OrderedDict()
        for line in open(args.data_filename+'.lengths'):
            l = line.split('\t')
            lengths[l[0]] = int(l[1])
        coords = np.genfromtxt(args.data_filename+'.coords')
    print('Drawing')
    import matplotlib.pyplot as plt
    # TODO: make some solution that does not crash on 6+ colours
    colours = ['red', 'blue', 'green', 'yellow', 'green']
    method_colours = {f: colours[d] for d, f in enumerate(lengths)}
    r = 0
    for x, d in enumerate(lengths):
        print(r, method_colours[d])
        plt.scatter(coords[r:r+lengths[d], 0],
                    coords[r:r+lengths[d], 1],
                    color=method_colours[d],
                    label=d)
        plt.legend(scatterpoints=1)
        r += lengths[d]
    plt.savefig(args.o)
