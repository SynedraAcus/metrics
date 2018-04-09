#! /usr/bin/env python3

from argparse import ArgumentParser
from glob import glob
from gmpy2 import mpz
from sklearn import manifold
from metrics import euclidean
import os
import numpy as np
import matplotlib.pyplot as plt


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
parser.add_argument('-d', type=str, nargs='*', help='Vector directories')
parser.add_argument('-o', type=str, help='Image output filename',
                    default='mds.svg')
parser.add_argument('-e', type=float, help='MDS epsilon',
		    default=1e-3)
args = parser.parse_args()

lengths = {}
vectors = []
for dir in args.d:
    if not os.path.exists(dir):
        raise ValueError('Inexistent directory {}'.format(dir))
    files = glob(dir+'/*.vector')
    lengths[dir] = len(files)
    for file in files:
        vectors.append(dict_from_file(file))

diss = np.ndarray(shape=(sum(lengths.values()), sum(lengths.values())),
                  dtype=np.float32)
for index, vector in enumerate(vectors):
    diss[index, index] = 0
    for index2 in range(index+1, len(vectors)):
        vector2 = vectors[index2]
        diss[index][index2] = np.float32(euclidean(vector, vector2))
        diss[index2][index] = diss[index, index2]
mds = manifold.MDS(dissimilarity='precomputed', eps=1e-9)
coords = mds.fit(diss).embedding_
# TODO: make some solution that does not crash on 6+ colours
colours = ['red', 'blue', 'green', 'yellow', 'green']
method_colours = {f: colours[d] for d, f in enumerate(args.d)}
r = 0
for x, d in enumerate(args.d):
    print(r, method_colours[d])
    plt.scatter(coords[r:r+len(vectors[x]), 0],
                coords[r:r+len(vectors[x]), 1],
                color=method_colours[d],
                label=d)
    plt.legend(scatterpoints=1)
    r += lengths[d]
plt.savefig(args.o)
