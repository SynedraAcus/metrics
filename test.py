#! /usr/bin/env python3.6

from metrics import get_rooted_vector, annotate_unrooted_tree,\
    get_unrooted_vector, vector_dict, euclidean
from dendropy import Tree
from time import time

t2 = Tree.get(data='((A, B), (C, D));', schema='newick')
t1 = Tree.get(data='((A, (B, C)), (D,E));', schema='newick')
start = time()
print(get_rooted_vector(t2))
print(time()-start)
