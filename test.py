#! /usr/bin/env python3.6

from metrics import get_rooted_vector, annotate_unrooted_tree,\
    get_unrooted_vector, vector_dict, euclidean
from dendropy import Tree
from time import time

t2 = Tree.get(data='((A, B), (C, D));', schema='newick')
t1 = Tree.get(data='((A, (B, C)), (D,E));', schema='newick')
start = time()
annotate_unrooted_tree(t1)
print(time()-start)
# print(euclidean(vector_dict(get_unrooted_vector(t1)),
#                 vector_dict(get_unrooted_vector(t2))))

