#! /usr/bin/env python3

from metrics import get_rooted_vector, annotate_unrooted_tree,\
    get_unrooted_vector, vector_dict, euclidean
from dendropy import Tree

t1 = Tree.get(data='((A, (B, C)), (D,E));', schema='newick')
t2 = Tree.get(data='((A, B), (C, D));', schema='newick')
annotate_unrooted_tree(t1)
# annotate_unrooted_tree(t2)
print(euclidean(vector_dict(get_unrooted_vector(t1)),
                vector_dict(get_unrooted_vector(t2))))

