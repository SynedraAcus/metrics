from metrics import get_rooted_vector, annotate_unrooted_tree,\
    get_unrooted_vector
from dendropy import Tree

t = Tree.get(data='((A, (B, C)), (D,E));', schema='newick')
triangle = Tree.get(data='(A, (B, C));', schema='newick')
# print(get_rooted_vector(t))
annotate_unrooted_tree(t)
print(get_unrooted_vector(t))
