"""
A work on generalizing tree topology metric to the unrooted case.
The basic idea is from "Metric on phylogenetic tree shapes". All the labels are
stored in node.
"""

from collections import defaultdict, deque
from dendropy import Tree
from gmpy2 import mpz
from math import sqrt
from sys import stderr


def label_parent(k, j):
    """
    Return a label for a node given labels for its children
    :return: 
    """
    if j > k:
        k, j = j, k
    return k * (k-1) // 2 + j + 1


def annotate_rooted_tree(tree):
    """
    Take a rooted phylogenetic tree and annotate it using the Colijn-Plazotta
    metric (precisely as in the paper)
    Modifies the tree object passed to it, returns nothing
    :param tree: 
    :return: 
    """
    assert isinstance(tree, Tree)
    for node in tree.postorder_node_iter():
        if not node.child_nodes():
            node.annotations['CP-label'].value = mpz(1)
        else:
            k, j = (x.annotations['CP-label'].value for x in node.child_nodes())
            label = label_parent(k, j)
            node.annotations['CP-label'].value = label


def get_root_label(tree):
    """
    Return the root value of a tree.
    If not annotated, annotate it first
    :param tree: 
    :return: 
    """
    r = tree.seed_node.annotations['CP-label'].value
    if r:
        return r
    else:
        annotate_rooted_tree(tree)
        return tree.seed_node.annotations['CP-label'].value


def get_rooted_vector(tree):
    """
    For an annotated rooted tree, collect labels into a vector
    :param tree: 
    :return: 
    """
    if not tree.seed_node.annotations['CP-label'].value:
        annotate_rooted_tree(tree)
    r = []
    for node in tree.postorder_node_iter():
        r.append(node.annotations['CP-label'].value)
    r = sorted(r)
    return r


# Unrooted version

def get_neighbours(node):
    """
    Get neighbours for a node
    :param node: 
    :return: 
    """
    r = []
    if node.child_nodes():
        r += [x for x in node.child_nodes()]
    if node.parent_node:
        r.append(node.parent_node)
    return r


def _process_node_wave(wave):
    """
    Use a bunch of nodes to try and give CPM-labels to their neighbours.

    Takes an iterable of nodes (assumes them to have the 'CPM-labels'
    annotation) and returns a tuple of nodes that got some new labels this run.
    :param wave:
    :return:
    """
    next_wave = set()
    for node in wave:
        for neighbour in node.annotations['CPM-labels'].value.keys():
            try:
                if neighbour.annotations['CPM-labels'].value[node] is None:
                    others = tuple(node.annotations['CPM-labels'].value[x]
                                   for x in node.annotations['CPM-labels'].value
                                   if x != neighbour and
                                   node.annotations['CPM-labels'].value[x])
                    if len(others) == 2:
                        v = label_parent(*others)
                        neighbour.annotations['CPM-labels'].value[node] = v
                        next_wave.add(neighbour)
            except KeyError as e:
                # There is a heisenbug when some nodes get incorrect keys for
                # 'CPM-labels' annotation. This is a poor man's debugger for it.
                # TODO: Remove this if the bug doesn't show up again.
                if neighbour.child_nodes():
                    print(neighbour.child_nodes())
                else:
                    print('No children')
                if neighbour.parent_node:
                    print(neighbour.parent_node)
                else:
                    print('No parent')
                print(neighbour.annotations['CPM-labels'])
                raise e
    return tuple(next_wave)
    

def annotate_unrooted_tree(tree):
    """
    Annotate a tree with three labels per node.
    :param tree: 
    :return: 
    """
    # Creating correct CPM-labels dicts for each node
    for node in tree.preorder_node_iter():
        if not node.annotations['CPM-labels'].value:
            node.annotations['CPM-labels'].value = {x: None for x in get_neighbours(node)}
    # A hack around a root node, which does not really exist, but which dendropy
    # creates anyway.
    root = tree.seed_node
    a, b = root.child_nodes()
    a.annotations['CPM-labels'].value[b] = a.annotations['CPM-labels'].value[root]
    del(a.annotations['CPM-labels'].value[root])
    b.annotations['CPM-labels'].value[a] = b.annotations['CPM-labels'].value[root]
    del(b.annotations['CPM-labels'].value[root])
    # Populating the initial CPM labels for leaf parents
    parents = set()
    for node in tree.leaf_node_iter():
        # Each node should have only one parent, but these being dicts I cannot
        # just address neighbours[0]
        for parent in node.annotations['CPM-labels'].value:
            parent.annotations['CPM-labels'].value[node] = mpz(1)
            parents.add(parent)
    # The basic idea is as follows: on every iteration, keep a list of nodes
    # that had at least one CPM-label set. On the next iteration, try to use
    # these nodes to give labels to some more nodes and keep a similar list.
    # Repeat until there are no nodes affected, which should be only when the
    # entire tree is marked up (not proven formally, but seems to work on the
    # smaller trees).
    # The initial set is composed from leaves' parents who get 1s from leaves.
    wave = tuple(parents)
    cont = True
    while cont:
        next_wave = _process_node_wave(wave)
        if len(next_wave) > 0:
            wave = next_wave
        else:
            cont = False
    # Root 'CPM-labels' is set to -1 (int) after annotation to indicate that the
    # tree was labeled for calling or not calling `annotate_unrooted_tree` from
    # `get_unrooted_vector`.
    root.annotations['CPM-labels'].value = -1


def get_unrooted_vector(tree):
    """
    Collect labels from an unrooted tree.
    Skips seed node as in the unrooted case it's not a real node, but implement
    :param tree: 
    :return: 
    """
    if not tree.seed_node.annotations['CPM-labels'].value == -1:
        annotate_unrooted_tree(tree)
    r = []
    for node in tree.preorder_node_iter():
        if node is not tree.seed_node:
            r += list(node.annotations['CPM-labels'].value.values())
    return list(sorted(r))


# Operations on vectors
def vector_dict(vector):
    """
    Return a vector dictionary.
    
    The dictionary is a defaultdict that, for each subtree label, returns the
    count of this label in vector. For absent labels, returns zero.
    :param vector:
    :return:
    """
    r = defaultdict(lambda: 0)
    for element in vector:
        r[element] += 1
    return r


def euclidean(d1, d2):
    """
    Return the euclidean distance between two vector dicts.
    
    Does not do any sort of normalization. That would probably make sense as the
    operation on vector(s) before any comparison between trees anyway.
    :param d1: vector dict for a tree
    :param d2: vector dict for a tree
    :return:
    """
    square_sum = 0
    s1 = set(d1.keys())
    key_set = s1.union(set(d2.keys()))
    for label in key_set:
        a = d1[label] if label in d1 else 0
        b = d2[label] if label in d2 else 0
        square_sum += (a - b)**2
    return sqrt(square_sum)
