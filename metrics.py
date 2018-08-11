"""
A work on generalizing tree topology metric to the unrooted case.
The basic idea is from "Metric on phylogenetic tree shapes". All the labels are
stored in node.
"""

from collections import defaultdict, deque
from dendropy import Tree
from gmpy2 import mpz, to_binary
from hashlib import md5
from math import sqrt


def label_parent(k, j):
    """
    Return a label for a node given labels for its children
    :return: 
    """
    if j > k:
        k, j = j, k
    return k * (k-1) // 2 + j + 1


def annotate_rooted_tree(tree, hashing=False):
    """
    Take a rooted phylogenetic tree and annotate it using the Colijn-Plazotta
    metric (precisely as in the paper)
    Modifies the tree object passed to it, returns nothing
    :param tree: A tree to be annotated
    :param hashing: Boolean. If True, the labels are set to MD5 hashes of the
    CP-labels instead of the labels themselves. This is useful for large trees
    where labels get prohibitively high.
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
            if hashing:
                # Hashing children values using MD5 if necessary
                # The children values are not gonna be necessary anymore if the
                # current node value was set
                for child in node.child_nodes():
                    child.annotations['CP-label'].value =\
                            md5(str(child.annotations['CP-label'].value).
                                encode(encoding='utf-8')).hexdigest()
    # Hashing root separately because it has no parent
    if hashing:
        tree.seed_node.annotations['CP-label'].value = \
            md5(str(tree.seed_node.annotations['CP-label'].value).encode(
                encoding='utf-8')).hexdigest()


def get_root_label(tree, hashing=False):
    """
    Return the root value of a tree.
    If not annotated, annotate it first
    :param tree: A tree whose label is to be produced
    :param hashing: if True, return MD5 hash of the label
    :return: 
    """
    r = tree.seed_node.annotations['CP-label'].value
    if r:
        return r
    else:
        annotate_rooted_tree(tree, hashing=hashing)
        return tree.seed_node.annotations['CP-label'].value


def get_rooted_vector(tree, hashing=False):
    """
    For an annotated rooted tree, collect labels into a vector
    :param tree: a tree whose label vector is to be produced
    :param hashing: if True, return MD5s of labels
    :return: 
    """
    if not tree.seed_node.annotations['CP-label'].value:
        annotate_rooted_tree(tree, hashing=hashing)
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


def _process_node_wave(wave, hashing=False):
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
            if neighbour.annotations['CPM-labels'].value[node] is None:
                others = tuple(node.annotations['CPM-labels'].value[x]
                               for x in node.annotations['CPM-labels'].value
                               if x != neighbour and
                               node.annotations['CPM-labels'].value[x])
                if len(others) == len(node.annotations['CPM-labels'].value) - 1:
                    v = label_parent(*others)
                    neighbour.annotations['CPM-labels'].value[node] = v
                    next_wave.add(neighbour)
        if hashing:
            # The value can be hashed iff the values are set for all the
            # neighbouring branches pointing to this node, except maybe the one
            # with the value to be hashed.
            # The latter appears to always be set in the tests, but I don't have
            # a formal proof
            for neighbour in node.annotations['CPM-labels'].value.keys():
                count = 0
                for neighbour2 in node.annotations['CPM-labels'].value.keys():
                    if neighbour2 == neighbour:
                        continue
                    if neighbour2.annotations['CPM-labels'].value[node]:
                        count += 1
                if count >= len(node.annotations['CPM-labels'].value.keys()) -1\
                        and not isinstance(
                        node.annotations['CPM-labels'].value[neighbour], str):
                    node.annotations['CPM-labels'].value[neighbour] = md5(str(
                        node.annotations['CPM-labels'].value[neighbour]).encode(
                        encoding='utf-8')).hexdigest()
    return tuple(next_wave)
    

def annotate_unrooted_tree(tree, hashing=False):
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
        next_wave = _process_node_wave(wave, hashing=hashing)
        if len(next_wave) > 0:
            wave = next_wave
        else:
            cont = False
    # Root 'CPM-labels' is set to -1 (int) after annotation to indicate that the
    # tree was labeled for calling or not calling `annotate_unrooted_tree` from
    # `get_unrooted_vector`.
    root.annotations['CPM-labels'].value = -1


def get_unrooted_vector(tree, hashing=False):
    """
    Collect labels from an unrooted tree.
    Skips seed node as in the unrooted case it's not a real node, but implement
    :param tree: a tree whose labels are to be returned
    :param hashing: if True, return MD5s of labels
    :return: 
    """
    if not tree.seed_node.annotations['CPM-labels'].value == -1:
        annotate_unrooted_tree(tree, hashing=hashing)
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


def euclidean(d1, d2, process_zeroes = True):
    """
    Return the euclidean distance between two vector dicts.
    
    Does not do any sort of normalization. That would probably make sense as the
    operation on vector(s) before any comparison between trees anyway.
    :param d1: vector dict for a tree
    :param d2: vector dict for a tree
    :param process_zeroes: if True, values present in one vector but not another
    are treated as zeroes for the latter vector. If False, they are ignored
    :return:
    """
    square_sum = 0
    if process_zeroes:
        s1 = set(d1.keys())
        key_set = s1.union(set(d2.keys()))
        for label in key_set:
            a = d1[label] if label in d1 else 0
            b = d2[label] if label in d2 else 0
            square_sum += (a - b)**2
    else:
        for label in d1:
            if label in d2:
                square_sum += (d1[label]-d2[label])**2
    return sqrt(square_sum)
