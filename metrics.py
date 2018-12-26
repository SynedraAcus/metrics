"""
A work on generalizing tree topology metric to the unrooted case.
The basic idea is from "Metric on phylogenetic tree shapes". All the labels are
stored in node.
"""

from collections import defaultdict
from gmpy2 import mpz
from hashlib import md5
from math import sqrt

from dendropy import Tree
from networkx import DiGraph, topological_sort


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


def get_unrooted_vector(tree, hashing=False, annotation_method='graph'):
    """
    Collect labels from an unrooted tree.
    If the tree was not labeled before, it is labelled using the method supplied
    in the `annotation_method` kwarg.
    :param tree: a tree whose labels are to be returned
    :param hashing: if True, return MD5s of labels
    :param annotation_method: either 'wave' or 'graph'. 'graph' is faster, but
    requires networkx.
    :return:
    """
    functions = {'graph': label_graph_annotation,
                 'wave': wave_traversal_annotation,
                 'leaf': leaf_enumeration_annotation}
    if not tree.seed_node.annotations['CPM-labels'].value == -1:
        functions[annotation_method](tree, hashing=hashing)
    r = []
    for node in tree.preorder_node_iter():
        if node is not tree.seed_node:
            r += list(node.annotations['CPM-labels'].value.values())
    print(r)
    return list(sorted(r))


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

#  Two functions that define the wave traversal unrooted labeling


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
    

def wave_traversal_annotation(tree, hashing=False):
    """
    Annotate a tree with three labels per node.
    Uses a wave traversal algorithm which is, frankly, an ugly bunch of hacks.
    Not recommended except when networkx is unavailable.
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


# Label graph-based unrooted labeling


class LabelGraphNode:
    """
    A node for label graph. Has a home node, a target node, and a value
    """
    def __init__(self, home_node, target_node, value=None):
        self.home_node = home_node
        self.target_node = target_node
        self.value = value

    def __hash__(self):
        # A tuple of nodes is good enough for a unique label graph node
        return hash((self.home_node, self.target_node))


def label_graph_annotation(tree, hashing = True):
    """
    Annotate the tree using label graph.
    It's a graph such that every CPM label corresponds to a node and
    the labels that require other labels to be built are their descendants.
    :param tree: a Tree that has CPM-labels markup
    :param hashing: if True, return MD5 hashes of labels instead of themselves
    :return:
    """
    ### Walk over a tree, hanging graph nodes on their corresponding tree nodes
    for node in tree.preorder_node_iter():
        if not node.annotations['CPM-nodes'].value:
            node.annotations['CPM-nodes'].value = {
                    x: LabelGraphNode(node, x)
                    for x in get_neighbours(node)}
    # A hack around a root node, which does not really exist, but which dendropy
    # creates anyway.
    root = tree.seed_node
    a, b = root.child_nodes()
    a.annotations['CPM-nodes'].value[b] = a.annotations['CPM-nodes'].value[root]
    a.annotations['CPM-nodes'].value[b].target_node = b
    del(a.annotations['CPM-nodes'].value[root])
    b.annotations['CPM-nodes'].value[a] = b.annotations['CPM-nodes'].value[root]
    b.annotations['CPM-nodes'].value[a].target_node = a
    del(b.annotations['CPM-nodes'].value[root])
    # Populating the initial CPM labels for leaf parents
    for node in tree.leaf_node_iter():
        # Each node should have only one parent, but these being dicts I cannot
        # just address neighbours[0]
        for parent in node.annotations['CPM-nodes'].value:
            parent.annotations['CPM-nodes'].value[node].value = mpz(1)

    ### Walk over the tree again, collecting nodes and assembling them to a graph
    ### Necessary so that all nodes do exist by the time a graph is assembled
    label_graph = DiGraph()
    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        for label_node in node.annotations['CPM-nodes'].value.values():
            label_graph.add_node(label_node)
            for other_node in label_node.target_node.annotations['CPM-nodes'].value.values():
                if other_node.target_node is not node:
                    label_graph.add_edge(other_node, label_node)

    ### Traverse the graph, calculating values in nodes
    for label_node in topological_sort(label_graph):
        # Try and set value for self. Check if parents can be hashed
        if label_node.value is None:
            parents = tuple(label_graph.predecessors(label_node))
            # There are, by definition, two parents for each non-parentless node
            label_node.value = label_parent(parents[0].value, parents[1].value)
            if hashing:
                # Try and hash node's parents
                for parent in parents:
                    can_hash = True
                    for child in label_graph.successors(parent):
                        if child.value is None:
                            can_hash = False
                    if can_hash:
                        parent.value = md5(str(parent.value).
                                           encode(encoding='utf-8')).hexdigest()
                # If the node is not a parent, hash itself
                if list(label_graph.successors(label_node)) == []:
                    label_node.value = md5(str(label_node.value).
                                           encode(encoding='utf-8')).hexdigest()

    ### Walk over the tree the third time, collecting values from nodes
    #TODO: Discard nodes on the tree for memory saving and general clarity
    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            node.annotations['CPM-labels'] = -1
        else:
            node.annotations['CPM-labels'] =\
                {x: node.annotations['CPM-nodes'].value[x].value
                 for x in node.annotations['CPM-nodes'].value}


def recursive_label(node, direction, hashing):
    """
    An internal function for leaf_enumeration_annotation
    Essentially a lazy label calculation
    :param node:
    :param direction:
    :param hashing:
    :return:
    """
    # print(direction)
    # print(node.annotations)
    if node.annotations['CPM-labels'].value[direction] is None:
        next_nodes = [x for x in direction.annotations['CPM-labels'].value if x is not node]
        node.annotations['CPM-labels'].value[direction] = \
            label_parent(recursive_label(direction,
                                         next_nodes[0],
                                         hashing),
                         recursive_label(direction,
                                         next_nodes[1],
                                         hashing))
        if hashing:
            #TODO: Work out the correct hashing check
            can_hash = True
            for node2 in next_nodes:
                if node2.annotations['CPM-labels'].value[direction] is None:
                    can_hash = False
            if can_hash:
                direction.annotations['CPM-labels'].value[node] = md5(str(
                        direction.annotations['CPM-labels'].value[node]).encode(
                        encoding='utf-8')).hexdigest()
    return node.annotations['CPM-labels'].value[direction]


def leaf_enumeration_annotation(tree, hashing=False):
    """
    Annotate the unrooted tree using leaf enumeration.
    This approach uses a recursive function, so it may fail if the distance
    between two furthest nodes in a tree is above recursion limit.
    :param tree: a Tree that has CPM-labels markup
    :param hashing: if True, return MD5 hashes of labels instead of themselves
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
    del (a.annotations['CPM-labels'].value[root])
    b.annotations['CPM-labels'].value[a] = b.annotations['CPM-labels'].value[root]
    del (b.annotations['CPM-labels'].value[root])
    # Populating the initial CPM labels for leaf parents
    for node in tree.leaf_node_iter():
        # Each node should have only one parent, but these being dicts I cannot
        # just address neighbours[0]
        for parent in node.annotations['CPM-labels'].value:
            parent.annotations['CPM-labels'].value[node] = mpz(1)
    # For each leaf, compute the label directed towards the rest of the tree
    for node in tree.leaf_node_iter():
        for parent in node.annotations['CPM-labels'].value:
            recursive_label(node, parent, hashing=hashing)
            node.annotations['CPM-labels'].value[parent] = md5(str(
                        node.annotations['CPM-labels'].value[parent]).encode(
                        encoding='utf-8')).hexdigest()




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
