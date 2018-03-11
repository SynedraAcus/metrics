"""
A work on generalizing tree topology metric to the unrooted case.
The basic idea is from "Metric on phylogenetic tree shapes". All the labels are
stored in node.
"""

from collections import defaultdict, deque
from dendropy import Tree
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
            node.annotations['CP-label'].value = 1
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
    # It only makes sense to process nodes that were changed
    # on a previous iteration
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
                        stderr.write('Setting value\n')
                        neighbour.annotations['CPM-labels'].value[node] = v
                        next_wave.add(neighbour)
                    else:
                        stderr.write('Skipping\n')
            except KeyError as e:
                # There is a heisenbug when some nodes get incorrect keys for
                # 'CPM-labels' annotation. This is a poor man's debugger for it.
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
    queue = deque()
    parents = set()
    for node in tree.preorder_node_iter():
        if not node.annotations['CPM-labels'].value:
            node.annotations['CPM-labels'].value = {x: None for x in get_neighbours(node)}
        if not node.child_nodes():
            node.parent_node.annotations['CPM-labels'].value[node] = 1
            parents.add(node.parent_node)
    # for node in parents:
    #     if len([x for x in node.annotations['CPM-labels'].value.values() if x]) == 2:
    #         queue.append(node)
    # A hack around a root node, which does not really exist, but which dendropy
    # creates anyway.
    root = tree.seed_node
    a, b = root.child_nodes()
    a.annotations['CPM-labels'].value[b] = a.annotations['CPM-labels'].value[root]
    del(a.annotations['CPM-labels'].value[root])
    b.annotations['CPM-labels'].value[a] = b.annotations['CPM-labels'].value[root]
    del(b.annotations['CPM-labels'].value[root])
    # if root in queue:
    #     queue.remove(root)
    if root in parents:
        parents.remove(root)
    # print(queue)
    wave = tuple(parents)
    cont = True
    while cont:
        print(f'Wave of {len(wave)} node(s)', file=stderr)
        print(wave, file=stderr)
        next_wave = _process_node_wave(wave)
        if len(next_wave) > 0:
            wave = next_wave
        else:
            cont = False
    # # After this point, any time the node needs a value from its neignbour, the
    # # neighbour has precisely two other neighbours.
    # while len(queue) > 0:
    #     node = queue.popleft()
    #     stderr.write('Processing\n')
    #     rm_permission = 0
    #     for neighbour in node.annotations['CPM-labels'].value.keys():
    #         if neighbour.annotations['CPM-labels'].value[node] is not None:
    #             # This node has already set the value for neighbour in question
    #             rm_permission += 1
    #         else:
    #             others = tuple(node.annotations['CPM-labels'].value[x]
    #                            for x in node.annotations['CPM-labels'].value
    #                            if x != neighbour and node.annotations['CPM-labels'].value[x])
    #             if len(others) == 2:
    #                 stderr.write('Setting value\n')
    #                 neighbour.annotations['CPM-labels'].value[node] = \
    #                         label_parent(*others)
    #                 if len([x for x in neighbour.annotations['CPM-labels'].value.values() if x]) == 2:
    #                     queue.append(neighbour)
    #                 rm_permission += 1
    #         # if node.annotations['CPM-labels'].value[neighbour]:
    #         #     # The value is already set
    #         #     # Check if that neighbour is complete
    #         #     if None not in neighbour.annotations['CPM-labels'].value.values():
    #         #         rm_permission += 1
    #         #     continue
    #         # # Node can set its value from a particular neignbour only if said
    #         # # neighbour has values from its other two neighbours already set.
    #         # cousins = neighbour.annotations['CPM-labels'].value
    #         # cousin_values = tuple(cousins[x] for x in cousins if x != node)
    #         # if cousin_values and all(cousin_values):
    #         #     node.annotations['CPM-labels'].value[neighbour] = \
    #         #         label_parent(*cousin_values)
    #     if rm_permission < 3:
    #             # or None in node.annotations['CPM-labels'].value.keys():
    #         # Either the node or one of its neighbours is incomplete. Cannot
    #         # leave the queue just yet
    #         # Probably there are some conditions when the node can be removed
    #         # earlier, but for this prototype it will leave the analysis iff
    #         # both it and all its neighbours are complete
    #         queue.append(node)
    #     else:
    #         stderr.write('Removing node from queue\n')
    # Root 'CPM-labels' is set to -1 (int) after annotation to indicate that the
    # tree was labeled for calling `annotate_unrooted_tree` from
    # `get_unrooted_vector` by default
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
    print(d1, d2, key_set)
    for label in key_set:
        square_sum += (d1[label] - d2[label])**2
    return sqrt(square_sum)
