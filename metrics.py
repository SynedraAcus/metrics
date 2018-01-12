"""
A work on generalizing tree topology metric to the unrooted case.
The basic idea is from "Metric on phylogenetic tree shapes". All the labels are
stored in node.
"""

from dendropy import Tree


def label_parent(k, j):
    """
    Return a label for a node given labels for its children
    :return: 
    """
    if j > k:
        k, j = j, k
    return int(k * (k-1) / 2 + j + 1)


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


def annotate_unrooted_tree(tree):
    """
    Annotate a tree with three labels per node.
    :param tree: 
    :return: 
    """
    queue = []
    for node in tree.preorder_node_iter():
        if not node.annotations['CPM-labels'].value:
            node.annotations['CPM-labels'].value = {x: None for x in get_neighbours(node)}
        if not node.child_nodes():
            node.parent_node.annotations['CPM-labels'].value[node] = 1
        queue.append(node)
    # A hack around a root node, which does not really exist, but which dendropy
    # creates anyway.
    root = tree.seed_node
    a, b = root.child_nodes()
    a.annotations['CPM-labels'].value[b] = None
    del(a.annotations['CPM-labels'].value[root])
    b.annotations['CPM-labels'].value[a] = None
    del(b.annotations['CPM-labels'].value[root])
    queue.remove(root)
    # After this point, any time the node needs a value from its neignbour, the
    # neighbour has precisely two other neighbours.
    while len(queue) > 0:
        node = queue.pop(0)
        rm_permission = 0
        for neighbour in node.annotations['CPM-labels'].value.keys():
            if node.annotations['CPM-labels'].value[neighbour]:
                # The value is already set
                # Check if that neighbour is complete
                if None not in neighbour.annotations['CPM-labels'].value.values():
                    rm_permission += 1
                continue
            # Node can set its value from a particular neignbour only if said
            # neighbour has values from its other two neighbours already set.
            cousins = neighbour.annotations['CPM-labels'].value
            cousin_values = tuple(cousins[x] for x in cousins if x != node)
            if cousin_values and all(cousin_values):
                node.annotations['CPM-labels'].value[neighbour] = \
                    label_parent(*cousin_values)
        if not rm_permission == len(node.annotations['CPM-labels'].value)\
                or None in node.annotations['CPM-labels'].value.keys():
            # Either the node or one of its neighbours is incomplete. Cannot
            # leave the queue just yet
            # Probably there are some conditions when the node can be removed
            # earlier, but for this prototype it will leave the analysis iff
            # both it and all its neighbours are complete
            queue.append(node)


def get_unrooted_vector(tree):
    """
    Collect labels from an unrooted tree.
    Skips seed node as in the unrooted case it's not a real node, but implement
    :param tree: 
    :return: 
    """
    r = []
    for node in tree.preorder_node_iter():
        if node is not tree.seed_node:
            r += list(node.annotations['CPM-labels'].value.values())
    return list(sorted(r))
