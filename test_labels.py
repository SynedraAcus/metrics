#! /usr/bin/env python3

# A pytest-compatible test suite for metrics.py

import pytest
from dendropy import Tree
from gmpy2 import mpz
from metrics import label_parent, get_rooted_vector, get_root_label, \
    get_unrooted_vector


@pytest.fixture
def tree():
    return Tree.get_from_string('(((A, B), (C, D)), ((E, F), (G, H)));',
                                 schema='newick')


def test_sums():
    """
    Tests the parent node calculation
    :return:
    """
    assert label_parent(1, 2) == 3
    assert label_parent (1, 4) == 8
    # Should ignore arg order
    assert label_parent(4, 1) == 8


def test_rooted_labels(tree):
    assert get_rooted_vector(tree) == [mpz(x) for x in [1, 1, 1, 1, 1, 1, 1, 1,
                                                       2, 2, 2, 2,
                                                       4, 4,
                                                       11]]
    assert get_root_label(tree) == 11


def test_unrooted_labels(tree):
    assert get_unrooted_vector(tree) == [mpz(x) for x in
                                         [1, 1, 1, 1, 1, 1, 1, 1,
                                         2, 2, 2, 2,
                                         4, 4,
                                         9, 9, 9, 9,
                                         38, 38, 38, 38, 38, 38, 38, 38]]
