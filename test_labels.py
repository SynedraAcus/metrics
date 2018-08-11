#! /usr/bin/env python3

# A pytest-compatible test suite for metrics.py

import pytest
from dendropy import Tree
from gmpy2 import mpz
from metrics import label_parent, get_rooted_vector, get_root_label, \
    get_unrooted_vector, vector_dict


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
    
    
def test_hashed_labels(tree):
    assert vector_dict(get_rooted_vector(tree, hashing=True)) ==\
           vector_dict(['c4ca4238a0b923820dcc509a6f75849b',
       'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
       'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
       'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
       'c4ca4238a0b923820dcc509a6f75849b', 'c81e728d9d4c2f636f067f89cc14862c',
       'c81e728d9d4c2f636f067f89cc14862c', 'c81e728d9d4c2f636f067f89cc14862c',
       'c81e728d9d4c2f636f067f89cc14862c', 'a87ff679a2f3e71d9181a67b7542122c',
       'a87ff679a2f3e71d9181a67b7542122c', '6512bd43d9caa6e02c990b0a82652dca'])
    assert get_root_label(tree, hashing=True) ==\
           '6512bd43d9caa6e02c990b0a82652dca'
    assert get_unrooted_vector(tree, hashing=True) ==\
           ['c4ca4238a0b923820dcc509a6f75849b',
        'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
        'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
        'c4ca4238a0b923820dcc509a6f75849b', 'c4ca4238a0b923820dcc509a6f75849b',
        'c4ca4238a0b923820dcc509a6f75849b', 'c81e728d9d4c2f636f067f89cc14862c',
        'c81e728d9d4c2f636f067f89cc14862c', 'c81e728d9d4c2f636f067f89cc14862c',
        'c81e728d9d4c2f636f067f89cc14862c', 'a87ff679a2f3e71d9181a67b7542122c',
        'a87ff679a2f3e71d9181a67b7542122c', '45c48cce2e2d7fbdea1afc51c7c6ad26',
        '45c48cce2e2d7fbdea1afc51c7c6ad26', '45c48cce2e2d7fbdea1afc51c7c6ad26',
        '45c48cce2e2d7fbdea1afc51c7c6ad26', 'a5771bce93e200c36f7cd9dfd0e5deaa',
        'a5771bce93e200c36f7cd9dfd0e5deaa', 'a5771bce93e200c36f7cd9dfd0e5deaa',
        'a5771bce93e200c36f7cd9dfd0e5deaa', 'a5771bce93e200c36f7cd9dfd0e5deaa',
        'a5771bce93e200c36f7cd9dfd0e5deaa', 'a5771bce93e200c36f7cd9dfd0e5deaa',
        'a5771bce93e200c36f7cd9dfd0e5deaa']
    


