## CP-metric for urooted trees

This repo is a Python3 reimplementation of the distance metric between tree
topologies proposed in
[Colijn and Plazotta, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5790134/),
as well as the (yet) unpublished extension of this metric to unrooted trees.

The hashing-based optimization is provided for both rooted and unrooted
versions. If used, the numeric labels are replaced with their md5 hashes.
Although this optimization comes with a risk of collision, it is highly
recommended for trees of more than about 100 leaves. Without it, the memory
requirements grow extremely quickly.

### Dependencies

[Dendropy](https://dendropy.org) for the tree implementation,
[gmpy2](https://pypi.org/project/gmpy2/) for the arbitrary precision ints.

### Included scripts

See `scriptname -h` for script parameters

#### process_tree_set.py

Takes a newick (multi) tree file and produces a vector file (actually a txt with a
single number per string) for each tree in it.

#### mds_vectors.py

Takes multiple directories of vector files and performs MDS using euclidean
distances. Additionally depends on `numpy` to run and `matplotlib.pyplot` to
draw the image.

#### metrics.py

The module including the tree labeling and distance calculation implementations.
####

### License and citing

The code is provided under the terms of MIT license. If you find the rooted
variant useful, please consider citing the Colijn & Plazotta paper above.
