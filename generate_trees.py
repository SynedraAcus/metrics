#! /usr/bin/env python3

from argparse import ArgumentParser
from dendropy import TreeList, TaxonNamespace
from dendropy.simulate import treesim
import os

parser = ArgumentParser('Generate trees of a given size with different algos')
parser.add_argument('-n', type=int, help='Tree size', default=100)
parser.add_argument('-d', type=str, help='Output directory')
args = parser.parse_args()

if not os.path.isdir(args.d):
    os.mkdir(args.d)
os.chdir(args.d)
bd2 = TreeList([treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5,
                                         num_extant_tips=args.n,
                                         repeat_until_success=True)
                for _ in range(100)])
bd2.write_to_path('birth_death2.nwk', schema='newick')
bd5 = TreeList([treesim.birth_death_tree(birth_rate=1.0, death_rate=0.2,
                                         num_extant_tips=args.n,
                                         repeat_until_success=True)
                for _ in range(100)])
bd5.write_to_path('birth_death5.nwk', schema='newick')
taxa = TaxonNamespace(['T{}'.format(x) for x in range(1, args.n+1)])
king = TreeList([treesim.pure_kingman_tree(taxon_namespace=taxa) for _ in range(100)])
king.write_to_path('kingman.nwk', schema='newick')
