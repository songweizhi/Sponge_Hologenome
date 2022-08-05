import os
from ete3 import Tree

tree_str_in = '/Users/songweizhi/Desktop/count_sister_taxa_wd/OG1968_iqtree.contree'

t_in = Tree(tree_str_in, format=1)
for leaf in t_in:
    leaf_name = leaf.name
    gnm_id = '_'.join(leaf_name.split('_')[:-1])
    leaf.name = gnm_id
    print(leaf_name)
    print(gnm_id)
print(t_in.write())






