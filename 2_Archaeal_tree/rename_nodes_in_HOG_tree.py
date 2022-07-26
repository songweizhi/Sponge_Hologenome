import os
import glob
from ete3 import Tree

# file in
tree_folder                             = '/Users/songweizhi/Desktop/OrthologousGroupsFasta_cov_85_iqtree_UFBoot_contree_files'
archaeal_mags_renamed_for_prokka_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka.txt'

# file out
tree_folder_renamed                     = '/Users/songweizhi/Desktop/OrthologousGroupsFasta_cov_85_iqtree_UFBoot_contree_files_renamed'

mag_rename_dict = {}
for each_mag in open(archaeal_mags_renamed_for_prokka_txt):
    each_mag_split = each_mag.strip().split('\t')
    before_rename = each_mag_split[0]
    after_rename = each_mag_split[1]
    mag_rename_dict[after_rename] = before_rename

tree_file_re = '%s/*_iqtree.contree' % tree_folder
tree_file_list = [os.path.basename(file_name) for file_name in glob.glob(tree_file_re)]

for tree_file in tree_file_list:

    pwd_tree_file     = '%s/%s' % (tree_folder, tree_file)
    pwd_tree_file_out = '%s/%s' % (tree_folder_renamed, tree_file)

    t = Tree(pwd_tree_file, format=1)
    for leaf in t:
        leaf_name_new = '_'.join(leaf.name.split('_')[:-1])
        leaf_name_new = mag_rename_dict.get(leaf_name_new, leaf_name_new)
        leaf.name = leaf_name_new
    t.write(format=1, outfile=pwd_tree_file_out)
