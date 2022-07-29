import os
import glob
from ete3 import Tree
from db_files import sponge_archaeal_MAG_tax_dict
from db_files import gtdb_ar_gnm_tax_dict
from common_functions import gtdb_tax_str_parser


'''
cd /Users/songweizhi/Desktop/iqtree2_wd

# get the tree and bootstrap files (takes about XXX mins)
/Users/songweizhi/Software/iqtree/iqtree-2.2.0-MacOSX/bin/iqtree2 --wbtl --bnni -st AA -m LG+C60+F -T 10 -B 1000 -alrt 1000 -s example2.phy -pre example2_iqtree

# run count_sister_taxa.py
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree.contree example2_iqtree.ufboot outputFile_example1
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree2.contree example2_iqtree2.ufboot outputFile_example2wrong
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree3.contree example2_iqtree3.ufboot outputFile_example3wrong

cd /Users/songweizhi/Desktop/count_sister_taxa_wd
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py OG1423_iqtree.contree OG1423_iqtree.ufboot OG1423_iqtree_count_sister_taxa_op

cd /Users/songweizhi/Desktop/count_sister_taxa_wd
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py OG1423_iqtree_renamed.contree OG1423_iqtree_renamed.contree OG1423_iqtree_count_sister_taxa_op

'''


archaeal_mags_renamed_for_prokka_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka.txt'
pwd_tree_file_original = '/Users/songweizhi/Desktop/count_sister_taxa_wd/tree_in/OG1423_iqtree.contree'
pwd_tree_file_renamed  = '/Users/songweizhi/Desktop/count_sister_taxa_wd/tree_in/OG1423_iqtree_renamed.contree'


mag_rename_dict = {}
for each_mag in open(archaeal_mags_renamed_for_prokka_txt):
    each_mag_split = each_mag.strip().split('\t')
    before_rename = each_mag_split[0]
    after_rename = each_mag_split[1]
    mag_rename_dict[after_rename] = before_rename


# rename the tree
t = Tree(pwd_tree_file_original, format=1)
for leaf in t:

    leaf_name_gene = leaf.name
    leaf_name_gnm = '_'.join(leaf.name.split('_')[:-1])
    leaf_name_gnm = mag_rename_dict.get(leaf_name_gnm, leaf_name_gnm)

    leaf_name_gnm_no_source = leaf_name_gnm
    if '.gtdb' in leaf_name_gnm_no_source:
        leaf_name_gnm_no_source = leaf_name_gnm[:-5]
    if '.ncbi' in leaf_name_gnm:
        leaf_name_gnm_no_source = leaf_name_gnm[:-5]

    # get mag_taxon_str
    gnm_taxon_str = 'NA'
    if leaf_name_gnm_no_source in sponge_archaeal_MAG_tax_dict:
        gnm_taxon_str = sponge_archaeal_MAG_tax_dict[leaf_name_gnm_no_source]
    if leaf_name_gnm_no_source in gtdb_ar_gnm_tax_dict:
        gnm_taxon_str = gtdb_ar_gnm_tax_dict[leaf_name_gnm_no_source]

    # get mag_taxon_str (GCA GCF things)
    if gnm_taxon_str == 'NA':
        mag_id_no_ext_no_source_GCF = leaf_name_gnm_no_source.replace('GCA', 'GCF')
        if mag_id_no_ext_no_source_GCF in gtdb_ar_gnm_tax_dict:
            gnm_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source_GCF]

    gnm_taxon_str_no_space = gnm_taxon_str.replace(' ', '_')
    gnm_taxon_str_no_space = gnm_taxon_str_no_space.replace(';', '|')
    leaf_name_new = '%s|%s|strain__%s' % ('cluster_0', gnm_taxon_str_no_space, '_'.join(leaf.name.split('_')[:-1]))
    leaf.name = leaf_name_new
    print(leaf_name_gnm)
    print(gnm_taxon_str)
    print(leaf_name_new)
    print()





    #leaf_name_new = '_'.join(leaf.name.split('_')[:-1])
    #leaf_name_new = mag_rename_dict.get(leaf_name_new, leaf_name_new)
    #leaf.name = leaf_name_new
t.write(format=1, outfile=pwd_tree_file_renamed)











