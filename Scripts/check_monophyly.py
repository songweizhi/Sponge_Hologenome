import os
import glob
from ete3 import Tree
from db_files import sponge_archaeal_MAG_tax_dict
from db_files import gtdb_ar_gnm_tax_dict
from common_functions import gtdb_tax_str_parser


archaeal_mags_renamed_for_prokka_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka.txt'
pwd_tree_file_original                  = '/Users/songweizhi/Desktop/OG13686_iqtree.contree'
pwd_tree_file_renamed                   = '/Users/songweizhi/Desktop/OG13686_iqtree_renamed.contree'
taxon_rank                              = 'g'

monophyletic_groups = ['g__Nitrosarchaeum',
                       'g__Nitrosopelagicus',
                       'g__Nitrosopumilus',
                       'g__Nitrosotalea',
                       'g__Nitrosotenuis',
                       'g__Cenarchaeum',
                       'g__WTHC01', 'g__TA-20','g__VHBM01', 'g__DRGT01', 'g__PXYB01', 'g__VYCS01', 'g__JACEMX01']

# monophyletic_groups = ['g__Nitrosarchaeum',
#                        'g__Nitrosopelagicus',
#                        'g__Nitrosopumilus',
#                        'g__Nitrosotalea',
#                        'g__Nitrosotenuis',
#                        'g__Cenarchaeum']

mag_rename_dict = {}
for each_mag in open(archaeal_mags_renamed_for_prokka_txt):
    each_mag_split = each_mag.strip().split('\t')
    before_rename = each_mag_split[0]
    after_rename = each_mag_split[1]
    mag_rename_dict[after_rename] = before_rename


tree_folder         = '/Users/songweizhi/Desktop/0_contree_files'
tree_folder_renamed = '/Users/songweizhi/Desktop/0_contree_files_renamed'

tree_file_re = '%s/*_iqtree.contree' % tree_folder
tree_file_list = [os.path.basename(file_name) for file_name in glob.glob(tree_file_re)]


for tree_file in tree_file_list:

    HOG_id = tree_file[:-15]
    pwd_tree_file_original  = '%s/%s' % (tree_folder, tree_file)
    pwd_tree_file_renamed   = '%s/%s' % (tree_folder_renamed, tree_file)

    # rename the tree
    t = Tree(pwd_tree_file_original, format=1)
    for leaf in t:
        leaf_name_new = '_'.join(leaf.name.split('_')[:-1])
        leaf_name_new = mag_rename_dict.get(leaf_name_new, leaf_name_new)
        leaf.name = leaf_name_new
    t.write(format=1, outfile=pwd_tree_file_renamed)

    taxon_to_tree_mag_dict = {}
    t2 = Tree(pwd_tree_file_renamed, format=1)
    for leaf in t2:

        # get mag_id, mag_id_no_ext and mag_id_no_ext_no_source (gtdb or ncbi)
        mag_id_no_ext = leaf.name
        mag_id_no_ext_no_source = mag_id_no_ext
        if '.gtdb' in mag_id_no_ext_no_source:
            mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]
        if '.ncbi' in mag_id_no_ext_no_source:
            mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]

        # get mag_taxon_str
        mag_taxon_str = 'NA'
        if mag_id_no_ext_no_source in sponge_archaeal_MAG_tax_dict:
            mag_taxon_str = sponge_archaeal_MAG_tax_dict[mag_id_no_ext_no_source]
        if mag_id_no_ext_no_source in gtdb_ar_gnm_tax_dict:
            mag_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source]

        # get mag_taxon_str (GCA GCF things)
        if mag_taxon_str == 'NA':
            mag_id_no_ext_no_source_GCF = mag_id_no_ext_no_source.replace('GCA', 'GCF')
            if mag_id_no_ext_no_source_GCF in gtdb_ar_gnm_tax_dict:
                mag_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source_GCF]

        # get needed_taxon
        needed_taxon = 'NA'
        if mag_taxon_str != 'NA':
            needed_taxon = gtdb_tax_str_parser(mag_taxon_str, taxon_rank)

        if needed_taxon not in taxon_to_tree_mag_dict:
            taxon_to_tree_mag_dict[needed_taxon] = [leaf.name]
        else:
            taxon_to_tree_mag_dict[needed_taxon].append(leaf.name)

    checked_group = 0
    monophyletic_grp_num = 0
    for each_grp in taxon_to_tree_mag_dict:
        if each_grp in monophyletic_groups:
            grp_member = taxon_to_tree_mag_dict[each_grp]
            tree_monophyly_tuple = t.check_monophyly(values=grp_member, target_attr="name")
            #print('%s\t%s(%s)\t%s' % (tree_monophyly_tuple[1], each_grp, len(grp_member), grp_member))
            if tree_monophyly_tuple[1] == 'monophyletic':
                monophyletic_grp_num += 1
            checked_group += 1

    monophyletic_grp_pct = monophyletic_grp_num*100/checked_group
    monophyletic_grp_pct = float("{0:.2f}".format(monophyletic_grp_pct))
    print('%s\t%s\t%s/%s' % (HOG_id, monophyletic_grp_pct, monophyletic_grp_num, checked_group))
    #print()






