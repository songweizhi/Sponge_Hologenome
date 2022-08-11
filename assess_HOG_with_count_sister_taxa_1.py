import os
from ete3 import Tree
import multiprocessing as mp


def get_rename_dict(tree_str_in, mag_rename_dict, sponge_mag_tax_dict, gtdb_gnm_tax_dict):

    # rename dict: {'old_name':'new_name'}

    leaf_rename_dict = {}
    for leaf in Tree(tree_str_in, format=1):

        leaf_name_gnm = '_'.join(leaf.name.split('_')[:-1])
        leaf_name_gnm = mag_rename_dict.get(leaf_name_gnm, leaf_name_gnm)

        leaf_name_gnm_no_source = leaf_name_gnm
        if '.gtdb' in leaf_name_gnm_no_source:
            leaf_name_gnm_no_source = leaf_name_gnm[:-5]
        if '.ncbi' in leaf_name_gnm:
            leaf_name_gnm_no_source = leaf_name_gnm[:-5]

        # get mag_taxon_str
        gnm_taxon_str = 'NA'
        if leaf_name_gnm_no_source in sponge_mag_tax_dict:
            gnm_taxon_str = sponge_mag_tax_dict[leaf_name_gnm_no_source]
        if leaf_name_gnm_no_source in gtdb_gnm_tax_dict:
            gnm_taxon_str = gtdb_gnm_tax_dict[leaf_name_gnm_no_source]

        # get mag_taxon_str (GCA GCF things)
        if gnm_taxon_str == 'NA':
            mag_id_no_ext_no_source_GCF = leaf_name_gnm_no_source.replace('GCA', 'GCF')
            if mag_id_no_ext_no_source_GCF in gtdb_gnm_tax_dict:
                gnm_taxon_str = gtdb_gnm_tax_dict[mag_id_no_ext_no_source_GCF]

        gnm_taxon_str_no_space = gnm_taxon_str.replace(' ', '_')
        gnm_taxon_str_no_space = gnm_taxon_str_no_space.replace(';', '|')
        leaf_name_new = '%s|%s|strain__%s' % ('cluster_0', gnm_taxon_str_no_space, '_'.join(leaf.name.split('_')[:-1]))

        leaf_rename_dict[leaf.name] = leaf_name_new

    return leaf_rename_dict


def rename_tree(tree_str_in, rename_dict):

    # rename dict: {'old_name':'new_name'}

    t_in = Tree(tree_str_in, format=1)
    for leaf in t_in:
        leaf_name = leaf.name
        leaf_name_new = rename_dict[leaf_name]
        leaf.name = leaf_name_new

    return t_in.write()


def count_sister_taxa_worker(arg_list):

    mag_rename_dict                 = arg_list[0]
    sponge_archaeal_MAG_tax_dict    = arg_list[1]
    gtdb_ar_gnm_tax_dict            = arg_list[2]
    tree_ml                         = arg_list[3]
    ufboot_file                     = arg_list[4]
    target_label                    = arg_list[5]
    tree_ml_renamed                 = arg_list[6]
    ufboot_file_renamed             = arg_list[7]
    count_sister_taxa_op_txt        = arg_list[8]

    # rename ml tree
    tree_ml_renamed_handle = open(tree_ml_renamed, 'w')
    current_tree_rename_dict = get_rename_dict(tree_ml, mag_rename_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict)
    tree_ml_str_renamed = rename_tree(tree_ml, current_tree_rename_dict)
    tree_ml_renamed_handle.write(tree_ml_str_renamed + '\n')
    tree_ml_renamed_handle.close()

    # rename ufboot trees
    ufboot_file_renamed_handle = open(ufboot_file_renamed, 'w')
    for each_tree in open(ufboot_file):
        tree_str = each_tree.strip()
        current_tree_rename_dict = get_rename_dict(tree_str, mag_rename_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict)
        tree_str_renamed = rename_tree(tree_str, current_tree_rename_dict)
        ufboot_file_renamed_handle.write(tree_str_renamed + '\n')
    ufboot_file_renamed_handle.close()

    # run count_sister_taxa.py
    count_sister_taxa_cmd = 'python3 %s -ml %s -bs %s -l %s -out %s' % (
    count_sister_taxa_py, tree_ml_renamed, ufboot_file_renamed, target_label, count_sister_taxa_op_txt)
    print(count_sister_taxa_cmd)
    os.system(count_sister_taxa_cmd)


####################################################### file in ########################################################

from db_files import sponge_archaeal_MAG_tax_dict
from db_files import gtdb_ar_gnm_tax_dict

hog_id_txt                              = '/Users/songweizhi/Desktop/HOG_id.txt'
contree_and_ufboot_dir                  = '/Users/songweizhi/Desktop/contree_and_ufboot_files'
count_sister_taxa_py                    = '~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py'
archaeal_mags_renamed_for_prokka_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka.txt'
target_label                            = 'genus'
num_threads                             = 10

output_dir                              = '/Users/songweizhi/Desktop/count_sister_taxa_op'

########################################################################################################################

renamed_contree_dir      = '%s/renamed_contree'      % output_dir
renamed_ufboot_dir       = '%s/renamed_ufboot'       % output_dir
count_sister_taxa_op_dir = '%s/count_sister_taxa_op' % output_dir
os.mkdir(output_dir)
os.mkdir(renamed_contree_dir)
os.mkdir(renamed_ufboot_dir)
os.mkdir(count_sister_taxa_op_dir)

hog_list = []
for each_hog in open(hog_id_txt):
    hog_list.append(each_hog.strip())

mag_rename_dict = {}
for each_mag in open(archaeal_mags_renamed_for_prokka_txt):
    each_mag_split = each_mag.strip().split('\t')
    before_rename = each_mag_split[0]
    after_rename = each_mag_split[1]
    mag_rename_dict[after_rename] = before_rename

argument_lol = []
for og_id in hog_list:

    # define file name
    tree_ml                  = '%s/%s_iqtree.contree'               % (contree_and_ufboot_dir, og_id)
    ufboot_file              = '%s/%s_iqtree.ufboot'                % (contree_and_ufboot_dir, og_id)
    tree_ml_renamed          = '%s/%s_iqtree_renamed.contree'       % (output_dir, og_id)
    ufboot_file_renamed      = '%s/%s_iqtree.renamed.ufboot'        % (output_dir, og_id)
    count_sister_taxa_op_txt = '%s/%s_iqtree_count_sister_taxa.txt' % (output_dir, og_id)

    current_arg_list = [mag_rename_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict,
                        tree_ml, ufboot_file, target_label,
                        tree_ml_renamed, ufboot_file_renamed, count_sister_taxa_op_txt]

    argument_lol.append(current_arg_list)

# run with multiprocessing
pool = mp.Pool(processes=num_threads)
pool.map(count_sister_taxa_worker, argument_lol)
pool.close()
pool.join()
