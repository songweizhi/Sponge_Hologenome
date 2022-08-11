import os
from ete3 import Tree


step_1_op_dir                   = '/Users/songweizhi/Desktop/count_sister_taxa_op'
hog_id_txt                      = '/Users/songweizhi/Desktop/HOG_id.txt'
TaxaCountStats_wd               = '/Users/songweizhi/Desktop/Input_folder_to_R'
TaxaCountStats_Rscript          = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/TaxaCountStats.R'

# define input files to R script
combined_contree_file           = '%s/combined.contree'                  % TaxaCountStats_wd
genes_to_remove_txt             = '%s/Genes_to_remove.txt'               % TaxaCountStats_wd
list_of_trees_txt               = '%s/List_of_trees.txt'                 % TaxaCountStats_wd
mapping_txt                     = '%s/mapping.txt'                       % TaxaCountStats_wd
marker_list_txt                 = '%s/MarkerList.txt'                    % TaxaCountStats_wd
combined_count_sister_taxa_op   = '%s/combined_count_sister_taxa_op.txt' % TaxaCountStats_wd
TaxaCountStats_op               = '%s/TaxaCountStats_output.txt'         % TaxaCountStats_wd


# get sorted hog list
hog_list = []
for each_hog in open(hog_id_txt):
    hog_list.append(each_hog.strip())
hog_list_sorted = sorted(hog_list)


cluster_to_domain_dict = {}
marker_list_txt_handle = open(marker_list_txt, 'w')
marker_list_txt_handle.write('MarkerID\n')
list_of_trees_txt_handle = open(list_of_trees_txt, 'w')
combined_contree_file_handle = open(combined_contree_file, 'w')
combined_count_sister_taxa_op_handle = open(combined_count_sister_taxa_op, 'w')
combined_count_sister_taxa_op_handle.write('MarkerID\tGroup_of_interest\tSister_taxa\tNormalized_sum_of_occurances\tsplits\tNormalized2_sum_of_occurances\tClusters\n')
for each_hog in hog_list_sorted:

    # write out to combined_count_sister_taxa_op
    pwd_count_sister_taxa_op_txt = '%s/count_sister_taxa_op/%s_iqtree_count_sister_taxa.txt' % (step_1_op_dir, each_hog)
    with open(pwd_count_sister_taxa_op_txt) as count_sister_taxa_op_txt_opened:
        for each_line in count_sister_taxa_op_txt_opened:
            combined_count_sister_taxa_op_handle.write('%s\t%s' % (each_hog, each_line))

    # write out to combined_contree_file
    pwd_renamed_contree_file = '%s/renamed_contree/%s_iqtree_renamed.contree' % (step_1_op_dir, each_hog)
    with open(pwd_renamed_contree_file, 'r') as pwd_renamed_contree_file_opened:
        combined_contree_file_handle.write(pwd_renamed_contree_file_opened.readline())

    # add to cluster_to_domain_dict
    t_in = Tree(pwd_renamed_contree_file, format=1)
    for leaf in t_in:
        leaf_name_split = leaf.name.split('|')
        print(leaf_name_split)
        cluster_to_domain_dict[leaf_name_split[0]] = leaf_name_split[1]

    # write out to marker_list_txt
    marker_list_txt_handle.write(each_hog + '\n')

    # write out to list_of_trees_txt
    list_of_trees_txt_handle.write(each_hog + '\n')

marker_list_txt_handle.close()
list_of_trees_txt_handle.close()
combined_contree_file_handle.close()
combined_count_sister_taxa_op_handle.close()


# prepare mapping_txt
mapping_txt_handle = open(mapping_txt, 'w')
mapping_txt_handle.write('Cluster\tDomain\n')
for each_cluster in cluster_to_domain_dict:
    mapping_txt_handle.write('%s\t%s\n' % (each_cluster, cluster_to_domain_dict[each_cluster]))
mapping_txt_handle.close()


# prepare genes_to_remove_txt
genes_to_remove_txt_handle = open(genes_to_remove_txt, 'w')
genes_to_remove_txt_handle.write('MarkerID\n')
genes_to_remove_txt_handle.close()


# run TaxaCountStats.R
get_TaxaCountStats_cmd = 'Rscript %s -t %s -l %s -g %s -x %s -s %s -r %s -o %s' % (TaxaCountStats_Rscript,
                                                                                   combined_contree_file,
                                                                                   list_of_trees_txt,
                                                                                   mapping_txt, marker_list_txt,
                                                                                   combined_count_sister_taxa_op,
                                                                                   genes_to_remove_txt,
                                                                                   TaxaCountStats_op)
print(get_TaxaCountStats_cmd)
os.system(get_TaxaCountStats_cmd)
