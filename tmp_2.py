import ast


file_in = '/Users/songweizhi/Desktop/marker_gene_stats_STY8taxa.tsv'

marker_set = set()
marker_num_dict_of_dict = {}
for each in open(file_in):

    each_split = each.strip().split('\t')

    bin_id = each_split[0]
    bin_marker_dict_str = each_split[1]
    bin_marker_dict = ast.literal_eval(bin_marker_dict_str)
    marker_num_dict = {}
    for each_ctg in bin_marker_dict:
        for each_marker in bin_marker_dict[each_ctg]:
            marker_set.add(each_marker)
            if each_marker not in marker_num_dict:
                marker_num_dict[each_marker] = 1
            else:
                marker_num_dict[each_marker] += 1

    # for each_marker in marker_num_dict:
    #     print('%s\t%s' % (marker_num_dict[each_marker], each_marker))
    marker_num_dict_of_dict[bin_id] = marker_num_dict


marker_set_sorted = sorted([i for i in marker_set])


maker_data_matrix = '/Users/songweizhi/Desktop/marker_gene_stats_STY8taxa_matrix.tsv'

maker_data_matrix_handle = open(maker_data_matrix, 'w')

maker_data_matrix_handle.write('%s\t%s\n' % ('MAG', '\t'.join(marker_set_sorted)))
for each_bin in marker_num_dict_of_dict:

    marker_num_list = []
    for each_marker in marker_set_sorted:
        marker_num = marker_num_dict_of_dict[each_bin].get(each_marker, 0)
        marker_num_list.append(marker_num)

    maker_data_matrix_handle.write('%s\t%s\n' % (each_bin, '\t'.join([str(i) for i in marker_num_list])))

maker_data_matrix_handle.close()