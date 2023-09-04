from db_files import GTDB_r202_ar122_marker_set
from db_files import GTDB_r207_ar53_marker_set
from db_files import hmm_profile_metadata_txt


marker_table_txt = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5_Markers/Marker_table.txt'


combined_r202_r207_marker_set = GTDB_r202_ar122_marker_set.union(GTDB_r207_ar53_marker_set)


hmm_no_version_to_label_dict = {}
hmm_no_version_to_product_name_dict = {}
for each_hmm in open(hmm_profile_metadata_txt):
    each_hmm_split = each_hmm.strip().split('\t')
    hmm_source_identifier = each_hmm_split[1]
    hmm_label             = each_hmm_split[2]
    hmm_product_name      = each_hmm_split[10]
    hmm_source_identifier_no_version = hmm_source_identifier
    if '.' in hmm_source_identifier:
        hmm_source_identifier_no_version = hmm_source_identifier.split('.')[0]
    hmm_no_version_to_label_dict[hmm_source_identifier_no_version] = hmm_label
    hmm_no_version_to_product_name_dict[hmm_source_identifier_no_version] = hmm_product_name


marker_table_txt_handle = open(marker_table_txt, 'w')
marker_table_txt_handle.write('Marker\tr202\tr207\tProduct\n')
for each_marker in combined_r202_r207_marker_set:

    each_marker_no_version = each_marker
    if '.' in each_marker:
        each_marker_no_version = each_marker.split('.')[0]

    current_hmm_label        = hmm_no_version_to_label_dict.get(each_marker_no_version, 'NA')
    current_hmm_product_name = hmm_no_version_to_product_name_dict.get(each_marker_no_version, 'NA')

    in_r202 = 0
    if each_marker in GTDB_r202_ar122_marker_set:
        in_r202 = 1

    in_r207 = 0
    if each_marker in GTDB_r207_ar53_marker_set:
        in_r207 = 1

    marker_table_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_marker, in_r202, in_r207, current_hmm_product_name))

marker_table_txt_handle.close()





