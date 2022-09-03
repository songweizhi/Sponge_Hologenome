import os
from Bio import SeqIO
from time import sleep
from Bio import AlignIO


def concatenated_msa(msa_file_list, msa_file_dir, concatenated_msa_file, concatenated_msa_loci_file):

    # get genome id list
    gnm_id_set = set()
    for each_msa_file in msa_file_list:
        pwd_each_msa_file = '%s/%s' % (msa_file_dir, each_msa_file)
        for each_seq in AlignIO.read(pwd_each_msa_file, "fasta"):
            gnm_id = '_'.join(each_seq.id.split('_')[:-1])
            gnm_id_set.add(gnm_id)
    gnm_id_list_sorted = sorted([i for i in gnm_id_set])

    # initialize concatenated_msa_dict
    concatenated_msa_dict = {}
    for each_gnm in gnm_id_set:
        concatenated_msa_dict[each_gnm] = ''

    # concatenate msa
    concatenated_msa_loci_dict = {}
    concatenated_seq_len = 0
    for each_msa_file in msa_file_list:
        pwd_each_msa_file = '%s/%s' % (msa_file_dir, each_msa_file)

        # read in current msa
        current_marker_id_to_seq_dict = {}
        current_marker_aln_len_set = set()
        for each_marker in AlignIO.read(pwd_each_msa_file, "fasta"):
            current_marker_gnm = '_'.join(each_marker.id.split('_')[:-1])
            marker_aln_seq = str(each_marker.seq)
            current_marker_aln_len_set.add(len(marker_aln_seq))
            current_marker_id_to_seq_dict[current_marker_gnm] = marker_aln_seq

        # get length of current msa
        if len(current_marker_aln_len_set) > 1:
            print('MSA not with the same length, program exited!')
            exit()

        # add info to concatenated_msa_loci_dict
        current_marker_aln_len = [i for i in current_marker_aln_len_set][0]
        concatenated_msa_loci_dict[each_msa_file] = [str(concatenated_seq_len + 1), str(concatenated_seq_len + current_marker_aln_len)]
        concatenated_seq_len += current_marker_aln_len

        # add current msa to concatenated_msa_dict
        for each_gnm in gnm_id_list_sorted:
            current_gnm_aln = current_marker_id_to_seq_dict.get(each_gnm, current_marker_aln_len*'-')
            concatenated_msa_dict[each_gnm] += current_gnm_aln

    # write out concatenated msa
    concatenated_msa_file_handle = open(concatenated_msa_file, 'w')
    for each_gnm in concatenated_msa_dict:
        concatenated_msa_file_handle.write('>%s\n' % each_gnm)
        concatenated_msa_file_handle.write('%s\n'  % concatenated_msa_dict[each_gnm])
    concatenated_msa_file_handle.close()

    # write out concatenated msa loci
    concatenated_msa_loci_file_handle = open(concatenated_msa_loci_file, 'w')
    for each_msa in msa_file_list:
        marker_loci = concatenated_msa_loci_dict[each_msa]
        marker_loci_str = '\t'.join(marker_loci)
        concatenated_msa_loci_file_handle.write('%s\t%s\n' % (each_msa, marker_loci_str))
    concatenated_msa_loci_file_handle.close()


########################################################################################################################

# # file in
# marker_txt                    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/Matrix_distance_Willis_subset_2_marker_id_89.tsv'
# marker_seq_dir                = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/OrthologousGroupsFasta_cov_85_no_empty_line'
#
# # file out
# wd                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/marker_89'
# aln_dir                       = '%s/marker_89_alignment_dir'                          % wd
# concatenated_msa_file         = '%s/marker_89_concatenated.aln'                       % wd
# concatenated_msa_loci_file    = '%s/marker_89_concatenated_loci.txt'                  % wd
# concatenated_msa_file_trimmed = '%s/marker_89_concatenated_trimmed.aln'               % wd

########################################################################################################################

# # file in
# marker_txt                    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Yang_70_markers/62_markers.txt'
# marker_seq_dir                = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Yang_70_markers/Yang_70_Markers'
#
# # file out
# wd                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Yang_70_markers/62_markers'
# aln_dir                       = '%s/Yang_62_marker_alignment_dir'                          % wd
# concatenated_msa_file         = '%s/Yang_62_marker_concatenated.aln'                       % wd
# concatenated_msa_loci_file    = '%s/Yang_62_marker_concatenated_loci.txt'                  % wd
# concatenated_msa_file_trimmed = '%s/Yang_62_marker_concatenated_trimmed.aln'               % wd

########################################################################################################################

# # file in
# marker_txt                    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_50perc/HOGs_best_50perc.txt'
# marker_seq_dir                = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/OrthologousGroupsFasta_cov_85_no_empty_line'
# gnm_id_txt                    = ''
#
# # file out
# op_dir                        = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_50perc/best_50perc'
# marker_seq_dir_subset_gnm     = '%s/marker_best_50perc_marker_seq_dir_subset_gnm'   % op_dir
# aln_dir                       = '%s/marker_best_50perc_alignment_dir'               % op_dir
# concatenated_msa_file         = '%s/marker_best_50perc_concatenated.aln'            % op_dir
# concatenated_msa_loci_file    = '%s/marker_best_50perc_concatenated_loci.txt'       % op_dir
# concatenated_msa_file_trimmed = '%s/marker_best_50perc_concatenated_trimmed.aln'    % op_dir
#
# '''
# cd /Users/songweizhi/Desktop
# BioSAK rename_leaves -tree marker_by_split_best_25perc.contree -txt /Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka_for_renaming_back.txt -out marker_by_split_best_25perc.renamed2.contree
# '''

########################################################################################################################

# # file in
# marker_txt                    = '/Users/songweizhi/Desktop/combined_390_genomes/marker_id_63.txt'
# marker_seq_dir                = '/Users/songweizhi/Desktop/combined_390_genomes/Marker_390_genomes_with_bacteria_cov80'
# gnm_id_txt                    = '/Users/songweizhi/Desktop/combined_390_genomes/gnm_id_193_with_bacteria.txt'
#
# # file out
# op_dir                        = '/Users/songweizhi/Desktop/combined_390_genomes/MSA_gnm_193_maker63'
# marker_seq_dir_subset_gnm     = '%s/Yang_markers_seq_dir_subset_gnm'                                    % op_dir
# aln_dir                       = '%s/Yang_markers_alignment_dir'                                         % op_dir
# concatenated_msa_file         = '%s/Yang_markers_concatenated.aln'                                      % op_dir
# concatenated_msa_loci_file    = '%s/Yang_markers_concatenated_loci.txt'                                 % op_dir
# concatenated_msa_file_trimmed = '%s/Yang_markers_concatenated_trimmed.aln'                              % op_dir
#
# '''
# cd /Users/songweizhi/Desktop
# BioSAK rename_leaves -tree marker_by_split_best_25perc.contree -txt /Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka_for_renaming_back.txt -out marker_by_split_best_25perc.renamed2.contree
# '''

########################################################################################################################

# file in
wd                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Yang_70_markers/combined_390_genomes'
marker_txt                    = '%s/Marker_328_genomes_with_bacteria_cov79_marker_id_64.txt'    % wd
marker_seq_dir                = '%s/Marker_328_genomes_with_bacteria_cov79'                     % wd
gnm_id_txt                    = '%s/gnm_id_193_with_bacteria.txt'                               % wd
gnm_id_txt                    = '%s/gnm_id_193_with_bacteria_after_rename.txt'                  % wd
gnms_to_ignore_txt            = '%s/gnms_to_ignore.txt'                                         % wd

# file out
op_dir                        = '%s/MSA_gnm_120_maker64'                                        % wd
marker_seq_dir_subset_gnm     = '%s/Yang_markers_seq_dir_subset_gnm'                            % op_dir
aln_dir                       = '%s/Yang_markers_alignment_dir'                                 % op_dir
concatenated_msa_file         = '%s/Yang_markers_concatenated.aln'                              % op_dir
concatenated_msa_loci_file    = '%s/Yang_markers_concatenated_loci.txt'                         % op_dir
concatenated_msa_file_trimmed = '%s/Yang_markers_concatenated_trimmed.aln'                      % op_dir

########################################################################################################################

# # get gnms_to_ignore_list
# gnms_to_ignore_set = set()
# for each_gnm in open(gnms_to_ignore_txt):
#     gnms_to_ignore_set.add(each_gnm.strip())
#
# # get genome list
# gnm_list = []
# for each_gnm in open(gnm_id_txt):
#     if each_gnm.strip() not in gnms_to_ignore_set:
#         gnm_list.append(each_gnm.strip())
# print('gnm_list:\t%s' % len(gnm_list))
#
# from db_files import Archaeal_mags_renamed_for_prokka_dict_raw2new
# gnm_list_renamed = [Archaeal_mags_renamed_for_prokka_dict_raw2new.get(i, i) for i in gnm_list]
#
# # get marker list
# marker_list = []
# for each_marker in open(marker_txt):
#     marker_list.append(each_marker.strip())

list_1 = ['GCF_000250675.2.yang', 'GCF_000092245.1.yang', 'GCF_000021565.1.yang', 'GCF_000019965.1.yang', 'GCF_000019625.1.yang', 'GCF_000092565.1.yang', 'GCF_000317475.1.yang', 'APA_bin_56', 'GCA_000200715.1', 'CC_bin.11', 'COS4_bin_17', 'COS36387_bin_15', 'CYMC_67496', 'GCA_000024305.1.yang', 'GCA_000246735.1.yang', 'GCA_000698785.1.yang', 'GCA_000802205.2.yang', 'GCA_001917865.1.yang', 'GCA_002506665.1.gtdb', 'GCA_003056285.1.yang', 'GCA_003702485.1.yang', 'GCA_007570915.1.gtdb', 'GCA_007571135.1.gtdb', 'GCA_007571225.1.gtdb', 'GCA_011087695.1.ncbi', 'GCA_011088205.1.ncbi', 'GCA_012271085.1.ncbi', 'GCA_013867245.1.gtdb', 'GCA_013911135.1.gtdb', 'GCA_014075315.1.gtdb', 'GCA_016125975.1.gtdb', 'GCA_016126025.1s', 'GCA_021296165.1', 'GCA_021296175.1', 'GCA_900143675.1.yang', 'GCA_900177045.1.yang', 'GCA_900620265.1.gtdb', 'GCF_000007005.1.yang', 'GCF_000007065.1.yang', 'GCF_000007185.1.yang', 'GCF_000007225.1.yang', 'GCF_000008265.1.yang', 'GCF_000011185.1.yang', 'GCF_000011205.1.yang', 'GCF_000012285.1.yang', 'GCF_000015205.1.yang', 'GCF_000015225.1.yang', 'GCF_000015765.1.yang', 'GCF_000016385.1.yang', 'GCF_000016605.1.yang', 'GCF_000018305.1.yang', 'GCF_000018465.1.yang', 'GCF_000019605.1.yang', 'GCF_000020905.1.yang', 'GCF_000022205.1.yang', 'GCF_000022365.1.yang', 'GCF_000023945.1.yang', 'GCF_000025625.1.yang', 'GCF_000025665.1.yang', 'GCF_000026045.1.yang', 'GCF_000063445.1.yang', 'GCF_000145985.1.yang', 'GCF_000148385.1.yang', 'GCF_000151105.2.yang', 'GCF_000152265.2.yang', 'GCF_000166095.1.yang', 'GCF_000172995.2.yang', 'GCF_000194625.1.yang', 'GCF_000195915.1.yang', 'GCF_000195935.2.yang', 'GCF_000213215.1.yang', 'GCF_000214415.1.yang', 'GCF_000215995.1.yang', 'GCF_000220175.1.yang', 'GCF_000241145.1.yang', 'GCF_000242875.2.yang', 'GCF_000243315.1.yang', 'GCF_000246985.2.yang', 'GCF_000253055.1.yang', 'GCF_000270325.1.yang', 'GCF_000299365.1.yang', 'GCF_000299395.1.yang', 'GCF_000303155.1.yang', 'GCF_000376445.1.yang', 'GCF_000730285.1.yang', 'GCF_000812185.1.yang', 'GCF_000875775.1.yang', 'GCF_000956175.1.yang', 'GCF_002906215.1.yang', 'GCF_900065925.1.yang', 'GCF_900167955.1.yang', 'GCF_900248165.1.yang', 'GCF_900696045.1.yang', 'IMG_2263082000.yang', 'IMG_2264867067.yang', 'IMG_2264867070.yang', 'IMG_2264867229.yang', 'IMG_2513237066.yang', 'IMG_2513237068.yang', 'IMG_2524023104.yang', 'IMG_2527291500.yang', 'IMG_2527291509.yang', 'IMG_2545555825.yang', 'IMG_2558309099.yang', 'IMG_2619618950.yang', 'IMG_2654587960.yang', 'IMG_2654588083.yang', 'IMG_2708742552.yang', 'IMG_2718217642.yang', 'IMG_2721755844.yang', 'GCA_001543015.1', 'RHO1_bin_37', 'SB0662_bin33_GCA_009840065.1', 'SB0663_bin5_GCA_009839185.1', 'SB0664_bin35_GCA_009838485.1', 'SB0666_bin15_GCA_009837245.1', 'SB0667_bin13_GCA_009836545.1', 'SB0677_bin16_GCA_009842575.1', 'Shan_2019_B06_GCA_003724255.1', 'Shan_2019_D6_GCA_003724325.1', 'Shan_2019_H8_GCA_003724235.1', 'Shan_2019_H13_GCA_003724285.1', 'Shan_2019_S13_GCA_003724275.1', 'Shan_2019_S14_GCA_003724215.1', 'Shan_2019_S15_GCA_003724175.1', 'GCA_001541925.1']
list_1_renamed = ['GCF_000250675.2.yang', 'GCF_000092245.1.yang', 'GCF_000021565.1.yang', 'GCF_000019965.1.yang', 'GCF_000019625.1.yang', 'GCF_000092565.1.yang', 'GCF_000317475.1.yang', 'APA_bin_56', 'GCA_000200715.1', 'CC_bin.11', 'COS4_bin_17', 'COS36387_bin_15', 'CYMC_67496', 'GCA_000024305.1.yang', 'GCA_000246735.1.yang', 'GCA_000698785.1.yang', 'GCA_000802205.2.yang', 'GCA_001917865.1.yang', 'GCA_002506665.1.gtdb', 'GCA_003056285.1.yang', 'GCA_003702485.1.yang', 'GCA_007570915.1.gtdb', 'GCA_007571135.1.gtdb', 'GCA_007571225.1.gtdb', 'GCA_011087695.1.ncbi', 'GCA_011088205.1.ncbi', 'GCA_012271085.1.ncbi', 'GCA_013867245.1.gtdb', 'GCA_013911135.1.gtdb', 'GCA_014075315.1.gtdb', 'GCA_016125975.1.gtdb', 'GCA_016126025.1s', 'GCA_021296165.1', 'GCA_021296175.1', 'GCA_900143675.1.yang', 'GCA_900177045.1.yang', 'GCA_900620265.1.gtdb', 'GCF_000007005.1.yang', 'GCF_000007065.1.yang', 'GCF_000007185.1.yang', 'GCF_000007225.1.yang', 'GCF_000008265.1.yang', 'GCF_000011185.1.yang', 'GCF_000011205.1.yang', 'GCF_000012285.1.yang', 'GCF_000015205.1.yang', 'GCF_000015225.1.yang', 'GCF_000015765.1.yang', 'GCF_000016385.1.yang', 'GCF_000016605.1.yang', 'GCF_000018305.1.yang', 'GCF_000018465.1.yang', 'GCF_000019605.1.yang', 'GCF_000020905.1.yang', 'GCF_000022205.1.yang', 'GCF_000022365.1.yang', 'GCF_000023945.1.yang', 'GCF_000025625.1.yang', 'GCF_000025665.1.yang', 'GCF_000026045.1.yang', 'GCF_000063445.1.yang', 'GCF_000145985.1.yang', 'GCF_000148385.1.yang', 'GCF_000151105.2.yang', 'GCF_000152265.2.yang', 'GCF_000166095.1.yang', 'GCF_000172995.2.yang', 'GCF_000194625.1.yang', 'GCF_000195915.1.yang', 'GCF_000195935.2.yang', 'GCF_000213215.1.yang', 'GCF_000214415.1.yang', 'GCF_000215995.1.yang', 'GCF_000220175.1.yang', 'GCF_000241145.1.yang', 'GCF_000242875.2.yang', 'GCF_000243315.1.yang', 'GCF_000246985.2.yang', 'GCF_000253055.1.yang', 'GCF_000270325.1.yang', 'GCF_000299365.1.yang', 'GCF_000299395.1.yang', 'GCF_000303155.1.yang', 'GCF_000376445.1.yang', 'GCF_000730285.1.yang', 'GCF_000812185.1.yang', 'GCF_000875775.1.yang', 'GCF_000956175.1.yang', 'GCF_002906215.1.yang', 'GCF_900065925.1.yang', 'GCF_900167955.1.yang', 'GCF_900248165.1.yang', 'GCF_900696045.1.yang', 'IMG_2263082000.yang', 'IMG_2264867067.yang', 'IMG_2264867070.yang', 'IMG_2264867229.yang', 'IMG_2513237066.yang', 'IMG_2513237068.yang', 'IMG_2524023104.yang', 'IMG_2527291500.yang', 'IMG_2527291509.yang', 'IMG_2545555825.yang', 'IMG_2558309099.yang', 'IMG_2619618950.yang', 'IMG_2654587960.yang', 'IMG_2654588083.yang', 'IMG_2708742552.yang', 'IMG_2718217642.yang', 'IMG_2721755844.yang', 'GCA_001543015.1', 'RHO1_bin_37', 'GCA_001541925.1', 'GCA_009840065.1', 'GCA_009839185.1', 'GCA_009836545.1', 'GCA_003724325.1', 'GCA_003724285.1', 'GCA_003724275.1', 'GCA_003724175.1']


for each in list_1_renamed:
    if each not in  list_1:
        print(each)


print(len(list_1))
print(len(list_1_renamed))
