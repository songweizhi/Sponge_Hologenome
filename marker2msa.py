import os
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

# file in
marker_txt                    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_50perc/HOGs_best_50perc.txt'
marker_seq_dir                = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/OrthologousGroupsFasta_cov_85_no_empty_line'

# file out
wd                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_50perc/best_50perc'
aln_dir                       = '%s/marker_best_50perc_alignment_dir'                          % wd
concatenated_msa_file         = '%s/marker_best_50perc_concatenated.aln'                       % wd
concatenated_msa_loci_file    = '%s/marker_best_50perc_concatenated_loci.txt'                  % wd
concatenated_msa_file_trimmed = '%s/marker_best_50perc_concatenated_trimmed.aln'               % wd

'''
cd /Users/songweizhi/Desktop
BioSAK rename_leaves -tree marker_by_split_best_25perc.contree -txt /Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Archaeal_mags_renamed_for_prokka_for_renaming_back.txt -out marker_by_split_best_25perc.renamed2.contree
'''

########################################################################################################################

if os.path.isdir(wd):
    sleep(1)
    os.system('rm -r %s' % wd)

sleep(1)
os.system('mkdir %s' % wd)
os.system('mkdir %s' % aln_dir)

########################################################################################################################

# get marker list
marker_list = []
for each_marker in open(marker_txt):
    marker_list.append(each_marker.strip())

# align marker separately and trim them
trimmed_aln_file_list = []
for each_marker in marker_list:

    # define file name
    pwd_marker_seq         = '%s/%s.fa'          % (marker_seq_dir, each_marker)
    pwd_marker_aln         = '%s/%s.aln'         % (aln_dir, each_marker)
    marker_aln_trimmed     = '%s_trimmed.aln'    % each_marker
    pwd_marker_aln_trimmed = '%s/%s_trimmed.aln' % (aln_dir, each_marker)
    trimmed_aln_file_list.append(marker_aln_trimmed)

    # run muscle
    muscle_cmd = 'muscle -align %s -output %s' % (pwd_marker_seq, pwd_marker_aln)
    print(muscle_cmd)
    os.system(muscle_cmd)

    # run trimal
    trimal_cmd = 'trimal -in %s -out %s -gappyout -keepheader -fasta' % (pwd_marker_aln, pwd_marker_aln_trimmed)
    print(trimal_cmd)
    os.system(trimal_cmd)

# concatenate alignment
print('concatenating MSAs')
concatenated_msa(trimmed_aln_file_list, aln_dir, concatenated_msa_file, concatenated_msa_loci_file)

# trim concatenated alignment
concatenated_msa_trimal_cmd = 'trimal -in %s -out %s -gappyout -keepheader -fasta' % (concatenated_msa_file, concatenated_msa_file_trimmed)
print(concatenated_msa_trimal_cmd)
os.system(concatenated_msa_trimal_cmd)
