import os
import glob
from Bio import SeqIO

########################################################################################################################
################################################### get psiblast cmds ##################################################
########################################################################################################################

# faa_id_txt          = '/Users/songweizhi/Desktop/Nitrosopumilaceae_MAGs_328_300_Prokka_faa_file_id.txt'
# cmd_file            = '/Users/songweizhi/Desktop/Nitrosopumilaceae_MAGs_328_300_psiblast_cmds.txt'
#
# cmd_file_handle = open(cmd_file, 'w')
# for each_faa in open(faa_id_txt):
#     gnm_id = each_faa.strip()
#     psiblast_op_txt          = '%s.psiblast.txt'         % gnm_id
#     psiblast_op_txt_best_hit = '%s.psiblast.BestHit.txt' % gnm_id
#     psiblast_cmd             = 'psiblast -query ../S1_Data_70_marker.fasta -subject %s.faa -out %s -outfmt 6 -evalue 1e-05 -num_threads 1' % (gnm_id, psiblast_op_txt)
#     get_best_hit_cmd         = 'BioSAK BestHit -i %s -o %s' % (psiblast_op_txt, psiblast_op_txt_best_hit)
#     cmds_in_one_line         = '%s; %s' % (psiblast_cmd, get_best_hit_cmd)
#     cmd_file_handle.write(cmds_in_one_line + '\n')
# cmd_file_handle.close()

########################################################################################################################
################################################ parse psiblast results ################################################
########################################################################################################################

##### input/output 1 #####
# combined_faa_file              = '/Users/songweizhi/Desktop/Yang_90_genomes_renamed_Prokka.faa'
# combined_psiblast_BestHit_file = '/Users/songweizhi/Desktop/Yang_90_genomes_renamed_Prokka_faa_psiblast_BestHit.txt'
# gnms_to_ignore_list            = ["GCF_000250675.2.yang", "GCF_000092245.1.yang", "GCF_000021565.1.yang", "GCF_000019965.1.yang", "GCF_000019625.1.yang", "GCF_000092565.1.yang", "GCF_000317475.1.yang"]
# marker_cov_cutoff              = 80
# output_dir                     = '/Users/songweizhi/Desktop/Marker_cov%s'            % marker_cov_cutoff

##### input/output #####
combined_faa_file              = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_Yang_70/combined_390_genomes/combined_390_genomes.faa'
combined_psiblast_BestHit_file = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_Yang_70/combined_390_genomes/combined_390_genomes_psiblast_BestHit.txt'
gnms_to_ignore_txt             = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_Yang_70/combined_390_genomes/gnms_to_ignore_77.txt'

marker_cov_cutoff              = 80
output_dir                     = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_Yang_70/combined_390_genomes/Marker_323_genomes_with_bacteria_cov%s' % marker_cov_cutoff

########################################################################################################################

os.mkdir(output_dir)

gnms_to_ignore_list = []
for each_gnm in open(gnms_to_ignore_txt):
    gnms_to_ignore_list.append(each_gnm.strip())

# read in faa file
faa_seq_dict = {}
for each_seq in SeqIO.parse(combined_faa_file, 'fasta'):
    faa_seq_dict[each_seq.id] = str(each_seq.seq)

gnm_set = set()
marker_to_gene_dict = {}
for each_hit in open(combined_psiblast_BestHit_file):
    each_hit_split = each_hit.strip().split('\t')
    query_id = each_hit_split[0]
    subject_id = each_hit_split[1]
    subject_gnm = '_'.join(subject_id.split('_')[:-1])
    if subject_gnm not in gnms_to_ignore_list:
        gnm_set.add(subject_gnm)
        if query_id not in marker_to_gene_dict:
            marker_to_gene_dict[query_id] = {subject_id}
        else:
            marker_to_gene_dict[query_id].add(subject_id)

print(gnm_set)
for each_marker in marker_to_gene_dict:
    current_marker_gene_set = marker_to_gene_dict[each_marker]
    marker_cov = len(current_marker_gene_set)*100/len(gnm_set)
    marker_cov = float("{0:.2f}".format(marker_cov))
    if marker_cov < marker_cov_cutoff:
        print('%s\t%s/%s\t%s' % (each_marker, len(current_marker_gene_set), len(gnm_set), marker_cov))
    else:
        current_marker_fasta = '%s/%s.fa' % (output_dir, each_marker)
        current_marker_fasta_handle = open(current_marker_fasta, 'w')
        for each_gene in current_marker_gene_set:
            current_marker_fasta_handle.write('>%s\n' % each_gene)
            current_marker_fasta_handle.write('%s\n' % faa_seq_dict[each_gene])
        current_marker_fasta_handle.close()

