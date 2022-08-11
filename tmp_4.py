from Bio import SeqIO


################################################### get psiblast cmds ##################################################

# faa_folder = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5_1/combined_Nitrosopumilaceae_MAGs_328_300_Prokka_faa'
# cmd_txt    = '/Users/songweizhi/Desktop/psiblast_cmds.txt'

# cmd_txt_handle = open(cmd_txt, 'w')
# for each_gnm in open('/Users/songweizhi/Desktop/gnm_id.txt'):
#     gnm_id = each_gnm.strip()
#     psiblast_cmd = 'psiblast -query S1_Data_70_marker.fasta -subject %s/%s.faa -out %s_psiblast.txt -outfmt 6 -evalue 1e-05 -num_threads 1' % (faa_folder, gnm_id, gnm_id)
#     print(psiblast_cmd)
#     cmd_txt_handle.write(psiblast_cmd + '\n')
# cmd_txt_handle.close()

################################################ parse psiblast results ################################################

# read in faa file
combined_faa_file = '/Users/songweizhi/Desktop/combined_Nitrosopumilaceae_MAGs_328_300_Prokka.faa'
faa_seq_dict = {}
for each_seq in SeqIO.parse(combined_faa_file, 'fasta'):
    faa_seq_dict[each_seq.id] = str(each_seq.seq)


marker_to_gene_dict = {}
for each_hit in open('/Users/songweizhi/Desktop/psiblast_op_BestHit.txt'):
    each_hit_split = each_hit.strip().split('\t')
    query_id = each_hit_split[0]
    subject_id = each_hit_split[1]
    if query_id not in marker_to_gene_dict:
        marker_to_gene_dict[query_id] = {subject_id}
    else:
        marker_to_gene_dict[query_id].add(subject_id)


for each_marker in marker_to_gene_dict:
    current_marker_gene_set = marker_to_gene_dict[each_marker]
    if len(current_marker_gene_set)/300 < 0.8:
        pass
        #print('ignore\t%s\t%s' % (each_marker, len(current_marker_gene_set)))
    else:
        current_marker_fasta = '/Users/songweizhi/Desktop/Yang_70_Markers/%s.fa' % each_marker
        current_marker_fasta_handle = open(current_marker_fasta, 'w')
        for each_gene in current_marker_gene_set:
            current_marker_fasta_handle.write('>%s\n' % each_gene)
            current_marker_fasta_handle.write('%s\n' % faa_seq_dict[each_gene])
        current_marker_fasta_handle.close()


