from Bio import SeqIO


fa_in       = '/Users/songweizhi/Desktop/assembly.fasta'
split_num   = 100
op_folder   = '/Users/songweizhi/Desktop/assembly_split'


# seq_num = 0
# for each_seq in SeqIO.parse(fa_in, 'fasta'):
#     seq_num += 1
# seq_num_per_subset = round(seq_num/split_num) + 1
#
# seq_index = 0
# for each_seq in SeqIO.parse(fa_in, 'fasta'):
#     seq_index += 1
#     op_file_index = seq_index//seq_num_per_subset + 1
#
#     pwd_subset_file = '%s/subset_%s.fa' % (op_folder, op_file_index)
#     with open(pwd_subset_file, 'a') as pwd_subset_file_handle:
#         pwd_subset_file_handle.write('>%s\n' % each_seq.description)
#         pwd_subset_file_handle.write('%s\n' % str(each_seq.seq))
#


for each_num in range(1, 101):

    cmd = 'blastn -db /data/bio/blastv5/nt -query subset_%s.fa -outfmt 6 -out /srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/assembly_split_tab/subset_%s.tab -num_threads 12' % (each_num, each_num)
    print(cmd)
