import os
import glob
from Bio import SeqIO


mag_folder                  = '/Users/songweizhi/Desktop/Sponges/MAGs/Engelberts_2020/prokaryotes_2022-03-19_genomes'
mag_metadata                = '/Users/songweizhi/Desktop/Sponges/MAGs/Engelberts_2020/Engelberts_2020_MAGs_259.txt'
mag_folder_renamed_selected = '/Users/songweizhi/Desktop/Sponges/MAGs/Engelberts_2020/prokaryotes_2022-03-19_genomes_renamed_selected'
mag_folder_renamed_not_used = '/Users/songweizhi/Desktop/Sponges/MAGs/Engelberts_2020/prokaryotes_2022-03-19_genomes_renamed_not_used'
mag_ext                     = 'fna'


mag_re = '%s/*.%s' % (mag_folder, mag_ext)
mag_list = [os.path.basename(file_name) for file_name in glob.glob(mag_re)]

selected_mag_list = []
for each_line in open(mag_metadata):
    if not each_line.startswith('Accession'):
        each_line_split = each_line.strip().split('\t')
        col2_split = each_line_split[1].split('_')
        bin_id = '%s_%s' % (col2_split[-4], col2_split[-3])
        selected_mag_list.append(bin_id)

found_mag_list = []
for each_mag in mag_list:
    pwd_mag = '%s/%s' % (mag_folder, each_mag)
    current_mag_index = ''
    for each_seq in SeqIO.parse(pwd_mag, 'fasta'):
        mag_des_split = each_seq.description.split(' ')
        for each_str in mag_des_split:
            if '_bin_' in each_str:
                current_mag_index = each_str

    current_mag_index = current_mag_index.replace('_bin_', '_bin')
    pwd_mag_renamed = '%s/%s_%s' % (mag_folder_renamed_not_used, current_mag_index, each_mag)
    if current_mag_index in selected_mag_list:
        pwd_mag_renamed = '%s/%s_%s' % (mag_folder_renamed_selected, current_mag_index, each_mag)
        found_mag_list.append(current_mag_index)
    os.system('cp %s %s' % (pwd_mag, pwd_mag_renamed))

print(found_mag_list)
print(selected_mag_list)
