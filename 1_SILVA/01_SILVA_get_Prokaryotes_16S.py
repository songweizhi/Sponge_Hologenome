import matplotlib.pyplot as plt


# file in
silva_SSURef_metadata = '/Users/songweizhi/Desktop/Sponges/1_SILVA/SILVA_138.1_SSURef.full_metadata.txt'
silva_SSUParc_metadata = '/Users/songweizhi/Desktop/Sponges/1_SILVA/SILVA_138.1_SSUParc.full_metadata.txt'
# silva_metadata = '/Users/songweizhi/Desktop/Sponges/SILVA/SILVA_138.1_SSURef.full_metadata.head10000.txt'
# silva_metadata = '/Users/songweizhi/Desktop/Sponges/SILVA/SILVA_138.1_SSUParc.full_metadata.head10000.txt'


# file out
silva_SSUParc_metadata_sponge_short = '/Users/songweizhi/Desktop/Sponges/1_SILVA/SILVA_138.1_SSUParc.full_metadata_sponge_subset_300-1199bp.txt'
silva_SSUParc_metadata_sponge_long  = '/Users/songweizhi/Desktop/Sponges/1_SILVA/SILVA_138.1_SSUParc.full_metadata_sponge_subset_1200bp.txt'


silva_SSUParc_metadata_sponge_short_handle = open(silva_SSUParc_metadata_sponge_short, 'w')
silva_SSUParc_metadata_sponge_long_handle = open(silva_SSUParc_metadata_sponge_long, 'w')
seq_id_list = []
seq_len_list = []
header_index_dict = {}
for each in open(silva_SSUParc_metadata):
    each_split = each.strip().split('\t')
    if each.startswith('acc	start	stop'):
        silva_SSUParc_metadata_sponge_short_handle.write(each)
        silva_SSUParc_metadata_sponge_long_handle.write(each)
        header_index_dict = {k: v for v, k in enumerate(each_split)}
    else:
        isolation_source = each_split[header_index_dict['isolation_source']]
        if 'sponge' in isolation_source:
            length = int(each_split[2]) - int(each_split[1]) + 1
            seq_id_list.append(each_split[0])
            seq_len_list.append(length)

            if length >= 1200:
                silva_SSUParc_metadata_sponge_long_handle.write(each)
            else:
                silva_SSUParc_metadata_sponge_short_handle.write(each)

silva_SSUParc_metadata_sponge_short_handle.close()
silva_SSUParc_metadata_sponge_long_handle.close()


print('seq_id_list %s' % len(seq_id_list))

# plt.hist(seq_len_list, bins=300)
# plt.show()



'''
cd /Users/songweizhi/Desktop/Sponges/SILVA
BioSAK select_seq -seq SILVA_138.1_SSUParc_tax_silva.fasta -id SILVA_138.1_SSUParc.full_metadata_sponge_subset_300-1199bp_id.txt -out SILVA_138.1_SSUParc.full_metadata_sponge_subset_300-1199bp.fa -option 1 -oneline 
BioSAK select_seq -seq SILVA_138.1_SSUParc_tax_silva.fasta -id SILVA_138.1_SSUParc.full_metadata_sponge_subset_1200bp_id.txt -out SILVA_138.1_SSUParc.full_metadata_sponge_subset_1200bp.fa -option 1 -oneline

'''

