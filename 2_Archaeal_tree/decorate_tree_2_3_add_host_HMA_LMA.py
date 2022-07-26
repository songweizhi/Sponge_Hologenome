from db_files import Sponge_HMA_LMA_status_dict

########################################################################################################################

wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5'

# file in
mag_name_file               = '%s/combined_Nitrosopumilaceae_MAGs_328_300_MAG_file_name.txt'    % wd
mag_host_cate_file          = '%s/Archaeal_tree_host_1_s.txt'                                   % wd

# file out

########################################################################################################################

for each in open(mag_host_cate_file):
    each_split = each.strip().split('\t')
    gnm_id = each_split[0]
    host_cate = each_split[1]

    if host_cate == 'nonsponge':
        host_HMA_LMA_status = 'nonsponge'
    elif host_cate == 'NA':
        host_HMA_LMA_status = 'NA'
    else:
        host_HMA_LMA_status = Sponge_HMA_LMA_status_dict.get(host_cate, 'unknown')
        print('%s\t%s' % (host_cate, host_HMA_LMA_status))



