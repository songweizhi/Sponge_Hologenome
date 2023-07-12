import os


########################################################################################################################

#################### r207 ####################

# OrthologousGroups_txt                   = '/Users/songweizhi/Desktop/OrthologousGroups.txt'
# pwd_OrthologousGroupsFasta              = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5/Nitrosopumilaceae_OMA_wd/Output/OrthologousGroupsFasta'
# pwd_OrthologousGroupsFasta_qualified    = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5/Nitrosopumilaceae_OMA_wd/Output/OrthologousGroupsFasta_cov_85'
# gnm_num                                 = 300
# cov_cutoff                              = 85
# cov_cutoff    90  85  80
# OMA_groups    9   118 307

#################### r214 ####################

# file in
OrthologousGroups_txt                   = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_OMA_wd_r214/Output/OrthologousGroups.txt'
pwd_OrthologousGroupsFasta              = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_OMA_wd_r214/Output/OrthologousGroupsFasta'
gnm_num                                 = 375
cov_cutoff                              = 90

# file out
pwd_OrthologousGroupsFasta_qualified    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_OMA_wd_r214/Output/OrthologousGroupsFasta_min%s_cov%s' % (gnm_num, cov_cutoff)

# cov_cutoff    90      85      80
# OMA_groups    12     132     310

########################################################################################################################

if os.path.isdir(pwd_OrthologousGroupsFasta_qualified) is True:
    os.system('rm -r %s' % pwd_OrthologousGroupsFasta_qualified)
os.system('mkdir %s' % pwd_OrthologousGroupsFasta_qualified)


qualified_grp_num = 0
for each_group in open(OrthologousGroups_txt):
    each_group_split = each_group.strip().split('\t')
    group_id = each_group_split[0]

    group_id_only_num = group_id.replace('OMA', '')
    while group_id_only_num[0] == '0':
        group_id_only_num = group_id_only_num[1:]

    group_member_list = each_group_split[1:]
    group_member_list_no_description = []
    group_member_list_no_description_colon = []
    group_member_gnm_set = set()
    for each_protein in group_member_list:
        protein_no_description = each_protein.split(' ')[0]
        protein_no_description_colon = each_protein.split(' ')[0].split(':')[1]
        group_member_gnm = '_'.join(each_protein.split(' ')[0].split('_')[:-1])
        group_member_gnm_set.add(group_member_gnm)
        group_member_list_no_description.append(protein_no_description)
        group_member_list_no_description_colon.append(protein_no_description_colon)

    cov = len(group_member_gnm_set)*100/gnm_num

    if cov >= cov_cutoff:
        qualified_grp_num += 1
        print('%s\t%s' % (group_id, group_member_list_no_description_colon))
        fasta_file_id = 'OG%s' % group_id_only_num
        pwd_fasta_file = '%s/OG%s.fa' % (pwd_OrthologousGroupsFasta, group_id_only_num)
        cp_cmd = 'cp %s %s/' % (pwd_fasta_file, pwd_OrthologousGroupsFasta_qualified)
        os.system(cp_cmd)

        muscle_cmd = 'muscle -in %s.fa -out %s.aln' % (fasta_file_id, fasta_file_id)
        trimal_cmd = 'trimal -in %s.aln -out %s_trimmed.phy -gappyout -keepheader -phylip' % (fasta_file_id, fasta_file_id)
        #print(muscle_cmd)
        #print(trimal_cmd)

        iqtree2_cmd = 'mkdir %s_iqtree_with_UFBoot; cd %s_iqtree_with_UFBoot; /srv/scratch/z5039045/Softwares/iqtree-2.2.0-Linux/bin/iqtree2 --wbtl --bnni -st AA -m LG+C60+F -T 16 -B 1000 -alrt 1000 -s ../%s_trimmed.phy -pre %s_iqtree' % (fasta_file_id, fasta_file_id, fasta_file_id, fasta_file_id)
        #print(iqtree2_cmd)

print('Orthologous groups with coverage >= %s%s: %s' % (cov_cutoff, '%', qualified_grp_num))
