import os
import glob


########################################################################################################################

# file/parameter in
gnm_folder          = '/Users/songweizhi/Desktop/DB'
prokka_op_folder    = '/Users/songweizhi/Desktop/DB_Prokka_op'
gnm_ext             = 'fa'
gnm_domain          = 'Archaea'
num_thread          = 12

# output file
cmd_file = '/Users/songweizhi/Desktop/DB_prokka_cmd.txt'

########################################################################################################################

# file/parameter in
gnm_folder          = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5/combined_Nitrosopumilaceae_MAGs_328_300'
prokka_op_folder    = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5/combined_Nitrosopumilaceae_MAGs_328_300_Prokka_op'
gnm_ext             = 'fna'
gnm_domain          = 'Archaea'
num_thread          = 12

# output file
cmd_file            = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5/combined_Nitrosopumilaceae_MAGs_328_300_Prokka_cmd.txt'

########################################################################################################################

# file/parameter in
gnm_folder          = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5_4_Marker_Yang_70/Yang_90_genomes'
prokka_op_folder    = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5_4_Marker_Yang_70/Yang_90_genomes_Prokka_op'
gnm_ext             = 'fna'
gnm_domain          = 'Archaea'
num_thread          = 12

# output file
cmd_file            = '/srv/scratch/z5039045/Sponge_hologenome/2_Archaeal_tree_50_5_4_Marker_Yang_70/Yang_90_genomes_Prokka_cmd.txt'

########################################################################################################################

gnm_file_re = '%s/*.%s' % (gnm_folder, gnm_ext)
gnm_file_list = [os.path.basename(file_name) for file_name in glob.glob(gnm_file_re)]

cmd_file_handle = open(cmd_file, 'w')
for each_gnm in gnm_file_list:
    gnm_no_ext = each_gnm[: -(len(gnm_ext) + 1)]
    pwd_gnm_file = '%s/%s' % (gnm_folder, each_gnm)
    current_gnm_op_folder = '%s/%s_prokka_op' % (prokka_op_folder, gnm_no_ext)
    cmd_prokka = 'prokka --compliant --cpus %s --kingdom %s --prefix %s --locustag %s --strain %s --outdir %s %s' % (num_thread, gnm_domain, gnm_no_ext, gnm_no_ext, gnm_no_ext, current_gnm_op_folder, pwd_gnm_file)
    cmd_file_handle.write(cmd_prokka + '\n')
cmd_file_handle.close()

