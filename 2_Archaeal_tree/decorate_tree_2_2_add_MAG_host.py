from db_files import sponge_MAG_host_dict
from db_files import manually_added_gnm_to_host_low_dict
from db_files import sponge_g_to_c_dict
from db_files import sponge_g_to_sc_dict
from db_files import sponge_g_to_o_dict
from db_files import sponge_g_to_f_dict
from db_files import gtdb_ar_gnm_tax_dict
from db_files import sponge_archaeal_MAG_tax_dict

description = '''  '''

########################################################################################################################

wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5'

# file in
#mag_name_file               = '%s/combined_archaeal_MAGs_234_209_MAG_file_name.txt'       % wd
#biosample_txt               = '%s/biosample_id.txt'                                       % wd
#biosample_source_txt        = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_MAGs_BioProjects/biosample_to_source.txt'

mag_name_file                = '%s/combined_Nitrosopumilaceae_MAGs_328_300_MAG_file_name.txt'   % wd
gnm_to_biosample_txt         = '%s/genome_to_biosample_id.txt'                                  % wd
biosample_source_txt         = '%s/combined_Biosample_metadata_reformatted.txt'                 % wd

# file out
mag_host_file               = '%s/Archaeal_tree_host_0_source.txt'                              % wd
mag_host_file_itol          = '%s/Archaeal_tree_host_0_source_iTOL.txt'                         % wd
mag_host_cate_file          = '%s/Archaeal_tree_host_1_s.txt'                                   % wd
mag_host_cate_file_itol     = '%s/Archaeal_tree_host_1_s_iTOL.txt'                              % wd
mag_host_txt_subset_g       = '%s/Archaeal_tree_host_2_g.txt'                                   % wd
mag_host_txt_subset_g_iTOL  = '%s/Archaeal_tree_host_2_g_iTOL.txt'                              % wd
mag_host_txt_subset_f       = '%s/Archaeal_tree_host_3_f.txt'                                   % wd
mag_host_txt_subset_f_iTOL  = '%s/Archaeal_tree_host_3_f_iTOL.txt'                              % wd
mag_host_txt_subset_o       = '%s/Archaeal_tree_host_4_o.txt'                                   % wd
mag_host_txt_subset_o_iTOL  = '%s/Archaeal_tree_host_4_o_iTOL.txt'                              % wd
mag_host_txt_subset_sc      = '%s/Archaeal_tree_host_5_sc.txt'                                  % wd
mag_host_txt_subset_sc_iTOL = '%s/Archaeal_tree_host_5_sc_iTOL.txt'                             % wd
mag_host_txt_subset_c       = '%s/Archaeal_tree_host_6_c.txt'                                   % wd
mag_host_txt_subset_c_iTOL  = '%s/Archaeal_tree_host_6_c_iTOL.txt'                              % wd
mag_host_txt_subset_scofg   = '%s/Archaeal_tree_taxon_and_host_scofg.txt'                       % wd

########################################################################################################################

# read in biosample_txt
assembly_id_to_biosample_dict = {}
for each_gnm in open(gnm_to_biosample_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    assembly_id_to_biosample_dict[each_gnm_split[0]] = each_gnm_split[1]

# read in biosample source
biosample_source_dict = {}
biosample_source_cate_dict = {}
for each_biosample in open(biosample_source_txt):
    each_biosample_split = each_biosample.strip().split('\t')
    biosample_source_dict[each_biosample_split[0]] = each_biosample_split[2]
    biosample_source_cate_dict[each_biosample_split[0]] = each_biosample_split[1]

mag_host_file_handle = open(mag_host_file, 'w')
mag_host_cate_file_handle = open(mag_host_cate_file, 'w')
for each_MAG in open(mag_name_file):

    # get mag_id, mag_id_no_ext and mag_id_no_ext_no_source (gtdb or ncbi)
    mag_id = each_MAG.strip()
    mag_id_no_ext = '.'.join(mag_id.split('.')[:-1])
    mag_id_no_ext_no_source = mag_id_no_ext
    if '.gtdb' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]
    if '.ncbi' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]

    mag_id_no_ext_no_source_no_suffix = mag_id_no_ext_no_source
    if mag_id_no_ext_no_source_no_suffix[-3:] in ['.1s', '.2s', '.3s', '.4s', '.5s']:
        mag_id_no_ext_no_source_no_suffix = mag_id_no_ext_no_source_no_suffix[:-1]

    mag_host = 'NA'
    mag_host_cate = 'NA'
    if mag_id_no_ext_no_source_no_suffix in sponge_MAG_host_dict:
        mag_host = sponge_MAG_host_dict[mag_id_no_ext_no_source_no_suffix]
        mag_host_cate = mag_host
    elif mag_id_no_ext_no_source_no_suffix in manually_added_gnm_to_host_low_dict:
        mag_host = manually_added_gnm_to_host_low_dict[mag_id_no_ext_no_source_no_suffix]
        mag_host_cate = mag_host
    else:
        gnm_biosample = assembly_id_to_biosample_dict[mag_id_no_ext_no_source_no_suffix]
        mag_host      = biosample_source_dict[gnm_biosample]
        mag_host_cate = biosample_source_cate_dict[gnm_biosample]

    # write out
    mag_host_file_handle.write('%s\t%s\n' % (mag_id_no_ext, mag_host))
    mag_host_cate_file_handle.write('%s\t%s\n' % (mag_id_no_ext, mag_host_cate))
mag_host_file_handle.close()
mag_host_cate_file_handle.close()


# get MAG host at higher taxon ranks
mag_host_txt_subset_g_handle = open(mag_host_txt_subset_g, 'w')
mag_host_txt_subset_f_handle = open(mag_host_txt_subset_f, 'w')
mag_host_txt_subset_o_handle = open(mag_host_txt_subset_o, 'w')
mag_host_txt_subset_sc_handle = open(mag_host_txt_subset_sc, 'w')
mag_host_txt_subset_c_handle = open(mag_host_txt_subset_c, 'w')
mag_host_txt_subset_scofg_handle = open(mag_host_txt_subset_scofg, 'w')
for each in open(mag_host_cate_file):
    each_split = each.strip().split('\t')
    gnm_id = each_split[0]

    # get mag_id, mag_id_no_ext and mag_id_no_ext_no_source (gtdb or ncbi)
    mag_id_no_ext_no_source = gnm_id
    if '.gtdb' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]
    if '.ncbi' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]

    mag_id_no_ext_no_source_no_suffix = mag_id_no_ext_no_source
    if mag_id_no_ext_no_source_no_suffix[-3:] in ['.1s', '.2s', '.3s', '.4s', '.5s']:
        mag_id_no_ext_no_source_no_suffix = mag_id_no_ext_no_source_no_suffix[:-1]

    # get mag_taxon_str
    mag_taxon_str = 'NA'
    if mag_id_no_ext_no_source_no_suffix in sponge_archaeal_MAG_tax_dict:
        mag_taxon_str = sponge_archaeal_MAG_tax_dict[mag_id_no_ext_no_source_no_suffix]
    elif mag_id_no_ext_no_source_no_suffix in gtdb_ar_gnm_tax_dict:
        mag_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source_no_suffix]

    sponge_species = each_split[1]
    if sponge_species == 'nonsponge':
        sponge_g  = 'nonsponge'
        sponge_f  = 'nonsponge'
        sponge_o  = 'nonsponge'
        sponge_sc = 'nonsponge'
        sponge_c  = 'nonsponge'
    elif sponge_species == 'NA':
        sponge_g  = 'NA'
        sponge_f  = 'NA'
        sponge_o  = 'NA'
        sponge_sc = 'NA'
        sponge_c  = 'NA'
    else:
        sponge_g = 'g__' + sponge_species.split('_')[0]
        if '(coral)' in sponge_species:
            sponge_g = 'g__%s(coral)' % sponge_species.split('_')[0]
        sponge_f  = sponge_g_to_f_dict.get(sponge_g, 'NA')
        sponge_o  = sponge_g_to_o_dict.get(sponge_g, 'NA')
        sponge_sc = sponge_g_to_sc_dict.get(sponge_g, 'NA')
        sponge_c  = sponge_g_to_c_dict.get(sponge_g, 'NA')

    mag_host_txt_subset_g_handle.write('%s\t%s\n' % (gnm_id, sponge_g))
    mag_host_txt_subset_f_handle.write('%s\t%s\n' % (gnm_id, sponge_f))
    mag_host_txt_subset_o_handle.write('%s\t%s\n' % (gnm_id, sponge_o))
    mag_host_txt_subset_sc_handle.write('%s\t%s\n' % (gnm_id, sponge_sc))
    mag_host_txt_subset_c_handle.write('%s\t%s\n' % (gnm_id, sponge_c))

    if sponge_sc not in ['NA', 'nonsponge']:
        mag_host_txt_subset_scofg_handle.write('%s\t%s\t%s;%s;%s;%s;s__%s\n' % (gnm_id, mag_taxon_str, sponge_sc, sponge_o, sponge_f, sponge_g, sponge_species))

mag_host_txt_subset_g_handle.close()
mag_host_txt_subset_f_handle.close()
mag_host_txt_subset_o_handle.close()
mag_host_txt_subset_sc_handle.close()
mag_host_txt_subset_c_handle.close()
mag_host_txt_subset_scofg_handle.close()


'''
cd /Users/songweizhi/Documents/Research/Sponge_Hologenome/2_Archaeal_tree_50_5
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_1_s.txt -lt Sponge_Species -gc group_color.txt -out Archaeal_tree_host_1_s_iTOL.txt
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_2_g.txt -lt Sponge_Genus -gc group_color.txt -out Archaeal_tree_host_2_g_iTOL.txt
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_3_f.txt -lt Sponge_Family -gc group_color.txt -out Archaeal_tree_host_3_f_iTOL.txt
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_4_o.txt -lt Sponge_Order -gc group_color.txt -out Archaeal_tree_host_4_o_iTOL.txt
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_5_sc.txt -lt Sponge_Subclass -gc group_color.txt -out Archaeal_tree_host_5_sc_iTOL.txt
BioSAK iTOL -ColorStrip -lg Archaeal_tree_host_6_c.txt -lt Sponge_Subclass -gc group_color.txt -out Archaeal_tree_host_6_c_iTOL.txt
'''
