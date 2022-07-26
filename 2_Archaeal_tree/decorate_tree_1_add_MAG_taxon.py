import os

from db_files import sponge_archaeal_MAG_tax_dict
from db_files import gtdb_ar_gnm_tax_dict
from common_functions import gtdb_tax_str_parser

description = '''

This script extracts the taxon of genomes in a tree.

'''

########################################################################################################################

#wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_Archaeal_tree'
wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5'

# file in
#mag_name_file       = '%s/combined_archaeal_MAGs_234_209_MAG_file_name.txt' % wd
mag_name_file       = '%s/combined_Nitrosopumilaceae_MAGs_328_300_MAG_file_name.txt' % wd
taxon_rank          = 'g'

# file out
mag_taxon_s         = '%s/Archaeal_tree_genome_taxon_s.txt'                 % wd
mag_taxon_g         = '%s/Archaeal_tree_genome_taxon_g.txt'                 % wd
mag_taxon_g_itol    = '%s/Archaeal_tree_genome_taxon_g_iTOL_ColorRange.txt' % wd
mag_taxon_s_itol    = '%s/Archaeal_tree_genome_taxon_s_iTOL_Labels.txt'     % wd

########################################################################################################################

mag_taxon_s_itol_handle = open(mag_taxon_s_itol, 'w')
mag_taxon_s_itol_handle.write('LABELS\nSEPARATOR TAB\nDATA\n\n')
mag_taxon_s_handle = open(mag_taxon_s, 'w')
mag_taxon_g_handle = open(mag_taxon_g, 'w')
for each_MAG in open(mag_name_file):

    # get mag_id, mag_id_no_ext and mag_id_no_ext_no_source (gtdb or ncbi)
    mag_id = each_MAG.strip()
    mag_id_no_ext = '.'.join(mag_id.split('.')[:-1])
    mag_id_no_ext_no_source = mag_id_no_ext
    if '.gtdb' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]
    if '.ncbi' in mag_id_no_ext_no_source:
        mag_id_no_ext_no_source = mag_id_no_ext_no_source[:-5]

    # get mag_taxon_str
    mag_taxon_str = 'NA'
    if mag_id_no_ext_no_source in sponge_archaeal_MAG_tax_dict:
        mag_taxon_str = sponge_archaeal_MAG_tax_dict[mag_id_no_ext_no_source]
    if mag_id_no_ext_no_source in gtdb_ar_gnm_tax_dict:
        mag_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source]

    # get mag_taxon_str (GCA GCF things)
    if mag_taxon_str == 'NA':
        mag_id_no_ext_no_source_GCF = mag_id_no_ext_no_source.replace('GCA', 'GCF')
        if mag_id_no_ext_no_source_GCF in gtdb_ar_gnm_tax_dict:
            mag_taxon_str = gtdb_ar_gnm_tax_dict[mag_id_no_ext_no_source_GCF]

    # get needed_taxon
    needed_taxon = 'NA'
    mag_species = 'NA'
    if mag_taxon_str != 'NA':
        needed_taxon = gtdb_tax_str_parser(mag_taxon_str, taxon_rank)
        mag_species = gtdb_tax_str_parser(mag_taxon_str, 's')

    # write out
    mag_taxon_s_handle.write('%s\t%s\n' % (mag_id_no_ext, mag_species))
    mag_taxon_g_handle.write('%s\t%s\n' % (mag_id_no_ext, needed_taxon))
    mag_taxon_s_itol_handle.write('%s\t%s__%s\n' % (mag_id_no_ext, mag_species, mag_id_no_ext))
mag_taxon_s_handle.close()
mag_taxon_g_handle.close()
mag_taxon_s_itol_handle.close()

itol_cmd = 'BioSAK iTOL -ColorRange -lg %s -lt Genus -out %s' % (mag_taxon_g, mag_taxon_g_itol)
print(itol_cmd)
os.system(itol_cmd)

'''
/usr/local/bin/python
/Library/Frameworks/Python.framework/Versions/3.7/bin/python3
'''
