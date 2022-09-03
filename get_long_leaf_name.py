from db_files import sponge_archaeal_MAG_tax_dict
from db_files import gtdb_ar_gnm_tax_dict
from common_functions import gtdb_tax_str_parser
from db_files import sponge_MAG_host_dict
from db_files import manually_added_gnm_to_host_low_dict
from db_files import gtdb_ar_gnm_tax_dict
from db_files import sponge_archaeal_MAG_tax_dict
from db_files import sponge_g_to_full_lineage_dict

########################################################################################################################

# file in
wd                      = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5'
gnm_to_biosample_txt    = '%s/genome_to_biosample_id.txt'                                  % wd
biosample_source_txt    = '%s/combined_Biosample_metadata_reformatted.txt'                 % wd
leaf_id_txt             = '/Users/songweizhi/Desktop/leaf_id.txt'

# file out
leaf_label_iTOL_txt     = '/Users/songweizhi/Desktop/leaf_label_iTOL.txt'

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

mag_to_species_str_dict = {}
mag_to_taxon_str_dict = {}
mag_to_sponge_lineage_dict = {}
mag_max_len = 0
taxon_str_max_len = 0
species_max_len = 0
sponge_lineage_max_len = 0
for each_MAG in open(leaf_id_txt):

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
        needed_taxon = gtdb_tax_str_parser(mag_taxon_str, 'g')
        mag_species = gtdb_tax_str_parser(mag_taxon_str, 's')

    if len(mag_id_no_ext) > mag_max_len:
        mag_max_len = len(mag_id_no_ext)
    if len(mag_taxon_str) > taxon_str_max_len:
        taxon_str_max_len = len(mag_taxon_str)
    if len(mag_species) > species_max_len:
        species_max_len = len(mag_species)

    mag_to_taxon_str_dict[mag_id_no_ext] = mag_taxon_str
    mag_to_species_str_dict[mag_id_no_ext] = mag_species

    # get host info
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

    sponge_genus = 'NA'
    if mag_host_cate not in ['nonsponge', 'NA']:
        sponge_genus = mag_host_cate.split('_')[0]
    sponge_full_lineage = sponge_g_to_full_lineage_dict.get(('g__'+ sponge_genus), 'NA')
    mag_to_sponge_lineage_dict[mag_id_no_ext] = sponge_full_lineage
    if len(sponge_full_lineage) > sponge_lineage_max_len:
        sponge_lineage_max_len = len(sponge_full_lineage)


leaf_label_iTOL_txt_handle = open(leaf_label_iTOL_txt, 'w')
leaf_label_iTOL_txt_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for each_mag in mag_to_taxon_str_dict:
    current_taxon_str = mag_to_taxon_str_dict[each_mag]
    current_taxon_str_polished = current_taxon_str + (taxon_str_max_len - len(current_taxon_str))*'_'
    mag_id_str_polished = each_mag + (mag_max_len - len(each_mag))*'_'
    mag_species_str = mag_to_species_str_dict[each_mag]
    mag_species_str_polished = mag_species_str + (species_max_len - len(mag_species_str))*'_'
    sponge_lineage = mag_to_sponge_lineage_dict[each_mag]
    print(sponge_lineage)
    leaf_label_iTOL_txt_handle.write('%s\t%s__%s__%s\n' % (each_mag, mag_species_str_polished.replace(' ', '_'), mag_id_str_polished, sponge_lineage))
leaf_label_iTOL_txt_handle.close()

