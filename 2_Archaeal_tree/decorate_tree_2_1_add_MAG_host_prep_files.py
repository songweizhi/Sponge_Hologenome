from db_files import sponge_MAG_host_dict
from db_files import manually_added_gnm_to_host_low_dict
from db_files import gtdb_archaeal_gnm_biosample_dict

description = '''

This script was wrote to extract the source/host of the GTDB/NCBI genomes in a tree.
You'll need to manually execute the commands in esearch_cmd_txt and summarise the results.
BioSAK exe_cmds -c esearch_cmd.txt -t 8 (not working)

Input files: 
mag_name_file

'''

########################################################################################################################

#wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree'
wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/4_Archaeal_tree_50_5'

# file in
#mag_name_file   = '%s/combined_archaeal_MAGs_234_209_MAG_file_name.txt' % wd
mag_name_file   = '%s/combined_Nitrosopumilaceae_MAGs_328_300_MAG_file_name.txt' % wd

# file out
biosample_txt   = '%s/genome_to_biosample_id.txt'   % wd
esearch_cmd_txt = '%s/esearch_cmd.txt'              % wd

########################################################################################################################

biosample_txt_handle = open(biosample_txt, 'w')
esearch_cmd_txt_handle = open(esearch_cmd_txt, 'w')
for each_MAG in open(mag_name_file):

    # get mag_id, mag_id_no_ext, mag_id_no_ext_no_source (gtdb or ncbi) and mag_id_no_ext_no_source_no_suffix (.1s, .2s, .3s, ...)
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
    if mag_id_no_ext_no_source in sponge_MAG_host_dict:
        mag_host = sponge_MAG_host_dict[mag_id_no_ext_no_source]
    if mag_id_no_ext_no_source in manually_added_gnm_to_host_low_dict:
        mag_host = manually_added_gnm_to_host_low_dict[mag_id_no_ext_no_source]

    if mag_host == 'NA':
        mag_biosample = gtdb_archaeal_gnm_biosample_dict.get(mag_id_no_ext_no_source_no_suffix, 'NA')
        if mag_biosample == 'NA':
            mag_id_no_ext_no_source_no_suffix_GCF = mag_id_no_ext_no_source_no_suffix.replace('GCA', 'GCF')
            if mag_id_no_ext_no_source_no_suffix_GCF in gtdb_archaeal_gnm_biosample_dict:
                mag_biosample = gtdb_archaeal_gnm_biosample_dict[mag_id_no_ext_no_source_no_suffix_GCF]

        esearch_cmd = 'esearch -db biosample -query %s | esummary | xtract -pattern DocumentSummary -element Accession -first Title -group Attribute -if Attribute@harmonized_name -equals "host" -element Attribute -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute > %s_metadata.txt' % (mag_biosample, mag_biosample)
        biosample_txt_handle.write('%s\t%s\n' % (mag_id_no_ext_no_source_no_suffix, mag_biosample))
        esearch_cmd_txt_handle.write(esearch_cmd + '\n')
biosample_txt_handle.close()
esearch_cmd_txt_handle.close()


print('cd %s' % wd)
print('BioSAK exe_cmds -c esearch_cmd.txt -t 8')