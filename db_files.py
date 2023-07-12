
def gtdb_gnm_metadata_parser(gtdb_genome_metadata):

    genome_to_taxon_dict = {}
    genome_to_completeness_dict = {}
    genome_to_contamination_dict = {}
    genome_to_biosample_dict = {}

    col_index = {}
    for each_ref in open(gtdb_genome_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession'):
            col_index = {key: i for i, key in enumerate(each_ref_split)}
        else:
            ref_accession = each_ref_split[0][3:]
            gnm_completeness = float(each_ref_split[2])
            gnm_contamination = float(each_ref_split[3])
            gtdb_taxon = each_ref_split[col_index['gtdb_taxonomy']]
            ncbi_biosample = each_ref_split[col_index['ncbi_biosample']]
            genome_to_taxon_dict[ref_accession] = gtdb_taxon
            genome_to_completeness_dict[ref_accession] = gnm_completeness
            genome_to_contamination_dict[ref_accession] = gnm_contamination
            genome_to_biosample_dict[ref_accession] = ncbi_biosample

    return genome_to_completeness_dict, genome_to_contamination_dict, genome_to_taxon_dict, genome_to_biosample_dict


########################################################################################################################

metadata_dir    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac'
GTDB_dir_Katana = '/srv/scratch/z5039045/DB/GTDB_r207'
GTDB_dir_Mac    = '/Users/songweizhi/DB/GTDB_r207'

########################################################################################################################
################################################### GTDB files (r207) ##################################################
########################################################################################################################

genome_metadata_ar53_r207_Katana    = '%s/ar53_metadata_r207.tsv'               % GTDB_dir_Katana
genome_metadata_bac120_r207_Katana  = '%s/bac120_metadata_r207.tsv'             % GTDB_dir_Katana
genome_taxonomy_ar53_r207_Katana    = '%s/ar53_taxonomy.tsv'                    % GTDB_dir_Katana
genome_taxonomy_bac120_r207_Katana  = '%s/bac120_taxonomy.tsv'                  % GTDB_dir_Katana
GTDB_genome_paths_r207_Katana       = '%s/release207/fastani/genome_paths.tsv'  % GTDB_dir_Katana
genome_metadata_ar53_r207_Mac       = '%s/ar53_metadata_r207.tsv'               % GTDB_dir_Mac
genome_metadata_bac120_r207_Mac     = '%s/bac120_metadata_r207.tsv'             % GTDB_dir_Mac
genome_taxonomy_ar53_r207_Mac       = '%s/ar53_taxonomy.tsv'                    % GTDB_dir_Mac
genome_taxonomy_bac120_r207_Mac     = '%s/bac120_taxonomy.tsv'                  % GTDB_dir_Mac
GTDB_genome_paths_r207_Mac          = '%s/genome_paths.tsv'                     % GTDB_dir_Mac
GTDB_r202_ar122_markers_txt         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac/GTDB_r202_ar122_markers.txt'
GTDB_r207_ar53_markers_txt          = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac/GTDB_r207_ar53_markers.txt'
hmm_profile_metadata_txt            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac/hmm_PGAP.tsv'

gtdb_archaeal_gnm_cpl_dict, gtdb_archaeal_gnm_ctm_dict, gtdb_ar_gnm_tax_dict, gtdb_archaeal_gnm_biosample_dict   = gtdb_gnm_metadata_parser(genome_metadata_ar53_r207_Mac)
gtdb_bacterial_gnm_cpl_dict, gtdb_bacterial_gnm_ctm_dict, gtdb_bacterial_gnm_tax_dict, gtdb_bacterial_gnm_biosample_dict = gtdb_gnm_metadata_parser(genome_metadata_bac120_r207_Mac)

gtdb_ar_cpl_dict,  gtdb_ar_ctm_dict,  gtdb_ar_tax_dict,  gtdb_ar_biosample_dict  = gtdb_gnm_metadata_parser(genome_metadata_ar53_r207_Mac)
gtdb_bac_cpl_dict, gtdb_bac_ctm_dict, gtdb_bac_tax_dict, gtdb_bac_biosample_dict = gtdb_gnm_metadata_parser(genome_metadata_bac120_r207_Mac)

GTDB_r202_ar122_marker_set = set()
for each_r202_marker in open(GTDB_r202_ar122_markers_txt):
    GTDB_r202_ar122_marker_set.add(each_r202_marker.strip())

GTDB_r207_ar53_marker_set = set()
for each_r207_marker in open(GTDB_r207_ar53_markers_txt):
    GTDB_r207_ar53_marker_set.add(each_r207_marker.strip())

########################################################################################################################
###################################################### sponge MAG ######################################################
########################################################################################################################

sponge_MAG_host_txt                  = '%s/Sponge_MAGs_1677_host.txt'            % metadata_dir
manually_added_gnm_to_host_txt       = '%s/manually_added_genome_to_host.csv'    % metadata_dir
sponge_MAG_CheckM_txt                = '%s/Sponge_MAGs_1677_CheckM.txt'          % metadata_dir
sponge_MAG_GTDB_archaea              = '%s/Sponge_MAGs_1677.ar53.summary.tsv'    % metadata_dir
sponge_MAG_GTDB_bacteria             = '%s/Sponge_MAGs_1677.bac120.summary.tsv'  % metadata_dir
Archaeal_mags_renamed_for_prokka_txt = '%s/Archaeal_mags_renamed_for_prokka.txt' % metadata_dir
# Sponge_MAG_id_75_5_1299_txt         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_MAGs/Sponge_MAGs_75_5_1299_genomes.txt'
# Sponge_MAG_id_80_5_1171_txt         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_MAGs/Sponge_MAGs_80_5_1171_genomes.txt'
# Sponge_MAG_id_85_5_1059_txt         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_MAGs/Sponge_MAGs_85_5_1059_genomes.txt'

sponge_archaeal_MAG_tax_dict = {}
for each in open(sponge_MAG_GTDB_archaea):
    if not each.startswith('user_genome'):
        each_split = each.strip().split('\t')
        sponge_archaeal_MAG_tax_dict[each_split[0]] = each_split[1]

sponge_bacterial_MAG_tax_dict = {}
for each in open(sponge_MAG_GTDB_bacteria):
    if not each.startswith('user_genome'):
        each_split = each.strip().split('\t')
        sponge_bacterial_MAG_tax_dict[each_split[0]] = each_split[1]

sponge_MAG_cpl_dict = {}
sponge_MAG_ctm_dict = {}
for each in open(sponge_MAG_CheckM_txt):
    if not each.startswith('Genome,Completeness,Contamination'):
        each_split = each.strip().split(',')
        gnm_id = each_split[0]
        gnm_cpl = float(each_split[1])
        gnm_ctm = float(each_split[2])
        sponge_MAG_cpl_dict[gnm_id] = gnm_cpl
        sponge_MAG_ctm_dict[gnm_id] = gnm_ctm

sponge_MAG_host_dict = {}
for each in open(sponge_MAG_host_txt):
    each_split = each.strip().split('\t')
    sponge_MAG_host_dict[each_split[0]] = each_split[1]

# read in manually added genome to host information
manually_added_gnm_to_host_low_dict = {}
manually_added_gnm_to_host_high_dict = {}
for each_gnm in open(manually_added_gnm_to_host_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    id_mag = each_gnm_split[0]
    host_low = each_gnm_split[1]
    host_high = each_gnm_split[2]
    manually_added_gnm_to_host_low_dict[id_mag] = host_low
    manually_added_gnm_to_host_high_dict[id_mag] = host_high

Archaeal_mags_renamed_for_prokka_dict_raw2new = {}
Archaeal_mags_renamed_for_prokka_dict_new2raw = {}
for each_renamed_gnm in open(Archaeal_mags_renamed_for_prokka_txt):
    each_renamed_gnm_split = each_renamed_gnm.strip().split('\t')
    name_raw = each_renamed_gnm_split[0]
    name_new = each_renamed_gnm_split[1]
    Archaeal_mags_renamed_for_prokka_dict_raw2new[name_raw] = name_new
    Archaeal_mags_renamed_for_prokka_dict_new2raw[name_new] = name_raw

########################################################################################################################
######################################################## sponge ########################################################
########################################################################################################################

sponge_full_taxon_txt       = '%s/sponge_full_lineage_GTDB_format.txt'  % metadata_dir
Sponge_HMA_LMA_status_txt   = '%s/Sponge_HMA_LMA_status.txt'            % metadata_dir

# read in Sponge_HMA_LMA_status_txt
Sponge_HMA_LMA_status_dict = {}
for each_sponge in open(Sponge_HMA_LMA_status_txt):
    each_sponge_split = each_sponge.strip().split('\t')
    Sponge_HMA_LMA_status_dict[each_sponge_split[0]] = each_sponge_split[1]

# read in full taxonomy of sponges
sponge_g_to_c_dict = {}
sponge_g_to_sc_dict = {}
sponge_g_to_o_dict = {}
sponge_g_to_f_dict = {}
sponge_g_to_full_lineage_dict = {}
for each in open(sponge_full_taxon_txt):
    each_split = each.strip().split(';')
    sponge_c = ''
    sponge_sc = ''
    sponge_o = ''
    sponge_f = ''
    sponge_g = ''
    for each_level in each_split:
        if each_level.startswith('c__'):
            sponge_c = each_level
        if each_level.startswith('sc__'):
            sponge_sc = each_level
        if each_level.startswith('o__'):
            sponge_o = each_level
        if each_level.startswith('f__'):
            sponge_f = each_level
        if each_level.startswith('g__'):
            sponge_g = each_level
    sponge_g_to_c_dict[sponge_g] = sponge_c
    sponge_g_to_sc_dict[sponge_g] = sponge_sc
    sponge_g_to_o_dict[sponge_g] = sponge_o
    sponge_g_to_f_dict[sponge_g] = sponge_f
    sponge_g_to_full_lineage_dict[sponge_g] = each.strip()

########################################################################################################################
########################################################################################################################
########################################################################################################################
