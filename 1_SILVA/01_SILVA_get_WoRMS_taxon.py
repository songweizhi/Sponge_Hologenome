import os


# file in
WoRMS_taxlist_txt                   = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/3_Sponge_diversity_and_evolution/combined_WoRMS_taxlist.txt'
SSURef_NR99_Sponge_18S_tax          = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_SILVA/SILVA_Porifera/SILVA_138.1_SSURef_NR99_tax_silva_Sponge_18S_taxon.txt'
LSURef_NR99_Sponge_28S_tax          = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_SILVA/SILVA_Porifera/SILVA_138.1_LSURef_NR99_tax_silva_Sponge_28S_taxon.txt'

# file out
SSURef_NR99_Sponge_18S_WoRMS_tax    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_SILVA/SILVA_Porifera/SILVA_138.1_SSURef_NR99_tax_silva_Sponge_18S_WoRMS_taxon.txt'
LSURef_NR99_Sponge_28S_WoRMS_tax    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_SILVA/SILVA_Porifera/SILVA_138.1_LSURef_NR99_tax_silva_Sponge_28S_WoRMS_taxon.txt'


WoRMS_taxon_str_set = set()
for each_line in open(WoRMS_taxlist_txt):
    each_line_split = each_line.strip().split('\t')
    if not each_line.startswith('AphiaID'):
        sponge__p = each_line_split[7]
        sponge__c = each_line_split[8]
        sponge__o = each_line_split[9]
        sponge__f = each_line_split[10]
        sponge__g = each_line_split[12]
        sponge__sg = each_line_split[13]
        sponge__s = each_line_split[14]
        sponge__ss = each_line_split[15]
        sponge__sn = each_line_split[4]
        WoRMS_taxon_str = 'p__%s;c__%s;o__%s;f__%s;g__%s;sg__%s;s__%s;ss__%s;sn__%s' % (sponge__p, sponge__c, sponge__o, sponge__f, sponge__g, sponge__sg, sponge__s, sponge__ss, sponge__sn)
        WoRMS_taxon_str_set.add(WoRMS_taxon_str)


n = 0
LSURef_NR99_Sponge_28S_WoRMS_tax_handle = open(LSURef_NR99_Sponge_28S_WoRMS_tax, 'w')
for each_line in open(LSURef_NR99_Sponge_28S_tax):
    each_line_split = each_line.strip().split('\t')
    seq_desc_split = each_line_split[1].split(';')

    # this is good
    genus_name = seq_desc_split[-1].split(' ')[0]
    if genus_name == 'Candidatus':
        genus_name = seq_desc_split[-1].split(' ')[1]

    matched_taxon_set = set()
    for each_taxon in WoRMS_taxon_str_set:
        if (genus_name) in each_taxon:
            matched_taxon_set.add(each_taxon)

    # get matched WoRMS taxon
    matched_WoRMS_taxon = ''
    if len(matched_taxon_set) == 0:
        matched_WoRMS_taxon = 'NA'
    elif len(matched_taxon_set) == 1:
        matched_WoRMS_taxon = list(matched_taxon_set)[0]
    else:
        # get lowest common taxon
        common_prefix = os.path.commonprefix(list(matched_taxon_set))
        while common_prefix.endswith('__'):
            common_prefix = ';'.join(common_prefix.split(';')[:-1])
        common_prefix_lowest_rank = common_prefix.split(';')[-1]
        if common_prefix_lowest_rank not in list(matched_taxon_set)[0].split(';'):
            common_prefix = ';'.join(common_prefix.split(';')[:-1])
        matched_WoRMS_taxon = common_prefix

    seq_id_formatted = each_line.strip().replace('\t', ' ').replace(' ', '_').replace('.', '_').replace(';', '_')
    seq_id_formatted_polished = seq_id_formatted + (80 - len(seq_id_formatted))*'_'

    LSURef_NR99_Sponge_28S_WoRMS_tax_handle.write('%s\t%s_WoRMS;%s\n' % (seq_id_formatted, seq_id_formatted_polished, matched_WoRMS_taxon))
LSURef_NR99_Sponge_28S_WoRMS_tax_handle.close()
