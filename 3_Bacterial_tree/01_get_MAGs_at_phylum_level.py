from db_files import sponge_bacterial_MAG_tax_dict


interested__p = 'p__Cyanobacteria'


interested__p_to_gnm_dict = {}
taxon_to_gnm_dict = {}
for each_gnm in open('/Users/songweizhi/Documents/Research/Sponge_Hologenome/2_MAGs_tree/Bacteria_tree_r207/Sponge_MAGs_80_5_1171_bacteria_1134_genomes.txt'):
    gnm_file = each_gnm.strip()
    gnm_id = gnm_file[:-4]
    gnm_taxon = sponge_bacterial_MAG_tax_dict[gnm_id]
    gnm_taxon_split = gnm_taxon.split(';')
    gnm__p = gnm_taxon_split[1]
    gnm__c = gnm_taxon_split[2]

    if gnm__p not in taxon_to_gnm_dict:
        taxon_to_gnm_dict[gnm__p] = [gnm_file]
    else:
        taxon_to_gnm_dict[gnm__p].append(gnm_file)

    if gnm__p == interested__p:
        if gnm__c not in interested__p_to_gnm_dict:
            interested__p_to_gnm_dict[gnm__c] = [gnm_file]
        else:
            interested__p_to_gnm_dict[gnm__c].append(gnm_file)


for each_tax in taxon_to_gnm_dict:
    print('%s\t%s' % (each_tax, len(taxon_to_gnm_dict[each_tax])))
print()


print('Interested taxon: %s (%s)' % (interested__p, len(taxon_to_gnm_dict.get(interested__p, []))))
for each_tax in interested__p_to_gnm_dict:
    print('%s\t%s' % (each_tax, len(interested__p_to_gnm_dict[each_tax])))
print()

