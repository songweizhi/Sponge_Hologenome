from db_files import sponge_MAG_cpl_dict
from db_files import sponge_MAG_ctm_dict
from db_files import sponge_archaeal_MAG_tax_dict

print(len(sponge_archaeal_MAG_tax_dict))

for each in sponge_archaeal_MAG_tax_dict:

    mag_id_no_suffix = each.strip()
    if mag_id_no_suffix[-3:] in ['.1s', '.2s', '.3s', '.4s', '.5s']:
        mag_id_no_suffix = mag_id_no_suffix[:-1]
    print('%s\t%s\t%s\t%s' % (each, sponge_MAG_cpl_dict.get(mag_id_no_suffix, 'na'), sponge_MAG_ctm_dict.get(mag_id_no_suffix, 'na'), sponge_archaeal_MAG_tax_dict[each]))


