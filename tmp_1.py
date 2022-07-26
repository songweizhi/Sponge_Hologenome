from db_files import sponge_MAG_GTDB_archaea
from db_files import sponge_MAG_GTDB_bacteria
from db_files import Sponge_MAG_id_80_5_1171_txt

aaa = 'aaa.txt'

archaeal_MAG_set = set()
for each in open(sponge_MAG_GTDB_archaea):
    if not each.startswith('user_genome'):
        archaeal_MAG_set.add(each.strip().split('\t')[0])

bacterial_MAG_set = set()
for each in open(sponge_MAG_GTDB_bacteria):
    if not each.startswith('user_genome'):
        bacterial_MAG_set.add(each.strip().split('\t')[0])

aaa_handle = open(aaa, 'w')
archaeal_MAG_num = 0
bacterial_MAG_num = 0
for each_mag in open(Sponge_MAG_id_80_5_1171_txt):
    mag_id = each_mag.strip()[:-4]
    if (mag_id in archaeal_MAG_set) or ((mag_id + 's') in archaeal_MAG_set):
        archaeal_MAG_num += 1
    elif (mag_id in bacterial_MAG_set) or ((mag_id + 's') in bacterial_MAG_set):
        bacterial_MAG_num += 1
        cp_md = 'cp Sponge_MAGs_80_5_1171/%s Sponge_MAGs_80_5_1171_bacteria_1134/\n' % each_mag.strip()
        print(cp_md)
        aaa_handle.write(cp_md)
    else:
        print(mag_id)

aaa_handle.close()

print('%s + %s = %s' % (bacterial_MAG_num, archaeal_MAG_num, (archaeal_MAG_num + bacterial_MAG_num)))