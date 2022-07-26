
id_file             = '/Users/songweizhi/Desktop/gapseq/pwy_ids.txt'
combined_db_file    = '/Users/songweizhi/Desktop/gapseq/combined_pwy.tbl'

interested_pwys_txt = '/Users/songweizhi/Desktop/gapseq/interested_pwys.txt'


pwy_id_set = set()
for each in open(id_file):
    pwy_id_set.add(each.strip())
print(pwy_id_set)

interested_pwys_txt_handle = open(interested_pwys_txt, 'w')

for each_pwy in open(combined_db_file):

    pwy_id = each_pwy.split('\t')[0][1:-1]
    if pwy_id in pwy_id_set:
        print(pwy_id)
        interested_pwys_txt_handle.write(each_pwy)

interested_pwys_txt_handle.close()