
done_ones_txt = '/Users/songweizhi/Desktop/done_ones.txt'
cmd_txt = '/Users/songweizhi/Desktop/OrthologousGroupsFasta_cov_85_iqtree2_cmds_with_UFBoot.txt'

done_list = []
for each in open(done_ones_txt):
    done_list.append(each.strip())



for each in open(cmd_txt):
    id = each.split('_iqtree_with_UFBoot; cd ')[0][6:]
    if id not in done_list:
        print(each.strip())
