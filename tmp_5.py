
hog_top_25_pct_txt      = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_25perc/HOGs_best_25perc.txt'
hog_top_50_pct_txt      = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers_by_split_best_50perc/HOGs_best_50perc.txt'
hog_48_by_cluster_txt   = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/Matrix_distance_Willis_subset_4_marker_id_48.tsv'


hog_top_25_pct_set = set()
for each in open(hog_top_25_pct_txt):
    hog_top_25_pct_set.add(each.strip())

hog_top_50_pct_set = set()
for each in open(hog_top_50_pct_txt):
    hog_top_50_pct_set.add(each.strip())

hog_48_by_cluster_set = set()
for each in open(hog_48_by_cluster_txt):
    hog_48_by_cluster_set.add(each.strip())


c = set(hog_top_25_pct_set).intersection(hog_48_by_cluster_set)
d = set(hog_top_50_pct_set).intersection(hog_48_by_cluster_set)

print('hog_top_25_pct_set:\t%s'                             % len(hog_top_25_pct_set))
print('hog_top_50_pct_set:\t%s'                             % len(hog_top_50_pct_set))
print('hog_48_by_cluster_set:\t%s'                          % len(hog_48_by_cluster_set))
print('hog_48_by_cluster_set and hog_top_25_pct_set:\t%s'   % len(c))
print('hog_48_by_cluster_set and hog_top_50_pct_set:\t%s'   % len(d))


