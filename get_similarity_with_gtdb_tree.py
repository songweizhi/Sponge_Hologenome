import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def boxplot_with_dots(num_lol, name_list, output_plot):

    data = [np.array(i) for i in num_lol]

    box_plot = plt.boxplot(data, labels=name_list, patch_artist=True,
                           whiskerprops=dict(color='lightblue', linewidth=2), capprops=dict(color='lightblue'))
    # set the color pf box
    for box in box_plot['boxes']:
        box.set(linewidth=0)
        box.set_facecolor('lightblue')

    # add dots, https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.plot.html
    col_index = 1
    for num_arrary in data:
        plt.plot(np.random.normal(col_index, 0.02, len(num_arrary)), num_arrary, '.', alpha=0.8, color='orange', markersize=6, markeredgewidth=0)
        col_index += 1

    # write out
    plt.tight_layout()
    plt.savefig(output_plot, bbox_inches='tight', dpi=300)
    plt.close()


marker_txt                    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/5_Archaeal_tree_50_5_Markers/Matrix_distance_Willis_subset_4_marker_id_48.tsv'
similarity_with_gtdb_tree_txt = '/Users/songweizhi/Desktop/Matrix_similarity_reformatted.txt'
output_plot                   = '/Users/songweizhi/Desktop/similarity_with_gtdb_tree_boxplot_with_dots.png'


# read in similarities
similarity_dict = {}
for each_value in open(similarity_with_gtdb_tree_txt):
    each_value_split = each_value.strip().split('\t')
    similarity_dict[each_value_split[0]] = float(each_value_split[1])


# get marker list
marker_list = []
for each_marker in open(marker_txt):
    marker_list.append(each_marker.strip())


select_marker_similarity_list = []
unselect_marker_similarity_list = []
for each in similarity_dict:
    if each in marker_list:
        select_marker_similarity_list.append(similarity_dict[each])
    else:
        unselect_marker_similarity_list.append(similarity_dict[each])

print('selected HOGs:\t%s +/- %s' % (np.mean(select_marker_similarity_list), np.std(select_marker_similarity_list)))
print('unselected HOGs:\t%s +/- %s' % (np.mean(unselect_marker_similarity_list), np.std(unselect_marker_similarity_list)))

# plot
boxplot_with_dots([select_marker_similarity_list, unselect_marker_similarity_list], ['selected HOGs', 'unselected HOGs'], output_plot)
