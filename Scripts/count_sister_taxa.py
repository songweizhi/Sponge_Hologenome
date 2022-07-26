from __future__ import print_function
import sys, re, numpy
from collections import defaultdict
from operator import itemgetter
from ete3 import Tree

# usage: python count_sister_taxa.py MLFile bootstrapFile outputFile

'''
Notes (Weizhi Song)

1. ete3 treats the input tree as rooted
2. leaf.name contain taxonomy information, format: to be added

cd /Users/songweizhi/Desktop/iqtree2_wd

# get the tree and bootstrap files (takes about XXX mins)
/Users/songweizhi/Software/iqtree/iqtree-2.2.0-MacOSX/bin/iqtree2 --wbtl --bnni -st AA -m LG+C60+F -T 10 -B 1000 -alrt 1000 -s example2.phy -pre example2_iqtree

# run count_sister_taxa.py
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree.contree example2_iqtree.ufboot outputFile_example1
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree2.contree example2_iqtree2.ufboot outputFile_example2wrong
python ~/PycharmProjects/Sponge_Hologenome/Scripts/count_sister_taxa.py example2_iqtree3.contree example2_iqtree3.ufboot outputFile_example3wrong

'''


def map_species_to_cluster(cluster_file):  # make a dict that links species name to the cluster to use for group compariosns
    spname_to_cluster = {}
    inh = open(cluster_file)
    for line in inh:
        if line.startswith("Names"):
            continue
        elements = re.split("\t", line)
        print(elements[1])
        e2 = re.split("\|", elements[1])

        species = e2[-2]
        cluster = e2[0]
        spname_to_cluster[species] = cluster
    return spname_to_cluster


def parse_taxonomy(taxon_name):  # given a taxon name, try to return whatever taxonomic info is available as a list starting with the highest level classification and going lower (or a map?)
    #name_elements = re.split("\|", taxon_name)
    name_elements = taxon_name.split('|')

    #print(name_elements)
    if (len(name_elements) < 8) or (len(name_elements) > 9):
        print(name_elements)
        print("Nonstandard!")
        quit()
    name_map = {}
    name_map['cluster'] = name_elements[0]
    name_map['domain'] = name_elements[1]
    name_map['phylum'] = name_elements[2]
    name_map['class'] = name_elements[3]
    name_map['order'] = name_elements[4]
    name_map['family'] = name_elements[5]
    name_map['genus'] = name_elements[6]
    name_map['species'] = name_elements[7]
    if len(name_elements) == 9:
        name_map['ncbi_id'] = name_elements[8]
    return name_map


def summarize_taxonomy(name_list, tax_level):  # take a list of names from a clade and summarize taxonomic info (labels and their frequencies)
    total_size = len(name_list)  # it perhaps makes sense to normalize by the size of the clade
    breakdown = {}
    for name in name_list:
        info = name_to_tax_info[name]
        if info[tax_level] in breakdown:
            breakdown[info[tax_level]] += 1.0 / float(total_size)
        else:
            breakdown[info[tax_level]] = 1.0 / float(total_size)
    return breakdown


# edit this to make the comparisons at a desired taxonomic level
target_label = 'cluster'
# compute the most frequent sister group of each (monophyletic?) group on the tree, to identify trends in gene transfers, "unstable" taxa, etc.

# read the ML tree, set up the taxonomy stuff, and calculate the number of clades per label, and the sizes of those clades (to report at the end)
labels = {}
name_to_tax_info = defaultdict(dict)
taxa_names = []
ml_tree = Tree(sys.argv[1])  # note that ete3 treats this input tree as rooted
for leaf in ml_tree:
    taxonomy = parse_taxonomy(leaf.name)
    name_to_tax_info[leaf.name] = taxonomy
    taxa_names.append(leaf.name)
    leaf.add_feature("tax", taxonomy[target_label])
    labels[taxonomy[target_label]] = 1
groups = labels.keys()

# print('labels')
# print(labels)
# print()
# print('groups')
# print(groups)
# print()
# print('name_to_tax_info')
# print(name_to_tax_info)
# print()
# print('taxa_names')
# print(taxa_names)
# print()

# compute the number of clades per label in the ML tree, and their sizes
ML_groups = defaultdict(list)  # the list is the size of each clade, len(list) is the number of clades for that label in the ML tree
for label in groups:
    for node in ml_tree.get_monophyletic(values=[label], target_attr="tax"):
        size_clade = 0
        for leaf in node:
            size_clade += 1
        ML_groups[label].append(size_clade)

summary = defaultdict(dict)
clades_per_group = defaultdict(list)
treeNum = -1
tree_sample_handle = open(sys.argv[2])
for line in tree_sample_handle:
    treeNum += 1
    tree = Tree(line.rstrip())
    for leaf in tree:


        tax = name_to_tax_info[leaf.name]  # this should set up taxonomy correctly...
        leaf.add_feature("tax", tax[target_label])  # this adds a feature called tax to the leaf, with the attribute of the phylum name
    for label in groups:
        clades_per_group[label].append(0.0)  # setup the clade counting for this particular tree
    tree.unroot()

    # iterate over groups that are monophyletic for the taxon label of choice. Choose the smallest sister branch for the comparison. (Assume root is within the larger sister clade)
    for label in groups:
        # print label
        for node in tree.get_monophyletic(values=[label], target_attr="tax"):
            clades_per_group[label][treeNum] += 1.0
            # print node.get_ascii()
            sisters = node.get_sisters()
            if node.is_root():
                continue
            elif len(sisters) == 1:  # not at the trifurcation. Do something a bit hacky to find the bigger sister clade
                taxa_in_OG = []
                taxa_in_sister = []
                taxa_in_group = []
                for leaf in sisters[0]:
                    taxa_in_sister.append(leaf.name)
                for leaf in node:
                    taxa_in_group.append(leaf.name)
                size_sister = len(taxa_in_sister)
                for leaf_name in taxa_names:
                    if leaf_name in taxa_in_sister:
                        continue
                    elif leaf_name in taxa_in_group:
                        continue
                    else:
                        taxa_in_OG.append(leaf_name)
                size_OG = len(taxa_in_OG)
                sister_tax = {}
                if size_OG > size_sister:
                    sister_tax = summarize_taxonomy(taxa_in_sister, target_label)
                else:
                    sister_tax = summarize_taxonomy(taxa_in_OG, target_label)
                # store the tax info of the sister group
                for element in sister_tax:
                    if element in summary[label]:
                        summary[label][element] += sister_tax[element]
                    else:
                        summary[label][element] = sister_tax[element]
            else:  # trifurcation in tree. Just treat the two sisters in the same way.
                taxa_in_sisters0 = []
                taxa_in_sisters1 = []

                for leaf in sisters[0]:
                    taxa_in_sisters0.append(leaf.name)
                size_s0 = len(taxa_in_sisters0)
                for leaf in sisters[1]:
                    taxa_in_sisters1.append(leaf.name)
                size_s1 = len(taxa_in_sisters1)

                sister_tax = {}
                if size_s0 > size_s1:
                    sister_tax = summarize_taxonomy(taxa_in_sisters1, target_label)
                else:
                    sister_tax = summarize_taxonomy(taxa_in_sisters0, target_label)

                for element in sister_tax:
                    if element in summary[label]:
                        summary[label][element] += sister_tax[element]
                    else:
                        summary[label][element] = sister_tax[element]

# now print out some kind of summary. For each label, the sorted list of sister taxa and their frequencies?

outh = open(sys.argv[3], "w")

print('summary')
print(summary)

for label in summary:
    num_groups = len(ML_groups[label])
    size_str = ''
    if num_groups == 1:
        size_str = ML_groups[label][0]
    else:
        size_str = ','.join(str(x) for x in (sorted(ML_groups[label], reverse=True)))
    avg_num_clades = numpy.mean(clades_per_group[label])
    total_num_clades = numpy.sum(clades_per_group[label])
    sorted_sisters = sorted(summary[label].items(), key=itemgetter(1), reverse=True)
    #print('sorted_sisters')
    #print(sorted_sisters)
    for tup in sorted_sisters:
        double_normalize = float(tup[1]) / float(total_num_clades)  # normalize the frequencies by the total number of clades, to account for different bootstrap numbers/MCMC sample numbers
        outh.write(label + "\t" + tup[0] + "\t" + str(tup[1]) + "\t" + str(avg_num_clades) + "\t" + str(double_normalize) + "\t" + str(size_str) + "\n")
outh.close()

print('Done!')
