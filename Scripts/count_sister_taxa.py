from __future__ import print_function
import argparse
import sys, re, numpy
from collections import defaultdict
from operator import itemgetter
from ete3 import Tree


# usage: python count_sister_taxa.py MLFile bootstrapFile outputFile

'''
Notes (Weizhi Song)

1. ete3 treats the input tree as rooted
2. leaf.name contain taxonomy information, format: to be added

'''


# def map_species_to_cluster(cluster_file):  # make a dict that links species name to the cluster to use for group compariosns
#     spname_to_cluster = {}
#     inh = open(cluster_file)
#     for line in inh:
#         if line.startswith("Names"):
#             continue
#         elements = re.split("\t", line)
#         print(elements[1])
#         e2 = re.split("\|", elements[1])
#
#         species = e2[-2]
#         cluster = e2[0]
#         spname_to_cluster[species] = cluster
#     return spname_to_cluster
#


def parse_taxonomy(taxon_name):  # given a taxon name, try to return whatever taxonomic info is available as a list starting with the highest level classification and going lower (or a map?)
    #name_elements = re.split("\|", taxon_name)
    name_elements = taxon_name.split('|')
    #print('name_elements')
    #print(name_elements)

    if (len(name_elements) < 8) or (len(name_elements) > 9):
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


def summarize_taxonomy(name_list, tax_level, name_to_tax_dict):  # take a list of names from a clade and summarize taxonomic info (labels and their frequencies)
    total_size = len(name_list)  # it perhaps makes sense to normalize by the size of the clade
    breakdown = {}
    for name in name_list:
        info = name_to_tax_dict[name]
        if info[tax_level] in breakdown:
            breakdown[info[tax_level]] += 1.0 / float(total_size)
        else:
            breakdown[info[tax_level]] = 1.0 / float(total_size)
    return breakdown


def count_sister_taxa(args):

    # read in argument
    target_label = args['l']
    tree_in_ml   = args['ml']
    tree_in_bs   = args['bs']
    output_file  = args['out']

    # edit target_label to make the comparisons at a desired taxonomic level

    # compute the most frequent sister group of each (monophyletic?) group on the tree, to identify trends in gene transfers, "unstable" taxa, etc.

    # read the ML tree, set up the taxonomy stuff, and calculate the number of clades per label, and the sizes of those clades (to report at the end)
    labels = {}
    name_to_tax_info = defaultdict(dict)
    all_tree_leaf_names = []
    ml_tree = Tree(tree_in_ml)  # note that ete3 treats this input tree as rooted
    for leaf in ml_tree:
        taxonomy = parse_taxonomy(leaf.name)
        name_to_tax_info[leaf.name] = taxonomy
        all_tree_leaf_names.append(leaf.name)
        leaf.add_feature("tax", taxonomy[target_label])
        labels[taxonomy[target_label]] = 1
    groups = labels.keys()

    # compute the number of clades (weizhi: monophyletic group) per label in the ML tree, and their sizes
    ML_groups = defaultdict(list)  # the list is the size of each clade, len(list) is the number of clades for that label in the ML tree
    # ML_groups: the number of leaves in each monophyletic groups of the corresponding target_label (e.g. genus)
    for label in groups:
        node_num = 0
        for monophyletic_clade in ml_tree.get_monophyletic(values=[label], target_attr="tax"):  # get monophyletic groups for each target_label (e.g. genus)
            # print('node')
            # print(node)
            size_clade = 0  # get the number of leaf (size_clade) in the monophyletic group
            for leaf in monophyletic_clade:
                size_clade += 1
            ML_groups[label].append(size_clade)
            node_num += 1
        # print('The number of monophyletic clade (node_num) of %s (label):\t%s' % (label, node_num))

    # print()
    # print('Dict holds the number of monophyletic clades per taxon, and their sizes')
    # print('ML_groups:\t %s' % ML_groups)
    # print()

    summary = defaultdict(dict)
    clades_per_group = defaultdict(list)
    treeNum = -1
    for line in open(tree_in_bs):  # read in each bootstrap tree

        treeNum += 1
        tree = Tree(line.rstrip())
        for leaf in tree:
            tax = name_to_tax_info[leaf.name]  # this should set up taxonomy correctly...
            leaf.add_feature("tax", tax[target_label])  # this adds a feature called tax to the leaf, with the attribute of the phylum name
        for label in groups:
            clades_per_group[label].append(0.0)  # setup the clade counting for this particular tree
        tree.unroot()  # Weizhi: why is this

        # iterate over groups that are monophyletic for the taxon label of choice.
        # Choose the smallest sister branch for the comparison. (Assume root is within the larger sister clade (Weizhi:why?))
        for label in groups:

            # print('tree.get_monophyletic(values=[label], target_attr="tax")')
            # print(tree.get_monophyletic(values=[label], target_attr="tax"))
            # print('---------------------------------------------------------------------------------------------------v')
            # print('label: %s' % label)
            monophyletic_clade_index = 1
            for monophyletic_clade in tree.get_monophyletic(values=[label], target_attr="tax"):  # node: monophyletic clade
                clades_per_group[label][treeNum] += 1.0
                # print node.get_ascii()
                sister_clades = monophyletic_clade.get_sisters()

                # print('--------------------v')
                # print('monophyletic clade %s in %s (label)' % (monophyletic_clade_index, label))
                monophyletic_clade_index += 1
                #print(monophyletic_clade)
                # print(len(sisters))
                # for leaf in sisters[0]:
                #     print(leaf.name)
                # print(sisters)
                # print('sisters of current monophyletic clade')
                sister_index = 1
                for each_sister in sister_clades:
                    current_sister_leaf_list = []
                    for leaf in each_sister:
                        current_sister_leaf_list.append(leaf.name)
                    # print('sister %s has %s leaves: %s' % (sister_index, len(current_sister_leaf_list), ','.join([])))
                    sister_index += 1

                if monophyletic_clade.is_root():  # monophyletic clade is root
                    continue

                # Weizhi: bifurcation
                elif len(sister_clades) == 1:  # not at the trifurcation. Do something a bit hacky to find the bigger sister clade

                    taxa_in_sister = []
                    for leaf in sister_clades[0]:
                        taxa_in_sister.append(leaf.name)

                    size_sister = len(taxa_in_sister)

                    taxa_in_group = []
                    for leaf in monophyletic_clade:
                        taxa_in_group.append(leaf.name)

                    taxa_in_other_groups = []  # what does OG mean? (other groups?)
                    for leaf_name in all_tree_leaf_names:
                        if leaf_name in taxa_in_sister:
                            continue
                        elif leaf_name in taxa_in_group:
                            continue
                        else:
                            taxa_in_other_groups.append(leaf_name)
                    size_other_groups = len(taxa_in_other_groups)

                    sister_tax = {}  # taxa in the smaller groups (either the sister group or the OG)
                    if size_other_groups > size_sister:
                        sister_tax = summarize_taxonomy(taxa_in_sister, target_label, name_to_tax_info)
                    else:
                        sister_tax = summarize_taxonomy(taxa_in_other_groups, target_label, name_to_tax_info)

                    # print('size_sister: %s' % size_sister)
                    # print('size_other_groups: %s' % size_other_groups)
                    # print('sister_tax (not really, actually taxa in the smaller one (either the sister group or the OG))')
                    # print(sister_tax)

                    # store the tax info of the sister group
                    for element in sister_tax:
                        # print('element: %s' % element)
                        #print('summary[label]: %s' % summary[label])
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                            #print('summary (in): %s' % summary)
                        else:
                            summary[label][element] = sister_tax[element]
                            #print('summary (not in): %s' % summary)

                else:  # trifurcation in tree. Just treat the two sisters in the same way.

                    taxa_in_sisters_1 = []
                    for leaf in sister_clades[0]:
                        taxa_in_sisters_1.append(leaf.name)

                    taxa_in_sisters_2 = []
                    for leaf in sister_clades[1]:
                        taxa_in_sisters_2.append(leaf.name)

                    # get the size of two sisters
                    size_s1 = len(taxa_in_sisters_1)
                    size_s2 = len(taxa_in_sisters_2)

                    # print('size_s1: %s' % size_s1)
                    # print('size_s2: %s' % size_s2)

                    # get taxa in the smaller sister group
                    sister_tax = {}
                    if size_s1 > size_s2:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_2, target_label, name_to_tax_info)
                    else:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_1, target_label, name_to_tax_info)

                    # print('sister_tax (taxa in the smaller sister group)')
                    # print(sister_tax)

                    for element in sister_tax:
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                        else:
                            summary[label][element] = sister_tax[element]

            #     print('--------------------^')
            #     print()
            # print('---------------------------------------------------------------------------------------------------^')

    # now print out some kind of summary. For each label, the sorted list of sister taxa and their frequencies?
    outh = open(output_file, "w")
    for label in summary:
        num_groups = len(ML_groups[label])
        size_str = ''
        if num_groups == 1:
            size_str = ML_groups[label][0]
        else:
            size_str = ','.join(str(x) for x in (sorted(ML_groups[label], reverse=True)))

        avg_num_clades = float("{0:.4f}".format(numpy.mean(clades_per_group[label])))
        total_num_clades = numpy.sum(clades_per_group[label])
        sorted_sisters = sorted(summary[label].items(), key=itemgetter(1), reverse=True)

        # if label == 'g__TA-20':
        #     print('ML_groups[label]:\t%s'               % ML_groups[label])
        #     print('clades_per_group[label]:\t%s'        % clades_per_group[label])
        #     print('avg_num_clades (mean of list):\t%s'  % numpy.mean(clades_per_group[label]))
        #     print('total_num_clades (sum of list):\t%s' % total_num_clades)
        #     print('summary[label]:\t%s'                 % summary[label])
        #     print('sorted_sisters\t%s'                  % sorted_sisters)

        for tup in sorted_sisters:
            double_normalize = float(tup[1]) / float(total_num_clades)  # normalize the frequencies by the total number of clades, to account for different bootstrap numbers/MCMC sample numbers
            double_normalize = float("{0:.4f}".format(double_normalize))

            str_to_write = '%s\t%s\t%s\t%s\t%s\t%s' % (label,                            # col 1
                                                       tup[0],                           # col 2
                                                       float("{0:.4f}".format(tup[1])),  # col 3
                                                       avg_num_clades,                   # col 4
                                                       double_normalize,                 # col 5
                                                       size_str)                         # col 6
            outh.write(str_to_write + '\n')
            # if label == 'g__TA-20':
            #     print(str_to_write)
    outh.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-ml',  required=True, help='ML tree')
    parser.add_argument('-bs',  required=True, help='bootstrap trees')
    parser.add_argument('-l',   required=True, help='label, choose from cluster, phylum,family, genus...')
    parser.add_argument('-out', required=True, help='output file')
    args = vars(parser.parse_args())
    count_sister_taxa(args)
