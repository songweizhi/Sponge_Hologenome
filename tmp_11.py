from collections import defaultdict


def def_value():
    return "Not Present"


clades_per_group = defaultdict(def_value)

print(clades_per_group['aaa'])

