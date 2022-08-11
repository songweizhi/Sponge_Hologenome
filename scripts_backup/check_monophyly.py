from ete3 import Tree

t = Tree("((((((a, e), i), o),h), u), ((f, g), j));")
print(t.check_monophyly(values=["a", "i", "o", "h"], target_attr="name"))
