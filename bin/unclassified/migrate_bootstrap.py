

t1 = ''
t2 = ''
from ete3 import Tree

tre1 = Tree(t1)
tre2 = Tree(t2)
tips2bootstrap = {}
for n in tre1.traverse():
    if not n.is_leaf():
        tips2bootstrap[tuple(sorted(n.get_leaf_names()))] = n.support

for n in tre2.traverse():
    if not n.is_leaf():
        n.support = tips2bootstrap.get(tuple(sorted(n.get_leaf_names())),None)
tre2.write(outfile='',format=3)


