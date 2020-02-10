
from ete3 import Tree
tree = Tree('./trees/final/83g_merged.newick',format=3)
tree = tree.write(format=0).replace(')1:0','):0')
align_file = '/mnt/home-backup/thliao/plancto/dating_for/beast_for/83g_files/alin.txt'

align_txt = open(align_file).read().split('\n')
align_txt = '\n'.join(align_txt[1:])
template=f"""#NEXUS
begin data;
dimensions ntax=83 nchar=6674;
Format datatype=protein gap=-;
matrix
{align_txt}
;
end;
begin assumptions;
    charset part1 = 1-136;
    charset part2 = 137-546;
    charset part3 = 547-762;
    charset part4 = 763-903;
    charset part5 = 904-1100;
    charset part6 = 1101-1978;
    charset part7 = 1979-2494;
    charset part8 = 2495-2686;
    charset part9 = 2687-2957;
    charset part10 = 2958-3137;
    charset part11 = 3138-3296;
    charset part12 = 3297-3449;
    charset part13 = 3450-3594;
    charset part14 = 3595-3717;
    charset part15 = 3718-3836;
    charset part16 = 3837-3945;
    charset part17 = 3946-4082;
    charset part18 = 4083-4417;
    charset part19 = 4418-4651;
    charset part20 = 4652-4829;
    charset part21 = 4830-5376;
    charset part22 = 5377-5549;
    charset part23 = 5550-5804;
    charset part24 = 5805-6104;
    charset part25 = 6105-6674;
end;

begin trees;
        tree TREE1 = [&R] {tree}
end;
"""

with open('./dating_for/beast_for/83g_files/test.nexus','w') as f1:
    f1.write(template)