import os

file = open('coils_list')

lines = file.readlines()

for l in lines:
    pdb = l.split()[0]
    if len(pdb) == 4:
        print pdb,
