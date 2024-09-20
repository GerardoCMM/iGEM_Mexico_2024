import pandas as pd
import math

import sys
import os

original_stdout = sys.stdout

up = pd.read_csv("upIn40.csv")

dif = pd.read_csv("difIn40.csv")

pos = pd.read_csv("gene_positions.tsv",sep="\t")

res = pd.read_csv("r_crudos.csv")

with open('genes_up.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    for gene in up["forty"]:
	    for i,genename in enumerate(list(pos.loc[:,"name"])):
	        if gene==genename:
	            print(str(pos.loc[i,"Chromosome"])+"\t"+str(int(pos.loc[i,"init"]))+"\t"+str(int(pos.loc[i,"end"])))
    sys.stdout = original_stdout # Reset the standard output to its original value


geneup = list(up["forty"])

genediff = list(dif["difforty"])

genenot = []

for gene in genediff:
    if gene not in geneup:
        genenot.append(gene)


with open('genes_down.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    for gene in genenot:
	    for i,genename in enumerate(list(pos.loc[:,"name"])):
	        if gene==genename:
	            print(str(pos.loc[i,"Chromosome"])+"\t"+str(int(pos.loc[i,"init"]))+"\t"+str(int(pos.loc[i,"end"])))
    sys.stdout = original_stdout # Reset the standard output to its original value

