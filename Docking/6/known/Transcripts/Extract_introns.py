from Bio import SeqIO
import pandas as pd

seqs = SeqIO.parse("chr_3.fa","fasta")

exons = pd.read_csv("Exonic_parts_info.csv")

for seq in seqs:
    for num,exon in enumerate(exons["Exonic part"]):
        print(">"+"Pn16.1237_"+str(exon))
        print(seq.seq[exons["Start"][num]-1:exons["End"][num]])