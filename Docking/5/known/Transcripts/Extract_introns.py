from Bio import SeqIO
import pandas as pd
import sys

seqs = SeqIO.parse("chr_9.fas","fasta")

exons = pd.read_csv("Exonic_parts_info.csv")

original_stdout = sys.stdout # Save a reference to the original standard output

with open('exons.fas', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.

    for seq in seqs:
        for num,exon in enumerate(exons["Exonic part"]):
            print(">"+"Pn3.4770_"+str(exon))
            print(seq.seq[exons["Start"][num]-1:exons["End"][num]])
    sys.stdout = original_stdout # Reset the standard output to its original value