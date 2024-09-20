import sys

from rdkit import Chem

from rdkit.Chem import AllChem

file = str(sys.argv[1])

m = Chem.MolFromMolFile(file)

m2 = Chem.AddHs(m)

AllChem.EmbedMolecule(m2)

AllChem.ComputeGasteigerCharges(m2)

AllChem.MMFFOptimizeMolecule(m2)

file2 = file.split(".")[0]+"_minimized.mol"

print(Chem.MolToMolBlock(m2),file=open(file2,'w+'))


