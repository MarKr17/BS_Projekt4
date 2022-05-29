#script for finding potencial pairs in PDB file
#Ogarnianie jak poruszac się pod strukturze pdb w Biopythonie
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
parser=PDBParser()
structure=parser.get_structure("rp01","rp/rp01.pdb")
print(structure)
model= structure.get_models()
print(model)
models=list(model)
print(models)
chains=list(models[0].get_chains())
print(chains)
residueA=list(chains[0].get_residues())
print(len(residueA))
residueB=list(chains[1].get_residues())
print(len(residueB))
atoms=list(residueA[0].get_atoms())
print(atoms)
vector=list(atoms[0].get_vector())
print(vector)
print(residueA)