#script for finding potencial pairs in PDB file
#Ogarnianie jak poruszac się pod strukturze pdb w Biopythonie
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import math
parser=PDBParser()
structure=parser.get_structure("rp02","rp/rp02.pdb")
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

G_C={'N1_N3': [2.67, 3.19],'N2_O2':[2.44, 3.20], 'O6_N4':[2.61, 3.31]}

print(residueA[0].get_resname())
print(chains[0].get_id())
G_list=[]
C_list=[]

wynik=[]
for model in structure:
	for chain in model:
		for residue in chain:
			name=residue.get_resname().strip()
			id=chain.get_id().strip()
			if name=='G':
				G_list.append(residue)
			if name=='C':
				C_list.append(residue)

print(G_list[0]["N1"].get_coord())
print(G_list[0].get_id())
print(C_list[0]["N3"].get_coord())
for G in G_list:
	for C in C_list:
		N1= G['N1'].get_coord()
		N3= C['N3'].get_coord()
		N1_N3= math.sqrt(pow(N1[0]-N3[0],2)+pow(N1[1]-N3[1],2)+pow(N1[2]-N3[2],2))
		if 2.67 <= N1_N3 <=3.19:
			N2=G["N2"].get_coord()
			O2=C["O2"].get_coord()
			N2_O2=math.sqrt(pow(N2[0]-O2[0],2)+pow(N2[1]-O2[1],2)+pow(N2[2]-O2[2],2))
			if 2.44 <= N2_O2 <= 3.20:
				O6=G["O6"].get_coord()
				N4=C["N4"].get_coord()
				O6_N4=math.sqrt(pow(O6[0]-N4[0],2)+pow(O6[1]-N4[1],2)+pow(O6[2]-N4[2],2))
				if 2.61 <= O6_N4 <= 3.31:
					s1="{}:{}{}".format(G.get_parent().get_id(),'G',G.get_id()[1])
					s2="{}:{}{}".format(C.get_parent().get_id(),"C",C.get_id()[1])
					wynik.append(sorted([s1,s2]))

wynik=sorted(wynik)
print(wynik)