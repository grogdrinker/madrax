import os,torch, string
import numpy as np
from madrax.sources.globalVariables import *
from madrax.sources import hashings
from madrax import dataStructures


def parsePDB(PDBFile,keep_only_chains=None,bb_only=False):

	"""
	function to parse pdb files. It can be used to parse a single file or all the pdb files in a folder. In case a folder is given, the coordinates are gonna be padded

	Parameters
	----------
	PDBFile : str
		path of the PDB file or of the folder containing multiple PDB files
	bb_only : bool
		if True ignores all the atoms but backbone N, C and CA
	keep_only_chains : str or None
		ignores all the chain but the one given. If None it keeps all chains
	keep_hetatm : bool
		if False it ignores heteroatoms
	Returns
	-------
	coords : torch.Tensor
		coordinates of the atoms in the pdb file(s). Shape ( batch, numberOfAtoms, 3)

	atomNames : list
		A list of the atom identifier. It encodes atom type, residue type, residue position and chain
		as an example GLN_39_N_B_0_0 refers to an atom N in a Glutamine, in position 39 of chain B. the last two zeros are used for the mutation engine and should be ignored

	pdbNames : list
		an ordered list of the structure names

	"""

	bbatoms = ["N", "CA", "C"]
	ring_residues = ["PHE","TYR"]
	R_C_atoms = ["CD1", "CE2"]
	pdbNames = []

	if not os.path.isdir(PDBFile):
		pdbNames=[PDBFile.replace(".pdb","")]
		fil = PDBFile
		atomNamesTMP = []
		coordsTMP = []
		RC_atoms_dict = {}
		cont = -1
		oldres = -999

		for line in open(fil).readlines():
			if line[:6] == "ENDMDL":
				break
			elif line[:4] == "ATOM":
				if keep_only_chains is not None and (not line[21] in keep_only_chains):
					continue
				if bb_only and not line[12:16].strip() in bbatoms:
					continue
				if oldres != int(line[22:26]):
					cont += 1
					oldres = int(line[22:26])
				if line[12:16].strip()[0] == "H":  # hydrogens are removed
					continue
				resnum = int(line[22:26])
				altern_conformation = "0"
				mutant = "0"

				atomNamesTMP += [line[17:20].strip() + "_" + str(resnum) + "_" + line[12:16].strip() + "_" + line[
					21] + "_" + mutant + "_" + altern_conformation]

				x = float(line[30:38])
				y = float(line[38:46])
				z = float(line[46:54])

				if line[17:20].strip() in ring_residues:
					if not resnum in RC_atoms_dict:
						name = line[17:20].strip() + "_" + str(resnum) + "_" + "RC" + "_" + line[
							21] + "_" + mutant + "_" + altern_conformation
						RC_atoms_dict[resnum] = ({}, name)
					if line[12:16].strip() in R_C_atoms:
						RC_atoms_dict[resnum][0][line[12:16].strip()] = [x, y, z]
				coordsTMP += [[x, y, z]]

		for r_c in RC_atoms_dict.keys():
			if ("CD1" in RC_atoms_dict[r_c][0] and "CE2" in RC_atoms_dict[r_c][0]):
				new_coords = [(RC_atoms_dict[r_c][0]["CD1"][0] + RC_atoms_dict[r_c][0]["CE2"][0]) / 2.0,
							  (RC_atoms_dict[r_c][0]["CD1"][1] + RC_atoms_dict[r_c][0]["CE2"][1]) / 2.0,
							  (RC_atoms_dict[r_c][0]["CD1"][2] + RC_atoms_dict[r_c][0]["CE2"][2]) / 2.0]
				coordsTMP += [new_coords]
				atomNamesTMP += [RC_atoms_dict[r_c][1]]

		return torch.tensor(coordsTMP).unsqueeze(0),[atomNamesTMP]

	else:
		coords = []
		atomNames = []

		for c,fil in enumerate(sorted(os.listdir(PDBFile))):
			#print(fil,c)
			pdbNames += [fil.replace(".pdb", "")]
			atomNamesTMP = []
			coordsTMP = []
			RC_atoms_dict = {}
			cont = -1
			oldres=-999
			for line in open(PDBFile+"/"+fil).readlines():
				if line[:6] == "ENDMDL":
					break
				elif line[:4] == "ATOM":
					if keep_only_chains is not None and (not line[21] in keep_only_chains):
						continue
					if bb_only and not line[12:16].strip() in bbatoms:
						continue
					if oldres != int(line[22:26]):
						cont += 1
						oldres = int(line[22:26])
					if line[12:16].strip()[0] == "H": # hydrogens are removed
						continue
					resnum = int(line[22:26])
					altern_conformation = "0"
					mutant = "0"

					atomNamesTMP += [line[17:20].strip()+"_"+str(resnum)+"_"+line[12:16].strip()+"_"+line[21]+"_"+mutant+"_"+ altern_conformation]

					x = float(line[30:38])
					y = float(line[38:46])
					z = float(line[46:54])


					if line[17:20].strip() in ring_residues:
						if not resnum in RC_atoms_dict:
							name = line[17:20].strip()+"_"+str(resnum)+"_"+"RC"+"_"+line[21]+"_"+mutant+"_"+ altern_conformation
							RC_atoms_dict[resnum] = ({ },name)
						if line[12:16].strip() in R_C_atoms:
							RC_atoms_dict[resnum][0][line[12:16].strip()] = [x,y,z]
					coordsTMP+=[[x,y,z]]


			for r_c in RC_atoms_dict.keys():
				if  ("CD1" in RC_atoms_dict[r_c][0] and "CE2" in RC_atoms_dict[r_c][0]):
					new_coords = [(RC_atoms_dict[r_c][0]["CD1"][0] + RC_atoms_dict[r_c][0]["CE2"][0] )/ 2.0,
								  (RC_atoms_dict[r_c][0]["CD1"][1] + RC_atoms_dict[r_c][0]["CE2"][1]) / 2.0,
								  (RC_atoms_dict[r_c][0]["CD1"][2] + RC_atoms_dict[r_c][0]["CE2"][2]) / 2.0]
					coordsTMP += [new_coords]
					atomNamesTMP += [RC_atoms_dict[r_c][1]]

			coords+=[torch.tensor(coordsTMP)]
			atomNames += [atomNamesTMP]


		return torch.nn.utils.rnn.pad_sequence(coords, batch_first=True, padding_value=PADDING_INDEX),atomNames,pdbNames

def writepdb(coords,atnames, pdb_names=None, output_folder="outpdb/"):

	"""
	function to write a pdb file from the coordinates and atom names

	Parameters
	----------
	coords : torch.Tensor
		shape: (Batch, nAtoms, 3)
		coordinates of the proteins. It can be generated using the madrax.utils.parsePDB function
	atName : list
		shape: list of lists
		A list of the atom identifier. It encodes atom type, residue type, residue position and chain
		as an example GLN_39_N_B_0_0 refers to an atom N in a Glutamine, in position 39 of chain B. the last two zeros are used for the mutation engine and should be ignored
		This list can be generated using the madrax.utils.parsePDB function
	pdb_names : list
		names of the PDBs. you can get them from the output of utils.parsePDB. If None is given, the proteins are named with an integer that represent their position in the batch
	output_folder : str
		output folder in which PDBs are written
	Returns
	-------
	"""

	chain_hashingTot=[]
	chain_hashing = {}
	for protN, prot in enumerate(atnames):
		for atn, atom in enumerate(prot):
			if not atom.split("_")[3] in chain_hashingTot:
				chain_hashingTot+=[atom.split("_")[3]]

	chain_hashingTot = sorted(chain_hashingTot)
	for i in range(len(chain_hashingTot)):
		chain_hashing[i] = chain_hashingTot[i]

	if pdb_names==None:
		pdb_names = []
		for i in range(len(atnames)):
			pdb_names+=[str(i)]

	info_tensors = dataStructures.create_info_tensors(atnames, device="cpu", verbose=False)

	_, atom_description, coordsIndexingAtom, _, _, altern_mask = info_tensors

	coords = coords[atom_description[:, 0].long(), coordsIndexingAtom]

	os.system("mkdir -p " + output_folder+"/")

	atname = atom_description[:, hashings.atom_description_hash["at_name"]]

	resname = atom_description[:, hashings.atom_description_hash["resname"]]

	batch_ind = atom_description[:, hashings.atom_description_hash["batch"]]

	chain_ind = atom_description[:, hashings.atom_description_hash["chain"]]

	res_num = atom_description[:, hashings.atom_description_hash["resnum"]]

	wtMask = altern_mask[:, 0]

	ca_mask = torch.zeros(coords.shape[0])
	for res in hashings.resi_hash.keys():
		ca_mask += (atname == hashings.atom_hash[res]["CA"])
	ca_mask = ca_mask.bool()

	for batch in range(len(pdb_names)):

		for alt in range(altern_mask.shape[-1]):

			n = pdb_names[batch]

			batchMask = batch_ind.eq(batch) #& altern_mask[:,alt]
			if ((altern_mask[:,alt] !=wtMask) & batchMask).sum()>0:
				### mutation ###

				changing = (altern_mask[:,alt] !=wtMask) & batchMask & ca_mask
				posWT = res_num[changing & wtMask].tolist()
				posMU = res_num[changing & altern_mask[:,alt]].tolist()

				chainWT = chain_ind[changing & wtMask].tolist()#.unique()
				chainMU = chain_ind[changing & altern_mask[:,alt]].tolist()#.unique()

				resnameWT = resname[changing & wtMask].tolist()#.unique()
				resnameMU = resname[changing & altern_mask[:,alt]].tolist()#.unique()

				#mutants = torch.cat([posWT,posMU,chainWT,chainMU,resnameWT,resnameMU],dim=-1)


				name=n
				for k in range(len(resnameMU)):
					nMut=hashings.resi_hash_inverse[resnameWT[k]]+chain_hashing[chainWT[k]]+str(posWT[k])+hashings.resi_hash_inverse[resnameMU[k]]
					name+="_"+nMut
				if len(name)>50:
					name = n+"_mutant"+str(alt)

			elif alt==0:
				name = n
			else:
				continue
			f = open(output_folder + name + ".pdb", "w")
			skip = ["RC","RE"]
			skip_mask = torch.zeros(coords.shape[0])
			for res in hashings.resi_hash.keys():
				for s in skip:
					if s in hashings.atom_hash[res]:
						skip_mask += (atname == hashings.atom_hash[res][s])
			skip_mask = ~skip_mask.bool()

			coords_mask = batchMask & altern_mask[:,alt] & skip_mask

			atnameLoc = atname[coords_mask].tolist()

			resnameLoc = resname[coords_mask].tolist()

			batch_indLoc = batch_ind[coords_mask].tolist()
			chain_indLoc =chain_ind[coords_mask].tolist()

			res_numLoc = res_num[coords_mask].tolist()

			coordsLoc = coords[coords_mask].tolist()
			for i in range(len(res_numLoc)):

				num = " " * (5 - len(str(i))) + str(i)
				at = hashings.atom_hash_inverse[atnameLoc[i]]
				if at == "tN":
					at = "N"
				if at == "OXT":
					at = "O"

				a_name = at + " " * (4 - len(at))
				numres = " " * (4 - len(str(res_numLoc[i]))) + str(res_numLoc[i])

				x = round(float(coordsLoc[i][0]), 3)
				sx = str(x)
				while len(sx.split(".")[1]) < 3:
					sx += "0"
				x = " " * (8 - len(sx)) + sx

				y = round(float(coordsLoc[i][1]), 3)
				sy = str(y)
				while len(sy.split(".")[1]) < 3:
					sy += "0"
				y = " " * (8 - len(sy)) + sy

				z = round(float(coordsLoc[i][2]), 3)
				sz = str(z)
				while len(sz.split(".")[1]) < 3:
					sz += "0"
				z = " " * (8 - len(sz)) + sz
				chain = " " * (2 - len(chain_hashing[chain_indLoc[i]])) + chain_hashing[chain_indLoc[i]]
				# "ATOM  "+num+"  "+a_name+""+resi_list[i]+" A"
				resNam = hashings.resi_hash_inverse[resnameLoc[i]]
				f.write("ATOM  " + num + "  " + a_name + "" + resNam + chain + numres + "    " + x + y + z + "  1.00 64.10           " + at[0] + "\n")
			f.close()

def atomName2Seq(atName):

	"""
	Function that returns the sequence of proteins given the atom names

	Parameters
	----------

	atName : list
		shape: list of lists
		A list of the atom identifier. It encodes atom type, residue type, residue position and chain
		as an example GLN_39_N_B_0_0 refers to an atom N in a Glutamine, in position 39 of chain B. the last two zeros are used for the mutation engine and should be ignored
		This list can be generated using the madrax.utils.parsePDB function

	Returns
	-------
	sequences :
		sequence of the PDBs
	"""

	letters = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ASN': 'N', 'PRO': 'P', 'THR': 'T',
			   'PHE': 'F', 'ALA': 'A', 'HIS': 'H', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
			   'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}  # gli aminoacidi, che male qui non fanno

	seqs=[]


	for batch in atName:
		ch=None
		minres =999
		maxres = -1
		rname = ""
		rnum = {}
		cont=0
		for a in batch:
			resname, resind, a, chains, mut, alt_conf = a.split("_")
			resind = int(resind)
			if not (ch==chains or ch is None):
				raise("this is a testing function, only single chian pdbs are allowed")
			if a=="CA":
				if resind>maxres:
					maxres = resind
				if resind < minres:
					minres = resind
				rname += letters[resname]
				rnum[int(resind)] = cont
				cont+=1
		s =""
		for k in range(0,maxres+1):#range(minres,maxres+1):
			if k in rnum:
				s+=rname[rnum[k]]
			else:
				s+="X"
		seqs+=[s]
	return seqs
