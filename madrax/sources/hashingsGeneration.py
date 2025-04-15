#   Copyright 2022 Gabriele Orlando
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


import torch,os,pickle
from madrax.sources.globalVariables import *

def generate_atom_hashing():
	atomsOfEachresi = {
		'LYS': ['C', 'CA', 'CB', 'CD', 'CE', 'CG', 'N', 'NZ', 'O', "tN", "OXT"],
		'ILE': ['C', 'CA', 'CB', 'CD1', 'CG1', 'CG2', 'N', 'O', "tN", "OXT"],
		'THR': ['C', 'CA', 'CB', 'CG2', 'N', 'O', 'OG1', "tN", "OXT"],
		'VAL': ['C', 'CA', 'CB', 'CG1', 'CG2', 'N', 'O', "tN", "OXT"],
		'HIS': ['C', 'CA', 'CB', 'CD2', 'CE1', 'CG', 'N', 'ND1', 'NE2', 'O',"ND1H1S","ND1H2S","NE2H1S","NE2H2S", "tN", "OXT"],
		'PHE': ['C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O',"RC","RE", "tN", "OXT"],
		'ASN': ['C', 'CA', 'CB', 'CG', 'N', 'ND2', 'O', 'OD1', "tN", "OXT"],
		'ARG': ['C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'N', 'NE', 'NH1', 'NH2', 'O',"ARG", "tN", "OXT"],
		'ASP': ['C', 'CA', 'CB', 'CG', 'N', 'O', 'OD1', 'OD2', "tN", "OXT"],
		'GLY': ['C', 'CA', 'N', 'O', "tN", "OXT"],
		'GLU': ['C', 'CA', 'CB', 'CD', 'CG', 'N', 'O', 'OE1', 'OE2', "tN", "OXT"],
		'LEU': ['C', 'CA', 'CB', 'CD1', 'CD2', 'CG', 'N', 'O', "tN", "OXT"],
		'SER': ['C', 'CA', 'CB', 'N', 'O', 'OG', "tN", "OXT"],
		'GLN': ['C', 'CA', 'CB', 'CD', 'CG', 'N', 'NE2', 'O', 'OE1', "tN", "OXT"],
		'ALA': ['C', 'CA', 'CB', 'N', 'O', "tN", "OXT"],
		'CYS': ['C', 'CA', 'CB', 'N', 'O', 'SG', "tN", "OXT"],
		'MET': ['C', 'CA', 'CB', 'CE', 'CG', 'N', 'O', 'SD', "tN", "OXT"],
		'TYR': ['C', 'CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CG', 'CZ', 'N', 'O', 'OH',"RC","RE", "tN", "OXT"],
		'PRO': ['C', 'CA', 'CB', 'CD', 'CG', 'N', 'O', "tN", "OXT"],
		'TRP': ['C', 'CA', 'CB', 'CD1', 'CD2', 'CE2', 'CE3', 'CG', 'CH2', 'CZ2', 'CZ3', 'N', 'NE1', 'O', "tN", "OXT"]}

	max_number_of_atoms = []
	for k in atomsOfEachresi.keys():
		max_number_of_atoms += [len(atomsOfEachresi[k])]

	max_number_of_atoms = max(max_number_of_atoms)

	atom_hash = {}
	bb_atoms = {}#{"C":0,"CA":1,"N":2,"O":3}
	atom_number_hash = len(bb_atoms.keys())

	#atom_to_coords_mapping = np.array()
	for res in atomsOfEachresi.keys():
		if not res in atom_hash:
			atom_hash[res]={}
		for atom in atomsOfEachresi[res]:

			if not atom in bb_atoms:
				atom_hash[res][atom] = atom_number_hash
				atom_number_hash+=1
			else:
				atom_hash[res][atom] = bb_atoms[atom]

	#atom_to_coords_mapping = np.zeros((atom_number_hash))
	#for res in atom_hash.keys():
	#	for cont,atom in enumerate(atom_hash[res]):
	#		atom_to_coords_mapping[atom_hash[res][atom]] = cont

	return 	atom_hash

def generateTensorMappings(atom_hash,property_hashings,property_hashingsFake):
	nH = 12
	nFO = 3
	maxAtomNumber = []
	for r in atom_hash.keys():
		maxAtomNumber+=[max(atom_hash[r].values())]
	maxAtomNumber = max(maxAtomNumber)+1

	path = os.path.dirname(os.path.abspath(__file__))
	
	solvenergy_props, hbond_params, partners, hydrogen_positions, freeOrb_positions = pickle.load(open(path+"/../parameters/parameters.pickle","rb"))

	### partners ###
	partnersT = torch.full((maxAtomNumber,4),PADDING_INDEX)
	for res in partners.keys():
		for at in partners[res].keys():
			at_numb = atom_hash[res][at]

			partnersT[at_numb,0] = atom_hash[res][partners[res][at][0][0]]
			partnersT[at_numb, 1] = partners[res][at][0][1]
			partnersT[at_numb, 2] = atom_hash[res][partners[res][at][1][0]]
			partnersT[at_numb, 3] = partners[res][at][1][1]


	### atom_properties ###
	numberOfProperties=0
	for k in property_hashings.keys():
		numberOfProperties += len(property_hashings[k].keys())

	numberOfFakeProperties = len(property_hashingsFake.keys())

	atom_properties = torch.full((maxAtomNumber,numberOfProperties),PADDING_INDEX)
	fake_atom_Properties = torch.full((maxAtomNumber,nH+nFO ,numberOfFakeProperties),PADDING_INDEX)

	virtual_atom_list = set(["NE2H1S","NE2H2S", "ND1H1S", "ND1H2S","Ncp","Ccp","RE","RC","ARG"])
	culo=[]
	for res in atom_hash.keys():
		for at in atom_hash[res].keys():

			for prop in property_hashings["solvenergy_props"].keys():

				#modifica solvenergy e hbond_params per far si che abbiano iresidui e atomi giusti, anche per i virtuali
				#aggiungi le proprietà extra, es virtuale e tutto quello riportato nel dataStructure vecchio
				if at in solvenergy_props[res]:
					at_numb = atom_hash[res][at]
					atom_properties[at_numb,property_hashings["solvenergy_props"][prop]] = solvenergy_props[res][at][prop]
				elif  at in solvenergy_props["XXX"]: ### virtual
					at_numb = atom_hash[res][at]
					atom_properties[at_numb,property_hashings["solvenergy_props"][prop]] = solvenergy_props["XXX"][at][prop]
				else:
					asd

			for prop in property_hashings["hbond_params"].keys():

				#if not res in hbond_params:
				#	continue
				if res in hbond_params and at in hbond_params[res]:
					at_numb = atom_hash[res][at]
					atom_properties[at_numb,property_hashings["hbond_params"][prop]] = hbond_params[res][at][prop]
				elif  at in hbond_params["XXX"]: ### virtual
					at_numb = atom_hash[res][at]
					atom_properties[at_numb,property_hashings["hbond_params"][prop]] = hbond_params["XXX"][at][prop]
				elif at in hbond_params["ooo"]:
					at_numb = atom_hash[res][at]
					atom_properties[at_numb,property_hashings["hbond_params"][prop]] = hbond_params["ooo"][at][prop]
				elif at == "ND1H1S" or at == "ND1H2S" or at == "NE2H1S" or at == ["NE2H2S"]:
					at_numb = atom_hash[res][at]
					r = at[3:]
					a = at[:3]
					atom_properties[at_numb,property_hashings["hbond_params"][prop]] = hbond_params[r][a][prop]
				else:
					culo+=[(res,at)]

			if at in virtual_atom_list:
				atom_properties[at_numb, property_hashings["other_params"]["virtual"]] = 1
			else:
				atom_properties[at_numb, property_hashings["other_params"]["virtual"]] = 0

			## close atoms ##
			if at == "O" or at == "OXT":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]] = 0
			elif at == "N" or at == "tN":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 2 
			elif at == "C":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 2
			elif at == "CA":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 3
			elif at[1] == "B":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 4
			elif at[1] == "G":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 5
			elif at[1] == "D":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 6
			elif at[1] == "E":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 7
			elif at[1] == "Z":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 8
			elif at == "OH" and res == "TYR":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 9
			elif (at == "NH2" or at == "NH1") and res == "ARG":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 9
			elif (at == "CH2" or at == "CH2") and res == "TRP":
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = 9
			elif at in virtual_atom_list:
				atom_properties[at_numb, property_hashings["other_params"]["close_atoms"]]  = PADDING_INDEX
			else:
				asd
			#####

			if (res in freeOrb_positions and at in freeOrb_positions[res]):

				for foNumber,fo in enumerate(freeOrb_positions[res][at].keys()):

					at_numb = atom_hash[res][at]
					rotation = int(fo.split("_")[1])
					positionInTensor = nH + foNumber

					assert positionInTensor>=12


					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionX"]] = freeOrb_positions[res][at][fo][property_hashingsFake["relative_positionX"]]
					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionY"]] = freeOrb_positions[res][at][fo][property_hashingsFake["relative_positionY"]]
					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionZ"]] = freeOrb_positions[res][at][fo][property_hashingsFake["relative_positionZ"]]
					fake_atom_Properties[at_numb, positionInTensor, property_hashingsFake["rotation"]] = rotation

			if (res in hydrogen_positions and at in hydrogen_positions[res]):

				for hNumber,h in enumerate(hydrogen_positions[res][at].keys()):

					at_numb = atom_hash[res][at]
					rotation = int(h.split("_")[1])
					positionInTensor = hNumber
					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionX"]] = hydrogen_positions[res][at][h][property_hashingsFake["relative_positionX"]]
					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionY"]] = hydrogen_positions[res][at][h][property_hashingsFake["relative_positionY"]]
					fake_atom_Properties[at_numb,positionInTensor,property_hashingsFake["relative_positionZ"]] = hydrogen_positions[res][at][h][property_hashingsFake["relative_positionZ"]]
					fake_atom_Properties[at_numb, positionInTensor, property_hashingsFake["rotation"]] = rotation

	return atom_properties,fake_atom_Properties,partners
