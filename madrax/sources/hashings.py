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

resi_hash={ 'CYS': 0, 'ASP': 1, 'SER': 2, 'GLN': 3, 'LYS': 4, 'ASN': 5, 'PRO': 6, 'THR': 7, 'PHE': 8, 'ALA': 9, 'HIS': 10, 'GLY': 11, 'ILE': 12, 'LEU': 13, 'ARG': 14, 'TRP': 15, 'VAL': 16, 'GLU': 17, 'TYR': 18, 'MET': 19}
hashing_hybrid = {"NO_HYBRID":0,"SP2_N_H1":1,"SP2_N_H2":2,"SP2_N_ORB1":3,"SP2_N_ORB2":4,"SP3_N_H3":5,"SP3_O_H1ORB2":6,"SP2_O_ORB2":7}

resi_hash_inverse = {}
for i in resi_hash.keys():
	resi_hash_inverse[resi_hash[i]] = i

atom_description_hash = {
						"batch": 0,
						"chain": 1,
						"resnum": 2,
						"resname": 3,
						"at_name": 4,
						"alternative": 5,
						"alternative_conf": 6
						}
property_hashings = {		  "solvenergy_props":
									{
									"radius":0,
									"Occmax":1,
									"volume":2,
									"Occ":3,
									"VdW":4,
									"enerG":5,
									"minradius":6,
                                    "level":7,
                                    "cycle":8,
									},
							  "hbond_params":
									{
									"hbond_plane_dist":9,
									"charge":10,
									"dipole":11,
									"hybridation":12,
									"charged":13,
									"donor":14,
									"acceptor":15,
									},
							  "other_params":
								  {
								  "virtual": 16,
								  "close_atoms": 17,
								  }
							  }

property_hashingsFake = {
						"relative_positionX":0,
						"relative_positionY":1,
						"relative_positionZ":2,
						"rotation":3,
						}
level_hashing = {
				"LEVEL_N":0,
				"LEVEL_A":1,
				"LEVEL_O":-1,
				"LEVEL_B":2,
				"LEVEL_G":3,
				"LEVEL_D":4,
				"LEVEL_E": 5,
				"LEVEL_Z": 6,
				"LEVEL_H": 7,
				"LEVEL_I": 8,
				"LEVEL_K": 9,
	}
#### angle calculation ####
chi_atoms = dict(
	chi1=dict(ARG=['N', 'CA', 'CB', 'CG'], ASN=['N', 'CA', 'CB', 'CG'], ASP=['N', 'CA', 'CB', 'CG'],
			  CYS=['N', 'CA', 'CB', 'SG'], GLN=['N', 'CA', 'CB', 'CG'], GLU=['N', 'CA', 'CB', 'CG'],
			  HIS=['N', 'CA', 'CB', 'CG'], ILE=['N', 'CA', 'CB', 'CG1'], LEU=['N', 'CA', 'CB', 'CG'],
			  LYS=['N', 'CA', 'CB', 'CG'], MET=['N', 'CA', 'CB', 'CG'], PHE=['N', 'CA', 'CB', 'CG'],
			  PRO=['N', 'CA', 'CB', 'CG'], SER=['N', 'CA', 'CB', 'OG'], THR=['N', 'CA', 'CB', 'OG1'],
			  TRP=['N', 'CA', 'CB', 'CG'], TYR=['N', 'CA', 'CB', 'CG'],
			  VAL=['N', 'CA', 'CB', 'CG1'], ),

	chi2=dict(ARG=['CA', 'CB', 'CG', 'CD'], ASN=['CA', 'CB', 'CG', 'OD1'],
			  ASP=['CA', 'CB', 'CG', 'OD1'],
			  GLN=['CA', 'CB', 'CG', 'CD'], GLU=['CA', 'CB', 'CG', 'CD'],
			  HIS=['CA', 'CB', 'CG', 'ND1'],
			  ILE=['CA', 'CB', 'CG1', 'CD1'], LEU=['CA', 'CB', 'CG', 'CD1'],
			  LYS=['CA', 'CB', 'CG', 'CD'],
			  MET=['CA', 'CB', 'CG', 'SD'], PHE=['CA', 'CB', 'CG', 'CD2'],
			  PRO=['CA', 'CB', 'CG', 'CD'],
			  TRP=['CA', 'CB', 'CG', 'CD1'], TYR=['CA', 'CB', 'CG', 'CD1'], ),

	chi3=dict(ARG=['CB', 'CG', 'CD', 'NE'], GLN=['CB', 'CG', 'CD', 'OE1'],
			  GLU=['CB', 'CG', 'CD', 'OE1'],
			  LYS=['CB', 'CG', 'CD', 'CE'], MET=['CB', 'CG', 'SD', 'CE'], ),

	chi4=dict(ARG=['CG', 'CD', 'NE', 'CZ'], LYS=['CG', 'CD', 'CE', 'NZ'], ),

	chi5=dict(ARG=['CD', 'NE', 'CZ', 'NH1'], ),

)

res_acid_point ={
				"PRO":[1125,474],
				"GLN":[1875,731],
				"VAL":[1125,512],
				"ASN":[1500,593],
				"THR":[1125,498],
				"ALA":[375,120],
				"ASP":[1500,580],
				"PHE":[2625,1072],
				"LEU":[1500,627],
				"SER":[750,271],
				"CYS":[750,266],
				"ILE":[1500,633],
				"TRP":[3750,1757],
				"ARG":[2625,1003],
				"LYS":[1875,744],
				"TYR":[3000,1263],
				"GLU":[1865,679],
				"MET":[1500,607],
				"HIS":[2250,971]
				}
res_acid_pointMC = { "LEU":[1500,834],
					 "ASN":[1500,759],
					 "CYS":[1500,746],
					 "ASP":[1500,745],
					 "GLY":[1500,595],
					 "GLN":[1500,780],
					 "VAL":[1500,804],
					 "SER":[1500,639],
					 "MET":[1500,799],
					 "LYS":[1500,797],
					 "PRO":[1500,789],
					 "CYS":[1500,746],
					 "LEU":[1500,834],
					 "THR":[1500,766],
					 "TYR":[1500,836],
					 "HIS":[1500,836],
					 "GLU":[1500,767],
					 "ALA":[1500,692],
					 "ILE":[1500,842],
					 "ARG":[1500,780],
					 "PHE":[1500,836],
					 "ILE":[1500,842],
					 "TRP":[1500,865]
					}


chi_angle_builder_hashing = {}
for chi in chi_atoms.keys():
    for res in chi_atoms[chi].keys():
        if not res in chi_angle_builder_hashing:
            chi_angle_builder_hashing[res] = {}
        for i,at in enumerate(chi_atoms[chi][res]):
            if not at in chi_angle_builder_hashing[res]:
                chi_angle_builder_hashing[res][at] = []
            chi_angle_builder_hashing[res][at] += [(chi, i)]

####


from madrax.sources import hashingsGeneration
atom_hash = hashingsGeneration.generate_atom_hashing()
atom_Properties, fake_atom_Properties,partners = hashingsGeneration.generateTensorMappings(atom_hash,property_hashings,property_hashingsFake)

atom_hash_inverse = {}
for i in atom_hash.keys():
	for k in atom_hash[i].keys():
		atom_hash_inverse[atom_hash[i][k]] = k

atom_hashTPL = {

	'CYS': {
		"N": 13,
		"O": 15,
		"C": 14,
		"SG": 1,
		"CB": 11,
		"CA": 12,

	},
	'ASP': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"OD": 8,
		"OD1": 8,
		"OD2": 8,
		"CG": 9,
		"CB": 11
	},
	'SER': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"OG": 7,
		"CB": 11,

	},
	'GLN': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"NE2": 2,
		"OE1": 6,
		"CD": 9,
		"CB": 11,
		"CD": 11,
		"CG": 11,

	},
	'LYS': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"NZ": 5,
		"CB": 11,
		"CG": 11,
		"CD": 11,
		"CE": 11,
	},
	'ASN': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"ND2": 2,
		"OD1": 6,
		"CG": 9,
		"CB": 11,

	},
	'PRO': {
		"N": 16,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
		"CG": 11,
		"CD": 11,

	},
	'THR': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"OG1": 7,
		"CB": 11,
		"CG2": 11,

	},
	'PHE': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CG": 10,
		"CD1": 10,
		"CD2": 10,
		"CE2": 10,
		"CE1": 10,
		"CZ": 10,
		"CB": 11,
	},
	'ALA': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
	},
	'HIS': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"ND1": 3,
		"ND2": 3,
		"NE1": 3,
		"NE2": 3,
		"CG": 10,
		"CD2": 10,
		"CE1": 10,
		"CB": 11,

	},
	'GLY': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
	},
	'ILE': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
		"CG2": 11,
		"CG1": 11,
		"CD1": 11
	},
	'LEU': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
		"CG": 11,
		"CD1": 11,
		"CD2": 11,
		"CD3": 11,
	},
	'ARG': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
		"CD": 11,
		"CG": 11,
		"NE": 4,
		"NH1": 4,
		"NH2": 4,
		"NH3": 4,
		"CZ": 9,

	},
	'TRP': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"NE1": 3,
		"CG": 10,
		"CD1": 10,
		"CD2": 10,
		"CE1": 10,
		"CE2": 10,
		"CE3": 10,
		"CZ1": 10,
		"CZ2": 10,
		"CZ3": 11,
		"CH2": 10,
		"CB": 11,

	},
	'VAL': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"CB": 11,
		"CG1": 11,
		"CG2": 11,
		"CG3": 11,
	},
	'GLU': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"OE1": 8,
		"OE2": 8,
		"CD": 9,
		"CB": 11,
		"CG": 11,

	},
	'TYR': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"OH": 7,
		"CG": 10,
		"CD1": 10,
		"CD2": 10,
		"CE1": 10,
		"CE2": 10,
		"CZ": 10,
		"CB": 11,

	},
	'MET': {
		"N": 13,
		"O": 15,
		"C": 14,
		"CA": 12,
		"SD": 1,
		"CB": 11,
		"CG": 11,
		"CE": 11,
	}
}