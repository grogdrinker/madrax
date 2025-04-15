import torch, multiprocessing,time

from madrax.sources.globalVariables import *
from madrax.sources import hashings
def create_info_tensors(atName,int_dtype=torch.short,device="cpu",verbose=False):
	"""
	function that calculates the actual energy of the protein(s) or complex(es)

	Parameters
	----------
	atName : list
		shape: list of lists
		A list of the atom identifier. It encodes atom type, residue type, residue position and chain
		as an example GLN_39_N_B_0_0 refers to an atom N in a Glutamine, in position 39 of chain B. the last two zeros are used for the mutation engine and should be ignored
		This list can be generated using the madrax.utils.parsePDB function
	int_dtype : torch.dtype
		experimental! it is supposed to be used to run MadraX with different data types
	device : str
		device on which the object should be created. Example: "cpu" for CPU calculation, "cuda" for GPU calculation, "cuda:0" to use GPU 0. Refer to pytorch devices semantic for more information

	Returns
	-------

	info_tensor : tuple
		some precalculated values used to run MadraX
	"""

	terminals={}
	alternative = []
	batch_index = []
	atom_number= []
	resname = []
	chain = []
	at_name = []
	res_num = []
	alternative_conf = []

	coordsIndexingAtom = []
	partnersTens = []
	atom_counter = 0

	partners = hashings.partners

	maxchain = 0
	maxres = 0
	atnToAtIndex ={}
	chain_hashingTot=[]
	chain_hashing = {}
	for protN, prot in enumerate(atName):
		for atn, atom in enumerate(prot):
			if not atom.split("_")[3] in chain_hashingTot:
				chain_hashingTot+=[atom.split("_")[3]]

	chain_hashingTot = sorted(chain_hashingTot)
	for i in range(len(chain_hashingTot)):
		chain_hashing[chain_hashingTot[i]] = i

	for protN, prot in enumerate(atName):
		chainhashNum = 0

		atnToAtIndex[protN] = {}
		indices_dict = {}

		torsionAngle_hashings = {
									"phi":0,
									"psi":1,
									"omega":2,
									"chi1":3,
									"chi2":4,
									"chi3":5,
									"chi4":6,
									"chi5":7
								}


		for atn,atom in enumerate(prot):
			r,resind,a,chainInd,mut,alt_conf = atom.split("_")

			resind = int(resind)

			if not r in hashings.resi_hash or not a in hashings.atom_hash[r]:
				continue

			res = hashings.resi_hash[r]

			at = hashings.atom_hash[r][a]
			at_name += [at]
			batch_index+=[protN]
			res_num += [resind]


			if not (chain_hashing[chainInd],resind) in indices_dict:
				indices_dict[(chain_hashing[chainInd],resind)]={}


			alternative += [int(mut)]
			alternative_conf += [int(alt_conf)]

			chain += [chain_hashing[chainInd]]
			resname += [res]
			atom_number+=[atom_counter]
			atnToAtIndex[protN][atn] = atom_counter
			current_atom_number = atom_counter

			maxchain = max(maxchain,chain[-1])

			maxres = max(maxres,res_num[-1])
			atom_counter += 1
			if atom_counter<0:
				raise("Overflow of the atom indexing! reduce the number of proteins in the batch or ask Gabriele to implement stuff in 64 bits (even if it is gonna slow down everything and Gabriele does not like it)")

			indices_dict[(chain_hashing[chainInd],resind)][a] = atn

			if not chain_hashing[chainInd] in terminals:
				terminals[chain_hashing[chainInd]] = [current_atom_number,current_atom_number,resind,resind]
			else:
				if  resind >= terminals[chain_hashing[chainInd]][3] and a=="O":
					terminals[chain_hashing[chainInd]][3] = resind
					terminals[chain_hashing[chainInd]][1] = current_atom_number
				elif resind <= terminals[chain_hashing[chainInd]][2] and a == "N":
					terminals[chain_hashing[chainInd]][2] = resind
					terminals[chain_hashing[chainInd]][0] = current_atom_number
			coordsIndexingAtom += [atn]

		### partner calculation! it can probably be implemented in tensorial form ###
		for ch in terminals.keys():
			tN_atom = terminals[ch][0]
			OXT_atom = terminals[ch][1]

			at_name[tN_atom] = hashings.atom_hash[hashings.resi_hash_inverse[resname[tN_atom]]]["tN"]
			at_name[OXT_atom] = hashings.atom_hash[hashings.resi_hash_inverse[resname[OXT_atom]]]["OXT"]

		for atn,atom in enumerate(prot):

			r,resind,a,chainInd,mut,alt_conf = atom.split("_")
			resind = int(resind)

			if not r in hashings.resi_hash or not a in hashings.atom_hash[r]:
				continue

			if not a in partners[r]:

				partnersTens += [[PADDING_INDEX, PADDING_INDEX]]
				continue

			partResNum1 = resind + partners[r][a][0][1]

			partAtom1 = partners[r][a][0][0]

			partResNum2 = resind + partners[r][a][1][1]
			partAtom2 = partners[r][a][1][0]

			if (chain_hashing[chainInd], partResNum1) in indices_dict and partAtom1 in indices_dict[(chain_hashing[chainInd], partResNum1)]:
				part1Ind = indices_dict[(chain_hashing[chainInd], partResNum1)][partAtom1]
			else:
				part1Ind = PADDING_INDEX

			if (chain_hashing[chainInd], partResNum2) in indices_dict and partAtom2 in indices_dict[(chain_hashing[chainInd], partResNum2)]:
				part2Ind = indices_dict[(chain_hashing[chainInd], partResNum2)][partAtom2]
			else:
				part2Ind = PADDING_INDEX

			partnersTens+=[[part1Ind,part2Ind]]

		chain_hashingTot += [chain_hashing]
			##### ANGLE STUFF ####

			# backbone #
			#phi 3,0,1,2
	if verbose:
		timeOld = time.time()

	if verbose:
		print("angle index time calculation: ",time.time()-timeOld)
	atom_description = torch.cat(
								[torch.tensor(batch_index,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(chain,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(res_num,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(resname,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(at_name,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(alternative,dtype=int_dtype).unsqueeze(-1),
								torch.tensor(alternative_conf,dtype=int_dtype).unsqueeze(-1),
								 ],dim=-1).to(device)

	alternative_mask = build_altern_mask(atom_description)
	coordsIndexingAtom = torch.tensor(coordsIndexingAtom,dtype=torch.long,device=device)
	partnersIndexingAtom = torch.tensor(partnersTens,dtype=torch.long,device=device)
	atom_number = torch.tensor(atom_number, dtype=torch.long,device=device)

	angle_indices = torch.full((len(atName),maxchain +1 ,maxres+1,alternative_mask.shape[-1], 8, 4),int(PADDING_INDEX),device=device)


	for protN, prot in enumerate(atName):

		for atn, atom in enumerate(prot):
			for alt in range(alternative_mask.shape[-1]):
				r, resind, a, chainInd, mut, alt_conf = atom.split("_")
				resind = int(resind)

				if not r in hashings.resi_hash or not a in hashings.atom_hash[r]:
					continue

				if not alternative_mask[atnToAtIndex[protN][atn],alt]:

					continue




					#asd
				#print("qui atn è deve tenere conto del fatto che potrebbero mancare robe del padding")
				if a=="N":
					#if protN==1 and chain_hashing[chainInd]==1 and resind==4:
					#	asd

					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["phi"], 1] = atnToAtIndex[protN][atn]

					angle_indices[protN, chain_hashing[chainInd], resind, alt,torsionAngle_hashings["psi"], 0] = atnToAtIndex[protN][atn]
					if resind-1>=0:#resind+1<angle_indices.shape[2]:

						angle_indices[protN,chain_hashing[chainInd], resind-1, alt, torsionAngle_hashings["psi"], 3] = atnToAtIndex[protN][atn]

						angle_indices[protN,chain_hashing[chainInd], resind-1, alt, torsionAngle_hashings["omega"], 2] = atnToAtIndex[protN][atn]


				elif a=="C":
					if resind+1<angle_indices.shape[2]:#resind-1>=0:

						angle_indices[protN,chain_hashing[chainInd], resind+1,alt, torsionAngle_hashings["phi"], 0] = atnToAtIndex[protN][atn]

					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["phi"], 3] = atnToAtIndex[protN][atn]
					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["psi"], 2] = atnToAtIndex[protN][atn]

					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["omega"], 1] =atnToAtIndex[protN][atn]


				elif a=="CA":
					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["phi"], 2] = atnToAtIndex[protN][atn]
					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["psi"], 1] = atnToAtIndex[protN][atn]

					angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings["omega"], 0] = atnToAtIndex[protN][atn]

					if resind-1>=0:#resind+1<angle_indices.shape[2]:

						angle_indices[protN, chain_hashing[chainInd], resind-1,alt, torsionAngle_hashings["omega"], 3] = atnToAtIndex[protN][atn]

				if r in hashings.chi_angle_builder_hashing and a in hashings.chi_angle_builder_hashing[r]:
					for sideAngle in hashings.chi_angle_builder_hashing[r][a]:
						angle_indices[protN, chain_hashing[chainInd], resind,alt, torsionAngle_hashings[sideAngle[0]], sideAngle[1]] = atnToAtIndex[protN][atn]

	#########################################

	#secStruct = generate_sec_struct(coords=)
	#getAngleIndices(atom_description)
	return atom_number,atom_description,coordsIndexingAtom,partnersIndexingAtom,angle_indices,alternative_mask

def getHelicalDipoles(atomPairs, hbondNet, atom_description, alternMask):

	print("WARNING: helix dipole is under developement. dipole is set at the C and N of the helix extremities, but it should be in the middle of the last turn")
	#a = atom_description[:, hashings.atom_description_hash["at_name"]]
	initial_atomPairs_shape = atomPairs.shape

	OatomsList = []
	NatomsList = []
	CatomsList = []
	for res in hashings.atom_hash.keys():
		OatomsList += [hashings.atom_hash[res]["O"]]
		OatomsList += [hashings.atom_hash[res]["OXT"]]
		NatomsList += [hashings.atom_hash[res]["N"]]
		NatomsList += [hashings.atom_hash[res]["tN"]]
		CatomsList += [hashings.atom_hash[res]["C"]]

	a = atom_description[:, hashings.atom_description_hash["at_name"]]
	for i, oatom in enumerate(OatomsList):
		if i == 0:
			oatom_mask = a.eq(oatom)
		else:
			oatom_mask += a.eq(oatom)

	for i, natom in enumerate(NatomsList):
		if i == 0:
			natom_mask = a.eq(natom)
		else:
			natom_mask += a.eq(natom)

	for i, catom in enumerate(CatomsList):
		if i == 0:
			catom_mask = a.eq(catom)
		else:
			catom_mask += a.eq(catom)


	maxRes = atom_description[:,hashings.atom_description_hash["resnum"]].max()+1
	maxBatch = atom_description[:,hashings.atom_description_hash["batch"]].max()+1
	maxChain = atom_description[:, hashings.atom_description_hash["chain"]].max() +1
	nalter = alternMask.shape[-1]

	O_N_atoms = natom_mask[atomPairs[:, 1]] & oatom_mask[atomPairs[:, 0]]
	atomPairsON = atomPairs[O_N_atoms]


	same_chain_mask = atom_description[atomPairsON[:,0],hashings.atom_description_hash["chain"]] == atom_description[atomPairsON[:,1],hashings.atom_description_hash["chain"]]
	helix3 =  torch.abs(atom_description[atomPairsON[:,0],hashings.atom_description_hash["resnum"]] -  atom_description[atomPairsON[:,1],hashings.atom_description_hash["resnum"]]).eq(3) & same_chain_mask & hbondNet[O_N_atoms]
	helix4 =  torch.abs(atom_description[atomPairsON[:,0],hashings.atom_description_hash["resnum"]] -  atom_description[atomPairsON[:,1],hashings.atom_description_hash["resnum"]]).eq(4) & same_chain_mask & hbondNet[O_N_atoms]
	helix5 =  torch.abs(atom_description[atomPairsON[:,0],hashings.atom_description_hash["resnum"]] -  atom_description[atomPairsON[:,1],hashings.atom_description_hash["resnum"]]).eq(5) & same_chain_mask & hbondNet[O_N_atoms]


	## forse è ogni 4 !!!as

	#atomPairs = atomPairs[residue_contiguous]
	helix = torch.zeros((maxBatch,maxChain,maxRes,nalter),dtype=torch.bool,device=atom_description.device)


	for alt in range(alternMask.shape[1]):

		atom_pairsH = atomPairsON[helix3]
		mask = alternMask[atom_pairsH[:,0],alt] & alternMask[atom_pairsH[:,1],alt]
		batchAlt = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["batch"]][mask].long()
		chainAlt = atom_description[atom_pairsH[:, 0], hashings.atom_description_hash["chain"]][mask].long()
		res1 = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["resnum"]][mask].long()
		res2 = atom_description[atom_pairsH[:, 1], hashings.atom_description_hash["resnum"]][mask].long()
		alt_index = torch.full(res1.shape,int(alt),device=res1.device)
		helix.index_put_((batchAlt,chainAlt,res1,alt_index),torch.ones(res1.shape,dtype=torch.bool,device=res1.device))
		helix.index_put_((batchAlt,chainAlt,res2,alt_index),torch.ones(res2.shape,dtype=torch.bool,device=res1.device))



		atom_pairsH = atomPairsON[helix4]
		mask = alternMask[atom_pairsH[:,0],alt] & alternMask[atom_pairsH[:,1],alt]
		batchAlt = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["batch"]][mask].long()
		chainAlt = atom_description[atom_pairsH[:, 0], hashings.atom_description_hash["chain"]][mask].long()
		res1 = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["resnum"]][mask].long()
		res2 = atom_description[atom_pairsH[:, 1], hashings.atom_description_hash["resnum"]][mask].long()
		alt_index = torch.full(res1.shape,int(alt),device=res1.device)
		helix.index_put_((batchAlt,chainAlt,res1,alt_index),torch.ones(res1.shape,dtype=torch.bool,device=res1.device))
		helix.index_put_((batchAlt,chainAlt,res2,alt_index),torch.ones(res2.shape,dtype=torch.bool,device=res1.device))




		atom_pairsH = atomPairsON[helix5]
		mask = alternMask[atom_pairsH[:,0],alt] & alternMask[atom_pairsH[:,1],alt]
		batchAlt = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["batch"]][mask].long()
		chainAlt = atom_description[atom_pairsH[:, 0], hashings.atom_description_hash["chain"]][mask].long()
		res1 = atom_description[atom_pairsH[:,0],hashings.atom_description_hash["resnum"]][mask].long()
		res2 = atom_description[atom_pairsH[:, 1], hashings.atom_description_hash["resnum"]][mask].long()
		alt_index = torch.full(res1.shape,int(alt),device=res1.device)
		helix.index_put_((batchAlt,chainAlt,res1,alt_index),torch.ones(res1.shape,dtype=torch.bool,device=res1.device))
		helix.index_put_((batchAlt,chainAlt,res2,alt_index),torch.ones(res2.shape,dtype=torch.bool,device=res1.device))


	### remove helices of size 2 and 1
	resnum = torch.arange(0,helix.shape[2]-3,device= helix.device)
	resnum1 = resnum + 1
	resnum2 = resnum + 2
	resnum3 = resnum + 3

	bad_two_helices = (~helix[:,:,resnum]) & (helix[:,:,resnum1]) & (helix[:,:,resnum2]) & (~helix[:,:,resnum3])
	#helix[:, :, resnum1][bad_two_helices] = False
	#helix[:, :, resnum2][bad_two_helices] = False

	helix[:, :, resnum1] *= ~bad_two_helices
	helix[:, :, resnum2] *= ~bad_two_helices

	resnum = torch.arange(0,helix.shape[2]-2,device= helix.device)
	resnum1 = resnum + 1
	resnum2 = resnum + 2
	bad_singleHelix = (~helix[:,:,resnum]) & (helix[:,:,resnum1]) & (~helix[:,:,resnum2])
	helix[:, :, resnum1] *= ~bad_singleHelix#] = False

	resnum = torch.arange(0,helix.shape[2]-1,device= helix.device)
	resnum1 = resnum + 1
	helixDipoledC = torch.cat([((helix[:,:,resnum]) & (~helix[:,:,resnum1])), (helix[:,:,[helix.shape[2]-1]])],dim = 2)

	resnum = torch.arange(1,helix.shape[2],device= helix.device)
	resnum1 = resnum - 1

	helixDipoledN =  torch.cat([(helix[:,:,[0]]),((~helix[:,:,resnum1]) & (helix[:,:,resnum]))],dim=2)

	#Ncp = atom_numbers[helixDipoledN]
	#Ccp = atom_numbers[helixDipoledC]
	#removing last
	NcpMask = torch.zeros(initial_atomPairs_shape,dtype=torch.bool,device=atomPairs.device)
	CcpMask = torch.zeros(initial_atomPairs_shape,dtype=torch.bool,device=atomPairs.device)


	resnum = atom_description[:,hashings.atom_description_hash["resnum"]].long()
	batch = atom_description[:,hashings.atom_description_hash["batch"]].long()
	chain = atom_description[:, hashings.atom_description_hash["chain"]].long()

	for alt in range(alternMask.shape[1]):

		mask = alternMask[atomPairs[:,0],alt] & alternMask[atomPairs[:,1],alt]
		resnumAlt = resnum[atomPairs[mask]]
		batchAlt = batch[atomPairs[mask]]
		chainAlt = chain[atomPairs[mask]]

		n_mask1 = mask & (natom_mask[atomPairs[:,0]]) & (helixDipoledN[batchAlt[:,0],chainAlt[:,0],resnumAlt[:,0],alt])
		n_mask2 = mask & (natom_mask[atomPairs[:,1]]) & (helixDipoledN[batchAlt[:,0],chainAlt[:,1],resnumAlt[:,1],alt])

		NcpMask[n_mask1, 0] = True
		NcpMask[n_mask2, 1] = True

		c_mask1 = mask & (catom_mask[atomPairs[:,0]]) & (helixDipoledC[batchAlt[:,0],chainAlt[:,0],resnumAlt[:,0],alt])
		c_mask2 = mask & (catom_mask[atomPairs[:,1]]) & (helixDipoledC[batchAlt[:,0],chainAlt[:,1],resnumAlt[:,1],alt])

		CcpMask[c_mask1, 0] = True
		CcpMask[c_mask2, 1] = True
	return NcpMask, CcpMask

def build_altern_mask(atom_description):

	alternIndices = atom_description[:, hashings.atom_description_hash["alternative"]].long()

	if alternIndices.max()==0:
		return 	torch.ones(atom_description.shape[0],device=atom_description.device,dtype=torch.bool).unsqueeze(-1)

	batch = atom_description[:, hashings.atom_description_hash["batch"]].long()
	chain = atom_description[:, hashings.atom_description_hash["chain"]].long()
	resnum = atom_description[:, hashings.atom_description_hash["resnum"]].long()

	#### mapping to residues ###

	max_altern = alternIndices.max() +1

	#for i in range(max_altern):

	#torch.zeros(alternIndices.shape[0],max_altern,device=alternIndices.device,dtype=torch.bool)

	mask_wild_type = alternIndices.eq(0)
	altern_mask = [ mask_wild_type.unsqueeze(-1)]

	for alt in range(1,max_altern):
		mask = alternIndices.eq(alt)

		alt_mask = mask_wild_type.clone()

		removeIndices = torch.cat([batch[mask].unsqueeze(-1),
								   chain[mask].unsqueeze(-1),
								   resnum[mask].unsqueeze(-1)
								   ],dim=-1).unique(dim=0)

		for mut in removeIndices:
			m = batch.eq(mut[0]) & chain.eq(mut[1]) & resnum.eq(mut[2]) & (alternIndices != alt)
			alt_mask[m] = False

		alt_mask[mask] = True
		altern_mask+=[alt_mask.unsqueeze(-1)]
	altern_mask = torch.cat(altern_mask,dim=-1)
	return altern_mask
