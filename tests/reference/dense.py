"""Reference implementations copied verbatim from upstream master (7aaee04).

Used only as oracles by the differential tests; never imported by the library.
"""
import torch
from madrax.sources import hashings
from madrax.sources.globalVariables import PADDING_INDEX


def dense_pairwise_representation(self,coords,atom_number,atom_description,distance_threshold = 7):
	L = atom_number.shape[0]
	atom_number_pw1 = atom_number.unsqueeze(1).expand(L, L).reshape(-1)
	atom_number_pw2 = atom_number.unsqueeze(0).expand(L, L).reshape(-1)

	triMask = torch.triu(torch.ones((L, L), dtype=torch.bool, device=self.device), diagonal=1).reshape(-1)
	padding_mask = ~(coords[atom_number_pw1 , 0].eq(PADDING_INDEX) | coords[atom_number_pw2 , 0].eq(PADDING_INDEX))
	#print("implementa la alt mask")
	altmask = True

	atom_number_pw1 = atom_number_pw1[padding_mask]
	atom_number_pw2 = atom_number_pw2[padding_mask]
	triMask = triMask[padding_mask]
	#full_mask = padding_mask
	pairwise_atom_name1 = atom_description[:,hashings.atom_description_hash["at_name"]][atom_number_pw1.long()]
	pairwise_atom_name2 = atom_description[:,hashings.atom_description_hash["at_name"]][atom_number_pw2.long()]

	for k, donor_atom in enumerate(self.donor):
		if k == 0:
			donor_mask = pairwise_atom_name2.eq(donor_atom)
			donor_maskOther = pairwise_atom_name1.eq(donor_atom)
		else:
			donor_mask += pairwise_atom_name2.eq(donor_atom)
			donor_maskOther += pairwise_atom_name1.eq(donor_atom)

	for k, acceptor_atom in enumerate(self.acceptor):
		if k == 0:
			acceptor_mask = pairwise_atom_name1.eq(acceptor_atom)
			acceptor_maskOther = pairwise_atom_name2.eq(acceptor_atom)
		else:
			acceptor_mask += pairwise_atom_name1.eq(acceptor_atom)
			acceptor_maskOther += pairwise_atom_name2.eq(acceptor_atom)

	for k, aromatic_atom in enumerate(self.aromatic):
		if k == 0:
			aromatic_mask = pairwise_atom_name1.eq(aromatic_atom)
			aromatic_maskOther = pairwise_atom_name2.eq(aromatic_atom)

		else:
			aromatic_mask += pairwise_atom_name1.eq(aromatic_atom)
			aromatic_maskOther += pairwise_atom_name2.eq(aromatic_atom)


	donor_acceptor_mask = (donor_mask & acceptor_mask)  # .unsqueeze(3)
	donor_aromatic_mask = (donor_mask & aromatic_mask)  # .unsqueeze(3)
	donor_acceptor_maskOther = (donor_maskOther & acceptor_maskOther)
	donor_aromatic_maskOther = (donor_maskOther & aromatic_maskOther)
	same_atom_mask = ~torch.eye(L, dtype=torch.bool, device=self.device).reshape(-1)[padding_mask]

	hbond_or_triangular_mask = (
									   (triMask & ((~donor_acceptor_mask) & (~donor_aromatic_mask)) & (
												   (~donor_acceptor_maskOther) & (~donor_aromatic_maskOther))) |
									   (triMask & (donor_acceptor_mask | donor_aromatic_mask)) |
									   (~triMask & (
										   ~(donor_acceptor_maskOther | donor_aromatic_maskOther)) & (
													(donor_acceptor_mask) | (donor_aromatic_mask)))
							   ) & same_atom_mask

	atom_number_pw1 = atom_number_pw1[hbond_or_triangular_mask]
	atom_number_pw2 = atom_number_pw2[hbond_or_triangular_mask]
	longMask = torch.nn.functional.pairwise_distance(coords[atom_number_pw1.long()],coords[atom_number_pw2.long()]).le(distance_threshold)
	atom_number_pw1 = atom_number_pw1[longMask]
	atom_number_pw2 = atom_number_pw2[longMask]

	return torch.cat([atom_number_pw1.unsqueeze(-1),atom_number_pw2.unsqueeze(-1)],dim=-1)


def reference_angle_indices(atName, alternative_mask, atom_description, device="cpu"):
	"""master's per-atom loop from create_info_tensors, unchanged apart from being a function."""
	chain_hashingTot = sorted({a.split("_")[3] for prot in atName for a in prot})
	chain_hashing = {c: i for i, c in enumerate(chain_hashingTot)}
	torsionAngle_hashings = {"phi": 0, "psi": 1, "omega": 2, "chi1": 3, "chi2": 4, "chi3": 5, "chi4": 6, "chi5": 7}
	maxchain = int(atom_description[:, hashings.atom_description_hash["chain"]].max())
	maxres = int(atom_description[:, hashings.atom_description_hash["resnum"]].max())
	atnToAtIndex = {}
	k = 0
	for protN, prot in enumerate(atName):
		atnToAtIndex[protN] = {}
		for atn, atom in enumerate(prot):
			r, resind, a = atom.split("_")[:3]
			if not r in hashings.resi_hash or not a in hashings.atom_hash[r]:
				continue
			atnToAtIndex[protN][atn] = k
			k += 1
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

	return angle_indices
