#!/usr/bin/env python
# -*- coding: utf-8 -*-

import torch,time,sys,os
from madrax.sources import fakeAtomsGeneration

from madrax.energies import BondLenConstrain

from madrax.energies import HBond_net as HBond_net
from madrax.energies import Electro_net as Electro_net

from madrax.energies import Clash_net

from madrax.energies import Disulfide_net
from madrax.energies import Vdw
from madrax.energies import Solvatation
from madrax.energies import EntropySC
from madrax.energies import AngleScorer

from madrax.energies import SolventAccessibility

import madrax.sources.hashings	as hashings
import madrax.sources.math_utils as math_utils
from madrax.sources.globalVariables import *


def score_min(x, dim, score): #gives you the tensor with the dim dimension as min
    _tmp=[1]*len(x.size())
    _tmp[dim] = x.size(dim)
    return torch.gather(x,dim,score.min(dim)[1].unsqueeze(dim).repeat(tuple(_tmp))).select(dim,0)

class ForceField(torch.nn.Module):


	def __init__(self,device='cpu'):
		"""
		It initializes the force field object

		Parameters
		----------
		device : str
			device on which the object should be created. Example: "cpu" for CPU calculation, "cuda" for GPU calculation, "cuda:0" to use GPU 0. Refer to pytorch devices semantic for more information
			default = 'cpu'
		Returns
		-------

		"""

		super(ForceField, self).__init__()
		hashings.atom_Properties=hashings.atom_Properties.to(device)
		hashings.fake_atom_Properties = hashings.fake_atom_Properties.to(device)
		self.device = device

		hbond_donors_sidechain = {
			"ARG": ["NE", "NH1", "NH2"],
			"ASN": ["ND2"],
			"GLN": ["NE2"],
			"HIS": ["ND1", "NE2", "ND1H1S", "NE2H2S"],
			"LYS": ["NZ"],
			"SER": ["OG"],
			"THR": ["OG1"],
			"TRP": ["NE1"],
			"TYR": ["OH"],
			"CYS": ["SG"]
		}
		hbond_acceptor_sidechain = {
			"ASN": ["OD1"],
			"ASP": ["OD1", "OD2"],
			"GLN": ["OE1"],
			"GLU": ["OE1", "OE2"],
			"HIS": ["ND1H2S", "NE2H1S"],
			"SER": ["OG"],
			"THR": ["OG1"],
			"TYR": ["OH"],
			"CYS": ["SG"]
		}
		hbond_aromatic = {
			"HIS": ["ND1", "NE2"],
			"PHE": ["RC"],
			"TYR": ["RC"]
			# "TRP":["NE1"]
		}

		backbone_atoms = ["N","C","O","CA","tN","OXT"]

		hbond_donor = set()
		for res in hbond_donors_sidechain.keys():
				for at in hbond_donors_sidechain[res]:
					hbond_donor.add(hashings.atom_hash[res][at])

		for res in hashings.resi_hash.keys():
			if res != "PRO":
				hbond_donor.add(hashings.atom_hash[res]["N"])
				hbond_donor.add(hashings.atom_hash[res]["tN"])

		hbond_acceptor = set()
		for res in hbond_acceptor_sidechain.keys():
			for at in hbond_acceptor_sidechain[res]:
				hbond_acceptor.add(hashings.atom_hash[res][at])

		for res in hashings.resi_hash.keys():
			hbond_acceptor.add(hashings.atom_hash[res]["O"])
			#hbond_acceptor.add(hashings.atom_hash[res]["OXT"])


		hbond_ar = set()
		for res in hbond_aromatic.keys():
			for at in hbond_aromatic[res]:
				hbond_ar.add(hashings.atom_hash[res][at])

		bb_atoms = set()
		for res in hashings.resi_hash.keys():
			for at in backbone_atoms:
				bb_atoms.add(hashings.atom_hash[res][at])

		self.donor = list(hbond_donor)
		self.acceptor = list(hbond_acceptor)
		self.aromatic = list(hbond_ar)
		self.backbone = list(bb_atoms)

		#### list of energies ####
		self.solventAccessibiliry = SolventAccessibility.SolventAccessibility(dev=device,donor=self.donor, acceptor = self.acceptor, backbone_atoms = self.backbone)
		self.hbond_net = HBond_net.HBond_net(dev = device,donor=self.donor, acceptor = self.acceptor, backbone_atoms = self.backbone, hbond_ar=self.aromatic)
		self.disulfide_net= Disulfide_net.Disulfide_net(dev = device)
		self.electro_net = Electro_net.Electro_net(dev=device,backbone_atoms = self.backbone)
		self.clash_net = Clash_net.Clash_net(dev = device, backbone_atoms = self.backbone, donor=self.donor, acceptor = self.acceptor)
		self.vdw = Vdw.Vdw(dev = device, backbone_atoms = self.backbone)
		self.solvatation = Solvatation.Solvatation(dev = device, backbone_atoms = self.backbone, acceptor = self.acceptor, donor = self.donor)
		self.EntropySC=EntropySC.EntropySC(dev = device)
		self.angleScorer = AngleScorer.AngleScorer(dev=device)
		self.bondLenConstraint = BondLenConstrain.BondLenConstrain(dev=device)
		self.load_state_dict(torch.load(os.path.realpath(os.path.dirname(__file__))+"/parameters/final_model.weights",map_location=device))
		
	def getPairwiseRepresentation(self,coords,atom_number,atom_description,distance_threshold = 7):
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

	def forward(self,coordinates,info_tensors,verbose=False,return_network=False):

		"""
		function that calculates the actual energy of the protein(s) or complex(es)

		Parameters
		----------
		coordinates : torch.Tensor
			shape: (Batch, nAtoms, 3)
			coordinates of the proteins. It can be generated using the madrax.utils.parsePDB function
		info_tensors : tuple
			a set of precalculated information tensors required by the forcefield. It can be created out of the box with the function madrax.dataStructures.create_info_tensors starting from the atom names (that can be obtained, along with the coordinates, by the  madrax.utils.parsePDB function)
		verbose : bool
			if you wanna see a lot of text everywhere in your terminal
		return_network : bool
			it return the interacting atoms. Used to build graphs

		Returns
		-------
		final_energy : torch.Tensor
		shape: shape (Batch, nChains, nResi, nMutants, 10)
		The Gibbs energy of the input proteins. The dimensions are organized as follow:

		The batch dimension refers to the protein number (same order of the one defined by the madrax.utils.parsePDB finction), chain refers to the chain index (sorted alphabetically), residue number refers to the residue position, nMutants refers to the mutants you might have implemented in the calculation (this dimension is 1 if no mutants have been added). The last dimension refers to the different types of energy:

		0: Disulfide bonds Energy
		1: Hydrogen Bonds Energy
		2: Electrostatics Energy
		3: Van der Waals Clashes
		4: Polar Solvation Energy
		5: Hydrophobic Solvation Energy
		6: Van der Waals Energy
		7: Backbone Entropy
		8: Side Chain Entropy
		9: Peptide Bond Violations
		10: Rotamer Violation
		"""

		atom_number, atom_description, coordsIndexingAtom, partnersIndexingAtom, angle_indices,alternativeMask = info_tensors

		partnersFinal1 = torch.full((atom_description.shape[0],3),float(PADDING_INDEX),device=self.device)
		partnersFinal2 = torch.full((atom_description.shape[0],3),float(PADDING_INDEX),device=self.device)

		padPart1 = partnersIndexingAtom[..., 0] != PADDING_INDEX
		coordsP1 = coordinates[atom_description[:, 0].long()[padPart1], partnersIndexingAtom[..., 0][padPart1]]
		partnersFinal1[padPart1] = coordsP1

		padPart2 = partnersIndexingAtom[..., 1] != PADDING_INDEX
		coordsP2 = coordinates[atom_description[:, 0].long()[padPart2], partnersIndexingAtom[..., 1][padPart2]]
		partnersFinal2[padPart2] = coordsP2

		del coordsP1,coordsP2,padPart1,padPart2

		coords = coordinates[atom_description[:, 0].long(), coordsIndexingAtom]

		partnersFinal = torch.cat([partnersFinal1.unsqueeze(1), partnersFinal2.unsqueeze(1)], dim=1)
		del coordsIndexingAtom
		#########################################

		fakeAtoms = fakeAtomsGeneration.generateFakeAtomTensor(coords, partnersFinal, atom_description,hashings.fake_atom_Properties.to(self.device))

		timeOld = time.time()
		independent_groups = atom_description[:,hashings.atom_description_hash["batch"]]
		atomPairs = []

		for batch in range(independent_groups.max()+1):

			independentMask = independent_groups.eq(batch)
			atomPairs+=[self.getPairwiseRepresentation( coords, atom_number[independentMask], atom_description)]
			torch.cuda.empty_cache()

		atomPairs = torch.cat(atomPairs,dim=0)

		if verbose:
			print("atom pairs calculation time for",coordinates.shape[0],"proteins", time.time()-timeOld)

		timeOld = time.time()

		if verbose:
			print("pairwise Dist calculation time for", coordinates.shape[0], "proteins", time.time() - timeOld)

		#### torsion angles calculation ####
		existentAnglesMask = (angle_indices!=PADDING_INDEX).prod(-1).bool()
		flat_angle_indices = angle_indices[existentAnglesMask]
		flat_angles,_ = math_utils.dihedral2dVectors(coords[flat_angle_indices[:,0]],coords[flat_angle_indices[:,1]],coords[flat_angle_indices[:,2]],coords[flat_angle_indices[:,3]])
		angles = torch.full(existentAnglesMask.shape,float(PADDING_INDEX),device=self.device)
		angles[existentAnglesMask] = flat_angles.squeeze(-1)


		#####################################
		timeOld = time.time()
		contRat, facc, contRatPol, contRatMC, contRatSC = self.solventAccessibiliry(coords,
																					atom_description,
																					atom_number,
																					atomPairs,
																					fakeAtoms,
																					alternativeMask,
																					)

		torch.cuda.empty_cache()
		if verbose:
			print("Solvent accessibility time:",time.time()-timeOld)

		timeOld = time.time()
		residueEnergyDisulf,atomDISULF,disulfideNetwork = self.disulfide_net(coords,
										atom_description,
										atom_number,
										atomPairs,
										alternativeMask,
										partnersFinal,
										facc
										)
		torch.cuda.empty_cache()
		if verbose:
			print("disulfide time:",time.time()-timeOld)

		timeOld = time.time()
		hbondMC,hbondSC,atomHB,hbondNetwork = self.hbond_net( coords,
																			atom_description,
																			atom_number,
																			atomPairs,
																			fakeAtoms,
																			alternativeMask,
																			disulfideNetwork,
																			partnersFinal,
																			facc
																			)

		torch.cuda.empty_cache()

		if verbose:
			print("hbond time:",time.time()-timeOld)

		timeOld = time.time()
		residueElectroMC,residueElectroSC,atomEnergyElectro,electroNetwork  =self.electro_net( coords,
															partnersFinal,
															atom_description,
															atom_number,
															atomPairs,
															alternativeMask,
															hbondNetwork,
															facc
															)
		torch.cuda.empty_cache()
		if verbose:
			print("electostatics time:",time.time()-timeOld)

		timeOld = time.time()
		residueClash, atomEnergy, clash_mask = self.clash_net(coords,
															atom_description,
															atom_number,
															atomPairs,
															alternativeMask,
															facc,
															hbondNetwork,
															disulfideNetwork
															)
		if verbose:
			print("clashes time:",time.time()-timeOld)

		tot = residueElectroSC + hbondSC + residueClash

		residueHbondSC = score_min(hbondSC, -1, tot)
		residueClash = score_min(residueClash, -1, tot)
		residueElectroSC = score_min(residueElectroSC, -1, tot)

		timeOld = time.time()

		vdwEnergyMC, vdwEnergySC = self.vdw( coords,
											atom_description,
											alternativeMask,
											facc)
		if verbose:
			print("VdW time:",time.time()-timeOld)

		timeOld = time.time()
		solvP, solvH = self.solvatation(atom_description, facc, contRat, contRatPol, atomHB,
										atomDISULF)
		if verbose:
			print("solvation time:",time.time()-timeOld)

		timeOld = time.time()

		entropySC = self.EntropySC(atom_description, contRatSC,residueHbondSC, vdwEnergySC, residueElectroSC,alternativeMask)
		if verbose:
			print("Entropy Side Chain time:", time.time() - timeOld)

		timeOld = time.time()
		peptideBondConstraints = self.bondLenConstraint(atom_description,coords,alternativeMask)
		if verbose:
			print("bond constrain time:", time.time() - timeOld)

		timeOld = time.time()
		entropyMC,rotamerViolation = self.angleScorer(atom_description, angles, alternativeMask)
		if verbose:
			print("angle_scorer:", time.time() - timeOld)

		hbonds = residueHbondSC + hbondMC
		electro = residueElectroMC + residueElectroSC
		VdW = vdwEnergyMC + vdwEnergySC

		if not return_network:
			final_energy = torch.cat([residueEnergyDisulf.unsqueeze(-1),
									  hbonds.unsqueeze(-1),
									  electro.unsqueeze(-1),
									  residueClash.unsqueeze(-1),
									  solvP.unsqueeze(-1),
									  solvH.unsqueeze(-1),
									  VdW.unsqueeze(-1),
									  entropyMC.unsqueeze(-1),
									  entropySC.unsqueeze(-1),
									  peptideBondConstraints.unsqueeze(-1),
									  rotamerViolation.unsqueeze(-1)
									  ],dim=-1)
			return final_energy

		else:
			mask = hbondNetwork
			actual_pairs = atomPairs[mask]
			return atom_description[actual_pairs,:5]