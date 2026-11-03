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
		self._donor_t = torch.tensor(self.donor)
		self._accar_t = torch.tensor(self.acceptor + self.aromatic)

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
		
	def getPairwiseRepresentation(self, coords, atom_number, atom_description, distance_threshold=7, block=2048):
		"""Same pairs/orientation/order as the original O(L^2)-memory version (reference copy in tests/reference/dense.py)."""
		atom_number = atom_number[coords[atom_number, 0] != PADDING_INDEX]
		L = atom_number.shape[0]
		if L < 2:
			return torch.zeros((0, 2), dtype=atom_number.dtype, device=coords.device)
		names = atom_description[atom_number.long(), hashings.atom_description_hash["at_name"]]
		donor = torch.isin(names, self._donor_t.to(names.device, names.dtype))
		acc = torch.isin(names, self._accar_t.to(names.device, names.dtype))
		x = coords[atom_number.long()]
		I, J = [], []
		with torch.no_grad():
			for s in range(0, L, block):
				e = min(s + block, L)
				ii, jj = torch.nonzero(torch.cdist(x[s:e], x[s:]) <= distance_threshold + 1e-3, as_tuple=True)
				ii = ii + s; jj = jj + s
				k = jj > ii
				I.append(ii[k]); J.append(jj[k])
		i = torch.cat(I); j = torch.cat(J)
		flip = (acc[j] & donor[i]) & ~(acc[i] & donor[j])
		a = torch.where(flip, j, i); b = torch.where(flip, i, j)
		pa, pb = atom_number[a], atom_number[b]
		ok = torch.nn.functional.pairwise_distance(coords[pa.long()], coords[pb.long()]).le(distance_threshold)
		a, b, pa, pb = a[ok], b[ok], pa[ok], pb[ok]
		order = torch.argsort(a * L + b)
		return torch.stack([pa[order], pb[order]], dim=-1)

	def _bond_len_per_alternative(self, atom_description, coords, alternativeMask):
		"""Peptide-bond term evaluated one alternative at a time.

		The compiled BondLenConstrain mixes alternatives: with a LYS mutant present, the wild-type row of
		that residue takes the mutant's value. On a single alternative it is exact, so each alternative is
		scored as its own structure and the columns are stacked.
		"""
		nalt = alternativeMask.shape[-1]
		if nalt == 1:
			return self.bondLenConstraint(atom_description, coords, alternativeMask)
		dims = [int(atom_description[:, hashings.atom_description_hash[k]].max()) + 1 for k in ("batch", "chain", "resnum")]
		alt_col = hashings.atom_description_hash["alternative"]
		cols = []
		for k in range(nalt):
			m = alternativeMask[:, k]
			sub = atom_description[m].clone()
			sub[:, alt_col] = 0
			e = self.bondLenConstraint(sub, coords[m], torch.ones((int(m.sum()), 1), dtype=torch.bool, device=coords.device))[..., 0]
			if list(e.shape) != dims:
				e = torch.nn.functional.pad(e, [0, dims[2] - e.shape[2], 0, dims[1] - e.shape[1], 0, dims[0] - e.shape[0]])
			cols.append(e)
		return torch.stack(cols, dim=-1)

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
		peptideBondConstraints = self._bond_len_per_alternative(atom_description, coords, alternativeMask)
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