#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  vectors_utils.py
#  
#  Copyright 2020 Gabriele Orlando <orlando.gabriele89@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import torch
from madrax.sources import math_utils
import numpy as np
from madrax.sources import hashings
from madrax.sources.globalVariables import *

def get_interaction_angles(acceptCoords,donorCoords,acceptorPartners1,donorPartners1,acceptorPartners2,donorPartners2,hydrogens, freeOrb):

	
		
	AccProton= acceptCoords-hydrogens
	
	DonProton= donorCoords-hydrogens
	
	#Part1Acceptor= acceptorPartners1 - acceptCoords
	
	#nanmask = torch.norm(AccProton- DonProton,dim=1).le(0.000000001).unsqueeze(1).repeat(1,3)

	freeProt = freeOrb - hydrogens
	AccFree= acceptCoords - freeOrb
	
	# angles involved ##

		
	freeProtDon,nanMaskFpd,_=math_utils.angle2dVectors(freeProt,DonProton)

	protFreeAcc,nanMaskPfa,test=math_utils.angle2dVectors(-freeProt, AccFree)

	dihed,nanMaskDihed,test = math_utils.dihedral2dVectors(donorCoords, hydrogens, freeOrb, acceptCoords,testing=True)
	test = hydrogens
	ang_planes,nanMaskPlane = math_utils.plane_angle(acceptCoords,acceptorPartners1,acceptorPartners2,donorCoords,donorPartners1,donorPartners2)

	#assert np.isnan(torch.sum(protFreeAcc).cpu().data.numpy())==False

	#non usi i plane angles, metti il nan separato
	fullNanMask =  nanMaskFpd * nanMaskPfa * nanMaskDihed

	return torch.cat([freeProtDon, protFreeAcc, dihed, ang_planes],dim=1),fullNanMask.squeeze(1),nanMaskPlane,test


def get_interaction_anglesStatisticalPotential(acceptCoords, donorCoords, acceptorPartners1, donorPartners1):

	DonProton = donorPartners1 - donorCoords

	freeProt = acceptCoords - donorCoords

	AccFree = acceptorPartners1 - acceptCoords

	# angles involved ##

	freeProtDon, nanMaskFpd, _ = math_utils.angle2dVectors(freeProt, DonProton)

	protFreeAcc, nanMaskPfa, test = math_utils.angle2dVectors(-freeProt, AccFree)

	dihed, nanMaskDihed, test = math_utils.dihedral2dVectors(donorPartners1, donorCoords, acceptCoords, acceptorPartners1, testing=True)

	assert np.isnan(torch.sum(protFreeAcc).cpu().data.numpy()) == False

	# non usi i plane angles, metti il nan separato
	fullNanMask = nanMaskFpd * nanMaskPfa * nanMaskDihed

	return torch.cat([freeProtDon, protFreeAcc, dihed], dim=1), fullNanMask.squeeze(1)


def get_standard_angles(donorCoords,acceptCoords,hydrogens,freeOrb):
	AccProton= acceptCoords-hydrogens
	DonProton= donorCoords-hydrogens
	freeProt = freeOrb - hydrogens
	AccFree = acceptCoords - freeOrb

	# angles involved ##
	freeProtDon,nanMaskFpd=math_utils.angle2dVectors(freeProt,DonProton)
	protFreeAcc,nanMaskPfa=math_utils.angle2dVectors(-freeProt, AccFree)
	dihed,nanMaskDihed = math_utils.dihedral2dVectors(donorCoords, hydrogens, freeOrb, acceptCoords)
	fullNanMask =  nanMaskFpd & nanMaskPfa & nanMaskDihed

	return torch.cat([freeProtDon, protFreeAcc, dihed],dim=1),fullNanMask.squeeze(1)

def define_partners_coords(coords,atom_hashing_position,res_num,atom_names,res_names, partners):
	partners_tens1=[]
	partners_tens2=[]
	mask_partners=[]
	batch_len = len(coords)#.shape[0]
	orig_tot=[]
	indexes1_tot = []
	indexes2_tot = []

	for batch in range(batch_len):
		p_tens1=torch.zeros(coords[batch].shape).type(coords[0].type())
		p_tens2=torch.zeros(coords[batch].shape).type(coords[0].type())
		indexes1=[]
		indexes2=[]
		orig=[]
		maxres = max(res_num[batch])
		for i in range(len(coords[batch])): # every atom

			respos=res_num[batch][i]
			atom_pos=i
			atom_n=atom_names[batch][i]
			res_n=res_names[batch][i]
			
			if atom_n in partners[res_n]:
				par1 = partners[res_n][atom_n][0][0]
				residue_correction1 = partners[res_n][atom_n][0][1] # atoms on the peptide bond define the plane using atoms from the preevious/successive residue
				par2 = partners[res_n][atom_n][1][0]
				residue_correction2 = partners[res_n][atom_n][1][1]
				#print(batch)

				if not (respos+residue_correction1<0 or respos+residue_correction1>=maxres or respos+residue_correction2<0 or respos+residue_correction2>=maxres) and par1 in atom_hashing_position[batch][respos+residue_correction1] and par2 in atom_hashing_position[batch][respos+residue_correction2]: # for missing atoms
					indexes1+=[atom_hashing_position[batch][respos+residue_correction1][par1]]
					indexes2+=[atom_hashing_position[batch][respos+residue_correction2][par2]]
					orig+=[i]

		#p_tens1[orig]=coords[batch][indexes1]
		#p_tens2[orig]=coords[batch][indexes2]
		#partners_tens1+=[p_tens1]
		#partners_tens2+=[p_tens2]
		indexes1_tot +=[indexes1]
		indexes2_tot += [indexes2]

		orig_tot += [orig]
		ners_tens = p_tens1.unsqueeze(0)
		mp = torch.zeros(coords[batch].shape[0])
		mp[orig] = 1
		mp = mp.ge(1).unsqueeze(1)
		mask_partners+=[mp]	
		
	return indexes1_tot,indexes2_tot,orig_tot,mask_partners#partners_tens1,partners_tens2,mask_partners

def calculateTorsionAngles(coords_tot,atom_description,angle_coords,float_type=torch.float):

	### bb ###
	c_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["C"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["C"])
	ca_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["CA"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["CA"])
	n_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["N"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["N"])

	batch = coords_tot.shape[0]
	L = coords_tot.shape[1]


	coords = coords_tot[:,:,hashings.property_hashings["coords_total"]["coords"],:]
	resiMinus1 = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]] -1
	resiPlus1 = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]] +1
	resi = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]]
	chain_ind = atom_description[:,:,hashings.property_hashings["atom_description"]["chain"]].long()

	maxResi = resi.max()+1
	max_chain_ind = chain_ind.max() + 1

	batch_ind = torch.arange(0,batch,device=coords_tot.device).unsqueeze(1).expand(batch,L)

	anglesFull = torch.full((batch,max_chain_ind,maxResi,9,3),PADDING_INDEX,device=coords_tot.device,dtype=float_type)

	padding_mask = ~(resi==PADDING_INDEX)

	nMask = n_hahsingMask &  padding_mask
	full_ind = (batch_ind[nMask].long(), chain_ind[nMask],resi[nMask].long(), torch.zeros(batch_ind.shape,device=coords_tot.device,dtype=torch.long)[nMask])
	anglesFull = anglesFull.index_put_(full_ind, coords[nMask])

	caMask = ca_hahsingMask &  padding_mask
	full_ind = (batch_ind[caMask].long(), chain_ind[caMask],resi[caMask].long(), torch.full(batch_ind.shape,1,device=coords_tot.device,dtype=torch.long)[caMask])
	anglesFull = anglesFull.index_put_(full_ind, coords[caMask])

	cMask = c_hahsingMask &  padding_mask
	full_ind = (batch_ind[cMask].long(),chain_ind[cMask], resi[cMask].long(), torch.full(batch_ind.shape,2,device=coords_tot.device,dtype=torch.long)[cMask])
	anglesFull = anglesFull.index_put_(full_ind, coords[cMask])

	c_plus1_Mask = c_hahsingMask &  padding_mask & (resiPlus1 < maxResi)
	full_ind = (batch_ind[c_plus1_Mask].long(), chain_ind[c_plus1_Mask],resiPlus1[c_plus1_Mask].long(), torch.full(c_plus1_Mask.shape,3,device=coords_tot.device,dtype=torch.long)[c_plus1_Mask])
	anglesFull = anglesFull.index_put_(full_ind, coords[c_plus1_Mask])

	n_minus_oneMaks = n_hahsingMask & (resiMinus1 >= 0) &  padding_mask
	full_ind = (batch_ind[n_minus_oneMaks].long(), chain_ind[n_minus_oneMaks], resiMinus1[n_minus_oneMaks].long(), torch.full(batch_ind.shape,4,device=coords_tot.device,dtype=torch.long)[n_minus_oneMaks])
	anglesFull = anglesFull.index_put_(full_ind, coords[n_minus_oneMaks])

	ca_minus_oneMaks = ca_hahsingMask & (resiMinus1 >= 0) &  padding_mask
	full_ind = (batch_ind[ca_minus_oneMaks].long(),chain_ind[ca_minus_oneMaks], resiMinus1[ca_minus_oneMaks].long(), torch.full(batch_ind.shape,5,device=coords_tot.device,dtype=torch.long)[ca_minus_oneMaks])
	anglesFull = anglesFull.index_put_(full_ind, coords[ca_minus_oneMaks])

	residue_padding_mask = ~(anglesFull[:,:,:,:,0] == PADDING_INDEX)

	everythingForPhi= residue_padding_mask[:,:,:,(3,0,1,2)].prod(dim=-1).bool()
	phi_angles= anglesFull[everythingForPhi]
	phi,_ = math_utils.dihedral2dVectors(phi_angles[:,3], phi_angles[:,0], phi_angles[:,1], phi_angles[:,2])

	#Sbagliatolecoordatomi

	everythingForPsi= residue_padding_mask[:,:,:,(0,1,2,4)].prod(dim=-1).bool()
	psi_angles= anglesFull[everythingForPsi]
	psi,_ = math_utils.dihedral2dVectors(psi_angles[:,0], psi_angles[:,1], psi_angles[:,2],psi_angles[:,4])

	everythingForOmega= residue_padding_mask[:,:,:,(1,2,4,5)].prod(dim=-1).bool()
	omega_angles= anglesFull[everythingForOmega]
	omega,_ = math_utils.dihedral2dVectors( omega_angles[:,1], omega_angles[:,2],omega_angles[:,4],omega_angles[:,5])

	#### full angles ####
	phiFull = torch.full((batch,max_chain_ind,maxResi),PADDING_INDEX,device=coords_tot.device)
	psiFull = torch.full((batch,max_chain_ind,maxResi),PADDING_INDEX,device=coords_tot.device)
	omegaFull = torch.full((batch,max_chain_ind,maxResi),PADDING_INDEX,device=coords_tot.device)

	phiFull[everythingForPhi] = phi.squeeze(-1)
	psiFull[everythingForPsi] = psi.squeeze(-1)
	omegaFull[everythingForOmega] = omega.squeeze(-1)


	sidechain_todo = (~(angle_coords[:,:,:,:,0].eq(PADDING_INDEX))).prod(-1).bool() # calculate angles for which all atoms are present
	todo = angle_coords[sidechain_todo]
	sidechain, _ = math_utils.dihedral2dVectors( todo[:,0], todo[:,1],todo[:,2],todo[:,3])

	fullanglesSC = torch.full((batch,max_chain_ind,maxResi,5),PADDING_INDEX).type_as(sidechain)
	fullanglesSC[sidechain_todo] = sidechain.squeeze(-1)
	fullangles = torch.cat([phiFull.unsqueeze(-1),psiFull.unsqueeze(-1),omegaFull.unsqueeze(-1),fullanglesSC],dim=3)

	return fullangles
	asd

def calculateTorsionAnglesBuildModels(coords_tot,atom_description,angle_coords,alternatives,alternative_resi,seqChain,seqNum,float_type=torch.float):

	### bb ###
	c_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["C"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["C"])
	ca_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["CA"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["CA"])
	n_hahsingMask = (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["ALA"]["N"]) | (atom_description[:,:,hashings.property_hashings["atom_description"]["atom_types"]] == hashings.atom_hash["PRO"]["N"])


	batch = coords_tot.shape[0]
	L = coords_tot.shape[1]
	max_alternatives = alternatives.shape[-1]
	nrots = coords_tot.shape[-3]

	coords = coords_tot[:,:,:,hashings.property_hashings["coords_total"]["coords"],:]
	resiMinus1 = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]] -1
	resiPlus1 = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]] +1
	resi = atom_description[:,:,hashings.property_hashings["atom_description"]["res_num"]]
	chain_ind = atom_description[:,:,hashings.property_hashings["atom_description"]["chain"]].long()

	maxResi = resi.max()+1
	max_chain_ind = chain_ind.max() + 1

	batch_ind = torch.arange(0,batch,device=coords_tot.device).unsqueeze(1).expand(batch,L)

	#batch_indResi = torch.arange(0,batch,device=coords_tot.device).unsqueeze(-1).unsqueeze(-1).unsqueeze(-1).expand(batch,max_chain_ind,maxResi,nrots)
	#chain_indResi = torch.arange(0,max_chain_ind,device=coords_tot.device).unsqueeze(0).unsqueeze(-1).unsqueeze(-1).expand(batch,max_chain_ind,maxResi,nrots)
	#resid_indResi = torch.arange(0,maxResi,device=coords_tot.device).unsqueeze(0).unsqueeze(0).unsqueeze(-1).expand(batch,max_chain_ind,maxResi,nrots)
	#rotat_indResi = torch.arange(0,nrots,device=coords_tot.device).unsqueeze(0).unsqueeze(0).unsqueeze(0).expand(batch,max_chain_ind,maxResi,nrots)


	anglesFull = torch.full((batch,max_chain_ind,maxResi,9,nrots,3),PADDING_INDEX,device=coords_tot.device,dtype=float_type)

	fullangles = []
	for alt in range(max_alternatives):
		padding_mask = ~(resi==PADDING_INDEX) & alternatives[:,:,alt]

		camask = atom_description[:, :, hashings.property_hashings["atom_description"]["atom_types"]].eq(
			hashings.atom_hash["ALA"]["CA"]) | atom_description[:, :,
											   hashings.property_hashings["atom_description"]["atom_types"]].eq(
			hashings.atom_hash["PRO"]["CA"])

		batch_indRes = batch_ind[camask & padding_mask].unsqueeze(1).unsqueeze(1).expand(-1,nrots,5)
		resnu_indRes = resi[camask & padding_mask].unsqueeze(1).unsqueeze(1).expand(-1,nrots,5)
		chain_indRes = chain_ind[camask & padding_mask].unsqueeze(1).unsqueeze(1).expand(-1,nrots,5)


		alternative_Resi_mask = alternative_resi[:,:,alt]
		nMask = n_hahsingMask &  padding_mask
		full_ind = (batch_ind[nMask].long(), chain_ind[nMask],resi[nMask].long(), torch.zeros(batch_ind.shape,device=coords_tot.device,dtype=torch.long)[nMask])
		anglesFull = anglesFull.index_put_(full_ind, coords[nMask])

		caMask = ca_hahsingMask &  padding_mask
		full_ind = (batch_ind[caMask].long(), chain_ind[caMask],resi[caMask].long(), torch.full(batch_ind.shape,1,device=coords_tot.device,dtype=torch.long)[caMask])
		anglesFull = anglesFull.index_put_(full_ind, coords[caMask])

		cMask = c_hahsingMask &  padding_mask
		full_ind = (batch_ind[cMask].long(),chain_ind[cMask], resi[cMask].long(), torch.full(batch_ind.shape,2,device=coords_tot.device,dtype=torch.long)[cMask])
		anglesFull = anglesFull.index_put_(full_ind, coords[cMask])

		c_plus1_Mask = c_hahsingMask &  padding_mask & (resiPlus1 < maxResi)
		full_ind = (batch_ind[c_plus1_Mask].long(), chain_ind[c_plus1_Mask],resiPlus1[c_plus1_Mask].long(), torch.full(c_plus1_Mask.shape,3,device=coords_tot.device,dtype=torch.long)[c_plus1_Mask])
		anglesFull = anglesFull.index_put_(full_ind, coords[c_plus1_Mask])

		n_minus_oneMaks = n_hahsingMask & (resiMinus1 >= 0) &  padding_mask
		full_ind = (batch_ind[n_minus_oneMaks].long(), chain_ind[n_minus_oneMaks], resiMinus1[n_minus_oneMaks].long(), torch.full(batch_ind.shape,4,device=coords_tot.device,dtype=torch.long)[n_minus_oneMaks])
		anglesFull = anglesFull.index_put_(full_ind, coords[n_minus_oneMaks])

		ca_minus_oneMaks = ca_hahsingMask & (resiMinus1 >= 0) &  padding_mask
		full_ind = (batch_ind[ca_minus_oneMaks].long(),chain_ind[ca_minus_oneMaks], resiMinus1[ca_minus_oneMaks].long(), torch.full(batch_ind.shape,5,device=coords_tot.device,dtype=torch.long)[ca_minus_oneMaks])
		anglesFull = anglesFull.index_put_(full_ind, coords[ca_minus_oneMaks])

		residue_padding_mask = ~(anglesFull[:,:,:,:,:,0] == PADDING_INDEX)


		everythingForPhi= residue_padding_mask[:,:,:,(3,0,1,2)].prod(dim=-2).bool()
		phi_angles= anglesFull.transpose(-2,-3)[everythingForPhi]
		phi,_ = math_utils.dihedral2dVectors(phi_angles[:,3], phi_angles[:,0], phi_angles[:,1], phi_angles[:,2])

		#Sbagliatolecoordatomi

		everythingForPsi= residue_padding_mask[:,:,:,(0,1,2,4)].prod(dim=-2).bool()
		psi_angles= anglesFull.transpose(-2,-3)[everythingForPsi]
		psi,_ = math_utils.dihedral2dVectors(psi_angles[:,0], psi_angles[:,1], psi_angles[:,2],psi_angles[:,4])

		everythingForOmega= residue_padding_mask[:,:,:,(1,2,4,5)].prod(dim=-2).bool()
		omega_angles= anglesFull.transpose(-2,-3)[everythingForOmega]
		omega,_ = math_utils.dihedral2dVectors( omega_angles[:,1], omega_angles[:,2],omega_angles[:,4],omega_angles[:,5])

		#### full angles ####
		phiFull = torch.full((batch,max_chain_ind,maxResi,nrots),PADDING_INDEX,device=coords_tot.device)
		psiFull = torch.full((batch,max_chain_ind,maxResi,nrots),PADDING_INDEX,device=coords_tot.device)
		omegaFull = torch.full((batch,max_chain_ind,maxResi,nrots),PADDING_INDEX,device=coords_tot.device)

		phiFull[everythingForPhi] = phi.squeeze(-1)
		psiFull[everythingForPsi] = psi.squeeze(-1)
		omegaFull[everythingForOmega] = omega.squeeze(-1)


		alterAngleCoords = angle_coords[alternative_Resi_mask]
		sidechain_todo = (~(alterAngleCoords[:,:,:,:,0].eq(PADDING_INDEX))).prod(-1).bool() # calculate angles for which all atoms are present
		todo = alterAngleCoords[sidechain_todo]
		sidechain, _ = math_utils.dihedral2dVectors( todo[:,0], todo[:,1],todo[:,2],todo[:,3])

		#batch_indRes = batch_indRes[sidechain_todo]
		#chain_indRes = chain_indRes[sidechain_todo]
		#resnu_indRes = resnu_indRes[sidechain_todo]


		seqshape =seqNum.shape[1]

		batch_SC = torch.arange(angle_coords.shape[0]).unsqueeze(-1).unsqueeze(-1).unsqueeze(-1).expand(batch,seqshape, nrots,5)[alternative_Resi_mask][sidechain_todo]
		chain_SC = seqChain.unsqueeze(-1).unsqueeze(-1).expand(batch,seqshape, nrots,5)[alternative_Resi_mask][sidechain_todo]
		resi_SC = seqNum.unsqueeze(-1).unsqueeze(-1).expand(batch,seqshape, nrots,5)[alternative_Resi_mask][sidechain_todo]
		numrot = torch.arange(  nrots                ).unsqueeze(0 ).unsqueeze(0 ).unsqueeze(-1).expand(batch,seqshape, nrots,5)[alternative_Resi_mask][sidechain_todo]
		numtor = torch.arange(  5                    ).unsqueeze(0 ).unsqueeze(0 ).unsqueeze(0 ).expand(batch,seqshape, nrots,5)[alternative_Resi_mask][sidechain_todo]
		indices = (batch_SC, chain_SC, resi_SC, numrot, numtor)


		'''
		numrot = torch.arange(0,nrots,device=coords_tot.device,dtype=torch.long).unsqueeze(0).unsqueeze(-1).expand(batch_indRes.shape[0],nrots,5)[sidechain_todo]
		numtor = torch.arange(0,5,device=coords_tot.device,dtype=torch.long).unsqueeze(0).unsqueeze(0).expand(batch_indRes.shape[0],nrots,5)[sidechain_todo]

		indices =	(batch_indRes[sidechain_todo].long(),
					 chain_indRes[sidechain_todo].long(),
					 resnu_indRes[sidechain_todo].long(),
					 numrot,
					 numtor
					)
		'''
		fullanglesSC = torch.full((batch,max_chain_ind,maxResi,nrots,5),PADDING_INDEX).type_as(sidechain)
		fullanglesSC.index_put_(indices,sidechain.squeeze(-1))

		fullangles += [torch.cat([phiFull.unsqueeze(-1),psiFull.unsqueeze(-1),omegaFull.unsqueeze(-1),fullanglesSC],dim=-1).unsqueeze(-3)]

	fullangles = torch.cat(fullangles,dim=-3)
	return fullangles

def main(args):
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
