from madrax.sources import geometry
import torch

from madrax.sources import hashings
from madrax.sources.globalVariables import *
def select_closest_hydrogen(x1_coords, hydrogen_p):
	
	dist=torch.norm(torch.add(x1_coords.unsqueeze(1),-hydrogen_p),dim=2)

	closest_H_dist, indices  = torch.min(dist,dim=1)

	mask = torch.arange(dist.size(1)).to(dist.device).reshape(1, -1) == indices.unsqueeze(1)
	mask = mask.unsqueeze(2)#.repeat(1,1,1)
	

	hcoords = torch.masked_select(hydrogen_p,mask).view(-1,3)
	
	return hcoords,closest_H_dist.view(-1,1)

def select_closest_freeOrb(freeorbs_pw, hydrogens_pairwise):
	
	hydrogens_pairwise =hydrogens_pairwise.unsqueeze(1)
	
	dist=torch.norm(torch.add(freeorbs_pw,-hydrogens_pairwise),dim=2)
	
	closest_orb_dist, indices  = torch.min(dist,dim=1)

	mask = torch.arange(dist.size(1)).to(dist.device).reshape(1, -1) == indices.unsqueeze(1)

	focoords = freeorbs_pw[mask]

	return focoords,closest_orb_dist.view(-1,1)	

def mask_angles(angles, nanmask, minAng, maxAng):
	
	freeProtDon_mask = torch.abs(angles[:,0]).ge(minAng[:,0]) * torch.abs(angles[:,0]).le(maxAng[:,0])
	ProtFreeAcc_mask = torch.abs(angles[:,1]).ge(minAng[:,1]) * torch.abs(angles[:,1]).le(maxAng[:,1])
	Dihed_mask = torch.abs(angles[:,2]).ge(minAng[:,2]) 	
	#print(nanmask.shape,freeProtDon_mask.shape,aromatic_mask)
	interaction_angles_fullmask = (freeProtDon_mask & ProtFreeAcc_mask & Dihed_mask)  & nanmask

	return interaction_angles_fullmask
	
def select_allVA(acceptor,donor,hydrogens_o,freeorbs_o):
	#print(donor.shape,hydrogens.shape,freeorbs.shape,aromatic_mask.shape)

	n_elements = freeorbs_o.shape[0]
	n_FO = freeorbs_o.shape[1]
	n_H  = hydrogens_o.shape[1]

	
	donor = donor.unsqueeze(1).repeat(1,n_H*n_FO,1)
	acceptor = acceptor.unsqueeze(1).repeat(1,n_H*n_FO,1)

	freeorbs = freeorbs_o.unsqueeze(2).expand(n_elements, n_FO,n_H, 3).permute(3,0,1,2).reshape(3,n_elements,n_H*n_FO).permute(1,2,0)
	hydrogens = hydrogens_o.unsqueeze(1).expand(n_elements,n_FO,n_H, 3).permute(3,0,1,2).reshape(3,n_elements,n_H*n_FO).permute(1,2,0)

	hydrogens =hydrogens.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	freeorbs =freeorbs.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	donor =donor.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	acceptor =acceptor.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	
	maksMissingH = ~hydrogens.eq(PADDING_INDEX)
	maskMissingFO= ~freeorbs.eq(PADDING_INDEX)

	maskMissing=(maksMissingH & maskMissingFO).sum(-1).type(maksMissingH.type())

	hydrogens = hydrogens[maskMissing]	
	freeorbs = freeorbs[maskMissing]
	donor = donor[maskMissing]
	acceptor = acceptor[maskMissing]
	
	angles, nanmask = geometry.get_standard_angles(donor, acceptor, hydrogens, freeorbs)

	ia_mask = mask_angles(angles, nanmask)

	#hydrogens=hydrogens.transpose(0,1).view(3,n_elements,n_H,n_FO).permute(1,2,3,0)
	#freeorbs = freeorbs.transpose(0,1).view(3,n_elements,n_H,n_FO).permute(1,2,3,0)
	#interaction_angles_fullmask = interaction_angles_fullmask.view(n_elements,n_FO*n_H)
	distH = torch.full((n_elements*n_H*n_FO,3), PADDING_INDEX).type(hydrogens.type())
	distFO = torch.full((n_elements*n_H*n_FO,3), PADDING_INDEX).type(hydrogens.type())
	interaction_angles_fullmask =torch.full([n_elements*n_H*n_FO], False).type(maksMissingH.type())
	
	distH[maskMissing]  = hydrogens
	distFO[maskMissing] = freeorbs
	interaction_angles_fullmask[maskMissing] = ia_mask

	dist = torch.norm(distH-distFO,dim=1)
	dist.masked_fill_(~interaction_angles_fullmask,NON_ACCEPTABLE_DISTANCE)
	dist = dist.view(n_elements,n_H*n_FO)	

	distH = distH.transpose(1,0).view(3,n_elements,n_H*n_FO).permute(1,2,0)
	distFO = distFO.transpose(1,0).view(3,n_elements,n_H*n_FO).permute(1,2,0)

	_, indices  = torch.min(dist,dim=1)
	
	mask_min = torch.arange(dist.size(1)).to(dist.device).reshape(1, -1) == indices.unsqueeze(1)

	hydrogens_o = distH[mask_min]
	freeorbs_o = distFO[mask_min]
	interaction_angles_fullmask= interaction_angles_fullmask.view(n_elements,n_H*n_FO).sum(dim=1).type(mask_min.type())
	
	return hydrogens_o,freeorbs_o,interaction_angles_fullmask
	
def do_pairwise_hbonds(acceptor,donor,hydrogens_o,freeorbs_o,acceptorPartners1,donorPartners1,acceptorPartners2,donorPartners2, minAng, maxAng):
	dev=acceptor.device
	n_elements = freeorbs_o.shape[0]
	n_FO = freeorbs_o.shape[1]
	n_H  = hydrogens_o.shape[1]
	
	donor = donor.unsqueeze(1).repeat(1,n_H*n_FO,1)
	acceptor = acceptor.unsqueeze(1).repeat(1,n_H*n_FO,1)
	acceptorPartners1 = acceptorPartners1.unsqueeze(1).repeat(1,n_H*n_FO,1)
	acceptorPartners2 = acceptorPartners2.unsqueeze(1).repeat(1,n_H*n_FO,1)
	donorPartners1 = donorPartners1.unsqueeze(1).repeat(1,n_H*n_FO,1)
	donorPartners2 = donorPartners2.unsqueeze(1).repeat(1,n_H*n_FO,1)

	#test = torch.load("test2.m").unsqueeze(1).repeat(1,n_H*n_FO)

	
	minAng = minAng.unsqueeze(1).repeat(1,n_H*n_FO,1)
	maxAng = maxAng.unsqueeze(1).repeat(1,n_H*n_FO,1)
	
	freeorbs = freeorbs_o.unsqueeze(2).expand(n_elements, n_FO,n_H, 3).permute(3,0,1,2).reshape(3,n_elements,n_H*n_FO).permute(1,2,0)
	hydrogens = hydrogens_o.unsqueeze(1).expand(n_elements,n_FO,n_H, 3).permute(3,0,1,2).reshape(3,n_elements,n_H*n_FO).permute(1,2,0)

	hydrogens =hydrogens.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	freeorbs =freeorbs.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	donor =donor.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	acceptor =acceptor.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	
	acceptorPartners1 =acceptorPartners1.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	acceptorPartners2 =acceptorPartners2.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	
	donorPartners1 =donorPartners1.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)
	donorPartners2 =donorPartners2.permute(2,0,1).view(3,n_elements*n_FO*n_H).transpose(0,1)

	maxAng =maxAng.permute(2,0,1).view(maxAng.shape[2],n_elements*n_FO*n_H).transpose(0,1)
	minAng =minAng.permute(2,0,1).view(minAng.shape[2],n_elements*n_FO*n_H).transpose(0,1)
	#[22.0475,  2.3447, 22.1966],[21.2329,  3.5961, 23.6807]

	maksMissingH = ~hydrogens.eq(PADDING_INDEX)
	maskMissingFO= ~freeorbs.eq(PADDING_INDEX)

	maskMissing=(maksMissingH & maskMissingFO).sum(-1).type(maksMissingH.type())

	hydrogens = hydrogens[maskMissing]
	#test = test.view(-1)[maskMissing]
	freeorbs = freeorbs[maskMissing]
	donor = donor[maskMissing]

	acceptor = acceptor[maskMissing]
	maxAng = maxAng[maskMissing]
	minAng = minAng[maskMissing]
	
	donorPartners1 = donorPartners1[maskMissing]
	donorPartners2 = donorPartners2[maskMissing]
	acceptorPartners1 = acceptorPartners1[maskMissing]
	acceptorPartners2 = acceptorPartners2[maskMissing]

	angles, nanmask,planeNanmask,testz = geometry.get_interaction_angles(acceptor, donor,acceptorPartners1,donorPartners1,acceptorPartners2,donorPartners2, hydrogens, freeorbs)

	ia_mask = mask_angles(angles[:,:3], nanmask, minAng, maxAng) # no plane angle required

	distH =  torch.full((n_elements*n_H*n_FO,3), MISSING_ATOM,dtype=torch.float,device=dev)

	distFO = torch.full((n_elements*n_H*n_FO,3), MISSING_ATOM,dtype=torch.float,device=dev)
	finAngs = torch.full((n_elements*n_H*n_FO,4), MISSING_ATOM,dtype=torch.float,device=dev)
	interaction_angles_fullmask =torch.full(tuple([n_elements*n_H*n_FO]), 0,dtype=torch.bool,device=dev)
	
	distH[maskMissing]  = hydrogens
	distFO[maskMissing] = freeorbs
	finAngs[maskMissing] = angles

	test = torch.full((n_elements*n_H*n_FO,3), MISSING_ATOM,dtype=torch.float,device=dev)
	test[maskMissing] = testz
	interaction_angles_fullmask[maskMissing] = ia_mask


	distH = distH.transpose(1,0).view(3,n_elements,n_H*n_FO).permute(1,2,0)
	distFO = distFO.transpose(1,0).view(3,n_elements,n_H*n_FO).permute(1,2,0)
	finAngs = finAngs.transpose(1,0).view(4,n_elements,n_H*n_FO).permute(1,2,0)

	interaction_angles_fullmask = interaction_angles_fullmask.view(n_elements,n_H*n_FO)

	return distH,distFO,finAngs,interaction_angles_fullmask,planeNanmask
