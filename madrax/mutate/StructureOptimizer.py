#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  learn_forcefield_main.py
#  
#  Copyright 2019 Gabriel Orlando <orlando.gabriele89@gmail.com>
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

letters={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ASN': 'N', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'HIS': 'H', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'} #gli aminoacidi, che male qui non fanno

import torch,math,time
from madrax.mutate import rotator
from madrax.sources import hashings
import madrax.utils as utils

def optimize(model, coords, info_tensors, epochs=5,verbose=False,learning_rate=0.1,backbone_rotation=False,atnames=None):
	dev=model.device
	rotator_obj = rotator.RotateStruct()

	start = time.time()
	atom_number, atom_description, coordsIndexingAtom, partnersIndexingAtom, angle_indices,alternativeMask = info_tensors

	batch = atom_description[:,hashings.atom_description_hash["batch"]].max()+1
	maxseq = atom_description[:,hashings.atom_description_hash["resnum"]].max()+1
	maxchain = atom_description[:,hashings.atom_description_hash["chain"]].max()+1

	naltern = alternativeMask.shape[-1]
	dummy_rotationSC = torch.zeros((batch, maxchain, maxseq,naltern ,8), device=dev,dtype=torch.float)

	dummy_translationBB = torch.zeros((batch, maxchain, maxseq,naltern, 3), device=dev,dtype=torch.float)
	dummy_rotationSC.requires_grad=True
	dummy_translationBB.requires_grad = True

	start=time.time()

	optimizer = torch.optim.Adam([
								{'params': dummy_translationBB, "lr":learning_rate},
								{'params': dummy_rotationSC, 'lr': learning_rate}
								],amsgrad=True,eps=0.1)

	scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.9, patience=10, verbose=False, threshold=0.0001, threshold_mode='rel', cooldown=0, min_lr=0, eps=1e-08)
	#early_stopping = EarlyStopping(patience=100)

	optimizer.zero_grad()

	startRotation = []
	giving_up = torch.zeros((batch,maxchain,maxseq),dtype=torch.long,device=dev)

	for e in range(epochs):
		start_time = time.time()

		dummy_rotationSCAngles = torch.sin(dummy_rotationSC) * math.pi

		start = time.time()
		#if e>20:
		#	backbone_rotation=True
		coords_local = rotator_obj(coords, info_tensors, dummy_rotationSCAngles,backbone_rotation=backbone_rotation)

		startRotation+=[start-time.time()]

		yp = model(coords_local,info_tensors)
		loss = yp.sum()
		loss.backward()
		optimizer.step()
		optimizer.zero_grad()
		scheduler.step(float(loss.data.cpu()))
		if verbose:
			print(" \t optimizing epoch:", e, "loss:", round(yp.sum().data.cpu().tolist(), 4), "time ",round((time.time() - start_time), 4))
		#utils.writepdbNew(coords_local.cpu().data, atnames, pdb_names=["testOptim" + str(e)])

	return yp,coords_local
