#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  math_utils.py
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
import numpy as np
from madrax.sources.globalVariables import *
def dot2dVectors(v1,v2):

	a = torch.bmm(v1.unsqueeze(1), v2.unsqueeze(2))

	return a.squeeze(2)


def angle2dVectors(v1,v2):
	#magic_number = 57.2958 just for radians --> degrees. I use radians
	'''
	double dih=0.0;
	c.Normalize();
	d.Normalize();
	double scal=c*d;
	if (fabs(scal+1.0)<EPS) dih=180.0;
	else if (fabs(scal-1.0)<EPS) dih=0.0;
	else
		dih=57.2958*acos(scal);
	return dih;	
	'''
	
	v1Norm = torch.norm(v1,dim=1).unsqueeze(1)
	nanmask1=v1Norm.le(EPS)
	#v1Norm[nanmask1]=1.0
	#v1Norm = torch.masked_fill(v1Norm, nanmask1, 1)
	v3=v1/v1Norm.clamp(min= 0.0 + EPS)


	#v2=v2/(torch.norm(v2,dim=1).unsqueeze(1))
	v2Norm = torch.norm(v2,dim=1).unsqueeze(1).clamp(min= 0.0 + EPS)
	nanmask2=v2Norm.le(EPS)

	#v2Norm[nanmask2]=1.0
	#v2Norm = torch.masked_fill(v2Norm,nanmask2,1)

	v4=v2/(v2Norm)
	#print(v2Norm.min())
	#v4.register_hook(lambda x: print(x.min()))
	#assert np.isnan(torch.sum(v2).cpu().data.numpy())==False
	#assert np.isnan(torch.sum(v1).cpu().data.numpy())==False

	scal = dot2dVectors(v3,v4)
	#ssert np.isnan(torch.sum(scal).cpu().data.numpy())==False

	mask1 = torch.abs(scal+1.0).le(EPS)
	mask2 = torch.abs(scal-1.0).le(EPS)
	scal.masked_fill_(mask1,-1+EPS)
	scal.masked_fill_(mask2,1-EPS)
	g = scal#.clamp(-0.999,0.999)
	dih = torch.acos(g)
	
	return dih,~nanmask1 & ~nanmask2,"scimmia"
	
def dihedral2dVectors(i,j,k,l,testing=False):
	'''
	double dih;
	Vector3D jk=k-j;
	Vector3D c=(i-j)%jk; //cross product
	c.Normalize();
	Vector3D d=(l-k)%jk;
	d.Normalize();
	double scal=c*d;
	if (fabs(scal+1.0)<EPS) dih=180.0;
	else if (fabs(scal-1.0)<EPS) dih=0.0;
	else
		dih=57.2958*acos(scal);
	double chiral=(c%d)*jk;
	if (chiral<0) dih=-dih;
	return dih;
	'''
	
	jk = k-j
	
	c = torch.cross(i-j,jk, dim=1)
	cNorm = torch.norm(c,dim=1).unsqueeze(1)
	nanmask1=cNorm.le(EPS)
	cNorm = torch.masked_fill(cNorm, nanmask1, 1)
	#cNorm[nanmask1]=1.0
	c = c/cNorm

	d = torch.cross(l-k,jk, dim=1)
	dNorm = torch.norm(d,dim=1).unsqueeze(1)	
	nanmask2=dNorm.le(EPS)
	#dNorm[nanmask2]=1.0
	#dNorm = dNorm.clamp(min= 0.0 + EPS)#torch.masked_fill(dNorm, nanmask1, 1)
	dNorm = dNorm.clamp(min=0.0 + EPS)  # torch.masked_fill(dNorm, nanmask1, 1)

	d = d/dNorm
	d[torch.isnan(d)] = 0

	scal = dot2dVectors(c,d)

	#mask1 = torch.abs(scal+1.0).le(EPS)
	#mask2 = torch.abs(scal-1.0).le(EPS)

	mask1 = torch.abs(scal + 1.0).le(EPS)
	mask2 = torch.abs(scal - 1.0).le(EPS)
	scal=torch.masked_fill(scal,mask1, -1 + EPS)
	scal=torch.masked_fill(scal,mask2, 1 - EPS)

	#qui
	dih = torch.acos(scal)
	test = dih
	#dih.masked_fill_(mask1,math.pi)
	#dih.masked_fill_(mask2,0)
	
	chiral = dot2dVectors((torch.cross(c,d,dim=1)), jk) < 0
	dih[chiral]*=-1
	if testing:
		return dih, ~nanmask1 * ~nanmask2,test
	else:
		return dih,~nanmask1*~nanmask2

def plane_angle(p1,p2,p3,r1,r2,r3):
	'''
	Vector3D n1,n2;
	n1 = (r2 - r1)%(r3-r1);
	n2 = (p2 - p1)%(p3-p1);
	n1.Normalize();
	n2.Normalize();
	return Angle(n1,n2);
	'''	
	n1 = torch.cross(r2 - r1, r3 - r1,dim=1)
	n2 = torch.cross(p2 - p1, p3 - p1,dim=1)
	#assert np.isnan(torch.sum(n2).cpu().data.numpy())==False
	#assert np.isnan(torch.sum(n1).cpu().data.numpy())==False
	#n1 = n1/torch.norm(n1,dim=1).unsqueeze(1)
	#n2 = n2/torch.norm(n2,dim=1).unsqueeze(1)
	an,mask,_ = angle2dVectors(n1, n2)
	#an[~mask] = 1
	return an,mask

def point_to_plane_dist(r1,r2,r3,r0):
	'''
	n = (r2 - r1)%(r3-r1);
	n.Normalize();
	return (fabs((r0-r1)*n));	
	'''

	n = torch.cross(r2-r1, r3-r1, dim=-1)
	n = n/torch.norm(n,dim=-1).unsqueeze(-1)

	return torch.abs(dot2dVectors(r0-r1,n))	



def rotation_matrix_from_vectors(vec1, vec2):
	""" Find the rotation matrix that aligns vec1 to vec2
	:param vec1: A 3d "source" vector
	:param vec2: A 3d "destination" vector
	:return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
	"""
	#a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
	nElements = vec1.shape[0]
	a = vec1/torch.norm(vec1,dim=-1).unsqueeze(-1)
	b = vec2/torch.norm(vec2,dim=-1).unsqueeze(-1)

	v = torch.cross(a, b,dim=-1)
	c = dot2dVectors(a, b).squeeze(1)
	s = torch.norm(v,dim=-1)
	s = torch.where(s>EPS,s,torch.tensor(1.0).type_as(s))
	#kmat = torch.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

	kmat1 = torch.cat([torch.zeros(nElements,1).type_as(vec1), -v[:,[2]],v[:,[1]]],dim=1).unsqueeze(-1)
	kmat2 = torch.cat([v[:, [2]], torch.zeros(nElements, 1).type_as(vec1), -v[:,[0]]],dim=1).unsqueeze(-1)
	kmat3 = torch.cat([-v[:, [1]], v[:, [0]], torch.zeros(nElements, 1).type_as(vec1)], dim=1).unsqueeze(-1)
	kmat = torch.cat([kmat1,kmat2,kmat3],dim=-1).transpose(1,2)
	rotation_matrix = torch.eye(3,device=kmat.device).unsqueeze(0) + kmat + torch.bmm(kmat,kmat) * ((1 - c) / (s ** 2)).unsqueeze(-1).unsqueeze(-1)

	#rotation_matrix = torch.eye(3).unsqueeze(0) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

	return rotation_matrix
