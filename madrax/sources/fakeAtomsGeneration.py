#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright 2021 Gabriele Orlando
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
import torch
from madrax.sources.globalVariables import *
from madrax.sources import hashings,math_utils

def generateFakeAtomTensor(coords,partnersCoords,atomDescription,fake_atom_properties):

    atom_names = atomDescription[:,hashings.atom_description_hash["at_name"]]

    fake_atom_relative_positions = fake_atom_properties[atom_names.long()][... , :3]

    hCoords = add_hydrogen(coords,partnersCoords[:,0],partnersCoords[:,1],fake_atom_relative_positions)

    return hCoords

def add_hydrogen(r1, r2, r3, hpos):

    center = r1

    padding_mask = (r1[...,0] != PADDING_INDEX) & (r2[...,0] != PADDING_INDEX) & (r2[...,0] != PADDING_INDEX) & (hpos[...,0] != PADDING_INDEX).sum(-1).bool()

    local_z_to_lab = (r2[padding_mask] - r1[padding_mask]) / torch.norm(r2[padding_mask] - r1[padding_mask],dim=-1).unsqueeze(-1).expand(-1,3)
    zx_vec = r3[padding_mask] - r1[padding_mask]  # )/np.linalg.norm(r3-r1)

    local_x_to_lab = zx_vec - math_utils.dot2dVectors(zx_vec, local_z_to_lab) * local_z_to_lab  # design a vector orthogonal to local_z

    locx_norm = torch.norm(local_x_to_lab,dim=-1).unsqueeze(-1)

    normExceptionValue = torch.tensor([1., 0., 0.]).type_as(r1)
    local_x_to_lab = torch.where(locx_norm > EPS, local_x_to_lab / locx_norm, normExceptionValue)


    x = torch.tensor( [[1.0, 0.0, 0.0]]).type_as(r1).expand(local_x_to_lab.shape)
    z = torch.tensor( [[0.0, 0.0, 1.0]]).type_as(r1).expand(local_x_to_lab.shape)

    mx = math_utils.rotation_matrix_from_vectors(x, local_x_to_lab)
    z = mx.bmm(z.unsqueeze(2).expand(-1,3,1)).squeeze(-1)

    mz = math_utils.rotation_matrix_from_vectors(z, local_z_to_lab)

    fakeAtomsPadding = (hpos[padding_mask][..., 0]!=PADDING_INDEX)


    expandedMX = mx.unsqueeze(1).expand(-1,hpos.shape[1],-1,-1)[fakeAtomsPadding]
    expandedMZ = mz.unsqueeze(1).expand(-1, hpos.shape[1], -1, -1)[fakeAtomsPadding]
    expandedCenter = center.unsqueeze(1).expand(-1,hpos.shape[1],-1)[padding_mask][fakeAtomsPadding]

    hcoords = expandedMX.bmm(hpos[padding_mask][fakeAtomsPadding].unsqueeze(2).expand(-1,3,1)).squeeze(-1)
    hcoords_fin = expandedMZ.bmm(hcoords.unsqueeze(2).expand(-1,3,1)).squeeze(-1) + expandedCenter

    full_mask = padding_mask.unsqueeze(1).expand(-1,hpos.shape[1]).clone()
    full_mask[padding_mask] = fakeAtomsPadding

    final_coords = torch.full((r1.shape[0],hpos.shape[1],3),PADDING_INDEX,device=r1.device)
    final_coords[full_mask] = hcoords_fin


    return final_coords