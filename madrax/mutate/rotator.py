#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  learn_forcefield_main.py
#
#  Copyright 2020 Gabriel Orlando <orlando.gabriele89@gmail.com>
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
from madrax.sources import hashings
import torch

from madrax.sources.globalVariables import *

class RotateStruct():
    def __init__(self):
        self.hashing = hashings.atom_hash

        self.torsion_hash ={

                                "chi1":0,
                                "chi2": 1,
                                "chi3": 2,
                                "chi4": 3,
                                "chi5": 4,
                             }



        self.torsion_axesHR ={
                            "chi1": {

                                #"PRO": ["CA", "CB"],
                                "GLN": ["CA", "CB"],
                                "VAL": ["CA", "CB"],
                                "ASN": ["CA", "CB"],
                                "THR": ["CA", "CB"],
                                #"ALA": ["CA", "CB"],
                                "ASP": ["CA", "CB"],
                                "PHE": ["CA", "CB"],
                                "LEU": ["CA", "CB"],
                                "SER": ["CA", "CB"],
                                "CYS": ["CA", "CB"],
                                "ILE": ["CA", "CB"],
                                "TRP": ["CA", "CB"],
                                "ARG": ["CA", "CB"],
                                "LYS": ["CA", "CB"],
                                "TYR": ["CA", "CB"],
                                "GLU": ["CA", "CB"],
                                "MET": ["CA", "CB"],
                                "HIS": ["CA", "CB"],
                                    },
                            "chi2": {
                                    "GLN": ["CB", "CG"],
                                    "ASN": ["CB", "CG"],
                                    "ASP": ["CB", "CG"],
                                    "PHE": ["CB", "CG"],
                                    "LEU": ["CB", "CG"],
                                    "TRP": ["CB", "CG"],
                                    "ARG": ["CB", "CG"],
                                    "LYS": ["CB", "CG"],
                                    "TYR": ["CB", "CG"],
                                    "GLU": ["CB", "CG"],
                                    "MET": ["CB", "CG"],
                                    "HIS": ["CB", "CG"]
                                    },

                            "chi3": {
                                "GLN": ["CG", "CD"],
                                "ARG": ["CG", "CD"],
                                "LYS": ["CG", "CD"],
                                "GLU": ["CG", "CD"],
                                "MET": ["CG", "SD"],
                                   },
                            "chi4": {

                                "ARG": ["CD", "NE"],
                                "LYS": ["CD", "CE"],
                            },
                            "chi5": {

                                "ARG": ["NE","CZ"],
                            },
                            "chi5": {

                                "ARG": ["NE", "CZ"],
                            },

        }

        torsion_axesAffectedAnglesHR = {
            "chi1": {
                "ILE": ["CB", "CG1", "CG2", "CD1"],
                # "PRO": ["CA", "CB"],
                "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
                "VAL": ["CB", "CG1", "CG2"],
                "ASN": ["CB", "ND2", "CG", "OD1"],
                "THR": ["CB", "OG1", "CG2"],
                "ASP": ["CB", "OD1", "OD2", "CG"],
                "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CZ", "CE2", "CE2", "RC", "RE"],
                "LEU": ["CB", "CG", "CD1", "CD2"],
                "SER": ["CB", "OG"],
                "CYS": ["CB", "SG"],
                "TRP": ["CB", "CG", "CD1", "NE1", "CE2", "CD1", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
                "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "ARG"],
                "LYS": ["CB", "CG", "CD", "CE", "NZ"],
                "TYR": ["CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2", "OH", "RC", "RE"],
                "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
                "MET": ["CB", "CG", "SD", "CE"],
                "HIS": ["CB", "CG", "CD2", "NE2", "CE1", "ND1", "NE2H1S", "NE2H2S", "ND1H1S", "ND1H2S"],

            },

            "chi2": {
                "GLN": ["CG", "CD", "OE1", "NE2"],
                "ASN": ["CG", "ND2", "OD1"],
                "ASP": ["CG", "OD1", "OD2"],
                "PHE": ["CG", "CD1", "CE1", "CD2", "CZ", "CE2", "CE2", "RC", "RE"],
                "LEU": ["CG", "CD1", "CD2"],
                "TRP": ["CG", "CD1", "NE1", "CE2", "CD1", "CD2", "CE3", "CZ2", "CZ3", "CH2"],
                "ARG": ["CG", "CD", "NE", "CZ", "NH1", "NH2", "ARG"],
                "LYS": ["CG", "CD", "CE", "NZ"],
                "TYR": ["CG", "CD1", "CE1", "CZ", "CE2", "CD2", "OH", "RC", "RE"],
                "GLU": ["CG", "CD", "OE1", "OE2"],
                "MET": ["CG", "SD", "CE"],
                "HIS": ["CG", "CD2", "NE2", "CE1", "ND1", "NE2H1S", "NE2H2S", "ND1H1S", "ND1H2S"],
            },

            "chi3": {
                "GLN": ["CD", "OE1", "NE2"],
                "ARG": ["CD", "NE", "CZ", "NH1", "NH2", "ARG"],
                "LYS": ["CD", "CE", "NZ"],
                "GLU": ["CD", "OE1", "OE2"],
                "MET": ["SD", "CE"],
            },

            "chi4": {

                "ARG": ["NE", "CZ", "NH1", "NH2", "ARG"],
                "LYS": ["CE", "NZ"],
            },

            "chi5": {

                "ARG": ["CZ", "NH1", "NH2", "ARG"],
            }
        }
        self.torsion_axesAffectedAnglesHR = {}
        for tors in torsion_axesAffectedAnglesHR.keys():
            self.torsion_axesAffectedAnglesHR[self.torsion_hash[tors]] = {}
            for res in torsion_axesAffectedAnglesHR[tors]:
                self.torsion_axesAffectedAnglesHR[self.torsion_hash[tors]][hashings.resi_hash[res]]=[]
                for at in torsion_axesAffectedAnglesHR[tors][res]:
                    self.torsion_axesAffectedAnglesHR[self.torsion_hash[tors]][hashings.resi_hash[res]] += [hashings.atom_hash[res][at]]
        #self.torsion_axesAffectedAnglesHR = list(set(self.torsion_axesAffectedAnglesHR))

        #self.residue_coords_hashing =

        self.torsion_axes = {}
        self.torsion_axesAffectedAngles = {}

        for i in self.torsion_axesHR.keys():
            self.torsion_axes[self.torsion_hash[i]] = {}

            for res in  self.torsion_axesHR[i]:

                self.torsion_axes[self.torsion_hash[i]][hashings.resi_hash[res]] = [self.hashing[res][self.torsion_axesHR[i][res][0]],self.hashing[res][self.torsion_axesHR[i][res][1]]]

        self.Natoms = []
        self.Catoms = []
        self.Oatoms = []
        self.CAatoms = []
        for res1 in hashings.resi_hash.keys():
            self.Natoms += [hashings.atom_hash[res1]["N"]]
            self.Natoms += [hashings.atom_hash[res1]["tN"]]
            self.Catoms += [hashings.atom_hash[res1]["C"]]
            self.Oatoms += [hashings.atom_hash[res1]["O"]]
            self.CAatoms += [hashings.atom_hash[res1]["CA"]]

        self.Natoms = list(set(self.Natoms))
        self.Catoms = list(set(self.Catoms))
        self.Oatoms = list(set(self.Oatoms))
        self.CAatoms = list(set(self.CAatoms))


    def generalRotationMatrix(self, angle): ## used for backbone cloud

        assert angle.shape[-1] == 3
        batch = angle.shape[0]
        cos = torch.cos(angle)
        sin = torch.sin(angle)

        a00=(cos[:,:,:,:,0] * cos[:,:,:,:,1]).unsqueeze(-1).unsqueeze(-1)
        a10=(sin[:,:,:,:,0] * cos[:,:,:,:,1]).unsqueeze(-1).unsqueeze(-1)
        a20 =(-sin[:,:,:,:,1]).unsqueeze(-1).unsqueeze(-1)

        a01 = (cos[:,:,:,:,0] * sin[:,:,:,:,1] * sin[:,:,:,:,2] - sin[:,:,:,:,0] * cos[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)
        a11 = (sin[:,:,:,:,0] * sin[:,:,:,:,1] * sin[:,:,:,:,2] + cos[:,:,:,:,0] * cos[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)
        a21 = (cos[:,:,:,:,1] * sin[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)

        a02 = (cos[:,:,:,:,0] * sin[:,:,:,:,1] * cos[:,:,:,:,2] + sin[:,:,:,:,0] * sin[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)
        a12 = (sin[:,:,:,:,0] * sin[:,:,:,:,1] * cos[:,:,:,:,2] - cos[:,:,:,:,0] * sin[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)
        a22 = (cos[:,:,:,:,1] * cos[:,:,:,:,2]).unsqueeze(-1).unsqueeze(-1)

        a0 = torch.cat([a00, a10, a20], dim=-2)
        a1 = torch.cat([a01, a11, a21], dim=-2)
        a2 = torch.cat([a02, a12, a22], dim=-2)
        mat = torch.cat([a0, a1, a2], dim=-1)
        return mat

    def rotationMatrix(self, angle, axis):
        batch = angle.shape[0]
        cos = torch.cos(angle)
        sin = torch.sin(angle)

        a00=(cos + axis[:, 0] ** 2 * (1 - cos)).unsqueeze(-1).unsqueeze(-1)
        a10=(axis[:, 0] * axis[:, 1] * (1 - cos) + axis[:, 2] * sin).unsqueeze(-1).unsqueeze(-1)
        a20 =(axis[:, 0] * axis[:, 2] * (1 - cos) - axis[:, 1] * sin).unsqueeze(-1).unsqueeze(-1)

        a01 = (axis[:, 0] * axis[:, 1] * (1 - cos) - axis[:, 2] * sin).unsqueeze(-1).unsqueeze(-1)
        a11 = (cos + axis[:, 1] ** 2 * (1 - cos)).unsqueeze(-1).unsqueeze(-1)
        a21 = (axis[:, 1] * axis[:, 2] * (1 - cos) + axis[:, 0] * sin).unsqueeze(-1).unsqueeze(-1)

        a02 = (axis[:, 0] * axis[:, 2] * (1 - cos) + axis[:, 1] * sin).unsqueeze(-1).unsqueeze(-1)
        a12  = (axis[:, 1] * axis[:, 2] * (1 - cos) - axis[:, 0] * sin).unsqueeze(-1).unsqueeze(-1)
        a22  = (cos + axis[:, 2] ** 2 * (1 - cos)).unsqueeze(-1).unsqueeze(-1)

        a0 = torch.cat([a00, a10, a20], dim=1)
        a1 = torch.cat([a01, a11, a21], dim=1)
        a2 = torch.cat([a02, a12, a22], dim=1)
        mat = torch.cat([a0, a1, a2], dim=-1)
        '''
        mat[:, 0, 0] = cos + axis[:, 0] ** 2 * (1 - cos)
        mat[:, 1, 0] = axis[:, 0] * axis[:, 1] * (1 - cos) + axis[:, 2] * sin
        mat[:, 2, 0] = axis[:, 0] * axis[:, 2] * (1 - cos) - axis[:, 1] * sin

        mat[:, 0, 1] = axis[:, 0] * axis[:, 1] * (1 - cos) - axis[:, 2] * sin
        mat[:, 1, 1] = cos + axis[:, 1] ** 2 * (1 - cos)
        mat[:, 2, 1] = axis[:, 1] * axis[:, 2] * (1 - cos) + axis[:, 0] * sin

        mat[:, 0, 2] = axis[:, 0] * axis[:, 2] * (1 - cos) + axis[:, 1] * sin
        mat[:, 1, 2] = axis[:, 1] * axis[:, 2] * (1 - cos) - axis[:, 0] * sin
        mat[:, 2, 2] = cos + axis[:, 2] ** 2 * (1 - cos)
        '''
        return mat

    def preOptimizationStep(self): #building structure-specific data strctures that are not gonna change during optimization


        self.torsion_axesAffectedAngles[self.torsion_hash[i]][hashings.resi_hash[res]] = atnme.eq(self.hashing[res][atname])

    def moveCloudBackBone(self, coords, fake_atoms, atom_description, rotation, translation, partial_optimization,alternatives):
        resnum = atom_description[:, :, hashings.property_hashings["atom_description"]["res_num"]]
        chain = atom_description[:, :, hashings.property_hashings["atom_description"]["chain"]]
        atname = atom_description[:, :, hashings.property_hashings["atom_description"]["atom_types"]]
        maxchain = chain.max() + 1
        maxres = resnum.max() + 1
        nbatch = coords.shape[0]
        natoms = resnum.shape[1]
        nrot =coords.shape[2]
        nalter = rotation.shape[3]
        nfakeatoms = fake_atoms.shape[3]

        existing_fake_atoms_mask = ~fake_atoms.eq(MISSING_ATOM)[:, :, :, 0]

        if partial_optimization is not None:
            padding_mask = ~resnum.eq(PADDING_INDEX) & partial_optimization
        else:
            padding_mask = ~resnum.eq(PADDING_INDEX)

        N = ((atname == hashings.atom_hash["ALA"]["N"]) | (atname == hashings.atom_hash["PRO"]["N"])).unsqueeze(-1)
        CA = ((atname == hashings.atom_hash["ALA"]["CA"]) | (atname == hashings.atom_hash["PRO"]["CA"])).unsqueeze(-1)
        C = ((atname == hashings.atom_hash["ALA"]["C"]) | (atname == hashings.atom_hash["PRO"]["C"])).unsqueeze(-1)

        rotated_coords = coords.clone()
        rotated_fake_atoms = fake_atoms.clone()

        batch_ind = torch.arange(nbatch,device = coords.device).unsqueeze(-1).expand(nbatch, natoms)
        rotInd = torch.arange(nrot,device = coords.device).unsqueeze(0).unsqueeze(0).expand(nbatch, natoms,nrot)

        batch_ind_resi = torch.arange(nbatch,device = coords.device).unsqueeze(-1).unsqueeze(-1).expand(nbatch, maxchain, maxres)
        chain_resi = torch.arange(maxchain,device = coords.device).unsqueeze(0).unsqueeze(-1).expand(nbatch, maxchain, maxres)

        rotInd = torch.arange(nrot,device = coords.device).unsqueeze(0).unsqueeze(0).expand(nbatch, natoms,nrot)

        for alt in range(nalter):
            centers = torch.full((nbatch,maxchain,maxres,nalter,nrot,3),PADDING_INDEX,device = coords.device)

            altMask =  alternatives[:,:,alt]
            altIndex = torch.full((nbatch,natoms),alt,device=self.dev)
            #gestisci le alter con roba da ciclo for


            rotation_matrix = self.generalRotationMatrix(rotation[:,:,:,alt])

            existing_rotations = ~coords[:, :, :, 0].eq(PADDING_INDEX)

            if partial_optimization is not None:
                padding = existing_rotations & partial_optimization.unsqueeze(-1) & altMask.unsqueeze(-1)
            else:
                padding = existing_rotations & altMask.unsqueeze(-1)

            centers = torch.full((nbatch,maxchain,maxres,nalter,nrot,3),PADDING_INDEX,device = coords.device)
            centers = centers.index_put((batch_ind.unsqueeze(-1).expand(-1,-1,nrot)[CA & padding], chain.unsqueeze(-1).expand(-1,-1,nrot)[CA & padding], resnum.unsqueeze(-1).expand(-1,-1,nrot)[CA & padding], altIndex.unsqueeze(-1).expand(-1,-1,nrot)[CA & padding],rotInd[CA & padding]), coords[CA & padding])

            centersAtom = centers[batch_ind.unsqueeze(-1).expand(-1,-1,nrot)[padding],chain.unsqueeze(-1).expand(-1,-1,nrot)[padding],resnum.unsqueeze(-1).expand(-1,-1,nrot)[padding],altIndex.unsqueeze(-1).expand(-1,-1,nrot)[padding],rotInd[padding]]

            assert not True in centersAtom.eq(PADDING_INDEX)
            rotationsAtoms = rotation_matrix[batch_ind.unsqueeze(-1).expand(-1,-1,nrot)[padding],chain.unsqueeze(-1).expand(-1,-1,nrot)[padding],resnum.unsqueeze(-1).expand(-1,-1,nrot)[padding],rotInd[padding]]
            translationAtoms = translation[batch_ind.unsqueeze(-1).expand(-1,-1,nrot)[padding],chain.unsqueeze(-1).expand(-1,-1,nrot)[padding],resnum.unsqueeze(-1).expand(-1,-1,nrot)[padding],alt,rotInd[padding]]





            rot_coords = torch.bmm(rotationsAtoms,(coords[padding] - centersAtom).unsqueeze(-1)).squeeze(-1) + centersAtom + translationAtoms

            # fakeatoms
            fakeatomsExisting_Mask = ((~fake_atoms[:,:,:,:,0].eq(PADDING_INDEX)) & padding.unsqueeze(-1).expand(fake_atoms.shape[:-1])) & altMask.unsqueeze(-1).unsqueeze(-1)
            fake_atnumber = torch.arange(0,nfakeatoms,device=coords.device).unsqueeze(0).unsqueeze(0).unsqueeze(0).expand(fakeatomsExisting_Mask.shape)[fakeatomsExisting_Mask]
            batch_indFake = batch_ind.unsqueeze(-1).unsqueeze(-1).expand(nbatch,natoms,nrot,nfakeatoms)[fakeatomsExisting_Mask]
            chainFake = chain.unsqueeze(-1).unsqueeze(-1).expand(nbatch,natoms,nrot,nfakeatoms)[fakeatomsExisting_Mask]
            resnumFake = resnum.unsqueeze(-1).unsqueeze(-1).expand(nbatch, natoms,nrot, nfakeatoms)[fakeatomsExisting_Mask]
            nrotFake = rotInd.unsqueeze(-1).expand(nbatch, natoms,nrot, nfakeatoms)[fakeatomsExisting_Mask]

            centersFake = centers[batch_indFake, chainFake, resnumFake, alt,nrotFake ]

            rotation_matrixFake = rotation_matrix.unsqueeze(-3).expand(nbatch,maxchain,maxres,nrot,nfakeatoms,3,3)[batch_indFake,chainFake,resnumFake,nrotFake,fake_atnumber]#[fakeatomsExisting_Mask]
            translation_Fake = translation.unsqueeze(-2).expand(nbatch,maxchain,maxres,nalter,nrot,nfakeatoms,3)[batch_indFake,chainFake,resnumFake,alt,nrotFake,fake_atnumber]#[fakeatomsExisting_Mask]

            assert not True in centersFake.eq(PADDING_INDEX)

            rot_fakeAtoms = torch.bmm(rotation_matrixFake,(fake_atoms[fakeatomsExisting_Mask] - centersFake).unsqueeze(-1)).squeeze(-1) + centersFake + translation_Fake

            fake_atoms[fakeatomsExisting_Mask] = rot_fakeAtoms
            coords[padding] = rot_coords

        return coords,fake_atoms
        #rot_fakeAtoms =
    def rotateBackBone(self, coords,atom_description,rotation,actually_moving):


        if atom_description[:, hashings.atom_description_hash["alternative"]].max() > 0:
            raise("rotation backbone is only valid for a protein with no mutations. Run mutate with keep_WT==False")

        batch_ind = atom_description[:, hashings.atom_description_hash["batch"]].long()
        resnum = atom_description[:, hashings.atom_description_hash["resnum"]].long()
        chain = atom_description[:, hashings.atom_description_hash["chain"]].long()
        atName = atom_description[:, hashings.atom_description_hash["at_name"]].long()

        batch = batch_ind.max() + 1
        chainMax = chain.max() + 1

        for i, natom in enumerate(self.Natoms):
            if i == 0:
                N = atName.eq(natom)
            else:
                N += atName.eq(natom)

        for i, catom in enumerate(self.Catoms):
            if i == 0:
                C = atName.eq(catom)
            else:
                C += atName.eq(catom)

        for i, catom in enumerate(self.CAatoms):
            if i == 0:
                CA = atName.eq(catom)
            else:
                CA += atName.eq(catom)


        rotated_coords = coords.clone()


        rotations_to_do = torch.cat(
                                    [
                                        batch_ind[actually_moving].unsqueeze(-1),
                                        chain[actually_moving].unsqueeze(-1),
                                        resnum[actually_moving].unsqueeze(-1),
                                    ],dim=-1).unique(dim=0)

        for rot in rotations_to_do:


            targetBatch,targetChain,targetResnum = rot
            targeting_mask = chain.eq(targetChain) & batch_ind.eq(targetBatch)

            ### rotating omega ###
            axis2 = ((resnum == targetResnum) & N) & targeting_mask
            axis1 = ((resnum == (targetResnum - 1)) & C) & targeting_mask

            if True in axis2 and True in axis1:
                N_curr = rotated_coords[axis2]
                C_prev = rotated_coords[axis1]
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum,0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                            (rotated_coords[targeting_mask][affected] - masked_translation.unsqueeze(0)).squeeze(1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(affected.unsqueeze(-1).repeat(1, 3), fat)

            ### rotating phi ###
            axis2 = ((resnum == targetResnum) & CA) & targeting_mask
            axis1 = ((resnum == targetResnum) & N) & targeting_mask

            if True in axis2 and True in axis1:
                N_curr = rotated_coords[axis2]
                C_prev = rotated_coords[axis1]
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum,0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                            (rotated_coords[targeting_mask][affected] - masked_translation.unsqueeze(0)).squeeze(1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(affected.unsqueeze(-1).repeat(1, 3), fat)

            ### rotating psi ###
            axis2 = ((resnum == targetResnum) & C) & targeting_mask
            axis1 = ((resnum == targetResnum) & CA) & targeting_mask

            if True in axis2 and True in axis1:
                N_curr = rotated_coords[axis2]
                C_prev = rotated_coords[axis1]
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum, 0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(
                    int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                                      (rotated_coords[targeting_mask][
                                                           affected] - masked_translation.unsqueeze(0)).squeeze(
                                                          1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(
                    affected.unsqueeze(-1).repeat(1, 3), fat)


        return rotated_coords





    def buildBackBone(self, coords,atom_description,rotation,actually_moving):


        if atom_description[:, hashings.atom_description_hash["alternative"]].max() > 0:
            raise("rotation backbone is only valid for a protein with no mutations. Run mutate with keep_WT==False")

        batch_ind = atom_description[:, hashings.atom_description_hash["batch"]].long()
        resnum = atom_description[:, hashings.atom_description_hash["resnum"]].long()
        chain = atom_description[:, hashings.atom_description_hash["chain"]].long()
        atName = atom_description[:, hashings.atom_description_hash["at_name"]].long()

        batch = batch_ind.max() + 1
        chainMax = chain.max() + 1

        for i, natom in enumerate(self.Natoms):
            if i == 0:
                N = atName.eq(natom)
            else:
                N += atName.eq(natom)

        for i, catom in enumerate(self.Catoms):
            if i == 0:
                C = atName.eq(catom)
            else:
                C += atName.eq(catom)

        for i, catom in enumerate(self.CAatoms):
            if i == 0:
                CA = atName.eq(catom)
            else:
                CA += atName.eq(catom)


        rotated_coords = coords.clone()


        rotations_to_do = torch.cat(
                                    [
                                        batch_ind[actually_moving].unsqueeze(-1),
                                        chain[actually_moving].unsqueeze(-1),
                                        resnum[actually_moving].unsqueeze(-1),
                                    ],dim=-1).unique(dim=0)

        for rot in rotations_to_do:


            targetBatch,targetChain,targetResnum = rot
            targeting_mask = chain.eq(targetChain) & batch_ind.eq(targetBatch)

            ### rotating omega ###
            axis2 = ((resnum == targetResnum) & N) & targeting_mask
            axis1 = ((resnum == (targetResnum - 1)) & C) & targeting_mask

            if True in axis2 and True in axis1:

                #phi
                C_prev = rotated_coords[((resnum == (targetResnum - 1)) & C) & targeting_mask]
                N_cur = rotated_coords[((resnum == (targetResnum )) & N) & targeting_mask]
                CA_cur = rotated_coords[((resnum == (targetResnum)) & CA) & targeting_mask]
                C_cur = rotated_coords[((resnum == (targetResnum)) & C) & targeting_mask]

                replacingAtom = C_cur

                from sources import math_utils
                bond_angles = math_utils.angle2dVectors(CA_cur-N_cur,C_cur-N_cur)
                distance
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum,0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                            (rotated_coords[targeting_mask][affected] - masked_translation.unsqueeze(0)).squeeze(1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(affected.unsqueeze(-1).repeat(1, 3), fat)

            ### rotating phi ###
            axis2 = ((resnum == targetResnum) & CA) & targeting_mask
            axis1 = ((resnum == targetResnum) & N) & targeting_mask

            if True in axis2 and True in axis1:
                N_curr = rotated_coords[axis2]
                C_prev = rotated_coords[axis1]
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum,0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                            (rotated_coords[targeting_mask][affected] - masked_translation.unsqueeze(0)).squeeze(1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(affected.unsqueeze(-1).repeat(1, 3), fat)

            ### rotating psi ###
            axis2 = ((resnum == targetResnum) & C) & targeting_mask
            axis1 = ((resnum == targetResnum) & CA) & targeting_mask

            if True in axis2 and True in axis1:
                N_curr = rotated_coords[axis2]
                C_prev = rotated_coords[axis1]
                atomAxis = torch.cat([C_prev, N_curr], dim=0)
                affected = resnum[targeting_mask].ge(targetResnum)

                norm = torch.norm(atomAxis[1, :] - atomAxis[0, :], dim=-1).unsqueeze(-1)
                masked_unitary_axis = (atomAxis[1, :] - atomAxis[0, :]) / norm
                masked_translation = atomAxis[0, :]

                masked_rotation = rotation[targetBatch, targetChain, targetResnum, 0, 2].unsqueeze(0)
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis.unsqueeze(0)).expand(
                    int(affected.sum()), 3, 3)
                fat = masked_translation + (torch.bmm(rotationMatrix,
                                                      (rotated_coords[targeting_mask][
                                                           affected] - masked_translation.unsqueeze(0)).squeeze(
                                                          1).unsqueeze(-1)).squeeze(-1))

                rotated_coords[targeting_mask] = rotated_coords[targeting_mask].masked_scatter(
                    affected.unsqueeze(-1).repeat(1, 3), fat)


        return rotated_coords




    def __call__(self,coordinates,info_tensors,rotation,rotation_mask=None,backbone_rotation=True):



        dev = coordinates.device
        atom_number, atom_description, coordsIndexingAtom, partnersIndexingAtom, angle_indices, alternativeMask = info_tensors

        coords = coordinates[atom_description[:, 0].long(), coordsIndexingAtom]

        rotated_coords = coords.clone()


        batch = rotation.shape[0]
        nchain = rotation.shape[1]
        nresi = rotation.shape[2] # resnum.max() + 1
        nalter = rotation.shape[3]
        ntorsions = 8

        if rotation_mask is  None:
            padding_mask = coords[..., 0] != (PADDING_INDEX)
        else:

            padding_mask = coords[..., 0] != (PADDING_INDEX)

            rotAtomsMask = rotation_mask[atom_description[:, hashings.atom_description_hash["batch"]].long(),
                                         atom_description[:, hashings.atom_description_hash["chain"]].long(),
                                         atom_description[:, hashings.atom_description_hash["resnum"]].long(),
                                         atom_description[:, hashings.atom_description_hash["alternative"]].long()
                                        ]
            padding_mask = padding_mask & rotAtomsMask

        L = coords.shape[0]

        atnumber = torch.arange(L,dtype=torch.long,device=dev)#.unsqueeze(0).expand(batch,L)



        axis_ind = torch.full((batch, nchain, nresi,nalter,ntorsions, 2), PADDING_INDEX,dtype=torch.long,device=dev)

        affectedAt = torch.zeros(L,nalter, ntorsions, dtype=torch.bool, device=dev)

        if backbone_rotation:

            rotated_coords = self.rotateBackBone(rotated_coords, atom_description, rotation[..., :3], padding_mask)

        resname = atom_description[:, hashings.atom_description_hash["resname"]].long()
        aname = atom_description[:, hashings.atom_description_hash["at_name"]].long()
        batch_ind = atom_description[:, hashings.atom_description_hash["batch"]].long()
        resnum = atom_description[:, hashings.atom_description_hash["resnum"]].long()
        chain = atom_description[:, hashings.atom_description_hash["chain"]].long()

        for torsion in self.torsion_axes.keys():
            #testmask = resname.eq(12)# ile only

            for res in self.torsion_axes[torsion].keys():
                at1 = self.torsion_axes[torsion][res][0]
                at2 = self.torsion_axes[torsion][res][1]
                for alt in range(nalter):

                    mask1 = resname.eq(res) & aname.eq(at1) & alternativeMask[:,alt]
                    mask2 = resname.eq(res) & aname.eq(at2) & alternativeMask[:,alt]

                    axis_ind[:,:,:,alt, torsion,0] = axis_ind[:, :, :,alt, torsion, 0].index_put([batch_ind[mask1], chain[mask1], resnum[mask1]], atnumber[mask1], accumulate = False)
                    axis_ind[:,:,:,alt, torsion,1] = axis_ind[:, :, :,alt, torsion, 1].index_put([batch_ind[mask2], chain[mask2], resnum[mask2]], atnumber[mask2], accumulate = False)

                    for i, at in enumerate(self.torsion_axesAffectedAnglesHR[torsion][res]):

                        affectedAt[: ,alt, torsion] += (aname.eq(at) & alternativeMask[:,alt])

                    pass
                    #self.torsion_axesAffectedAngles[self.torsion_hash[i]][hashings.resi_hash[res]] = atnme.eq(self.hashing[res][atname])
        #calculating axis
        atomAxis_ind = axis_ind[batch_ind[padding_mask],chain[padding_mask], resnum[padding_mask]]#.permute(1, 2, 0).reshape(nrotating_atoms, 2, batch, L).permute(2, 3, 0, 1)

        rotationAtom = rotation[batch_ind[padding_mask], chain[padding_mask], resnum[padding_mask],:,3:]
        existing_mask = ~atomAxis_ind.eq(PADDING_INDEX)[:,:,:,0] # it tells me if the torsion of a residue exists, not a single coord or atom
        #masked_axis = atomAxis_ind[existing_mask]
        #masked_rotation = rotationAtom[existing_mask]

        affectedAt = affectedAt[padding_mask]
        batch_ind = batch_ind[padding_mask]

        for torsion in self.torsion_axes.keys():
            for alt in range(nalter):
                todoMask = existing_mask[:,alt, torsion] & affectedAt[:,alt, torsion]
                masked_axes_ind = atomAxis_ind[:,alt,torsion,:][todoMask]

                atomAxis = torch.cat([rotated_coords[masked_axes_ind[:,0]].unsqueeze(1),
                                      rotated_coords[masked_axes_ind[:,1]].unsqueeze(1)],dim=1)

                existingRotation = ~(atomAxis[:,0,0].eq(PADDING_INDEX)  | atomAxis[:,1,0].eq(PADDING_INDEX))

                norm = torch.norm(atomAxis[:, 1, :][existingRotation] - atomAxis[:, 0, :][existingRotation], dim=-1).unsqueeze(-1)

                masked_unitary_axis = (atomAxis[:, 1, :][existingRotation]  - atomAxis[:, 0, :][existingRotation] ) / norm

                masked_translation = atomAxis[:, 0, :][existingRotation]

                masked_rotation = rotationAtom[:,alt,torsion][todoMask][existingRotation]


                #masked_unitary_axis = unitary_axis[affectedAt[:,:,torsion]]
                #masked_translation = translation[affectedAt[:,:,torsion]]
                rotationMatrix = self.rotationMatrix(masked_rotation, masked_unitary_axis)
                fat = masked_translation + (torch.bmm(rotationMatrix,(rotated_coords[padding_mask][todoMask][existingRotation]-masked_translation).unsqueeze(-1)).squeeze(-1))


                #rotated_coords[todoMask] = fat

                rotate_coord_mask = padding_mask.clone()
                rotate_coord_mask[rotate_coord_mask == True]=todoMask

                rotated_coords = rotated_coords.masked_scatter(rotate_coord_mask.unsqueeze(-1).repeat(1,3),fat)

        torch.arange(rotated_coords.shape[0])
        rotated_coords_final = torch.full(coordinates.shape,float(PADDING_INDEX),device=coordinates.device)

        rotated_coords_final[ atom_description[:, hashings.atom_description_hash["batch"]].long(),coordsIndexingAtom] = rotated_coords
        return rotated_coords_final
