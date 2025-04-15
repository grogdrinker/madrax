#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  parse_foldx_params.py
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
import os,string
import numpy as np

letters={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ASN': 'N', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ALA': 'A', 'HIS': 'H', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'} #gli aminoacidi, che male qui non fanno
hashing_hybrid = {"NO_HYBRID":0,"SP2_N_H1":1,"SP2_N_H2":2,"SP2_N_ORB1":3,"SP2_N_ORB2":4,"SP3_N_H3":5,"SP3_O_H1ORB2":6,"SP2_O_ORB2":7}

path = os.path.dirname(os.path.abspath(__file__))

def read_h_positions(fil=path+"/../parameters/hbond_coords_params.txt"): ### ADD BACKBONE PARTNERS MANUALLY!!! ###
	diz={}
	for i in open(fil).readlines():
		l=i.split("\t")
		
		if not l[0] in letters and l[0]!="ooo":
			continue
			
		aa=l[0]
		if l[6]=="_" or l[6][0]!="H":
			continue
			
		if not aa in diz:
			diz[aa]={}
			
		if not l[1] in diz[aa]:
			diz[aa][l[1]]={}
			
		if  l[6]+"_0" in diz[aa][l[1]]:
			#continue

			#hlabel_addition[0]
			cont=1
			nameH = l[6] + "_" + str(cont)
			while nameH in diz[aa][l[1]]:
				nameH = l[6]+"_"+str(cont)
				cont+=1
			l[6]=nameH
		else:
			l[6] = l[6]+"_0"
		diz[aa][l[1]][l[6]]=(float(l[7]),float(l[8]),float(l[9]))
	
	for aa in letters.keys():
	
		if aa!="PRO":
			diz[aa]["N"]={}
			diz[aa]["N"]["HN_0"] = diz["ooo"]["N"]["HN_0"]


	del diz["ooo"]

			
		
	return diz


def read_hbond_params(fil=path+"/../parameters/hbond_params.txt"):
	level_hashing = {
					"LEVEL_N": 0,
					"LEVEL_A": 1,
					"LEVEL_O": -1,
					"LEVEL_B": 2,
					"LEVEL_G": 3,
					"LEVEL_D": 4,
					"LEVEL_E": 5,
					"LEVEL_Z": 6,
					"LEVEL_H": 7,
					"LEVEL_I": 8,
					"LEVEL_K": 9,
					}

	header=[]
	fin={}
	for l in open(fil,"r").readlines():
		if l[0]=="#":
			header+=[l.strip().replace("#","")]
		elif l.strip()=="":
			continue
		else:
			a=l.strip().replace('"',"").replace(' ',"").replace('\t',"").strip(",").split(",")

			assert len(a)==len(header)
			diz={}
			for k in range(len(a)):
				if header[k]=="aa":
					aa_ind = k

				elif header[k]=="atom":
					at_ind = k
				elif header[k] == "level":
					diz[header[k]] = level_hashing[a[k]]
				elif header[k]!="hybridation":
					try:
						diz[header[k]]=float(a[k])
					except:
						diz[header[k]]=a[k]

				else:
					diz[header[k]]=hashing_hybrid[a[k]]
					
			
			if not a[0] in fin:
				fin[a[aa_ind]]={}
			'''
			if a[aa_ind]=="H1S" and a[at_ind]=="NE2":
				print("H1S NE2",diz["donor"],diz["acceptor"])
			if a[aa_ind] == "H1S" and a[at_ind] == "ND1":
				print("H1S ND1",diz["donor"],diz["acceptor"])
			if a[aa_ind]=="H2S" and a[at_ind]=="NE2":
				print("asdasd", diz["donor"],diz["acceptor"])
			if a[aa_ind] == "H2S" and a[at_ind] == "ND1":
				print("H2S ND1",diz["donor"],diz["acceptor"])
			if a[aa_ind]=="HIS" and a[at_ind]=="NE2":
				print("hisNH2", diz["donor"],diz["acceptor"])
			if a[aa_ind]=="HIS" and a[at_ind]=="ND1":
				print("hisND1", diz["donor"],diz["acceptor"])
			'''
			fin[a[aa_ind]][a[at_ind]]=diz

	return fin
	
def read_hbond_partners(fil=path+"/../parameters/hbond_coords_params.txt"):
	diz={}
	for i in open(fil).readlines():
		l=i.split("\t")
		
		if not l[0] in letters:
			continue
			
		aa=l[0]
		if l[6]=="_":
			continue
			
		if not aa in diz:
			diz[aa]={}

		diz[aa][l[1]]=[[l[2],0],[l[4],0]]
		
		#### backbone manually ####

	for aa in letters.keys():

		if aa!="PRO":
			diz[aa]["N"] = [["C",-1],["O",-1]]
		diz[aa]["C"] = [["CA",0],["N",0]]
		diz[aa]["O"] = [["C",0],["CA",0]]

		#### virtual manually ####


	diz["ARG"]["ARG"] = [["NE",0],["NH1",0]]
	diz["ARG"]["CZ"] = [["NE", 0], ["NH1", 0]]
	diz["PHE"]["RC"] = [["CG", 0], ["CD1", 0]]
	diz["TYR"]["RC"] = [["CG", 0], ["CD1", 0]]

	return diz
	
def read_FO_positions(fil=path+"/../parameters/hbond_coords_params.txt"):
	diz={}
	hlabel_addition=string.ascii_lowercase
	for i in open(fil).readlines():
		l=i.split("\t")
		
		if not (l[0] in letters or l[0]=="ooo"):
			continue
			
		aa=l[0]

		if len(l[6])<2 or l[6][:2]!="FO":
			continue



		if not aa in diz:
			diz[aa]={}
			
		if not l[1] in diz[aa]:
			diz[aa][l[1]]={}
		l[6]=l[6].replace("FO","X")	# names are too long for pdb
		#if aa=="HIS" and l[1]=="ND1H2S":
		#	ads
		if l[6]+"_0" in diz[aa][l[1]]:
			#continue
			

			cont=1
			nameFO = l[6] + "_" + str(cont)
			while nameFO in diz[aa][l[1]]:
				nameFO = l[6]+"_"+str(cont)
				cont+=1
			l[6]=nameFO
		else:
			l[6] = l[6] +"_0"

		diz[aa][l[1]][l[6]]=(float(l[7]),float(l[8]),float(l[9]))
		#if aa=="HIS" and l[1]=="ND1H2S":

	for aa in letters.keys():
		if not aa in diz:
			diz[aa]={}
		
		diz[aa]["O"]={}
		diz[aa]["O"]["X1_0"] = diz["ooo"]["O"]["X1_0"]
		diz[aa]["O"]["X2_0"] = diz["ooo"]["O"]["X2_0"]
	

	del diz["ooo"]

	return diz

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
