from madrax.ForceField import ForceField # the main MadraX module
from madrax import utils,dataStructures # the MadraX utility module
import time,os,urllib # some standard python modules we are going to us

os.mkdir('exampleStructures')
urllib.request.urlretrieve('http://files.rcsb.org/download/101M.pdb', 'exampleStructures/101m.pdb')
urllib.request.urlretrieve('http://files.rcsb.org/download/5BMZ.pdb', 'exampleStructures/5bmz.pdb')
urllib.request.urlretrieve('http://files.rcsb.org/download/5BOX.pdb', 'exampleStructures/5box.pdb')




device = "cpu"
coords, atnames, pdbNames = utils.parsePDB("exampleStructures/")
info_tensors = dataStructures.create_info_tensors(atnames,device=device)
forceField_Obj = ForceField(device=device)
energy = forceField_Obj(coords.to(device), info_tensors)
print("en",energy.sum())
