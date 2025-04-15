==========
Quickstart
==========
This quickstart example shows how to easily run MadraX on a protein.

Importing libraries 
^^^^^^^^^^^^^^^^^^^
You can use the following code to import all the necessary libraries used in this example:

.. code-block:: python

   from madrax.ForceField import ForceField # the main MadraX module
   from madrax import utils,dataStructures # the MadraX utility module
   import time,os,urllib # some standard python modules we are going to use

Collecting the data
^^^^^^^^^^^^^^^^^^^

We now need some protein structures to test. You can download them from the PDB website (https://www.rcsb.org/). In this case, we will get them using urllib, a standard and preinstalled library of python. This will allow a painless data collection with a simple copy-paste of the code. Lets fetch a couple of structures:

.. code-block:: python

   os.mkdir('exampleStructures')
   urllib.request.urlretrieve('http://files.rcsb.org/download/101M.pdb', 'exampleStructures/101m.pdb')
   urllib.request.urlretrieve('http://files.rcsb.org/download/5BMZ.pdb', 'exampleStructures/5bmz.pdb')
   urllib.request.urlretrieve('http://files.rcsb.org/download/5BOX.pdb', 'exampleStructures/5box.pdb')

Parsing the structures
^^^^^^^^^^^^^^^^^^^^^^

We first set the device we want to run MadraX on. This follows the typical pytorch way to deal with devices.
For simplicity, lets assign the device we want to use to a variable called "device"

.. code-block:: python

   device = "cpu" #or if you wanna use GPU you can write device = "cuda"


Next step is to parse the information of the structures. In order to do so, we can use the utility module of MadraX:

.. code-block:: python

   coords, atnames, pdbNames = utils.parsePDB("exampleStructures/") # get coordinates and atom names

We also need to "tensorize" the atom names. This step organizes the data and it needs to be run once only. If the **order,type and number** of atom does not change you don't need to run it again.

.. code-block:: python

   info_tensors = dataStructures.create_info_tensors(atnames,device=device)

Running MadraX
^^^^^^^^^^^^^^

then we can create the main MadraX object

.. code-block:: python

   forceField_Obj = ForceField(device=device)

And we can calculate the energy of the proteins simply providing the coordinates and the previously calculated info_tensors

.. code-block:: python

   energy = forceField_Obj(coords.to(device), info_tensors)

The output energy is structured as following:
.. code-block:: python

   dim 0: batch --> protein number
   dim 1: chain --> chain index
   dim 3: residue --> residue index
   dim 4: mutant --> mutant index. Ignore it for now, we will use it later
   dim 5: energy types structured as following:

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
