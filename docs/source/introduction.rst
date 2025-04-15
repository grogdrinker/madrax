Introduction
============

What is a force field?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A force field is  physical description of the forces involved in (de)stabilizing the conformation of a macromolecule. A forcefield, therefore, is a tool that takes as input a set of atom coordinates, alogn with their atom types and characteristics, and provides he Gibbs energy (\Delta G) of the complex. There are a lot of forcefield currently available for the scientific community, such as Rosetta or FoldX
 
Machine learning in structural biology
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Machine learning in structural biology made a massive step forward with the introduction of modern artificial intelligence.
These techniques extract all the information required to solve a problem from data, not requiring any a-priori model.

The Problem
^^^^^^^^^^^^

While learning everything from data is definetly an advantage for systems which are way to difficult to handle analytically (i.e. protein folding), this complete dependance by experimental information can become a problem in topics with scarce data. In this case, the artificial intelligence don't have enough data to learn even the known physical rules that we find in the system.


Our Solution
^^^^^^^^^^^^

The solution we propose is to hard code the physical rules inside the neural network in order to release it from the taks to learn from the data something we are already able to model analytically

What is MadraX?
^^^^^^^^^^^^^^^

MadraX is a differentiable force field implemented as a pytorch module inspired by FoldX.
Given a set of atoms, it therefore defines the energies that stabilize of destabilize the protein or complex conformation.

Why is Differentiability Important?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MadraX is a fully differentiable pytorch module. This means the gradient tracked by pytorch can flow through it. And what does it mean? It means that MadraX can be included INSIDE a neural network, and pytorch will not complain about it. The network will be trained normally, but now, at every backpropagation step, it will also get an energy estimation.
Let's take for example Alphafold: MatdraX could be included in the training loop, as an additional loss function.

To our knowledge, Madrax is the first and only forcefield that can do such a thing.

Long story short, if you want a force field to interact with a network, affecting the actual evolution of it, you need differentiability

Serendipity of Differentiability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A very welcome side effect of differentiablity is that pytorch, while mainly used for neural networks developement, can also be used to calculate gradients respect an arbitary tensor.

If we calculate the gradient of the energy of a protein (the sum of the energies provided by MadraX) respect to a set of transformation matrices that are applied on the protein itself (i.e. rotation matrices that are applied on torsion angles), we get how we should rotate EVERY torsion angle in order to reduce the energy.

We can use pytorch optimizers that are usually used for neural network training (such as Adam), to relax a protein.
Additionally, we turned protein optimization into a standard gradient minimization task and every future update on pytorch optimizers, such as the developement of new and more effective optimizers, will be directly applicable to proteins as well thanks to MadraX.

Isn't it nice? For once is the University that makes use of the developements coming from companies!

What does MadraX mean?
^^^^^^^^^^^^^^^^^^^^^^

Madrac (singular for madracs, read mardax) is the Friulan name (a Rhaeto-Romance language) for hierophis viridiflavus carbonarius, a very common snake in Friuli, a mountainous region in north-east of Italy full of very grumpy people.
