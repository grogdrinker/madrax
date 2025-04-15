![logo](docs/_static/logo_small.png)

---

## About MadraX

MadraX is a Pytorch-based Force Field designed to calculate the stability of a protein or complex in a differentiable way. 

MadraX is capable of interacting in an end-to-end way with neural networks and it can record the gradient via the autograd function of Pytorch

If you use MadraX in your research, please consider citing:


## Installation

For installation, refer to the official documentation: https://madrax.readthedocs.io/install.html

We recommend using MadraX with Python 3.7 3.8, 3.9, 3.10, 3.11 or 3.12. 
Package installation should only take a few minutes with any of these methods (conda, pip, source).

### Installing MadraX with pip:

We suggest to install pytorch separately following the instructions from https://pytorch.org/get-started/locally/

```sh
python -m pip install "madrax @git+https://github.com/grogdrinker/madrax/"
```

### Installing MadraX from source:

If you want to install MadraX from this repository, you need to install the dependencies first.
First, install [PyTorch](https://pytorch.org/get-started/locally/). The library is currently compatible with PyTorch versions between 1.8 and 1.13. We will continue to update MadraX to be compatible with the latest version of PyTorch.

Finally, you can clone the repository with the following command:

```sh
git clone https://github.com/grogdrinker/madrax/
```

## Documentation

The documentation for MadraX is available https://madrax.readthedocs.io

## Quickstart

We provide a quickstart example to show how to easily run madrax on PDB structures.
The quickstart is available https://madrax.readthedocs.io/quickstart.html#


## Tutorials

Tutorials for madrax are available https://madrax.readthedocs.io/tutorials/index.html

We recommend the tutorials to be run on a machine with a GPU, as they will take longer when run on a CPU machine.

## Help

For bug reports, features addition and technical questions please contact gabriele.orlando@umontpellier.fr
