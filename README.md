<!-- # CPPE -->
# <img src="https://gist.githubusercontent.com/maxscheurer/43b3dd040ea09ab06546bc6c2c771f56/raw/ced0c420e4840faf203dbca4d719f90cd66ca3fb/cppe_logo.png" height=150>


[![Build Status](https://travis-ci.org/maxscheurer/cppe.svg?branch=master)](https://travis-ci.org/maxscheurer/cppe)
[![Documentation Status](https://readthedocs.org/projects/cppe/badge/?version=latest)](https://cppe.readthedocs.io/en/latest/?badge=latest)
![](https://img.shields.io/badge/ChemRxiv-%20preprint-critical?link=https://chemrxiv.org/articles/CPPE_An_Open-Source_C_and_Python_Library_for_Polarizable_Embedding/8949101)


CPPE is an open-source, light-weight **C**-**P**lus-**P**lus library for **P**olarizable **E**mbedding (PE)<sup>1,2</sup>
calculations.
It provides an easy-to-use API to implement PE for ground-state self-consistent field
calculations and post-SCF methods. A convenient Python interface is also available.

<!-- CPPE is currently implemented in the Q-Chem program package for PE-SCF
and PE-ADC calculations <sup>3</sup>, and the open-source
packages [Psi4](http://psicode.org) and [pyscf](https://github.com/pyscf/pyscf).
The latter implementation makes use of the Python interface. -->

## Build
If used in a C++ package, CPPE needs to be built from sources, i.e.,
```
git clone https://github.com/maxscheurer/cppe
cd cppe; mkdir build; cd build
cmake ..
make -j4
```

### Python interface
If the Python interface should be built, specify the CMake option
`-DENABLE_PYTHON_INTERFACE=ON`. If `pybind11` is not installed, CMake
will automatically download `pybind11` and install it locally.

### Dependencies
- C++ 14 compiler
- Python 3.6 (interpreter and development packages)

<!-- ## Installation via pip
will be possible in the future. -->

## Citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3345696.svg)](https://doi.org/10.5281/zenodo.3345696)

A preprint of the paper describing the CPPE library can be found [here](https://chemrxiv.org/articles/CPPE_An_Open-Source_C_and_Python_Library_for_Polarizable_Embedding/8949101).




## Literature
<sup>1</sup> Olsen, J. M., Aidas, K., & Kongsted, J. (2010). Excited States in Solution through Polarizable Embedding. _J. Chem. Theory Comput._, **6** (12), 3721–3734. https://doi.org/10.1021/ct1003803

<sup>2</sup> Olsen, J. M. H., & Kongsted, J. (2011). Molecular Properties through Polarizable Embedding. _Advances in Quantum Chemistry_ (Vol. 61). https://doi.org/10.1016/B978-0-12-386013-2.00003-6
