<!-- # CPPE -->
# <img src="https://gist.githubusercontent.com/maxscheurer/43b3dd040ea09ab06546bc6c2c771f56/raw/ced0c420e4840faf203dbca4d719f90cd66ca3fb/cppe_logo.png" height=150>


[![Build Status](https://travis-ci.com/maxscheurer/cppe.svg?branch=master)](https://travis-ci.com/maxscheurer/cppe)
[![Documentation Status](https://readthedocs.org/projects/cppe/badge/?version=latest)](https://cppe.readthedocs.io/en/latest/?badge=latest)


CPPE is an open-source, light-weight C++ and Python library for Polarizable Embedding (PE)<sup>1,2</sup>
calculations.
It provides an easy-to-use API to implement PE for ground-state self-consistent field (SCF)
calculations and post-SCF methods. A convenient Python interface is also available.

CPPE enables PE calculations in the following programs:
- [PySCF](https://github.com/pyscf/pyscf)
- [Psi4](https://github.com/psi4/psi4)
- [Q-Chem](https://www.q-chem.com)
- [VeloxChem](https://veloxchem.org)

__Examples__ for the open-source Python-driven programs can be found [here](https://github.com/maxscheurer/cppe_examples).

## Installation
The easiest way to install CPPE is via `pip`:
```
pip install cppe
```
Note that CPPE will be built from source and a C++14 compatible compiler is required (see below).
If the compiler supports OpenMP, it will also be enabled for the CPPE build.

Manual builds can be done using CMake by running
```
git clone https://github.com/maxscheurer/cppe
cd cppe; mkdir build; cd build
cmake ..
make -j4
```

Alternatively, CPPE can be installed using the `setup.py` script with
```
git clone https://github.com/maxscheurer/cppe
cd cppe
python setup.py install
```

### Python interface
If the Python interface should be built, specify the CMake option
`-DENABLE_PYTHON_INTERFACE=ON`. If `pybind11` is not installed, CMake
will automatically download `pybind11` and install it locally.
Installing through `setup.py` will always build the Python interface.

### Dependencies
- C++ 14 compiler
- Python >= 3.6 (interpreter and development packages)

### Tests
The tests can be run with
```
python setup.py build_ext -i; python setup.py test
```
for the `setup.py` build, or
```
source setup_environment.sh; py.test
```
for the CMake build.


## Citation
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3345696.svg)](https://doi.org/10.5281/zenodo.3345696)

The journal article describing CPPE can be found [here](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00758).


**CPPE: An Open-Source C++ and Python Library for Polarizable Embedding**</br>
Maximilian Scheurer, Peter Reinholdt, Erik Rosendahl Kjellgren, Jógvan Magnus Haugaard Olsen, Andreas Dreuw, and Jacob Kongsted;
_Journal of Chemical Theory and Computation_ 2019 15 (11), 6154-6163,
DOI: 10.1021/acs.jctc.9b00758



## Literature
<sup>1</sup> Olsen, J. M., Aidas, K., & Kongsted, J. (2010). Excited States in Solution through Polarizable Embedding. _J. Chem. Theory Comput._, **6** (12), 3721–3734. https://doi.org/10.1021/ct1003803

<sup>2</sup> Olsen, J. M. H., & Kongsted, J. (2011). Molecular Properties through Polarizable Embedding. _Advances in Quantum Chemistry_ (Vol. 61). https://doi.org/10.1016/B978-0-12-386013-2.00003-6
