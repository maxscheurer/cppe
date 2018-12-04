# CPPE

CPPE is an open-source, light-weight **C**-**P**lus-**P**lus library for **P**olarizable **E**mbedding (PE)<sup>1,2</sup>
calculations.
It provides an easy-to-use API to implement PE for ground-state self-consistent field
calculations and post-SCF methods. A convenient Python interface is also available.

For testing purposes, CPPE is also interfaced to the original PE implementation, namely
[PElib](https://gitlab.com/pe-software/pelib-public), used in the DALTON program package.

<!-- CPPE is currently implemented in the Q-Chem program package for PE-SCF
and PE-ADC calculations <sup>3</sup>, and the open-source
packages [Psi4](http://psicode.org) and [pyscf](https://github.com/pyscf/pyscf).
The latter implementation makes use of the Python interface. -->

## Installation
CPPE needs to be built from sources, i.e.,
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

### PElib interface
To build PElib (together with gen1int and the appropriate interface), set the following option:
`-DENABLE_PELIB=ON`.


### Dependencies
- C++ 11 compiler
- [armadillo](http://arma.sourceforge.net/docs.html) (CMake will at least try to install armadillo from sources if not installed)
- Python 3.6 (interpreter and development packages)


## Literature
<sup>1</sup> Olsen, J. M., Aidas, K., & Kongsted, J. (2010). Excited States in Solution through Polarizable Embedding. _J. Chem. Theory Comput._, **6** (12), 3721–3734.

<sup>2</sup> Olsen, J. M. H., & Kongsted, J. (2011). Molecular Properties through Polarizable Embedding. _Advances in Quantum Chemistry_ (Vol. 61). https://doi.org/10.1016/B978-0-12-386013-2.00003-6