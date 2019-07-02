# The CPPE C++ layer

This page contains a reference of the *CPPE* C++ library
and its classes and functions generated automatically
from the *CPPE* source code.

## Quantum region molecule

### Atom container
```eval_rst
This class represents an atom in the quantum region molecule

.. doxygenstruct:: libcppe::Atom
   :members:
   :undoc-members:

```

### Molecule container
```eval_rst
This class contains data of the quantum region molecule.
It is simply a decorated :cpp:class:`std::vector` of :cpp:class:`libcppe::Atom`.

.. doxygenstruct:: libcppe::Molecule
   :members:
   :undoc-members:

```
<!-- 
```eval_rst
This category lists the *adccore* functionality,
which imports the data from the :cpp:class:`adcc::HartreeFockSolution_i`
interface into the :cpp:class:`adcc::ReferenceState`
for internal use by the library.
Important classes in the process are :cpp:class:`adcc::MoSpaces`,
which collects information about the occupied and virtual
orbital spaces, and :cpp:class:`adcc::MoIndexTranslation`,
which maps orbitals indices between the ordering used by adccore
and the one used by the SCF program.

.. doxygengroup:: ReferenceObjects
   :members:
   :content-only:

``` -->

