# CPPE

`cppe` is a C++ and Python library for the Polarizable Embedding (PE) model.
It implements the necessary routines to implement the PE model in any
quantum-chemical program package (host program) _via_ a light-weight, minimal interface.

The interface can be implemented either on the C++ level or through Python.
In general, several components of the host program need to be accessed in the
interface to `cppe`, though `cppe` as such is completely host-program-agnostic.

From a technical point of view, the only requirement to a quantum-chemical
program to run PE calculations with `cppe` are integrals of the form

```eval_rst
.. math::

   t_{p q}^{(k)}\left(\underline{\mathbf{R}}_{s}\right)=-\int \phi_{p}^{*}(\underline{\mathbf{r}})\left(\frac{\partial^{|k|}}{\partial \underline{\mathbf{r}}^{k}} \frac{1}{\left|\underline{\mathbf{r}}-\underline{\mathbf{R}}_{s}\right|}\right) \phi_{q}(\underline{\mathbf{r}}) \mathrm{d} \underline{\mathbf{r}}
   
evaluated at the fixed coordinates :math:`\underline{\mathbf{R}}_{s}` in the environment.
```

## Contents
```eval_rst
.. toctree::
   :maxdepth: 2

   cppecore
   publications
   license

* :ref:`genindex`
```
