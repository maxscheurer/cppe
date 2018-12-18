#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/math.hh"

#include "pybind_arma.h"

namespace py = pybind11;

using namespace libcppe;

void export_math(py::module &m) {
  m.def("smat_vec", &smat_vec);
  m.def("multipole_derivative", &multipole_derivative);
  m.def("T", &T);
  m.def("Tk_tensor", &Tk_tensor);
  m.def("Tk_coefficients", &Tk_coefficients);
  m.def("xyz2idx", &xyz2idx);
  m.def("factorial", &factorial);
  m.def("make_df", &make_df);
  m.def("trinom", &trinom);
  m.def("symmetry_factors", &symmetry_factors);
  m.def("prefactors", py::overload_cast<unsigned>(&prefactors),
        "Prefactors for the multipole components (alphabetical order)");
  m.def("prefactors_nuclei", &prefactors_nuclei);
  m.def("multipole_components", &multipole_components);
}