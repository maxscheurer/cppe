#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/math.hh"

namespace py = pybind11;

using namespace libcppe;

void export_math(py::module& m) {
  m.def("triangle_to_mat", &triangle_to_mat<std::vector<double>>);
  m.def("mat_to_triangle", &mat_to_triangle);
  m.def("factorial", &factorial);
  m.def("make_df", &make_df);
  m.def("trinom", &trinom);
  m.def("symmetry_factors", &symmetry_factors);
  m.def("prefactors", py::overload_cast<unsigned>(&prefactors),
        "Prefactors for the multipole components (alphabetical order)");
  m.def("prefactors_nuclei", &prefactors_nuclei);
  m.def("multipole_components", &multipole_components);
}
