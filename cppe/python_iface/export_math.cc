#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/math.hh"

namespace py = pybind11;

using namespace libcppe;

void export_math(py::module& m) {
  // m.def("smat_vec", &smat_vec);
  m.def("triangle_to_mat", &triangle_to_mat<std::vector<double>>);
  m.def("mat_to_triangle", &mat_to_triangle);
  m.def("multipole_derivative", &multipole_derivative);
  m.def("T", &T, py::arg("Rij"), py::arg("x"), py::arg("y"), py::arg("z"),
        py::arg("Cijn"), py::arg("damping_factor") = 0.0, py::arg("alpha_i") = 0.0,
        py::arg("alpha_j") = 0.0);
  m.def("Tk_tensor", &Tk_tensor, py::arg("k"), py::arg("Rij"), py::arg("Tk_coeffs"),
        py::arg("damping_factor") = 0.0, py::arg("alpha_i") = 0.0,
        py::arg("alpha_j") = 0.0);
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
