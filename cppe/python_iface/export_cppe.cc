#include "../metadata.hh"
#include <pybind11/pybind11.h>

namespace py = pybind11;

void export_molecule(py::module&);
void export_multipole(py::module&);
void export_utils(py::module&);
void export_state(py::module&);
void export_fields(py::module&);
void export_math(py::module&);
void export_tensors(py::module&);

PYBIND11_MODULE(cppe, cppe) {
  cppe.doc()                    = "Python interface for CPPE";
  cppe.attr("__version__")      = libcppe::version::version_string();
  cppe.attr("__build_type__")   = libcppe::version::is_debug() ? "Debug" : "Release";
  cppe.attr("__authors__")      = libcppe::__authors__();
  cppe.attr("__contributors__") = libcppe::__contributors__();
  cppe.attr("__email__")        = libcppe::__email__();

  export_molecule(cppe);
  export_multipole(cppe);
  export_utils(cppe);
  export_state(cppe);
  export_fields(cppe);
  export_math(cppe);
  export_tensors(cppe);
}