#include <pybind11/pybind11.h>

namespace py = pybind11;

void export_molecule(py::module &);
void export_multipole(py::module &);
void export_options(py::module &);
void export_utils(py::module &);
void export_state(py::module &);
void export_fields(py::module &);
void export_math(py::module &);

PYBIND11_MODULE(cppe, cppe) {
  cppe.doc() = "Python interface for CPPE";
  export_molecule(cppe);
  export_multipole(cppe);
  export_options(cppe);
  export_utils(cppe);
  export_state(cppe);
  export_fields(cppe);
  export_math(cppe);
}