#include <pybind11/pybind11.h>

namespace py = pybind11;

void export_molecule(py::module &);
void export_multipole(py::module &);
void export_utils(py::module &);

PYBIND11_MODULE(pycppe, pycppe) {
  pycppe.doc() = "Python interface for CPPE";
  export_molecule(pycppe);
  export_multipole(pycppe);
  export_utils(pycppe);
}