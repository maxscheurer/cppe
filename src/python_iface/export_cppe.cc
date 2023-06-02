#include <pybind11/pybind11.h>

namespace py = pybind11;

void export_molecule(py::module&);
void export_multipole(py::module&);
void export_utils(py::module&);
void export_state(py::module&);
void export_fields(py::module&);
void export_math(py::module&);
void export_tensors(py::module&);
void export_fmm(py::module&);

PYBIND11_MODULE(pycppe, cppe) {
  cppe.doc()                  = "Python interface for CPPE";
  cppe.attr("__build_type__") = [] {
#ifdef NDEBUG
    return "Release";
#else
    return "Debug";
#endif
  }();
  cppe.attr("__parallel__") = [] {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
  }();

  export_molecule(cppe);
  export_multipole(cppe);
  export_utils(cppe);
  export_state(cppe);
  export_fields(cppe);
  export_math(cppe);
  export_tensors(cppe);
  export_fmm(cppe);
}
