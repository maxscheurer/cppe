#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/pe_options.hh"

namespace py = pybind11;

void export_options(py::module &m) {
  // libcppe::BorderOptions
  py::class_<libcppe::BorderOptions> pe_border_options(
      m, "PeBorderOptions", "Border Options for CPPE library");
  py::enum_<libcppe::BorderType>(pe_border_options, "BorderType")
      .value("rem", libcppe::BorderType::rem)
      .value("redist", libcppe::BorderType::redist);
  pe_border_options.def(py::init<>())
      .def_readwrite("border_type", &libcppe::BorderOptions::border_type)
      .def_readwrite("rmin", &libcppe::BorderOptions::rmin)
      .def_readwrite("nredist", &libcppe::BorderOptions::nredist)
      .def_readwrite("redist_order", &libcppe::BorderOptions::redist_order)
      .def_readwrite("redist_pol", &libcppe::BorderOptions::redist_pol);

  // libcppe::PeOptions
  py::class_<libcppe::PeOptions> pe_options(m, "PeOptions",
                                            "Options for CPPE library");
  pe_options.def(py::init<>())
      .def_readwrite("potfile", &libcppe::PeOptions::potfile)
      .def_readwrite("iso_pol", &libcppe::PeOptions::iso_pol)

      .def_readwrite("induced_thresh", &libcppe::PeOptions::induced_thresh)
      .def_readwrite("do_diis", &libcppe::PeOptions::do_diis)
      .def_readwrite("diis_start_norm", &libcppe::PeOptions::diis_start_norm)
      .def_readwrite("maxiter", &libcppe::PeOptions::maxiter)

      .def_readwrite("pe_border", &libcppe::PeOptions::pe_border)
      .def_readwrite("border_options", &libcppe::PeOptions::border_options);
}