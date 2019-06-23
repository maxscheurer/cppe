#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <variant>
#include <unordered_map>

#include "../core/pe_options.hh"

using option_t = std::variant<int, double, bool, std::string>;
using OptionMap = std::unordered_map<std::string, option_t>;

PYBIND11_MAKE_OPAQUE(OptionMap);

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
  py::bind_map<OptionMap>(m, "OptionMap");
  py::class_<libcppe::PeOptions> pe_options(m, "PeOptions",
                                            "Options for CPPE library");
  pe_options.def(py::init<>())
      .def_readwrite("options", &libcppe::PeOptions::options)
      .def_readwrite("potfile", &libcppe::PeOptions::potfile)
      .def_readwrite("iso_pol", &libcppe::PeOptions::iso_pol)

      .def_readwrite("induced_thresh", &libcppe::PeOptions::induced_thresh)
      .def_readwrite("do_diis", &libcppe::PeOptions::do_diis)
      .def_readwrite("diis_maxiter", &libcppe::PeOptions::diis_maxiter)

      .def_readwrite("pe_border", &libcppe::PeOptions::pe_border)
      .def_readwrite("border_options", &libcppe::PeOptions::border_options);
}