#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "pybind_arma.h"

#include "../core/math.hh"
#include "../core/multipole.hh"

namespace py = pybind11;

void export_multipole(py::module &m) {
  m.def("prefactors", py::overload_cast<unsigned>(&libcppe::prefactors),
        "Prefactors for the multipole components (alphabetical order)");
  // libcppe::Multipole
  py::class_<libcppe::Multipole> mul(m, "Multipole");
  mul.def(py::init<unsigned>())
      .def_readwrite("k", &libcppe::Multipole::m_k)
      .def("get_values", &libcppe::Multipole::get_values_vec)
      .def("remove_trace", &libcppe::Multipole::remove_trace)
      .def("add_value", &libcppe::Multipole::add_value);

  // libcppe::Polarizability
  py::class_<libcppe::Polarizability> polarizability(m, "Polarizability");
  polarizability.def(py::init<>())
      .def("get_values", &libcppe::Polarizability::get_values_vec)
      .def("add_value", &libcppe::Polarizability::add_value);

  // libcppe::Potential
  py::class_<libcppe::Potential> pot(m, "Potential", "Potential (Site)");
  pot.def(py::init<double, double, double, int>())
      .def_readwrite("x", &libcppe::Potential::m_x, "x coordinate")
      .def_readwrite("y", &libcppe::Potential::m_y, "y coordinate")
      .def_readwrite("z", &libcppe::Potential::m_z, "z coordinate")
      .def("is_polarizable", &libcppe::Potential::is_polarizable)
      .def("get_site_position", &libcppe::Potential::get_site_position)
      .def("get_multipoles", &libcppe::Potential::get_multipoles)
      .def("get_polarizabilities", &libcppe::Potential::get_polarizabilities)
      .def_readwrite("index", &libcppe::Potential::index, "site index");
}