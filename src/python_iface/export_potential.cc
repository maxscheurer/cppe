#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../core/math.hh"
#include "../core/potential.hh"

namespace py = pybind11;

void export_multipole(py::module& m) {
  // libcppe::Multipole
  py::class_<libcppe::Multipole> mul(m, "Multipole");
  mul.def(py::init<unsigned>())
        .def_readwrite("k", &libcppe::Multipole::m_k)
        .def_property_readonly("values", &libcppe::Multipole::get_values_vec)
        .def("remove_trace", &libcppe::Multipole::remove_trace)
        .def("add_value", &libcppe::Multipole::add_value);

  // libcppe::Polarizability
  py::class_<libcppe::Polarizability> polarizability(m, "Polarizability");
  polarizability.def(py::init<>())
        .def(py::init<std::vector<double>>())
        .def_property_readonly("values", &libcppe::Polarizability::get_values_vec)
        .def_property_readonly("isotropic_value",
                               &libcppe::Polarizability::get_isotropic_value)
        .def("make_isotropic", &libcppe::Polarizability::make_isotropic);

  // libcppe::Potential
  py::class_<libcppe::Potential> pot(m, "Potential", "Potential (Site)");
  pot.def(py::init<double, double, double, std::string, int>())
        .def_readwrite("x", &libcppe::Potential::m_x, "x coordinate")
        .def_readwrite("y", &libcppe::Potential::m_y, "y coordinate")
        .def_readwrite("z", &libcppe::Potential::m_z, "z coordinate")
        .def_readwrite("element", &libcppe::Potential::m_element, "element")
        .def("add_multipole", &libcppe::Potential::add_multipole)
        .def("add_exclusion", &libcppe::Potential::add_exclusion)
        .def("set_polarizability", &libcppe::Potential::set_polarizability)
        .def_property_readonly("is_polarizable", &libcppe::Potential::is_polarizable)
        .def("excludes_site", &libcppe::Potential::excludes_site)
        .def_property_readonly("exclusions", &libcppe::Potential::get_exclusions)
        .def_property_readonly("position", &libcppe::Potential::get_site_position)
        .def_property_readonly("multipoles", &libcppe::Potential::get_multipoles)
        .def_property_readonly("polarizability", &libcppe::Potential::get_polarizability)
        .def_readwrite("index", &libcppe::Potential::index, "site index");
}
