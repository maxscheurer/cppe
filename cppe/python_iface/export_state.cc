#include <string>
#include <unordered_map>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../core/cppe_state.hh"
#include "../core/pe_options.hh"

namespace py = pybind11;

using UnorderedMapStringDouble = std::unordered_map<std::string, double>;
using UnorderedMapMapStringDouble =
      std::unordered_map<std::string, std::unordered_map<std::string, double>>;

PYBIND11_MAKE_OPAQUE(UnorderedMapStringDouble);
PYBIND11_MAKE_OPAQUE(UnorderedMapMapStringDouble);

void export_state(py::module& m) {
  py::bind_map<UnorderedMapStringDouble>(m, "UnorderedMapStringDouble");
  py::bind_map<UnorderedMapMapStringDouble>(m, "UnorderedMapMapStringDouble");

  py::class_<libcppe::CppeState> cppe_state(m, "CppeState");
  cppe_state
        .def(py::init<libcppe::PeOptions, libcppe::Molecule, libcppe::PrintCallback>(),
             "Constructor", py::arg("options") = libcppe::PeOptions(),
             py::arg("molecule") = libcppe::Molecule(),
             py::arg("printer")  = libcppe::default_printer)
        .def("set_potentials", &libcppe::CppeState::set_potentials)
        .def("calculate_static_energies_and_fields",
             &libcppe::CppeState::calculate_static_energies_and_fields)
        .def("get_induced_moments", &libcppe::CppeState::get_induced_moments)
        .def_readwrite("energies", &libcppe::CppeState::m_pe_energy)
        .def_property_readonly("total_energy", &libcppe::CppeState::get_total_energy)
        .def("update_induced_moments", &libcppe::CppeState::update_induced_moments)
        .def("get_static_fields", &libcppe::CppeState::get_static_fields)
        .def_property_readonly("summary_string",
                               &libcppe::CppeState::get_energy_summary_string)
        .def_property_readonly("potentials", &libcppe::CppeState::get_potentials)
        .def("get_polarizable_site_number",
             &libcppe::CppeState::get_polarizable_site_number);
}