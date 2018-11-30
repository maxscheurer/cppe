#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cppe_state.hh"
#include "../core/pe_options.hh"

#include "pybind_arma.h"

namespace py = pybind11;

void export_state(py::module &m) {
  py::class_<libcppe::CppeState> cppe_state(m, "CppeState");
  cppe_state.def(py::init<libcppe::PeOptions, libcppe::Molecule>())
      .def("set_es_operator", &libcppe::CppeState::set_es_operator)
      .def("get_es_operator", &libcppe::CppeState::es_operator_copy)
      .def("update_energies", &libcppe::CppeState::update_energies,
           "Updates the PE electrostatic-electronic energy given a density",
           py::arg("D"))
      .def("set_potentials", &libcppe::CppeState::set_potentials)
      .def("get_current_energies", &libcppe::CppeState::get_current_energies)
      .def("calculate_static_energies_and_fields",
           &libcppe::CppeState::calculate_static_energies_and_fields)
      .def("get_induced_moments", &libcppe::CppeState::get_induced_moments)
      .def("update_induced_moments",
           &libcppe::CppeState::update_induced_moments)
      .def("get_static_fields", &libcppe::CppeState::get_static_fields)
      .def("print_summary", &libcppe::CppeState::print_summary)
      .def("get_potentials", &libcppe::CppeState::get_potentials)
      .def("get_polarizable_site_number",
           &libcppe::CppeState::get_polarizable_site_number);
}