#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/cppe_state.hh"
#include "../core/pe_energies.hh"
#include "../core/pe_options.hh"

namespace py = pybind11;

void export_state(py::module &m) {
  py::class_<libcppe::CppeState> cppe_state(m, "CppeState");
  cppe_state.def(py::init<libcppe::PeOptions, libcppe::Molecule>())
      .def("set_potentials", &libcppe::CppeState::set_potentials)
      .def("get_energies", &libcppe::CppeState::get_energies)
      .def("set_energies", &libcppe::CppeState::set_energies)
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

  py::class_<libcppe::PeEnergy> pe_energy(m, "PeEnergy");
  pe_energy.def(py::init<>())
      .def("get_total_energy", &libcppe::PeEnergy::get_total_energy)
      .def("set", &libcppe::PeEnergy::set)
      .def("get", &libcppe::PeEnergy::get);
}