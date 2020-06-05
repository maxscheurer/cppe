#include <string>
#include <unordered_map>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../core/cppe_state.hh"

namespace py = pybind11;

using UnorderedMapStringDouble = std::unordered_map<std::string, double>;
using UnorderedMapMapStringDouble =
      std::unordered_map<std::string, std::unordered_map<std::string, double>>;

PYBIND11_MAKE_OPAQUE(UnorderedMapStringDouble);
PYBIND11_MAKE_OPAQUE(UnorderedMapMapStringDouble);

libcppe::PeOptions _dict_to_options(py::dict py_options) {
  libcppe::PeOptions options;
  py::list valid_keys = py::cast(libcppe::valid_option_keys);
  for (const auto& entry : py_options) {
    std::string key = entry.first.cast<std::string>();
    auto value      = entry.second;
    if (!valid_keys.contains(entry.first)) {
      throw std::invalid_argument("Option key \'" + key + "\' is invalid.");
    }
    if (key == "potfile") {
      options.potfile = value.cast<std::string>();
    } else if (key == "iso_pol") {
      options.iso_pol = value.cast<bool>();
    } else if (key == "induced_thresh") {
      options.induced_thresh = value.cast<double>();
    } else if (key == "maxiter") {
      options.maxiter = value.cast<int>();
    } else if (key == "damp_induced") {
      options.damp_induced = value.cast<bool>();
    } else if (key == "damp_multipole") {
      options.damp_multipole = value.cast<bool>();
    } else if (key == "damping_factor_induced") {
      options.damping_factor_induced = value.cast<double>();
    } else if (key == "damping_factor_multipole") {
      options.damping_factor_multipole = value.cast<double>();
    } else if (key == "pe_border") {
      options.pe_border = value.cast<bool>();
    } else if (key == "border_type") {
      options.border_type = value.cast<std::string>();
    } else if (key == "border_rmin") {
      options.border_rmin = value.cast<double>();
    } else if (key == "border_nredist") {
      options.border_nredist = value.cast<int>();
    } else if (key == "border_redist_order") {
      options.border_redist_order = value.cast<int>();
    } else if (key == "border_redist_pol") {
      options.border_redist_pol = value.cast<bool>();
    }
  }
  return options;
}

static std::shared_ptr<libcppe::CppeState> _init_state(py::dict py_options,
                                                       libcppe::Molecule mol,
                                                       libcppe::PrintCallback callback) {
  libcppe::PeOptions options = _dict_to_options(py_options);
  return std::make_shared<libcppe::CppeState>(options, mol, callback);
}

static py::dict _options_to_dict(libcppe::CppeState state) {
  libcppe::PeOptions self = state.get_options();
  py::list valid_keys     = py::cast(libcppe::valid_option_keys);
  py::dict ret;
  for (const auto& pykey : valid_keys) {
    std::string key = pykey.cast<std::string>();
    if (key == "potfile") {
      ret[pykey] = py::cast(self.potfile);
    } else if (key == "iso_pol") {
      ret[pykey] = py::cast(self.iso_pol);
    } else if (key == "induced_thresh") {
      ret[pykey] = py::cast(self.induced_thresh);
    } else if (key == "maxiter") {
      ret[pykey] = py::cast(self.maxiter);
    } else if (key == "damp_induced") {
      ret[pykey] = py::cast(self.damp_induced);
    } else if (key == "damp_multipole") {
      ret[pykey] = py::cast(self.damp_multipole);
    } else if (key == "damping_factor_induced") {
      ret[pykey] = py::cast(self.damping_factor_induced);
    } else if (key == "damping_factor_multipole") {
      ret[pykey] = py::cast(self.damping_factor_multipole);
    } else if (key == "pe_border") {
      ret[pykey] = py::cast(self.pe_border);
    } else if (key == "border_type") {
      ret[pykey] = py::cast(self.border_type);
    } else if (key == "border_rmin") {
      ret[pykey] = py::cast(self.border_rmin);
    } else if (key == "border_nredist") {
      ret[pykey] = py::cast(self.border_nredist);
    } else if (key == "border_redist_order") {
      ret[pykey] = py::cast(self.border_redist_order);
    } else if (key == "border_redist_pol") {
      ret[pykey] = py::cast(self.border_redist_pol);
    }
  }
  return ret;
}

void export_state(py::module& m) {
  py::bind_map<UnorderedMapStringDouble>(m, "UnorderedMapStringDouble");
  py::bind_map<UnorderedMapMapStringDouble>(m, "UnorderedMapMapStringDouble");

  py::class_<libcppe::CppeState, std::shared_ptr<libcppe::CppeState>> cppe_state(
        m, "CppeState");
  cppe_state
        .def(py::init(&_init_state),
             "Create a CppeState using a dictionary with options, molecule, and print "
             "callback",
             py::arg("options") = py::dict(), py::arg("molecule") = libcppe::Molecule(),
             py::arg("printer") = libcppe::default_printer)
        .def("set_potentials", &libcppe::CppeState::set_potentials)
        .def("calculate_static_energies_and_fields",
             &libcppe::CppeState::calculate_static_energies_and_fields)
        .def("get_induced_moments", &libcppe::CppeState::get_induced_moments)
        .def_readwrite("energies", &libcppe::CppeState::m_pe_energy)
        .def_property_readonly("total_energy", &libcppe::CppeState::get_total_energy)
        .def("update_induced_moments", &libcppe::CppeState::update_induced_moments)
        .def_property_readonly("static_fields", &libcppe::CppeState::get_static_fields)
        .def_property_readonly("nuclear_fields", &libcppe::CppeState::get_nuclear_fields)
        .def_property_readonly("multipole_fields",
                               &libcppe::CppeState::get_multipole_fields)
        .def_property_readonly("summary_string",
                               &libcppe::CppeState::get_energy_summary_string)
        .def_property_readonly("potentials", &libcppe::CppeState::get_potentials)
        .def("get_polarizable_site_number",
             &libcppe::CppeState::get_polarizable_site_number)
        .def_property_readonly("options", &_options_to_dict);
  m.attr("valid_option_keys") = py::cast(libcppe::valid_option_keys);
}