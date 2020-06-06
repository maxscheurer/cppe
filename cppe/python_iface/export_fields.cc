#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/electric_fields.hh"
#include "../core/molecule.hh"
#include "../core/pe_options.hh"

#include "../core/bmatrix.hh"

namespace py = pybind11;

libcppe::PeOptions _dict_to_options(py::dict);

static std::shared_ptr<libcppe::BMatrix> _init_bmatrix(
      std::vector<libcppe::Potential> potentials, py::dict py_options) {
  libcppe::PeOptions options = _dict_to_options(py_options);
  return std::make_shared<libcppe::BMatrix>(potentials, options);
}

static std::shared_ptr<libcppe::InducedMoments> _init_indmom(
      std::vector<libcppe::Potential> potentials, py::dict py_options) {
  libcppe::PeOptions options = _dict_to_options(py_options);
  return std::make_shared<libcppe::InducedMoments>(potentials, options);
}

void export_fields(py::module& m) {
  py::class_<libcppe::NuclearFields> nuc_fields(m, "NuclearFields",
                                                "Electric fields created by nuclei");
  nuc_fields.def(py::init<libcppe::Molecule, std::vector<libcppe::Potential>>())
        .def("compute", &libcppe::NuclearFields::compute);

  py::class_<libcppe::MultipoleFields> mul_fields(
        m, "MultipoleFields", "Electric fields created by multipoles");
  mul_fields.def(py::init<std::vector<libcppe::Potential>, libcppe::PeOptions>())
        .def("compute", &libcppe::MultipoleFields::compute);

  py::class_<libcppe::InducedMoments, std::shared_ptr<libcppe::InducedMoments>>
        ind_moments(m, "InducedMoments");
  ind_moments.def(py::init(&_init_indmom))
        .def("compute",
             py::overload_cast<Eigen::VectorXd&, bool>(&libcppe::InducedMoments::compute),
             "Compute the induced moments solving the classical response "
             "equation",
             py::arg("total_fields"), py::arg("make_guess"))
        .def("compute_cg", &libcppe::InducedMoments::compute_cg);

  m.def("get_polarizable_sites", &libcppe::get_polarizable_sites);

  py::class_<libcppe::BMatrix, std::shared_ptr<libcppe::BMatrix>> bmatrix(m, "BMatrix");
  bmatrix.def(py::init(&_init_bmatrix))
        .def("direct_inverse", &libcppe::BMatrix::direct_inverse)
        .def("compute_apply", &libcppe::BMatrix::compute_apply)
        .def("compute_apply_slice", &libcppe::BMatrix::compute_apply_slice);

  m.def("multipole_derivative", &libcppe::multipole_derivative);
}
