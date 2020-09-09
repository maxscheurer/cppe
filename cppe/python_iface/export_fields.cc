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

static std::shared_ptr<libcppe::MultipoleFields> _init_multipole_fields(
      std::vector<libcppe::Potential> potentials, py::dict py_options) {
  libcppe::PeOptions options = _dict_to_options(py_options);
  return std::make_shared<libcppe::MultipoleFields>(potentials, options);
}

void export_fields(py::module& m) {
  py::class_<libcppe::NuclearFields> nuc_fields(m, "NuclearFields",
                                                "Electric fields created by nuclei");
  nuc_fields.def(py::init<libcppe::Molecule, std::vector<libcppe::Potential>>())
        .def("compute", &libcppe::NuclearFields::compute);

  py::class_<libcppe::MultipoleFields, std::shared_ptr<libcppe::MultipoleFields>>
        mul_fields(m, "MultipoleFields", "Electric fields created by multipoles");
  mul_fields.def(py::init(&_init_multipole_fields))
        .def("compute", &libcppe::MultipoleFields::compute)
        .def("compute_legacy", &libcppe::MultipoleFields::compute_legacy)
        .def("compute_tree", &libcppe::MultipoleFields::compute_tree);

  py::class_<libcppe::InducedMoments, std::shared_ptr<libcppe::InducedMoments>>
        ind_moments(m, "InducedMoments");
  ind_moments.def(py::init(&_init_indmom))
        .def("compute",
             py::overload_cast<Eigen::VectorXd&, bool>(&libcppe::InducedMoments::compute),
             "Compute the induced moments solving the classical response "
             "equation",
             py::arg("total_fields"), py::arg("make_guess"))
        .def("compute_cg",
             py::overload_cast<Eigen::VectorXd&, bool>(
                   &libcppe::InducedMoments::compute_cg),
             "Compute the induced moments solving the classical response "
             "equation using a CG solver",
             py::arg("total_fields"), py::arg("make_guess"));

  m.def("get_polarizable_sites", &libcppe::get_polarizable_sites);

  py::class_<libcppe::BMatrix, std::shared_ptr<libcppe::BMatrix>> bmatrix(m, "BMatrix");
  bmatrix.def(py::init(&_init_bmatrix))
        .def("direct_inverse", &libcppe::BMatrix::direct_inverse)
        .def("to_dense_matrix", &libcppe::BMatrix::to_dense_matrix)
        .def("apply", &libcppe::BMatrix::apply)
        .def("apply_direct", &libcppe::BMatrix::apply_direct)
        .def("apply_fast_summation", &libcppe::BMatrix::apply_fast_summation)
        .def("apply_diagonal", &libcppe::BMatrix::apply_diagonal)
        .def("apply_diagonal_inverse", &libcppe::BMatrix::apply_diagonal_inverse);
  bmatrix.def_property_readonly("exclusions", &libcppe::BMatrix::get_exclusions);

  m.def("multipole_derivative", &libcppe::multipole_derivative);
}
