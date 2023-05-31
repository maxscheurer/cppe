#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../core/tensors.hh"

namespace py = pybind11;

using namespace libcppe::tensors_recursive;
using namespace libcppe::tensors;

using Tfun      = std::function<Eigen::VectorXd(const Eigen::Vector3d&)>;
using Tfun_damp = std::function<Eigen::VectorXd(const Eigen::Vector3d&, double)>;

void export_tensors(py::module& m) {
  py::module ten = m.def_submodule("tensors", "Tensor submodule");
  ten.def("T_recursive", &T_recursive, py::arg("k"), py::arg("Rij"),
          py::arg("damping_factor") = 0.0);
  ten.def("xyz2idx", &xyz2idx);

  py::bind_vector<std::vector<Tfun>>(ten, "Tensor functions.");
  py::bind_vector<std::vector<Tfun_damp>>(ten, "Damped tensor functions.");
  ten.def("T0", &T0);
  ten.def("T1", &T1);
  ten.def("T2", &T2);
  ten.def("T3", &T3);
  ten.def("T4", &T4);
  ten.def("T5", &T5);
  ten.def("T6", &T6);
  ten.attr("T") = libcppe::tensors::T;
  ten.def("T0_damp_thole", &T0_damp_thole);
  ten.def("T0_damp_thole", &T1_damp_thole);
  ten.def("T2_damp_thole", &T2_damp_thole);
  ten.def("T3_damp_thole", &T3_damp_thole);
  ten.def("T4_damp_thole", &T4_damp_thole);
  ten.def("T5_damp_thole", &T5_damp_thole);
  ten.def("T6_damp_thole", &T6_damp_thole);
  ten.attr("T_damp_thole") = libcppe::tensors::T_damp_thole;
}