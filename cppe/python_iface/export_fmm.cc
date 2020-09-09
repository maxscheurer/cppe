#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/fmm/calculate.hh"
#include "../core/fmm/tree.hh"

namespace py = pybind11;

static std::shared_ptr<Tree<1, 3>> __build_tree(
      py::array_t<double, py::array::c_style> positions,
      py::array_t<double, py::array::c_style> sources, size_t nparticles, size_t ncrit,
      size_t order, double theta, std::vector<std::vector<int>> exclusion_lists) {
  return build_shared_tree<1, 3>(positions.mutable_data(), sources.mutable_data(),
                                 nparticles, ncrit, order, theta, exclusion_lists);
}

static void __compute_field_fmm(const std::shared_ptr<Tree<1, 3>>& self,
                                py::array_t<double, py::array::c_style> field) {
  self->compute_field_fmm(field.mutable_data());
}

static void __compute_field_exact(const std::shared_ptr<Tree<1, 3>>& self,
                                  py::array_t<double, py::array::c_style> field) {
  self->compute_field_exact(field.mutable_data());
}

static void __set_sources(const std::shared_ptr<Tree<1, 3>>& self,
                          py::array_t<double, py::array::c_style> sources) {
  self->set_sources(sources.mutable_data());
}

void export_fmm(py::module& m) {
  py::module fmm = m.def_submodule("fmm", "Module for fast multipole methods.");
  py::class_<Tree<1, 3>, std::shared_ptr<Tree<1, 3>>> tree(fmm, "Tree",
                                                           "Tree class for FMM/BH.");
  tree.def("compute_field_fmm", &__compute_field_fmm)
        .def("compute_field_exact", &__compute_field_exact)
        .def("set_sources", &__set_sources);
  fmm.def("build_tree", &__build_tree);
}
