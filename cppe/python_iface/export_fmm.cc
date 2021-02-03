#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/fmm/calculate.hh"
#include "../core/fmm/tree.hh"

namespace py = pybind11;

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
}
