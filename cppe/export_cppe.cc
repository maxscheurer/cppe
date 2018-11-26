#include <pybind11/pybind11.h>


namespace py = pybind11;

PYBIND11_MODULE(example,m)
{
  m.doc() = "pybind11 example plugin";

  m.def("mul", &mul, py::call_guard<py::scoped_ostream_redirect,
                     py::scoped_estream_redirect>());
  m.def("set_element", &set_element);
  m.def("matmul", &matmul);
}