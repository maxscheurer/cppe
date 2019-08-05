#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../core/molecule.hh"
#include "../core/pot_manipulation.hh"
#include "../core/potential.hh"
#include "../core/potfile_reader.hh"

namespace py = pybind11;

void export_utils(py::module& m) {
  py::class_<libcppe::PotfileReader> potfile_reader(m, "PotfileReader",
                                                    "Reader for potential files");
  potfile_reader.def(py::init<std::string>())
        .def("read", &libcppe::PotfileReader::read, "Read the potential file");

  py::class_<libcppe::PotManipulator> pot_manip(m, "PotManipulator",
                                                "Manipulator for potentials");
  pot_manip.def(py::init<std::vector<libcppe::Potential>, libcppe::Molecule>());
}