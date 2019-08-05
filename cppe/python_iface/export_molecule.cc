#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "../core/molecule.hh"

namespace py = pybind11;

void export_molecule(py::module& m) {
  // libcppe::Atom
  py::class_<libcppe::Atom> atom(m, "Atom", "Atom class");
  atom.def(py::init<int, double, double, double>())
        .def_readwrite("atomic_number", &libcppe::Atom::atomic_number)
        .def_readwrite("charge", &libcppe::Atom::charge)
        .def_readwrite("x", &libcppe::Atom::m_x)
        .def_readwrite("y", &libcppe::Atom::m_y)
        .def_readwrite("z", &libcppe::Atom::m_z)
        .def_property_readonly("position", &libcppe::Atom::get_pos);

  // libcppe::Molecule
  py::bind_vector<libcppe::Molecule>(m, "Molecule");
}