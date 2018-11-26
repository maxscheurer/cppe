#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "pybind_arma.h"

#include "../core/cppe_state.hh"
#include "../core/molecule.hh"

// namespace pybind11 { namespace detail {
    // template<>
    // struct type_caster<std::vector<libcppe::Atom, std::allocator<libcppe::Atom> >>
    //  : list_caster<std::vector<libcppe::Atom, std::allocator<libcppe::Atom>>, libcppe::Atom> { };
//     template <> struct type_caster<std::vector<libcppe::Atom>> : list_caster<std::vector<libcppe::Atom>, libcppe::Atom> { };
// }
// }

namespace py = pybind11;

void export_core(py::module& m) {
    py::class_<libcppe::Atom> atom(m, "Atom", "Atom class");
    atom.def(py::init<int, double, double, double>())
        .def_readwrite("atomic_number", &libcppe::Atom::atomic_number)
        .def_readwrite("charge", &libcppe::Atom::charge)
        .def_readwrite("x", &libcppe::Atom::m_x)
        .def_readwrite("y", &libcppe::Atom::m_y)
        .def_readwrite("z", &libcppe::Atom::m_z)
        .def("get_position", &libcppe::Atom::get_pos);
    
    // py::class_<std::vector<libcppe::Atom>> atomlist(m, "AtomList");
    // py::bind_vector<std::vector<libcppe::Atom>>(m, "AtomList");
    // py::class_<libcppe::Molecule, std::vector<libcppe::Atom>> mol(m, "Molecule", "Molecule class");
    py::bind_vector<libcppe::Molecule>(m, "Molecule");
    
    py::class_<libcppe::Potential> pot(m, "Potential", "Potential (Site)");
    pot.def(py::init<double, double, double, int>())
        .def_readwrite("x", &libcppe::Potential::m_x, "x coordinate")
        .def_readwrite("y", &libcppe::Potential::m_y, "y coordinate")
        .def_readwrite("z", &libcppe::Potential::m_z, "z coordinate")
        .def_readwrite("index", &libcppe::Potential::index, "site index");
}