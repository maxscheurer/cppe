#pragma once

#include <vector>

#include <Eigen/Core>

namespace libcppe {

struct Atom {
  int atomic_number;
  int charge;
  double m_x, m_y, m_z;
  Atom(int an) : atomic_number(an) { charge = an; }
  Atom(int an, double x, double y, double z) : atomic_number(an), m_x(x), m_y(y), m_z(z) {
    charge = an;
  }

  Eigen::Vector3d get_pos() { return Eigen::Vector3d(m_x, m_y, m_z); }
};

// Molecule is a slightly decorated std::vector
struct Molecule : std::vector<Atom> {
  Eigen::Vector3d get_atom_position(int atom) {
    if (this->size() <= atom) {
      throw std::out_of_range("Not enough atoms in Molecule.");
    }
    return (*this)[atom].get_pos();
  }

  ~Molecule()       = default;
  Molecule& operator=(const Molecule&) = default;
};

}  // namespace libcppe
