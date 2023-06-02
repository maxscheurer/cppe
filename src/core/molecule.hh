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

  Eigen::Vector3d get_position() { return Eigen::Vector3d(m_x, m_y, m_z); }
};

using Molecule = std::vector<Atom>;

}  // namespace libcppe
