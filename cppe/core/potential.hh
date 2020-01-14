#pragma once

#include <Eigen/Core>

#include "math.hh"
#include <iostream>
#include <vector>

namespace libcppe {

class Multipole {
 private:
  std::vector<double> m_values;

 public:
  Multipole(unsigned k) : m_k(k){};
  ~Multipole(){};

  // TODO: add check if too many values
  void add_value(double val) { m_values.push_back(val); }

  // TODO: check if required number of values in m_values
  void remove_trace() {
    double trace;
    if (m_k == 2) {
      trace = (m_values[0] + m_values[3] + m_values[5]) / 3.0;
      m_values[0] -= trace;
      m_values[3] -= trace;
      m_values[5] -= trace;
    } else if (m_k > 2) {
      throw std::runtime_error(
            "remove_trace() not implemented for multipoles of order > 2");
    }
  }

  std::vector<double>& get_values() { return m_values; }
  Eigen::VectorXd get_values_vec() {
    return Eigen::Map<Eigen::VectorXd>(m_values.data(), m_values.size());
    ;
  }
  unsigned m_k;
};

// Currently, only dipole-dipole polarizability are used
class Polarizability {
 private:
  Eigen::Matrix3d m_values;

 public:
  Polarizability() = default;
  explicit Polarizability(std::vector<double> values) {
    if (values.size() != 6) {
      throw std::runtime_error(
            "Invalid number of elements to construct polarizability. Vector must have "
            "exactly 6 elements.");
    }
    m_values = triangle_to_mat(values);
  };

  void make_isotropic() {
    Eigen::Matrix3d new_values = Eigen::Matrix3d::Zero();
    double trace;
    trace            = get_isotropic_value();
    new_values(0, 0) = new_values(1, 1) = new_values(2, 2) = trace;
    m_values                                               = new_values;
  }

  double get_isotropic_value() { return m_values.trace() / 3.0; }

  Eigen::VectorXd get_values_vec() { return mat_to_triangle(m_values); }
  Eigen::Matrix3d& get_matrix() { return m_values; }
};

class Potential {
 private:
  std::vector<Multipole> m_multipoles;
  Polarizability m_polarizability;
  bool m_is_polarizable = false;
  // sites to exclude, 0-based index
  std::vector<int> m_exclusions;

 public:
  Potential(double x, double y, double z, int idx) : m_x(x), m_y(y), m_z(z), index(idx){};
  ~Potential(){};

  double m_x, m_y, m_z;
  int index;

  void add_multipole(Multipole mul) { m_multipoles.push_back(mul); }

  // 0-based!!!
  void add_exclusion(int excl) { m_exclusions.push_back(excl); }

  bool excludes_site(int other_site) {
    return (std::find(m_exclusions.begin(), m_exclusions.end(), other_site) !=
            m_exclusions.end());
  }

  std::vector<int>& get_exclusions() { return m_exclusions; }

  std::vector<Multipole>& get_multipoles() { return m_multipoles; }

  void set_polarizability(Polarizability pol) {
    m_polarizability = pol;
    m_is_polarizable = true;
  }

  Polarizability& get_polarizability() { return m_polarizability; }

  bool is_polarizable() { return m_is_polarizable; }

  Eigen::Vector3d get_site_position() { return Eigen::Vector3d(m_x, m_y, m_z); }
};

}  // namespace libcppe
