#pragma once

#include <Eigen/Core>

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
      std::cout << "remove_trace() not implemented for multipoles of order "
                << m_k << std::endl;
      std::cout << "Results from integral calculations could be wrong."
                << std::endl;
    }
  }

  std::vector<double> &get_values() { return m_values; }
  Eigen::VectorXd get_values_vec() {
    return Eigen::Map<Eigen::VectorXd>(m_values.data(), m_values.size());
    ;
  }
  unsigned m_k;
};

// Currently, only dipole-dipole polarizabilities are used
class Polarizability {
 private:
  std::vector<double> m_values;

 public:
  Polarizability(){};
  ~Polarizability(){};

  void add_value(double val) { m_values.push_back(val); }

  Eigen::VectorXd get_values_vec() {
    return Eigen::Map<Eigen::VectorXd>(m_values.data(), m_values.size());
  }
  std::vector<double> &get_values() { return m_values; }
};

class Potential {
 private:
  std::vector<Multipole> m_multipoles;
  std::vector<Polarizability> m_polarizabilities;
  // sites to exclude, 0-based index
  std::vector<int> m_exclusions;

 public:
  Potential(double x, double y, double z, int idx)
      : m_x(x), m_y(y), m_z(z), index(idx){};
  ~Potential(){};

  double m_x, m_y, m_z;
  int index;

  void add_multipole(Multipole mul) { m_multipoles.push_back(mul); }

  void add_polarizability(Polarizability pol) {
    m_polarizabilities.push_back(pol);
  }

  // 0-based!!!
  void add_exclusion(int excl) { m_exclusions.push_back(excl); }

  bool excludes_site(int other_site) {
    return (std::find(m_exclusions.begin(), m_exclusions.end(), other_site) !=
            m_exclusions.end());
  }

  std::vector<int> &get_exclusions() { return m_exclusions; }

  std::vector<Multipole> &get_multipoles() { return m_multipoles; }

  std::vector<Polarizability> &get_polarizabilities() {
    return m_polarizabilities;
  }

  bool is_polarizable() { return (m_polarizabilities.size() > 0); }

  Eigen::Vector3d get_site_position() { return Eigen::Vector3d(m_x, m_y, m_z); }
};

}  // namespace libcppe
