#ifndef LIBCPPE_CORE_MULTIPOLE_H
#define LIBCPPE_CORE_MULTIPOLE_H

#include <vector>
#include <iostream>
#include <armadillo>

namespace libcppe {

class Multipole;
using Potential = std::vector<Multipole>;
  
class Multipole {
private:
  std::vector<double> m_values;

public:
  Multipole (unsigned k, double x, double y, double z) :
    m_k(k), m_x(x), m_y(y), m_z(z) {  };
  
  ~Multipole () {};
  
  // TODO: add check if too many values
  void add_value(double val) {
    m_values.push_back(val);
  }
  
  // TODO: check if required number of values in m_values
  void remove_trace() {
    double trace;
    if (m_k == 2) {
      trace = (m_values[0] + m_values[3] + m_values[5]) / 3.0;
      m_values[0] -= trace;
      m_values[3] -= trace;
      m_values[5] -= trace;
    } else if (m_k > 2) {
      std::cout << "remove_trace() not implemented for multipoles of order " << m_k << std::endl;
      std::cout << "Results from integral calculations could be wrong." << std::endl;
    }
  }
  
  arma::vec get_site_position() {
    arma::vec pos(3);
    pos[0] = m_x;
    pos[1] = m_y;
    pos[2] = m_z;
    return pos;
  }
  
  std::vector<double>& get_values() {return m_values;}
  arma::vec get_values_vec() {return arma::vec(m_values.data(), m_values.size());}
  unsigned m_k;
  double m_x, m_y, m_z;
  
};

}

#endif // LIBCPPE_CORE_MULTIPOLE_H