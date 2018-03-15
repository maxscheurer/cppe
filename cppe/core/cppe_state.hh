#ifndef INCLUDE_CPPE_STATE_H
#define INCLUDE_CPPE_STATE_H

#include <armadillo>
#include "pe_energies.hh"

namespace libcppe {
  
class CppeState {
private:
  unsigned m_scf_iteration;
  int m_nbas;
  arma::mat m_es_operator;
  arma::mat m_pol_operator;
  
  PeEnergy m_pe_energy;
  

public:
  // TODO: extend constructor
  CppeState (int nbas) { m_nbas = nbas; }
  ~CppeState () {};
  
  arma::mat& pol_operator() { return m_pol_operator; }
  arma::mat& es_operator() { return m_es_operator; }
  
  arma::mat pol_operator_copy() const { return m_pol_operator; }
  arma::mat es_operator_copy() const { return m_es_operator; }
  
  void set_es_operator(const arma::mat& es_operator) { m_es_operator = es_operator; }
  void set_pol_operator(const arma::mat& pol_operator) { m_pol_operator = pol_operator; }
  
  void next_iteration() { m_scf_iteration++; }
  unsigned get_current_iteration() { return m_scf_iteration; }
  
};

} /* libcppe */

#endif //INCLUDE_CPPE_STATE_H