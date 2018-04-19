#ifndef INCLUDE_CPPE_STATE_H
#define INCLUDE_CPPE_STATE_H

#include <armadillo>

#include "pe_energies.hh"
#include "molecule.hh"
#include "multipole.hh"

namespace libcppe {

class CppeState {
private:
  // Operators
  arma::mat m_es_operator;
  // arma::mat m_pol_operator;

  // PE Energy Container
  PeEnergy m_pe_energy;

  // Molecule and Potentials
  Molecule m_mol;
  std::vector<Potential> m_potentials;

  size_t m_polarizable_sites;
  // Static Fields
  arma::vec m_nuc_fields;
  arma::vec m_multipole_fields;

  // Induced Moments
  arma::vec m_induced_moments;


public:
  // TODO: extend constructor
  CppeState();
  ~CppeState() {};

  // arma::mat pol_operator_copy() const { return m_pol_operator; }
  arma::mat es_operator_copy() const { return m_es_operator; }

  void set_es_operator(arma::mat es_operator) { m_es_operator = es_operator; }
  // void set_pol_operator(arma::mat pol_operator) { m_pol_operator = pol_operator; }

  void update_energies(arma::mat P);

  void set_molecule(Molecule mol) { m_mol = mol; }
  void set_potentials(std::vector<Potential> potentials);

  PeEnergy get_current_energies() const { return m_pe_energy; }
  void calculate_static_energies_and_fields();

  arma::vec get_induced_moments() const { return m_induced_moments; }
  void update_induced_moments(arma::vec elec_fields, int iteration, bool elec_only = false);

  size_t get_polarizable_site_number() { return m_polarizable_sites; }

  void print_summary();

};

} /* libcppe */

#endif //INCLUDE_CPPE_STATE_H
