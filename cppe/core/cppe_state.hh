#ifndef INCLUDE_CPPE_STATE_H
#define INCLUDE_CPPE_STATE_H

#include <armadillo>

#include "pe_energies.hh"
#include "molecule.hh"
#include "multipole.hh"
#include "pe_options.hh"

namespace libcppe {

class CppeState {
private:
  arma::mat m_es_operator; //!< PE electrostatics operator

  PeEnergy m_pe_energy; //!< PE Energy Container

  // Molecule and Potentials
  Molecule m_mol; //!< core region molecule
  std::vector<Potential> m_potentials; //!< vector with all site potentials

  size_t m_polarizable_sites; //!< number of polarizable sites
  // Static Fields
  arma::vec m_nuc_fields; //!< electric fields from nuclei
  arma::vec m_multipole_fields; //!< electric fields from multipole moments

  arma::vec m_induced_moments; //!< Vector with induced moments

  PeOptions m_options;


public:
  // TODO: extend constructor
  CppeState(PeOptions options);
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

  arma::vec get_static_fields() { return (m_nuc_fields + m_multipole_fields); }

  void print_summary(std::ostream& output_stream = std::cout);

};

} /* libcppe */

#endif //INCLUDE_CPPE_STATE_H
