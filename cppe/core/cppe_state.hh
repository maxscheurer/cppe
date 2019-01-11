#ifndef INCLUDE_CPPE_STATE_H
#define INCLUDE_CPPE_STATE_H

#include <iostream>

#include <Eigen/Core>

#include "molecule.hh"
#include "multipole.hh"
#include "pe_energies.hh"
#include "pe_options.hh"

namespace libcppe {

class CppeState {
 private:
  Eigen::MatrixXd m_es_operator;  //!< PE electrostatics operator

  PeEnergy m_pe_energy;  //!< PE Energy Container

  // Molecule and Potentials
  Molecule m_mol;                       //!< core region molecule
  std::vector<Potential> m_potentials;  //!< vector with all site potentials

  size_t m_polarizable_sites;  //!< number of polarizable sites
  // Static Fields
  Eigen::VectorXd m_nuc_fields;  //!< electric fields from nuclei
  Eigen::VectorXd
      m_multipole_fields;  //!< electric fields from multipole moments

  Eigen::VectorXd m_induced_moments;  //!< Vector with induced moments

  PeOptions m_options;

  std::ostream& m_output_stream = std::cout;  //!< Output stream for printing

  bool m_make_guess = true;

 public:
  CppeState(){};
  CppeState(PeOptions options, Molecule mol, std::ostream& = std::cout);
  ~CppeState(){};

  void set_options(PeOptions options) { m_options = options; }
  void set_molecule(Molecule mol) { m_mol = mol; }

  void set_potentials(std::vector<Potential> potentials);
  std::vector<Potential> get_potentials() { return m_potentials; }

  PeEnergy& get_energies() { return m_pe_energy; }
  void set_energies(PeEnergy energy) { m_pe_energy = energy; }
  void calculate_static_energies_and_fields();

  std::vector<double> get_induced_moments() const {
    return std::vector<double>(
        m_induced_moments.data(),
        m_induced_moments.data() + m_induced_moments.size());
  }

  Eigen::VectorXd get_induced_moments_vec() const { return m_induced_moments; }

  void update_induced_moments(Eigen::VectorXd elec_fields,
                              bool elec_only = false);

  size_t get_polarizable_site_number() { return m_polarizable_sites; }

  std::vector<double> get_static_fields() {
    Eigen::VectorXd static_fields = m_nuc_fields + m_multipole_fields;
    return std::vector<double>(static_fields.data(),
                               static_fields.data() + static_fields.size());
  }

  // TODO: summary as string
  void print_summary();
};

}  // namespace libcppe

#endif  // INCLUDE_CPPE_STATE_H
