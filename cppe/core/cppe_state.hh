#pragma once

#include <iomanip>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

#include "math.hh"
#include "molecule.hh"
#include "pe_options.hh"
#include "potential.hh"

namespace libcppe {

using PeEnergy = std::unordered_map<std::string, std::unordered_map<std::string, double>>;

using PrintCallback = std::function<void(std::string)>;

const PrintCallback default_printer = [](std::string print_string) {
  std::cout << std::setprecision(12) << print_string << std::endl;
};

class CppeState {
 private:
  // Molecule and Potentials
  Molecule m_mol;                       //!< core region molecule
  std::vector<Potential> m_potentials;  //!< vector with all site potentials

  size_t m_polarizable_sites;  //!< number of polarizable sites
  // Static Fields
  Eigen::VectorXd m_nuc_fields;        //!< electric fields from nuclei
  Eigen::VectorXd m_multipole_fields;  //!< electric fields from multipole moments

  Eigen::VectorXd m_induced_moments;  //!< Vector with induced moments

  PeOptions m_options;

  bool m_make_guess = true;
  PrintCallback m_printer;

 public:
  CppeState(){};
  explicit CppeState(PeOptions options, Molecule mol,
                     PrintCallback printer = default_printer);
  ~CppeState(){};

  void set_options(PeOptions options) { m_options = options; }
  PeOptions get_options() { return m_options; }
  void set_molecule(Molecule mol) { m_mol = mol; }

  void set_potentials(std::vector<Potential> potentials);
  std::vector<Potential> get_potentials() { return m_potentials; }

  void calculate_static_energies_and_fields();

  std::vector<double> get_induced_moments() const {
    return std::vector<double>(m_induced_moments.data(),
                               m_induced_moments.data() + m_induced_moments.size());
  }

  PeEnergy m_pe_energy;  //!< PE Energy Container
  double get_total_energy_for_category(std::string);
  double get_total_energy();

  Eigen::VectorXd get_induced_moments_vec() const { return m_induced_moments; }

  void update_induced_moments(Eigen::VectorXd elec_fields, bool elec_only = false);

  size_t get_polarizable_site_number() { return m_polarizable_sites; }

  Eigen::VectorXd get_static_fields() { return m_nuc_fields + m_multipole_fields; }

  Eigen::VectorXd get_nuclear_fields() { return m_nuc_fields; }

  Eigen::VectorXd get_multipole_fields() { return m_multipole_fields; }

  std::string get_energy_summary_string();
};

}  // namespace libcppe
