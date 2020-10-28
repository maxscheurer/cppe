#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>

#include "cppe_state.hh"
#include "electric_fields.hh"
#include "multipole_expansion.hh"
#include "pot_manipulation.hh"
#include "potfile_reader.hh"

namespace libcppe {

CppeState::CppeState(const PeOptions& options, Molecule mol, PrintCallback printer)
      : m_options(options), m_mol(mol), m_printer(printer) {
  m_pe_energy                       = PeEnergy{};
  std::vector<Potential> potentials = PotfileReader(m_options.potfile).read();
  auto manip                        = PotManipulator(potentials, m_mol);
  manip.set_print_callback(m_printer);
  potentials = manip.manipulate(m_options);
  set_potentials(std::move(potentials));
  // create empty energy container
  std::map<std::string, std::vector<std::string>> energy_terms{
        {"Electrostatic", {"Electronic", "Nuclear", "Multipoles"}},
        {"Polarization", {"Electronic", "Nuclear", "Multipoles"}}};
  for (auto it : energy_terms) {
    for (auto nm : it.second) {
      m_pe_energy[it.first][nm] = 0.0;
    }
  }
}

void CppeState::set_potentials(std::vector<Potential> potentials) {
  m_potentials        = potentials;
  m_polarizable_sites = std::count_if(m_potentials.begin(), m_potentials.end(),
                                      [](Potential p) { return p.is_polarizable(); });
  m_induced_moments   = Eigen::VectorXd::Zero(m_polarizable_sites * 3);

  m_positions             = Eigen::MatrixXd::Zero(m_potentials.size(), 3);
  m_positions_polarizable = Eigen::MatrixXd::Zero(m_polarizable_sites, 3);

  for (int i = 0; i < m_potentials.size(); ++i) {
    m_positions(i, 0) = m_potentials[i].m_x;
    m_positions(i, 1) = m_potentials[i].m_y;
    m_positions(i, 2) = m_potentials[i].m_z;
  }

  auto m_potentials_polarizable = get_polarizable_sites(m_potentials);
  for (int i = 0; i < m_polarizable_sites; ++i) {
    m_positions_polarizable(i, 0) = m_potentials_polarizable[i].m_x;
    m_positions_polarizable(i, 1) = m_potentials_polarizable[i].m_y;
    m_positions_polarizable(i, 2) = m_potentials_polarizable[i].m_z;
  }
}

void CppeState::calculate_static_energies_and_fields() {
  // Electrostatic energy (nuclei-multipoles)
  MultipoleExpansion mexp(m_mol, m_potentials);
  double nuc_mul_energy                   = mexp.calculate_interaction_energy();
  m_pe_energy["Electrostatic"]["Nuclear"] = nuc_mul_energy;

  // Calculate static fields

  // Nuclear fields
  NuclearFields nfields(m_mol, m_potentials);
  m_nuc_fields = nfields.compute();
  // Multipole fields
  MultipoleFields mul_fields(m_potentials, m_options);
  m_multipole_fields = mul_fields.compute();
}

void CppeState::update_induced_moments(Eigen::VectorXd elec_fields, bool elec_only) {
  Eigen::VectorXd tmp_total_fields = Eigen::VectorXd::Zero(m_polarizable_sites * 3);
  if (elec_only) {
    tmp_total_fields = elec_fields;
  } else {
    tmp_total_fields = elec_fields + m_nuc_fields + m_multipole_fields;
  }
  InducedMoments ind(m_potentials, m_options);
  ind.set_print_callback(m_printer);
  if (m_options.summation_induced_fields == "direct") {
    ind.compute(tmp_total_fields, m_induced_moments, m_make_guess);
  } else {
    m_induced_moments = ind.compute_cg(tmp_total_fields, m_induced_moments, m_make_guess);
  }
  if (m_make_guess) {
    m_make_guess = false;
  }

  if (elec_only) {
    double epol_elec                          = -0.5 * m_induced_moments.dot(elec_fields);
    m_pe_energy["Polarization"]["Electronic"] = epol_elec;
  } else {
    double epol_elec       = -0.5 * m_induced_moments.dot(elec_fields);
    double epol_nuclear    = -0.5 * m_induced_moments.dot(m_nuc_fields);
    double epol_multipoles = -0.5 * m_induced_moments.dot(m_multipole_fields);

    m_pe_energy["Polarization"]["Electronic"] = epol_elec;
    m_pe_energy["Polarization"]["Nuclear"]    = epol_nuclear;
    m_pe_energy["Polarization"]["Multipoles"] = epol_multipoles;
  }
}

Eigen::MatrixXd CppeState::induced_moments_eef() {
  InducedMoments ind(m_potentials, m_options);
  ind.set_print_callback(m_printer);

  Eigen::MatrixXd ret = Eigen::MatrixXd::Zero(m_polarizable_sites * 3, 3);
  Eigen::MatrixXd Fdn = Eigen::MatrixXd::Zero(m_polarizable_sites * 3, 3);
  for (int s = 0; s < m_polarizable_sites; ++s) {
    int l         = 3 * s;
    Fdn(l, 0)     = 1;
    Fdn(l + 1, 1) = 1;
    Fdn(l + 2, 2) = 1;
  }
  for (int a = 0; a < 3; ++a) {
    Eigen::VectorXd ind_mom = Eigen::VectorXd::Zero(m_polarizable_sites * 3);
    if (m_options.summation_induced_fields == "direct") {
      ind.compute(Fdn.col(a), ind_mom, true);
      ret.col(a) = ind_mom;
    } else {
      ret.col(a) = ind.compute_cg(Fdn.col(a), ind_mom, true);
    }
  }
  return ret;
}

double CppeState::get_total_energy_for_category(std::string category) {
  auto energy_cat   = m_pe_energy[category];
  double acc_energy = std::accumulate(
        energy_cat.begin(), energy_cat.end(), 0.0,
        [](double value, const std::unordered_map<std::string, double>::value_type& p) {
          return value + p.second;
        });
  return acc_energy;
}

double CppeState::get_total_energy() {
  return get_total_energy_for_category("Electrostatic") +
         get_total_energy_for_category("Polarization");
}

std::string CppeState::get_energy_summary_string() {
  size_t w = 30;
  std::string off(3, ' ');
  std::string off2(6, ' ');
  std::stringstream sstream;
  sstream << std::string(2 * w + 10, '-') << std::endl;

  sstream << "Polarizable Embedding Summary:";
  sstream << std::setprecision(12) << std::endl << std::endl;
  sstream << off << "Electrostatics:" << std::endl;
  sstream << off2 << "Electronic:" << std::string(w - 11, ' ')
          << m_pe_energy["Electrostatic"]["Electronic"] << std::endl;
  sstream << off2 << "Nuclear:" << std::string(w - 8, ' ')
          << m_pe_energy["Electrostatic"]["Nuclear"] << std::endl;
  sstream << off2 << "Multipole:" << std::string(w - 10, ' ')
          << m_pe_energy["Electrostatic"]["Multipoles"] << std::endl;
  double electrostatic_energy = get_total_energy_for_category("Electrostatic");
  sstream << off2 << "Total:" << std::string(w - 6, ' ') << electrostatic_energy
          << std::endl;
  sstream << std::endl;
  sstream << off << "Polarization:" << std::endl;
  sstream << off2 << "Electronic:" << std::string(w - 11, ' ')
          << m_pe_energy["Polarization"]["Electronic"] << std::endl;
  sstream << off2 << "Nuclear:" << std::string(w - 8, ' ')
          << m_pe_energy["Polarization"]["Nuclear"] << std::endl;
  sstream << off2 << "Multipole:" << std::string(w - 10, ' ')
          << m_pe_energy["Polarization"]["Multipoles"] << std::endl;
  double polarization_energy = get_total_energy_for_category("Polarization");
  sstream << off2 << "Total:" << std::string(w - 6, ' ')
          << get_total_energy_for_category("Polarization") << std::endl;
  sstream << std::endl;
  sstream << off << "Total Energy:" << std::string(w - 10, ' ')
          << electrostatic_energy + polarization_energy << std::endl;
  sstream << std::string(2 * w + 10, '-') << std::endl << std::endl;
  return sstream.str();
}

}  // namespace libcppe
