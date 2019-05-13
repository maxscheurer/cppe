#include <algorithm>
#include <iomanip>
#include <iostream>

#include "cppe_state.hh"

#include "electric_fields.hh"
#include "molecule.hh"
#include "multipole.hh"
#include "multipole_expansion.hh"
#include "pe_energies.hh"
#include "pe_options.hh"

#include "pot_manipulation.hh"
#include "potfile_reader.hh"

namespace libcppe {

CppeState::CppeState(PeOptions options, Molecule mol,
                     std::ostream &output_stream)
    : m_options(options), m_output_stream(output_stream), m_mol(mol) {
  m_pe_energy = PeEnergy{};
  std::vector<Potential> potentials = PotfileReader(m_options.potfile).read();
  potentials = PotManipulator(potentials, m_mol).manipulate(m_options);
  set_potentials(potentials);
}

void CppeState::set_potentials(std::vector<Potential> potentials) {
  m_potentials = potentials;
  m_polarizable_sites =
      std::count_if(m_potentials.begin(), m_potentials.end(),
                    [](Potential p) { return p.is_polarizable(); });
  m_induced_moments = Eigen::VectorXd::Zero(m_polarizable_sites * 3);
}

void CppeState::calculate_static_energies_and_fields() {
  // Electrostatic energy (nuclei-multipoles)
  MultipoleExpansion mexp(m_mol, m_potentials);
  double nuc_mul_energy = mexp.calculate_interaction_energy();
  m_pe_energy.set("Electrostatic/Nuclear", nuc_mul_energy);

  // Calculate static fields

  // Nuclear fields
  NuclearFields nfields(m_mol, m_potentials);
  m_nuc_fields = nfields.compute();
  // Multipole fields
  MultipoleFields mul_fields(m_potentials);
  m_multipole_fields = mul_fields.compute();
}

void CppeState::update_induced_moments(Eigen::VectorXd elec_fields,
                                       bool elec_only) {
  Eigen::VectorXd tmp_total_fields =
      Eigen::VectorXd::Zero(m_polarizable_sites * 3);
  if (elec_only) {
    tmp_total_fields = elec_fields;
  } else {
    tmp_total_fields = elec_fields + m_nuc_fields + m_multipole_fields;
  }
  InducedMoments ind(m_potentials, m_options);
  ind.compute(tmp_total_fields, m_induced_moments, m_make_guess,
              m_output_stream);
  if (m_make_guess) {
    m_make_guess = false;
  }

  if (elec_only) {
    double epol_elec = -0.5 * m_induced_moments.dot(elec_fields);
    // std::cout << std::setprecision(15) << "epol_elec = " << epol_elec <<
    // std::endl;
    m_pe_energy.set("Polarization/Electronic", epol_elec);
  } else {
    double epol_elec = -0.5 * m_induced_moments.dot(elec_fields);
    double epol_nuclear = -0.5 * m_induced_moments.dot(m_nuc_fields);
    double epol_multipoles = -0.5 * m_induced_moments.dot(m_multipole_fields);

    // std::cout << "epol_elec = " << epol_elec << std::endl;
    // std::cout << "epol_nuclear = " << epol_nuclear << std::endl;
    // std::cout << "epol_multipoles = " << epol_multipoles << std::endl;

    m_pe_energy.set("Polarization/Electronic", epol_elec);
    m_pe_energy.set("Polarization/Nuclear", epol_nuclear);
    m_pe_energy.set("Polarization/Multipoles", epol_multipoles);
  }
}

void CppeState::print_summary() {
  size_t w = 30;
  std::string off(3, ' ');
  std::string off2(6, ' ');
  // std::cout << off << "Corrected Excitation Energy (" << energy_name << "): "
  // << std::string(w - 52, ' '); std::cout << std::setw(11) << (state.energy +
  // ene) * conversion::au2ev << " eV" << std::endl;
  m_output_stream << std::string(2 * w + 10, '-') << std::endl;
  // m_output_stream << "__________      .__               .__              ___.
  // .__          \n"
  // "\\______   \\____ |  | _____ _______|__|____________ \\_ |__ |  |   ____
  // \n" " |     ___/  _ \\|  | \\__  \\\\_  __ \\  \\___   /\\__  \\ | __ \\|
  // | _/ __ \\ \n" " |    |  (  <_> )  |__/ __ \\|  | \\/  |/    /  / __ \\|
  // \\_\\ \\  |_\\  ___/ \n" " |____|   \\____/|____(____  /__|  |__/_____
  // \\(____  /___  /____/\\___  >\n"
  // "___________      ___.     \\/       .___  .___.__  \\/    \\/          \\/
  // \n"
  // "\\_   _____/ _____\\_ |__   ____   __| _/__| _/|__| ____    ____ \n" " |
  // __)_ /     \\| __ \\_/ __ \\ / __ |/ __ | |  |/    \\  / ___\\        \n"
  // " |        \\  Y Y  \\ \\_\\ \\  ___// /_/ / /_/ | |  |   |  \\/ /_/  > \n"
  // "/_______  /__|_|  /___  /\\___  >____ \\____ | |__|___|  /\\___  / \n" "
  // _________     \\/    \\/     \\/     \\/    \\/         \\//_____/ \n" " /
  // _____/__ __  _____   _____ _____ _______ ___.__.                  \n" "
  // \\_____  \\|  |  \\/     \\ /     \\\\__  \\\\_  __ <   |  | \n" " / \\  |
  // /  Y Y  \\  Y Y  \\/ __ \\|  | \\/\\___  |                  \n"
  // "/_______  /____/|__|_|  /__|_|  (____  /__|   / ____|                  \n"
  // "        \\/            \\/      \\/     \\/       \\/ " << std::endl;
  m_output_stream << "Polarizable Embedding Summary:";
  m_output_stream << std::setprecision(12) << std::endl << std::endl;
  m_output_stream << off << "Electrostatics:" << std::endl;
  m_output_stream << off2 << "Electronic:" << std::string(w - 11, ' ')
                  << m_pe_energy.get("Electrostatic/Electronic") << std::endl;
  m_output_stream << off2 << "Nuclear:" << std::string(w - 8, ' ')
                  << m_pe_energy.get("Electrostatic/Nuclear") << std::endl;
  m_output_stream << off2 << "Multipole:" << std::string(w - 10, ' ')
                  << m_pe_energy.get("Electrostatic/Multipoles") << std::endl;
  m_output_stream << off2 << "Total:" << std::string(w - 6, ' ')
                  << m_pe_energy.get("Electrostatic") << std::endl;
  m_output_stream << std::endl;
  m_output_stream << off << "Polarization:" << std::endl;
  m_output_stream << off2 << "Electronic:" << std::string(w - 11, ' ')
                  << m_pe_energy.get("Polarization/Electronic") << std::endl;
  m_output_stream << off2 << "Nuclear:" << std::string(w - 8, ' ')
                  << m_pe_energy.get("Polarization/Nuclear") << std::endl;
  m_output_stream << off2 << "Multipole:" << std::string(w - 10, ' ')
                  << m_pe_energy.get("Polarization/Multipoles") << std::endl;
  m_output_stream << off2 << "Total:" << std::string(w - 6, ' ')
                  << m_pe_energy.get("Polarization") << std::endl;
  m_output_stream << std::endl;
  m_output_stream << off << "Total Energy:" << std::string(w - 10, ' ')
                  << m_pe_energy.get_total_energy() << std::endl;
  m_output_stream << std::string(2 * w + 10, '-') << std::endl << std::endl;
}

} // namespace libcppe
