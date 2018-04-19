#include <iomanip>

#include "cppe_state.hh"

#include "multipole_expansion.hh"
#include "electric_fields.hh"

namespace libcppe {

CppeState::CppeState() {
  m_pe_energy = PeEnergy{};
}

void CppeState::set_potentials(std::vector<Potential> potentials) {
  m_potentials = potentials;
  m_polarizable_sites = 0;
  for (auto& site : m_potentials) {
    if (site.is_polarizable()) {
      m_polarizable_sites++;
    }
  }
  m_induced_moments = arma::vec(m_polarizable_sites*3, arma::fill::zeros);
}

void CppeState::update_energies(arma::mat P) {
  // Electrostatic energy (electrons-multipoles)
  double ees_elec = arma::dot(P, m_es_operator);
  m_pe_energy.set("Electrostatic/Electronic", ees_elec);

}

void CppeState::calculate_static_energies_and_fields() {
  // Electrostatic energy (nuclei-multipoles)
  MultipoleExpansion mexp(m_mol, m_potentials);
  double nuc_mul_energy = mexp.calculate_interaction_energy();
  m_pe_energy.set("Electrostatic/Nuclear", nuc_mul_energy);


  // Calculate static fields

  // Nuclear fields
  NuclearFields nfields(m_mol, m_potentials);
  m_nuc_fields = arma::vec(m_polarizable_sites*3, arma::fill::zeros);
  nfields.compute(m_nuc_fields, false);
  // Multipole fields
  MultipoleFields mul_fields(m_potentials);
  m_multipole_fields = arma::vec(m_polarizable_sites*3, arma::fill::zeros);
  mul_fields.compute(m_multipole_fields, false);


}

void CppeState::update_induced_moments(arma::vec elec_fields, int iteration, bool elec_only) {
  arma::vec tmp_total_fields(m_polarizable_sites*3, arma::fill::zeros);
  if (elec_only) {
    tmp_total_fields = elec_fields;
  } else {
    tmp_total_fields = elec_fields + m_nuc_fields + m_multipole_fields;

    // elec_fields.save("elec_fields.txt", arma::raw_ascii);
    // m_nuc_fields.save("nuc_fields.txt", arma::raw_ascii);
    // m_multipole_fields.save("multipole_fields.txt", arma::raw_ascii);
  }

  bool make_guess = true;
  if (iteration > 0) {
    make_guess = false;
  }
  InducedMoments ind(m_potentials);
  ind.compute(tmp_total_fields, m_induced_moments, make_guess);

  if (elec_only) {
    double epol_elec = -0.5*arma::dot(m_induced_moments, elec_fields);
    std::cout << std::setprecision(15) << "epol_elec = " << epol_elec << std::endl;
    m_pe_energy.set("Polarization/Electronic", epol_elec);
  } else {
    double epol_elec = -0.5*arma::dot(m_induced_moments, elec_fields);
    double epol_nuclear = -0.5*arma::dot(m_induced_moments, m_nuc_fields);
    double epol_multipoles = -0.5*arma::dot(m_induced_moments, m_multipole_fields);

    std::cout << "epol_elec = " << epol_elec << std::endl;
    std::cout << "epol_nuclear = " << epol_nuclear << std::endl;
    std::cout << "epol_multipoles = " << epol_multipoles << std::endl;

    m_pe_energy.set("Polarization/Electronic", epol_elec);
    m_pe_energy.set("Polarization/Nuclear", epol_nuclear);
    m_pe_energy.set("Polarization/Multipoles", epol_multipoles);
  }

  // m_induced_moments.save("induced_moments.txt", arma::raw_ascii);
}

void CppeState::print_summary() {
  size_t w = 30;
  std::string off(3, ' ');
  std::string off2(6, ' ');
  // std::cout << off << "Corrected Excitation Energy (" << energy_name << "): " << std::string(w - 52, ' ');
  // std::cout << std::setw(11) << (state.energy + ene) * conversion::au2ev << " eV" << std::endl;
  std::cout << std::string(2*w+10, '-') << std::endl;
  // std::cout << "__________      .__               .__              ___.   .__          \n"
  // "\\______   \\____ |  | _____ _______|__|____________ \\_ |__ |  |   ____  \n"
  // " |     ___/  _ \\|  | \\__  \\\\_  __ \\  \\___   /\\__  \\ | __ \\|  | _/ __ \\ \n"
  // " |    |  (  <_> )  |__/ __ \\|  | \\/  |/    /  / __ \\| \\_\\ \\  |_\\  ___/ \n"
  // " |____|   \\____/|____(____  /__|  |__/_____ \\(____  /___  /____/\\___  >\n"
  // "___________      ___.     \\/       .___  .___.__  \\/    \\/          \\/ \n"
  // "\\_   _____/ _____\\_ |__   ____   __| _/__| _/|__| ____    ____         \n"
  // " |    __)_ /     \\| __ \\_/ __ \\ / __ |/ __ | |  |/    \\  / ___\\        \n"
  // " |        \\  Y Y  \\ \\_\\ \\  ___// /_/ / /_/ | |  |   |  \\/ /_/  >       \n"
  // "/_______  /__|_|  /___  /\\___  >____ \\____ | |__|___|  /\\___  /        \n"
  // "  _________     \\/    \\/     \\/     \\/    \\/         \\//_____/         \n"
  // " /   _____/__ __  _____   _____ _____ _______ ___.__.                  \n"
  // " \\_____  \\|  |  \\/     \\ /     \\\\__  \\\\_  __ <   |  |                  \n"
  // " /        \\  |  /  Y Y  \\  Y Y  \\/ __ \\|  | \\/\\___  |                  \n"
  // "/_______  /____/|__|_|  /__|_|  (____  /__|   / ____|                  \n"
  // "        \\/            \\/      \\/     \\/       \\/                       " << std::endl;
  std::cout << "Polarizable Embedding Summary:";
  std::cout << std::endl << std::endl;
  std::cout << off << "Electrostatics:" << std::endl;
  std::cout << off2 << "Electronic:" << std::string(w-11, ' ') << m_pe_energy.get("Electrostatic/Electronic") << std::endl;
  std::cout << off2 << "Nuclear:" << std::string(w-8, ' ') << m_pe_energy.get("Electrostatic/Nuclear") << std::endl;
  std::cout << off2 << "Multipole:" << std::string(w-10, ' ') << m_pe_energy.get("Electrostatic/Multipoles") << std::endl;
  std::cout << off2 << "Total:" << std::string(w-6, ' ') << m_pe_energy.get("Electrostatic") << std::endl;
  std::cout << std::endl;
  std::cout << off << "Polarization:" << std::endl;
  std::cout << off2 << "Electronic:" << std::string(w-11, ' ') << m_pe_energy.get("Polarization/Electronic") << std::endl;
  std::cout << off2 << "Nuclear:" << std::string(w-8, ' ') << m_pe_energy.get("Polarization/Nuclear") << std::endl;
  std::cout << off2 << "Multipole:" << std::string(w-10, ' ') << m_pe_energy.get("Polarization/Multipoles") << std::endl;
  std::cout << off2 << "Total:" << std::string(w-6, ' ') << m_pe_energy.get("Polarization") << std::endl;
  std::cout << std::endl;
  std::cout << off << "Total Energy:" << std::string(w-10, ' ') << m_pe_energy.get_total_energy() << std::endl;
  std::cout << std::string(2*w+10, '-') << std::endl << std::endl;
}

} /* libcppe */
