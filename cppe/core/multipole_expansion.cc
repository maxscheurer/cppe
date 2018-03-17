#include "multipole_expansion.hh"

namespace libcppe {
  
std::vector<int> mul_v{1, 3, 6, 10, 15, 21};

// calculates the multipole-nuclei interaction energy through the given order
double MultipoleExpansion::calculate_interaction_energy() {
  double total_energy = 0.0;
  arma::Cube<int> Tk_coeffs = Tk_coefficients(5);
  for (auto& potential : m_potentials) {
    for (auto& multipole : potential.get_multipoles()) {
      std::vector<double> pref(mul_v[multipole.m_k]);
      prefactors_nuclei(multipole.m_k, pref);
      arma::vec pref_v(pref.data(), pref.size());
      arma::vec mul_v = multipole.get_values_vec();
      arma::vec site_position = multipole.get_site_position();
      for (auto& atom : m_mol) {
        arma::vec core_position = atom.get_pos();
        arma::vec diff = core_position-site_position;
        arma::vec Tsm = Tk_tensor(multipole.m_k, diff, Tk_coeffs);
        total_energy +=  arma::dot(pref_v, mul_v % Tsm) * atom.charge;
      }
    }
  }
  return total_energy;
}

} /* libcppe */