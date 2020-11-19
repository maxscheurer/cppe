#include <Eigen/Dense>

#include "multipole_expansion.hh"
#include "electric_fields.hh"
#include "tensors.hh"

namespace libcppe {

// calculates the multipole-nuclei interaction energy through the given order
double MultipoleExpansion::interaction_energy() {
  double total_energy = 0.0;
  int npots           = m_potentials.size();

#pragma omp parallel for reduction(+ : total_energy)
  for (size_t i = 0; i < npots; i++) {
    Potential& potential          = m_potentials[i];
    Eigen::Vector3d site_position = potential.get_site_position();
    for (auto& multipole : potential.get_multipoles()) {
      // TODO: refactor in math.cc
      std::vector<double> pref = prefactors_nuclei(multipole.m_k);
      Eigen::VectorXd pref_v =
            Eigen::Map<Eigen::VectorXd>(std::move(pref.data()), pref.size());
      Eigen::VectorXd mul_v = multipole.get_values_vec();
      for (auto& atom : m_mol) {
        Eigen::Vector3d core_position = atom.get_position();
        Eigen::Vector3d diff          = core_position - site_position;
        Eigen::VectorXd Tsm           = tensors::T[multipole.m_k](diff);
        total_energy += pref_v.dot(mul_v.cwiseProduct(Tsm)) * atom.charge;
      }
    }
  }
  return total_energy;
}

Eigen::MatrixXd MultipoleExpansion::nuclear_gradient() {
  int natoms = m_mol.size();
  int npots           = m_potentials.size();
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(natoms, 3);

  #pragma omp parallel for
  for (size_t i = 0; i < npots; i++) {
    Potential& potential          = m_potentials[i];
    Eigen::Vector3d site_position = potential.get_site_position();
    for (auto& multipole : potential.get_multipoles()) {
      // TODO: refactor in math.cc
      std::vector<double> pref = prefactors_nuclei(multipole.m_k);
      Eigen::VectorXd pref_v =
            Eigen::Map<Eigen::VectorXd>(std::move(pref.data()), pref.size());
      Eigen::VectorXd mul_v = multipole.get_values_vec();
      for (int ai = 0; ai < natoms; ++ai) {
        auto& atom = m_mol[ai];
        Eigen::Vector3d core_position = atom.get_position();
        Eigen::Vector3d diff          = core_position - site_position;
        Eigen::Vector3d tmp_grad = multipole_derivative(multipole.m_k, 1, diff, multipole.get_values_vec());
        #pragma omp critical
        {
          grad.row(ai) += -tmp_grad * atom.charge;
        }
      }
    }
  }
  return grad;
}

}  // namespace libcppe
