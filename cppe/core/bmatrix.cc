#include "bmatrix.hh"
#include "math.hh"
#include "tensors.hh"

namespace libcppe {

Eigen::VectorXd BMatrix::compute_apply(Eigen::VectorXd induced_moments) {
  return compute_apply_slice(induced_moments, 0, m_n_polsites);
}

Eigen::VectorXd BMatrix::compute_apply_slice(Eigen::VectorXd induced_moments, int start,
                                             int stop) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(induced_moments.size());
  if (start < 0 || stop > m_n_polsites) {
    throw std::runtime_error("Invalid range in compute_apply_slice.");
  }

#pragma omp parallel for
  for (int i = start; i < stop; ++i) {
    int l           = i * 3;
    Potential& pot1 = m_polsites[i];
    for (auto j : m_polmask[i]) {
      int m                = j * 3;
      Potential& pot2      = m_polsites[j];
      Eigen::Vector3d diff = pot2.get_site_position() - pot1.get_site_position();
      Eigen::VectorXd T2;
      if (m_options.damp_induced) {
        Polarizability& alpha_i = pot1.get_polarizability();
        Polarizability& alpha_j = pot2.get_polarizability();
        double v                = m_options.damping_factor_induced /
                   std::pow(alpha_i.get_isotropic_value() * alpha_j.get_isotropic_value(),
                            1.0 / 6.0);
        T2 = tensors::T2_damp_thole(diff, v);
      } else {
        T2 = tensors::T2(diff);
      }
      Eigen::Matrix3d T2m = triangle_to_mat(T2);
      ret.segment<3>(l) -= T2m * induced_moments.segment<3>(m);
    }
    ret.segment<3>(l) += m_alpha_inverse[i] * induced_moments.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::compute_apply_diagonal(Eigen::VectorXd in) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(in.size());

#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    int l                   = i * 3;
    Polarizability& alpha_i = m_polsites[i].get_polarizability();
    ret.segment<3>(l)       = alpha_i.get_matrix() * in.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::compute_gauss_seidel_update(
      Eigen::VectorXd induced_moments, const Eigen::VectorXd& total_fields) {
#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    Eigen::Vector3d Ftmp = Eigen::Vector3d::Zero();
    int l                = i * 3;
    Potential& pot1      = m_polsites[i];
    for (int j = 0; j < m_n_polsites; ++j) {
      int m           = 3 * j;
      Potential& pot2 = m_polsites[j];
      if (pot1.excludes_site(pot2.index) || i == j) {
        continue;
      }
      Eigen::Vector3d diff = pot2.get_site_position() - pot1.get_site_position();
      Eigen::VectorXd T2;
      if (m_options.damp_induced) {
        Polarizability& alpha_i = pot1.get_polarizability();
        Polarizability& alpha_j = pot2.get_polarizability();
        double v                = m_options.damping_factor_induced /
                   std::pow(alpha_i.get_isotropic_value() * alpha_j.get_isotropic_value(),
                            1.0 / 6.0);
        T2 = tensors::T2_damp_thole(diff, v);
      } else {
        T2 = tensors::T2(diff);
      }
      Eigen::Matrix3d T2m = triangle_to_mat(T2);
      Ftmp += T2m * induced_moments.segment<3>(m);
    }
    Ftmp += total_fields.segment<3>(l);
    induced_moments.segment<3>(l) = pot1.get_polarizability().get_matrix() * Ftmp;
  }
  return induced_moments;
}

Eigen::MatrixXd BMatrix::direct_inverse() {
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_n_polsites * 3, m_n_polsites * 3);
  for (int i = 0; i < m_n_polsites; ++i) {
    int l           = i * 3;
    Potential& pot1 = m_polsites[i];
    for (auto j : m_polmask[i]) {
      int m                = j * 3;
      Potential& pot2      = m_polsites[j];
      Eigen::Vector3d diff = pot2.get_site_position() - pot1.get_site_position();
      Eigen::VectorXd T2;
      if (m_options.damp_induced) {
        Polarizability& alpha_i = pot1.get_polarizability();
        Polarizability& alpha_j = pot2.get_polarizability();
        double v                = m_options.damping_factor_induced /
                   std::pow(alpha_i.get_isotropic_value() * alpha_j.get_isotropic_value(),
                            1.0 / 6.0);
        T2 = tensors::T2_damp_thole(diff, v);
      } else {
        T2 = tensors::T2(diff);
      }
      Eigen::Matrix3d T2m = triangle_to_mat(T2);
      B.block<3, 3>(l, m) = -T2m;
    }
    B.block<3, 3>(l, l) = m_alpha_inverse[i];
  }
  return B.inverse();
}

}  // namespace libcppe
