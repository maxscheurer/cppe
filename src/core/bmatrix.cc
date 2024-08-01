#include "bmatrix.hh"
#include "math.hh"
#include "tensors.hh"

#include "fmm/tree.hh"

namespace libcppe {

Eigen::VectorXd BMatrix::apply(const Eigen::VectorXd& induced_moments) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(induced_moments.size());

  int n_crit   = m_options.tree_ncrit;
  int order    = m_options.tree_expansion_order;
  double theta = m_options.theta;
  std::vector<double> S(3 * m_n_polsites);
  for (decltype(m_n_polsites) i = 0; i < m_n_polsites; ++i) {
    int l             = i * 3;
    Eigen::Vector3d s = induced_moments.segment<3>(l);
    S[i * 3 + 0]      = s(0);
    S[i * 3 + 1]      = s(1);
    S[i * 3 + 2]      = s(2);
  }

  double damping = 0.0;
  if (m_options.damp_induced) {
    damping = m_options.damping_factor_induced;
  }
  std::shared_ptr<Tree<1, 3>> tree =
        build_shared_tree<1, 3>(m_polsites, S.data(), n_crit, order, theta, damping);
  std::vector<double> induced_fields_v(3 * m_n_polsites);
  auto scheme = m_options.summation_induced_fields;
  if (scheme == "fmm") {
    tree->compute_field_fmm(induced_fields_v.data());
  } else {
    tree->compute_field_exact(induced_fields_v.data());
  }
  Eigen::VectorXd induced_fields =
        Eigen::Map<Eigen::VectorXd>(induced_fields_v.data(), induced_fields_v.size());
  return induced_fields + apply_diagonal(induced_moments);
}

Eigen::VectorXd BMatrix::apply_diagonal_inverse(const Eigen::VectorXd& in) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(3 * m_n_polsites);
#pragma omp parallel for
  for (decltype(m_n_polsites) i = 0; i < m_n_polsites; ++i) {
    int l                   = i * 3;
    Polarizability& alpha_i = m_polsites[i].get_polarizability();
    ret.segment<3>(l)       = alpha_i.get_matrix() * in.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::apply_diagonal(const Eigen::VectorXd& in) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(3 * m_n_polsites);

#pragma omp parallel for
  for (decltype(m_n_polsites) i = 0; i < m_n_polsites; ++i) {
    int l             = i * 3;
    ret.segment<3>(l) = m_alpha_inverse[i] * in.segment<3>(l);
  }
  return ret;
}

Eigen::MatrixXd BMatrix::to_dense_matrix() {
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_n_polsites * 3, m_n_polsites * 3);
#pragma omp parallel for
  for (decltype(m_n_polsites) i = 0; i < m_n_polsites; ++i) {
    int l           = i * 3;
    Potential& pot1 = m_polsites[i];
    for (decltype(m_n_polsites) j = 0; j < m_n_polsites; ++j) {
      int m           = j * 3;
      Potential& pot2 = m_polsites[j];
      if (pot1.excludes_site(pot2.index) || i == j) continue;
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
  return B;
}

Eigen::MatrixXd BMatrix::direct_inverse() {
  Eigen::MatrixXd B = to_dense_matrix();
  return B.inverse();
}

}  // namespace libcppe
