#include "bmatrix.hh"
#include "math.hh"
#include "tensors.hh"

#include "fmm/tree.hh"

namespace libcppe {

Eigen::VectorXd BMatrix::apply_direct(const Eigen::VectorXd& induced_moments) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(induced_moments.size());
#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    int l           = i * 3;
    Potential& pot1 = m_polsites[i];
    double* retx    = &ret[l + 0];
    double* rety    = &ret[l + 1];
    double* retz    = &ret[l + 2];
    for (int j = 0; j < m_n_polsites; ++j) {
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
      double fx = induced_moments[m + 0];
      double fy = induced_moments[m + 1];
      double fz = induced_moments[m + 2];
      // inline matrix-vector product, faster than unfolding and multiplying
      *retx -= T2[0] * fx + T2[1] * fy + T2[2] * fz;
      *rety -= T2[1] * fx + T2[3] * fy + T2[4] * fz;
      *retz -= T2[2] * fx + T2[4] * fy + T2[5] * fz;
    }
    ret.segment<3>(l) += m_alpha_inverse[i] * induced_moments.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::apply_fast_summation(const Eigen::VectorXd& induced_moments,
                                              std::string scheme) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(induced_moments.size());

  int n_crit   = m_options.tree_ncrit;
  int order    = m_options.tree_expansion_order;
  double theta = m_options.theta;
  std::vector<double> S(3 * m_n_polsites);
  for (int i = 0; i < m_n_polsites; ++i) {
    int l             = i * 3;
    Eigen::Vector3d s = induced_moments.segment<3>(l);
    S[i * 3 + 0]      = s(0);
    S[i * 3 + 1]      = s(1);
    S[i * 3 + 2]      = s(2);
  }
  std::shared_ptr<Tree<1, 3>> tree = build_shared_tree<1, 3>(
        m_positions.data(), S.data(), m_n_polsites, n_crit, order, theta, m_exclusions);
  std::vector<double> induced_fields_v(3 * m_n_polsites);
  if (scheme == "fmm") {
    tree->compute_field_fmm(induced_fields_v.data());
  } else {
    throw std::runtime_error("No such summation scheme " + scheme);
  }
  Eigen::VectorXd induced_fields =
        Eigen::Map<Eigen::VectorXd>(induced_fields_v.data(), induced_fields_v.size());
  return induced_fields + apply_diagonal(induced_moments);
}

Eigen::VectorXd BMatrix::apply(const Eigen::VectorXd& induced_moments) {
  auto scheme = m_options.summation_induced_fields;
  if (scheme == "direct") {
    return apply_direct(induced_moments);
  } else if (scheme == "fmm" or scheme == "bh") {
    return apply_fast_summation(induced_moments, scheme);
  } else {
    throw std::invalid_argument("Invalid summation scheme for induced fields provided.");
  }
}

Eigen::VectorXd BMatrix::apply_diagonal_inverse(const Eigen::VectorXd& in) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(3 * m_n_polsites);

#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    int l                   = i * 3;
    Polarizability& alpha_i = m_polsites[i].get_polarizability();
    ret.segment<3>(l)       = alpha_i.get_matrix() * in.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::apply_diagonal(const Eigen::VectorXd& in) {
  Eigen::VectorXd ret = Eigen::VectorXd::Zero(3 * m_n_polsites);

#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    int l             = i * 3;
    ret.segment<3>(l) = m_alpha_inverse[i] * in.segment<3>(l);
  }
  return ret;
}

Eigen::VectorXd BMatrix::gauss_seidel_update(Eigen::VectorXd induced_moments,
                                             const Eigen::VectorXd& total_fields) {
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
      double fx = induced_moments[m + 0];
      double fy = induced_moments[m + 1];
      double fz = induced_moments[m + 2];
      // inline matrix-vector product, faster than unfolding and multiplying
      Ftmp[0] += T2[0] * fx + T2[1] * fy + T2[2] * fz;
      Ftmp[1] += T2[1] * fx + T2[3] * fy + T2[4] * fz;
      Ftmp[2] += T2[2] * fx + T2[4] * fy + T2[5] * fz;
    }
    Ftmp += total_fields.segment<3>(l);
    induced_moments.segment<3>(l) = pot1.get_polarizability().get_matrix() * Ftmp;
  }
  return induced_moments;
}

Eigen::MatrixXd BMatrix::to_dense_matrix() {
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m_n_polsites * 3, m_n_polsites * 3);
#pragma omp parallel for
  for (int i = 0; i < m_n_polsites; ++i) {
    int l           = i * 3;
    Potential& pot1 = m_polsites[i];
    for (int j = 0; j < m_n_polsites; ++j) {
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
