#pragma once

#include <Eigen/Dense>

#include "pe_options.hh"
#include "potential.hh"

namespace libcppe {

class BMatrix {
 private:
  std::vector<Potential> m_polsites;  //!< vector with all potentials of polarizable sites
  size_t m_n_polsites;                //!< number of polarizable sites
  PeOptions m_options;
  std::vector<Eigen::Matrix3d>
        m_alpha_inverse;            //!< List with inverse polarizability tensors
  std::vector<double> m_positions;  //! list of all polarizabile site positions
 public:
  BMatrix(std::vector<Potential> polsites, const PeOptions& options)
        : m_polsites(polsites), m_options(options) {
    m_n_polsites = polsites.size();
    m_alpha_inverse.clear();
    m_positions = std::vector<double>(3 * m_n_polsites);
    std::transform(m_polsites.begin(), m_polsites.end(),
                   std::back_inserter(m_alpha_inverse),
                   [](Potential& p) -> Eigen::Matrix3d {
                     return p.get_polarizability().get_matrix().inverse();
                   });
#pragma omp parallel for
    for (int i = 0; i < m_n_polsites; ++i) {
      Potential& pot1        = m_polsites[i];
      m_positions[i * 3 + 0] = pot1.m_x;
      m_positions[i * 3 + 1] = pot1.m_y;
      m_positions[i * 3 + 2] = pot1.m_z;
    }
  }

  Eigen::VectorXd apply(const Eigen::VectorXd& induced_moments);
  Eigen::VectorXd apply_diagonal_inverse(const Eigen::VectorXd& in);
  Eigen::VectorXd apply_diagonal(const Eigen::VectorXd& in);
  Eigen::MatrixXd to_dense_matrix();
  Eigen::MatrixXd direct_inverse();
};

}  // namespace libcppe
