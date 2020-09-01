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
        m_alpha_inverse;  //!< List with inverse polarizability tensors
  std::vector<std::vector<int>>
        m_exclusions;  //!< List for each polarizable sites with sites that are excluded
  std::vector<double> m_positions;  //! list of all polarizabile site positions
 public:
  BMatrix(std::vector<Potential> polsites, const PeOptions& options)
        : m_polsites(polsites), m_options(options) {
    m_n_polsites = polsites.size();
    m_alpha_inverse.clear();
    m_exclusions.clear();
    m_positions = std::vector<double>(3 * m_n_polsites);
    std::transform(m_polsites.begin(), m_polsites.end(),
                   std::back_inserter(m_alpha_inverse),
                   [](Potential& p) -> Eigen::Matrix3d {
                     return p.get_polarizability().get_matrix().inverse();
                   });

    m_exclusions = std::vector<std::vector<int>>(m_n_polsites);
#pragma omp parallel for
    for (int i = 0; i < m_n_polsites; ++i) {
      Potential& pot1        = m_polsites[i];
      m_positions[i * 3 + 0] = pot1.m_x;
      m_positions[i * 3 + 1] = pot1.m_y;
      m_positions[i * 3 + 2] = pot1.m_z;
      std::vector<int> pot_excludes;
      for (int j = 0; j < m_n_polsites; ++j) {
        Potential& pot2 = m_polsites[j];
        if (pot1.excludes_site(pot2.index)) {
          pot_excludes.push_back(j);
        }
      }
      m_exclusions[i] = pot_excludes;
    }
  }

  Eigen::VectorXd apply(const Eigen::VectorXd& induced_moments);
  Eigen::VectorXd apply_direct(const Eigen::VectorXd& induced_moments);
  Eigen::VectorXd apply_fast_summation(const Eigen::VectorXd& induced_moments,
                                       std::string scheme);
  Eigen::VectorXd apply_diagonal_inverse(const Eigen::VectorXd& in);
  Eigen::VectorXd apply_diagonal(const Eigen::VectorXd& in);
  Eigen::VectorXd gauss_seidel_update(Eigen::VectorXd induced_moments,
                                      const Eigen::VectorXd& total_fields);
  Eigen::MatrixXd to_dense_matrix();
  Eigen::MatrixXd direct_inverse();
  std::vector<std::vector<int>> get_exclusions() { return m_exclusions; }
};

}  // namespace libcppe
