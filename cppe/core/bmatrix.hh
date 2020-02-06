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
        m_polmask;  //!< List for each polarizable sites with sites that are not excluded
 public:
  BMatrix(std::vector<Potential> polsites, PeOptions options)
        : m_polsites(polsites), m_options(options) {
    m_n_polsites = polsites.size();
    m_alpha_inverse.clear();
    m_polmask.clear();
    std::transform(m_polsites.begin(), m_polsites.end(),
                   std::back_inserter(m_alpha_inverse),
                   [](Potential& p) -> Eigen::Matrix3d {
                     return p.get_polarizability().get_matrix().inverse();
                   });

    for (int i = 0; i < m_n_polsites; ++i) {
      Potential& pot1 = m_polsites[i];
      std::vector<int> pot_pols;
      for (int j = 0; j < m_n_polsites; ++j) {
        Potential& pot2 = m_polsites[j];
        if (pot1.excludes_site(pot2.index) || i == j) {
          continue;
        } else {
          pot_pols.push_back(j);
        }
      }
      m_polmask.push_back(pot_pols);
    }
  }

  Eigen::VectorXd compute_apply(Eigen::VectorXd induced_moments);
  Eigen::VectorXd compute_apply_slice(Eigen::VectorXd induced_moments, int start,
                                      int stop);
  Eigen::VectorXd compute_apply_diagonal(Eigen::VectorXd in);
  Eigen::VectorXd compute_gauss_seidel_update(Eigen::VectorXd induced_moments,
                                              const Eigen::VectorXd& total_fields);
  Eigen::MatrixXd direct_inverse();
};

}  // namespace libcppe