#pragma once

#include "molecule.hh"
#include "multipole_expansion.hh"
#include "pe_options.hh"

namespace libcppe {

Eigen::VectorXd multipole_derivative(int k, int l, const Eigen::Vector3d& Rji,
                                     Eigen::VectorXd Mkj, double damping_factor = 0.0);

std::vector<Potential> get_polarizable_sites(std::vector<Potential> potentials);

class ElectricFields {
 protected:
  std::vector<Potential> m_potentials;  //!< vector with all site potentials
  std::vector<Potential> m_polsites;  //!< vector with all potentials of polarizable sites
  size_t m_n_polsites;                //!< number of polarizable sites

 public:
  ElectricFields(std::vector<Potential> potentials) : m_potentials(potentials) {
    m_polsites   = get_polarizable_sites(m_potentials);
    m_n_polsites = m_polsites.size();
  };
};

class NuclearFields : public ElectricFields {
 private:
  Molecule m_mol;  //!< core region molecule

 public:
  NuclearFields(Molecule mol, std::vector<Potential> potentials)
        : ElectricFields(potentials), m_mol(mol){};
  Eigen::VectorXd compute();
};

class MultipoleFields : public ElectricFields {
 public:
  MultipoleFields(std::vector<Potential> potentials, const PeOptions& options)
        : ElectricFields(potentials), m_options(options) {
    m_exclusions.clear();
    int n_sites  = m_potentials.size();
    m_positions  = std::vector<double>(3 * n_sites);
    m_exclusions = std::vector<std::vector<int>>(n_sites);
#pragma omp parallel for
    for (int i = 0; i < n_sites; ++i) {
      Potential& pot1        = m_potentials[i];
      m_positions[i * 3 + 0] = pot1.m_x;
      m_positions[i * 3 + 1] = pot1.m_y;
      m_positions[i * 3 + 2] = pot1.m_z;
      m_exclusions[i]        = pot1.get_exclusions();
    }
  };
  Eigen::VectorXd compute();
  Eigen::VectorXd compute_legacy();
  Eigen::VectorXd compute_tree();

 private:
  PeOptions m_options;
  std::vector<std::vector<int>>
        m_exclusions;               //!< List for each site with sites that are excluded
  std::vector<double> m_positions;  //! list of all polarizabile site positions
};

class InducedMoments {
 private:
  std::vector<Potential> m_potentials;  //!< vector with all site potentials
  std::vector<Potential> m_polsites;  //!< vector with all potentials of polarizable sites
  size_t m_n_polsites;                //!< number of polarizable sites
  PeOptions m_options;
  std::function<void(std::string)> m_printer = [](std::string str) {
    std::cout << str << std::endl;
  };

 public:
  InducedMoments(std::vector<Potential> potentials, const PeOptions& options)
        : m_potentials(potentials), m_options(options) {
    m_polsites   = get_polarizable_sites(m_potentials);
    m_n_polsites = m_polsites.size();
  };
  void set_print_callback(std::function<void(std::string)> printer) {
    m_printer = printer;
  }
  void compute(const Eigen::VectorXd& total_fields, Eigen::VectorXd& induced_moments,
               bool make_guess);
  Eigen::VectorXd compute_cg(const Eigen::VectorXd& rhs, Eigen::VectorXd guess,
                             bool make_guess);
  /**
      overloads the compute method for induced moments and returns
     a copy of the induced moments vector
  */
  Eigen::VectorXd compute(Eigen::VectorXd& total_fields, bool make_guess) {
    Eigen::VectorXd induced_moments = Eigen::VectorXd::Zero(total_fields.size());
    compute(total_fields, induced_moments, make_guess);
    return induced_moments;
  }
  Eigen::VectorXd compute_cg(Eigen::VectorXd& total_fields, bool make_guess) {
    Eigen::VectorXd induced_moments = Eigen::VectorXd::Zero(total_fields.size());
    return compute_cg(total_fields, induced_moments, make_guess);
  }
};

}  // namespace libcppe
