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
  ~ElectricFields(){};
  virtual Eigen::VectorXd compute() = 0;
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
  MultipoleFields(std::vector<Potential> potentials, PeOptions options)
        : ElectricFields(potentials), m_options(options){};
  Eigen::VectorXd compute();

 private:
  PeOptions m_options;
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
  InducedMoments(std::vector<Potential> potentials, PeOptions options)
        : m_potentials(potentials), m_options(options) {
    m_polsites   = get_polarizable_sites(m_potentials);
    m_n_polsites = m_polsites.size();
  };
  ~InducedMoments(){};
  void set_print_callback(std::function<void(std::string)> printer) {
    m_printer = printer;
  }
  void compute(const Eigen::VectorXd& total_fields, Eigen::VectorXd& induced_moments,
               bool make_guess);
  Eigen::VectorXd compute_cg(const Eigen::VectorXd& rhs);
  /**
      overloads the compute method for induced moments and returns
     a copy of the induced moments vector
  */
  Eigen::VectorXd compute(Eigen::VectorXd& total_fields, bool make_guess) {
    Eigen::VectorXd induced_moments = Eigen::VectorXd::Zero(total_fields.size());
    compute(total_fields, induced_moments, make_guess);
    return induced_moments;
  }
};

}  // namespace libcppe
