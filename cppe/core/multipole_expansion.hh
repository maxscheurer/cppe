#pragma once

#include "molecule.hh"
#include "potential.hh"

namespace libcppe {

class MultipoleExpansion {
 private:
  Molecule m_mol;                       //!< core region molecule
  std::vector<Potential> m_potentials;  //!< vector with all site potentials

 public:
  MultipoleExpansion(Molecule core, std::vector<Potential> potentials)
        : m_mol(core),
          m_potentials(potentials){

          };

  double interaction_energy();
  Eigen::MatrixXd nuclear_gradient();
};

}  // namespace libcppe
