#ifndef CPPE_UTILS_POT_MANIPULATION
#define CPPE_UTILS_POT_MANIPULATION

#include "../core/multipole.hh"
#include "../core/pe_options.hh"
#include "../core/molecule.hh"

namespace libcppe {

class PotManipulator {
private:
  std::vector<Potential> m_potentials;
  Molecule m_mol;

public:
  PotManipulator (std::vector<Potential> potentials, Molecule mol) :
    m_potentials(potentials), m_mol(mol)
    {
    };
  ~PotManipulator () {};
  std::vector<Potential> manipulate(PeOptions& pe_options);
};


} /* libcppe */



#endif // CPPE_UTILS_POT_MANIPULATION
