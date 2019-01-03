#ifndef CPPE_UTILS_POT_MANIPULATION
#define CPPE_UTILS_POT_MANIPULATION

#include <iostream>

#include "../core/molecule.hh"
#include "../core/multipole.hh"
#include "../core/pe_options.hh"

namespace libcppe {

class PotManipulator {
 private:
  std::vector<Potential> m_potentials;
  Molecule m_mol;
  std::ostream &m_output_stream;

 public:
  PotManipulator(std::vector<Potential> potentials, Molecule mol,
                 std::ostream &output_stream = std::cout)
      : m_potentials(potentials), m_mol(mol), m_output_stream(output_stream){};
  ~PotManipulator(){};
  std::vector<Potential> manipulate(PeOptions &pe_options);
};

}  // namespace libcppe

#endif  // CPPE_UTILS_POT_MANIPULATION
