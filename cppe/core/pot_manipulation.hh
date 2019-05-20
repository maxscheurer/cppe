#pragma once

#include <iostream>

#include "molecule.hh"
#include "potential.hh"
#include "pe_options.hh"

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
  std::vector<Potential> manipulate(const PeOptions &pe_options);
  std::vector<Potential> manipulate_border(const PeOptions &pe_options);
};

}  // namespace libcppe
