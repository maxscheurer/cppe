#pragma once

#include <functional>

#include "molecule.hh"
#include "pe_options.hh"
#include "potential.hh"

namespace libcppe {

class PotManipulator {
 private:
  std::vector<Potential> m_potentials;
  Molecule m_mol;
  std::function<void(std::string)> m_printer;

 public:
  PotManipulator(std::vector<Potential> potentials, Molecule mol)
        : m_potentials(potentials), m_mol(mol){};
  ~PotManipulator(){};
  void set_print_callback(std::function<void(std::string)> printer) {
    m_printer = printer;
  }
  std::vector<Potential> manipulate(const PeOptions& pe_options);
  std::vector<Potential> manipulate_border(const PeOptions& pe_options);
};

}  // namespace libcppe
