#ifndef LIBCPPE_CPPE_CORE_ELECTRIC_FIELDS_H
#define LIBCPPE_CPPE_CORE_ELECTRIC_FIELDS_H

#include "molecule.hh"
#include "multipole_expansion.hh"

// TODO: create Field class which others inherit

namespace libcppe {
  

class ElectronFiels {
private:
  /* data */

public:
  ElectronFiels ();
  ~ElectronFiels ();
};

class NuclearFields {
private:
  std::vector<Potential> m_potentials;
  Molecule m_mol;

public:
  NuclearFields(Molecule mol, std::vector<Potential> potentials) : 
    m_mol(mol), m_potentials(potentials) {};
  ~NuclearFields() {};
  // nuc_fields has length of 3*n_polarizable_sites
  void compute(arma::vec& nuc_fields, bool damp_core);
};


class MultipoleFields {
private:
  std::vector<Potential> m_potentials;

public:
  MultipoleFields(std::vector<Potential> potentials) : 
    m_potentials(potentials) {};
  ~MultipoleFields() {};
  // nuc_fields has length of 3*n_polarizable_sites
  void compute(arma::vec& mult_fields, bool damp);
};


class InducedMoments {
private:
  std::vector<Potential> m_potentials;

public:
  InducedMoments(std::vector<Potential> potentials) : m_potentials(potentials) {};
  ~InducedMoments() {};
  void compute(arma::vec& total_fields, arma::vec& induced_moments, bool make_guess); 
};



} // namespace libcppe

#endif // LIBCPPE_CPPE_CORE_ELECTRIC_FIELDS_H