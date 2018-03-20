#include <iostream>
#include <cassert>

#include "interface.h"

#include "libcppe.hh"

#include "utils/potfile_reader.hh"
#include "core/multipole_expansion.hh"
#include "core/electric_fields.hh"

namespace libcppe {

CPPE::CPPE() {
  // constructor for CPPE class
}

CPPE::~CPPE() {
  // destructor for CPPE class
  // needs to destroy pelib and gen1int
}

void CPPE::initialize_gen1int(int natoms, int nshells, const double* coords, const double *charges) {
  gen1int_api_initialize(natoms, nshells, coords, charges);
  m_gen1int_initialized = true;
}

void CPPE::gen1int_add_shell(int spher_gto, int idx_center, const double* coord_center,
  int ang_num, int num_prim, const double* exponents, int num_contr, const double* contr_coef) {
  gen1int_api_add_shell(spher_gto, idx_center, coord_center, ang_num, num_prim, exponents, num_contr,
  contr_coef);
}

void CPPE::gen1int_print_shells() {
  if (!m_gen1int_initialized) {
    throw std::runtime_error("Gen1int needs to be initialized before calling print_shells");
  }
  print_shells();
}

void CPPE::initialize_pelib(std::string potfile, int natoms, int nbas, const double* coords, const double* charges) {
  if (!m_gen1int_initialized) {
    throw std::runtime_error("Gen1int needs to be initialized before calling PElib");
  }
  m_natoms = natoms;
  m_nbas = nbas;
  m_nnbas = nbas*(nbas+1)/2;
  assert(potfile.size() < 80);
  pe_set_potfile(potfile.c_str(), static_cast<int>(potfile.size()));
  pe_interface_init(natoms, coords, charges);
  m_pe_initialized = true;
  
}

double CPPE::calculate_nulcei_multipole_interaction(Molecule& mol, std::vector<Potential>& potentials) {
  std::cout << "calculate_nulcei_multipole_interaction" << std::endl;
  MultipoleExpansion mexp(mol,potentials);
  double nuc_mul_energy = mexp.calculate_interaction_energy();
  return nuc_mul_energy;
}

arma::vec CPPE::calculate_nulcear_fields(Molecule& mol, std::vector<Potential>& potentials, size_t polarizable_sites) {
  std::cout << "calculate nuclear fields" << std::endl;
  NuclearFields nfields(mol, potentials);
  arma::vec result(polarizable_sites*3, arma::fill::zeros);
  nfields.compute(result, false);
  return result;
}

arma::vec CPPE::calculate_multipole_fields(std::vector<Potential>& potentials, size_t polarizable_sites) {
  MultipoleFields mul_fields(potentials);
  arma::vec result(polarizable_sites*3, arma::fill::zeros);
  mul_fields.compute(result, false);
  return result;
}

arma::vec CPPE::calculate_induced_moments(std::vector<Potential>& potentials, arma::vec& total_fields) {
  InducedMoments ind(potentials);
  return ind.compute(total_fields);
}

void CPPE::call_pe_energy(const double* densmat) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_energy(densmat, m_nbas, m_nnbas);
}

void CPPE::call_pe_pol_energy(const double* densmat, double *energy) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_pol_energy(densmat, m_nbas, m_nnbas, energy);
}

void CPPE::call_full_fock(const double* densmat, double* fockmat, double* energy) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_fock(densmat, m_nbas, m_nnbas, fockmat, energy);
  // pe_interface_energy(densmat, m_nbas, m_nnbas);
}

void CPPE::call_dynamic_response(const double* densmat, double* fockmat, double* energy) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_response(densmat, m_nbas, m_nnbas, fockmat, energy);
}


std::vector<Potential> CPPE::read_potfile(std::string potfile_name) {
  PotfileReader reader(potfile_name);
  std::vector<Potential> result = reader.read();
  assert(result.size() > 0);
  return result;
}

} // namespace libcppe
