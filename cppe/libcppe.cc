#include <iostream>
#include <cassert>

#include "interface.h"

#include "libcppe.hh"

namespace libcppe {

CPPE::CPPE() {
  // constructor for CPPE class
  std::cout << "Constructor CPPE" << std::endl;
}

CPPE::~CPPE() {
  // destructor for CPPE class
}

void CPPE::initialize_gen1int(int natoms, int nshells, const double* coords, const double *charges) {
  gen1int_api_initialize(natoms, nshells, coords, charges);
}

void CPPE::gen1int_add_shell(int spher_gto, int idx_center, const double* coord_center,
  int ang_num, int num_prim, const double* exponents, int num_contr, const double* contr_coef) {
  gen1int_api_add_shell(spher_gto, idx_center, coord_center, ang_num, num_prim, exponents, num_contr,
  contr_coef);
  m_gen1int_initialized = true;
}

void CPPE::gen1int_print_shells() {
  print_shells();
}

void CPPE::initialize_gen1int(int ntypes, int natoms, const double* coords, const double* charges,
                    const int* num_sym_atom, const int* shells_per_type, const int* max_l_per_type,
                    const int* num_cgto, const int* num_prim,
                    const int* num_contr,
                    const double* exponents,
                    const double* ucontr_coefs, int pure) {
  gen1int_api_create(ntypes, natoms, coords, charges, num_sym_atom, max_l_per_type, shells_per_type, num_cgto,
    num_prim, num_contr, exponents, ucontr_coefs, pure);
  m_gen1int_initialized = true;
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

void CPPE::call_pe_energy(const double* densmat) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_energy(densmat, m_nbas, m_nnbas);
}

void CPPE::call_full_fock(const double* densmat, double* fockmat, double* energy) {
  if (!m_gen1int_initialized || !m_pe_initialized) {
    throw std::runtime_error("Gen1int and PElib need to be initialized before calling PElib");
  }
  pe_interface_fock(densmat, m_nbas, m_nnbas, fockmat, energy);
  // pe_interface_energy(densmat, m_nbas, m_nnbas);
}

} // namespace libcppe
