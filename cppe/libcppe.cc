#include <cassert>
#include <iostream>

// #include "interface.h"

#include "libcppe.hh"

#include "core/electric_fields.hh"
#include "core/multipole_expansion.hh"
#include "utils/potfile_reader.hh"

namespace libcppe {

CPPE::CPPE() {
  m_gen1int_initialized = false;
  m_pe_initialized = false;
}

CPPE::~CPPE() {
  // if (m_gen1int_initialized && m_pe_initialized) {
  //   finalize_all();
  // }
}

// void CPPE::initialize_gen1int(int natoms, int nshells, const double *coords,
//                               const double *charges) {
//   gen1int_api_initialize(natoms, nshells, coords, charges);
//   m_gen1int_initialized = true;
// }
//
// void CPPE::gen1int_add_shell(int spher_gto, int idx_center,
//                              const double *coord_center, int ang_num,
//                              int num_prim, const double *exponents,
//                              int num_contr, const double *contr_coef) {
//   gen1int_api_add_shell(spher_gto, idx_center, coord_center, ang_num,
//   num_prim,
//                         exponents, num_contr, contr_coef);
// }
//
// void CPPE::gen1int_print_shells() {
//   if (!m_gen1int_initialized) {
//     throw std::runtime_error("Gen1int needs to be initialized before"
//                              "calling print_shells");
//   }
//   print_shells();
// }
//
// void CPPE::initialize_pelib(PeOptions &options, int natoms, int nbas,
//                             const double *coords, const double *charges) {
//   if (!m_gen1int_initialized) {
//     throw std::runtime_error("Gen1int needs to be initialized"
//                              "before calling PElib");
//   }
//   m_natoms = natoms;
//   m_nbas = nbas;
//   m_nnbas = nbas * (nbas + 1) / 2;
//   assert(options.potfile.size() < 80);
//   pe_set_potfile(options.potfile.c_str(),
//                  static_cast<int>(options.potfile.size()));
//   pe_set_border_options(options.pe_border, options.border_options.rmin,
//                         options.border_options.border_type,
//                         options.border_options.redist_pol
//                             ? -options.border_options.redist_order
//                             : options.border_options.redist_order,
//                         options.border_options.nredist);
//   pe_interface_init(natoms, coords, charges);
//   m_pe_initialized = true;
// }
//
// void CPPE::call_pe_energy(const double *densmat) {
//   if (!m_gen1int_initialized || !m_pe_initialized) {
//     throw std::runtime_error(
//         "Gen1int and PElib need to be initialized before calling PElib");
//   }
//   pe_interface_energy(densmat, m_nbas, m_nnbas);
// }
//
// void CPPE::call_pe_pol_energy(const double *densmat, double *energy) {
//   if (!m_gen1int_initialized || !m_pe_initialized) {
//     throw std::runtime_error(
//         "Gen1int and PElib need to be initialized before calling PElib");
//   }
//   pe_interface_pol_energy(densmat, m_nbas, m_nnbas, energy);
// }
//
// void CPPE::call_full_fock(const double *densmat, double *fockmat,
//                           double *energy) {
//   if (!m_gen1int_initialized || !m_pe_initialized) {
//     throw std::runtime_error(
//         "Gen1int and PElib need to be initialized before calling PElib");
//   }
//   pe_interface_fock(densmat, m_nbas, m_nnbas, fockmat, energy);
// }

std::vector<Potential> CPPE::read_potfile(std::string potfile_name) {
  PotfileReader reader(potfile_name);
  std::vector<Potential> result = reader.read();
  assert(result.size() > 0);
  return result;
}

}  // namespace libcppe
