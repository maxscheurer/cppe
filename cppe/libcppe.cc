#include <iostream>

#include "interface.h"

#include "libcppe.hh"
// Dalton variables for gen1int API
#define KATOM 500
#define MXCENT 500
#define KANG 17
#define KBLOCK 18000
#define KPRIM 35

namespace libcppe {

CPPE::CPPE() {
  // constructor for CPPE class
  std::cout << "Constructor CPPE" << std::endl;
  // test
}

CPPE::~CPPE() {
  // destructor for CPPE class
}

void CPPE::call_pe_energy(const double* densmat) {
  pe_interface_energy(densmat, m_nbas, m_nnbas);
}

} // namespace libcppe
