#ifndef INTERFACE_H
#define INTERFACE_H

// extern "C" avoids C++ mangling of names.
extern "C" void pe_interface_init_(int*, const double*, const double*);
extern "C" void pe_set_potfile_(const char*, int);
extern "C" void pe_interface_energy_(const double*, int*, int*);
extern "C" void pe_interface_pol_energy_(const double*, int*, int*, double*);
extern "C" void pe_interface_fock_(const double*, int*, int*, const double*, double*);
extern "C" void pe_interface_response_(const double*, int*, int*, const double*, double*);
extern "C" int gen1int_api_initialized_();
extern "C" void gen1int_api_create_(int*, int*, const double* coords, const double* charges, const int*, const int*, const int*, const int*, const int*, const int*, const double*, const double*, int*);
extern "C" void set_coord_nuclei_(int*, double*, double*);
extern "C" void print_shells_();

extern "C" void gen1int_api_initialize_(int *natoms, int *num_shells, const double *coords, const double *charges);
extern "C" void gen1int_create_shell_(int* spher_gto, int* idx_center, const double* coord_center,
  int* ang_num, int* num_prim, const double* exponents, int* num_contr, const double* contr_coef);
// As usually FORTRAN passes by reference => use pointers


void pe_interface_init(int n, const double* coords, const double* charges) { return pe_interface_init_(&n, coords, charges); }
void pe_set_potfile(const char* c, int n) { return pe_set_potfile_(c, n); }
void pe_interface_energy(const double *densmat, int ndim, int nnbas) {
  return pe_interface_energy_(densmat, &ndim, &nnbas);
}

void pe_interface_fock(const double *densmat, int ndim, int nnbas, const double *fckmat, double* energy) {
  return pe_interface_fock_(densmat, &ndim, &nnbas, fckmat, energy);
}

void pe_interface_pol_energy(const double *densmat, int ndim, int nnbas, double* energy) {
  return pe_interface_pol_energy_(densmat, &ndim, &nnbas, energy);
}

void pe_interface_response(const double *densmat, int ndim, int nnbas, const double *fckmat, double* energy) {
  return pe_interface_response_(densmat, &ndim, &nnbas, fckmat, energy);
}

int gen1int_api_initialized() { return gen1int_api_initialized_(); }

void gen1int_api_create(int n, int natoms, const double* coords, const double* charges, const int* num_sym_atom, const int* max_l_per_type, const int* shells_per_type, const int* num_cgto, const int* num_prim, const int* num_contr, const double* exponents, const double* ucontr_coefs, int pure) {
  return gen1int_api_create_(&n, &natoms, coords, charges, num_sym_atom, max_l_per_type, shells_per_type, num_cgto, num_prim, num_contr, exponents, ucontr_coefs, &pure);
}

void set_coord_nuclei(int natoms, double* coords, double* charges) {
  return set_coord_nuclei_(&natoms, coords, charges);
}

void gen1int_api_initialize(int natoms, int num_shells, const double *coords, const double *charges) {
  return gen1int_api_initialize_(&natoms, &num_shells, coords, charges);
}

void gen1int_api_add_shell(int spher_gto, int idx_center, const double* coord_center,
  int ang_num, int num_prim, const double* exponents, int num_contr, const double* contr_coef) {
  return gen1int_create_shell_(&spher_gto, &idx_center, coord_center, &ang_num, &num_prim, exponents, &num_contr, contr_coef);
}

void print_shells() {
  return print_shells_();
}

#endif
