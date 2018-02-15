#ifndef INTERFACE_H
#define INTERFACE_H

// extern "C" avoids C++ mangling of names.
extern "C" void pe_interface_init_(int*, double*, double*);
extern "C" void pe_set_potfile_(const char*, int);
extern "C" void pe_interface_energy_(const double*, int*, int*);
extern "C" void pe_interface_fock_(double*, int*, int*, double*, double*);
extern "C" int gen1int_api_initialized_();
extern "C" void gen1int_api_create_(int*, int*, int*, int*, int*, int*, int*, double*, double*, int*);
extern "C" void set_coord_nuclei_(int*, double*, double*);
// As usually FORTRAN passes by reference => use pointers


void pe_interface_init(int n, double* coords, double* charges) { return pe_interface_init_(&n, coords, charges); }
void pe_set_potfile(const char* c, int n) { return pe_set_potfile_(c, n); }
void pe_interface_energy(const double *densmat, int ndim, int nnbas) {
  return pe_interface_energy_(densmat, &ndim, &nnbas);
}

void pe_interface_fock(double *densmat, int ndim, int nnbas, double *fckmat, double energy) {
  return pe_interface_fock_(densmat, &ndim, &nnbas, fckmat, &energy);
}

int gen1int_api_initialized() { return gen1int_api_initialized_(); }

void gen1int_api_create(int n, int* num_sym_atom, int* max_l_per_type, int* shells_per_type, int* num_cgto, int* num_prim, int* num_contr, double* exponents, double* ucontr_coefs, int pure) {
  return gen1int_api_create_(&n, num_sym_atom, max_l_per_type, shells_per_type, num_cgto, num_prim, num_contr, exponents, ucontr_coefs, &pure);
}

void set_coord_nuclei(int natoms, double* coords, double* charges) {
  return set_coord_nuclei_(&natoms, coords, charges);
}

#endif
