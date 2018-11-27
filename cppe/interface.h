#ifndef CPPE_INTERFACE_H
#define CPPE_INTERFACE_H

// extern "C" avoids C++ mangling of names.
extern "C" void pe_interface_init_(int *, const double *, const double *);
extern "C" void pe_set_potfile_(const char *, int);
extern "C" void pe_interface_energy_(const double *, int *, int *);
extern "C" void pe_interface_pol_energy_(const double *, int *, int *,
                                         double *);
extern "C" void pe_interface_fock_(const double *, int *, int *, const double *,
                                   double *);
// extern "C" void pe_interface_response_(const double*, int*, int*, const
// double*, double*);
extern "C" int gen1int_api_initialized_();
extern "C" void print_shells_();
extern "C" void finalize_all_();

extern "C" void gen1int_api_initialize_(int *natoms, int *num_shells,
                                        const double *coords,
                                        const double *charges);
extern "C" void gen1int_create_shell_(int *spher_gto, int *idx_center,
                                      const double *coord_center, int *ang_num,
                                      int *num_prim, const double *exponents,
                                      int *num_contr, const double *contr_coef);

extern "C" void pe_set_border_options_(int *m_pe_border, double *m_rmin,
                                       int *type_flag, int *redist_order,
                                       int *nredist);

void pe_interface_init(int n, const double *coords, const double *charges) {
  return pe_interface_init_(&n, coords, charges);
}
void pe_set_potfile(const char *c, int n) { return pe_set_potfile_(c, n); }
void pe_interface_energy(const double *densmat, int ndim, int nnbas) {
  return pe_interface_energy_(densmat, &ndim, &nnbas);
}

void pe_interface_fock(const double *densmat, int ndim, int nnbas,
                       const double *fckmat, double *energy) {
  return pe_interface_fock_(densmat, &ndim, &nnbas, fckmat, energy);
}

void pe_interface_pol_energy(const double *densmat, int ndim, int nnbas,
                             double *energy) {
  return pe_interface_pol_energy_(densmat, &ndim, &nnbas, energy);
}

// void pe_interface_response(const double *densmat, int ndim, int nnbas, const
// double *fckmat, double* energy) {
//   return pe_interface_response_(densmat, &ndim, &nnbas, fckmat, energy);
// }

int gen1int_api_initialized() { return gen1int_api_initialized_(); }

void gen1int_api_initialize(int natoms, int num_shells, const double *coords,
                            const double *charges) {
  return gen1int_api_initialize_(&natoms, &num_shells, coords, charges);
}

void gen1int_api_add_shell(int spher_gto, int idx_center,
                           const double *coord_center, int ang_num,
                           int num_prim, const double *exponents, int num_contr,
                           const double *contr_coef) {
  return gen1int_create_shell_(&spher_gto, &idx_center, coord_center, &ang_num,
                               &num_prim, exponents, &num_contr, contr_coef);
}

void print_shells() { return print_shells_(); }

void finalize_all() { return finalize_all_(); }

void pe_set_border_options(int m_pe_border, double m_rmin, int type_flag,
                           int redist_order, int nredist) {
  return pe_set_border_options_(&m_pe_border, &m_rmin, &type_flag,
                                &redist_order, &nredist);
}

#endif // CPPE_INTERFACE_H
