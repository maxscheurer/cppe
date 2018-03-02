// Dalton variables for gen1int API
#define KATOM 500
#define MXCENT 500
#define KANG 17
#define KBLOCK 18000
#define KPRIM 35


namespace libcppe {

  class CPPE {
    public:
      CPPE();
      ~CPPE();
      void call_pe_energy(const double *densmat);
      void call_full_fock(const double* densmat, double* fockmat, double* energy);
      bool gen1int_initialized() { return m_gen1int_initialized; }
      void initialize_gen1int(int ntypes, int natoms, const double* coords, const double* charges,
                          const int* num_sym_atom, const int* shells_per_type, const int* max_l_per_type,
                          const int* num_cgto, const int* num_prim,
                          const int* num_contr,
                          const double* exponents,
                          const double* ucontr_coefs, int pure);
      void initialize_gen1int(int natoms, int nshells, const double* coords, const double *charges);
      void gen1int_add_shell(int spher_gto, int idx_center, const double* coord_center,
        int ang_num, int num_prim, const double* exponents, int num_contr, const double* contr_coef);
      void gen1int_print_shells();
      void initialize_pelib(std::string potfile, int natoms, int nbas,
                            const double* coords, const double* charges);

    private:
      bool m_gen1int_initialized;
      // TODO: build a struct/class with the basis set info?
      int m_nbas;  // number of basis functions
      int m_nnbas; // nbas*(nbas+1)/2
      int m_natoms;
      bool m_pe_initialized;
  };

}
