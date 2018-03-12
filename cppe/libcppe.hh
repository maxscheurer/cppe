#include <string>
#include <vector>

#include "utils/potfile_reader.hh"
#include "core/multipole.hh"


namespace libcppe {

  class CPPE {
    public:
      CPPE();
      ~CPPE();
      void call_pe_energy(const double *densmat);
      void call_pe_pol_energy(const double *densmat, double* energy);
      void call_full_fock(const double* densmat, double* fockmat, double* energy);
      void call_dynamic_response(const double* densmat, double* fockmat, double* energy);
      bool gen1int_initialized() { return m_gen1int_initialized; }
      void initialize_gen1int(int natoms, int nshells, const double* coords, const double *charges);
      void gen1int_add_shell(int spher_gto, int idx_center, const double* coord_center,
        int ang_num, int num_prim, const double* exponents, int num_contr, const double* contr_coef);
      void gen1int_print_shells();
      void initialize_pelib(std::string potfile, int natoms, int nbas,
                            const double* coords, const double* charges);
      
      std::vector<Potential> read_potfile(std::string potfile_name);

    private:
      bool m_gen1int_initialized;
      // TODO: build a struct/class with the basis set info?
      int m_nbas;  // number of basis functions
      int m_nnbas; // nbas*(nbas+1)/2
      int m_natoms;
      bool m_pe_initialized;
  };

}
