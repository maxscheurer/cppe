#ifndef LIBCPPE_INCLUDE_H
#define LIBCPPE_INCLUDE_H

#include <string>
#include <vector>

#include <armadillo>

#include "core/multipole.hh"
#include "core/molecule.hh"


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
      double calculate_nulcei_multipole_interaction(Molecule& mol, std::vector<Potential>& );
      arma::vec calculate_nulcear_fields(Molecule& mol, std::vector<Potential>& potentials, size_t polarizable_sites);
      arma::vec calculate_multipole_fields(std::vector<Potential>& potentials, size_t polarizable_sites);
      arma::vec calculate_induced_moments(std::vector<Potential>& potentials, arma::vec& total_fields);
      
    private:
      bool m_gen1int_initialized;
      bool m_pe_initialized;
      // TODO: build a struct/class with the basis set info?
      int m_nbas;  // number of basis functions
      int m_nnbas; // nbas*(nbas+1)/2
      int m_natoms;
  };

}

#endif // LIBCPPE_INCLUDE_H