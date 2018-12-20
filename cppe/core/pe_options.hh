#ifndef CPPE_CORE_PE_OPTIONS_H
#define CPPE_CORE_PE_OPTIONS_H

#include <string>
namespace libcppe {

enum BorderType { rem, redist };

struct BorderOptions {
  BorderType border_type{rem};
  // rmin is in AU!
  double rmin = 2.2;

  // number of sites that parameters are redistributed to, i.e., the
  // nredist nearest neighbors
  int nredist = 1;

  // order to which multipole moments are redistributed
  // moments > redist_order will be removed
  int redist_order = 1;
  bool redist_pol = false;
};

struct PeOptions {
  std::string potfile{"potential.pot"};

  int print_level = 1;

  bool damp_induced = false;
  bool damp_multipoles = false;
  bool damp_core = false;

  double damp_coeff_ind = 2.1304;
  double damp_coeff_mult = 2.1304;
  double damp_coeff_core = 2.1304;

  bool zero_pol = false;
  bool zero_mul = false;
  int zero_mul_order = 1;

  int induced_thresh = 8;
  bool do_diis = true;
  int diis_maxiter = 50;
  double diis_start_norm = 1.0;

  bool pe_border = false;
  BorderOptions border_options{};
};

}  // namespace libcppe
#endif  // CPPE_CORE_PE_OPTIONS_H
