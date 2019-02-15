#ifndef CPPE_CORE_PE_OPTIONS_H
#define CPPE_CORE_PE_OPTIONS_H

#include <string>
namespace libcppe {

enum BorderType { rem, redist };

struct BorderOptions {
  BorderOptions() : border_type(rem),
                    rmin(2.2), nredist(1), redist_order(1), redist_pol(false) {};
  //BorderType border_type{rem};
  BorderType border_type;
  // rmin is in AU!
  double rmin;

  // number of sites that parameters are redistributed to, i.e., the
  // nredist nearest neighbors
  int nredist;

  // order to which multipole moments are redistributed
  // moments > redist_order will be removed
  int redist_order;
  bool redist_pol;
};

struct PeOptions {
  PeOptions() : potfile("potential.pot"), print_level(1), damp_induced(false),
                damp_multipoles(false), damp_core(false), damp_coeff_ind(2.1304),
                damp_coeff_mult(2.1304), damp_coeff_core(2.1304), zero_pol(false),
                zero_mul(false), zero_mul_order(1), induced_thresh(8),
                do_diis(true), diis_maxiter(50), diis_start_norm(1.0), pe_border(false),
                border_options(BorderOptions()) {};
  std::string potfile;

  int print_level;

  bool damp_induced;
  bool damp_multipoles;
  bool damp_core;

  double damp_coeff_ind;
  double damp_coeff_mult;
  double damp_coeff_core;

  bool zero_pol;
  bool zero_mul;
  int zero_mul_order;

  int induced_thresh;
  bool do_diis;
  int diis_maxiter;
  double diis_start_norm;

  bool pe_border;
  BorderOptions border_options;
};

}  // namespace libcppe
#endif  // CPPE_CORE_PE_OPTIONS_H
