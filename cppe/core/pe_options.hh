#pragma once

#include <string>

namespace libcppe {

enum BorderType { rem, redist };

struct BorderOptions {
  BorderType border_type{rem};  ///< border type, either remove or redistribute
                                ///< moments/polarizabilities
  double rmin = 2.2;            ///< minimum radius from QM atoms to MM sites to be taken
                                ///< into account for removal/redistribution (in AU)

  int nredist = -1;  ///< number of neighbor sites to redistribute to. The default (-1)
                     ///< redistributes to all sites which are not in the border region

  // order to which multipole moments are redistributed
  // moments > redist_order will be removed
  int redist_order = 1;     ///< order from which moments are removed, e.g., if set
                            ///< to 1 (default), only charges are redistributed and
                            ///< all higher order moments are removed
  bool redist_pol = false;  ///< redistribute polarizabilities? If false,
                            ///< polarizabilities are removed (default)
};

struct PeOptions {
  std::string potfile{"potential.pot"};  ///< Name of the potential file

  bool iso_pol = false;  ///< Make polarizabilities isotropic

  // bool zero_pol = false;
  // bool zero_mul = false;
  // int zero_mul_order = 1;

  double induced_thresh = 1e-8;  ///< Threshold for induced moments convergence
  bool do_diis          = true;  ///< Use DIIS acceleration to obtain induced moments
  int maxiter           = 50;    ///< Maximum number of iterations for induced moments
  double diis_start_norm =
        1.0;  ///< maximal residual norm for which DIIS is being enabled

  bool pe_border = false;          ///< Activate border options for sites in proximity
                                   ///< to the QM/MM border
  BorderOptions border_options{};  ///< Options for QM/MM border
};

}  // namespace libcppe
