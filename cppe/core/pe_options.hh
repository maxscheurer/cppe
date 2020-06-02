#pragma once

#include <string>
#include <vector>

namespace libcppe {

struct PeOptions {
  std::string potfile{"potential.pot"};  ///< Name of the potential file

  bool iso_pol = false;  ///< Make polarizabilities isotropic

  // bool zero_pol = false;
  // bool zero_mul = false;
  // int zero_mul_order = 1;

  double induced_thresh = 1e-8;  ///< Threshold for induced moments convergence
  int maxiter           = 50;    ///< Maximum number of iterations for induced moments

  bool damp_induced             = false;   ///< Enable Thole damping for induced moments
  double damping_factor_induced = 2.1304;  ///< damping factor for induced moments
  bool damp_multipole =
        false;  ///< Enable Thole damping for electric fields created by multipole moments
  double damping_factor_multipole = 2.1304;  ///< damping factor for electric
                                             ///< fields created by multipole moments

  bool pe_border = false;             ///< Activate border options for sites in proximity
                                      ///< to the QM/MM border
  std::string border_type{"remove"};  ///< border type, either "remove" or "redist"
                                      ///< moments/polarizabilities
  double border_rmin = 2.2;  ///< minimum radius from QM atoms to MM sites to be taken
                             ///< into account for removal/redistribution (in AU)

  int border_nredist = -1;  ///< number of neighbor sites to redistribute to. The default
                            ///< (-1) redistributes to all sites which are not in the
                            ///< border region order to which multipole moments are
                            ///< redistributed moments > redist_order will be removed
  int border_redist_order = 1;     ///< order from which moments are removed, e.g., if set
                                   ///< to 1 (default), only charges are redistributed and
                                   ///< all higher order moments are removed
  bool border_redist_pol = false;  ///< redistribute polarizabilities? If false,
  ///< polarizabilities are removed (default)

  // TODO: remove these options completely
  bool do_diis = true;  ///< Use DIIS acceleration to obtain induced moments
  double diis_start_norm =
        1.0;  ///< maximal residual norm for which DIIS is being enabled
};

static const std::vector<std::string> valid_option_keys{"potfile",
                                                        "iso_pol",
                                                        "induced_thresh",
                                                        "maxiter",
                                                        "damp_induced",
                                                        "damping_factor_induced",
                                                        "damp_multipole",
                                                        "damping_factor_multipole",
                                                        "pe_border",
                                                        "border_type",
                                                        "border_rmin",
                                                        "border_nredist",
                                                        "border_redist_order",
                                                        "border_redist_pol"};

}  // namespace libcppe
