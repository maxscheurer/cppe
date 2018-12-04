!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
!...  This file contains the data used to generate cube files.
!
!...  2012-03-09, Bin Gao
!...  * first version

#include "gen1int_host.h"

!> \brief data module of cube files
!> \author Bin Gao
!> \date 2012-03-09
module gen1int_cube

  implicit none

  logical, save, public :: do_density_cube = .false.                !electron density cube file
  logical, save, public :: do_homo_cube = .false.                   !HOMO cube file
  logical, save, public :: do_lumo_cube = .false.                   !LUMO cube file
  logical, save, public :: do_mo_cube = .false.                     !MO cube files
  integer, save, public :: num_cube_mo = 0                          !number of MOs to generate
  integer, save, allocatable, public :: idx_cube_mo(:)              !indices of MOs
  character(MAX_LEN_STR), save, public :: cube_format = "GAUSSIAN"  !format of cube file
  real(REALK), save, public :: cube_origin(3) = 0.0_REALK           !origin of cube file
  real(REALK), save, public :: cube_increment(3,3) = 0.0_REALK      !increments of cube file
  integer, save, public :: cube_num_inc(3) = 0                      !number of increments of cube file

end module gen1int_cube
