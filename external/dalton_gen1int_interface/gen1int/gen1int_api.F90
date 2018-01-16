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
!...  This file is the module of API of Gen1Int interface.
!
!...  2013-05-16, Bin Gao
!...  * add a subroutine to update the information of molecule in this interface
!
!...  2013-05-03, Bin Gao
!...  * fix the bug of returning wrong partial geometric derivatives
!
!...  2012-05-09, Radovan Bast
!...  * implements large and small components
!
!...  2012-01-10, Bin Gao
!...  * first version

#include "gen1int_host.h"

!> \brief module of API of Gen1Int interface
!> \author Bin Gao and Radovan Bast
!> \date 2012-01-10
module gen1int_api

  ! Fortran 90 module of Gen1Int library
  use gen1int_geom, Gen1IntAPINaryTreeDestroy => NaryTreeDestroy, &
                    Gen1IntAPINaryTreeView => NaryTreeView
  ! AO sub-shells
  use gen1int_shell
  ! matrix module
  use gen1int_matrix

  implicit none

  ! one-electron operator with non-zero components
  type, public :: prop_comp_t
    private
    ! one-electron operator created by Gen1Int library
    type(one_prop_t) one_prop
    ! non-zero components, the first dimension is 2 for bra and ket sub-shells
    integer, allocatable :: nnz_comp(:,:)
  end type

  ! if Gen1Int interface initialized
  logical, save, private :: api_inited = .false.

  ! large and small components (the latter is only used in Dirac by -DPRG_DIRAC)
#ifdef PRG_DIRAC
  integer, parameter, public :: NUM_COMPONENTS = 2
  integer, parameter, public :: LARGE_COMP = 1
  integer, parameter, public :: SMALL_COMP = 2
#else
  integer, parameter, public :: NUM_COMPONENTS = 1
  integer, parameter, public :: LARGE_COMP = 1
#endif

  ! number of AO sub-shells from Dalton/DIRAC
  ! Dalton only uses num_sub_shells(1), DIRAC uses either the first or both
  integer, save, private :: num_sub_shells(NUM_COMPONENTS) = 0
  ! AO sub-shells from Dalton/DIRAC
  ! Dalton only uses sub_shells(:,1), DIRAC uses either only sub_shells(:,1) or both sub_shells(:,1:2)
  type(sub_shell_t), save, allocatable, private :: sub_shells(:,:)

  ! molecule information used when creating operator and generating cube files
  integer, save, private :: api_num_atoms = 0                      !number of atoms
  real(REALK), save, allocatable, private :: api_coord_atoms(:,:)  !coordinates of atoms
  real(REALK), save, allocatable, private :: api_charge_atoms(:)   !charges of atoms
  real(REALK), save, private :: api_dipole_origin(3) = 0.0_REALK   !coordinates of dipole origin
  real(REALK), save, private :: api_gauge_origin(3) = 0.0_REALK    !coordinates of gauge origin
  real(REALK), save, private :: api_origin_LPF(3) = 0.0_REALK      !coordinates of origin of London phase factor

  ! \fn(Gen1IntAPICreate) might be the only program specific subroutine (depends on common blocks)
  public :: Gen1IntAPICreate
#if defined(VAR_MPI)
  public :: Gen1IntAPIBcast
#endif
  public :: Gen1IntAPIInited
  public :: Gen1IntAPIUpdateMolecule
  public :: Gen1IntAPIShellView
  public :: Gen1IntAPIGetNumAO
  public :: Gen1IntAPIGetMO
  public :: Gen1IntAPIDestroy

  public :: Gen1IntAPIPropCreate
  public :: Gen1IntAPIPropView
  public :: Gen1IntAPIPropGetNumProp
  public :: Gen1IntAPIPropGetSymmetry
  public :: Gen1IntAPIPropGetIntExpt
  public :: Gen1IntAPIPropGetFunExpt
  public :: Gen1IntAPIPropDestroy

  public :: Gen1IntAPINaryTreeCreate
  public :: Gen1IntAPINaryTreeView
  public :: Gen1IntAPINaryTreeDestroy

  public :: Gen1IntOnePropGetIntExpt

  contains

  !> \brief initializes Gen1Int interface, for instance, creates the AO sub-shells
  !>        of host program (based on \fn(ORBPRO) subroutine by getting the unnormalized
  !>        contraction coefficients); should be called before any calculation
  !> \author Bin Gao and Radovan Bast
  !> \date 2010-12-06
  !> \param num_comp is the number of components
  !> \param num_atom_type is the number of atomic types
  !> \param num_sym_atom contains the number of symmetry independent centers of atomic types
  !> \param ang_numbers contains the angular momentum (1=s, 2=p, 3=d, ...)
  !> \param num_cgto contains the number of CGTOs in the AO blocks for an angular momentum
  !> \param num_prim contains the number of uncontracted functions
  !> \param num_contr contains the number of contracted functions
  !> \param exponents contains the exponents of primitive shells
  !> \param ucontr_coefs contains the unnormalized contraction coefficients
  !> \note this subroutine is program specific; please add the meaning of other
  !>       arguments if you know, thanks!
  subroutine Gen1IntAPICreate(num_comp, num_atom_type, KATOM, num_sym_atom, &
                              ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                              num_prim, num_contr, KPRIM, exponents, ucontr_coefs)
    integer, intent(in) :: num_comp
    integer, intent(in) :: num_atom_type
    integer, intent(in) :: KATOM
    integer, intent(in) :: num_sym_atom(KATOM)
    integer, intent(in) :: ang_numbers(KATOM,num_comp)
    integer, intent(in) :: NBLCK(KATOM,num_comp)
    integer, intent(in) :: KANG
    integer, intent(in) :: num_cgto(KANG,KATOM,num_comp)
    integer, intent(in) :: KBLOCK
    integer, intent(in) :: num_prim(KBLOCK,num_comp)
    integer, intent(in) :: num_contr(KBLOCK,num_comp)
    integer, intent(in) :: KPRIM
    real(REALK), intent(in) :: exponents(KPRIM,KBLOCK,num_comp)
    real(REALK), intent(in) :: ucontr_coefs(KPRIM,KPRIM,KBLOCK,num_comp)
#include "mxcent.h"
#include "maxaqn.h"
#include "ccom.h"
#include "nuclei.h"
#include "orgcom.h"
    integer icomp      !incremental recorder over components
    integer IDX_CENT   !index of symmetry independent center
    integer IDX_BLOCK  !
    integer ITYP       !incremental recorder over number of atomic types
    integer ICENT      !incremental recorder over number of symmetry independent centers
    integer IANG       !incremental recorder over angular momentum
    integer KBCH       !
    logical spher_gto  !if SGTOs
    integer ishell     !incremental recorder over AO sub-shells
    integer ang_num    !angular number
    real(REALK), allocatable :: contr_coef(:,:)  !contraction coefficients
    integer icontr, iprim                        !incremental recorder over contractions
    integer ierr                                 !error information
    integer iang_sub

    if (num_comp>NUM_COMPONENTS) then
      call quit("Gen1IntAPICreate>> too many components!")
    end if
    ! terminates previous created interface
    if (api_inited) call Gen1IntAPIDestroy()
    ! gets the number of AO sub-shells and kind of GTOs
    num_sub_shells = 0
    spher_gto = .false.
    ! loops over components of basis sets
    do icomp = 1, num_comp
      ! number of atomic types
      do ITYP = 1, num_atom_type
        ! number of symmetry independent centers of this type
        do ICENT = 1, num_sym_atom(ITYP)
          ! angular momentum 1=s, 2=p, 3=d, etc.
          do IANG = 1, ang_numbers(ITYP,icomp)
             ! radovan: basis does not have to start with s
             !          the blocks do not have to be s p d f
             !          they can be s p p d f
             num_sub_shells(icomp) = num_sub_shells(icomp) &
                                   + num_cgto(IANG, ITYP, icomp)
             if (SPH(IANG)) spher_gto = .true.  !Dalton always marks s and p sub-shells as CGTOs
          end do
         end do
      end do
    end do
    ! initializes the AO sub-shells
    allocate(sub_shells(maxval(num_sub_shells),NUM_COMPONENTS), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntAPICreate>> failed to allocate sub_shells!")
    end if
    ! loops over components of basis sets
    do icomp = 1, num_comp
      ishell = 0
      IDX_BLOCK = 0
      IDX_CENT = 0
      ! number of atomic types
      do ITYP = 1, num_atom_type
        ! number of symmetry independent centers of this type
        do ICENT = 1, num_sym_atom(ITYP)
          IDX_CENT = IDX_CENT+1
          KBCH = IDX_BLOCK
          ! angular momentum 1=s, 2=p, 3=d, etc.
          do IANG = 1, ang_numbers(ITYP,icomp)

            ! radovan: basis does not have to start with s
            !          and there can be two consecutive blocks
            !          with the same IANG

            do iang_sub = 1, num_cgto(IANG, ITYP, icomp)
               ! next block
               KBCH = KBCH+1
               ! gets the contraction coefficients
               allocate(contr_coef(num_contr(KBCH,icomp),num_prim(KBCH,icomp)), stat=ierr)
               if (ierr/=0) then
                 call quit("Gen1IntAPICreate>> failed to allocate contr_coef!")
               end if
               do iprim = 1, num_prim(KBCH,icomp)
                 do icontr = 1, num_contr(KBCH,icomp)
                   contr_coef(icontr,iprim) = ucontr_coefs(iprim,icontr,KBCH,icomp)
                 end do
               end do
               ! normalizes the contraction coefficients
               ang_num = IANG-1
               ! Dalton/Dirac do not use mixed CGTOs and SGTOs
               !if (SPH(IANG).neqv.spher_gto) then
               !  stop "Gen1IntAPICreate>> mixed CGTOs and SGTOs not supported!"
               !end if
               ! Gen1Int library uses HGTOs in recurrence relations, while Dalton uses
               ! CGTOs, so that we need to normalize the contraction coefficients of
               ! SGTOs using Gen1Int subroutines
               if (spher_gto) then
                 call norm_contr_sgto(ang_num, num_prim(KBCH,icomp),                &
                                      exponents(1:num_prim(KBCH,icomp),KBCH,icomp), &
                                      num_contr(KBCH,icomp), contr_coef)
               else
                 call norm_contr_cgto(ang_num, num_prim(KBCH,icomp),                &
                                      exponents(1:num_prim(KBCH,icomp),KBCH,icomp), &
                                      num_contr(KBCH,icomp), contr_coef)
               end if
               !-ISTBNU(IDX_CENT)  !stabiliser: basic sym. op. that do not move center
               ishell = ishell+1
               if (ishell>1) then
                 call Gen1IntShellCreate(spher_gto=spher_gto,                        &
                                         idx_cent=IDX_CENT,                          &
                                         coord_cent=CORD(1:3,IDX_CENT),              &
                                         ang_num=ang_num,                            &
                                         num_prim=num_prim(KBCH,icomp),              &
                                         exponents=exponents(1:num_prim(KBCH,icomp), &
                                                             KBCH,icomp),            &
                                         num_contr=num_contr(KBCH,icomp),            &
                                         contr_coef=contr_coef,                      &
                                         last_shell=sub_shells(ishell-1,icomp),      &
                                         sub_shell=sub_shells(ishell,icomp))
               ! sets the first AO sub-shell
               else
#ifdef PRG_DIRAC
                 select case (icomp)
                 case (LARGE_COMP)
#endif         
                   call Gen1IntShellCreate(spher_gto=spher_gto,                        &
                                           idx_cent=IDX_CENT,                          &
                                           coord_cent=CORD(1:3,IDX_CENT),              &
                                           ang_num=ang_num,                            &
                                           num_prim=num_prim(KBCH,icomp),              &
                                           exponents=exponents(1:num_prim(KBCH,icomp), &
                                                               KBCH,icomp),            &
                                           num_contr=num_contr(KBCH,icomp),            &
                                           contr_coef=contr_coef,                      &
                                           sub_shell=sub_shells(ishell,icomp))
#ifdef PRG_DIRAC
                 case (SMALL_COMP)
                   ! put small component block behind large component block
                   call Gen1IntShellCreate(spher_gto=spher_gto,                              &
                                           idx_cent=IDX_CENT,                                &
                                           coord_cent=CORD(1:3,IDX_CENT),                    &
                                           ang_num=ang_num,                                  &
                                           num_prim=num_prim(KBCH,icomp),                    &
                                           exponents=exponents(1:num_prim(KBCH,icomp),       &
                                                               KBCH,icomp),                  &
                                           num_contr=num_contr(KBCH,icomp),                  &
                                           contr_coef=contr_coef,                            &
                                           last_shell=sub_shells(num_sub_shells(LARGE_COMP), &
                                                                 LARGE_COMP),                &
                                           sub_shell=sub_shells(ishell,icomp))
                 end select
#endif         
               end if
               deallocate(contr_coef)
            end do
          end do
        end do
        IDX_BLOCK = IDX_BLOCK+NBLCK(ITYP,icomp)
      end do
    end do
    ! number of atoms
    api_num_atoms = NUCDEP
    ! coordinates of atoms
    allocate(api_coord_atoms(3,api_num_atoms), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntAPICreate>> failed to allocate api_coord_atoms!")
    end if
    api_coord_atoms = CORD(:,1:NUCDEP)
    ! charges of atoms
    allocate(api_charge_atoms(api_num_atoms), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntAPICreate>> failed to allocate api_charge_atoms!")
    end if
    api_charge_atoms = -CHARGE(1:NUCDEP)
    ! coordinates of origins
    api_dipole_origin = DIPORG
    api_gauge_origin = GAGORG
    api_origin_LPF = ORIGIN
    api_inited = .true.
  end subroutine Gen1IntAPICreate

#if defined(VAR_MPI)
  !> \brief broadcasts AO sub-shells
  !> \author Bin Gao
  !> \date 2012-05-13
  !> \param root is the root processor which broadcasts the AO sub-shells
  !> \param api_comm is the MPI communicator
  subroutine Gen1IntAPIBcast(root, api_comm)
    integer, intent(in) :: root
    integer, intent(in) :: api_comm
#include "mpif.h"
    integer rank_proc  !rank of processor
    integer icomp      !incremental recorder over components
    integer ierr       !error information
    ! gets the rank of processor
    call MPI_Comm_rank(api_comm, rank_proc, ierr)
    if (rank_proc==root) then
      ! stops if the interface is not initialized
      if (.not.api_inited) stop "Gen1IntAPIBcast>> interface is not initialized!"
      ! broadcasts the number of components
      call MPI_Bcast(num_sub_shells, NUM_COMPONENTS, MPI_INTEGER, root, api_comm, ierr)
    else
      ! terminates previous created interface
      if (api_inited) call Gen1IntAPIDestroy()
      ! gets the number of components
      call MPI_Bcast(num_sub_shells, NUM_COMPONENTS, MPI_INTEGER, root, api_comm, ierr)
      ! allocates memory for sub-shells
      allocate(sub_shells(maxval(num_sub_shells),NUM_COMPONENTS), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntAPIBcast>> failed to allocate sub_shells!")
      end if
    end if
    ! broadcasts sub-shells
    do icomp = 1, NUM_COMPONENTS
      if (num_sub_shells(icomp)>0)                               &
        call Gen1IntShellBcast(num_shells=num_sub_shells(icomp), &
                               sub_shells=sub_shells(:,icomp),   &
                               root=root,                        &
                               api_comm=api_comm)
    end do
    ! number of atoms
    call MPI_Bcast(api_num_atoms, 1, MPI_INTEGER, root, api_comm, ierr)
    ! coordinates and charges of atoms
    if (rank_proc==root) then
      call MPI_Bcast(api_coord_atoms, 3*api_num_atoms, MPI_REALK, root, api_comm, ierr)
      call MPI_Bcast(api_charge_atoms, api_num_atoms, MPI_REALK, root, api_comm, ierr)
    else
      allocate(api_coord_atoms(3,api_num_atoms), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntAPIBcast>> failed to allocate api_coord_atoms!")
      end if
      call MPI_Bcast(api_coord_atoms, 3*api_num_atoms, MPI_REALK, root, api_comm, ierr)
      allocate(api_charge_atoms(api_num_atoms), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntAPIBcast>> failed to allocate api_charge_atoms!")
      end if
      call MPI_Bcast(api_charge_atoms, api_num_atoms, MPI_REALK, root, api_comm, ierr)
    end if
    ! coordinates of origins
    call MPI_Bcast(api_dipole_origin, 3, MPI_REALK, root, api_comm, ierr)
    call MPI_Bcast(api_gauge_origin, 3, MPI_REALK, root, api_comm, ierr)
    call MPI_Bcast(api_origin_LPF, 3, MPI_REALK, root, api_comm, ierr)
    api_inited = .true.
  end subroutine Gen1IntAPIBcast
#endif

  !> \brief returns if the interface is initialized or not
  !> \author Bin Gao
  !> \date 2012-03-09
  function Gen1IntAPIInited() result(api_status)
    logical api_status
    api_status = api_inited
  end function Gen1IntAPIInited

  !> \brief visualizes the information of AO sub-shells in host programs
  !> \author Bin Gao and Radovan Bast
  !> \date 2012-03-09
  !> \param io_viewer is the logical unit number of the viewer
  subroutine Gen1IntAPIShellView(io_viewer)
    integer, intent(in) :: io_viewer
    integer icomp  !incremental recorder over components
    if (api_inited) then
      do icomp = 1, NUM_COMPONENTS
        if (num_sub_shells(icomp)>0)                              &
          call Gen1IntShellView(num_shells=num_sub_shells(icomp), &
                                sub_shells=sub_shells(:,icomp),   &
                                io_viewer=io_viewer)
      end do
    end if
  end subroutine Gen1IntAPIShellView

  !> \brief updates the information of molecule
  !> \author Bin Gao
  !> \date 2013-05-16
  !> \param num_atoms is the number of atoms to update
  !> \param idx_atoms contains the indices of atoms to update
  !> \param charge_atoms contains the charges of atoms
  !> \param coord_atoms contains the coordinates of atoms
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param gauge_origin contains the coordinates of gauge origin of the magnetic vector potential
  !> \param origin_LPF contains the coordinates of origin of the London phase factor
  subroutine Gen1IntAPIUpdateMolecule(num_atoms, idx_atoms, charge_atoms, coord_atoms, &
                                      dipole_origin, gauge_origin, origin_LPF)
    integer, intent(in) :: num_atoms
    integer, optional, intent(in) :: idx_atoms(num_atoms)
    real(REALK), optional, intent(in) :: charge_atoms(num_atoms)
    real(REALK), optional, intent(in) :: coord_atoms(3,num_atoms)
    real(REALK), optional, intent(in) :: dipole_origin(3)
    real(REALK), optional, intent(in) :: gauge_origin(3)
    real(REALK), optional, intent(in) :: origin_LPF(3)
    integer iatom  !incremental recorder over atoms
    if (.not.api_inited) stop "Gen1IntAPIUpdateMolecule>> interface is not initialized!"
    if (present(idx_atoms)) then
      if (present(charge_atoms)) then
        do iatom = 1, num_atoms
          if (idx_atoms(iatom)>=1 .and. idx_atoms(iatom)<=api_num_atoms) then
            api_charge_atoms(idx_atoms(iatom)) = charge_atoms(iatom)
          else
            call quit("Gen1IntAPIUpdateMolecule>> wrong index when updating charges!")
          end if
        end do
      end if
      if (present(coord_atoms)) then
        do iatom = 1, num_atoms
          if (idx_atoms(iatom)>=1 .and. idx_atoms(iatom)<=api_num_atoms) then
            api_coord_atoms(:,idx_atoms(iatom)) = coord_atoms(:,iatom)
          else
            call quit("Gen1IntAPIUpdateMolecule>> wrong index when updating coordinates!")
          end if
        end do
      end if
    ! the first \var(num_atoms) atoms will update
    else
      if (num_atoms>api_num_atoms) then
        call quit("Gen1IntAPIUpdateMolecule>> too many atoms to update!")
      end if
      if (present(charge_atoms)) then
        do iatom = 1, num_atoms
          api_charge_atoms(iatom) = charge_atoms(iatom)
        end do
      end if
      if (present(coord_atoms)) then
        do iatom = 1, num_atoms
          api_coord_atoms(:,iatom) = coord_atoms(:,iatom)
        end do
      end if
    end if
    if (present(dipole_origin)) api_dipole_origin = dipole_origin
    if (present(gauge_origin)) api_gauge_origin = gauge_origin
    if (present(origin_LPF)) api_origin_LPF = origin_LPF
  end subroutine Gen1IntAPIUpdateMolecule

  !> \brief gets the number of atomic orbitals in host programs
  !> \author Bin Gao and Radovan Bast
  !> \date 2012-03-09
  !> \return num_ao is the number of atomic orbitals
  subroutine Gen1IntAPIGetNumAO(num_ao)
    integer, intent(out) :: num_ao
    integer icomp      !incremental recorder over components
    integer idx_first  !index of the first orbital in the last AO sub-shell
    integer idx_last   !index of the last orbital in the last AO sub-shell
    num_ao = 0
    if (api_inited) then
      do icomp = 1, NUM_COMPONENTS
        if (num_sub_shells(icomp)>0) then
          call Gen1IntShellGetRangeAO(sub_shell=sub_shells(num_sub_shells(icomp),icomp), &
                                      idx_first=idx_first, idx_last=idx_last)
          num_ao = idx_last
        end if
      end do
    else
      call quit("Gen1IntAPIGetNumAO>> sub-shells are not created!")
    end if
  end subroutine Gen1IntAPIGetNumAO

  !> \brief calculates molecular orbitals at grid points
  !> \author Bin Gao and Radovan Bast
  !> \date 2012-03-11
  !> \param comp_shell contains the components
  !> \param mo_coef contains the molecular orbital coefficients
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_derv is the number of derivatives
  !> \param num_mo is the number of molecular orbitals
  !> \param api_comm is the MPI communicator
  !> \param gto_type specifies the type of GTOs, should be either NON_LAO (non London atomic
  !>        orbital), LONDON (London atomic orbital, LAO), or ROT_LAO (rotational LAO), only
  !>        NON_LAO implemented
  !> \param order_mag is the order of magnetic derivatives
  !> \param order_ram is the order of derivatives w.r.t. total rotational angular momentum
  !> \param order_geo is the order of geometric derivatives
  !> \return val_mo contains the value of molecular orbitals at grid points
  subroutine Gen1IntAPIGetMO(comp_shell, mo_coef, num_points, grid_points, &
                             num_derv, num_mo, val_mo, api_comm, gto_type, &
                             order_mag, order_ram, order_geo)
    integer, intent(in) :: comp_shell(:)
    type(matrix), intent(in) :: mo_coef
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_derv
    integer, intent(in) :: num_mo
    real(REALK), intent(out) :: val_mo(num_points*num_derv,num_mo)
    integer, optional, intent(in) :: api_comm
    integer, optional, intent(in) :: gto_type
    integer, optional, intent(in) :: order_mag
    integer, optional, intent(in) :: order_ram
    integer, optional, intent(in) :: order_geo
    integer icomp, jcomp  !incremental recorders over components
    ! initializes
    val_mo = 0.0_REALK
    ! loops over components
    do icomp = 1, size(comp_shell)
      jcomp = comp_shell(icomp)
      if (jcomp>0 .and. jcomp<=NUM_COMPONENTS) then
        if (num_sub_shells(jcomp)>0) then
          call Gen1IntShellGetMO(num_shells=num_sub_shells(jcomp), &
                                 sub_shells=sub_shells(:,jcomp),   &
                                 mo_coef=mo_coef,                  &
                                 num_points=num_points,            &
                                 grid_points=grid_points,          &
                                 num_derv=num_derv,                &
                                 num_mo=num_mo,                    &
                                 val_mo=val_mo,                    &
                                 api_comm=api_comm,                &
                                 gto_type=gto_type,                &
                                 order_mag=order_mag,              &
                                 order_ram=order_ram,              &
                                 order_geo=order_geo)
        end if
      else
        call quit("Gen1IntAPIGetMO>> invalid component!")
      end if
    end do
  end subroutine Gen1IntAPIGetMO

  !> \brief terminates Gen1Int interface after all calculations
  !> \author Bin Gao and Radovan Bast
  !> \date 2011-10-02
  subroutine Gen1IntAPIDestroy
    integer icomp  !incremental recorder over components
    if (api_inited) then
      do icomp = 1, NUM_COMPONENTS
        call Gen1IntShellDestroy(num_shells=num_sub_shells(icomp), &
                                 sub_shells=sub_shells(:,icomp))
      end do
      deallocate(sub_shells)
      num_sub_shells = 0
      api_num_atoms = 0              !number of atoms
      deallocate(api_coord_atoms)    !coordinates of atoms
      deallocate(api_charge_atoms)   !charges of atoms
      api_dipole_origin = 0.0_REALK  !coordinates of dipole origin
      api_gauge_origin = 0.0_REALK   !coordinates of gauge origin
      api_origin_LPF = 0.0_REALK     !coordinates of origin of the London phase factor
    end if
    api_inited = .false.
  end subroutine Gen1IntAPIDestroy

  !> \brief creates the operator of property integrals with non-zero components
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param gto_type specifies the type of GTOs, should be either NON_LAO (non London atomic
  !>        orbital), LONDON (London atomic orbital, LAO), or ROT_LAO (rotational LAO), only
  !>        NON_LAO implemented
  !> \param prop_name is the name of property integrals
  !> \param order_mom is the order of multipole moments
  !> \param order_mag_bra is the order of partial magnetic derivatives on bra center, not implemented
  !> \param order_mag_ket is the order of partial magnetic derivatives on ket center, not implemented
  !> \param order_mag_total is the order of total magnetic derivatives, not implemented
  !> \param order_ram_bra is the order of partial derivatives w.r.t. the total rotational
  !>        angular momentum on bra center, not implemented
  !> \param order_ram_ket is the order of partial derivatives w.r.t. the total rotational
  !>        angular momentum on ket center, not implemented
  !> \param order_ram_total is the order of total derivatives w.r.t. the total rotational
  !>        angular momentum, not implemented
  !> \param add_sr is for scalar-relativistic (SR) correction, not implemented
  !> \param add_so is for spin-orbit (SO) correction, not implemented
  !> \param add_london transforms the operator by the LAO type gauge-including projector, not implemented
  !> \return prop_comp is the operator of property integrals
  subroutine Gen1IntAPIPropCreate(gto_type, prop_name, order_mom, &
                                  order_elec,                     &
                                  order_mag_bra, order_mag_ket,   &
                                  order_mag_total,                &
                                  order_ram_bra, order_ram_ket,   &
                                  order_ram_total,                &
                                  add_sr, add_so, add_london,     &
                                  nr_active_blocks,               &
                                  active_component_pairs,         &
                                  prop_comp)

    use london_ao
    integer,           intent(in)    :: gto_type
    character*(*),     intent(in)    :: prop_name
    integer,           intent(in)    :: order_mom
    integer,           intent(in)    :: order_elec
    integer,           intent(in)    :: order_mag_bra
    integer,           intent(in)    :: order_mag_ket
    integer,           intent(in)    :: order_mag_total
    integer,           intent(in)    :: order_ram_bra
    integer,           intent(in)    :: order_ram_ket
    integer,           intent(in)    :: order_ram_total
    logical,           intent(in)    :: add_sr
    logical,           intent(in)    :: add_so
    logical,           intent(in)    :: add_london
    integer,           intent(in)    :: nr_active_blocks
    integer,           intent(in)    :: active_component_pairs(*)
    type(prop_comp_t), intent(inout) :: prop_comp

    integer ierr  !error information

    ! in Dalton active_component_pairs is (/1, 1/)
    ! in DIRAC active_component_pairs can be (/1, 1, 2, 2/)
    !                                        (/1, 2, 2, 1/)
    !                                        (/1, 1/)
    !                                        (/2, 2/)
    !                                        (/1, 2/)
    !                                        (/2, 1/)
    allocate(prop_comp%nnz_comp(2, nr_active_blocks), stat=ierr)
    if (ierr /= 0) then
       call quit("Gen1IntAPIPropCreate>> failed to allocate nnz_comp!")
    end if
    prop_comp%nnz_comp = reshape(active_component_pairs(1:nr_active_blocks*2), (/2, nr_active_blocks/))

    select case (trim(prop_name))
    ! one-electron Hamiltonian
    case (INT_ONE_HAMIL)
      call OnePropCreate(prop_name=INT_ONE_HAMIL,      &
                         one_prop=prop_comp%one_prop,  &
                         info_prop=ierr,               &
                         coord_nuclei=api_coord_atoms, &
                         charge_nuclei=api_charge_atoms)
    ! Cartesian multipole moments
    case (INT_CART_MULTIPOLE)
      call OnePropCreate(prop_name=INT_CART_MULTIPOLE,    &
                         one_prop=prop_comp%one_prop,     &
                         info_prop=ierr,                  &
                         dipole_origin=api_dipole_origin, &
                         order_mom=order_mom,             &
                         order_elec=order_elec)
    ! overlap integrals
    case (INT_OVERLAP)
      call OnePropCreate(prop_name=INT_OVERLAP,       &
                         one_prop=prop_comp%one_prop, &
                         info_prop=ierr)
    ! kinetic energy integrals
    case (INT_KIN_ENERGY)
      call OnePropCreate(prop_name=INT_KIN_ENERGY,    &
                         one_prop=prop_comp%one_prop, &
                         info_prop=ierr)
    ! one-electron potential energy integrals
    case (INT_POT_ENERGY)
      call OnePropCreate(prop_name=INT_POT_ENERGY,     &
                         one_prop=prop_comp%one_prop,  &
                         info_prop=ierr,               &
                         coord_nuclei=api_coord_atoms, &
                         charge_nuclei=api_charge_atoms)
    ! angular momentum integrals
    case (INT_ANGMOM)
      call OnePropCreate(prop_name=INT_ANGMOM,        &
                         one_prop=prop_comp%one_prop, &
                         info_prop=ierr,              &
                         dipole_origin=api_dipole_origin)
    case default
      write(STDOUT,999) "unknown property "//trim(prop_name)//"!"
      call quit('unknown property')
    end select
    if (ierr/=0) then
      write(STDOUT,999) "failed to create operator of "//trim(prop_name)//"!"
      call quit('failed to create propert operator')
    end if
    ! sets magnetic derivatives
    call OnePropSetMag(one_prop=prop_comp%one_prop, &
                       order_mag=order_mag_total,   &
                       order_mag_bra=order_mag_bra, &
                       order_mag_ket=order_mag_ket)
    ! sets derivatives w.r.t. total rotational angular momentum
    call OnePropSetRAM(one_prop=prop_comp%one_prop, &
                       order_ram=order_ram_total,   &
                       order_ram_bra=order_ram_bra, &
                       order_ram_ket=order_ram_ket)
    ! sets the information of London atomic orbitals
    if (gto_type/=NON_LAO .and.                                           &
        (order_mag_total>0 .or. order_mag_bra>0 .or. order_mag_ket>0 .or. &
         order_ram_total>0 .or. order_ram_bra>0 .or. order_ram_ket>0)) then
      call OnePropSetLAO(one_prop=prop_comp%one_prop,   &
                         gauge_origin=api_gauge_origin, &
                         origin_London_PF=api_origin_LPF)
    end if
    if (ierr/=0) then
      call quit("Gen1IntAPIPropCreate>> invalid type of GTOs!")
    end if
999 format("Gen1IntAPIPropCreate>> ",A)
  end subroutine Gen1IntAPIPropCreate

  !> \brief visualizes the operator of property integrals with non-zero components
  !> \author Bin Gao
  !> \date 2012-05-14
  !> \param prop_comp is the operator of property integrals
  !> \param io_viewer is the logical unit number of the viewer
  subroutine Gen1IntAPIPropView(prop_comp, io_viewer)
    type(prop_comp_t), intent(in) :: prop_comp
    integer, intent(in) :: io_viewer
    integer icomp  !incremental recorder over components
    call OnePropView(one_prop=prop_comp%one_prop, io_viewer=io_viewer)
    do icomp = 1, size(prop_comp%nnz_comp,2)
      write(io_viewer,100) "non-zero components", prop_comp%nnz_comp(:,icomp)
    end do
100 format("Gen1IntAPIPropView>> ",A,2I4)
  end subroutine Gen1IntAPIPropView

  !> \brief returns the number of property integral matrices for given one-electron property integrals
  !> \author Bin Gao
  !> \date 2012-05-15
  !> \param prop_comp is the operator of property integrals
  !> \return num_prop is the number of property integral matrices
  subroutine Gen1IntAPIPropGetNumProp(prop_comp, num_prop)
    type(prop_comp_t), intent(in) :: prop_comp
    integer, intent(out) :: num_prop
    call OnePropGetNumProp(one_prop=prop_comp%one_prop, num_prop=num_prop)
  end subroutine Gen1IntAPIPropGetNumProp

  !> \brief returns the symmetry of property integral matrices for given one-electron property integrals
  !> \author Bin Gao
  !> \date 2012-01-12
  !> \param prop_comp is the operator of property integrals
  !> \return prop_sym indicates the symmetry of property integral matrices (SYMM_INT_MAT,
  !>         ANTI_INT_MAT, or SQUARE_INT_MAT)
  subroutine Gen1IntAPIPropGetSymmetry(prop_comp, prop_sym)
    type(prop_comp_t), intent(in) :: prop_comp
    integer, intent(out) :: prop_sym
    call OnePropGetSymmetry(one_prop=prop_comp%one_prop, prop_sym=prop_sym)
  end subroutine Gen1IntAPIPropGetSymmetry

  !> \brief evaluates the integral matrices and/or expectation values
  !> \author Bin Gao and Radovan Bast
  !> \date 2011-01-11
  !> \param prop_comp contains the information of one-electron property integrals and non-zero components
  !> \param nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \param nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \param nary_tree_total is the N-ary tree for total geometric derivatives
  !> \param api_comm is the MPI communicator
  !> \param num_ints is the number of property integral matrices including various derivatives
  !> \param write_ints indicates if writing integral matrices on file
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param write_expt indicates if writing expectation values on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note the arrangement of var(val_ints) and \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntAPIPropGetIntExpt(prop_comp, nary_tree_bra, nary_tree_ket, &
                                      nary_tree_total, api_comm,               &
                                      num_ints, val_ints, write_ints,          &
                                      num_dens, ao_dens, val_expt, write_expt, &
                                      io_viewer, level_print)
    type(prop_comp_t), intent(in) :: prop_comp
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer, optional, intent(in) :: api_comm
    integer, intent(in) :: num_ints
    type(matrix), optional, intent(inout) :: val_ints(num_ints)
    logical, optional, intent(in) :: write_ints
    integer, intent(in) :: num_dens
    type(matrix), optional, intent(in) :: ao_dens(num_dens)
    real(REALK), optional, intent(inout) :: val_expt(num_ints*num_dens)
    logical, optional, intent(in) :: write_expt
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    ! checks if the AO sub-shells are created
    if (.not.api_inited) then
      call quit("Gen1IntAPIPropGetIntExpt>> sub-shells are not created!")
    end if
    call Gen1IntOnePropGetIntExpt(nnz_comp=prop_comp%nnz_comp,     &
                                  one_prop=prop_comp%one_prop,     &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
                                  api_comm=api_comm,               &
                                  num_ints=num_ints,               &
                                  val_ints=val_ints,               &
                                  write_ints=write_ints,           &
                                  num_dens=num_dens,               &
                                  ao_dens=ao_dens,                 &
                                  val_expt=val_expt,               &
                                  write_expt=write_expt,           &
                                  io_viewer=io_viewer,             &
                                  level_print=level_print)
  end subroutine Gen1IntAPIPropGetIntExpt

  !> \brief evaluates the one-electron property intergrands contracted with AO density matrices
  !> \author Bin Gao
  !> \date 2012-05-15
  !> \param prop_comp contains the information of one-electron property integrals and non-zero components
  !> \param nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \param nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \param nary_tree_total is the N-ary tree for total geometric derivatives
  !> \param api_comm is the MPI communicator
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param num_ints is the number of property integral matrices including various derivatives
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_expt contains the one-electron property intergrands contracted with AO density matrices
  !> \note the arrangement of \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntAPIPropGetFunExpt(prop_comp, nary_tree_bra, nary_tree_ket, &
                                      nary_tree_total, api_comm,               &
                                      num_points, grid_points,                 &
                                      num_dens, ao_dens, num_ints, val_expt,   &
                                      io_viewer, level_print)
    type(prop_comp_t), intent(in) :: prop_comp
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer, optional, intent(in) :: api_comm
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_dens
    type(matrix), intent(in) :: ao_dens(num_dens)
    integer, intent(in) :: num_ints
    real(REALK), intent(inout) :: val_expt(num_points,num_ints,num_dens)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer num_nnz_comp         !number of non-zero components
    integer num_prop             !number of property integrals
    integer num_geo_bra          !number of geometric derivatives on bra center
    integer num_geo_ket          !number of geometric derivatives on ket center
    integer num_geo_total        !number of total geometric derivatives
    integer idx_path_bra         !index of current path of N-ary tree for geometric derivatives on bra center
    integer num_paths_bra        !total number of different paths of N-ary tree for geometric derivatives on bra center
    integer idx_path_ket         !index of current path of N-ary tree for geometric derivatives on ket center
    integer num_paths_ket        !total number of different paths of N-ary tree for geometric derivatives on ket center
    integer idx_path_total       !index of current path of N-ary tree for total geometric derivatives
    integer num_paths_total      !total number of different paths of N-ary tree for total geometric derivatives
    integer num_cent_bra         !number of differentiated centers for geometric derivatives on bra center
    integer num_cent_ket         !number of differentiated centers for geometric derivatives on ket center
    integer ipath, jpath, kpath  !incremental recorder over different paths
    integer icomp                !incremental recorder over components
    integer comp_bra             !which component of AO sub-shells on bra center
    integer comp_ket             !which component of AO sub-shells on ket center
    logical same_braket          !if the AO sub-shells are the same on bra and ket centers
    ! checks if the AO sub-shells are created
    if (.not.api_inited) then
      call quit("Gen1IntAPIPropGetFunExpt>> sub-shells are not created!")
    end if
    ! gets the number of non-zero components
    num_nnz_comp = size(prop_comp%nnz_comp,2)
    ! dumps AO sub-shells
    if (level_print>=15) then
      do icomp = 1, num_nnz_comp
        comp_bra = prop_comp%nnz_comp(1,icomp)
        comp_ket = prop_comp%nnz_comp(2,icomp)
        same_braket = comp_bra==comp_ket
        write(io_viewer,100) "AO sub-shells on bra center", comp_bra
        call Gen1IntShellView(num_shells=num_sub_shells(comp_bra), &
                              sub_shells=sub_shells(:,comp_bra),   &
                              io_viewer=io_viewer)
        if (.not.same_braket) then
          write(io_viewer,100) "AO sub-shells on ket center", comp_ket
          call Gen1IntShellView(num_shells=num_sub_shells(comp_ket), &
                                sub_shells=sub_shells(:,comp_ket),   &
                                io_viewer=io_viewer)
        end if
      end do
    end if
    ! dumps other information to check
    if (level_print>=10) then
      ! the information of one-electron property integrals
      call OnePropView(one_prop=prop_comp%one_prop, io_viewer=io_viewer)
      write(io_viewer,100) "information of partial geometric derivatives on bra center"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
      write(io_viewer,100) "information of partial geometric derivatives on ket center"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
      write(io_viewer,100) "information of total geometric derivatives"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
      ! MPI communicator
      if (present(api_comm)) write(io_viewer,100) "MPI communicator provided"
      ! arguments related to AO density matrices
      write(io_viewer,100) "number of AO density matrices", num_dens
      if (level_print>=20) then
        do icomp = 1, num_dens
          write(io_viewer,"()")
          write(io_viewer,100) "AO density matrix", icomp
          call MatView(A=ao_dens(icomp), io_viewer=io_viewer)
          write(io_viewer,"()")
        end do
      end if
    end if
    ! gets the number of property integrals
    call OnePropGetNumProp(one_prop=prop_comp%one_prop, num_prop=num_prop)
    ! gets the numbers of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)
    ! checks the size of input arguments
    if (num_prop*num_geo_bra*num_geo_ket*num_geo_total>num_ints) then
      write(io_viewer,100) "input size", num_ints
      write(io_viewer,100) "required size", num_prop*num_geo_bra*num_geo_ket*num_geo_total
      call quit("Gen1IntOnePropGetIntExpt>> input array not enough!")
    end if
    ! gets the index of current path and total number of different paths
    call NaryTreePathGetIndex(nary_tree=nary_tree_bra, idx_path=idx_path_bra)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_bra, num_paths=num_paths_bra)
    call NaryTreePathGetIndex(nary_tree=nary_tree_ket, idx_path=idx_path_ket)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_ket, num_paths=num_paths_ket)
    call NaryTreePathGetIndex(nary_tree=nary_tree_total, idx_path=idx_path_total)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_total, num_paths=num_paths_total)
    ! loops over different paths
    do kpath = idx_path_total, num_paths_total
      ! dumps the information of current path
      if (level_print>=20) then
        call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
      end if
      do jpath = idx_path_ket, num_paths_ket
        ! dumps the information of current path
        if (level_print>=20) then
          call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
        end if
        ! gets the number of differentiated centers on ket center
        call NaryTreePathGetNumCenters(nary_tree=nary_tree_ket, num_centers=num_cent_ket)
        if (num_cent_ket<=1) then
          do ipath = idx_path_bra, num_paths_bra
            ! dumps the information of current path
            if (level_print>=20) then
              call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
            end if
            ! gets the number of differentiated centers on bra center
            call NaryTreePathGetNumCenters(nary_tree=nary_tree_bra, num_centers=num_cent_bra)
            if (num_cent_bra<=1) then
              ! calculates the property intergrands of current path
              do icomp = 1, num_nnz_comp
                comp_bra = prop_comp%nnz_comp(1,icomp)
                comp_ket = prop_comp%nnz_comp(2,icomp)
                same_braket = comp_bra==comp_ket
                call Gen1IntShellGetFunExpt(num_shells_bra=num_sub_shells(comp_bra), &
                                            sub_shells_bra=sub_shells(:,comp_bra),   &
                                            num_shells_ket=num_sub_shells(comp_ket), &
                                            sub_shells_ket=sub_shells(:,comp_ket),   &
                                            same_braket=same_braket,                 &
                                            one_prop=prop_comp%one_prop,             &
                                            nary_tree_bra=nary_tree_bra,             &
                                            nary_tree_ket=nary_tree_ket,             &
                                            nary_tree_total=nary_tree_total,         &
                                            api_comm=api_comm,                       &
                                            num_points=num_points,                   &
                                            grid_points=grid_points,                 &
                                            num_dens=num_dens,                       &
                                            ao_dens=ao_dens,                         &
                                            num_prop=num_prop,                       &
                                            num_geo_bra=num_geo_bra,                 &
                                            num_geo_ket=num_geo_ket,                 &
                                            num_geo_total=num_geo_total,             &
                                            val_expt=val_expt)
              end do
            end if
            ! generates the differentiated centers and their orders of the next path
            call NaryTreeSearch(nary_tree=nary_tree_bra)
          end do
        end if
        ! generates the differentiated centers and their orders of the next path
        call NaryTreeSearch(nary_tree=nary_tree_ket)
      end do
      ! generates the differentiated centers and their orders of the next path
      call NaryTreeSearch(nary_tree=nary_tree_total)
    end do
100 format("Gen1IntAPIPropGetFunExpt>> ",A,I8)
  end subroutine Gen1IntAPIPropGetFunExpt

  !> \brief frees space taken by the operator of property integrals with non-zero components
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param prop_comp is the operator of property integrals
  subroutine Gen1IntAPIPropDestroy(prop_comp)
    type(prop_comp_t), intent(inout) :: prop_comp
    call OnePropDestroy(one_prop=prop_comp%one_prop)
    deallocate(prop_comp%nnz_comp)
  end subroutine Gen1IntAPIPropDestroy

  !> \brief creates N-ary tree for geometric derivatives
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param max_num_cent is the maximum number of differentiated centers for total
  !>        geometric derivatives
  !> \param order_geo is the order of geometric derivatives
  !> \param num_geo_atoms is the number of selected atoms which might be chosen as the
  !>        differentiated centers, <1 means using all atoms, not tested
  !> \param idx_geo_atoms contains the indices of the selected atoms, will not be used
  !>        if \var(num_geo_atoms) is <1, not tested
  !> \return nary_tree is the N-ary tree for total geometric derivatives
  subroutine Gen1IntAPINaryTreeCreate(max_num_cent, order_geo, num_geo_atoms, &
                                      idx_geo_atoms, nary_tree)
    integer, intent(in) :: max_num_cent
    integer, intent(in) :: order_geo
    integer, intent(in) :: num_geo_atoms
    integer, intent(in) :: idx_geo_atoms(*)
    type(nary_tree_t), intent(inout) :: nary_tree
    integer ierr  !error information
    if (order_geo>=0) then
      if (num_geo_atoms>0) then
        call NaryTreeCreate(num_atoms=num_geo_atoms,   &
                            order_geo=order_geo,       &
                            max_num_cent=max_num_cent, &
                            nary_tree=nary_tree,       &
                            info_geom=ierr)
        if (ierr/=0) then
          call quit("Gen1IntAPINaryTreeCreate>> error occurred when calling NaryTreeCreate!")
        end if
        call NaryTreeSetAtoms(num_atoms=num_geo_atoms, &
                              idx_atoms=idx_geo_atoms, &
                              nary_tree=nary_tree,     &
                              info_geom=ierr)
        if (ierr/=0) then
          call quit("Gen1IntAPINaryTreeCreate>> error occurred when calling NaryTreeSetAtoms!")
        end if
      else
        call NaryTreeCreate(num_atoms=api_num_atoms,   &
                            order_geo=order_geo,       &
                            max_num_cent=max_num_cent, &
                            nary_tree=nary_tree,       &
                            info_geom=ierr)
        if (ierr/=0) then
          call quit("Gen1IntAPINaryTreeCreate>> error occurred when calling NaryTreeCreate!")
        end if
      end if
    else
      call quit("Gen1IntAPINaryTreeCreate>> negative order of geometric derivatives!")
    end if
  end subroutine Gen1IntAPINaryTreeCreate

  !> \brief evaluates the integral matrices and/or expectation values
  !> \author Bin Gao and Radovan Bast
  !> \date 2011-01-11
  !> \param nnz_comp contains the non-zero components
  !> \param one_prop contains the information of one-electron property integrals
  !> \param nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \param nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \param nary_tree_total is the N-ary tree for total geometric derivatives
  !> \param api_comm is the MPI communicator
  !> \param num_ints is the number of property integral matrices including various derivatives
  !> \param write_ints indicates if writing integral matrices on file
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param write_expt indicates if writing expectation values on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note the arrangement of var(val_ints) and \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntOnePropGetIntExpt(nnz_comp, one_prop,                            &
                                      nary_tree_bra, nary_tree_ket, nary_tree_total, &
                                      api_comm, num_ints, val_ints, write_ints,      &
                                      num_dens, ao_dens, val_expt, write_expt,       &
                                      io_viewer, level_print)
    integer, intent(in) :: nnz_comp(:,:)
    type(one_prop_t), intent(in) :: one_prop
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer, optional, intent(in) :: api_comm
    integer, intent(in) :: num_ints
    type(matrix), optional, intent(inout) :: val_ints(num_ints)
    logical, optional, intent(in) :: write_ints
    integer, intent(in) :: num_dens
    type(matrix), optional, intent(in) :: ao_dens(num_dens)
    real(REALK), optional, intent(inout) :: val_expt(num_ints*num_dens)
    logical, optional, intent(in) :: write_expt
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer num_nnz_comp         !number of non-zero components
    integer num_prop             !number of property integrals
    integer num_geo_bra          !number of geometric derivatives on bra center
    integer num_geo_ket          !number of geometric derivatives on ket center
    integer num_geo_total        !number of total geometric derivatives
    integer idx_path_bra         !index of current path of N-ary tree for geometric derivatives on bra center
    integer num_paths_bra        !total number of different paths of N-ary tree for geometric derivatives on bra center
    integer idx_path_ket         !index of current path of N-ary tree for geometric derivatives on ket center
    integer num_paths_ket        !total number of different paths of N-ary tree for geometric derivatives on ket center
    integer idx_path_total       !index of current path of N-ary tree for total geometric derivatives
    integer num_paths_total      !total number of different paths of N-ary tree for total geometric derivatives
    integer num_cent_bra         !number of differentiated centers for geometric derivatives on bra center
    integer num_cent_ket         !number of differentiated centers for geometric derivatives on ket center
    integer ipath, jpath, kpath  !incremental recorder over different paths
    integer icomp                !incremental recorder over components
    integer comp_bra             !which component of AO sub-shells on bra center
    integer comp_ket             !which component of AO sub-shells on ket center
    logical same_braket          !if the AO sub-shells are the same on bra and ket centers
    ! checks if the AO sub-shells are created
    if (.not.api_inited) then
      call quit("Gen1IntOnePropGetIntExpt>> sub-shells are not created!")
    end if
    ! gets the number of non-zero components
    if (size(nnz_comp,1)/=2) then
      call quit("Gen1IntOnePropGetIntExpt>> needs non-zero components on bra and ket centers!")
    end if
    num_nnz_comp = size(nnz_comp,2)
    if (num_nnz_comp>NUM_COMPONENTS) then
      call quit("Gen1IntOnePropGetIntExpt>> too many components!")
    end if
    ! dumps AO sub-shells
    if (level_print>=15) then
      do icomp = 1, num_nnz_comp
        comp_bra = nnz_comp(1,icomp)
        comp_ket = nnz_comp(2,icomp)
        same_braket = comp_bra==comp_ket
        write(io_viewer,100) "AO sub-shells on bra center", comp_bra
        call Gen1IntShellView(num_shells=num_sub_shells(comp_bra), &
                              sub_shells=sub_shells(:,comp_bra),   &
                              io_viewer=io_viewer)
        if (.not.same_braket) then
          write(io_viewer,100) "AO sub-shells on ket center", comp_ket
          call Gen1IntShellView(num_shells=num_sub_shells(comp_ket), &
                                sub_shells=sub_shells(:,comp_ket),   &
                                io_viewer=io_viewer)
        end if
      end do
    end if
    ! dumps other information to check
    if (level_print>=10) then
      ! the information of one-electron property integrals
      call OnePropView(one_prop=one_prop, io_viewer=io_viewer)
      write(io_viewer,100) "information of partial geometric derivatives on bra center"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
      write(io_viewer,100) "information of partial geometric derivatives on ket center"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
      write(io_viewer,100) "information of total geometric derivatives"
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
      ! MPI communicator
      if (present(api_comm)) write(io_viewer,100) "MPI communicator provided"
      ! arguments related to integral matrices
      if (present(val_ints)) write(io_viewer,100) "integral matrices returned"
      if (present(write_ints)) then
        if (write_ints) write(io_viewer,100) "integral matrices will be written on file"
      end if
      ! arguments related to expectation values
      if (present(ao_dens)) then
        write(io_viewer,100) "number of AO density matrices", num_dens
        if (level_print>=20) then
          do icomp = 1, num_dens
            write(io_viewer,"()")
            write(io_viewer,100) "AO density matrix", icomp
            call MatView(A=ao_dens(icomp), io_viewer=io_viewer)
            write(io_viewer,"()")
          end do
        end if
        if (present(val_expt)) write(io_viewer,100) "expectation values returned"
        if (present(write_expt)) then
          if (write_expt) write(io_viewer,100) "expectation values will be written on file"
        end if
      end if
    end if
    ! gets the number of property integrals
    call OnePropGetNumProp(one_prop=one_prop, num_prop=num_prop)
    ! gets the numbers of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)
    ! checks the size of input arguments
    if (present(val_ints) .or. present(val_expt)) then
      if (num_prop*num_geo_bra*num_geo_ket*num_geo_total>num_ints) then
        write(io_viewer,100) "input size", num_ints
        write(io_viewer,100) "required size", num_prop*num_geo_bra*num_geo_ket*num_geo_total
        call quit("Gen1IntOnePropGetIntExpt>> input array not enough!")
      end if
    end if
    ! gets the index of current path and total number of different paths
    call NaryTreePathGetIndex(nary_tree=nary_tree_bra, idx_path=idx_path_bra)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_bra, num_paths=num_paths_bra)
    call NaryTreePathGetIndex(nary_tree=nary_tree_ket, idx_path=idx_path_ket)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_ket, num_paths=num_paths_ket)
    call NaryTreePathGetIndex(nary_tree=nary_tree_total, idx_path=idx_path_total)
    call NaryTreeGetNumPaths(nary_tree=nary_tree_total, num_paths=num_paths_total)
    ! loops over different paths
    do kpath = idx_path_total, num_paths_total
      ! dumps the information of current path
      if (level_print>=20) then
        call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
      end if
      do jpath = idx_path_ket, num_paths_ket
        ! dumps the information of current path
        if (level_print>=20) then
          call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
        end if
        ! gets the number of differentiated centers on ket center
        call NaryTreePathGetNumCenters(nary_tree=nary_tree_ket, num_centers=num_cent_ket)
        if (num_cent_ket<=1) then
          do ipath = idx_path_bra, num_paths_bra
            ! dumps the information of current path
            if (level_print>=20) then
              call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
            end if
            ! gets the number of differentiated centers on bra center
            call NaryTreePathGetNumCenters(nary_tree=nary_tree_bra, num_centers=num_cent_bra)
            if (num_cent_bra<=1) then
              ! calculates the property integrals of current path
              do icomp = 1, num_nnz_comp
                comp_bra = nnz_comp(1,icomp)
                comp_ket = nnz_comp(2,icomp)
                same_braket = comp_bra==comp_ket
                call Gen1IntShellGetIntExpt(num_shells_bra=num_sub_shells(comp_bra), &
                                            sub_shells_bra=sub_shells(:,comp_bra),   &
                                            num_shells_ket=num_sub_shells(comp_ket), &
                                            sub_shells_ket=sub_shells(:,comp_ket),   &
                                            same_braket=same_braket,                 &
                                            one_prop=one_prop,                       &
                                            nary_tree_bra=nary_tree_bra,             &
                                            nary_tree_ket=nary_tree_ket,             &
                                            nary_tree_total=nary_tree_total,         &
                                            api_comm=api_comm,                       &
                                            num_prop=num_prop,                       &
                                            num_geo_bra=num_geo_bra,                 &
                                            num_geo_ket=num_geo_ket,                 &
                                            num_geo_total=num_geo_total,             &
                                            val_ints=val_ints,                       &
                                            write_ints=write_ints,                   &
                                            num_dens=num_dens,                       &
                                            ao_dens=ao_dens,                         &
                                            val_expt=val_expt,                       & 
                                            write_expt=write_expt)
              end do
            end if
            ! generates the differentiated centers and their orders of the next path
            call NaryTreeSearch(nary_tree=nary_tree_bra)
          end do
        end if
        ! generates the differentiated centers and their orders of the next path
        call NaryTreeSearch(nary_tree=nary_tree_ket)
      end do
      ! generates the differentiated centers and their orders of the next path
      call NaryTreeSearch(nary_tree=nary_tree_total)
    end do
100 format("Gen1IntOnePropGetIntExpt>> ",A,I8)
  end subroutine Gen1IntOnePropGetIntExpt

end module gen1int_api
