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
!...  This file contains most host program specific subroutines for Gen1Int interface.
!
!...  2013-05-02, Bin Gao
!...  * fix the bug of returning wrong partial geometric derivatives
!
!...  2012-05-09, Bin Gao
!...  * first version

#include "gen1int_host.h"

  ! public subroutines

  !> \brief initializes Gen1Int interface, for instance, creates the AO sub-shells
  !>        of host program (based on \fn(ORBPRO) subroutine by getting the unnormalized
  !>        contraction coefficients); should be called before any calculation
  !> \author Bin Gao
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
  subroutine gen1int_host_init(num_comp, num_atom_type, KATOM, num_sym_atom, &
                               ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                               num_prim, num_contr, KPRIM, exponents,        &
                               ucontr_coefs)
    use gen1int_api
    implicit none
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
    integer :: ierr  !error information
    logical :: mpi_sync = .true.
#if defined(VAR_MPI)
#include "mpif.h"
#include "iprtyp.h"
#endif
    ! in case of initializing the interface multiple times
    if (Gen1IntAPIInited()) mpi_sync = .false.
    ! initializes API of Gen1Int
    call Gen1IntAPICreate(num_comp, num_atom_type, KATOM, num_sym_atom, &
                          ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
                          num_prim, num_contr, KPRIM, exponents, ucontr_coefs)
#if defined(VAR_MPI)
    if (mpi_sync) then
      ! wakes up workers
      call MPI_Bcast(GEN1INT_INIT, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
      ! broadcasts level of print
      call MPI_Bcast(0, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
      ! broadcasts information of API of Gen1Int
      call Gen1IntAPIBcast(MANAGER, MPI_COMM_WORLD)
    end if
#endif
    return
  end subroutine gen1int_host_init

#if defined(VAR_MPI)
  !> \brief initializes Gen1Int interface, for instance, creates the AO sub-shells
  !>        of host program (based on \fn(ORBPRO) subroutine by getting the unnormalized
  !>        contraction coefficients); should be called before any calculation
  !> \author Bin Gao
  !> \date 2012-05-13
  subroutine gen1int_worker_init
    use gen1int_api
    implicit none
#include "mpif.h"
    ! initializes API of Gen1Int by information from manager processor
    call Gen1IntAPIBcast(MANAGER, MPI_COMM_WORLD)
    return
  end subroutine gen1int_worker_init
#endif

  !> \brief evaluates the one-electron property integral matrices on manager processor
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
  !> \param max_ncent_bra is the maximum number of differentiated centers on bra center
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param num_atoms_bra is the number of selected atoms for differentiated centers
  !>        on bra center, <1 means using all atoms
  !> \param idx_atoms_bra contains the indices of the selected atoms on bra center,
  !>        will not be used if \var(num_atoms_bra) is <1
  !> \param max_ncent_ket is the maximum number of differentiated centers on ket center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param num_atoms_ket is the number of selected atoms for differentiated centers
  !>        on ket center, <1 means using all atoms
  !> \param idx_atoms_ket contains the indices of the selected atoms on ket center,
  !>        will not be used if \var(num_atoms_ket) is <1
  !> \param max_num_cent is the maximum number of differentiated centers for total
  !>        geometric derivatives
  !> \param order_geo_total is the order of total geometric derivatives
  !> \param num_geo_atoms is the number of selected atoms chosen as the differentiated
  !>        centers, <1 means using all atoms
  !> \param idx_geo_atoms contains the indices of the selected atoms, will not be used
  !>        if \var(num_geo_atoms) is <1
  !> \param add_sr is for scalar-relativistic (SR) correction, not implemented
  !> \param add_so is for spin-orbit (SO) correction, not implemented
  !> \param add_london transforms the operator by the LAO type gauge-including projector, not implemented
  !> \param num_ints is the number of property integral matrices including various derivatives
  !> \param write_ints indicates if writing integral matrices on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_ints contains the property integral matrices
  !> \note the arrangement of var(val_ints) will be in the order of \var(order_mom),
  !>       \var(order_mag_bra), ..., \var(order_geo_total), and each of them is arranged
  !>       in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz), see
  !>       Gen1Int library manual, for instance Section 2.2.
  !>       Moreover, the "triangular" total geometric derivatives could be obtained
  !>       from the unique total geometric derivatives by giving \var(max_num_cent)=\var(order_geo_total)
  subroutine gen1int_host_get_int(gto_type, prop_name, order_mom, &
                                  order_elec,                     &
                                  order_mag_bra, order_mag_ket,   &
                                  order_mag_total,                &
                                  order_ram_bra, order_ram_ket,   &
                                  order_ram_total,                &
                                  max_ncent_bra, order_geo_bra,   &
                                  num_atoms_bra, idx_atoms_bra,   &
                                  max_ncent_ket, order_geo_ket,   &
                                  num_atoms_ket, idx_atoms_ket,   &
                                  max_num_cent, order_geo_total,  &
                                  num_geo_atoms, idx_geo_atoms,   &
                                  add_sr, add_so, add_london,     &
                                  num_ints, val_ints, write_ints, &
                                  nr_active_blocks,               &
                                  active_component_pairs,         &
                                  io_viewer, level_print)
    use gen1int_api
    implicit none
    integer,       intent(in)    :: gto_type
    character*(*), intent(in)    :: prop_name
    integer,       intent(in)    :: order_mom
    integer,       intent(in)    :: order_elec
    integer,       intent(in)    :: order_mag_bra
    integer,       intent(in)    :: order_mag_ket
    integer,       intent(in)    :: order_mag_total
    integer,       intent(in)    :: order_ram_bra
    integer,       intent(in)    :: order_ram_ket
    integer,       intent(in)    :: order_ram_total
    integer,       intent(in)    :: max_ncent_bra
    integer,       intent(in)    :: order_geo_bra
    integer,       intent(in)    :: num_atoms_bra
    integer,       intent(in)    :: idx_atoms_bra(*)
    integer,       intent(in)    :: max_ncent_ket
    integer,       intent(in)    :: order_geo_ket
    integer,       intent(in)    :: num_atoms_ket
    integer,       intent(in)    :: idx_atoms_ket(*)
    integer,       intent(in)    :: max_num_cent
    integer,       intent(in)    :: order_geo_total
    integer,       intent(in)    :: num_geo_atoms
    integer,       intent(in)    :: idx_geo_atoms(*)
    logical,       intent(in)    :: add_sr
    logical,       intent(in)    :: add_so
    logical,       intent(in)    :: add_london
    integer,       intent(in)    :: num_ints
    type(matrix),  intent(inout) :: val_ints(*)
    logical,       intent(in)    :: write_ints
    integer,       intent(in)    :: nr_active_blocks
    integer,       intent(in)    :: active_component_pairs(*)
    integer,       intent(in)    :: io_viewer
    integer,       intent(in)    :: level_print

    real(REALK) start_time             !start time
    type(prop_comp_t) prop_comp        !operator of property integrals with non-zero components
    type(nary_tree_t) nary_tree_bra    !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) nary_tree_ket    !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) nary_tree_total  !N-ary tree for total geometric derivatives
    integer ierr                       !error information
#if defined(VAR_MPI)
#include "mpif.h"
#include "iprtyp.h"
#endif
    ! gets the start time
    call xtimer_set(start_time)
#if defined(VAR_MPI)
    ! wakes up workers
    call MPI_Bcast(GEN1INT_GET_INT, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts level of print
    call MPI_Bcast(level_print, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts number of property integral matrices including various derivatives
    call MPI_Bcast(num_ints, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
#endif
    ! creates the operator of property integrals and N-ary tree for total
    ! geometric derivatives on manager processor, and broadcasts other input
    ! arguments to worker processors
    call gen1int_host_prop_create(gto_type, prop_name, order_mom, &
                                  order_elec,                     &
                                  order_mag_bra, order_mag_ket,   &
                                  order_mag_total,                &
                                  order_ram_bra, order_ram_ket,   &
                                  order_ram_total,                &
                                  add_sr, add_so, add_london,     &
                                  io_viewer, level_print,         &
                                  nr_active_blocks,               &
                                  active_component_pairs,         &
                                  prop_comp)
    call gen1int_host_geom_create(max_ncent_bra, order_geo_bra,   &
                                  num_atoms_bra, idx_atoms_bra,   &
                                  max_ncent_ket, order_geo_ket,   &
                                  num_atoms_ket, idx_atoms_ket,   &
                                  max_num_cent, order_geo_total,  &
                                  num_geo_atoms, idx_geo_atoms,   &
                                  io_viewer, level_print,         &
                                  nary_tree_bra, nary_tree_ket, nary_tree_total)
    ! performs calculations
    call Gen1IntAPIPropGetIntExpt(prop_comp=prop_comp,             &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
#if defined(VAR_MPI)
                                  api_comm=MPI_COMM_WORLD,         &
#endif
                                  num_ints=num_ints,               &
                                  val_ints=val_ints,               &
                                  write_ints=write_ints,           &
                                  num_dens=0,                      &
                                  io_viewer=io_viewer,             &
                                  level_print=level_print)
    ! free spaces
    call Gen1IntAPIPropDestroy(prop_comp=prop_comp)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)
#if defined(VAR_MPI)
    ! blocks until all processors have finished
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    ! prints the CPU elapsed time
    call xtimer_view(start_time, trim(prop_name)//"@gen1int_host_get_int", io_viewer)
    return
  end subroutine gen1int_host_get_int

#if defined(VAR_MPI)
  !> \brief evaluates the one-electron property integral matrices on worker processors
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  subroutine gen1int_worker_get_int(len_work, wrk_space, io_viewer, level_print)
    use gen1int_api
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer order_geo_bra              !order of geometric derivatives on bra center
    integer order_geo_ket              !order of geometric derivatives on ket center
    integer order_geo_total            !order of total geometric derivatives
    integer num_ints                   !number of property integral matrices including various derivatives
    type(prop_comp_t) prop_comp        !operator of property integrals with non-zero components
    type(nary_tree_t) nary_tree_bra    !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) nary_tree_ket    !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) nary_tree_total  !N-ary tree for total geometric derivatives
    integer ierr                       !error information
#include "mpif.h"
    ! gets the number of property integral matrices including various derivatives
    call MPI_Bcast(num_ints, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! creates the operator of property integrals and N-ary trees for geometric
    ! derivatives on worker processors by the arguments from manager
    call gen1int_worker_prop_create(io_viewer, level_print, prop_comp)
    call gen1int_worker_geom_create(io_viewer, level_print, nary_tree_bra, nary_tree_ket, nary_tree_total)
    ! performs calculations
    call Gen1IntAPIPropGetIntExpt(prop_comp=prop_comp,             &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
                                  api_comm=MPI_COMM_WORLD,         &
                                  num_ints=num_ints,               &
                                  num_dens=0,                      &
                                  io_viewer=io_viewer,             &
                                  level_print=level_print)
    ! free spaces
    call Gen1IntAPIPropDestroy(prop_comp=prop_comp)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)
    ! blocks until all processors have finished
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    return
  end subroutine gen1int_worker_get_int
#endif

  !> \brief evaluates the one-electron property expectation values on manager processor
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param write_expt indicates if writing expectation values on file
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_expt contains the expectation values
  !> \note please see the comments of \fn(gen1int_host_get_int) of other arguments,
  !>       \var(val_expt) should be zero by users before calculations
  subroutine gen1int_host_get_expt(gto_type, prop_name, order_mom, &
                                   order_elec,                     &
                                   order_mag_bra, order_mag_ket,   &
                                   order_mag_total,                &
                                   order_ram_bra, order_ram_ket,   &
                                   order_ram_total,                &
                                   max_ncent_bra, order_geo_bra,   &
                                   num_atoms_bra, idx_atoms_bra,   &
                                   max_ncent_ket, order_geo_ket,   &
                                   num_atoms_ket, idx_atoms_ket,   &
                                   max_num_cent, order_geo_total,  &
                                   num_geo_atoms, idx_geo_atoms,   &
                                   add_sr, add_so, add_london,     &
                                   num_dens, ao_dens,              &
                                   num_ints, val_expt, write_expt, &
                                   nr_active_blocks,               &
                                   active_component_pairs,         &
                                   io_viewer, level_print)
    use gen1int_api
    implicit none
    integer,       intent(in)    :: gto_type
    character*(*), intent(in)    :: prop_name
    integer,       intent(in)    :: order_mom
    integer,       intent(in)    :: order_elec
    integer,       intent(in)    :: order_mag_bra
    integer,       intent(in)    :: order_mag_ket
    integer,       intent(in)    :: order_mag_total
    integer,       intent(in)    :: order_ram_bra
    integer,       intent(in)    :: order_ram_ket
    integer,       intent(in)    :: order_ram_total
    integer,       intent(in)    :: max_ncent_bra
    integer,       intent(in)    :: order_geo_bra
    integer,       intent(in)    :: num_atoms_bra
    integer,       intent(in)    :: idx_atoms_bra(*)
    integer,       intent(in)    :: max_ncent_ket
    integer,       intent(in)    :: order_geo_ket
    integer,       intent(in)    :: num_atoms_ket
    integer,       intent(in)    :: idx_atoms_ket(*)
    integer,       intent(in)    :: max_num_cent
    integer,       intent(in)    :: order_geo_total
    integer,       intent(in)    :: num_geo_atoms
    integer,       intent(in)    :: idx_geo_atoms(*)
    logical,       intent(in)    :: add_sr
    logical,       intent(in)    :: add_so
    logical,       intent(in)    :: add_london
    integer,       intent(in)    :: num_dens
    type(matrix),  intent(inout) :: ao_dens(num_dens)
    integer,       intent(in)    :: num_ints
    real(REALK),   intent(inout) :: val_expt(num_ints*num_dens)
    logical,       intent(in)    :: write_expt
    integer,       intent(in)    :: nr_active_blocks
    integer,       intent(in)    :: active_component_pairs(*)
    integer,       intent(in)    :: io_viewer
    integer,       intent(in)    :: level_print

    real(REALK) start_time             !start time
    type(prop_comp_t) prop_comp        !operator of property integrals with non-zero components
    type(nary_tree_t) nary_tree_bra    !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) nary_tree_ket    !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) nary_tree_total  !N-ary tree for total geometric derivatives
    integer idens                      !incremental recorder over AO density matrices
    integer ierr                       !error information
#if defined(VAR_MPI)
#include "mpif.h"
#include "iprtyp.h"
#endif
    ! gets the start time
    call xtimer_set(start_time)
#if defined(VAR_MPI)
    ! wakes up workers
    call MPI_Bcast(GEN1INT_GET_EXPT, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts level of print
    call MPI_Bcast(level_print, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts the number of AO density matrices
    call MPI_Bcast(num_dens, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts AO density matrices
    do idens = 1, num_dens
      call MatBcast(ao_dens(idens), MANAGER, MPI_COMM_WORLD)
    end do
    ! broadcasts the number of property integral matrices including various derivatives
    call MPI_Bcast(num_ints, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
#endif
    ! creates the operator of property integrals and N-ary tree for total
    ! geometric derivatives on manager processor, and broadcasts other input
    ! arguments to worker processors
    call gen1int_host_prop_create(gto_type, prop_name, order_mom, &
                                  order_elec,                     &
                                  order_mag_bra, order_mag_ket,   &
                                  order_mag_total,                &
                                  order_ram_bra, order_ram_ket,   &
                                  order_ram_total,                &
                                  add_sr, add_so, add_london,     &
                                  io_viewer, level_print,         &
                                  nr_active_blocks,               &
                                  active_component_pairs,         &
                                  prop_comp)
    call gen1int_host_geom_create(max_ncent_bra, order_geo_bra,   &
                                  num_atoms_bra, idx_atoms_bra,   &
                                  max_ncent_ket, order_geo_ket,   &
                                  num_atoms_ket, idx_atoms_ket,   &
                                  max_num_cent, order_geo_total,  &
                                  num_geo_atoms, idx_geo_atoms,   &
                                  io_viewer, level_print,         &
                                  nary_tree_bra, nary_tree_ket, nary_tree_total)
    ! performs calculations
    call Gen1IntAPIPropGetIntExpt(prop_comp=prop_comp,             &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
#if defined(VAR_MPI)
                                  api_comm=MPI_COMM_WORLD,         &
#endif
                                  num_dens=num_dens,               &
                                  ao_dens=ao_dens,                 &
                                  num_ints=num_ints,               &
                                  val_expt=val_expt,               &
                                  write_expt=write_expt,           &
                                  io_viewer=io_viewer,             &
                                  level_print=level_print)
    ! free spaces
    call Gen1IntAPIPropDestroy(prop_comp=prop_comp)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)
#if defined(VAR_MPI)
    ! blocks until all processors have finished
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    ! prints the CPU elapsed time
    call xtimer_view(start_time, trim(prop_name)//"@gen1int_host_get_expt", io_viewer)
    return
  end subroutine gen1int_host_get_expt

#if defined(VAR_MPI)
  !> \brief evaluates the one-electron property expectation values on worker processors
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  subroutine gen1int_worker_get_expt(len_work, wrk_space, io_viewer, level_print)
    use gen1int_api
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    integer order_geo_total                  !order of total geometric derivatives
    integer num_ints                         !number of property integral matrices including various derivatives
    integer num_dens                         !number of AO density matrices
    type(matrix), allocatable :: ao_dens(:)  !AO density matrices
    integer idens                            !incremental recorder over AO density matrices
    type(prop_comp_t) prop_comp              !operator of property integrals with non-zero components
    type(nary_tree_t) nary_tree_bra          !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) nary_tree_ket          !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) nary_tree_total        !N-ary tree for total geometric derivatives
    integer ierr                             !error information
#include "mpif.h"
    ! gets the number of AO density matrices
    call MPI_Bcast(num_dens, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! allocates memory for AO density matrices
    allocate(ao_dens(num_dens), stat=ierr)
    if (ierr/=0) then
      call quit("gen1int_worker_get_expt>> failed to allocate ao_dens!")
    end if
    ! gets AO density matrices
    do idens = 1, num_dens
      call MatBcast(ao_dens(idens), MANAGER, MPI_COMM_WORLD)
    end do
    ! gets the number of property integral matrices including various derivatives
    call MPI_Bcast(num_ints, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! creates the operator of property integrals and N-ary trees for geometric
    ! derivatives on worker processors by the arguments from manager
    call gen1int_worker_prop_create(io_viewer, level_print, prop_comp)
    call gen1int_worker_geom_create(io_viewer, level_print, nary_tree_bra, nary_tree_ket, nary_tree_total)
    ! performs calculations
    call Gen1IntAPIPropGetIntExpt(prop_comp=prop_comp,             &
                                  nary_tree_bra=nary_tree_bra,     &
                                  nary_tree_ket=nary_tree_ket,     &
                                  nary_tree_total=nary_tree_total, &
                                  api_comm=MPI_COMM_WORLD,         &
                                  num_dens=num_dens,               &
                                  ao_dens=ao_dens,                 &
                                  num_ints=num_ints,               &
                                  io_viewer=io_viewer,             &
                                  level_print=level_print)
    ! free spaces
    call Gen1IntAPIPropDestroy(prop_comp=prop_comp)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)
    do idens = 1, num_dens
      call MatDestroy(A=ao_dens(idens))
    end do
    deallocate(ao_dens)
    ! blocks until all processors have finished
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    return
  end subroutine gen1int_worker_get_expt
#endif

  !> \brief reads the information of cube file and initializes the information
  !>        for cube files; the following keywords should be added into **WAVE FUNCTIONS
  !>        in DALTON.INP/DIRAC.INP:
  !>        *CUBE
  !>        .DENSITY
  !>        .HOMO
  !>        .LUMO
  !>        .MO
  !>        1,5-9,12
  !>        .FORMAT
  !>        GAUSSIAN
  !>        .ORIGIN
  !>        -2.0  -4.0   -5.0
  !>        .INCREMENT
  !>        100    0.0    0.0    0.1  #of increments in the slowest running direction
  !>         80    0.0    0.1    0.0
  !>         40    0.1    0.0    0.0  #of increments in the fastest running directio
  !>        where .DENSITY, .HOMO, .LUMO, and .MO (followed by the indices of molecular orbtials)
  !>        are not all necessarily needed
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param io_input is the logical unit number of standard input
  !> \param io_viewer is the logical unit number of the viewer
  !> \param word is the keyword from input file
  subroutine gen1int_host_cube_init(io_input, io_viewer, word)
    ! module of generating cube files
    use gen1int_cube
    ! to decode the string of MOs (in Gen1Int library)
    use str_decode
    implicit none
    integer, intent(in) :: io_input
    integer, intent(in) :: io_viewer
    character*(*), intent(inout) :: word
! uses MXCORB
#include "maxorb.h"
! uses DO_CUBE
#include "infinp.h"
    character(MAX_LEN_STR) key_word  !key words read from standard input
    type(decode_str_t) str_idx_mo    !string of indices of MOs
    integer ixyz                     !incremental recorder along XYZ directions
    integer ierr                     !error information
    ! reads in the first line after *CUBE
    read(io_input,"(A)",err=999,end=999) key_word
    do while(key_word(1:1)/="*")
      select case(trim(key_word))
      ! electron density
      case(".DENSITY")
        do_density_cube = .true.
      ! HOMO
      case(".HOMO")
        do_homo_cube = .true.
      ! LUMO
      case(".LUMO")
        do_lumo_cube = .true.
      ! MOs
      case(".MO")
        read(io_input,"(A)",err=999,end=999) key_word
        call StrDecodeCreate(the_str=trim(key_word), &
                             conn_char="-",          &
                             sep_char=",",           &
                             num_ints=num_cube_mo,   &
                             decode_str=str_idx_mo)
        allocate(idx_cube_mo(num_cube_mo), stat=ierr)
        if (ierr/=0) then
          call quit("gen1int_host_cube_init>> failed to allocate idx_cube_mo!")
        end if
        call StrDecodeGetInts(decode_str=str_idx_mo, &
                              num_ints=num_cube_mo,  &
                              convert_ints=idx_cube_mo)
        call StrDecodeDestroy(decode_str=str_idx_mo)
        do_mo_cube = .true.
      ! format of cube file
      case(".FORMAT")
        read(io_input,"(A)",err=999,end=999) cube_format
      ! origin
      case(".ORIGIN")
        read(io_input,*,err=999,end=999) cube_origin
      ! increments
      case(".INCREMENT")
        do ixyz = 1, 3
          ! cube_increment(:,1) is the increment for X
          ! cube_increment(:,2) is the increment for Y
          ! cube_increment(:,3) is the increment for Z
          ! reads N, X, Y, Z
          read(io_input,*,err=999,end=999) cube_num_inc(ixyz), cube_increment(ixyz,:)
        end do
      ! comments or illegal keyword
      case default
        if (key_word(1:1)/="#" .and. key_word(1:1)/="!") then
          write(io_viewer,100) "keyword """//trim(key_word)// &
                               """ is not recognized in *CUBE!"
          call quit('unknown keyword in *CUBE')
        end if
      end select
      ! reads next line
      read(io_input,"(A)",err=999,end=999) key_word
    end do
    ! writes input information to check
    if (do_density_cube) &
      write(io_viewer,100) "generates cube file of electron density"
    if (do_homo_cube) &
      write(io_viewer,100) "generates cube file of HOMO"
    if (do_lumo_cube) &
      write(io_viewer,100) "generates cube file of LUMO"
    if (do_mo_cube) then
      write(io_viewer,100) "generates cube file of MOs"
      write(io_viewer,110) idx_cube_mo
    end if
    write(io_viewer,100) "format: "//trim(cube_format)
    write(io_viewer,120) "origin:", cube_origin
    do ixyz = 1, 3
      write(io_viewer,130) "increments:", cube_num_inc(ixyz), cube_increment(ixyz,:)
    end do
    ! returns the last read keyword back
    word = trim(key_word)
    ! generates cube files later on
    DO_CUBE = do_density_cube .or. do_homo_cube .or. do_lumo_cube .or. do_mo_cube
    return
999 write(io_viewer,100) "failed to process input after reading "// &
                         trim(key_word)//"!"
    call quit('failed to process input')
100 format("gen1int_host_cube_init>> ",A,I8)
110 format("gen1int_host_cube_init>> ",10I5)
120 format("gen1int_host_cube_init>> ",A,3F16.8)
130 format("gen1int_host_cube_init>> ",A,I8,3F16.8)
  end subroutine gen1int_host_cube_init

  !> \brief generates the cube file of the electron density and/or molecular orbitals
  !>        using Gen1Int library
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  subroutine gen1int_host_get_cube(len_work, wrk_space, io_viewer, level_print)
    use london_ao
    use gen1int_matrix
    use gen1int_api
    use gen1int_cube
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
! uses MXCENT, etc.
#include "mxcent.h"
! uses NUCDEP, CHARGE, CORD
#include "nuclei.h"
! uses NCMOT
#include "inforb.h"
! uses LUSIFC
#include "inftap.h"
    integer num_points                                !number of points in cube file
    real(REALK), allocatable :: grid_points(:,:)      !XYZ coordinates of points in cube file
    integer ipoint                                    !incremental recorder over points
    integer ixyz                                      !incremental recorder along XYZ directions
    integer ic, jc, kc                                !incremental recorder over cube increments
    type(matrix) ao_dens(1)                           !AO density matrix
    type(matrix) mo_coef                              !MO coefficients
    real(REALK) start_time                            !start time
    type(prop_comp_t) prop_comp                       !operator of property integrals with non-zero components
    type(nary_tree_t) nary_tree_bra                   !N-ary tree for partial geometric derivatives on bra center
    type(nary_tree_t) nary_tree_ket                   !N-ary tree for partial geometric derivatives on ket center
    type(nary_tree_t) nary_tree_total                 !N-ary tree for total geometric derivatives
    real(REALK), allocatable :: cube_values(:,:,:,:)  !values of points in cube file
    integer io_cube                                   !logical unit number of cube file
    logical close_sirifc                              !if closing SIRIFC afterwards
    integer ierr                                      !error information
    logical found                                     !if found required data from SIRIFC
    integer start_ao, end_ao                          !start and end addresses of AOs
    integer imo                                       !incremental recorder over MOs
#if defined(VAR_MPI)
#include "mpif.h"
#include "iprtyp.h"
    ! wakes up workers
    call MPI_Bcast(GEN1INT_GET_CUBE, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    ! broadcasts level of print
    call MPI_Bcast(level_print, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
#endif
    ! sets the XYZ coordinates of points in cube file
    num_points = product(cube_num_inc)
    allocate(grid_points(3,num_points), stat=ierr)
    if (ierr/=0) then
      call quit("gen1int_host_get_cube>> failed to allocate grid_points!")
    end if
    ipoint = 0
    do ic = 1, cube_num_inc(1)
      do jc = 1, cube_num_inc(2)
        do kc = 1, cube_num_inc(3)
          ! if the origin is (X0,Y0,Z0), and the increment is (X1,Y1,Z1),
          ! then point (I1,I2,I3) has the coordinates:
          ! X-coordinate: X0+(I1-1)*X1+(I2-1)*X2+(I3-1)*X3
          ! Y-coordinate: Y0+(I1-1)*Y1+(I2-1)*Y2+(I3-1)*Y3
          ! Z-coordinate: Z0+(I1-1)*Z1+(I2-1)*Z2+(I3-1)*Z3
          ipoint = ipoint+1
          do ixyz = 1, 3
            grid_points(ixyz,ipoint) = cube_origin(ixyz) &
              + real(ic-1,REALK)*cube_increment(1,ixyz)  &
              + real(jc-1,REALK)*cube_increment(2,ixyz)  &
              + real(kc-1,REALK)*cube_increment(3,ixyz)
          end do
        end do
      end do
    end do
    ! gets the AO density matrix
    if (do_density_cube) then
      call gen1int_host_get_dens(ao_dens(1), len_work, wrk_space)
      ! writes matrix to check
      if (level_print>10) then
        write(io_viewer,"()")
        write(io_viewer,100) "AO density matrix"
        call MatView(A=ao_dens(1), io_viewer=io_viewer)
      end if
      ! creates the operator of overlap integrals
      call gen1int_host_prop_create(NON_LAO, INT_OVERLAP,      &
                                    0, 0,                      &
                                    0, 0, 0,                   &
                                    0, 0, 0,                   &
                                    .false., .false., .false., &
                                    io_viewer, level_print,    &
                                    1, (/1, 1/),               &   !hardcoded for Dalton
                                    prop_comp)
      call gen1int_host_geom_create(0, 0, 0, (/0/),         &
                                    0, 0, 0, (/0/),         &
                                    0, 0, 0, (/0/),         &
                                    io_viewer, level_print, &
                                    nary_tree_bra, nary_tree_ket, nary_tree_total)
      ! evaluates the electron density at points of cube file
      allocate(cube_values(cube_num_inc(3),cube_num_inc(2),cube_num_inc(1),1), &
               stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_get_cube>> failed to allocate cube_values!")
      end if
      cube_values = 0.0_REALK  !necessary to zero
      call Gen1IntAPIPropGetFunExpt(prop_comp=prop_comp,             &
                                    nary_tree_bra=nary_tree_bra,     &
                                    nary_tree_ket=nary_tree_ket,     &
                                    nary_tree_total=nary_tree_total, &
!FIXME: be parallel
                                    !api_comm=MPI_COMM_WORLD,         &
                                    num_points=num_points,           &
                                    grid_points=grid_points,         &
                                    num_dens=1,                      &
                                    ao_dens=ao_dens,                 &
                                    num_ints=1,                      &
                                    val_expt=cube_values,            &
                                    io_viewer=io_viewer,             &
                                    level_print=level_print)
      ! free spaces
      call MatDestroy(A=ao_dens(1))
      call Gen1IntAPIPropDestroy(prop_comp=prop_comp)
      call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
      call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
      call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)
      ! writes cube file
      write(io_viewer,100) "writes cube file of electron density"
      io_cube = -1
      call GPOPEN(io_cube, "density.cube", "unknown", " ", "formatted", &
                  ierr, .false.)
      write(io_cube,"(1X,A)") "molecule density=scf"
      write(io_cube,"(1X,A)") "electron density from total SCF density"
      write(io_cube,"(I5,3F12.6)") NUCDEP, cube_origin
      do ixyz = 1, 3
        write(io_cube,"(I5,3F12.6)") cube_num_inc(ixyz), cube_increment(ixyz,:)
      end do
      do ipoint = 1, NUCDEP
        write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                     CORD(:,ipoint)
      end do
      do ic = 1, cube_num_inc(1)
        do jc = 1, cube_num_inc(2)
          write(io_cube,"(6Es13.5)") cube_values(:,jc,ic,1)
        end do
      end do
      call GPCLOSE(io_cube, "KEEP")
      ! free spaces
      deallocate(cube_values)
    end if
    ! gets the MO coefficients
    if (do_homo_cube .or. do_lumo_cube .or. do_mo_cube) then
      ! gets molecular orbital coefficienct from SIRIFC file
      close_sirifc = LUSIFC<=0
      if (close_sirifc) &
        call GPOPEN(LUSIFC, "SIRIFC", "OLD", " ", "UNFORMATTED", ierr, .false.)
      rewind(LUSIFC)
      ! reads the molecular orbital coefficients
      wrk_space(1:NCMOT) = 0.0_REALK
#ifdef PRG_DIRAC
      print *, 
      call quit('error: RD_SIRIFC not available in DIRAC')
#else
      call RD_SIRIFC("CMO", found, wrk_space(1), wrk_space(NCMOT+1), &
                     len_work-NCMOT)
#endif
      if (.not.found) then
        call quit("gen1int_host_get_cube>> CMO is not found on SIRIFC!")
      end if
      if (close_sirifc) call GPCLOSE(LUSIFC, "KEEP")
      ! gets the number of required MOs
      if (do_homo_cube) num_cube_mo = num_cube_mo+1
      if (do_lumo_cube) num_cube_mo = num_cube_mo+1
      ! sets MO coefficients
      call MatCreate(A=mo_coef, num_row=NBAST, num_col=num_cube_mo, info_mat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_get_cube>> failed to create mo_coef!")
      end if
      if (do_mo_cube) then
        do imo = 1, size(idx_cube_mo)
          if (idx_cube_mo(imo)>NCMOT .or. idx_cube_mo(imo)<1) then
            call quit("gen1int_host_get_cube>> incorrect MOs!")
          end if
          start_ao = (idx_cube_mo(imo)-1)*NBAST
          end_ao = start_ao+NBAST
          start_ao = start_ao+1
          call MatSetValues(mo_coef, 1, NBAST, imo, imo, &
                            wrk_space(start_ao:end_ao), trans=.false.)
        end do
        imo = size(idx_cube_mo)
      else
        imo = 0
      end if
      ! MO coefficients of HOMO
      if (do_homo_cube) then
        start_ao = (NOCCT-1)*NBAST
        end_ao = start_ao+NBAST
        start_ao = start_ao+1
        imo = imo+1
        call MatSetValues(mo_coef, 1, NBAST, imo, imo, &
                          wrk_space(start_ao:end_ao), trans=.false.)
      end if
      ! MO coefficients of LUMO
      if (do_lumo_cube) then
        start_ao = NOCCT*NBAST
        end_ao = start_ao+NBAST
        start_ao = start_ao+1
        imo = imo+1
        call MatSetValues(mo_coef, 1, NBAST, imo, imo, & 
                          values=wrk_space(start_ao:end_ao), trans=.false.)
      end if
      ! evaluates MOs at points in cube file
      allocate(cube_values(cube_num_inc(3),cube_num_inc(2), &
                           cube_num_inc(1),num_cube_mo),    &
               stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_get_cube>> failed to allocate cube_values!")
      end if
#ifdef PRG_DIRAC
      call Gen1IntAPIGetMO(comp_shell=(/1,2/),      &
#else
      call Gen1IntAPIGetMO(comp_shell=(/1/),        &
#endif
                           mo_coef=mo_coef,         &
                           num_points=num_points,   &
                           grid_points=grid_points, &
                           num_derv=1,              &
                           num_mo=num_cube_mo,      &
                           val_mo=cube_values,      &
                           !api_comm=MPI_COMM_WORLD, &
                           gto_type=NON_LAO,        &
                           order_mag=0,             &
                           order_ram=0,             &
                           order_geo=0)
      ! free spaces
      call MatDestroy(A=mo_coef)
      ! resets the number of indices of MOs in cube file
      if (do_mo_cube) then
        num_cube_mo = size(idx_cube_mo)
      else
        num_cube_mo = 0
      end if
      ! writes cube file of MOs
      write(io_viewer,100) "writes cube file of HOMO, LUMO and/or MOs"
      if (do_mo_cube) then
        io_cube = -1
        call GPOPEN(io_cube, "mo.cube", "unknown", " ", "formatted", &
                    ierr, .false.)
        write(io_cube,"(1X,A)") "molecule mo=selected"
        write(io_cube,"(1X,A)") "MO coefficients"
        write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
        do ixyz = 1, 3
          write(io_cube,"(I5,3F12.6)") cube_num_inc(ixyz), cube_increment(ixyz,:)
        end do
        do ipoint = 1, NUCDEP
          write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                       CORD(:,ipoint)
        end do
        write(io_cube,"(10I5)") num_cube_mo, idx_cube_mo
        do ic = 1, cube_num_inc(1)
          do jc = 1, cube_num_inc(2)
!FIXME: the memory access here is awful, might need to change
            write(io_cube,"(6Es13.5)") &
              ((cube_values(kc,jc,ic,imo), imo=1,num_cube_mo), kc=1,cube_num_inc(3))
          end do
        end do
        call GPCLOSE(io_cube, "KEEP")
        imo = num_cube_mo
      else
        imo = 0
      end if
      ! writes cube file of HOMO
      if (do_homo_cube) then
        io_cube = -1
        call GPOPEN(io_cube, "homo.cube", "unknown", " ", "formatted", &
                    ierr, .false.)
        write(io_cube,"(1X,A)") "molecule mo=HOMO"
        write(io_cube,"(1X,A)") "MO coefficients"
        write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
        do ixyz = 1, 3
          write(io_cube,"(I5,3F12.6)") cube_num_inc(ixyz), cube_increment(ixyz,:)
        end do
        do ipoint = 1, NUCDEP
          write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                       CORD(:,ipoint)
        end do
        write(io_cube,"(10I5)") 1, NOCCT
        imo = imo+1
        do ic = 1, cube_num_inc(1)
          do jc = 1, cube_num_inc(2)
            write(io_cube,"(6Es13.5)") cube_values(:,jc,ic,imo)
          end do
        end do
        call GPCLOSE(io_cube, "KEEP")
      end if
      ! writes cube file of LUMO
      if (do_lumo_cube) then
        io_cube = -1
        call GPOPEN(io_cube, "lumo.cube", "unknown", " ", "formatted", &
                    ierr, .false.)
        write(io_cube,"(1X,A)") "molecule mo=LUMO"
        write(io_cube,"(1X,A)") "MO coefficients"
        write(io_cube,"(I5,3F12.6)") -NUCDEP, cube_origin
        do ixyz = 1, 3
          write(io_cube,"(I5,3F12.6)") cube_num_inc(ixyz), cube_increment(ixyz,:)
        end do
        do ipoint = 1, NUCDEP
          write(io_cube,"(I5,4F12.6)") int(CHARGE(ipoint)), CHARGE(ipoint), &
                                       CORD(:,ipoint)
        end do
        write(io_cube,"(10I5)") 1, NOCCT+1
        imo = imo+1
        do ic = 1, cube_num_inc(1)
          do jc = 1, cube_num_inc(2)
            write(io_cube,"(6Es13.5)") cube_values(:,jc,ic,imo)
          end do
        end do
        call GPCLOSE(io_cube, "KEEP")
      end if
      ! free spaces
      deallocate(cube_values)
    end if
    ! frees spaces
    deallocate(grid_points)
    return
100 format("gen1int_host_get_cube>> ",A,I8)
  end subroutine gen1int_host_get_cube

#if defined(VAR_MPI)
  subroutine gen1int_worker_get_cube(len_work, wrk_space, io_viewer, level_print)
    use gen1int_matrix
    use gen1int_api
    use gen1int_cube
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
#include "mpif.h"
    return
  end subroutine gen1int_worker_get_cube
#endif

  !> \brief frees the space taken by the cube files
  !> \author Bin Gao
  !> \date 2012-05-14
  subroutine gen1int_host_cube_finalize()
    use gen1int_cube
    implicit none
    do_density_cube = .false.
    do_homo_cube = .false.
    do_lumo_cube = .false.
    do_mo_cube = .false.
    num_cube_mo = 0
    if (allocated(idx_cube_mo)) deallocate(idx_cube_mo)
    cube_format = "GAUSSIAN"
    cube_origin = 0.0_REALK
    cube_increment = 0.0_REALK
    cube_num_inc = 0
    return
  end subroutine gen1int_host_cube_finalize

  !> \brief terminates Gen1Int interface after all calculations
  !> \author Bin Gao
  !> \date 2011-10-02
  subroutine gen1int_host_finalize
    use gen1int_api
    implicit none
    call Gen1IntAPIDestroy()
    return
  end subroutine gen1int_host_finalize

  !> \brief test suite of Gen1Int interface, enabled by adding the following lines
  !>        in DALTON.INP/DIRAC.INP
  !>        **INTEGRAL
  !>        .GENINT
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  subroutine gen1int_host_test(len_work, wrk_space, io_viewer, level_print)
    use london_ao
    use gen1int_api
    implicit none
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    real(REALK), parameter :: ERR_THRSH = 10.0_REALK**(-8)    !threshold of error
    real(REALK), parameter :: RATIO_THRSH = 10.0_REALK**(-6)  !threshold of ratio to the referenced result
    logical test_failed                                   !indicator if the test failed
    integer num_ao                                        !number of orbitals
    integer, parameter :: NUM_TEST = 11                   !number of tests
    character*20, parameter :: PROP_NAME(NUM_TEST) = &    !names of testing property integrals,
      (/"INT_KIN_ENERGY    ", "INT_OVERLAP       ",  &    !see Gen1int library src/gen1int.F90
        "INT_POT_ENERGY    ", "INT_CART_MULTIPOLE",  &
        "INT_ANGMOM        ", "INT_ONE_HAMIL     ",  &
        "INT_ONE_HAMIL     ", "INT_OVERLAP       ",  &
        "INT_OVERLAP       ", "INT_OVERLAP       ",  &
        "INT_CART_MULTIPOLE"/)
    character*8, parameter :: HERM_PROP(NUM_TEST) = &     !names of property integrals,
      (/"KINENERG", "SQHDOR  ", "POTENERG",         &     !see \fn(PR1IN1) in abacus/her1pro.F
        "CARMOM  ", "ANGMOM  ", "MAGMOM  ",         &
        "ANGMOM  ", "S1MAG   ", "S1MAGL  ",         &
        "S1MAGR  ", "CM1     "/)
    integer, parameter :: GTO_TYPE(NUM_TEST) = &          !
      (/NON_LAO, NON_LAO, NON_LAO, NON_LAO,    &
        NON_LAO, LONDON, NON_LAO, LONDON,      &
        LONDON, LONDON, LONDON/)
    integer, parameter :: ORDER_MOM(NUM_TEST) = &         !order of Cartesian multipole moments
      (/0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1/)
    integer, parameter :: ORDER_MAG_BRA(NUM_TEST) = &     !
      (/0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0/)
    integer, parameter :: ORDER_MAG_KET(NUM_TEST) = &     !
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0/)
    integer, parameter :: ORDER_MAG_TOTAL(NUM_TEST) = &   !
      (/0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1/)
    integer, parameter :: ORDER_RAM_BRA(NUM_TEST) = &     !
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: ORDER_RAM_KET(NUM_TEST) = &     !
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: ORDER_RAM_TOTAL(NUM_TEST) = &   !
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: ORDER_GEO_BRA(NUM_TEST) = &     !
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: ORDER_GEO_KET(NUM_TEST) = &     !
      (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: MAX_NUM_CENT(NUM_TEST) = &      !maximum number of differentiated centers
      (/0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    integer, parameter :: ORDER_GEO_TOTAL(NUM_TEST) = &   !order of total geometric derivatives
      (/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
    logical, parameter :: ADD_SR(NUM_TEST) = &            !
      (/.false., .false., .false., .false.,  &
        .false., .false., .false., .false.,  &
        .false., .false., .false./)
    logical, parameter :: ADD_SO(NUM_TEST) = &            !
      (/.false., .false., .false., .false.,  &
        .false., .false., .false., .false.,  &
        .false., .false., .false./)
    logical, parameter :: ADD_LONDON(NUM_TEST) = &        !
      (/.false., .false., .false., .false.,      &
        .false., .false., .false., .false.,      &
        .false., .false., .false./)
    logical, parameter :: TRIANG(NUM_TEST) = &            !integral matrices are triangularized or squared
      (/.true., .false., .true., .true.,     &
        .true., .true., .true., .true.,      &
        .false., .false., .true./)
    logical, parameter :: SYMMETRIC(NUM_TEST) = &         !integral matrices are symmetric or anti-symmetric
      (/.true., .false., .true., .true.,        &
        .false., .false., .false., .false.,     &
        .false., .false., .false./)
    integer NUM_INTS(NUM_TEST)                            !number of integral matrices
    type(matrix), allocatable :: val_ints(:)              !integral matrices
    logical, parameter :: WRITE_INTS = .false.            !if writing integrals on file
    logical, parameter :: WRITE_EXPT = .false.            !if writing expectation values on file
    integer itest                                         !incremental recorder over different tests
    integer imat, jmat                                    !incremental recorders over matrices
    integer ierr                                          !error information
    ! variables related to \fn(PR1IN1) ...
! uses \var(MXCENT)
#include "mxcent.h"
! number of atomic centers
#include "nuclei.h"
! uses \var(MXCORB)
#include "maxorb.h"
! uses \var(MXQN)
#include "maxaqn.h"
! uses \var(MAXREP)
#include "symmet.h"
! uses FIELD1
#include "efield.h"
!FIXME: having problem of calling \fn(PR1IN1) with \var(GET_EXPT)=.true., which gives wrong integrals
    integer size_int                                      !size of property integrals
    integer start_herm_int, end_herm_int                  !addresses for integrals from \fn(PR1IN1)
    logical, parameter :: GET_EXPT = .false.              !if getting expectation values back
    integer max_typ
    integer, allocatable :: int_rep(:)
    integer lint_ad
    integer, allocatable :: int_adr(:)
    character*8, allocatable :: lb_int(:)
    integer NCOMP
    logical, parameter :: TOFILE = .false.                !if writing integrals on file
    character*6, parameter :: MTFORM = "TRIANG"
    logical, parameter :: DOINT(4) = (/.true., .false., .false., .false./)
    logical, parameter :: PROP_PRINT = .false.            !if printing referenced property integrals
    integer, parameter :: NUM_PQUAD = 40                  !number of integration points for DSO integrals
    integer len_free                                      !length of free Dalton/Dirac workspace
#if !defined (PRG_DIRAC)
    integer base_free                                     !base address of free Dalton/Dirac workspace
#endif
    logical almost_equal                                  !indicates if the results from Gen1Int are
                                                          !almost equal to those from \fn(PR1IN1)
#if !defined (BUILD_OPENRSP)
    ! test suite of matrix module
    call MatTestSuite(test_failed=test_failed, io_viewer=io_viewer, &
                      level_print=level_print, threshold=ERR_THRSH)
#endif
    ! gets the number of atomic orbitals
    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    write(io_viewer,100) "number of orbitals", num_ao
    write(io_viewer,110) "threshold of error", ERR_THRSH
    write(io_viewer,110) "threshold of ratio to the referenced result", RATIO_THRSH
    ! sets the number of integral matrices
    NUM_INTS(1) = 1
    NUM_INTS(2) = 3*NUCDEP
    NUM_INTS(3) = 1
    NUM_INTS(4) = 3
    NUM_INTS(5) = 3
    NUM_INTS(6) = 3
    NUM_INTS(7) = 3
    NUM_INTS(8) = 3
    NUM_INTS(9) = 3
    NUM_INTS(10) = 3
    NUM_INTS(11) = 9
    ! loops over different tests
    do itest = 1, NUM_TEST
      test_failed = .false.
#if defined(PRG_DIRAC)
      if (HERM_PROP(itest)=="POTENERG") cycle
#endif
      ! allocates integral matrices
      allocate(val_ints(NUM_INTS(itest)), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_test>> failed to allocate val_ints!")
      end if
      do imat = 1, NUM_INTS(itest)
        call MatCreate(A=val_ints(imat), num_row=num_ao, info_mat=ierr, &
                       triangular=TRIANG(itest), symmetric=SYMMETRIC(itest))
        if (ierr/=0) then
          call quit("gen1int_host_test>> failed to creates integral matrices!")
        end if
      end do
      ! calculates integrals using Gen1Int
      call gen1int_host_get_int(GTO_TYPE(itest),                                       &
                                trim(PROP_NAME(itest)), ORDER_MOM(itest),              &
                                0,                                                     &
                                ORDER_MAG_BRA(itest), ORDER_MAG_KET(itest),            &
                                ORDER_MAG_TOTAL(itest),                                &
                                ORDER_RAM_BRA(itest), ORDER_RAM_KET(itest),            &
                                ORDER_RAM_TOTAL(itest),                                &
                                MAX_NUM_CENT(itest), ORDER_GEO_BRA(itest), 0, (/0/),   &
                                MAX_NUM_CENT(itest), ORDER_GEO_KET(itest), 0, (/0/),   &
                                MAX_NUM_CENT(itest), ORDER_GEO_TOTAL(itest), 0, (/0/), &
                                ADD_SR(itest), ADD_SO(itest), ADD_LONDON(itest),       &
                                NUM_INTS(itest), val_ints, WRITE_INTS,                 &
                                1, (/1, 1/),                                           & !hardcoded for Dalton
                                io_viewer, level_print)
      ! gets the referenced results from HERMIT
!FIXME: \var(FORQM3)
      if (trim(PROP_NAME(itest))=="DSO") then
        max_typ = (3*NUCDEP)**2
      else
        max_typ = 3*MXCOOR
      end if
      allocate(int_rep(max_typ), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_test>> failed to allocate int_rep!")
      end if
      int_rep = 0
      if (trim(PROP_NAME(itest))=="ELFGRDC" .or. trim(PROP_NAME(itest))=="ELFGRDS") then
         lint_ad = 9*NUCIND*(MAXREP+1)
      else
         lint_ad = max_typ
      end if
      allocate(int_adr(lint_ad), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_test>> failed to allocate int_adr!")
      end if
      int_adr = 0
      allocate(lb_int(max_typ), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_host_test>> failed to allocate lb_int!")
      end if
#if !defined(PRG_DIRAC)
      base_free = 1
#endif
      if (TRIANG(itest)) then
        size_int = num_ao*(num_ao+1)/2
      else
        size_int = num_ao*num_ao
      end if
      end_herm_int = size_int*NUM_INTS(itest)
      len_free = len_work-end_herm_int
      write(io_viewer,100) "gets the referenced results from HERMIT ..."
      ! not equals to 0 so that \fn(PR1IN1) will copy the results back when first calling it
      NCOMP = -1
      if (TRIANG(itest)) then
#if !defined(PRG_DIRAC)
!FIXME: the last argument is symmetric AO density matrix
        if (HERM_PROP(itest)(1:7)=="CM1    ") then
          FIELD1 = 'X-FIELD'
          call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep, & 
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest), & 
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,        & 
                      wrk_space, NCOMP, TOFILE, MTFORM, DOINT,                  & 
                      wrk_space(end_herm_int+1:), GET_EXPT, wrk_space(end_herm_int+1:))
          FIELD1 = 'Y-FIELD'
          call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep,  &
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest),  &
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,         &
                      wrk_space(end_herm_int/3+1), NCOMP, TOFILE, MTFORM, DOINT, &
                      wrk_space(end_herm_int+1:), GET_EXPT, wrk_space(end_herm_int+1:))
          FIELD1 = 'Z-FIELD'
          call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep,    &
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest),    &
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,           &
                      wrk_space(end_herm_int/3*2+1), NCOMP, TOFILE, MTFORM, DOINT, &
                      wrk_space(end_herm_int+1:), GET_EXPT, wrk_space(end_herm_int+1:))
        else
          call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep, &
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest), &
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,        &
                      wrk_space(1:end_herm_int), NCOMP, TOFILE, MTFORM, DOINT,  &
                      wrk_space(end_herm_int+1:), GET_EXPT, wrk_space(end_herm_int+1:))
        end if
#else
        call PR1INT_1(wrk_space(end_herm_int+1:), len_free, int_rep,            &
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest), &
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,        &
                      wrk_space(1:end_herm_int), NCOMP, TOFILE, MTFORM, DOINT)
#endif
      else
#if !defined(PRG_DIRAC)
!FIXME: the last argument is square AO density matrix
        call PR1IN1(wrk_space(end_herm_int+1:), base_free, len_free, int_rep, &
                    int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest), &
                    NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,        &
                    wrk_space(1:end_herm_int), NCOMP, TOFILE, MTFORM, DOINT,  &
                    wrk_space(end_herm_int+1:), GET_EXPT, wrk_space(end_herm_int+1:))
#else
        call PR1INT_1(wrk_space(end_herm_int+1:), len_free, int_rep,            &
                      int_adr, lb_int, HERM_PROP(itest)(1:7), ORDER_MOM(itest), &
                      NUM_PQUAD, TRIANG(itest), PROP_PRINT, level_print,        &
                      wrk_space(1:end_herm_int), NCOMP, TOFILE, MTFORM, DOINT)
#endif
      end if
      deallocate(int_rep)
      deallocate(int_adr)
      deallocate(lb_int)
      ! magnetic derivatives of one-electron Hamiltonian using non-LAOs
      if (PROP_NAME(itest)=="INT_ONE_HAMIL     " .and. &
          GTO_TYPE(itest)==NON_LAO .and.               &
          ORDER_MAG_TOTAL(itest)==1) then
        wrk_space(1:end_herm_int) = -0.5_REALK*wrk_space(1:end_herm_int)
      ! HERMIT uses different sign for the following integrals
      else if (HERM_PROP(itest)=="DPLGRA  " .or. HERM_PROP(itest)=="POTENERG" .or. &
          HERM_PROP(itest)=="NUCSLO  " .or. HERM_PROP(itest)=="PSO     " .or. &
          HERM_PROP(itest)=="ANGMOM  " .or. HERM_PROP(itest)=="SQHDOR  " .or. &
          HERM_PROP(itest)=="MAGMOM  " .or. HERM_PROP(itest)=="S1MAG   " .or. &
          HERM_PROP(itest)=="CM1     ") then
        wrk_space(1:end_herm_int) = -wrk_space(1:end_herm_int)
      end if
      ! checks the results
      write(io_viewer,100) "checks the results of "//trim(PROP_NAME(itest))
      start_herm_int = 0
      do imat = 1, NUM_INTS(itest)
        end_herm_int = start_herm_int+size_int
        start_herm_int = start_herm_int+1
        ! DALTON (i,j): (i-1)*3+j
        ! 1 2 3
        ! 4 5 6
        ! 7 8 9
        ! Gen1Int (i,j): (j-1)*3+i
        ! 1 4 7
        ! 2 5 8
        ! 3 6 9
        if (HERM_PROP(itest)=="CM1     ") then
            jmat = mod(imat,3)
            if (jmat==0) then
                jmat = 6+imat/3
            else
                jmat = (jmat-1)*3+imat/3+1
            end if
            call MatArrayAlmostEqual(A=val_ints(jmat),                              &
                                     values=wrk_space(start_herm_int:end_herm_int), &
                                     io_viewer=io_viewer,                           &
                                     almost_equal=almost_equal,                     &
                                     triangular=TRIANG(itest),                      &
                                     symmetric=SYMMETRIC(itest),                    &
                                     threshold=ERR_THRSH,                           &
                                     ratio_thrsh=RATIO_THRSH)
            call MatDestroy(A=val_ints(jmat))
        else
            call MatArrayAlmostEqual(A=val_ints(imat),                              &
                                     values=wrk_space(start_herm_int:end_herm_int), &
                                     io_viewer=io_viewer,                           &
                                     almost_equal=almost_equal,                     &
                                     triangular=TRIANG(itest),                      &
                                     symmetric=SYMMETRIC(itest),                    &
                                     threshold=ERR_THRSH,                           &
                                     ratio_thrsh=RATIO_THRSH)
            call MatDestroy(A=val_ints(imat))
        end if
        if (.not. test_failed) test_failed = .not.almost_equal
        start_herm_int = end_herm_int
      end do
      deallocate(val_ints)
      if (test_failed) then
        write(io_viewer,100) "test of "//trim(PROP_NAME(itest))//".Gen1Int failed!"
      else
        write(io_viewer,100) "test of "//trim(PROP_NAME(itest))//".Gen1Int passed!"
      end if
    end do
    return
100 format("gen1int_host_test>> ",A,I8)
110 format("gen1int_host_test>> ",A,Es16.6)
  end subroutine gen1int_host_test

  ! private subroutines

  !> \brief creates the operator of property integrals, and broadcasts input
  !>        arguments to worker processors
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \return prop_comp is the operator of property integrals with non-zero components
  !> \note see \fn(gen1int_host_get_int) for the explanation of other arguments
  subroutine gen1int_host_prop_create(gto_type, prop_name, order_mom, &
                                      order_elec,                     &
                                      order_mag_bra, order_mag_ket,   &
                                      order_mag_total,                &
                                      order_ram_bra, order_ram_ket,   &
                                      order_ram_total,                &
                                      add_sr, add_so, add_london,     &
                                      io_viewer, level_print,         &
                                      nr_active_blocks,               &
                                      active_component_pairs,         &
                                      prop_comp)
    use gen1int_api
    implicit none
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
    integer,           intent(in)    :: io_viewer
    integer,           intent(in)    :: level_print
    integer,           intent(in)    :: nr_active_blocks
    integer,           intent(in)    :: active_component_pairs(*)
    type(prop_comp_t), intent(inout) :: prop_comp
    integer len_name  !length of property name
    integer ierr      !error information
#if defined(VAR_MPI)
#include "mpif.h"
    ! broadcasts input arguments
    call MPI_Bcast(gto_type, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    len_name = len_trim(prop_name)
    if (len_name>MAX_LEN_STR) then
      call quit("gen1int_host_prop_create>> too long property name!")
    end if
    call MPI_Bcast(len_name, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(prop_name(1:len_name), len_name, MPI_CHARACTER, &
                   MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mom,       1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_elec,      1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_bra,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_ket,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_bra,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_ket,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nr_active_blocks,       1,                  MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(active_component_pairs, nr_active_blocks*2, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_sr,     1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_so,     1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_london, 1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
#endif
    ! creates the operator of property integrals
    call Gen1IntAPIPropCreate(gto_type, prop_name, order_mom, &
                              order_elec,                     &
                              order_mag_bra, order_mag_ket,   &
                              order_mag_total,                &
                              order_ram_bra, order_ram_ket,   &
                              order_ram_total,                &
                              add_sr, add_so, add_london,     &
                              nr_active_blocks,               &
                              active_component_pairs,         &
                              prop_comp)
    if (level_print>=5) then
      call Gen1IntAPIPropView(prop_comp=prop_comp, io_viewer=io_viewer)
    end if
    return
  end subroutine gen1int_host_prop_create

#if defined(VAR_MPI)
  !> \brief creates the operator of property integrals on worker processors by the arguments from manager
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return prop_comp is the operator of property integrals with non-zero components
  subroutine gen1int_worker_prop_create(io_viewer, level_print, prop_comp)
    use gen1int_api
    implicit none
    integer,           intent(in)    :: io_viewer
    integer,           intent(in)    :: level_print
    type(prop_comp_t), intent(inout) :: prop_comp
    ! local variables, see the explanation of input arguments in \fn(gen1int_host_get_int)
    integer gto_type
    character(MAX_LEN_STR) prop_name
    integer order_mom
    integer order_elec
    integer order_mag_bra
    integer order_mag_ket
    integer order_mag_total
    integer order_ram_bra
    integer order_ram_ket
    integer order_ram_total
    integer              :: nr_active_blocks
    integer, allocatable :: active_component_pairs(:)
    logical add_sr
    logical add_so
    logical add_london
    integer len_name  !length of property name
    integer ierr      !error information
#include "mpif.h"
    ! gets input arguments from manager
    call MPI_Bcast(gto_type, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(len_name, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (len_name>MAX_LEN_STR) then
      call quit("gen1int_worker_prop_create>> too long property name!")
    end if
    prop_name = ""
    call MPI_Bcast(prop_name(1:len_name), len_name, MPI_CHARACTER, &
                   MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mom,       1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_elec,      1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_bra,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_ket,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_mag_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_bra,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_ket,   1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_ram_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nr_active_blocks,       1,                  MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    allocate(active_component_pairs(nr_active_blocks*2), stat=ierr)
    if (ierr /= 0) then
      call quit("gen1int_worker_prop_create>> failed to allocate active_component_pairs!")
    end if
    call MPI_Bcast(active_component_pairs, nr_active_blocks*2, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_sr,     1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_so,     1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(add_london, 1, MPI_LOGICAL, MANAGER, MPI_COMM_WORLD, ierr)
    ! creates the operator of property integrals
    call Gen1IntAPIPropCreate(gto_type, prop_name, order_mom, &
                              order_elec,                     &
                              order_mag_bra, order_mag_ket,   &
                              order_mag_total,                &
                              order_ram_bra, order_ram_ket,   &
                              order_ram_total,                &
                              add_sr, add_so, add_london,     &
                              nr_active_blocks,               &
                              active_component_pairs,         &
                              prop_comp)
    if (level_print>=5) then
      call Gen1IntAPIPropView(prop_comp=prop_comp, io_viewer=io_viewer)
    end if
    return
  end subroutine gen1int_worker_prop_create
#endif

  !> \brief creates N-ary tree for geometric derivatives on manager processor, and broadcasts
  !>        input arguments to worker processors
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \return nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \return nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \return nary_tree_total is the N-ary tree for total geometric derivatives
  !> \note see \fn(gen1int_host_get_int) for the explanation of other arguments
  subroutine gen1int_host_geom_create(max_ncent_bra, order_geo_bra,  &
                                      num_atoms_bra, idx_atoms_bra,  &
                                      max_ncent_ket, order_geo_ket,  &
                                      num_atoms_ket, idx_atoms_ket,  &
                                      max_num_cent, order_geo_total, &
                                      num_geo_atoms, idx_geo_atoms,  &
                                      io_viewer, level_print,        &
                                      nary_tree_bra, nary_tree_ket, nary_tree_total)
    use gen1int_api
    implicit none
    integer,           intent(in)    :: max_ncent_bra
    integer,           intent(in)    :: order_geo_bra
    integer,           intent(in)    :: num_atoms_bra
    integer,           intent(in)    :: idx_atoms_bra(*)
    integer,           intent(in)    :: max_ncent_ket
    integer,           intent(in)    :: order_geo_ket
    integer,           intent(in)    :: num_atoms_ket
    integer,           intent(in)    :: idx_atoms_ket(*)
    integer,           intent(in)    :: max_num_cent
    integer,           intent(in)    :: order_geo_total
    integer,           intent(in)    :: num_geo_atoms
    integer,           intent(in)    :: idx_geo_atoms(*)
    integer,           intent(in)    :: io_viewer
    integer,           intent(in)    :: level_print
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer ierr  !error information
#if defined(VAR_MPI)
#include "mpif.h"
    ! broadcasts information of N-ary trees
    call MPI_Bcast(max_ncent_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_atoms_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_atoms_bra>0) then
      call MPI_Bcast(idx_atoms_bra, num_atoms_bra, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
    call MPI_Bcast(max_ncent_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_atoms_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_atoms_ket>0) then
      call MPI_Bcast(idx_atoms_ket, num_atoms_ket, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
    call MPI_Bcast(max_num_cent, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_geo_atoms, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_geo_atoms>0) then
      call MPI_Bcast(idx_geo_atoms, num_geo_atoms, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
#endif
    ! creates N-ary tree for geometric derivatives
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_ncent_bra,  &
                                  order_geo=order_geo_bra,     &
                                  num_geo_atoms=num_atoms_bra, &
                                  idx_geo_atoms=idx_atoms_bra, &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_ncent_ket,  &
                                  order_geo=order_geo_ket,     &
                                  num_geo_atoms=num_atoms_ket, &
                                  idx_geo_atoms=idx_atoms_ket, &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_num_cent,   &
                                  order_geo=order_geo_total,   &
                                  num_geo_atoms=num_geo_atoms, &
                                  idx_geo_atoms=idx_geo_atoms, &
                                  nary_tree=nary_tree_total)
    if (level_print>=5) then
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
    end if
    return
  end subroutine gen1int_host_geom_create

#if defined(VAR_MPI)
  !> \brief creates N-ary tree for geometric derivatives on worker processors by the arguments from manager
  !> \author Bin Gao
  !> \date 2012-05-09
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \return nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \return nary_tree_total is the N-ary tree for total geometric derivatives
  subroutine gen1int_worker_geom_create(io_viewer, level_print, nary_tree_bra, nary_tree_ket, nary_tree_total)
    use gen1int_api
    implicit none
    integer,           intent(in)    :: io_viewer
    integer,           intent(in)    :: level_print
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    ! local variables, see the explanation of input arguments in \fn(gen1int_host_get_int)
    integer max_ncent_bra
    integer order_geo_bra
    integer num_atoms_bra
    integer, allocatable :: idx_atoms_bra(:)
    integer max_ncent_ket
    integer order_geo_ket
    integer num_atoms_ket
    integer, allocatable :: idx_atoms_ket(:)
    integer max_num_cent
    integer order_geo_total
    integer num_geo_atoms
    integer, allocatable :: idx_geo_atoms(:)
    integer ierr  !error information
#include "mpif.h"
    ! gets information of N-ary trees from manager node
    call MPI_Bcast(max_ncent_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_atoms_bra, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_atoms_bra>0) then
      allocate(idx_atoms_bra(num_atoms_bra), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_worker_geom_create>> failed to allocate idx_atoms_bra!")
      end if
      call MPI_Bcast(idx_atoms_bra, num_atoms_bra, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
    call MPI_Bcast(max_ncent_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_atoms_ket, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_atoms_ket>0) then
      allocate(idx_atoms_ket(num_atoms_ket), stat=ierr)
      if (ierr/=0) then          
        call quit("gen1int_worker_geom_create>> failed to allocate idx_atoms_ket!")
      end if
      call MPI_Bcast(idx_atoms_ket, num_atoms_ket, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
    call MPI_Bcast(max_num_cent, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(order_geo_total, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(num_geo_atoms, 1, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    if (num_geo_atoms>0) then
      allocate(idx_geo_atoms(num_geo_atoms), stat=ierr)
      if (ierr/=0) then
        call quit("gen1int_worker_geom_create>> failed to allocate idx_geo_atoms!")
      end if
      call MPI_Bcast(idx_geo_atoms, num_geo_atoms, MPI_INTEGER, MANAGER, MPI_COMM_WORLD, ierr)
    end if
    ! creates N-ary trees for geometric derivatives
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_ncent_bra,  &
                                  order_geo=order_geo_bra,     &
                                  num_geo_atoms=num_atoms_bra, &
                                  idx_geo_atoms=idx_atoms_bra, &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_ncent_ket,  &
                                  order_geo=order_geo_ket,     &
                                  num_geo_atoms=num_atoms_ket, &
                                  idx_geo_atoms=idx_atoms_ket, &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=max_num_cent,   &
                                  order_geo=order_geo_total,   &
                                  num_geo_atoms=num_geo_atoms, &
                                  idx_geo_atoms=idx_geo_atoms, &
                                  nary_tree=nary_tree_total)
    if (level_print>=5) then
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_bra, io_viewer=io_viewer)
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_ket, io_viewer=io_viewer)
      call Gen1IntAPINaryTreeView(nary_tree=nary_tree_total, io_viewer=io_viewer)
    end if
    return
  end subroutine gen1int_worker_geom_create
#endif

  !> \brief gets the atomic orbital (AO) density matrix
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param len_work is the length of Dalton/Dirac workspace
  !> \param wrk_space is the Dalton/Dirac workspace
  !> \return ao_dens is the AO density matrix
  subroutine gen1int_host_get_dens(ao_dens, len_work, wrk_space)
    use gen1int_matrix
    implicit none
    type(matrix), intent(inout) :: ao_dens
    integer, intent(in) :: len_work
    real(REALK), intent(inout) :: wrk_space(len_work)
! uses NCMOT, NBAST, NNASHX, N2BASX
#include "inforb.h"
! uses LUSIFC
#include "inftap.h"
    integer start_dv_mo     !start of active part of one-electron density matrix (MO) in the workspace
    integer start_dv_ao     !start of active part of one-electron density matrix (AO) in the workspace
    integer start_ao_dens   !start of AO density matrix
    integer start_left_wrk  !start of left workspace
    integer len_left_wrk    !length of left workspace
    logical get_dc          !if calculating DCAO
    logical get_dv          !if calculating DVAO
    logical close_sirifc    !if closing SIRIFC afterwards
    integer ierr            !error information
    logical found           !if found required data from SIRIFC
    ! start of active part of one-electron density matrix (MO)
    start_dv_mo = 1+NCMOT
    ! start of active part of one-electron density matrix (AO)
    start_dv_ao = start_dv_mo+1
    ! start of AO density matrix
    start_ao_dens = start_dv_ao+1
    ! start of left workspace
    start_left_wrk = start_ao_dens+N2BASX+1
    ! only calculates DCAO
    get_dc = .true.
    get_dv = .false.
    ! lenght of the left workspace
    len_left_wrk = len_work-start_left_wrk+1
    ! checks if the workspace is enough
    if (len_left_wrk<0) &
      call STOPIT("GEN1INT", "gen1int_host_get_dens", start_left_wrk-1, len_work)
    ! opens file SIRIFC
    close_sirifc = LUSIFC<=0
    if (close_sirifc) &
      call GPOPEN(LUSIFC, "SIRIFC", "OLD", " ", "UNFORMATTED", ierr, .false.)
    rewind(LUSIFC)
    ! reads the molecular orbital coefficients
    wrk_space(1:NCMOT) = 0.0_REALK
#ifdef PRG_DIRAC
      print *, 'error: RD_SIRIFC not available in DIRAC'
      call quit('error: RD_SIRIFC not available in DIRAC')
#else
    call RD_SIRIFC("CMO", found, wrk_space(1), wrk_space(start_left_wrk), &
                   len_left_wrk)
#endif
    if (.not.found) then
      call quit("gen1int_host_get_dens>> CMO is not found on SIRIFC!")
    end if
    ! reads active part of one-electron density matrix (MO)
    if (get_dv) then
      wrk_space(start_dv_mo:start_dv_mo+NNASHX-1) = 0.0_REALK
#ifdef PRG_DIRAC
      print *,'error: RD_SIRIFC not available in DIRAC' 
      call quit('error: RD_SIRIFC not available in DIRAC')
#else
      call RD_SIRIFC("DV", found, wrk_space(start_dv_mo), &
                     wrk_space(start_left_wrk), len_left_wrk)
#endif
      if (.not.found) then
        call quit("gen1int_host_get_dens>> DV is not found on SIRIFC!")
      end if
      wrk_space(start_dv_ao:start_dv_ao+N2BASX-1) = 0.0_REALK
    end if
    if (close_sirifc) call GPCLOSE(LUSIFC, "KEEP")
    ! gets the AO density matrix, using
    !
    ! FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LFRSAV)
    ! Input:
    !   GETDC   if true calculate DCAO
    !   GETDV   if true calculate DVAO
    !   CMO(*)  molecular orbital coefficients
    !   DV(*)   active part of one-electron density matrix (over MO's)
    ! Scratch:
    !   WRK(LFRSAV)
    call FCKDEN(get_dc, get_dv, wrk_space(start_ao_dens), wrk_space(start_dv_ao), &
                wrk_space(1), wrk_space(start_dv_mo), wrk_space(start_left_wrk),  &
                len_left_wrk)
    ! sums DCAO and DVAO
    if (get_dv) then
      start_dv_ao = start_dv_ao-1
      do ierr = 0, N2BASX-1
        wrk_space(start_ao_dens+ierr) = wrk_space(start_ao_dens+ierr) &
                                      + wrk_space(start_dv_ao+ierr+1)
      end do
    end if
    ! sets AO density matrix
    call MatCreate(A=ao_dens, num_row=NBAST, info_mat=ierr)
    if (ierr/=0) then
      call quit("gen1int_host_get_dens>> failed to create ao_dens!")
    end if
    call MatSetValues(ao_dens, 1, NBAST, 1, NBAST, &
                      values=wrk_space(start_ao_dens), trans=.false.)
    return
  end subroutine gen1int_host_get_dens
