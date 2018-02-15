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
!...  This file takes care the atomic orbital (AO) sub-shell used in Gen1Int interface.
!
!...  2013-05-16, Bin Gao
!...  * fix the bug of calculating expectation values in parallel
!
!...  2013-05-04, Bin Gao
!...  * fix the bug of returning wrong partial geometric derivatives
!
!...  2012-05-10, Bin Gao
!...  * implements manager/worker parallelization scheme in Gen1IntShellGetIntExpt and
!...    Gen1IntShellGetMO
!
!...  2012-05-09, Radovan Bast
!...  * implements large and small components
!
!...  2012-03-11, Bin Gao
!...  * adds subroutine Gen1IntShellGetMO to calculate molecular orbitals at grid points
!
!...  2011-10-02, Bin Gao
!...  * first version

#include "gen1int_host.h"

!> \brief defines the AO sub-shell used in Gen1Int interface and corresponding subroutines
!> \author Bin Gao
!> \date 2011-10-02
module gen1int_shell

  ! Fortran 90 module of Gen1Int library
  use gen1int

  implicit none

  ! contracted Cartesian or spherical GTOs, or an atomic orbital (AO) sub-shell
  ! here, an AO sub-shell contains the orbitals with the same exponents, and the
  ! same/different contraction coefficients; for instance 2p and 3p sub-shells
  type, public :: sub_shell_t
    private
    ! if spherical GTOs
    logical spher_gto
    ! index of atomic center where this AO sub-shell locates
    integer idx_cent
    ! coordinates of the atomic center
    real(REALK) coord_cent(3)
    ! angular number
    integer ang_num
    ! number of primitive GTOs
    integer num_prim
    ! exponents of primitive GTOs
    real(REALK), allocatable :: exponents(:)
    ! number of contractions
    integer num_contr
    ! contraction coefficients
    real(REALK), allocatable :: contr_coef(:,:)
    ! number of atomic orbitals
    integer num_ao
    !-Dalton uses the same sequence of SGTOs as those in Gen1Int except for p-shell
    !-which could be done by hand coding; Dirac does not use SGTOs
    !-! magnetic numbers if spherical GTOs
    !-integer, allocatable :: mag_num(:)
    ! Cartesian powers if Cartesian GTOs
    integer, allocatable :: powers(:,:)
    ! base index of the atomic orbitals in this sub-shell, i.e., the index of the
    ! first atomic orbital in this shell minus 1
    integer base_idx
  end type sub_shell_t

  public :: Gen1IntShellCreate
#if defined(VAR_MPI)
  public :: Gen1IntShellBcast
#endif
  public :: Gen1IntShellView
  public :: Gen1IntShellGetNumAO
  public :: Gen1IntShellGetNumContr
  public :: Gen1IntShellGetRangeAO
  public :: Gen1IntShellGetIntExpt
  public :: Gen1IntShellGetFunExpt
  public :: Gen1IntShellGetMO
  public :: Gen1IntShellDestroy

  contains

  !> \brief creates an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param spher_gto indicates if spherical GTOs
  !> \param idx_cent is the index of atomic center where this AO sub-shell locates
  !> \param coord_cent contains the coordinates of the atomic center
  !> \param ang_num is the angular number
  !> \param num_prim is the number of primitive GTOs
  !> \param exponents contains the exponents of primitive GTOs
  !> \param num_contr is the number of contractions
  !> \param contr_coef contains the contraction coefficients
  !> \param powers contains the Cartesian powers if Cartesian GTOs
  !> \param last_shell is the last sub-shell before this AO sub-shell
  !> \return sub_shell is the initialized AO sub-shell
  subroutine Gen1IntShellCreate(spher_gto, idx_cent, coord_cent, ang_num,   &
                                num_prim, exponents, num_contr, contr_coef, &
                                powers, last_shell, sub_shell)
    logical, intent(in) :: spher_gto
    integer, intent(in) :: idx_cent
    real(REALK), intent(in) :: coord_cent(3)
    integer, intent(in) :: ang_num
    integer, intent(in) :: num_prim
    real(REALK), intent(in) :: exponents(num_prim)
    integer, intent(in) :: num_contr
    real(REALK), intent(in) :: contr_coef(num_contr,num_prim)
    integer, optional, intent(in) :: powers(3,(ang_num+1)*(ang_num+2)/2)
    type(sub_shell_t), optional, intent(in) :: last_shell
    type(sub_shell_t), intent(inout) :: sub_shell
    integer iao         !incremental recorder
    integer xpow, ypow  !incremental recorders over xy powers
    integer ierr        !error information
    sub_shell%spher_gto = spher_gto
    sub_shell%idx_cent = idx_cent
    sub_shell%coord_cent = coord_cent
    sub_shell%ang_num = ang_num
    sub_shell%num_prim = num_prim
    allocate(sub_shell%exponents(num_prim), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntShellCreate>> failed to allocate exponents!")
    end if
    sub_shell%exponents = exponents
    sub_shell%num_contr = num_contr
    allocate(sub_shell%contr_coef(num_contr,num_prim), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntShellCreate>> failed to allocate contr_coef!")
    end if
    sub_shell%contr_coef = contr_coef
    ! spherical GTOs
    if (spher_gto) then
      sub_shell%num_ao = 2*ang_num+1
      !-allocate(sub_shell%mag_num(sub_shell%num_ao), stat=ierr)
      !-if (ierr/=0) then
      !-  stop "Gen1IntShellCreate>> failed to allocate mag_num!"
      !-end if
      !-if (present(mag_num)) then
      !-  sub_shell%mag_num = mag_num
      !-else
      !-  ! Dalton's order of SGTOs
      !-  if (ang_num==1) then
      !-    sub_shell%mag_num(1) = 1   !px
      !-    sub_shell%mag_num(2) = -1  !py
      !-    sub_shell%mag_num(3) = 0   !pz
      !-  else
      !-    do iao = -ang_num, ang_num
      !-      sub_shell%mag_num(ang_num+iao+1) = iao
      !-    end do
      !-  end if
      !-end if
    ! Cartesian GTOs
    else
      sub_shell%num_ao = (ang_num+1)*(ang_num+2)/2
      allocate(sub_shell%powers(3,sub_shell%num_ao), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntShellCreate>> failed to allocate powers!")
      end if
      if (present(powers)) then
        sub_shell%powers = powers
      ! Dalton/Dirac's order of CGTOs, for instance, dxx, dxy, dxz, dyy, dyz, dzz
      else
        iao = 0
        do xpow = ang_num, 0, -1
          do ypow = ang_num-xpow, 0, -1
            iao = iao+1
            sub_shell%powers(1,iao) = xpow
            sub_shell%powers(2,iao) = ypow
            sub_shell%powers(3,iao) = ang_num-(xpow+ypow)
          end do
        end do
      end if  
    end if
    if (present(last_shell)) then
      sub_shell%base_idx = last_shell%base_idx &
                         + last_shell%num_ao*last_shell%num_contr
    ! this is the first sub-shell
    else
      sub_shell%base_idx = 0
    end if
  end subroutine Gen1IntShellCreate

#if defined(VAR_MPI)
  !> \brief broadcasts AO sub-shells
  !> \author Bin Gao
  !> \date 2012-05-13
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  !> \param root is the root processor which broadcasts the AO sub-shells
  !> \param api_comm is the MPI communicator
  subroutine Gen1IntShellBcast(num_shells, sub_shells, root, api_comm)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(inout) :: sub_shells(num_shells)
    integer, intent(in) :: root
    integer, intent(in) :: api_comm
#include "mpif.h"
    integer rank_proc         !rank of processor
    logical bcast_mag_powers  !if broadcasting magnetic numbers or Cartesian powers
    integer ishell            !incremental recorder sub-shells
    integer ierr              !error information
    ! gets the rank of processor
    call MPI_Comm_rank(api_comm, rank_proc, ierr)
    do ishell = 1, num_shells
      ! broadcasts
      call MPI_Bcast(sub_shells(ishell)%spher_gto, 1, MPI_LOGICAL, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%idx_cent, 1, MPI_INTEGER, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%coord_cent, 3, MPI_REALK, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%ang_num, 1, MPI_INTEGER, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%num_prim, 1, MPI_INTEGER, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%num_contr, 1, MPI_INTEGER, root, api_comm, ierr)
      if (rank_proc==root) then
        ! spherical GTOs
        if (sub_shells(ishell)%spher_gto) then
          !-bcast_mag_powers = allocated(sub_shells(ishell)%mag_num)
          !-call MPI_Bcast(bcast_mag_powers, 1, MPI_LOGICAL, root, api_comm, ierr)
          !-if (bcast_mag_powers) then
          !-  call MPI_Bcast(sub_shells(ishell)%mag_num, sub_shells(ishell)%num_ao, &
          !-                 MPI_INTEGER, root, api_comm, ierr)
          !-end if
        ! Cartesian GTOs
        else
          bcast_mag_powers = allocated(sub_shells(ishell)%powers)
          call MPI_Bcast(bcast_mag_powers, 1, MPI_LOGICAL, root, api_comm, ierr)
          if (bcast_mag_powers) then
            call MPI_Bcast(sub_shells(ishell)%powers, 3*sub_shells(ishell)%num_ao, &
                           MPI_INTEGER, root, api_comm, ierr)
          end if
        end if
      else
        ! spherical GTOs
        if (sub_shells(ishell)%spher_gto) then
          sub_shells(ishell)%num_ao = 2*sub_shells(ishell)%ang_num+1
          !-call MPI_Bcast(bcast_mag_powers, 1, MPI_LOGICAL, root, api_comm, ierr)
          !-if (bcast_mag_powers) then
          !-  allocate(sub_shells(ishell)%mag_num(sub_shells(ishell)%num_ao), stat=ierr)
          !-  if (ierr/=0) then
          !-    stop "Gen1IntShellBcast>> failed to allocate mag_num!"
          !-  end if
          !-  call MPI_Bcast(sub_shells(ishell)%mag_num, sub_shells(ishell)%num_ao, &
          !-                 MPI_INTEGER, root, api_comm, ierr)
          !-end if
        ! Cartesian GTOs
        else
          sub_shells(ishell)%num_ao = (sub_shells(ishell)%ang_num+1) &
                                    * (sub_shells(ishell)%ang_num+2)/2
          call MPI_Bcast(bcast_mag_powers, 1, MPI_LOGICAL, root, api_comm, ierr)
          if (bcast_mag_powers) then
            allocate(sub_shells(ishell)%powers(3,sub_shells(ishell)%num_ao), stat=ierr)
            if (ierr/=0) then
              call quit("Gen1IntShellBcast>> failed to allocate powers!")
            end if
            call MPI_Bcast(sub_shells(ishell)%powers, 3*sub_shells(ishell)%num_ao, &
                           MPI_INTEGER, root, api_comm, ierr)
          end if
        end if
        ! allocates memory for exponents and contraction coefficients
        allocate(sub_shells(ishell)%exponents(sub_shells(ishell)%num_prim), stat=ierr)
        if (ierr/=0) then
          call quit("Gen1IntShellBcast>> failed to allocate exponents!")
        end if
        allocate(sub_shells(ishell)%contr_coef(sub_shells(ishell)%num_contr, &
                                               sub_shells(ishell)%num_prim), stat=ierr)
        if (ierr/=0) then
          call quit("Gen1IntShellBcast>> failed to allocate contr_coef!")
        end if
      end if
      call MPI_Bcast(sub_shells(ishell)%exponents, sub_shells(ishell)%num_prim, &
                     MPI_REALK, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%contr_coef,                            &
                     sub_shells(ishell)%num_contr*sub_shells(ishell)%num_prim, &
                     MPI_REALK, root, api_comm, ierr)
      call MPI_Bcast(sub_shells(ishell)%base_idx, 1, MPI_INTEGER, root, api_comm, ierr)
    end do
  end subroutine Gen1IntShellBcast
#endif

  !> \brief visualizes the information of several AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  !> \param io_viewer is the logical unit number of the viewer
  subroutine Gen1IntShellView(num_shells, sub_shells, io_viewer)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(in) :: sub_shells(num_shells)
    integer, intent(in) :: io_viewer
    integer ishell, iprim  !incremental recorders
    write(io_viewer,100) "number of AO sub-shells", num_shells
    do ishell = 1, num_shells
      write(io_viewer,100) "sub-shell", ishell
      write(io_viewer,100) "index of atomic center of the sub-shell", &
                           sub_shells(ishell)%idx_cent
      write(io_viewer,110) "coordinates of the atomic center", &
                           sub_shells(ishell)%coord_cent
      write(io_viewer,100) "angular number", sub_shells(ishell)%ang_num
      write(io_viewer,100) "number of primitive GTOs", sub_shells(ishell)%num_prim
      write(io_viewer,100) "number of contractions", sub_shells(ishell)%num_contr
      write(io_viewer,100) "    exponents      coefficients"
      do iprim = 1, sub_shells(ishell)%num_prim
        write(io_viewer,120)                   &
          sub_shells(ishell)%exponents(iprim), &
          sub_shells(ishell)%contr_coef(1:min(4,sub_shells(ishell)%num_contr),iprim)
        write(io_viewer,130) &
          sub_shells(ishell)%contr_coef(5:sub_shells(ishell)%num_contr,iprim)
      end do
      if (sub_shells(ishell)%spher_gto) then
        write(io_viewer,100) "SGTOs used"
        !-write(io_viewer,140) sub_shells(ishell)%mag_num
      else
        write(io_viewer,100) "CGTOs used"
        write(io_viewer,150) sub_shells(ishell)%powers
      end if
      write(io_viewer,100) "base index of the orbitals", sub_shells(ishell)%base_idx
    end do
100 format("Gen1IntShellView>> "A,2I8)
110 format("Gen1IntShellView>> "A,3F16.8)
120 format("Gen1IntShellView>> ",5Es16.8)
130 format("Gen1IntShellView>> ",16X,4Es16.8)
140 format("Gen1IntShellView>> magnetic numbers>> ",8I4)
150 format("Gen1IntShellView>> Cartesian powers>> ",4(3I4,2X))
  end subroutine Gen1IntShellView

  !> \brief gets the number of atomic orbitals of an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param sub_shell is the AO sub-shell
  !> \return num_ao is the number of atomic orbitals
  subroutine Gen1IntShellGetNumAO(sub_shell, num_ao)
    type(sub_shell_t), intent(in) :: sub_shell
    integer, intent(out) :: num_ao
    num_ao = sub_shell%num_ao
  end subroutine Gen1IntShellGetNumAO

  !> \brief gets the number of contractions of an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-04
  !> \param sub_shell is the AO sub-shell
  !> \return num_contr is the number of contractions
  subroutine Gen1IntShellGetNumContr(sub_shell, num_contr)
    type(sub_shell_t), intent(in) :: sub_shell
    integer, intent(out) :: num_contr
    num_contr = sub_shell%num_contr
  end subroutine Gen1IntShellGetNumContr

  !> \brief gets the indices of the first and last orbtials of an AO sub-shell
  !> \author Bin Gao
  !> \date 2011-10-06
  !> \param sub_shell is the AO sub-shell
  !> \return idx_first is the index of the first orbital
  !> \return idx_last is the index of the last orbital
  subroutine Gen1IntShellGetRangeAO(sub_shell, idx_first, idx_last)
    type(sub_shell_t), intent(in) :: sub_shell
    integer, intent(out) :: idx_first
    integer, intent(out) :: idx_last
    idx_first = sub_shell%base_idx+1
    idx_last = sub_shell%base_idx+sub_shell%num_ao*sub_shell%num_contr
  end subroutine Gen1IntShellGetRangeAO

  !> \brief calculates property integral matrices and/or expectation values for
  !>        given AO sub-shells on bra and ket centers
  !> \author Bin Gao and Radovan Bast
  !> \date 2011-10-07
  !> \param num_shells_bra is the number of AO sub-shells on bra center
  !> \param sub_shells_bra are the AO sub-shells on bra center
  !> \param num_shells_ket is the number of AO sub-shells on ket center
  !> \param sub_shells_ket are the AO sub-shells on ket center
  !> \param same_braket indicates if the AO sub-shells are the same on bra and ket centers
  !> \param one_prop contains the information of one-electron property integrals
  !> \param nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \param nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \param nary_tree_total is the N-ary tree for total geometric derivatives
  !> \param api_comm is the MPI communicator
  !> \param num_prop is the number of property integrals including various derivatives (except
  !>        geometric derivatives)
  !> \param num_geo_bra is number of geometric derivatives on bra center
  !> \param num_geo_ket is number of geometric derivatives on ket center
  !> \param num_geo_total is number of total geometric derivatives
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \return val_ints contains the integral matrices
  !> \return val_expt contains the expectation values
  !> \note the arrangement of var(val_ints) and \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntShellGetIntExpt(num_shells_bra, sub_shells_bra, &
                                    num_shells_ket, sub_shells_ket, &
                                    same_braket, one_prop,          &
                                    nary_tree_bra, nary_tree_ket,   &
                                    nary_tree_total, api_comm,      &
                                    num_prop, num_geo_bra,          &
                                    num_geo_ket, num_geo_total,     &
                                    val_ints, write_ints,           &
                                    num_dens, ao_dens, val_expt, write_expt)
    ! matrix module
    use gen1int_matrix
    integer, intent(in) :: num_shells_bra
    type(sub_shell_t), intent(in) :: sub_shells_bra(*)
    integer, intent(in) :: num_shells_ket
    type(sub_shell_t), intent(in) :: sub_shells_ket(*)
    logical, intent(in) :: same_braket
    type(one_prop_t), intent(in) :: one_prop
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer, optional, intent(in) :: api_comm
    integer, intent(in) :: num_prop
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_geo_ket
    integer, intent(in) :: num_geo_total
    type(matrix), optional, intent(inout) :: val_ints(num_prop,    &
                                                      num_geo_bra, &
                                                      num_geo_ket, &
                                                      num_geo_total)
    logical, optional, intent(in) :: write_ints
    integer, intent(in) :: num_dens
    type(matrix), optional, intent(in) :: ao_dens(num_dens)
    real(REALK), optional, target, intent(inout) :: val_expt(num_prop,      &
                                                             num_geo_bra,   &
                                                             num_geo_ket,   &
                                                             num_geo_total, &
                                                             num_dens)
    logical, optional, intent(in) :: write_expt
    logical spher_gto                               !if using spherical GTOs
    integer prop_sym                                !symmetry of property integrals
    integer num_cent_bra                            !number of differentiated centers on bra center
    integer num_cent_ket                            !number of differentiated centers on ket center
    integer idx_cent_bra(1)                         !indices of differentiated centers on bra center
    integer idx_cent_ket(1)                         !indices of differentiated centers on ket center
    integer order_cent_bra(1)                       !orders of differentiated centers on bra center
    integer order_cent_ket(1)                       !orders of differentiated centers on ket center
    integer min_shell_bra, max_shell_bra            !minimum and maximum of sub-shells on bra center for calculations
    integer min_shell_ket, max_shell_ket            !minimum and maximum of sub-shells on ket center for calculations
    integer max_ao_bra                              !maximum number of AOs in a sub-shell on bra center
    integer max_ao_ket                              !maximum number of AOs in a sub-shell on ket center
    integer size_ao                                 !size of AOs
    integer num_pairs                               !number of AO sub-shell pairs to calculate
    integer path_num_bgeo                           !number of geometric derivatives on bra center
    integer path_num_kgeo                           !number of geometric derivatives on ket center
    integer path_num_tgeo                           !number of unique total geometric derivatives of current path
    integer path_offset_bgeo                        !offset of geometric derivatives on bra center
    integer path_offset_kgeo                        !offset of geometric derivatives on ket center
    integer path_offset_tgeo                        !offset of unique total geometric derivatives of current path
    integer num_matrices                            !number of integral matrices
    integer size_ints                               !size of contracted integrals
    real(REALK), allocatable :: contr_ints(:)       !contracted integrals between two AO sub-shells (pair)
    logical p_write_ints                            !if writing integral matrices on file
    logical p_write_expt                            !if writing expectation values on file
    logical do_integral                             !if returning and/or writing integral matrices
    logical do_expectation                          !if calculating or writing expectaion values on file
    real(REALK), pointer :: unique_expt(:,:,:,:,:)  !expectation values with unique geometric derivatives
    integer remaining_jobs                          !number of remaining jobs
    integer shell_pair(2)                           !indices of AO sub-shell pair on bra and ket centers
    integer min_row_idx, max_row_idx                !minimum and maximum indices of rows (on bra center)
    integer min_col_idx, max_col_idx                !minimum and maximum indices of columns (on ket center)
    integer start_ints, end_ints                    !start and end addresses of contracted integrals
    integer ishell, idens, bgeo, kgeo, tgeo, iprop  !incremental recorders
    integer ierr                                    !error information
#if defined(VAR_MPI)
#include "mpif.h"
    integer rank_proc                               !rank of processor
    integer num_proc                                !number of processors
    integer worker_request(3)                       !request from a worker, in which the first two elements
                                                    !are either \var(REQUEST_WORK) or the AO sub-shell pair
                                                    !to send back, the third is the rank of the worker
    integer msg_tag                                 !message tag
    integer mpi_status(MPI_STATUS_SIZE)             !MPI status
#endif
    ! since we do not use mixed CGTOs and SGTOs, we get this information only from
    ! the first sub-shell on bra center
    spher_gto = sub_shells_bra(1)%spher_gto
    ! gets the number of differentiated centers on bra and ket centers
    call NaryTreePathGetNumCenters(nary_tree=nary_tree_bra, num_centers=num_cent_bra)
    call NaryTreePathGetNumCenters(nary_tree=nary_tree_ket, num_centers=num_cent_ket)
    if (num_cent_bra>1 .or. num_cent_ket>1) then
      return
    ! one-center partial geometric derivatives
    else if (num_cent_bra==1 .and. num_cent_ket==1) then
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_bra, idx_cent=idx_cent_bra)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_bra, order_cent=order_cent_bra)
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_ket, idx_cent=idx_cent_ket)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_ket, order_cent=order_cent_ket)
      ! finds the range of shells for calculations
      min_shell_bra = num_shells_bra+1
      max_shell_bra = 0
      do ishell = 1, num_shells_bra
        if (sub_shells_bra(ishell)%idx_cent==idx_cent_bra(1)) then
          if (ishell<min_shell_bra) min_shell_bra = ishell
          if (ishell>max_shell_bra) max_shell_bra = ishell
        end if
      end do
      if (min_shell_bra>max_shell_bra) return
      min_shell_ket = num_shells_ket+1
      max_shell_ket = 0
      do ishell = 1, num_shells_ket
        if (sub_shells_ket(ishell)%idx_cent==idx_cent_ket(1)) then
          if (ishell<min_shell_ket) min_shell_ket = ishell
          if (ishell>max_shell_ket) max_shell_ket = ishell
        end if
      end do
      if (min_shell_ket>max_shell_ket) return
    ! one-center partial geometric derivatives on bra center
    else if (num_cent_bra==1 .and. num_cent_ket==0) then
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_bra, idx_cent=idx_cent_bra)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_bra, order_cent=order_cent_bra)
      order_cent_ket = 0
      ! finds the range of shells for calculations
      min_shell_bra = num_shells_bra+1
      max_shell_bra = 0
      do ishell = 1, num_shells_bra
        if (sub_shells_bra(ishell)%idx_cent==idx_cent_bra(1)) then
          if (ishell<min_shell_bra) min_shell_bra = ishell
          if (ishell>max_shell_bra) max_shell_bra = ishell
        end if
      end do
      if (min_shell_bra>max_shell_bra) return
      min_shell_ket = 1
      max_shell_ket = num_shells_ket
    ! one-center partial geometric derivatives on ket center
    else if (num_cent_bra==0 .and. num_cent_ket==1) then
      order_cent_bra = 0
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_ket, idx_cent=idx_cent_ket)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_ket, order_cent=order_cent_ket)
      ! finds the range of shells for calculations
      min_shell_bra = 1
      max_shell_bra = num_shells_bra
      min_shell_ket = num_shells_ket+1
      max_shell_ket = 0
      do ishell = 1, num_shells_ket
        if (sub_shells_ket(ishell)%idx_cent==idx_cent_ket(1)) then
          if (ishell<min_shell_ket) min_shell_ket = ishell
          if (ishell>max_shell_ket) max_shell_ket = ishell
        end if
      end do
      if (min_shell_ket>max_shell_ket) return
    ! no partial geometric derivatives
    else
      order_cent_bra = 0
      order_cent_ket = 0
      min_shell_bra = 1
      max_shell_bra = num_shells_bra
      min_shell_ket = 1
      max_shell_ket = num_shells_ket
    end if
    ! sets the maximum number of AOs in a sub-shell on bra center
    max_ao_bra = 0
    do ishell = min_shell_bra, max_shell_bra
      size_ao = sub_shells_bra(ishell)%num_ao*sub_shells_bra(ishell)%num_contr
      if (max_ao_bra<size_ao) max_ao_bra = size_ao
    end do
    ! sets the number of AO sub-shell pairs to calculate
    if (same_braket .and. num_shells_bra==num_shells_ket .and. &
        num_cent_bra==0 .and. num_cent_ket==0) then
      ! gets the symmetry of property integrals
      call OnePropGetSymmetry(one_prop=one_prop, prop_sym=prop_sym)
      select case (prop_sym)
      case (SYMM_INT_MAT,ANTI_INT_MAT)
        num_pairs = num_shells_bra*(num_shells_bra+1)/2
      case default
        num_pairs = num_shells_bra*num_shells_ket
      end select
      ! sets the maximum number of AOs in a sub-shell on ket center
      max_ao_ket = max_ao_bra
    else
      ! we need to calculate all the sub-shell pairs even for symmetric and anti-symmetric integrals
      prop_sym = SQUARE_INT_MAT
      num_pairs = (max_shell_bra-min_shell_bra+1)*(max_shell_ket-min_shell_ket+1)
      ! sets the maximum number of AOs in a sub-shell on ket center
      max_ao_ket = 0
      do ishell = min_shell_ket, max_shell_ket
        size_ao = sub_shells_ket(ishell)%num_ao*sub_shells_ket(ishell)%num_contr
        if (max_ao_ket<size_ao) max_ao_ket = size_ao
      end do
    end if
    ! gets the numbers of unique geometric derivatives of current path
    path_num_bgeo = (order_cent_bra(1)+1)*(order_cent_bra(1)+2)/2
    path_num_kgeo = (order_cent_ket(1)+1)*(order_cent_ket(1)+2)/2
    call NaryTreePathGetNumUnique(nary_tree=nary_tree_total, path_num_unique=path_num_tgeo)
    ! gets the offsets of geometric derivatives of current path
    call NaryTreePathGetOffset(nary_tree=nary_tree_bra, path_offset=path_offset_bgeo)
    call NaryTreePathGetOffset(nary_tree=nary_tree_ket, path_offset=path_offset_kgeo)
    call NaryTreePathGetOffset(nary_tree=nary_tree_total, path_offset=path_offset_tgeo)
    ! sets the number of integral matrices
    num_matrices = num_prop*path_num_bgeo*path_num_kgeo*path_num_tgeo
    ! allocates cache of contracted integrals between two AO sub-shells
    size_ints = max_ao_bra*max_ao_ket*num_matrices
    allocate(contr_ints(size_ints), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntShellGetIntExpt>> failed to allocate contr_ints!")
    end if
#if defined(VAR_MPI)
    ! gets the rank of this processor and the number of processors
    if (present(api_comm)) then
      call MPI_Comm_rank(api_comm, rank_proc, ierr)
      call MPI_Comm_size(api_comm, num_proc, ierr)
    else
      rank_proc = MANAGER
      num_proc = 1
    end if
    ! manager processor has arguments \var(val_ints) and/or \var(val_expt)
    if (rank_proc==MANAGER) then
#endif
      ! if writing integral matrices on file
      if (present(write_ints)) then
        p_write_ints = write_ints
      else
        p_write_ints = .false.
      end if
      ! if calculating contracted integrals (only needed for parallel mode)
      do_integral = present(val_ints) .or. p_write_ints
      ! if writing expectation values on file
      if (present(write_expt)) then
        p_write_expt = write_expt
      else
        p_write_expt = .false.
      end if
      ! if calculating expectation values
      do_expectation = present(ao_dens) .and. (present(val_expt) .or. p_write_expt)
      ! sets the array of expectation values
      if (do_expectation) then
        if (.not. present(val_expt)) then
          allocate(unique_expt(num_prop,                                          &
                               path_offset_bgeo+1:path_offset_bgeo+path_num_bgeo, &
                               path_offset_kgeo+1:path_offset_kgeo+path_num_kgeo, &
                               path_offset_tgeo+1:path_offset_tgeo+path_num_tgeo, &
                               num_dens), stat=ierr)
          if (ierr/=0) then
            call quit("Gen1IntShellGetIntExpt>> failed to allocate unique_expt!")
          end if
          unique_expt = 0.0_REALK
        else
          unique_expt => val_expt
        end if
      end if
#if defined(VAR_MPI)
      ! broadcasts what kind of jobs to do
      if (present(api_comm)) then
        call MPI_Bcast(do_integral, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
        call MPI_Bcast(do_expectation, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
      end if
    ! worker processors do not have arguments \var(val_ints) and/or \var(val_expt)
    else
      call MPI_Bcast(do_integral, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
      call MPI_Bcast(do_expectation, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
      if (do_expectation) then
        if (.not.present(ao_dens)) then
          call quit("Gen1IntShellGetIntExpt>> worker processor does not have ao_dens!")
        end if
        allocate(unique_expt(num_prop,                                          &
                             path_offset_bgeo+1:path_offset_bgeo+path_num_bgeo, &
                             path_offset_kgeo+1:path_offset_kgeo+path_num_kgeo, &
                             path_offset_tgeo+1:path_offset_tgeo+path_num_tgeo, &
                             num_dens), stat=ierr)
        if (ierr/=0) then
          call quit("Gen1IntShellGetIntExpt>> failed to allocate unique_expt on worker processor!")
        end if
        unique_expt = 0.0_REALK
      end if
    end if
    ! calculations with more than one processor
    if (num_proc>1) then
      ! sends the message tag
      msg_tag = 1
      ! manager code
      if (rank_proc==MANAGER) then
        ! initializes the number of remaining jobs (the manager needs to send "finish" signal
        ! to other worker processors)
        remaining_jobs = num_proc-1
        ! intializes the indices of AO sub-shell pair on bra and ket centers
        shell_pair(1) = min_shell_bra-1
        shell_pair(2) = min_shell_ket
        do while (remaining_jobs>0)
          ! receives a request from a woker
          call MPI_Recv(worker_request, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, api_comm, mpi_status, ierr)
          ! the worker requests new work
          if (worker_request(1)==REQUEST_WORK) then
            ! no more sub-shell pair to calculate
            if (shell_pair(1)>=max_shell_bra .and. shell_pair(2)>=max_shell_ket) then
              call MPI_Send((/NO_MORE_WORK,NO_MORE_WORK/), 2, MPI_INTEGER, &
                            worker_request(3), msg_tag, api_comm, ierr)
              ! decreases the number of remaining jobs
              remaining_jobs = remaining_jobs-1
            else
              ! prepares the next sub-shell pair to calculate
              select case (prop_sym)
              ! for symmetric and anti-symmetric matrices, we only calculate the upper
              ! and diagonal parts explicitly
              case (SYMM_INT_MAT,ANTI_INT_MAT)
                if (shell_pair(1)<shell_pair(2)) then
                  shell_pair(1) = shell_pair(1)+1
                else
                  shell_pair(2) = shell_pair(2)+1
                  shell_pair(1) = 1
                end if
              case default
                if (shell_pair(1)<max_shell_bra) then
                  shell_pair(1) = shell_pair(1)+1
                else
                  shell_pair(2) = shell_pair(2)+1
                  shell_pair(1) = min_shell_bra
                end if
              end select
              ! sends the next sub-shell pair to the worker
              call MPI_Send(shell_pair, 2, MPI_INTEGER, worker_request(3), &
                            msg_tag, api_comm, ierr)
            end if
          ! the worker wants to send the contracted integrals back
          else
            size_ao = sub_shells_bra(worker_request(1))%num_ao    &
                    * sub_shells_bra(worker_request(1))%num_contr &
                    * sub_shells_ket(worker_request(2))%num_ao    &
                    * sub_shells_ket(worker_request(2))%num_contr
            size_ints =  size_ao*num_matrices
            ! receives results from the worker
            call MPI_Recv(contr_ints(1:size_ints), size_ints, MPI_REALK, &
                          worker_request(3), MPI_ANY_TAG, api_comm,      &
                          mpi_status, ierr)
            ! sets the minimum and maximum of indices of rows of the integral matrices
            min_row_idx = sub_shells_bra(worker_request(1))%base_idx+1
            max_row_idx = sub_shells_bra(worker_request(1))%base_idx &
                        + sub_shells_bra(worker_request(1))%num_ao   &
                        * sub_shells_bra(worker_request(1))%num_contr
            ! sets the minimum and maximum of indices of columns of the integral matrices
            min_col_idx = sub_shells_ket(worker_request(2))%base_idx+1
            max_col_idx = sub_shells_ket(worker_request(2))%base_idx &
                        + sub_shells_ket(worker_request(2))%num_ao   &
                        * sub_shells_ket(worker_request(2))%num_contr
            ! assigns the returned integrals
            if (present(val_ints)) then
              start_ints = 0
              do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
                do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
                  do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                    do iprop = 1, num_prop
                      end_ints = start_ints+size_ao
                      start_ints = start_ints+1
                      call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                        min_row_idx, max_row_idx,       &
                                        min_col_idx, max_col_idx,       &
                                        contr_ints(start_ints:end_ints), .false.)
                      if (prop_sym==SYMM_INT_MAT .and. &
                          worker_request(1)/=worker_request(2)) then
                        call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                          min_row_idx, max_row_idx,       &
                                          min_col_idx, max_col_idx,       &
                                          contr_ints(start_ints:end_ints), .true.)
                      else if (prop_sym==ANTI_INT_MAT .and. &
                               worker_request(1)/=worker_request(2)) then
                        call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                          min_row_idx, max_row_idx,       &
                                          min_col_idx, max_col_idx,       &
                                          -contr_ints(start_ints:end_ints), .true.)
                      end if
                      start_ints = end_ints
                    end do
                  end do
                end do
              end do
            end if
!FIXME
            ! writes the integrals on file
            if (p_write_ints) then
            end if
          end if
        end do
        ! receives expectation values from worker processors
        if (do_expectation) then
          call MPI_Reduce(MPI_IN_PLACE,                                                     &
                          unique_expt(:,path_offset_bgeo+1:path_offset_bgeo+path_num_bgeo,  &
                                      path_offset_kgeo+1:path_offset_kgeo+path_num_kgeo,    &
                                      path_offset_tgeo+1:path_offset_tgeo+path_num_tgeo,:), &
                          num_matrices*num_dens, MPI_REALK, MPI_SUM, MANAGER, api_comm, ierr)
        end if
      ! worker code
      else
        worker_request(3) = rank_proc
        do while (.true.)
          ! sends request for a new work to manager
          worker_request(1) = REQUEST_WORK
          call MPI_Send(worker_request, 3, MPI_INTEGER, MANAGER, &
                        msg_tag, api_comm, ierr)
          ! receives the next sub-shell pair or "finished" signal from manager
          call MPI_Recv(shell_pair, 2, MPI_INTEGER, MANAGER, MPI_ANY_TAG, &
                        api_comm, mpi_status, ierr)
          if (shell_pair(1)==NO_MORE_WORK) then
            exit
          else
            ! calculates the contracted integrals
            size_ao = sub_shells_bra(shell_pair(1))%num_ao    &
                    * sub_shells_bra(shell_pair(1))%num_contr &
                    * sub_shells_ket(shell_pair(2))%num_ao    &
                    * sub_shells_ket(shell_pair(2))%num_contr
            size_ints =  size_ao*num_matrices
            ! spherical GTOs
            if (spher_gto) then
              ! calls Gen1Int subroutines to evaluate property integrals
              call OnePropGetIntegral(                                        &
                     idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                     coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                     angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                     num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                     exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                     num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                     contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                     idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                     coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                     angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                     num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                     exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                     num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                     contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                     spher_gto=spher_gto,                                     &
                     one_prop=one_prop,                                       &
                     order_geo_bra=order_cent_bra(1),                         &
                     order_geo_ket=order_cent_ket(1),                         &
                     nary_tree_total=nary_tree_total,                         &
                     num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                     num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                     num_opt=num_matrices, contr_ints=contr_ints)
              ! reorders p-shell spherical GTOs, since Dalton uses x(+1), y(-1) and z(0),
              ! while Gen1Int uses y(-1), z(0) and x(+1)
              if (sub_shells_bra(shell_pair(1))%ang_num==1)                              &
                call gen1int_reorder_p_sgto(1, sub_shells_bra(shell_pair(1))%num_contr   &
                                               *sub_shells_ket(shell_pair(2))%num_ao     &
                                               *sub_shells_ket(shell_pair(2))%num_contr, &
                                            num_matrices, contr_ints)
              if (sub_shells_ket(shell_pair(2))%ang_num==1)                           &
                call gen1int_reorder_p_sgto(sub_shells_bra(shell_pair(1))%num_ao      &
                                            *sub_shells_bra(shell_pair(1))%num_contr, &
                                            sub_shells_ket(shell_pair(2))%num_contr,  &
                                            num_matrices, contr_ints)
            ! Cartesian GTOs
            else
              ! calls Gen1Int subroutines to evaluate property integrals, and reorders
              ! the integrals according to Cartesian powers
              call OnePropGetIntegral(                                        &
                     idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                     coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                     angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                     num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                     exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                     num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                     contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                     idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                     coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                     angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                     num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                     exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                     num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                     contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                     spher_gto=spher_gto,                                     &
                     one_prop=one_prop,                                       &
                     order_geo_bra=order_cent_bra(1),                         &
                     order_geo_ket=order_cent_ket(1),                         &
                     nary_tree_total=nary_tree_total,                         &
                     num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                     num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                     num_opt=num_matrices, contr_ints=contr_ints,             &
                     powers_bra=sub_shells_bra(shell_pair(1))%powers,         &
                     powers_ket=sub_shells_ket(shell_pair(2))%powers)
            end if
            ! sends the contracted integrals to manager
            if (do_integral) then
              worker_request(1:2) = shell_pair
              call MPI_Send(worker_request, 3, MPI_INTEGER, MANAGER, &
                            msg_tag, api_comm, ierr)
              call MPI_Send(contr_ints(1:size_ints), size_ints, MPI_REALK, &
                            MANAGER, msg_tag, api_comm, ierr)
            end if
            ! calculates the expectation values with unique total geometric derivatives
            if (do_expectation) then
              ! sets the minimum and maximum of indices of rows of the integral matrices
              min_row_idx = sub_shells_bra(shell_pair(1))%base_idx+1
              max_row_idx = sub_shells_bra(shell_pair(1))%base_idx &
                          + sub_shells_bra(shell_pair(1))%num_ao   &
                          * sub_shells_bra(shell_pair(1))%num_contr
              ! sets the minimum and maximum of indices of columns of the integral matrices
              min_col_idx = sub_shells_ket(shell_pair(2))%base_idx+1
              max_col_idx = sub_shells_ket(shell_pair(2))%base_idx &
                          + sub_shells_ket(shell_pair(2))%num_ao   &
                          * sub_shells_ket(shell_pair(2))%num_contr
              do idens = 1, num_dens
                start_ints = 0
                do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
                  do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
                    do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                      do iprop = 1, num_prop
                        end_ints = start_ints+size_ao
                        start_ints = start_ints+1
                        call MatMultBlockedTrace(ao_dens(idens),                          &
                                                 min_row_idx, max_row_idx,                &
                                                 min_col_idx, max_col_idx,                &
                                                 contr_ints(start_ints:end_ints),         &
                                                 unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                                 .false.)
                        if (prop_sym==SYMM_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                          call MatMultBlockedTrace(ao_dens(idens),                          &
                                                   min_row_idx, max_row_idx,                &
                                                   min_col_idx, max_col_idx,                &
                                                   contr_ints(start_ints:end_ints),         &
                                                   unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                                   .true.)
                        else if (prop_sym==ANTI_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                          call MatMultBlockedTrace(ao_dens(idens),                          &
                                                   min_row_idx, max_row_idx,                &
                                                   min_col_idx, max_col_idx,                &
                                                   -contr_ints(start_ints:end_ints),        &
                                                   unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                                   .true.)
                        end if
                        start_ints = end_ints
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end if
        end do
        ! sends expectation values to manager processor
        if (do_expectation) then
          call MPI_Reduce(unique_expt, unique_expt, num_matrices*num_dens, &
                          MPI_REALK, MPI_SUM, MANAGER, api_comm, ierr)
        end if
      end if
    ! calculations with only one processor
    else
#endif
      ! initializes the number of remaining pairs to calculate
      remaining_jobs = num_pairs
      ! intializes the indices of AO sub-shell pair on bra and ket centers
      shell_pair(1) = min_shell_bra-1
      shell_pair(2) = min_shell_ket
      do while (remaining_jobs>0)
        ! prepares the next sub-shell pair to calculate
        select case (prop_sym)
        ! for symmetric and anti-symmetric matrices, we only calculate the upper
        ! and diagonal parts explicitly
        case (SYMM_INT_MAT,ANTI_INT_MAT)
          if (shell_pair(1)<shell_pair(2)) then
            shell_pair(1) = shell_pair(1)+1
          else
            shell_pair(2) = shell_pair(2)+1
            shell_pair(1) = 1
          end if
        case default
          if (shell_pair(1)<max_shell_bra) then
            shell_pair(1) = shell_pair(1)+1
          else
            shell_pair(2) = shell_pair(2)+1
            shell_pair(1) = min_shell_bra
          end if
        end select
        size_ao = sub_shells_bra(shell_pair(1))%num_ao    &
                * sub_shells_bra(shell_pair(1))%num_contr &
                * sub_shells_ket(shell_pair(2))%num_ao    &
                * sub_shells_ket(shell_pair(2))%num_contr
        size_ints =  size_ao*num_matrices
        ! spherical GTOs
        if (spher_gto) then
          ! calls Gen1Int subroutines to evaluate property integrals
          call OnePropGetIntegral(                                        &
                 idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                 coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                 angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                 num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                 exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                 num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                 contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                 idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                 coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                 angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                 num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                 exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                 num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                 contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                 spher_gto=spher_gto,                                     &
                 one_prop=one_prop,                                       &
                 order_geo_bra=order_cent_bra(1),                         &
                 order_geo_ket=order_cent_ket(1),                         &
                 nary_tree_total=nary_tree_total,                         &
                 num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                 num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                 num_opt=num_matrices, contr_ints=contr_ints)
          ! reorders p-shell spherical GTOs, since Dalton uses x(+1), y(-1) and z(0),
          ! while Gen1Int uses y(-1), z(0) and x(+1)
          if (sub_shells_bra(shell_pair(1))%ang_num==1)                              &
            call gen1int_reorder_p_sgto(1, sub_shells_bra(shell_pair(1))%num_contr   &
                                           *sub_shells_ket(shell_pair(2))%num_ao     &
                                           *sub_shells_ket(shell_pair(2))%num_contr, &
                                        num_matrices, contr_ints)
          if (sub_shells_ket(shell_pair(2))%ang_num==1)                           &
            call gen1int_reorder_p_sgto(sub_shells_bra(shell_pair(1))%num_ao      &
                                        *sub_shells_bra(shell_pair(1))%num_contr, &
                                        sub_shells_ket(shell_pair(2))%num_contr,  &
                                        num_matrices, contr_ints)
        ! Cartesian GTOs
        else
          ! calls Gen1Int subroutines to evaluate property integrals, and reorders
          ! the integrals according to Cartesian powers
          call OnePropGetIntegral(                                        &
                 idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                 coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                 angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                 num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                 exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                 num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                 contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                 idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                 coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                 angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                 num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                 exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                 num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                 contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                 spher_gto=spher_gto,                                     &
                 one_prop=one_prop,                                       &
                 order_geo_bra=order_cent_bra(1),                         &
                 order_geo_ket=order_cent_ket(1),                         &
                 nary_tree_total=nary_tree_total,                         &
                 num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                 num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                 num_opt=num_matrices, contr_ints=contr_ints,             &
                 powers_bra=sub_shells_bra(shell_pair(1))%powers,         &
                 powers_ket=sub_shells_ket(shell_pair(2))%powers)
        end if
        ! sets the minimum and maximum of indices of rows of the integral matrices
        min_row_idx = sub_shells_bra(shell_pair(1))%base_idx+1
        max_row_idx = sub_shells_bra(shell_pair(1))%base_idx &
                    + sub_shells_bra(shell_pair(1))%num_ao   &
                    * sub_shells_bra(shell_pair(1))%num_contr
        ! sets the minimum and maximum of indices of columns of the integral matrices
        min_col_idx = sub_shells_ket(shell_pair(2))%base_idx+1
        max_col_idx = sub_shells_ket(shell_pair(2))%base_idx &
                    + sub_shells_ket(shell_pair(2))%num_ao   &
                    * sub_shells_ket(shell_pair(2))%num_contr
        ! assigns the returned integrals
        if (present(val_ints)) then
          start_ints = 0
          do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
            do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
              do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                do iprop = 1, num_prop
                  end_ints = start_ints+size_ao
                  start_ints = start_ints+1
                  call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                    min_row_idx, max_row_idx,       &
                                    min_col_idx, max_col_idx,       &
                                    contr_ints(start_ints:end_ints), .false.)
                  if (prop_sym==SYMM_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                    call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                      min_row_idx, max_row_idx,       &
                                      min_col_idx, max_col_idx,       &
                                      contr_ints(start_ints:end_ints), .true.)
                  else if (prop_sym==ANTI_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                    call MatSetValues(val_ints(iprop,bgeo,kgeo,tgeo), &
                                      min_row_idx, max_row_idx,       &
                                      min_col_idx, max_col_idx,       &
                                      -contr_ints(start_ints:end_ints), .true.)
                  end if
                  start_ints = end_ints
                end do
              end do
            end do
          end do
        end if
!FIXME
        ! writes the integrals on file
        if (p_write_ints) then
        end if
        ! calculates the expectation values with unique total geometric derivatives
        if (do_expectation) then
          do idens = 1, num_dens
            start_ints = 0
            do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
              do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
                do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                  do iprop = 1, num_prop
                    end_ints = start_ints+size_ao
                    start_ints = start_ints+1
                    call MatMultBlockedTrace(ao_dens(idens),                          &
                                             min_row_idx, max_row_idx,                &
                                             min_col_idx, max_col_idx,                &
                                             contr_ints(start_ints:end_ints),         &
                                             unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                             .false.)
                    if (prop_sym==SYMM_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                      call MatMultBlockedTrace(ao_dens(idens),                          &
                                               min_row_idx, max_row_idx,                &
                                               min_col_idx, max_col_idx,                &
                                               contr_ints(start_ints:end_ints),         &
                                               unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                               .true.)
                    else if (prop_sym==ANTI_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                      call MatMultBlockedTrace(ao_dens(idens),                          &
                                               min_row_idx, max_row_idx,                &
                                               min_col_idx, max_col_idx,                &
                                               -contr_ints(start_ints:end_ints),        &
                                               unique_expt(iprop,bgeo,kgeo,tgeo,idens), &
                                               .true.)
                    end if
                    start_ints = end_ints
                  end do
                end do
              end do
            end do
          end do
        end if 
        ! decreases the number of remaining pairs
        remaining_jobs = remaining_jobs-1
      end do
#if defined(VAR_MPI)
    end if
#endif
    ! frees spaces
    deallocate(contr_ints)
    if (do_expectation) then
!FIXME
      ! writes expectation values \var(unique_expt) on file
#if defined(VAR_MPI)
      if (rank_proc==MANAGER) then
#endif
        if (p_write_expt) then
        end if
#if defined(VAR_MPI)
      end if
#endif
      if (.not.present(val_expt)) then
        deallocate(unique_expt)
      end if
      nullify(unique_expt)
    end if
#if defined(VAR_MPI)
    ! blocks until all processors have finished
    if (present(api_comm)) call MPI_Barrier(api_comm, ierr)
#endif
  end subroutine Gen1IntShellGetIntExpt

  !> \brief calculates property integrand matrices and/or expectation values for
  !>        given AO sub-shells on bra and ket centers
  !> \author Bin Gao and Radovan Bast
  !> \date 2012-05-15
  !> \param num_shells_bra is the number of AO sub-shells on bra center
  !> \param sub_shells_bra are the AO sub-shells on bra center
  !> \param num_shells_ket is the number of AO sub-shells on ket center
  !> \param sub_shells_ket are the AO sub-shells on ket center
  !> \param same_braket indicates if the AO sub-shells are the same on bra and ket centers
  !> \param one_prop contains the information of one-electron property integrals
  !> \param nary_tree_bra is the N-ary tree for geometric derivatives on bra center
  !> \param nary_tree_ket is the N-ary tree for geometric derivatives on ket center
  !> \param nary_tree_total is the N-ary tree for total geometric derivatives
  !> \param api_comm is the MPI communicator
  !> \param num_points is the number of grid points
  !> \param grid_points contains the coordinates of grid points
  !> \param num_dens is the number of AO density matrices
  !> \param ao_dens contains the AO density matrices
  !> \param num_prop is the number of property integrals including various derivatives (except
  !>        geometric derivatives)
  !> \param num_geo_bra is number of geometric derivatives on bra center
  !> \param num_geo_ket is number of geometric derivatives on ket center
  !> \param num_geo_total is number of total geometric derivatives
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \return val_expt contains the one-electron property intergrands contracted with AO density matrices
  !> \note the arrangement of \var(val_expt) will be in the order of
  !>       \var(order_mom), \var(order_mag_bra), ..., \var(order_geo_total), and each of
  !>       them is arranged in the order of (xx,xy,yy,xz,yz,zz) or (xx,yx,zx,xy,yy,zy,xz,yz,zz),
  !>       see Gen1Int library manual, for instance Section 2.2;
  !>       \var(val_expt) should be zero by users before calculations
  subroutine Gen1IntShellGetFunExpt(num_shells_bra, sub_shells_bra, &
                                    num_shells_ket, sub_shells_ket, &
                                    same_braket, one_prop,          &
                                    nary_tree_bra, nary_tree_ket,   &
                                    nary_tree_total, api_comm,      &
                                    num_points, grid_points,        &
                                    num_dens, ao_dens, num_prop,    &
                                    num_geo_bra, num_geo_ket,       &
                                    num_geo_total, val_expt)
    ! matrix module
    use gen1int_matrix
    integer, intent(in) :: num_shells_bra
    type(sub_shell_t), intent(in) :: sub_shells_bra(*)
    integer, intent(in) :: num_shells_ket
    type(sub_shell_t), intent(in) :: sub_shells_ket(*)
    logical, intent(in) :: same_braket
    type(one_prop_t), intent(in) :: one_prop
    type(nary_tree_t), intent(inout) :: nary_tree_bra
    type(nary_tree_t), intent(inout) :: nary_tree_ket
    type(nary_tree_t), intent(inout) :: nary_tree_total
    integer, optional, intent(in) :: api_comm
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_dens
    type(matrix), intent(in) :: ao_dens(num_dens)
    integer, intent(in) :: num_prop
    integer, intent(in) :: num_geo_bra
    integer, intent(in) :: num_geo_ket
    integer, intent(in) :: num_geo_total
    real(REALK), optional, target, intent(inout) :: val_expt(num_points,    &
                                                             num_prop,      &
                                                             num_geo_bra,   &
                                                             num_geo_ket,   &
                                                             num_geo_total, &
                                                             num_dens)
    logical spher_gto                                 !if using spherical GTOs
    integer num_cent_bra                              !number of differentiated centers on bra center
    integer num_cent_ket                              !number of differentiated centers on ket center
    integer idx_cent_bra(1)                           !indices of differentiated centers on bra center
    integer idx_cent_ket(1)                           !indices of differentiated centers on ket center
    integer order_cent_bra(1)                         !orders of differentiated centers on bra center
    integer order_cent_ket(1)                         !orders of differentiated centers on ket center
    integer min_shell_bra, max_shell_bra              !minimum and maximum of sub-shells on bra center for calculations
    integer min_shell_ket, max_shell_ket              !minimum and maximum of sub-shells on ket center for calculations
    integer prop_sym                                  !symmetry of property integrals
    integer max_ao_bra                                !maximum number of AOs in a sub-shell on bra center
    integer max_ao_ket                                !maximum number of AOs in a sub-shell on ket center
    integer size_ao                                   !size of AOs
    integer num_pairs                                 !number of AO sub-shell pairs to calculate
    integer path_num_bgeo                             !number of geometric derivatives on bra center
    integer path_num_kgeo                             !number of geometric derivatives on ket center
    integer path_num_tgeo                             !number of unique total geometric derivatives of current path
    integer path_offset_bgeo                          !offset of geometric derivatives on bra center
    integer path_offset_kgeo                          !offset of geometric derivatives on ket center
    integer path_offset_tgeo                          !offset of unique total geometric derivatives of current path
    integer num_matrices                              !number of integral matrices
    integer size_ints                                 !size of contracted integrals
    real(REALK), allocatable :: contr_ints(:)         !contracted integrals between two AO sub-shells (pair)
    logical do_expectation                            !if calculating or writing expectaion values on file
    real(REALK), pointer :: unique_expt(:,:,:,:,:,:)  !expectation values with unique geometric derivatives
    integer remaining_jobs                            !number of remaining jobs
    integer shell_pair(2)                             !indices of AO sub-shell pair on bra and ket centers
    integer min_row_idx, max_row_idx                  !minimum and maximum indices of rows (on bra center)
    integer min_col_idx, max_col_idx                  !minimum and maximum indices of columns (on ket center)
    integer start_ints, end_ints                      !start and end addresses of contracted integrals
    integer ishell, idens, bgeo, kgeo, tgeo, &        !incremental recorders
            iprop, ipoint
    integer ierr                                      !error information
#if defined(VAR_MPI)
#include "mpif.h"
    integer rank_proc                                 !rank of processor
    integer num_proc                                  !number of processors
    integer worker_request(3)                         !request from a worker, in which the first two elements
                                                      !are either \var(REQUEST_WORK) or the AO sub-shell pair
                                                      !to send back, the third is the rank of the worker
    integer msg_tag                                   !message tag
    integer mpi_status(MPI_STATUS_SIZE)               !MPI status
#endif
    ! since we do not use mixed CGTOs and SGTOs, we get this information only from
    ! the first sub-shell on bra center
    spher_gto = sub_shells_bra(1)%spher_gto
    ! gets the number of differentiated centers on bra and ket centers
    call NaryTreePathGetNumCenters(nary_tree=nary_tree_bra, num_centers=num_cent_bra)
    call NaryTreePathGetNumCenters(nary_tree=nary_tree_ket, num_centers=num_cent_ket)
    if (num_cent_bra>1 .or. num_cent_ket>1) then
      return
    ! one-center partial geometric derivatives
    else if (num_cent_bra==1 .and. num_cent_ket==1) then
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_bra, idx_cent=idx_cent_bra)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_bra, order_cent=order_cent_bra)
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_ket, idx_cent=idx_cent_ket)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_ket, order_cent=order_cent_ket)
      ! finds the range of shells for calculations
      min_shell_bra = num_shells_bra+1
      max_shell_bra = 0
      do ishell = 1, num_shells_bra
        if (sub_shells_bra(ishell)%idx_cent==idx_cent_bra(1)) then
          if (ishell<min_shell_bra) min_shell_bra = ishell
          if (ishell>max_shell_bra) max_shell_bra = ishell
        end if
      end do
      if (min_shell_bra>max_shell_bra) return
      min_shell_ket = num_shells_ket+1
      max_shell_ket = 0
      do ishell = 1, num_shells_ket
        if (sub_shells_ket(ishell)%idx_cent==idx_cent_ket(1)) then
          if (ishell<min_shell_ket) min_shell_ket = ishell
          if (ishell>max_shell_ket) max_shell_ket = ishell
        end if
      end do
      if (min_shell_ket>max_shell_ket) return
    ! one-center partial geometric derivatives on bra center
    else if (num_cent_bra==1 .and. num_cent_ket==0) then
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_bra, idx_cent=idx_cent_bra)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_bra, order_cent=order_cent_bra)
      order_cent_ket = 0
      ! finds the range of shells for calculations
      min_shell_bra = num_shells_bra+1
      max_shell_bra = 0
      do ishell = 1, num_shells_bra
        if (sub_shells_bra(ishell)%idx_cent==idx_cent_bra(1)) then
          if (ishell<min_shell_bra) min_shell_bra = ishell
          if (ishell>max_shell_bra) max_shell_bra = ishell
        end if
      end do
      if (min_shell_bra>max_shell_bra) return
      min_shell_ket = 1
      max_shell_ket = num_shells_ket
    ! one-center partial geometric derivatives on ket center
    else if (num_cent_bra==0 .and. num_cent_ket==1) then
      order_cent_bra = 0
      ! gets the indices and orders of differentitated centers
      call NaryTreePathGetIdxCent(nary_tree=nary_tree_ket, idx_cent=idx_cent_ket)
      call NaryTreePathGetOrderCent(nary_tree=nary_tree_ket, order_cent=order_cent_ket)
      ! finds the range of shells for calculations
      min_shell_bra = 1
      max_shell_bra = num_shells_bra
      min_shell_ket = num_shells_ket+1
      max_shell_ket = 0
      do ishell = 1, num_shells_ket
        if (sub_shells_ket(ishell)%idx_cent==idx_cent_ket(1)) then
          if (ishell<min_shell_ket) min_shell_ket = ishell
          if (ishell>max_shell_ket) max_shell_ket = ishell
        end if
      end do
      if (min_shell_ket>max_shell_ket) return
    ! no partial geometric derivatives
    else
      order_cent_bra = 0
      order_cent_ket = 0
      min_shell_bra = 1
      max_shell_bra = num_shells_bra
      min_shell_ket = 1
      max_shell_ket = num_shells_ket
    end if
    ! sets the maximum number of AOs in a sub-shell on bra center
    max_ao_bra = 0
    do ishell = 1, num_shells_bra
      size_ao = sub_shells_bra(ishell)%num_ao*sub_shells_bra(ishell)%num_contr
      if (max_ao_bra<size_ao) max_ao_bra = size_ao
    end do
    ! sets the number of AO sub-shell pairs to calculate
    if (same_braket .and. num_shells_bra==num_shells_ket .and. &
        num_cent_bra==0 .and. num_cent_ket==0) then
      ! gets the symmetry of property integrals
      call OnePropGetSymmetry(one_prop=one_prop, prop_sym=prop_sym)
      select case (prop_sym)
      case (SYMM_INT_MAT,ANTI_INT_MAT)
        num_pairs = num_shells_bra*(num_shells_bra+1)/2
      case default
        num_pairs = num_shells_bra*num_shells_ket
      end select
      ! sets the maximum number of AOs in a sub-shell on ket center
      max_ao_ket = max_ao_bra
    else
      ! we need to calculate all the sub-shell pairs even for symmetric and anti-symmetric integrals
      prop_sym = SQUARE_INT_MAT
      num_pairs = (max_shell_bra-min_shell_bra+1)*(max_shell_ket-min_shell_ket+1)
      ! sets the maximum number of AOs in a sub-shell on ket center
      max_ao_ket = 0
      do ishell = min_shell_ket, max_shell_ket
        size_ao = sub_shells_ket(ishell)%num_ao*sub_shells_ket(ishell)%num_contr
        if (max_ao_ket<size_ao) max_ao_ket = size_ao
      end do
    end if
    ! gets the numbers of unique geometric derivatives of current path
    path_num_bgeo = (order_cent_bra(1)+1)*(order_cent_bra(1)+2)/2
    path_num_kgeo = (order_cent_ket(1)+1)*(order_cent_ket(1)+2)/2
    call NaryTreePathGetNumUnique(nary_tree=nary_tree_total, path_num_unique=path_num_tgeo)
    ! gets the offsets of geometric derivatives of current path
    call NaryTreePathGetOffset(nary_tree=nary_tree_bra, path_offset=path_offset_bgeo)
    call NaryTreePathGetOffset(nary_tree=nary_tree_ket, path_offset=path_offset_kgeo)
    call NaryTreePathGetOffset(nary_tree=nary_tree_total, path_offset=path_offset_tgeo)
    ! sets the number of integral matrices
    num_matrices = num_prop*path_num_bgeo*path_num_kgeo*path_num_tgeo
    ! allocates cache of contracted integrals between two AO sub-shells
    size_ints = max_ao_bra*max_ao_ket*num_points*num_matrices
    allocate(contr_ints(size_ints), stat=ierr)
    if (ierr/=0) then
      call quit("Gen1IntShellGetFunExpt>> failed to allocate contr_ints!")
    end if
#if defined(VAR_MPI)
    ! gets the rank of this processor and the number of processors
    if (present(api_comm)) then
      call MPI_Comm_rank(api_comm, rank_proc, ierr)
      call MPI_Comm_size(api_comm, num_proc, ierr)
    else
      rank_proc = MANAGER
      num_proc = 1
    end if
    ! manager processor has argument \var(val_expt)
    if (rank_proc==MANAGER) then
#endif
      ! if calculating expectation values
      do_expectation = present(val_expt)
      ! sets the array of expectation values
      if (do_expectation) then
        unique_expt => val_expt
      end if
#if defined(VAR_MPI)
      ! broadcasts what kind of jobs to do
      if (present(api_comm)) then
        call MPI_Bcast(do_expectation, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
      end if
    ! worker processors do not have argument \var(val_expt)
    else
      call MPI_Bcast(do_expectation, 1, MPI_LOGICAL, MANAGER, api_comm, ierr)
      if (do_expectation) then
        allocate(unique_expt(num_points,                                        &
                             num_prop,                                          &
                             path_offset_bgeo+1:path_offset_bgeo+path_num_bgeo, &
                             path_offset_kgeo+1:path_offset_kgeo+path_num_kgeo, &
                             path_offset_tgeo+1:path_offset_tgeo+path_num_tgeo, &
                             num_dens), stat=ierr)
        if (ierr/=0) then
          call quit("Gen1IntShellGetFunExpt>> failed to allocate unique_expt on worker processor!")
        end if
        unique_expt = 0.0_REALK
      end if
    end if
    ! calculations with more than one processor
    if (num_proc>1) then
      ! sends the message tag
      msg_tag = 1
      ! manager code
      if (rank_proc==MANAGER) then
        ! initializes the number of remaining jobs (the manager needs to send "finish" signal
        ! to other worker processors)
        remaining_jobs = num_proc-1
        ! intializes the indices of AO sub-shell pair on bra and ket centers
        shell_pair(1) = min_shell_bra-1
        shell_pair(2) = min_shell_ket
        do while (remaining_jobs>0)
          ! receives a request from a woker
          call MPI_Recv(worker_request, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                        MPI_ANY_TAG, api_comm, mpi_status, ierr)
          ! the worker requests new work
          if (worker_request(1)==REQUEST_WORK) then
            ! no more sub-shell pair to calculate
            if (shell_pair(1)>=max_shell_bra .and. shell_pair(2)>=max_shell_ket) then
              call MPI_Send((/NO_MORE_WORK,NO_MORE_WORK/), 2, MPI_INTEGER, &
                            worker_request(3), msg_tag, api_comm, ierr)
              ! decreases the number of remaining jobs
              remaining_jobs = remaining_jobs-1
            else
              ! prepares the next sub-shell pair to calculate
              select case (prop_sym)
              ! for symmetric and anti-symmetric matrices, we only calculate the upper
              ! and diagonal parts explicitly
              case (SYMM_INT_MAT,ANTI_INT_MAT)
                if (shell_pair(1)<shell_pair(2)) then
                  shell_pair(1) = shell_pair(1)+1
                else
                  shell_pair(2) = shell_pair(2)+1
                  shell_pair(1) = 1
                end if
              case default
                if (shell_pair(1)<max_shell_bra) then
                  shell_pair(1) = shell_pair(1)+1
                else
                  shell_pair(2) = shell_pair(2)+1
                  shell_pair(1) = min_shell_bra
                end if
              end select
              ! sends the next sub-shell pair to the worker
              call MPI_Send(shell_pair, 2, MPI_INTEGER, worker_request(3), &
                            msg_tag, api_comm, ierr)
            end if
          ! the worker wants to send the contracted integrals back
          else
            call quit("Gen1IntShellGetFunExpt>> weird error!")
          end if
        end do
        ! receives expectation values from worker processors
        if (do_expectation) then
          call MPI_Reduce(MPI_IN_PLACE,                                                     &
                          unique_expt(:,:,                                                  &
                                      path_offset_bgeo+1:path_offset_bgeo+path_num_bgeo,    &
                                      path_offset_kgeo+1:path_offset_kgeo+path_num_kgeo,    &
                                      path_offset_tgeo+1:path_offset_tgeo+path_num_tgeo,:), &
                          num_points*num_matrices*num_dens, MPI_REALK, MPI_SUM, MANAGER,    &
                          api_comm, ierr)
        end if
      ! worker code
      else
        worker_request(3) = rank_proc
        do while (.true.)
          ! sends request for a new work to manager
          worker_request(1) = REQUEST_WORK
          call MPI_Send(worker_request, 3, MPI_INTEGER, MANAGER, &
                        msg_tag, api_comm, ierr)
          ! receives the next sub-shell pair or "finished" signal from manager
          call MPI_Recv(shell_pair, 2, MPI_INTEGER, MANAGER, MPI_ANY_TAG, &
                        api_comm, mpi_status, ierr)
          if (shell_pair(1)==NO_MORE_WORK) then
            exit
          else
            ! calculates the contracted integrals
            size_ao = sub_shells_bra(shell_pair(1))%num_ao    &
                    * sub_shells_bra(shell_pair(1))%num_contr &
                    * sub_shells_ket(shell_pair(2))%num_ao    &
                    * sub_shells_ket(shell_pair(2))%num_contr
            size_ints =  size_ao*num_points*num_matrices
            ! spherical GTOs
            if (spher_gto) then
              ! calls Gen1Int subroutines to evaluate property integrands
              call OnePropGetFunction(                                        &
                     idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                     coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                     angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                     num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                     exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                     num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                     contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                     idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                     coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                     angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                     num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                     exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                     num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                     contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                     spher_gto=spher_gto,                                     &
                     one_prop=one_prop,                                       &
                     order_geo_bra=order_cent_bra(1),                         &
                     order_geo_ket=order_cent_ket(1),                         &
                     nary_tree_total=nary_tree_total,                         &
                     num_points=num_points,                                   &
                     grid_points=grid_points,                                 &
                     num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                     num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                     num_opt=num_matrices, contr_ints=contr_ints)
              ! reorders p-shell spherical GTOs, since Dalton uses x(+1), y(-1) and z(0),
              ! while Gen1Int uses y(-1), z(0) and x(+1)
              if (sub_shells_bra(shell_pair(1))%ang_num==1)                              &
                call gen1int_reorder_p_sgto(1, sub_shells_bra(shell_pair(1))%num_contr   &
                                               *sub_shells_ket(shell_pair(2))%num_ao     &
                                               *sub_shells_ket(shell_pair(2))%num_contr, &
                                            num_points*num_matrices, contr_ints)
              if (sub_shells_ket(shell_pair(2))%ang_num==1)                           &
                call gen1int_reorder_p_sgto(sub_shells_bra(shell_pair(1))%num_ao      &
                                            *sub_shells_bra(shell_pair(1))%num_contr, &
                                            sub_shells_ket(shell_pair(2))%num_contr,  &
                                            num_points*num_matrices, contr_ints)
            ! Cartesian GTOs
            else
              ! calls Gen1Int subroutines to evaluate property integrands, and reorders
              ! the integrands according to Cartesian powers
              call OnePropGetFunction(                                        &
                     idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                     coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                     angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                     num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                     exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                     num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                     contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                     idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                     coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                     angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                     num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                     exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                     num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                     contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                     spher_gto=spher_gto,                                     &
                     one_prop=one_prop,                                       &
                     order_geo_bra=order_cent_bra(1),                         &
                     order_geo_ket=order_cent_ket(1),                         &
                     nary_tree_total=nary_tree_total,                         &
                     num_points=num_points,                                   &
                     grid_points=grid_points,                                 &
                     num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                     num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                     num_opt=num_matrices, contr_ints=contr_ints,             &
                     powers_bra=sub_shells_bra(shell_pair(1))%powers,         &
                     powers_ket=sub_shells_ket(shell_pair(2))%powers)
            end if
            ! calculates the expectation values with unique total geometric derivatives
            if (do_expectation) then
              ! sets the minimum and maximum of indices of rows of the integral matrices
              min_row_idx = sub_shells_bra(shell_pair(1))%base_idx+1
              max_row_idx = sub_shells_bra(shell_pair(1))%base_idx &
                          + sub_shells_bra(shell_pair(1))%num_ao   &
                          * sub_shells_bra(shell_pair(1))%num_contr
              ! sets the minimum and maximum of indices of columns of the integral matrices
              min_col_idx = sub_shells_ket(shell_pair(2))%base_idx+1
              max_col_idx = sub_shells_ket(shell_pair(2))%base_idx &
                          + sub_shells_ket(shell_pair(2))%num_ao   &
                          * sub_shells_ket(shell_pair(2))%num_contr
              do idens = 1, num_dens
                start_ints = 0
                do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
                  do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
                    do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                      do iprop = 1, num_prop
                        do ipoint = 1, num_points
                          end_ints = start_ints+size_ao
                          start_ints = start_ints+1
                          call MatMultBlockedTrace(ao_dens(idens),                  &
                                                   min_row_idx, max_row_idx,        &
                                                   min_col_idx, max_col_idx,        &
                                                   contr_ints(start_ints:end_ints), &
                                                   unique_expt(ipoint,              &
                                                               iprop,               &
                                                               bgeo,                &
                                                               kgeo,                &
                                                               tgeo,                &
                                                               idens),              &
                                                   .false.)
                          if (prop_sym==SYMM_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                            call MatMultBlockedTrace(ao_dens(idens),                  &
                                                     min_row_idx, max_row_idx,        &
                                                     min_col_idx, max_col_idx,        &
                                                     contr_ints(start_ints:end_ints), &
                                                     unique_expt(ipoint,              &
                                                                 iprop,               &
                                                                 bgeo,                &
                                                                 kgeo,                &
                                                                 tgeo,                &
                                                                 idens),              &
                                                     .true.)
                          else if (prop_sym==ANTI_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                            call MatMultBlockedTrace(ao_dens(idens),                   &
                                                     min_row_idx, max_row_idx,         &
                                                     min_col_idx, max_col_idx,         &
                                                     -contr_ints(start_ints:end_ints), &
                                                     unique_expt(ipoint,               &
                                                                 iprop,                &
                                                                 bgeo,                 &
                                                                 kgeo,                 &
                                                                 tgeo,                 &
                                                                 idens),               &
                                                     .true.)
                          end if
                          start_ints = end_ints
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end if
          end if
        end do
        ! sends expectation values to manager processor
        if (do_expectation) then
          call MPI_Reduce(unique_expt, unique_expt, num_points*num_matrices*num_dens, &
                          MPI_REALK, MPI_SUM, MANAGER, api_comm, ierr)
        end if
      end if
    ! calculations with only one processor
    else
#endif
      ! initializes the number of remaining pairs to calculate
      remaining_jobs = num_pairs
      ! intializes the indices of AO sub-shell pair on bra and ket centers
      shell_pair(1) = min_shell_bra-1
      shell_pair(2) = min_shell_ket
      do while (remaining_jobs>0)
        ! prepares the next sub-shell pair to calculate
        select case (prop_sym)
        ! for symmetric and anti-symmetric matrices, we only calculate the upper
        ! and diagonal parts explicitly
        case (SYMM_INT_MAT,ANTI_INT_MAT)
          if (shell_pair(1)<shell_pair(2)) then
            shell_pair(1) = shell_pair(1)+1
          else
            shell_pair(2) = shell_pair(2)+1
            shell_pair(1) = 1
          end if
        case default
          if (shell_pair(1)<max_shell_bra) then
            shell_pair(1) = shell_pair(1)+1
          else
            shell_pair(2) = shell_pair(2)+1
            shell_pair(1) = min_shell_bra
          end if
        end select
        size_ao = sub_shells_bra(shell_pair(1))%num_ao    &
                * sub_shells_bra(shell_pair(1))%num_contr &
                * sub_shells_ket(shell_pair(2))%num_ao    &
                * sub_shells_ket(shell_pair(2))%num_contr
        size_ints =  size_ao*num_points*num_matrices
        ! spherical GTOs
        if (spher_gto) then
          ! calls Gen1Int subroutines to evaluate property integrands
          call OnePropGetFunction(                                        &
                 idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                 coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                 angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                 num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                 exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                 num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                 contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                 idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                 coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                 angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                 num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                 exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                 num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                 contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                 spher_gto=spher_gto,                                     &
                 one_prop=one_prop,                                       &
                 order_geo_bra=order_cent_bra(1),                         &
                 order_geo_ket=order_cent_ket(1),                         &
                 nary_tree_total=nary_tree_total,                         &
                 num_points=num_points,                                   &
                 grid_points=grid_points,                                 &
                 num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                 num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                 num_opt=num_matrices, contr_ints=contr_ints)
          ! reorders p-shell spherical GTOs, since Dalton uses x(+1), y(-1) and z(0),
          ! while Gen1Int uses y(-1), z(0) and x(+1)
          if (sub_shells_bra(shell_pair(1))%ang_num==1)                              &
            call gen1int_reorder_p_sgto(1, sub_shells_bra(shell_pair(1))%num_contr   &
                                           *sub_shells_ket(shell_pair(2))%num_ao     &
                                           *sub_shells_ket(shell_pair(2))%num_contr, &
                                        num_points*num_matrices, contr_ints)
          if (sub_shells_ket(shell_pair(2))%ang_num==1)                           &
            call gen1int_reorder_p_sgto(sub_shells_bra(shell_pair(1))%num_ao      &
                                        *sub_shells_bra(shell_pair(1))%num_contr, &
                                        sub_shells_ket(shell_pair(2))%num_contr,  &
                                        num_points*num_matrices, contr_ints)
        ! Cartesian GTOs
        else
          ! calls Gen1Int subroutines to evaluate property integrands, and reorders
          ! the integrands according to Cartesian powers
          call OnePropGetFunction(                                        &
                 idx_bra=sub_shells_bra(shell_pair(1))%idx_cent,          &
                 coord_bra=sub_shells_bra(shell_pair(1))%coord_cent,      &
                 angular_bra=sub_shells_bra(shell_pair(1))%ang_num,       &
                 num_prim_bra=sub_shells_bra(shell_pair(1))%num_prim,     &
                 exponent_bra=sub_shells_bra(shell_pair(1))%exponents,    &
                 num_contr_bra=sub_shells_bra(shell_pair(1))%num_contr,   &
                 contr_coef_bra=sub_shells_bra(shell_pair(1))%contr_coef, &
                 idx_ket=sub_shells_ket(shell_pair(2))%idx_cent,          &
                 coord_ket=sub_shells_ket(shell_pair(2))%coord_cent,      &
                 angular_ket=sub_shells_ket(shell_pair(2))%ang_num,       &
                 num_prim_ket=sub_shells_ket(shell_pair(2))%num_prim,     &
                 exponent_ket=sub_shells_ket(shell_pair(2))%exponents,    &
                 num_contr_ket=sub_shells_ket(shell_pair(2))%num_contr,   &
                 contr_coef_ket=sub_shells_ket(shell_pair(2))%contr_coef, &
                 spher_gto=spher_gto,                                     &
                 one_prop=one_prop,                                       &
                 order_geo_bra=order_cent_bra(1),                         &
                 order_geo_ket=order_cent_ket(1),                         &
                 nary_tree_total=nary_tree_total,                         &
                 num_points=num_points,                                   &
                 grid_points=grid_points,                                 &
                 num_gto_bra=sub_shells_bra(shell_pair(1))%num_ao,        &
                 num_gto_ket=sub_shells_ket(shell_pair(2))%num_ao,        &
                 num_opt=num_matrices, contr_ints=contr_ints,             &
                 powers_bra=sub_shells_bra(shell_pair(1))%powers,         &
                 powers_ket=sub_shells_ket(shell_pair(2))%powers)
        end if
        ! sets the minimum and maximum of indices of rows of the integral matrices
        min_row_idx = sub_shells_bra(shell_pair(1))%base_idx+1
        max_row_idx = sub_shells_bra(shell_pair(1))%base_idx &
                    + sub_shells_bra(shell_pair(1))%num_ao   &
                    * sub_shells_bra(shell_pair(1))%num_contr
        ! sets the minimum and maximum of indices of columns of the integral matrices
        min_col_idx = sub_shells_ket(shell_pair(2))%base_idx+1
        max_col_idx = sub_shells_ket(shell_pair(2))%base_idx &
                    + sub_shells_ket(shell_pair(2))%num_ao   &
                    * sub_shells_ket(shell_pair(2))%num_contr
        ! calculates the expectation values with unique total geometric derivatives
        if (do_expectation) then
          do idens = 1, num_dens
            start_ints = 0
            do tgeo = path_offset_tgeo+1, path_offset_tgeo+path_num_tgeo
              do kgeo = path_offset_kgeo+1, path_offset_kgeo+path_num_kgeo
                do bgeo = path_offset_bgeo+1, path_offset_bgeo+path_num_bgeo
                  do iprop = 1, num_prop
                    do ipoint = 1, num_points
                      end_ints = start_ints+size_ao
                      start_ints = start_ints+1
                      call MatMultBlockedTrace(ao_dens(idens),                  &
                                               min_row_idx, max_row_idx,        &
                                               min_col_idx, max_col_idx,        &
                                               contr_ints(start_ints:end_ints), &
                                               unique_expt(ipoint,              &
                                                           iprop,               &
                                                           bgeo,                &
                                                           kgeo,                &
                                                           tgeo,                &
                                                           idens),              &
                                               .false.)
                      if (prop_sym==SYMM_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                        call MatMultBlockedTrace(ao_dens(idens),                  &
                                                 min_row_idx, max_row_idx,        &
                                                 min_col_idx, max_col_idx,        &
                                                 contr_ints(start_ints:end_ints), &
                                                 unique_expt(ipoint,              &
                                                             iprop,               &
                                                             bgeo,                &
                                                             kgeo,                &
                                                             tgeo,                &
                                                             idens),              &
                                                 .true.)
                      else if (prop_sym==ANTI_INT_MAT .and. shell_pair(1)/=shell_pair(2)) then
                        call MatMultBlockedTrace(ao_dens(idens),                   &
                                                 min_row_idx, max_row_idx,         &
                                                 min_col_idx, max_col_idx,         &
                                                 -contr_ints(start_ints:end_ints), &
                                                 unique_expt(ipoint,               &
                                                             iprop,                &
                                                             bgeo,                 &
                                                             kgeo,                 &
                                                             tgeo,                 &
                                                             idens),               &
                                                 .true.)
                      end if
                      start_ints = end_ints
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if 
        ! decreases the number of remaining pairs
        remaining_jobs = remaining_jobs-1
      end do
#if defined(VAR_MPI)
    end if
#endif
    ! frees spaces
    deallocate(contr_ints)
    if (do_expectation) then
      nullify(unique_expt)
    end if
#if defined(VAR_MPI)
    ! blocks until all processors have finished
    if (present(api_comm)) call MPI_Barrier(api_comm, ierr)
#endif
  end subroutine Gen1IntShellGetFunExpt

  !> \brief calculates molecular orbitals at grid points
  !> \author Bin Gao
  !> \date 2012-03-11
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
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
!FIXME: check the implementation of \var(order_geo)
  !> \param order_geo is the order of geometric derivatives
  !> \return val_mo contains the value of molecular orbitals at grid points
  !> \note val_mo should be zero before
  subroutine Gen1IntShellGetMO(num_shells, sub_shells, mo_coef, &
                               num_points, grid_points,         &
                               num_derv, num_mo, val_mo,        &
                               api_comm, gto_type,              &
                               order_mag, order_ram, order_geo)
    use gen1int_matrix
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(in) :: sub_shells(num_shells)
    type(matrix), intent(in) :: mo_coef
    integer, intent(in) :: num_points
    real(REALK), intent(in) :: grid_points(3,num_points)
    integer, intent(in) :: num_derv
    integer, intent(in) :: num_mo
    real(REALK), intent(inout) :: val_mo(num_points*num_derv,num_mo)
    integer, optional, intent(in) :: api_comm
    integer, optional, intent(in) :: gto_type
    integer, optional, intent(in) :: order_mag
    integer, optional, intent(in) :: order_ram
    integer, optional, intent(in) :: order_geo
    integer p_gto_type   !type of GTOs
    integer p_order_mag  !order of magnetic derivatives
    integer p_order_ram  !order of derivatives w.r.t. total rotational angular momentum
    integer p_order_geo  !order of geometric derivatives
    integer num_opt      !number of operators including grid points and derivatives
    integer ishell       !incremental recorders over AO sub-shells
    logical spher_gto    !if spherical GTOs
    real(REALK), allocatable :: contr_value(:,:,:)  !contracted GTOs at grid points
    real(REALK), allocatable :: mo_value(:,:,:)     !a block of values of the molecular orbital coefficients
    integer min_row_idx, max_row_idx                !minimum and maximum indices of AOs
    integer imo, iopt, icontr, iao                  !incremental recorders
    integer ierr                                    !error information
#if defined(VAR_MPI)
#include "mpif.h"
    integer rank_proc                               !rank of processor
    integer num_proc                                !number of processors
    integer worker_request(2)                       !request from a worker, in which the first element
                                                    !is either \var(REQUEST_WORK) or the AO sub-shell
                                                    !to send back, the second is the rank of the worker
    integer msg_tag                                 !message tag
    integer mpi_status(MPI_STATUS_SIZE)             !MPI status
#endif
    ! gets the type of GTOs
    if (present(gto_type)) then
      if (gto_type/=NON_LAO) then
        call quit("Gen1IntShellGetMO>> only non-LAO implemented!")
      end if
      p_gto_type = gto_type
    else
      p_gto_type = NON_LAO
    end if
    ! gets the number of derivatives
    num_opt = 1
    if (present(order_mag)) then
      p_order_mag = order_mag
      num_opt = num_opt*(order_mag+1)*(order_mag+2)/2
    else
      p_order_mag = 0
    end if
    if (present(order_ram)) then
      p_order_ram = order_ram
      num_opt = num_opt*(order_ram+1)*(order_ram+2)/2
    else
      p_order_ram = 0
    end if
    if (p_order_mag/=0 .or. p_order_ram/=0) then
      call quit("Gen1IntShellGetMO>> only geometric derivatives implemented!")
    end if
    if (present(order_geo)) then
      p_order_geo = order_geo
      num_opt = num_opt*(order_geo+1)*(order_geo+2)/2
    else
      p_order_geo = 0
    end if
    if (num_opt/=num_derv) then
      call quit("Gen1IntShellGetMO>> incorrect number of derivatives!")
    end if
    num_opt = num_opt*num_points
    ! loops over AO sub-shells
    do ishell = 1, num_shells
      spher_gto = sub_shells(ishell)%spher_gto
      allocate(contr_value(sub_shells(ishell)%num_ao,    &
                           sub_shells(ishell)%num_contr, &
                           num_opt), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntShellGetMO>> failed to allocate contr_value!")
      end if
      ! spherical GTOs
      if (spher_gto) then
        call contr_sgto_value(sub_shells(ishell)%coord_cent,        &
                              sub_shells(ishell)%ang_num,           &
                              sub_shells(ishell)%num_prim,          &
                              sub_shells(ishell)%exponents,         &
                              sub_shells(ishell)%num_contr,         &
                              sub_shells(ishell)%contr_coef,        &
                              p_order_geo, num_points, grid_points, &
                              sub_shells(ishell)%num_ao, num_derv,  &
                              contr_value)
        ! reorders p-shell spherical GTOs
        if (sub_shells(ishell)%ang_num==1)                             &
          call gen1int_reorder_p_sgto(1, sub_shells(ishell)%num_contr, &
                                      num_opt, contr_value)
      ! Cartesian GTOs
      else
        if (sub_shells(ishell)%ang_num==0) then
          call contr_cgto_value(sub_shells(ishell)%coord_cent,        &
                                sub_shells(ishell)%ang_num,           &
                                sub_shells(ishell)%num_prim,          &
                                sub_shells(ishell)%exponents,         &
                                sub_shells(ishell)%num_contr,         &
                                sub_shells(ishell)%contr_coef,        &
                                p_order_geo, num_points, grid_points, &
                                sub_shells(ishell)%num_ao, num_derv,  &
                                contr_value)
        else
          allocate(mo_value(sub_shells(ishell)%num_ao,    &
                            sub_shells(ishell)%num_contr, &
                            num_opt), stat=ierr)
          if (ierr/=0) then
            call quit("Gen1IntShellGetMO>> failed to allocate mo_value!")
          end if
          call contr_cgto_value(sub_shells(ishell)%coord_cent,        &
                                sub_shells(ishell)%ang_num,           &
                                sub_shells(ishell)%num_prim,          &
                                sub_shells(ishell)%exponents,         &
                                sub_shells(ishell)%num_contr,         &
                                sub_shells(ishell)%contr_coef,        &
                                p_order_geo, num_points, grid_points, &
                                sub_shells(ishell)%num_ao, num_derv,  &
                                mo_value)
          ! reorders the results according to Cartesian powers
          call reorder_cgtos(sub_shells(ishell)%ang_num,   &
                             sub_shells(ishell)%num_ao,    &
                             sub_shells(ishell)%powers, 1, &
                             sub_shells(ishell)%num_contr, &
                             num_opt, mo_value, contr_value)
          deallocate(mo_value)
        end if
      end if
      ! sets the minimum and maximum of indices of AOs
      min_row_idx = sub_shells(ishell)%base_idx+1
      max_row_idx = sub_shells(ishell)%base_idx &
                  + sub_shells(ishell)%num_ao   &
                  * sub_shells(ishell)%num_contr
      ! gets a block of values of the molecular orbital coefficients
      allocate(mo_value(sub_shells(ishell)%num_ao,    &
                        sub_shells(ishell)%num_contr, &
                        num_mo), stat=ierr)
      if (ierr/=0) then
        call quit("Gen1IntShellGetMO>> failed to allocate mo_value!")
      end if
      call MatGetValues(mo_coef, min_row_idx, max_row_idx, 1, num_mo, mo_value)
! mo_value(:,:,:) = reshape(mo_coef%elms(min_row_idx : max_row_idx, &
!                                                  1 : num_mo, 1), shape(mo_value))
      ! gets the value of MOs at grid points
      do imo = 1, num_mo
        do iopt = 1, num_opt
          do icontr = 1, sub_shells(ishell)%num_contr
            do iao = 1, sub_shells(ishell)%num_ao
              val_mo(iopt,imo) = val_mo(iopt,imo) &
                               + mo_value(iao,icontr,imo)*contr_value(iao,icontr,iopt)
            end do
          end do
        end do
      end do
      deallocate(mo_value)
      deallocate(contr_value)
    end do
  end subroutine Gen1IntShellGetMO

  !> \brief frees space taken by the AO sub-shells
  !> \author Bin Gao
  !> \date 2011-10-03
  !> \param num_shells is the number of AO sub-shells
  !> \param sub_shells are the AO sub-shells
  subroutine Gen1IntShellDestroy(num_shells, sub_shells)
    integer, intent(in) :: num_shells
    type(sub_shell_t), intent(inout) :: sub_shells(num_shells)
    integer ishell  !incremental recorder
    do ishell = 1, num_shells
      if (allocated(sub_shells(ishell)%exponents)) &
        deallocate(sub_shells(ishell)%exponents)
      if (allocated(sub_shells(ishell)%contr_coef)) &
        deallocate(sub_shells(ishell)%contr_coef)
      !-if (allocated(sub_shells(ishell)%mag_num)) &
      !-  deallocate(sub_shells(ishell)%mag_num)
      if (allocated(sub_shells(ishell)%powers)) &
        deallocate(sub_shells(ishell)%powers)
    end do
  end subroutine Gen1IntShellDestroy

  !> \brief reorders the p-shell contracted real solid-harmonic Gaussians in Dalton's order
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param dim_sgto_bra is the dimension of SGTOs on bra center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param num_opt is the number of operators
  !> \param gen_ints contains the contracted integrals from Gen1Int library
  subroutine gen1int_reorder_p_sgto(dim_sgto_bra, num_contr_ket, num_opt, gen_ints)
    integer, intent(in) :: dim_sgto_bra
    integer, intent(in) :: num_contr_ket
    integer, intent(in) :: num_opt
    real(REALK), intent(inout) :: gen_ints(dim_sgto_bra,3,num_contr_ket,num_opt)
    real(REALK), allocatable :: tmp_ints(:)  !temporary integrals
    integer icontr, iopt                     !incremental recorders
    integer ierr                             !error information
    ! Dalton's order of SGTOs: px(1), py(-1), pz(0),
    ! while those in Gen1Int is: py(-1), pz(0), px(1)
    allocate(tmp_ints(dim_sgto_bra), stat=ierr)
    if (ierr/=0) then
      call quit("gen1int_reorder_p_sgto>> failed to allocate tmp_ints!")
    end if
    do iopt = 1, num_opt
      do icontr = 1, num_contr_ket
        tmp_ints = gen_ints(:,3,icontr,iopt)
        gen_ints(:,3,icontr,iopt) = gen_ints(:,2,icontr,iopt)
        gen_ints(:,2,icontr,iopt) = gen_ints(:,1,icontr,iopt)
        gen_ints(:,1,icontr,iopt) = tmp_ints
      end do
    end do
    deallocate(tmp_ints)
  end subroutine gen1int_reorder_p_sgto

end module gen1int_shell
