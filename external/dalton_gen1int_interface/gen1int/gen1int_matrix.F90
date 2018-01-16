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
!...  This file defines the matrix type and its corresponding functions used in Gen1Int interface.
!
!...  2011-12-10, Bin Gao
!...  * first version

#include "gen1int_host.h"

!> \brief matrix module used in Gen1Int interface
!> \note the reason of introducing matrix type is due to other codes (like OpenRSP) ask
!>       for a matrix type output, in order to maintain only one Gen1Int interface, I
!>       have therefore introduced such matrix module (or wrapper); users who would like
!>       to have their results in an array could use \fn(MatAssociate) and \fn(MatNullify)
!>       subroutines
!> \author Bin Gao
!> \date 2011-12-10
module gen1int_matrix

#ifdef BUILD_OPENRSP
  ! matrix module in OpenRSP
  use matrix_lowlevel, &
           MatIsClosedShell    => mat_is_closed_shell, &
           MatDestroy          => mat_remove, &
           MatGetValues        => mat_get_block, &
           MatSetValues        => mat_set_block, &
           MatBcast            => mat_mpi_bcast, &
           matrix              => openrsp_matrix

  implicit none

  public :: MatView
  public :: MatCreate
  public :: MatMultBlockedTrace
  public :: MatArrayAlmostEqual

  contains

  !> \brief visualizes the information of a matrix
  !> \author Bin Gao
  !> \date 2012-05-14
  !> \param A is the matrix
  !> \param io_viewer is the logical unit number of the viewer
  subroutine MatView(A, io_viewer)
    type(matrix), intent(in) :: A
    integer, intent(in) :: io_viewer
    call mat_print(A, unit=io_viewer)
  end subroutine MatView

  !> \brief creates a matrix
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \param num_row is the number of rows
  !> \param num_col is the number of columns
  !> \param triangular indicates if the matrix element storage is in triangular format
  !> \param symmetric indicates if the matrix is symmetric or anti-symmetric in triangular format
  !> \return A is the matrix
  !> \return info_mat (/=0) indicates error happened when creating the matrix
  subroutine MatCreate(A, num_row, info_mat, num_col, triangular, symmetric)
    type(matrix), intent(inout) :: A
    integer, intent(in) :: num_row
    integer, intent(out) :: info_mat
    integer, optional, intent(in) :: num_col
    logical, optional, intent(in) :: triangular
    logical, optional, intent(in) :: symmetric
    integer nrow, ncol
    nrow  = num_row
    if (present(num_col)) then
      ncol = num_col
    else
      ncol = num_row
    end if
    if (present(triangular)) then
      if (triangular .and. (nrow/=ncol)) then
        info_mat = ncol !ajt FIXME ncol could be zero, indicating OK
        return
      end if
    end if
    if (present(symmetric)) then
      if (symmetric .and. (nrow/=ncol)) then
        info_mat = A%ncol !ajt FIXME ncol could be zero, indicating OK
        return
      end if
    end if
    call mat_init(A, nrow, ncol)
    info_mat = 0
  end subroutine MatCreate


  !> \brief gets the trace of the product of a block of values and the
  !>        corresponding part of a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param A is the matrix
  !> \param min_row_idx is the minimum row index
  !> \param max_row_idx is the maximum row index
  !> \param min_col_idx is the minimum column index
  !> \param max_col_idx is the maximum column index
  !> \param values contains the block of values
  !> \param trans indicates if transposing the values
  !> \return trace is the updated trace
  subroutine MatMultBlockedTrace(A, min_row_idx, max_row_idx, &
                                 min_col_idx, max_col_idx,    &
                                 values, trace, trans)
    type(matrix), intent(in) :: A
    integer, intent(in) :: min_row_idx
    integer, intent(in) :: max_row_idx
    integer, intent(in) :: min_col_idx
    integer, intent(in) :: max_col_idx
    real(8), intent(in) :: values(min_row_idx:max_row_idx,min_col_idx:max_col_idx)
    real(8), intent(inout) :: trace
    logical, optional, intent(in) :: trans
    real(8) blocked_trace  !temporary result for getting the trace
    logical p_trans            !indicates if transposing the values (private)
    integer irow, icol         !incremental recorders over rows and columns
    if (present(trans)) then
      p_trans = trans
    else
      p_trans = .false.
    end if
    ! mat_dot_block does the same, except 'trans' carries the opposite meaning
    trace = trace + mat_dot_block(A, min_row_idx, max_row_idx, &
                                  min_col_idx, max_col_idx, &
                                  values, .not.p_trans)
  end subroutine MatMultBlockedTrace


  !> \brief checks if the elements of a matrix are almost equal to those of a given array
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \param A is the matrix
  !> \param values is the array
  !> \param io_viewer is the logical unit number of the viewer
  !> \param triangular indicates if the array is in triangular format
  !> \param symmetric indicates if the array is symmetric or anti-symmetric in triangular format
  !> \param threshold is the threshold for comparison
  !> \param ratio_thrsh is the threshold of ratio between the matrix element and array
  !> \return almost_equal indicates if the elements of the matrix are almost equal to
  !>         those of the array
  subroutine MatArrayAlmostEqual(A, values, io_viewer, almost_equal, &
                                 triangular, symmetric, threshold,   &
                                 ratio_thrsh)
    type(matrix), intent(in) :: A
    real(8), intent(in) :: values(:)
    integer, intent(in) :: io_viewer
    logical, intent(out) :: almost_equal
    logical, optional, intent(in) :: triangular
    logical, optional, intent(in) :: symmetric
    real(8), optional, intent(in) :: threshold
    real(8), optional, intent(in) :: ratio_thrsh
    real(8) p_threshold       !threshold for comparison (private)
    real(8) p_ratio_thrsh(2)  !threshold of ratio between the matrix element and array(private)
    logical p_triangular      !if the array is in triangular format (private)
    logical p_symmetric       !if the array is symmetric or anti-symmetric in triangular format (private)
    integer nrow, ncol        !dimensions
    integer irow, icol        !incremental recorders over rows and columns
    integer addr_val          !address of the array value
    real(8) mat_to_val        !ratio between the matrix element and array
    real(8) tmp_val           !temporary value for comparision
    real(8), pointer :: elms(:,:)
    almost_equal = .true.
    if (present(threshold)) then
      p_threshold = threshold
    else
      p_threshold = 10.0d0**(-8)
    end if
    if (present(ratio_thrsh)) then
      p_ratio_thrsh(1) = 1.0d0-ratio_thrsh
      p_ratio_thrsh(2) = 1.0d0+ratio_thrsh
    else
      p_ratio_thrsh(1) = 1.0d0-10.0d0**(-6)
      p_ratio_thrsh(2) = 1.0d0+10.0d0**(-6)
    end if
    if (present(triangular)) then
      p_triangular = triangular
    else
      p_triangular = .false.
    end if
    if (p_triangular) then
      if (present(symmetric)) then
        p_symmetric = symmetric
      else
        p_symmetric = .true.
      end if
    else
      p_symmetric = .true.
    end if
    ! obtain dimensions and pointer to elements
    nrow = mat_get_nrow(A)
    ncol = mat_get_ncol(A)
    elms => mat_acquire_block(A, 1, nrow, 1, ncol)
    ! array is in triangular format
    if (p_triangular) then
      if (nrow/=ncol) then
        write(io_viewer,100) "matrix is a non-sqaure matrix"
        almost_equal = .false.
        return
      end if
      if (nrow*(nrow+1)/2/=size(values)) then
        write(io_viewer,100) "matrix and array have different sizes", &
                             nrow*(nrow+1)/2, size(values)
        almost_equal = .false.
        return
      end if
      addr_val = 0
      do icol = 1, ncol
        do irow = 1, nrow
          addr_val = addr_val+1
          if (irow<=icol) then
            tmp_val = values(irow+icol*(icol-1)/2)
          else
            if (p_symmetric) then
              tmp_val = values(icol+irow*(irow-1)/2)
            else
              tmp_val = -values(icol+irow*(irow-1)/2)
            end if
          end if
          ! we check the ratio for absolutely greater values
          if (abs(tmp_val)>p_threshold) then
            mat_to_val = elms(irow,icol)/tmp_val
            if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
              write(io_viewer,110) irow, icol, elms(irow,icol), tmp_val
              almost_equal = .false.
            end if
          ! checks the difference for smaller values
          else
            if (abs(elms(irow,icol)-tmp_val)>p_threshold) then
              write(io_viewer,110) irow, icol, elms(irow,icol), tmp_val
              almost_equal = .false.
            end if
          end if
        end do
      end do
    ! array is in square format
    else
      if (size(elms)/=size(values)) then
        write(io_viewer,100) "matrix and array have different sizes", &
                             size(elms), size(values)
        almost_equal = .false.
        return
      end if
      addr_val = 0
      do icol = 1, ncol
        do irow = 1, nrow
          addr_val = addr_val+1
          ! we check the ratio for absolutely greater values
          if (abs(values(addr_val))>p_threshold) then
            mat_to_val = elms(irow,icol)/values(addr_val)
            if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
              write(io_viewer,110) irow, icol, elms(irow,icol), values(addr_val)
              almost_equal = .false.
            end if
          ! checks the difference for smaller values
          else
            if (abs(elms(irow,icol)-values(addr_val))>p_threshold) then
              write(io_viewer,110) irow, icol, elms(irow,icol), values(addr_val)
              almost_equal = .false.
            end if
          end if
        end do
      end do
    end if
    call mat_unacquire_block(A, elms)
100 format("MatArrayAlmostEqual>> ",A,I8,",",I8)
110 format("MatArrayAlmostEqual>> element(",I4,",",I4,")",Es16.8," (Mat),",Es16.8," (array)")
  end subroutine MatArrayAlmostEqual

#else
  implicit none

  ! definition of matrix
  type, public :: matrix
    private
    ! if matrix element storage in triangular format
    logical :: triangular = .false.
    ! if symmetric or anti-symmetric matrix in triangular format
    logical :: symmetric = .true.
    ! number of rows
    integer :: num_row = huge(1)
    ! number of columns
    integer :: num_col = huge(1)
    ! pointer to matrix element storage alpha part
    real(REALK), pointer :: elms_alpha(:)
    !-! pointer to matrix element storage for beta part
    !-real(REALK), pointer :: elms_beta(:)
  end type matrix

  ! public subroutines which may be used before calling Gen1Int interface
  public :: MatAssociate
  public :: MatNullify
  public :: MatView

  ! public subroutines used in Gen1Int interface
  public :: MatSetValues
  public :: MatGetValues
  public :: MatMultBlockedTrace

  ! public subroutines for test suite of Gen1Int interface
  public :: MatCreate
#if defined(VAR_MPI)
  public :: MatBcast
#endif
  public :: MatDestroy
  public :: MatArrayAlmostEqual

  ! test suite of matrix module
  public :: MatTestSuite

  contains

  !> \brief associates the matrix element storage to Dalton/dirac workspace
  !> \author Bin Gao
  !> \date 2011-12-10
  !> \param work_alpha is the Dalton/Dirac workspace
  !> \param num_row is the number of rows
  !> \param triangular indicates if the matrix element storage is in triangular format
  !> \param symmetric indicates if the matrix is symmetric or anti-symmetric in triangular format
  !> \return A is the matrix
  !> \return info_mat (/=0) indicates error happened when associating the matrix elements
  subroutine MatAssociate(work_alpha, num_row, A, info_mat, triangular, symmetric)
    real(REALK), target, intent(in) :: work_alpha(:)
    integer, intent(in) :: num_row
    type(matrix), intent(inout) :: A
    integer, intent(out) :: info_mat
    logical, optional, intent(in) :: triangular
    logical, optional, intent(in) :: symmetric
    integer num_elms  !number of elements
    if (present(triangular)) A%triangular = triangular
    num_elms = size(work_alpha)
    ! matrix element storage in triangular format
    if (A%triangular) then
      info_mat = num_elms-num_row*(num_row+1)/2
      if (info_mat==0) then
        if (present(symmetric)) A%symmetric = symmetric
        A%num_row = num_row
        A%num_col = num_row
        A%elms_alpha => work_alpha
      ! incorrect number of elements
      else
        A%triangular = .false.
      end if
    ! matrix element storage in square format
    else
      info_mat = mod(num_elms,num_row)
      if (info_mat==0) then
        A%num_row = num_row
        A%num_col = num_elms/num_row
        A%elms_alpha => work_alpha
      end if
    end if
  end subroutine MatAssociate

  !> \brief deassociates the matrix element storage to Dalton/dirac workspace
  !> \author Bin Gao
  !> \date 2011-12-10
  !> \param A is the matrix
  subroutine MatNullify(A)
    type(matrix), intent(inout) :: A
    nullify(A%elms_alpha)
    A%triangular = .false.
    A%symmetric = .true.
    A%num_row = huge(1)
    A%num_col = huge(1)
  end subroutine MatNullify

  !> \brief visualizes the information of a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param A is the matrix
  !> \param io_viewer is the logical unit number of the viewer
  subroutine MatView(A, io_viewer)
    type(matrix), intent(in) :: A
    integer, intent(in) :: io_viewer
! number of columns per row
#if !defined(NCOL_PER_ROW)
#define NCOL_PER_ROW 5
#endif
! width of numbers
#if !defined(WIDTH_NUMBER)
#define WIDTH_NUMBER 16
#define ES_FORMAT Es16.6
#endif
    !integer, parameter :: NCOL_PER_ROW = 5  !number of columns per row
    integer strt_col, end_col  !incremental recorders for dumping columns
    integer irow, icol, jcol   !incremental recorders over rows and columns
    character(NCOL_PER_ROW*WIDTH_NUMBER), parameter :: &
      LEADING_BLANKS = " "     !leading blanks for each row
    write(io_viewer,100) "number of rows", A%num_row
    write(io_viewer,100) "number of columns", A%num_col
    ! matrix element storage in triangular format
    if (A%triangular) then
      write(io_viewer,100) "matrix element storage in triangular format"
      if (A%symmetric) then
        write(io_viewer,100) "matrix is symmetric"
      else
        write(io_viewer,100) "matrix is anti-symmetric"
      end if
      if (associated(A%elms_alpha)) then
        ! initializes the start and end columns to dump
        strt_col = 1
        end_col = min(A%num_col,NCOL_PER_ROW)
        ! loops over columns
        do while (strt_col<=A%num_col)
          ! dumps the column information
          write(io_viewer,110) ("Column", icol, icol=strt_col,end_col)
          ! loops over rows
          do irow = 1, min(A%num_row,end_col)
            ! the element \f$(i,j)\f$ of an upper triangular matrix in
            ! the column-major order is \f$\frac{j(j-1)}{2}+i\f$
            jcol = max(irow,strt_col)
            write(io_viewer,120) irow,                          &
              LEADING_BLANKS(1:WIDTH_NUMBER*(jcol-strt_col)+1), &
              (A%elms_alpha(icol*(icol-1)/2+irow), icol=jcol,end_col)
          end do
          write(io_viewer,"()")
          ! increases the column recorders
          strt_col = strt_col+NCOL_PER_ROW
          end_col = min(end_col+NCOL_PER_ROW,A%num_col)
        end do
      end if
    ! matrix element storage in square format
    else
      write(io_viewer,100) "matrix element storage in square format"
      if (associated(A%elms_alpha)) then
        ! the index of ending column
        end_col = min(NCOL_PER_ROW,A%num_col)
        ! loops over columns
        do strt_col = 1, A%num_col, NCOL_PER_ROW
          ! dumps the column information
          write(io_viewer,110) ("Column", icol, icol=strt_col,end_col)
          ! loops over rows
          do irow = 1, A%num_row
            write(io_viewer,120) irow, " ", &
             (A%elms_alpha((icol-1)*A%num_row+irow), icol=strt_col,end_col)
          end do
          ! increases the index of ending column
          end_col = min(end_col+NCOL_PER_ROW,A%num_col)
        end do
        write(io_viewer,"()")
      end if
    end if
100 format("MatView>> ",A,I6)
110 format(6X,NCOL_PER_ROW(4X,A,I6))
120 format(I5,A,NCOL_PER_ROW(ES_FORMAT))
  end subroutine MatView

  !> \brief inserts a block of values into a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param min_row_idx is the minimum row index
  !> \param max_row_idx is the maximum row index
  !> \param min_col_idx is the minimum column index
  !> \param max_col_idx is the maximum column index
  !> \param values contains the block of values
  !> \param trans indicates if transposing the values
  !> \return A is the matrix
  subroutine MatSetValues(A, min_row_idx, max_row_idx, &
                          min_col_idx, max_col_idx,    &
                          values, trans)
    type(matrix), intent(inout) :: A
    integer, intent(in) :: min_row_idx
    integer, intent(in) :: max_row_idx
    integer, intent(in) :: min_col_idx
    integer, intent(in) :: max_col_idx
    real(REALK), intent(in) :: values(min_row_idx:max_row_idx,min_col_idx:max_col_idx)
    logical, optional, intent(in) :: trans
    logical p_trans     !indicates if transposing the values (private)
    integer irow, icol  !incremental recorders over rows and columns
    integer base_row    !base address of a specific row
    integer base_col    !base address of a specific column
    if (present(trans)) then
      p_trans = trans
    else
      p_trans = .false.
    end if
    ! matrix element storage in triangular format
    if (A%triangular) then
      ! diagonal parts
      if (min_row_idx==min_col_idx .and. max_row_idx==max_col_idx) then
        do icol = min_col_idx, max_col_idx
          base_col = icol*(icol-1)/2
          do irow = min_col_idx, icol
            A%elms_alpha(base_col+irow) = values(irow,icol)
          end do
        end do
      ! off-diagonal parts (only upper part stored)
      else
        if (.not.p_trans) then
          do icol = min_col_idx, max_col_idx
            base_col = icol*(icol-1)/2
            do irow = min_row_idx, max_row_idx
              A%elms_alpha(base_col+irow) = values(irow,icol)
            end do
          end do
        end if
      end if
    ! matrix element storage in square format
    else
      if (p_trans) then
        do irow = min_row_idx, max_row_idx
          base_row = (irow-1)*A%num_row
          do icol = min_col_idx, max_col_idx
            A%elms_alpha(icol+base_row) = values(irow,icol)
          end do
        end do
      else
        do icol = min_col_idx, max_col_idx
          base_col = (icol-1)*A%num_row
          do irow = min_row_idx, max_row_idx
            A%elms_alpha(irow+base_col) = values(irow,icol)
          end do
        end do
      end if
    end if
  end subroutine MatSetValues

  !> \brief gets a block of values from a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param A is the matrix
  !> \param min_row_idx is the minimum row index
  !> \param max_row_idx is the maximum row index
  !> \param min_col_idx is the minimum column index
  !> \param max_col_idx is the maximum column index
  !> \return values contains the block of values
  subroutine MatGetValues(A, min_row_idx, max_row_idx, &
                          min_col_idx, max_col_idx,    &
                          values)
    type(matrix), intent(in) :: A
    integer, intent(in) :: min_row_idx
    integer, intent(in) :: max_row_idx
    integer, intent(in) :: min_col_idx
    integer, intent(in) :: max_col_idx
    real(REALK), intent(out) :: values(min_row_idx:max_row_idx,min_col_idx:max_col_idx)
    integer irow, icol  !incremental recorders over rows and columns
    integer base_col    !base address of a specific column
    integer addr_elms   !address of element in matrix
    ! matrix element storage in triangular format
    if (A%triangular) then
      ! symmetric matrix
      if (A%symmetric) then
        do icol = min_col_idx, max_col_idx
          do irow = min_row_idx, max_row_idx
            if (icol>=irow) then
              addr_elms = icol*(icol-1)/2+irow
            else
              addr_elms = irow*(irow-1)/2+icol
            end if
            values(irow,icol) = A%elms_alpha(addr_elms)
          end do
        end do
      ! anti-symmetric matrix
      else
        do icol = min_col_idx, max_col_idx
          do irow = min_row_idx, max_row_idx
            if (icol>=irow) then
              addr_elms = icol*(icol-1)/2+irow
              values(irow,icol) = A%elms_alpha(addr_elms)
            else
              addr_elms = irow*(irow-1)/2+icol
              values(irow,icol) = -A%elms_alpha(addr_elms)
            end if
          end do
        end do
      end if
    ! matrix element storage in square format
    else
      do icol = min_col_idx, max_col_idx
        base_col = (icol-1)*A%num_row
        do irow = min_row_idx, max_row_idx
          values(irow,icol) = A%elms_alpha(irow+base_col)
        end do
      end do
    end if
  end subroutine MatGetValues

  !> \brief gets the trace of the product of a block of values and the
  !>        corresponding part of a matrix
  !> \author Bin Gao
  !> \date 2012-01-17
  !> \param A is the matrix
  !> \param min_row_idx is the minimum row index
  !> \param max_row_idx is the maximum row index
  !> \param min_col_idx is the minimum column index
  !> \param max_col_idx is the maximum column index
  !> \param values contains the block of values
  !> \param trans indicates if transposing the values
  !> \return trace is the updated trace
  subroutine MatMultBlockedTrace(A, min_row_idx, max_row_idx, &
                                 min_col_idx, max_col_idx,    &
                                 values, trace, trans)
    type(matrix), intent(in) :: A
    integer, intent(in) :: min_row_idx
    integer, intent(in) :: max_row_idx
    integer, intent(in) :: min_col_idx
    integer, intent(in) :: max_col_idx
    real(REALK), intent(in) :: values(min_row_idx:max_row_idx,min_col_idx:max_col_idx)
    real(REALK), intent(inout) :: trace
    logical, optional, intent(in) :: trans
    logical p_trans            !indicates if transposing the values (private)
    real(REALK) blocked_trace  !temporary result for getting the trace
    integer irow, icol         !incremental recorders over rows and columns
    integer base_row           !base address of a specific row
    integer base_col           !base address of a specific column
    if (present(trans)) then
      p_trans = trans
    else
      p_trans = .false.
    end if
    blocked_trace = 0.0_REALK
    ! matrix element storage in triangular format (only upper and diagonal parts stored)
    if (A%triangular) then
      ! transposes the \var(values)
      if (p_trans) then
        do icol = min_col_idx, max_col_idx
          do irow = min_row_idx, max_row_idx
            if (irow<=icol) then
              blocked_trace = blocked_trace &
                            + A%elms_alpha(irow+icol*(icol-1)/2)*values(irow,icol)
            ! lower part
            else
              if (A%symmetric) then
                blocked_trace = blocked_trace &
                              + A%elms_alpha(icol+irow*(irow-1)/2)*values(irow,icol)
              else
                blocked_trace = blocked_trace &
                              - A%elms_alpha(icol+irow*(irow-1)/2)*values(irow,icol)
              end if
            end if
          end do
        end do
      else
        do irow = min_row_idx, max_row_idx
          do icol = min_col_idx, max_col_idx
            if (icol<=irow) then
              blocked_trace = blocked_trace &
                            + A%elms_alpha(icol+irow*(irow-1)/2)*values(irow,icol)
            ! lower part
            else
              if (A%symmetric) then
                blocked_trace = blocked_trace &
                              + A%elms_alpha(irow+icol*(icol-1)/2)*values(irow,icol)
              else
                blocked_trace = blocked_trace &
                              - A%elms_alpha(irow+icol*(icol-1)/2)*values(irow,icol)
              end if
            end if
          end do
        end do
      end if
    ! matrix element storage in square format
    else
      ! transposes the \var(values)
      if (p_trans) then
        do icol = min_col_idx, max_col_idx
          base_col = (icol-1)*A%num_col
          do irow = min_row_idx, max_row_idx
            blocked_trace = blocked_trace+A%elms_alpha(irow+base_col)*values(irow,icol)
          end do
        end do
      else
        do irow = min_row_idx, max_row_idx
          base_row = (irow-1)*A%num_row
          do icol = min_col_idx, max_col_idx
            blocked_trace = blocked_trace+A%elms_alpha(icol+base_row)*values(irow,icol)
          end do
        end do
      end if
    end if
    ! closed-shell, then additional factor 2
    trace = trace+2.0_REALK*blocked_trace
  end subroutine MatMultBlockedTrace

  !> \brief creates a matrix
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \param num_row is the number of rows
  !> \param num_col is the number of columns
  !> \param triangular indicates if the matrix element storage is in triangular format
  !> \param symmetric indicates if the matrix is symmetric or anti-symmetric in triangular format
  !> \return A is the matrix
  !> \return info_mat (/=0) indicates error happened when creating the matrix
  subroutine MatCreate(A, num_row, info_mat, num_col, triangular, symmetric)
    type(matrix), intent(inout) :: A
    integer, intent(in) :: num_row
    integer, intent(out) :: info_mat
    integer, optional, intent(in) :: num_col
    logical, optional, intent(in) :: triangular
    logical, optional, intent(in) :: symmetric
    A%num_row = num_row
    if (present(num_col)) then
      A%num_col = num_col
    else
      A%num_col = num_row
    end if
    if (present(triangular)) A%triangular = triangular
    ! matrix element storage in triangular format
    if (A%triangular) then
      if (A%num_row==A%num_col) then
        allocate(A%elms_alpha(A%num_row*(A%num_row+1)/2), stat=info_mat)
        A%elms_alpha = 0.0d0
        if (info_mat/=0) then
          A%num_row = huge(1)
          A%num_col = huge(1)
          A%triangular = .false.
        else
          if (present(symmetric)) A%symmetric = symmetric
        end if
      ! we can not use triangular format for non-square matrix
      else
        A%num_row = huge(1)
        A%num_col = huge(1)
        A%triangular = .false.
        info_mat = A%num_col
      end if
    ! matrix element storage in square format
    else
      allocate(A%elms_alpha(A%num_row*A%num_col), stat=info_mat)
      A%elms_alpha = 0.0d0
      if (info_mat/=0) then
        A%num_row = huge(1)
        A%num_col = huge(1)
        A%triangular = .false.
      end if
    end if
    A%elms_alpha = 0.0
  end subroutine MatCreate

#if defined(VAR_MPI)
  !> \brief broadcasts matrix
  !> \author Bin Gao
  !> \date 2012-05-13
  !> \param A is the matrix
  !> \param root is the root processor which broadcasts the matrix
  !> \param mat_comm is the MPI communicator
  subroutine MatBcast(A, root, mat_comm)
    type(matrix), intent(inout) :: A
    integer, intent(in) :: root
    integer, intent(in) :: mat_comm
#include "mpif.h"
    integer rank_proc  !rank of processor
    integer ierr       !error information
    ! broadcasts information of the matrix
    call MPI_Bcast(A%num_row, 1, MPI_INTEGER, root, mat_comm, ierr)
    call MPI_Bcast(A%num_col, 1, MPI_INTEGER, root, mat_comm, ierr)
    call MPI_Bcast(A%triangular, 1, MPI_LOGICAL, root, mat_comm, ierr)
    call MPI_Bcast(A%symmetric, 1, MPI_LOGICAL, root, mat_comm, ierr)
    ! allocates memory for the matrix on other processors
    call MPI_Comm_rank(mat_comm, rank_proc, ierr)
    if (rank_proc==root) then
      call MPI_Bcast(A%elms_alpha, size(A%elms_alpha), MPI_REALK, root, mat_comm, ierr)
    else
      call MatCreate(A=A, num_row=A%num_row, info_mat=ierr, num_col=A%num_col, &
                     triangular=A%triangular, symmetric=A%symmetric)
      call MPI_Bcast(A%elms_alpha, size(A%elms_alpha), MPI_REALK, root, mat_comm, ierr)
    end if
  end subroutine MatBcast
#endif

  !> \brief frees space taken by a matrix
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \return A is the matrix
  subroutine MatDestroy(A)
    type(matrix), intent(inout) :: A
    deallocate(A%elms_alpha)
    nullify(A%elms_alpha)
    A%triangular = .false.
    A%symmetric = .true.
    A%num_row = huge(1)
    A%num_col = huge(1)
  end subroutine MatDestroy

  !> \brief checks if the elements of a matrix are almost equal to those of a given array
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \param A is the matrix
  !> \param values is the array
  !> \param io_viewer is the logical unit number of the viewer
  !> \param triangular indicates if the array is in triangular format
  !> \param symmetric indicates if the array is symmetric or anti-symmetric in triangular format
  !> \param threshold is the threshold for comparison
  !> \param ratio_thrsh is the threshold of ratio between the matrix element and array
  !> \return almost_equal indicates if the elements of the matrix are almost equal to
  !>         those of the array
  subroutine MatArrayAlmostEqual(A, values, io_viewer, almost_equal, &
                                 triangular, symmetric, threshold,   &
                                 ratio_thrsh)
    type(matrix), intent(in) :: A
    real(REALK), intent(in) :: values(:)
    integer, intent(in) :: io_viewer
    logical, intent(out) :: almost_equal
    logical, optional, intent(in) :: triangular
    logical, optional, intent(in) :: symmetric
    real(REALK), optional, intent(in) :: threshold
    real(REALK), optional, intent(in) :: ratio_thrsh
    real(REALK) p_threshold       !threshold for comparison (private)
    real(REALK) p_ratio_thrsh(2)  !threshold of ratio between the matrix element and array(private)
    logical p_triangular          !if the array is in triangular format (private)
    logical p_symmetric           !if the array is symmetric or anti-symmetric in triangular format (private)
    integer irow, icol            !incremental recorders over rows and columns
    integer addr_elm              !address of element
    real(REALK) mat_to_val        !ratio between the matrix element and array
    real(REALK) tmp_val           !temporary value for comparision
    almost_equal = .true.
    if (present(threshold)) then
      p_threshold = threshold
    else
      p_threshold = 10.0_REALK**(-8)
    end if
    if (present(ratio_thrsh)) then
      p_ratio_thrsh(1) = 1.0_REALK-ratio_thrsh
      p_ratio_thrsh(2) = 1.0_REALK+ratio_thrsh
    else
      p_ratio_thrsh(1) = 1.0_REALK-10.0_REALK**(-6)
      p_ratio_thrsh(2) = 1.0_REALK+10.0_REALK**(-6)
    end if
    if (present(triangular)) then
      p_triangular = triangular
    else
      p_triangular = .false.
    end if
    if (p_triangular) then
      if (present(symmetric)) then
        p_symmetric = symmetric
      else
        p_symmetric = .true.
      end if
    else
      p_symmetric = .true.
    end if
    if (A%triangular) then
      ! matrix and array are in triangular format
      if (p_triangular) then
        if (A%symmetric.neqv.p_symmetric) then
          write(io_viewer,100) "matrix and array have different symmetries", &
                               A%symmetric, p_symmetric
          almost_equal = .false.
          return
        end if
        if (size(A%elms_alpha)/=size(values)) then
          write(io_viewer,110) "matrix and array have different sizes", &
                               size(A%elms_alpha), size(values)
          almost_equal = .false.
          return
        end if
        addr_elm = 0
        do icol = 1, A%num_col
          do irow = 1, icol
            addr_elm = addr_elm+1
            ! we check the ratio for absolutely greater values
            if (abs(values(addr_elm))>p_threshold) then
              mat_to_val = A%elms_alpha(addr_elm)/values(addr_elm)
              if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), values(addr_elm)
                almost_equal = .false.
              end if
            ! checks the difference for smaller values
            else
              if (abs(A%elms_alpha(addr_elm)-values(addr_elm))>p_threshold) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), values(addr_elm)
                almost_equal = .false.
              end if
            end if
          end do
        end do
      ! matrix is triangular, array is square
      else
        if (A%num_row*A%num_col/=size(values)) then
          write(io_viewer,110) "matrix and array have different sizes", &
                               A%num_row*A%num_col, size(values)
          almost_equal = .false.
          return
        end if
        addr_elm = 0
        do icol = 1, A%num_col
          do irow = 1, A%num_row
            addr_elm = addr_elm+1
            if (irow<=icol) then
              tmp_val = A%elms_alpha(irow+icol*(icol-1)/2)
            else
              if (A%symmetric) then
                tmp_val = A%elms_alpha(icol+irow*(irow-1)/2)
              else
                tmp_val = -A%elms_alpha(icol+irow*(irow-1)/2)
              end if
            end if
            ! we check the ratio for absolutely greater values
            if (abs(values(addr_elm))>p_threshold) then
              mat_to_val = tmp_val/values(addr_elm)
              if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
                write(io_viewer,120) irow, icol, &
                                     tmp_val, values(addr_elm)
                almost_equal = .false.
              end if
            ! checks the difference for smaller values
            else
              if (abs(tmp_val-values(addr_elm))>p_threshold) then
                write(io_viewer,120) irow, icol, &
                                     tmp_val, values(addr_elm)
                almost_equal = .false.
              end if
            end if
          end do
        end do
      end if
    else
      ! matrix is square, array is triangular
      if (p_triangular) then
        if (A%num_row/=A%num_col) then
          write(io_viewer,100) "matrix is a non-sqaure matrix"
          almost_equal = .false.
          return
        end if
        if (A%num_row*(A%num_row+1)/2/=size(values)) then
          write(io_viewer,110) "matrix and array have different sizes", &
                               A%num_row*(A%num_row+1)/2, size(values)
          almost_equal = .false.
          return
        end if
        addr_elm = 0
        do icol = 1, A%num_col
          do irow = 1, A%num_row
            addr_elm = addr_elm+1
            if (irow<=icol) then
              tmp_val = values(irow+icol*(icol-1)/2)
            else
              if (p_symmetric) then
                tmp_val = values(icol+irow*(irow-1)/2)
              else
                tmp_val = -values(icol+irow*(irow-1)/2)
              end if
            end if
            ! we check the ratio for absolutely greater values
            if (abs(tmp_val)>p_threshold) then
              mat_to_val = A%elms_alpha(addr_elm)/tmp_val
              if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), tmp_val
                almost_equal = .false.
              end if
            ! checks the difference for smaller values
            else
              if (abs(A%elms_alpha(addr_elm)-tmp_val)>p_threshold) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), tmp_val
                almost_equal = .false.
              end if
            end if
          end do
        end do
      ! matrix and array are in square format
      else
        if (size(A%elms_alpha)/=size(values)) then
          write(io_viewer,110) "matrix and array have different sizes", &
                               size(A%elms_alpha), size(values)
          almost_equal = .false.
          return
        end if
        addr_elm = 0
        do icol = 1, A%num_col
          do irow = 1, A%num_row
            addr_elm = addr_elm+1
            ! we check the ratio for absolutely greater values
            if (abs(values(addr_elm))>p_threshold) then
              mat_to_val = A%elms_alpha(addr_elm)/values(addr_elm)
              if (mat_to_val<p_ratio_thrsh(1) .or. mat_to_val>p_ratio_thrsh(2)) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), values(addr_elm)
                almost_equal = .false.
              end if
            ! checks the difference for smaller values
            else
              if (abs(A%elms_alpha(addr_elm)-values(addr_elm))>p_threshold) then
                write(io_viewer,120) irow, icol, &
                                     A%elms_alpha(addr_elm), values(addr_elm)
                almost_equal = .false.
              end if
            end if
          end do
        end do
      end if
    end if
100 format("MatArrayAlmostEqual>> ",A,L4,",",L4)
110 format("MatArrayAlmostEqual>> ",A,I8,",",I8)
120 format("MatArrayAlmostEqual>> element(",I4,",",I4,")",Es16.8," (Mat),",Es16.8," (array)")
  end subroutine MatArrayAlmostEqual

  !> \brief test suite of matrix module
  !> \author Bin Gao
  !> \date 2012-01-20
  !> \param test_failed indicates if tests failed
  !> \param io_viewer is the logical unit number of the viewer
  !> \param level_print is the level of print
  !> \param threshold is the threshold for comparison
  subroutine MatTestSuite(test_failed, io_viewer, level_print, threshold)
    logical, intent(out) :: test_failed
    integer, intent(in) :: io_viewer
    integer, intent(in) :: level_print
    real(REALK), optional, intent(in) :: threshold
! number of blocks of matrix
#if !defined(NUM_BLOCK)
#define NUM_BLOCK 2
#endif
! number of rows in a block
#if !defined(NUM_BLOCK_ROW)
#define NUM_BLOCK_ROW 4
#endif
    integer, parameter :: NUM_MAT_ROW = NUM_BLOCK*NUM_BLOCK_ROW
                                                  !number of rows of matrix
    integer, parameter :: DIM_TR_MAT = NUM_MAT_ROW*(NUM_MAT_ROW+1)/2
                                                  !dimension of matrix in triangular format
    integer, parameter :: DIM_SQ_MAT = NUM_MAT_ROW*NUM_MAT_ROW
                                                  !dimension of matrix in square format
    real(REALK) p_threshold                       !threshold for comparison (private)
    real(REALK), allocatable :: ref_tr_sym(:)     !ref. elements of symmetric matrix in triangular format
    real(REALK), allocatable :: ref_tr_anti(:)    !ref. elements of anti-symmetric matrix in triangular format
    real(REALK), allocatable :: ref_sq_sym(:,:)   !ref. elements of symmetric matrix in square format
    real(REALK), allocatable :: ref_sq_anti(:,:)  !ref. elements of anti-symmetric matrix in square format
    integer ierr                                  !error information
    integer irow, icol                            !incremental recorders over rows and columns
    integer addr_elm                              !address of element
    real(REALK), allocatable :: tr_sym_elms(:)    !element storage of symmetric matrix in triangular format
    real(REALK), allocatable :: tr_anti_elms(:)   !element storage of anti-symmetric matrix in triangular format
    real(REALK), allocatable :: sq_sym_elms(:)    !element storage of symmetric matrix in square format
    real(REALK), allocatable :: sq_anti_elms(:)   !element storage of anti-symmetric matrix in square format
    type(matrix) TrSym                            !symmetric matrix in triangular format
    type(matrix) TrAnti                           !anti-symmetric matrix in triangular format
    type(matrix) SqSym                            !symmetric matrix in square format
    type(matrix) SqAnti                           !anti-symmetric matrix in square format
    integer min_row_idx, max_row_idx              !minimum and maximum indices of rows
    integer min_col_idx, max_col_idx              !minimum and maximum indices of columns
    real(REALK) trace_tr_ss                       !trace of product of symmetric-symmetric matrices
    real(REALK) trace_tr_sa                       !trace of product of symmetric-anti-symmetric matrices
    real(REALK) trace_tr_as                       !trace of product of anti-symmetric-symmetric matrices
    real(REALK) trace_tr_aa                       !trace of product of anti-symmetric-anti-symmetric matrices
    real(REALK) trace_sq_ss                       !trace of product of symmetric-symmetric matrices
    real(REALK) trace_sq_sa                       !trace of product of symmetric-anti-symmetric matrices
    real(REALK) trace_sq_as                       !trace of product of anti-symmetric-symmetric matrices
    real(REALK) trace_sq_aa                       !trace of product of anti-symmetric-anti-symmetric matrices
    real(REALK) ref_trace_ss                      !ref. trace of product of symmetric-symmetric matrices
    real(REALK) ref_trace_aa                      !ref. trace of product of anti-symmetric-anti-symmetric matrices
    test_failed = .false.
    if (present(threshold)) then
      p_threshold = threshold
    else
      p_threshold = 10.0_REALK**(-8)
    end if
    ! sets up referenced elements
    allocate(ref_tr_sym(DIM_TR_MAT), stat=ierr)
    if (ierr/=0) then
      test_failed = .true.
      return
    end if
    allocate(ref_tr_anti(DIM_TR_MAT), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      test_failed = .true.
      return
    end if
    allocate(ref_sq_sym(NUM_MAT_ROW,NUM_MAT_ROW), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      test_failed = .true.
      return
    end if
    allocate(ref_sq_anti(NUM_MAT_ROW,NUM_MAT_ROW), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      test_failed = .true.
      return
    end if
    call random_number(ref_tr_sym)
    ! we uses the upper and diagonal parts in triangular format
    addr_elm = 0
    do icol = 1, NUM_MAT_ROW
      do irow = 1, icol-1
        addr_elm = addr_elm+1
        ref_tr_sym(addr_elm) = 0.1_REALK*ref_tr_sym(addr_elm)
        ref_tr_anti(addr_elm) = ref_tr_sym(addr_elm)
        ref_sq_sym(irow,icol) = ref_tr_sym(addr_elm)
        ref_sq_sym(icol,irow) = ref_tr_sym(addr_elm)
        ref_sq_anti(irow,icol) = ref_tr_sym(addr_elm)
        ref_sq_anti(icol,irow) = -ref_tr_sym(addr_elm)
      end do
      addr_elm = addr_elm+1
      ref_tr_sym(addr_elm) = 0.1_REALK*ref_tr_sym(addr_elm)
      ref_tr_anti(addr_elm) = 0.0_REALK
      ref_sq_sym(icol,icol) = ref_tr_sym(addr_elm)
      ref_sq_anti(icol,icol) = 0.0_REALK
    end do
    ! sets up test matrices
    allocate(tr_sym_elms(DIM_TR_MAT), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      test_failed = .true.
      return
    end if
    tr_sym_elms = 0.0_REALK
    call MatAssociate(work_alpha=tr_sym_elms, num_row=NUM_MAT_ROW, A=TrSym, &
                      info_mat=ierr, triangular=.true., symmetric=.true.)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      test_failed = .true.
      return
    end if
    allocate(tr_anti_elms(DIM_TR_MAT), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      test_failed = .true.
      return
    end if
    tr_anti_elms = 0.0_REALK
    call MatAssociate(work_alpha=tr_anti_elms, num_row=NUM_MAT_ROW, A=TrAnti, &
                      info_mat=ierr, triangular=.true., symmetric=.false.)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      deallocate(tr_anti_elms)
      test_failed = .true.
      return
    end if
    allocate(sq_sym_elms(DIM_SQ_MAT), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      deallocate(tr_anti_elms)
      test_failed = .true.
      return
    end if
    sq_sym_elms = 0.0_REALK
    call MatAssociate(work_alpha=sq_sym_elms, num_row=NUM_MAT_ROW, &
                      A=SqSym, info_mat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      deallocate(tr_anti_elms)
      deallocate(sq_sym_elms)
      test_failed = .true.
      return
    end if
    allocate(sq_anti_elms(DIM_SQ_MAT), stat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      deallocate(tr_anti_elms)
      deallocate(sq_sym_elms)
      test_failed = .true.
      return
    end if
    sq_anti_elms = 0.0_REALK
    call MatAssociate(work_alpha=sq_anti_elms, num_row=NUM_MAT_ROW, A=SqAnti, &
                      info_mat=ierr)
    if (ierr/=0) then
      deallocate(ref_tr_sym)
      deallocate(ref_tr_anti)
      deallocate(ref_sq_sym)
      deallocate(ref_sq_anti)
      deallocate(tr_sym_elms)
      deallocate(tr_anti_elms)
      deallocate(sq_sym_elms)
      deallocate(sq_anti_elms)
      test_failed = .true.
      return
    end if
    ! assigns elements of matrices
    do icol = 1, NUM_BLOCK
      ! sets the minimum and maximum of indices of columns of matrices
      min_col_idx = (icol-1)*NUM_BLOCK_ROW+1
      max_col_idx = icol*NUM_BLOCK_ROW
      do irow = 1, icol
        ! sets the minimum and maximum of indices of rows of matrices
        min_row_idx = (irow-1)*NUM_BLOCK_ROW+1
        max_row_idx = irow*NUM_BLOCK_ROW
        ! symmetric matrix in triangular format
        call MatSetValues(A=TrSym,                                                  &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trans=.false.)
        ! anti-symmetric matrix in triangular format
        call MatSetValues(A=TrAnti,                                                  &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trans=.false.)
        ! symmetric matrix in square format
        call MatSetValues(A=SqSym,                                                  &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trans=.false.)
        if (irow/=icol)                                                               &
          call MatSetValues(A=SqSym,                                                  &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                  values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trans=.true.)
        ! anti-symmetric matrix in square format
        call MatSetValues(A=SqAnti,                                                  &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trans=.false.)
        if (irow/=icol)                                                                 &
          call MatSetValues(A=SqAnti,                                                   &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                     &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                     &
                  values=-ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trans=.true.)
      end do
    end do
    ! visualizes the information of matrices
    if (level_print>=10) then
      call MatView(A=TrSym, io_viewer=io_viewer)
      call MatView(A=TrAnti, io_viewer=io_viewer)
      call MatView(A=SqSym, io_viewer=io_viewer)
      call MatView(A=SqAnti, io_viewer=io_viewer)
    end if
    ! checks results
    addr_elm = 0
    do icol = 1, NUM_MAT_ROW
      do irow = 1, icol-1
        addr_elm = addr_elm+1
        if (tr_sym_elms(addr_elm)/=ref_tr_sym(addr_elm)) then
          write(io_viewer,100) "TR-SYMM", irow, icol, &
                               tr_sym_elms(addr_elm), &
                               ref_tr_sym(addr_elm)
          test_failed = .true.
        end if
        if (tr_anti_elms(addr_elm)/=ref_tr_anti(addr_elm)) then
          write(io_viewer,100) "TR-ANTI", irow, icol,  &
                               tr_anti_elms(addr_elm), &
                               ref_tr_anti(addr_elm)
          test_failed = .true.
        end if
        if (sq_sym_elms(irow+(icol-1)*NUM_MAT_ROW)/=ref_sq_sym(irow,icol)) then
          write(io_viewer,100) "SQ-SYMM", irow, icol,                  &
                               sq_sym_elms(irow+(icol-1)*NUM_MAT_ROW), &
                               ref_sq_sym(irow,icol)
          test_failed = .true.
        end if
        if (sq_sym_elms(icol+(irow-1)*NUM_MAT_ROW)/=ref_sq_sym(icol,irow)) then
          write(io_viewer,100) "SQ-SYMM", icol, irow,                  &
                               sq_sym_elms(icol+(irow-1)*NUM_MAT_ROW), &
                               ref_sq_sym(icol,irow)
          test_failed = .true.
        end if
        if (sq_anti_elms(irow+(icol-1)*NUM_MAT_ROW)/=ref_sq_anti(irow,icol)) then
          write(io_viewer,100) "SQ-ANTI", irow, icol,                   &
                               sq_anti_elms(irow+(icol-1)*NUM_MAT_ROW), &
                               ref_sq_anti(irow,icol)
          test_failed = .true.
        end if
        if (sq_anti_elms(icol+(irow-1)*NUM_MAT_ROW)/=ref_sq_anti(icol,irow)) then
          write(io_viewer,100) "SQ-ANTI", icol, irow,                   &
                               sq_anti_elms(icol+(irow-1)*NUM_MAT_ROW), &
                               ref_sq_anti(icol,irow)
          test_failed = .true.
        end if
      end do
      addr_elm = addr_elm+1
      if (tr_sym_elms(addr_elm)/=ref_tr_sym(addr_elm)) then
        write(io_viewer,100) "TR-SYMM", icol, icol, &
                             tr_sym_elms(addr_elm), &
                             ref_tr_sym(addr_elm)
        test_failed = .true.
      end if
      if (tr_anti_elms(addr_elm)/=ref_tr_anti(addr_elm)) then
        write(io_viewer,100) "TR-ANTI", icol, icol,  &
                             tr_anti_elms(addr_elm), &
                             ref_tr_anti(addr_elm)
        test_failed = .true.
      end if
      if (sq_sym_elms(icol+(icol-1)*NUM_MAT_ROW)/=ref_sq_sym(icol,icol)) then
        write(io_viewer,100) "SQ-SYMM", icol, icol,                  &
                             sq_sym_elms(icol+(icol-1)*NUM_MAT_ROW), &
                             ref_sq_sym(icol,icol)
        test_failed = .true.
      end if
      if (sq_anti_elms(icol+(icol-1)*NUM_MAT_ROW)/=ref_sq_anti(icol,icol)) then
        write(io_viewer,100) "SQ-ANTI", icol, icol,                   &
                             sq_anti_elms(icol+(icol-1)*NUM_MAT_ROW), &
                             ref_sq_anti(icol,icol)
        test_failed = .true.
      end if
    end do
    ! evaluates expectation values
    trace_tr_ss = 0.0_REALK
    trace_tr_sa = 0.0_REALK
    trace_tr_as = 0.0_REALK
    trace_tr_aa = 0.0_REALK
    trace_sq_ss = 0.0_REALK
    trace_sq_sa = 0.0_REALK
    trace_sq_as = 0.0_REALK
    trace_sq_aa = 0.0_REALK
    do icol = 1, NUM_BLOCK
      ! sets the minimum and maximum of indices of columns of matrices
      min_col_idx = (icol-1)*NUM_BLOCK_ROW+1
      max_col_idx = icol*NUM_BLOCK_ROW
      do irow = 1, icol
        ! sets the minimum and maximum of indices of rows of matrices
        min_row_idx = (irow-1)*NUM_BLOCK_ROW+1
        max_row_idx = irow*NUM_BLOCK_ROW
        ! symmetric-symmetric matrix product in triangular format
        call MatMultBlockedTrace(A=TrSym,                                           &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_tr_ss, trans=.false.)
        if (irow/=icol)                                                               &
          call MatMultBlockedTrace(A=TrSym,                                           &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                  values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_tr_ss, trans=.true.)
        ! symmetric-anti-symmetric matrix product in triangular format
        call MatMultBlockedTrace(A=TrSym,                                            &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_tr_sa, trans=.false.)
        if (irow/=icol)                                                                 &
          call MatMultBlockedTrace(A=TrSym,                                             &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                     &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                     &
                  values=-ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_tr_sa, trans=.true.)
        ! anti-symmetric-symmetric matrix product in triangular format
        call MatMultBlockedTrace(A=TrAnti,                                          &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_tr_as, trans=.false.)
        if (irow/=icol)                                                               &
          call MatMultBlockedTrace(A=TrAnti,                                          &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                  values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_tr_as, trans=.true.)
        ! anti-symmetric-anti-symmetric matrix product in triangular format
        call MatMultBlockedTrace(A=TrAnti,                                           &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_tr_aa, trans=.false.)
        if (irow/=icol)                                                                 &
          call MatMultBlockedTrace(A=TrAnti,                                            &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                     &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                     &
                  values=-ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_tr_aa, trans=.true.)
        ! symmetric-symmetric matrix product in square format
        call MatMultBlockedTrace(A=SqSym,                                           &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_sq_ss, trans=.false.)
        if (irow/=icol)                                                               &
          call MatMultBlockedTrace(A=SqSym,                                           &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                  values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_sq_ss, trans=.true.)
        ! symmetric-anti-symmetric matrix product in square format
        call MatMultBlockedTrace(A=SqSym,                                            &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_sq_sa, trans=.false.)
        if (irow/=icol)                                                                 &
          call MatMultBlockedTrace(A=SqSym,                                             &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                     &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                     &
                  values=-ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_sq_sa, trans=.true.)
        ! anti-symmetric-symmetric matrix product in square format
        call MatMultBlockedTrace(A=SqAnti,                                          &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_sq_as, trans=.false.)
        if (irow/=icol)                                                               &
          call MatMultBlockedTrace(A=SqAnti,                                          &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                   &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                   &
                  values=ref_sq_sym(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_sq_as, trans=.true.)
        ! anti-symmetric-anti-symmetric matrix product in square format
        call MatMultBlockedTrace(A=SqAnti,                                           &
                min_row_idx=min_row_idx, max_row_idx=max_row_idx,                    &
                min_col_idx=min_col_idx, max_col_idx=max_col_idx,                    &
                values=ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                trace=trace_sq_aa, trans=.false.)
        if (irow/=icol)                                                                 &
          call MatMultBlockedTrace(A=SqAnti,                                            &
                  min_row_idx=min_row_idx, max_row_idx=max_row_idx,                     &
                  min_col_idx=min_col_idx, max_col_idx=max_col_idx,                     &
                  values=-ref_sq_anti(min_row_idx:max_row_idx,min_col_idx:max_col_idx), &
                  trace=trace_sq_aa, trans=.true.)
      end do
    end do
    ! referenced expectation values
    ref_trace_ss = 0.0_REALK
    ref_trace_aa = 0.0_REALK
    do irow = 1, NUM_MAT_ROW
      do icol = 1, NUM_MAT_ROW
        ref_trace_ss = ref_trace_ss+ref_sq_sym(icol,irow)*ref_sq_sym(irow,icol)
        ref_trace_aa = ref_trace_aa+ref_sq_anti(icol,irow)*ref_sq_anti(irow,icol)
      end do
    end do
    ! closed-shell, then additional factor 2
    ref_trace_ss = 2.0_REALK*ref_trace_ss
    ref_trace_aa = 2.0_REALK*ref_trace_aa
    ! checks expectation values
    if (abs(trace_tr_ss-ref_trace_ss)>p_threshold) then
      write(io_viewer,110) "TRACE-TR-SS", trace_tr_ss, ref_trace_ss
      test_failed = .true.
    end if
    if (abs(trace_tr_sa)>p_threshold) then
      write(io_viewer,110) "TRACE-TR-SA", trace_tr_sa, 0.0_REALK
      test_failed = .true.
    end if
    if (abs(trace_tr_as)>p_threshold) then
      write(io_viewer,110) "TRACE-TR-AS", trace_tr_as, 0.0_REALK
      test_failed = .true.
    end if
    if (abs(trace_tr_aa-ref_trace_aa)>p_threshold) then
      write(io_viewer,110) "TRACE-TR-AA", trace_tr_aa, ref_trace_aa
      test_failed = .true.
    end if
    if (abs(trace_sq_ss-ref_trace_ss)>p_threshold) then
      write(io_viewer,110) "TRACE-SQ-SS", trace_sq_ss, ref_trace_ss
      test_failed = .true.
    end if
    if (abs(trace_sq_sa)>p_threshold) then
      write(io_viewer,110) "TRACE-SQ-SA", trace_sq_sa, 0.0_REALK
      test_failed = .true.
    end if
    if (abs(trace_sq_as)>p_threshold) then
      write(io_viewer,110) "TRACE-SQ-AS", trace_sq_as, 0.0_REALK
      test_failed = .true.
    end if
    if (abs(trace_sq_aa-ref_trace_aa)>p_threshold) then
      write(io_viewer,110) "TRACE-SQ-AA", trace_sq_aa, ref_trace_aa
      test_failed = .true.
    end if
    ! frees space
    call MatNullify(A=TrSym)
    call MatNullify(A=TrAnti)
    call MatNullify(A=SqSym)
    call MatNullify(A=SqAnti)
    deallocate(tr_sym_elms)
    deallocate(tr_anti_elms)
    deallocate(sq_sym_elms)
    deallocate(sq_anti_elms)
    deallocate(ref_tr_sym)
    deallocate(ref_tr_anti)
    deallocate(ref_sq_sym)
    deallocate(ref_sq_anti)
100 format("MatTestSuite>> ",A," row",I4," col",I4,2Es14.6)
110 format("MatTestSuite>> ",A,2Es14.6)
  end subroutine MatTestSuite
#endif

end module gen1int_matrix
