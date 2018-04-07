subroutine gen1int_api_initialize(natoms, num_shells, coords, charges)
  use gen1int_api
  integer, intent(in) :: natoms
  integer, intent(in) :: num_shells
  real(REALK), intent(in) :: coords(3, natoms)
  real(REALK), intent(in) :: charges(natoms)
  call Gen1IntAPICreateEasy(natoms, num_shells, coords, charges)

end subroutine gen1int_api_initialize

subroutine gen1int_create_shell(spher_gto, idx_center, coord_center, &
  ang_num, num_prim, exponents, num_contr, contr_coef)
  use gen1int_api
  logical, intent(in) :: spher_gto
  integer, intent(in) :: idx_center
  real(REALK), intent(in) :: coord_center(3)
  integer, intent(in) :: ang_num
  integer, intent(in) :: num_prim
  real(REALK), intent(in) :: exponents(num_prim)
  integer, intent(in) :: num_contr
  real(REALK), intent(in) :: contr_coef(num_prim)
  call Gen1IntAPICreateShell(spher_gto, idx_center, coord_center, &
    ang_num, num_prim, exponents, num_contr, contr_coef)

end subroutine gen1int_create_shell

function gen1int_api_initialized()
  use gen1int_api, only: test_init => Gen1IntAPIInited
  logical :: gen1int_api_initialized

  write (*,*) "gen1int initialized: ", test_init()
  ! gen1int_api_initialized() = test_init()
end function gen1int_api_initialized


subroutine pe_set_potfile(potfileName)
  use pe_variables
  character(len=80) :: potfileName
  potfile = potfileName
end subroutine

subroutine pe_interface_init(nAtoms, coords, charges)
  use polarizable_embedding, only: pe_init
  use pe_variables
  !implicit none

  integer,intent(in) :: nAtoms
  real(dp), intent(in), optional :: charges(nAtoms)
  real(dp), intent(in), optional :: coords(3, nAtoms)

  !double precision :: parseArray
  ! set to 6 for stdout
  !double precision,dimension(no),intent(in) :: x
  ! lupri, ..., ...
  call pe_init(6, coords, charges)
end subroutine

subroutine print_shells()
  use gen1int_api
  call Gen1IntPrintShells()
end subroutine print_shells

subroutine pe_interface_energy(densmatrix, ndim, nnbas)
  use polarizable_embedding, only: pe_master
  integer, intent(in) :: ndim
  integer, intent(in) :: nnbas
  real(8), dimension(nnbas), intent(in) :: densmatrix

  call pe_master("print_energy", .true., ndim, 1, densmatrix)
  ! call pe_master(runtype, triang, ndim, nmats, denmats, fckmats, expvals)
end subroutine pe_interface_energy

subroutine pe_interface_fock(densmatrix, ndim, nnbas, fckmatrix, energy)
  use polarizable_embedding, only: pe_master
  use gen1int_api
  integer, intent(in) :: ndim
  integer, intent(in) :: nnbas
  real(8), dimension(nnbas), intent(in) :: densmatrix
  ! give storage for fckmatrix
  real(8), dimension(nnbas), intent(out) :: fckmatrix
  real(8), intent(out) :: energy
  real(8), dimension(1) :: energies

  call pe_master("full_fock", .true., ndim, 1, densmatrix, fckmatrix, energies)
  ! call pe_master("print_energy", .true., ndim, 1, densmatrix)
  energy = energies(1)
end subroutine pe_interface_fock

subroutine pe_interface_pol_energy(densmatrix, ndim, nnbas, energy)
  use polarizable_embedding, only: pe_master
  use gen1int_api
  integer, intent(in) :: ndim
  integer, intent(in) :: nnbas
  real(8), dimension(nnbas), intent(in) :: densmatrix
  real(8), intent(out) :: energy
  real(8), dimension(1) :: energies
  call pe_master(runtype="get_pol_energy", triang=.true., ndim=ndim, nmats=1, denmats=densmatrix, expvals=energies)
  energy = energies(1)
end subroutine pe_interface_pol_energy


subroutine pe_set_border_options(m_pe_border, m_rmin, type_flag)
  use pe_variables
  logical, intent(in) :: m_pe_border
  real(dp), intent(in) :: m_rmin
  integer, intent(in) :: type_flag

  if ( type_flag == 0 ) then
    border_type = "REMOVE"
  elseif ( type_flag == 1 ) then
    border_type = "REDIST"
  endif

  pe_border = m_pe_border
  Rmin = m_rmin
end subroutine pe_set_border_options

subroutine finalize_all()
  use polarizable_embedding
  use gen1int_api

  call pe_finalize()
  call Gen1IntAPIDestroy()

end subroutine finalize_all

! subroutine pe_interface_response(densmatrix, ndim, nnbas, fckmatrix, energy)
!   use polarizable_embedding, only: pe_master
!   use gen1int_api
!   integer, intent(in) :: ndim
!   integer, intent(in) :: nnbas
!   real(8), dimension(nnbas), intent(in) :: densmatrix
!   ! give storage for fckmatrix
!   real(8), dimension(nnbas), intent(out) :: fckmatrix
!   real(8), intent(out) :: energy
!
!   ! call pe_master("print_energy", .true., ndim, 1, densmatrix)
!   call pe_master("dynamic_response", .true., ndim, 1, densmatrix, fckmatrix)
! end subroutine pe_interface_response
