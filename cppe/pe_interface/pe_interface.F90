subroutine set_coord_nuclei(number_atoms, coordinates, charges)
  use gen1int
! #include "mxcent.h"
! #include "maxaqn.h"
! #include "ccom.h"
! #include "nuclei.h"
! #include "orgcom.h"
  integer, intent(in) :: number_atoms
  real(REALK), intent(in) :: coordinates(3,MXCENT)
  real(REALK), intent(in) :: charges(MXCENT)
  ! NUCDEP = number_atoms
  ! NATOMS = number_atoms
  ! CORD = coordinates
  ! CHARGE = charges
end subroutine set_coord_nuclei

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

subroutine gen1int_api_create(num_atom_types, natoms, coords, charges, num_sym_atom, ang_numbers, NBLCK, num_cgto, num_prim, &
  num_contr, exponents, ucontr_coefs, pure)
  use gen1int_api
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
  ! integer, intent(in) :: num_comp
  integer :: KATOM = 500
  integer :: KANG = 17
  integer :: KBLOCK = 18000
  integer :: KPRIM = 35
  integer, intent(in) :: num_atom_types
  integer, intent(in) :: natoms
  real(REALK), intent(in) :: coords(3,num_atom_type)
  real(REALK), intent(in) :: charges(num_atom_type)
  ! integer, intent(in) :: KATOM ! number of atoms
  integer, intent(in) :: num_sym_atom(KATOM)
  integer, intent(in) :: ang_numbers(KATOM,num_comp) ! at the moment, we are going to pass this as NBLCK
  integer, intent(in) :: NBLCK(KATOM,num_comp)
  ! integer, intent(in) :: KANG
  integer, intent(in) :: num_cgto(KANG,KATOM,num_comp)
  ! integer, intent(in) :: KBLOCK
  integer, intent(in) :: num_prim(KBLOCK,num_comp)
  integer, intent(in) :: num_contr(KBLOCK,num_comp)
  ! integer, intent(in) :: KPRIM
  real(REALK), intent(in) :: exponents(KPRIM,KBLOCK,num_comp)
  real(REALK), intent(in) :: ucontr_coefs(KPRIM,KPRIM,KBLOCK,num_comp)
  logical, intent(in) :: pure
  ! num_comp, num_atom_type, KATOM, num_sym_atom, &
  !                             ang_numbers, NBLCK, KANG, num_cgto, KBLOCK,   &
  !                             num_prim, num_contr, KPRIM, exponents, ucontr_coefs
  call Gen1IntAPICreate(1, num_atom_types, natoms, coords, charges, KATOM, num_sym_atom, ang_numbers, &
                        NBLCK, KANG, num_cgto, KBLOCK, num_prim, num_contr, KPRIM, &
                        exponents, ucontr_coefs, .TRUE., pure)
  ! write (*,*) num_sym_atom(1:num_atom_types)
  ! write (*,*) ang_numbers(1:num_atom_types,1)
  ! write (*,*) num_cgto(1:3,1:3,1)
  ! write (*,*) exponents(1:3,1,1)
  ! write (*,*) ucontr_coefs(1:3,1,1,1)
  ! write (*,*) NBLCK(1:3,1)

end subroutine gen1int_api_create

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

  ! call pe_master("print_energy", .true., ndim, 1, densmatrix)
  call pe_master("full_fock", .true., ndim, 1, densmatrix, fckmatrix, energies)
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

subroutine pe_interface_response(densmatrix, ndim, nnbas, fckmatrix, energy)
  use polarizable_embedding, only: pe_master
  use gen1int_api
  integer, intent(in) :: ndim
  integer, intent(in) :: nnbas
  real(8), dimension(nnbas), intent(in) :: densmatrix
  ! give storage for fckmatrix
  real(8), dimension(nnbas), intent(out) :: fckmatrix
  real(8), intent(out) :: energy

  ! call pe_master("print_energy", .true., ndim, 1, densmatrix)
  call pe_master("dynamic_response", .true., ndim, 1, densmatrix, fckmatrix)
end subroutine pe_interface_response
