!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "config_F90.h"

module hamiltonian
use global
use liboct
use spline
use fft
use units
use system
use specie
use hartree
use xc
use lasers

implicit none

type hamiltonian_type
  integer :: ispin ! How to handle spin (duplicated in states_type)
  integer :: np  ! number of points (duplicated in mesh)

  logical :: soc ! spin-orbit coupling for anyone?

  integer :: vpsl_space           ! How should the local potential be calculated
  integer :: vnl_space            ! How should the nl    potential be calculated
  integer :: nextra               ! extra points for the interpolation method(s)
  real(r8), pointer :: Vpsl(:)    ! the external potential
  real(r8), pointer :: Vhartree(:)! the hartree potential
  real(r8), pointer :: Vxc(:,:)   ! xc potential

  ! For non-collinear spin (ispin=3)
  real(r8), pointer :: dVxc_off(:)
  complex(r8), pointer :: zVxc_off(:)

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec, epot

  ! System under the independent particle approximation, or not.
  logical :: ip_app

  ! should we include the classical point charges
  logical :: classic_pot
  real(r8), pointer :: Vclassic(:)! potential created by classic point charges

  ! hartree potential structure
  type(hartree_type) :: hart
  type(xc_type) :: xc

end type hamiltonian_type

contains

subroutine hamiltonian_init(h, sys)
  type(hamiltonian_type), intent(out) :: h
  type(system_type), intent(inout) :: sys

  sub_name = 'hamiltonian_init'; call push_sub()

  ! Duplicate this two variables
  h%ispin = sys%st%ispin
  h%np  = sys%m%np

  ! allocate potentials and density of the cores
  allocate(h%Vpsl(h%np), h%Vhartree(h%np), h%Vxc(h%np, sys%st%nspin))
  h%Vhartree = 0._r8; h%Vxc = 0._r8
  if(h%ispin == 3) then
    allocate(h%R_FUNC(Vxc_off) (h%np))
    h%R_FUNC(Vxc_off) = R_TOTYPE(0._r8)
  else
    nullify(h%dVxc_off, h%zVxc_off)
  end if

  if(sys%ncatoms > 0) then
    call oct_parse_logical(C_string("ClassicPotential"), .false., h%classic_pot)
    if(h%classic_pot) allocate(h%Vclassic(h%np))
  end if

  ! should we calculate the local pseudopotentials in Fourier space?
  call oct_parse_int(C_string('LocalPotentialSpace'), 1, h%vpsl_space)
  if(h%vpsl_space < 0 .or. h%vpsl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vpsl_space, &
         "' is not a valid LocalPotentialSpace"
    message(2) = '(LocalPotentialSpace = 0 | 1)'
    call write_fatal(2)
  end if

  call oct_parse_int(C_string('NonLocalPotentialSpace'), 0, h%vnl_space)
  if(h%vnl_space < 0 .or. h%vnl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vnl_space, &
         "' is not a valid NonLocalPotentialSpace"
    message(2) = '(NonLocalPotentialSpace = 0 | 1)'
    call write_fatal(2)
  end if

  if(h%vpsl_space == 1) then
    call mesh_alloc_ffts(sys%m, 2)
    call specie_local_fourier_init(sys%nspecies, sys%specie, sys%m, sys%st%nlcc)
  end if

  if(h%vnl_space == 1) then
    call oct_parse_int(C_string('GridRefinement'), 3, h%nextra)
    if(h%nextra < 0) then
      write(message(1), '(a,i5,a)') "Input: '", h%nextra, &
           "' is not a valid GridRefinement"
      message(2) = '(GridRefinement >= 0)'
      call write_fatal(2)
    end if

    call specie_nl_fourier_init(sys%nspecies, sys%specie, sys%m, h%nextra)
  end if

#ifdef COMPLEX_WFNS
  call oct_parse_logical(C_string("NonInteractingElectrons"), .false., h%ip_app)
  if(h%ip_app .and. h%ispin.ne.3) then
    message(1) = "Spin-orbit coupling can only be included with SpinComponents = 3"
    call write_fatal(1)
  end if
#endif

  ! Should we treat the particles as independent?
  call oct_parse_logical(C_string("NonInteractingElectrons"), .false., h%ip_app)
  if(h%ip_app) then
    message(1) = 'Info: Treating the electrons as non-interacting'
    call write_info(1)
  else
    ! initilize hartree and xc modules
    call hartree_init(h%hart, sys%m)
    call xc_init(h%xc, sys%m, h%ispin)
    message(1) = "Info: Exchange and correlation"
    call write_info(1)
    if(conf%verbose > 20) call xc_write_info(h%xc, stdout)
  end if

  call pop_sub()
end subroutine hamiltonian_init

subroutine hamiltonian_end(h)
  type(hamiltonian_type) :: h

  sub_name = 'hamiltonian_end'; call push_sub()

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhartree, h%Vxc)
    nullify(h%Vpsl, h%Vhartree, h%Vxc)
  end if

  if(h%ispin == 3 .and. associated(h%dVxc_off)) then
    deallocate(h%dVxc_off); nullify(h%dVxc_off)
  end if
  if(h%ispin == 3 .and. associated(h%zVxc_off)) then
    deallocate(h%zVxc_off); nullify(h%zVxc_off)
  end if

  if(h%classic_pot .and. associated(h%Vclassic)) then
    deallocate(h%Vclassic); nullify(h%Vclassic)
  end if

  if(.not.h%ip_app) then
    call hartree_end(h%hart)
    call xc_end(h%xc)
  end if

  call pop_sub()
end subroutine hamiltonian_end

! This subroutine calculates the total energy of the system. Basically, it
! adds up the KS eigenvalues, and then it substracts the whatever double
! counts exist (see TDDFT theory for details).
subroutine hamiltonian_energy(h, sys, iunit, reduce)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys
  integer, intent(in) :: iunit
  logical, intent(in), optional :: reduce

  integer :: ik
  real(r8) :: s, e
#ifdef HAVE_MPI
  integer :: ierr
#endif 

  sub_name = 'hamiltonian_energy'; call push_sub()

  e = 0
  do ik = 1, sys%st%nik
    e = e + sys%st%kweights(ik) * sum(sys%st%occ(sys%st%st_start:sys%st%st_end, ik)* &
         sys%st%eigenval(sys%st%st_start:sys%st%st_end, ik))
  end do
#ifdef HAVE_MPI
  if(present(reduce) .and. reduce) then
    call MPI_ALLREDUCE(e, s, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    e = s
  end if
#endif
  
  h%etot = e + sys%eii + h%epot + h%ex + h%ec
  
  if(iunit>0) then
    write(message(1), '(6x,a, f15.8)')'Ion-ion     = ', sys%eii / units_out%energy%factor
    write(message(2), '(6x,a, f15.8)')'Eigenvalues = ', e       / units_out%energy%factor
    write(message(3), '(6x,a, f15.8)')'Potentials  = ', h%epot  / units_out%energy%factor
    write(message(4), '(6x,a, f15.8)')'Exchange    = ', h%ex    / units_out%energy%factor
    write(message(5), '(6x,a, f15.8)')'Correlation = ', h%ec    / units_out%energy%factor
    write(message(6), '(6x,a, f15.8)')'Total       = ', h%etot  / units_out%energy%factor
    call write_info(6, iunit)
  end if

  call pop_sub()
  return
end subroutine hamiltonian_energy

#include "h_external_pot.F90"

#include "undef.F90"
#include "real.F90"
#include "h_inc.F90"
#include "h_forces.F90"

#include "undef.F90"
#include "complex.F90"
#include "h_inc.F90"
#include "h_forces.F90"

end module hamiltonian
