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

#include "global.h"

module hamiltonian
use external_pot
use poisson
use xc

implicit none

type hamiltonian_type
  integer :: spin_channels ! How to handle spin (duplicated in states_type)
  integer :: nspin
  integer :: ispin

  integer :: reltype ! type of relativistic correction to use

  real(r8), pointer :: Vpsl(:)     ! the external potential
  real(r8), pointer :: Vhxc(:,:)   ! xc potential + hartree potential

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec, epot

  ! System under the independent particle approximation, or not.
  logical :: ip_app

  type(epot_type) :: ep  ! handles the external potential
  type(xc_type)   :: xc  ! handles the xc potential

  ! gauge
  integer :: gauge ! in which gauge shall we work in
                   ! 1 = length gauge
                   ! 2 = velocity gauge

  ! absorbing boundaries
  integer  :: ab         ! do we have absorbing boundaries?
  real(r8) :: ab_width   ! width of the absorbing boundary
  real(r8) :: ab_height  ! height of the absorbing boundary
  real(r8), pointer :: ab_pot(:) ! where we store the ab potential

  ! Spectral range.
  real(r8) :: spectral_middle_point
  real(r8) :: spectral_half_span

  ! Kinetic Cutoff
  real(r8) :: cutoff
end type hamiltonian_type

integer, parameter :: NOREL      = 0, &
                      SPIN_ORBIT = 1, &
                      APP_ZORA   = 2, &
                      ZORA       = 3

integer, parameter :: NO_ABSORBING        = 0, &
                      IMAGINARY_ABSORBING = 1, &
                      MASK_ABSORBING      = 2
contains

subroutine hamiltonian_init(h, sys)
  type(hamiltonian_type), intent(out) :: h
  type(system_type), intent(inout) :: sys

  integer :: i, j, n
  real(r8) :: d(3), r, x(3)

  call push_sub('hamiltonian_init')

  ! Duplicate this two variables
  h%ispin         = sys%st%ispin
  h%spin_channels = sys%st%spin_channels
  h%nspin         = sys%st%nspin

  ! allocate potentials and density of the cores
  ! In the case of spinors, vxc_11 = h%vxc(:, 1), vxc_22 = h%vxc(:, 2), Re(vxc_12) = h%vxc(:. 3);
  ! Im(vxc_12) = h%vxc(:, 4)
  allocate(h%Vpsl(sys%m%np), h%Vhxc(sys%m%np, sys%st%nspin))
  h%vpsl = M_ZERO
  h%Vhxc = M_ZERO

  call epot_init(h%ep, sys)

  call oct_parse_int("RelativisticCorrection", NOREL, h%reltype)
#ifdef COMPLEX_WFNS
  select case(h%reltype)
    case(NOREL);      message(1) = 'Info: Rel. correction: No relativistic corrections.'
    case(SPIN_ORBIT); message(1) = 'Info: Rel. corrections: Spin-Orbit.'
    case(APP_ZORA);   message(1) = 'Info: Rel. correction: Approximated ZORA.'
    case(ZORA);       message(1) = 'Info: Rel. correction: ZORA.'
    case default
      message(1) = "Relativistic corrections must be: 0 (none),"
      message(2) = "                                  1 (spin-orbit coupling)"
      message(3) = "                                  2 (approximated ZORA)"
      message(4) = "                                  3 (ZORA)"
      call write_fatal(4)
  end select
  call write_info(1)
  ! This is temporary...
  if(h%reltype > SPIN_ORBIT) then
    message(1) = 'Error: ZORA corrections not working yet. Visit us soon.'
    call write_fatal(1)
  endif
#else
  if(h%reltype .ne. NOREL) then
    message(1) = "Cannot apply relativistic corrections with an executable compiled"
    message(2) = "for real wavefunctions."
    call write_fatal(2)
  endif
#endif

  ! Should we treat the particles as independent?
  call oct_parse_logical("NonInteractingElectrons", .false., h%ip_app)
  if(h%ip_app) then
    message(1) = 'Info: Treating the electrons as non-interacting'
    call write_info(1)
  else
    ! initilize hartree and xc modules
    call poisson_init(sys%m)
    call xc_init(h%xc, sys%nlcc)
    message(1) = "Info: Exchange and correlation"
    call write_info(1)
    if(conf%verbose > 20) call xc_write_info(h%xc, stdout)
  end if

  ! gauge
  call oct_parse_int("TDGauge", 1, h%gauge)
  if (h%gauge < 1 .or. h%gauge > 2) then
    write(message(1), '(a,i6,a)') "Input: '", h%gauge, "' is not a valid TDGauge"
    message(2) = 'Accepted values are:'
    message(3) = '   1 = length gauge'
    message(4) = '   2 = velocity gauge'
    call write_fatal(4)
  end if

  ! absorbing boundaries
  call oct_parse_int("AbsorbingBoundaries", NO_ABSORBING, h%ab)
  nullify(h%ab_pot)

  absorbing_boundaries: if(h%ab.ne.NO_ABSORBING) then
    call oct_parse_double("ABWidth", 4._r8/units_inp%length%factor, h%ab_width)
    h%ab_width  = h%ab_width * units_inp%length%factor
    if(h%ab == 1) then
      call oct_parse_double("ABHeight", -0.2_r8/units_inp%energy%factor, h%ab_height)
      h%ab_height = h%ab_height * units_inp%energy%factor
    else
      h%ab_height = M_ONE
    end if
    
    ! generate boundary potential...
    allocate(h%ab_pot(sys%m%np))
    h%ab_pot = M_ZERO
    do i = 1, sys%m%np
         call mesh_inborder(sys%m, i, n, d, h%ab_width)
         if(n>0) then
           do j = 1, n
              h%ab_pot(i) = h%ab_pot(i) + h%ab_height * sin(d(j)*M_PI/(M_TWO*h%ab_width))**2
           enddo
         endif
         if(abs(h%ab_pot(i)) > abs(h%ab_height)) h%ab_pot(i) = h%ab_height
    enddo

  end if absorbing_boundaries
  
  ! Cutoff applied to the kinetic term.
  ! it is used *both* in the calculation of the derivatives and in the split operator method
  call oct_parse_double("KineticCutoff", -M_ONE, h%cutoff)
  if(h%cutoff > M_ZERO) then
    h%cutoff = h%cutoff * units_inp%energy%factor
    write(message(1),'(a,f7.2,a)') 'Info: The kinetic operator will have a cutoff of',&
                                  h%cutoff/units_out%energy%factor, units_out%energy%abbrev
    write(message(2),'(a)')        '      (only if DerivativesSpace = 1 is set)'
    call write_info(2)
  endif

  call pop_sub()
end subroutine hamiltonian_init

subroutine hamiltonian_end(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(in) :: sys
  
  call push_sub('hamiltonian_end')

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhxc)
    nullify(h%Vpsl, h%Vhxc)
  end if

  call epot_end(h%ep, sys)
  if(.not.h%ip_app) then
    call poisson_end()
    call xc_end(h%xc)
  end if

  if(associated(h%ab_pot)) then
    deallocate(h%ab_pot); nullify(h%ab_pot)
  end if

  call pop_sub()
end subroutine hamiltonian_end

! This subroutine calculates the total energy of the system. Basically, it
! adds up the KS eigenvalues, and then it substracts the whatever double
! counts exist (see TDDFT theory for details).
subroutine hamiltonian_energy(h, st, eii, iunit, reduce)
  type(hamiltonian_type), intent(inout) :: h
  type(states_type), intent(in) :: st
  real(r8), intent(in) :: eii
  integer, intent(in) :: iunit
  logical, intent(in), optional :: reduce

  integer :: ik
  real(r8) :: s, e
#ifdef HAVE_MPI
  integer :: ierr
#endif 

  call push_sub('hamiltonian_energy')

  e = 0
  do ik = 1, st%nik
    e = e + st%kweights(ik) * sum(st%occ(st%st_start:st%st_end, ik)* &
         st%eigenval(st%st_start:st%st_end, ik))
  end do
#ifdef HAVE_MPI
  if(present(reduce)) then
    if(reduce) then
      call MPI_ALLREDUCE(e, s, 1, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      e = s
    end if
  end if
#endif
  
  h%etot = e + eii + h%epot + h%ex + h%ec
  
  if(iunit>0) then
    write(message(1), '(6x,a, f15.8)')'Ion-ion     = ', eii    / units_out%energy%factor
    write(message(2), '(6x,a, f15.8)')'Eigenvalues = ', e      / units_out%energy%factor
    write(message(3), '(6x,a, f15.8)')'Potentials  = ', h%epot / units_out%energy%factor
    write(message(4), '(6x,a, f15.8)')'Exchange    = ', h%ex   / units_out%energy%factor
    write(message(5), '(6x,a, f15.8)')'Correlation = ', h%ec   / units_out%energy%factor
    write(message(6), '(6x,a, f15.8)')'Total       = ', h%etot / units_out%energy%factor
    call write_info(6, iunit)
  end if

  call pop_sub()
end subroutine hamiltonian_energy

subroutine hamiltonian_span(h, delta, emin)
  type(hamiltonian_type), intent(inout) :: h
  real(r8), intent(in) :: delta, emin

  call push_sub('hamiltonian_span')

  h%spectral_middle_point = ((M_Pi**2/(2*delta**2)) + emin)/2._r8
  h%spectral_half_span    = ((M_Pi**2/(2*delta**2)) - emin)/2._r8

  call pop_sub()
end subroutine hamiltonian_span

subroutine hamiltonian_output(h, m, dir, outp)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  character(len=*), intent(IN) :: dir
  type(output_type), intent(IN) :: outp

  integer :: is
  character(len=80) :: fname  
  real(r8) :: u

  call push_sub('hamiltonian_output')

  u = units_out%energy%factor
  if(outp%what(output_potential)) then
    call doutput_function(outp%how, dir, "v0", m, h%Vpsl, u)

    if(h%ep%classic_pot > 0) then
      call doutput_function(outp%how, dir, "vc", m, h%ep%Vclassic, u)
    end if

    if(.not.h%ip_app) then
      do is = 1, min(h%ispin, 2)
        write(fname, '(a,i1)') 'vhxc-', is
        call doutput_function(outp%how, dir, fname, m, h%Vhxc(:, is), u)
      end do
    end if
  end if

  call pop_sub()
end subroutine hamiltonian_output

#include "undef.F90"
#include "real.F90"
#include "h_inc.F90"
#include "h_xc_OEP.F90"

#include "undef.F90"
#include "complex.F90"
#include "h_inc.F90"
#include "h_xc_OEP.F90"

#if defined(COMPLEX_WFNS)
#include "h_so.F90"
#endif

end module hamiltonian
