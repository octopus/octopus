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
use output
use system
use specie
use hartree
use xc
use lasers

implicit none

type hamiltonian_type
  integer :: spin_channels ! How to handle spin (duplicated in states_type)
  integer :: nspin
  integer :: ispin

  integer :: np  ! number of points (duplicated in mesh)

  integer :: reltype ! type of relativistic correction to use

  integer :: vpsl_space           ! How should the local potential be calculated
  integer :: vnl_space            ! How should the nl    potential be calculated
  integer :: nextra               ! extra points for the interpolation method(s)
  real(r8), pointer :: Vpsl(:)    ! the external potential
  real(r8), pointer :: Vhartree(:)! the hartree potential
  real(r8), pointer :: Vxc(:,:)   ! xc potential

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec, epot

  ! System under the independent particle approximation, or not.
  logical :: ip_app

  ! should we include the classical point charges
  integer :: classic_pot
  real(r8), pointer :: Vclassic(:)! potential created by classic point charges

  ! hartree potential structure
  type(hartree_type) :: hart
  type(xc_type) :: xc

  ! gauge
  integer :: gauge ! in which gauge shall we work in
                   ! 1 = length gauge
                   ! 2 = velocity gauge

  ! lasers stuff
  integer :: no_lasers ! number of laser pulses used
  logical :: output_laser ! write laser field
  type(laser_type), pointer :: lasers(:)

  ! absorbing boundaries
  integer  :: ab         ! do we have absorbing boundaries?
  real(r8) :: ab_width   ! width of the absorbing boundary
  real(r8) :: ab_height  ! height of the absorbing boundary
  real(r8), pointer :: ab_pot(:) ! where we store the ab potential

  ! Spectral range.
  real(r8) :: spectral_middle_point
  real(r8) :: spectral_half_span
end type hamiltonian_type

integer, parameter :: NOREL      = 0, &
                      SPIN_ORBIT = 1, &
                      APP_ZORA   = 2, &
                      ZORA       = 3

contains

subroutine hamiltonian_init(h, sys)
  type(hamiltonian_type), intent(out) :: h
  type(system_type), intent(inout) :: sys

  integer :: i, j, dummy
  real(r8) :: d, r, x(3)

  sub_name = 'hamiltonian_init'; call push_sub()

  ! Duplicate this two variables
  h%ispin = sys%st%ispin
  h%spin_channels = sys%st%spin_channels
  h%nspin         = sys%st%nspin
  h%np  = sys%m%np

  ! allocate potentials and density of the cores
  allocate(h%Vpsl(h%np),              &
           h%Vhartree(h%np),          &
           h%Vxc(h%np, sys%st%nspin))
  ! In the case of spinors, vxc_11 = h%vxc(:, 1), vxc_22 = h%vxc(:, 2), Re(vxc_12) = h%vxc(:. 3);
  ! Im(vxc_12) = h%vxc(:, 4)
  h%vpsl = ZERO; h%vhartree = ZERO; h%vxc = ZERO

  if(sys%ncatoms > 0) then
    call oct_parse_int(C_string("ClassicPotential"), 0, h%classic_pot)
    if(h%classic_pot > 0) allocate(h%Vclassic(h%np))
  end if

  ! should we calculate the local pseudopotentials in Fourier space?
  call oct_parse_int(C_string('LocalPotentialSpace'), RECIPROCAL_SPACE, h%vpsl_space)
  select case(h%vpsl_space)
    case(RECIPROCAL_SPACE)
      message(1) = 'Info: Local Potential in Reciprocal Space.'
    case(REAL_SPACE)
      message(1) = 'Info: Local Potential in Real Space.'
    case default
      write(message(1), '(a,i5,a)') "Input: '", h%vpsl_space, &
           "' is not a valid LocalPotentialSpace"
      message(2) = '(LocalPotentialSpace = 0 | 1)'
      call write_fatal(2)
  end select
  call write_info(1)

  call oct_parse_int(C_string('NonLocalPotentialSpace'), REAL_SPACE, h%vnl_space)
  if(h%vnl_space < 0 .or. h%vnl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vnl_space, &
         "' is not a valid NonLocalPotentialSpace"
    message(2) = '(NonLocalPotentialSpace = 0 | 1)'
    call write_fatal(2)
  end if

  if(h%vpsl_space == RECIPROCAL_SPACE) then
    call mesh_alloc_ffts(sys%m, 2)
    call specie_local_fourier_init(sys%nspecies, sys%specie, sys%m, sys%st%nlcc)
  end if

  if(h%vnl_space == RECIPROCAL_SPACE) then
    call oct_parse_int(C_string('GridRefinement'), 3, h%nextra)
    if(h%nextra < 0) then
      write(message(1), '(a,i5,a)') "Input: '", h%nextra, &
           "' is not a valid GridRefinement"
      message(2) = '(GridRefinement >= 0)'
      call write_fatal(2)
    end if

    call specie_nl_fourier_init(sys%nspecies, sys%specie, sys%m, h%nextra)
  end if

  call oct_parse_int(C_string("RelativisticCorrection"), NOREL, h%reltype)
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
  call oct_parse_logical(C_string("NonInteractingElectrons"), .false., h%ip_app)
  if(h%ip_app) then
    message(1) = 'Info: Treating the electrons as non-interacting'
    call write_info(1)
  else
    ! initilize hartree and xc modules
    call hartree_init(h%hart, sys%m)
    call xc_init(h%xc)
    message(1) = "Info: Exchange and correlation"
    call write_info(1)
    if(conf%verbose > 20) call xc_write_info(h%xc, stdout)
  end if

  ! Temporarily spinors and 4-component densities do not work unless LDA.
  if(h%ispin == SPINORS .and. (h%xc%x_family > 1 .or. h%xc%c_family > 1)) then
    message(1) = 'Currently the spinors only work well with LDA functionals.'
    call write_fatal(1)
  endif

  ! gauge
  call oct_parse_int(C_string("TDGauge"), 1, h%gauge)
  if (h%gauge < 1 .or. h%gauge > 2) then
    write(message(1), '(a,i6,a)') "Input: '", h%gauge, "' is not a valid TDGauge"
    message(2) = 'Accepted values are:'
    message(3) = '   1 = length gauge'
    message(4) = '   2 = velocity gauge'
    call write_fatal(4)
  end if

  ! lasers
  call laser_init(sys%m, h%no_lasers, h%lasers)
  if(h%no_lasers>0 ) then
      message(1) = 'Info: Lasers'
      call write_info(1)
      if(conf%verbose > 20 .and. mpiv%node == 0) then
        call laser_write_info(h%no_lasers, h%lasers, stdout)
      end if
      call oct_parse_logical(C_string("TDOutputLaser"), .false., h%output_laser)  
  end if

  ! absorbing boundaries
  call oct_parse_int(C_string("TDAbsorbingBoundaries"), 0, dummy)
  nullify(h%ab_pot)

  absorbing_boundaries: if(dummy .eq. 1 .or. dummy .eq. 2) then

  h%ab = dummy
  call oct_parse_double(C_string("TDABWidth"), 4._r8/units_inp%length%factor, h%ab_width)
  h%ab_width  = h%ab_width * units_inp%length%factor
  if(h%ab == 1) then
    call oct_parse_double(C_string("TDABHeight"), -0.2_r8/units_inp%energy%factor, h%ab_height)
    h%ab_height = h%ab_height * units_inp%energy%factor
  else
    h%ab_height = 1._r8
  end if

 ! generate boundary potential...
  allocate(h%ab_pot(sys%m%np))
  h%ab_pot = 0._r8
  pot: do i = 1, sys%m%np
     call mesh_r(sys%m, i, r, x=x)

     select case(sys%m%box_shape)
     case(SPHERE)
          d = r - (sys%m%rsize - h%ab_width)
          if(d.gt.0._r8) then
            h%ab_pot(i) = h%ab_height * sin(d*M_PI/(2._r8*h%ab_width))**2
          end if

#if defined(THREE_D)
     case(CYLINDER)
          d = sqrt(x(1)**2 + x(2)**2) - (sys%m%rsize - h%ab_width)
          if(d.gt.0._r8)  &
               h%ab_pot(i) = h%ab_height * sin(d*M_PI/(2._r8*h%ab_width))**2
          d = abs(x(3)) - (sys%m%zsize - h%ab_width)
          if(d.gt.0._r8)  &
               h%ab_pot(i) = h%ab_pot(i) + h%ab_height * sin(d*M_PI/(2._r8*h%ab_width))**2

     case(PARALLELEPIPED)
          do j = 1, 3
            d = x(j) - (sys%m%lsize(j)/2._r8 - h%ab_width)
            if(d.gt.0._r8) then
              h%ab_pot(i) = h%ab_pot(i) + h%ab_height * sin(d*M_PI/(2._r8*h%ab_width))**2
            end if
          end do
#endif

     case default
          message(1) = "Absorbing boundaries are not implemented for"
          message(2) = "Box_shape = 3"
          call write_warning(2)
          exit pot
     end select

     if(abs(h%ab_pot(i)) > abs(h%ab_height)) h%ab_pot(i) = h%ab_height
  end do pot

  end if absorbing_boundaries

  call pop_sub(); return
end subroutine hamiltonian_init

subroutine hamiltonian_end(h)
  type(hamiltonian_type) :: h

  sub_name = 'hamiltonian_end'; call push_sub()

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhartree, h%Vxc)
    nullify(h%Vpsl, h%Vhartree, h%Vxc)
  end if

  if(h%classic_pot > 0 .and. associated(h%Vclassic)) then
    deallocate(h%Vclassic); nullify(h%Vclassic)
  end if

  if(.not.h%ip_app) then
    call hartree_end(h%hart)
    call xc_end(h%xc)
  end if

  if(associated(h%ab_pot)) then
    deallocate(h%ab_pot); nullify(h%ab_pot)
  end if

  call laser_end(h%no_lasers, h%lasers)

  call pop_sub(); return
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

  call pop_sub(); return
end subroutine hamiltonian_energy

subroutine hamiltonian_span(h, delta, emin)
  type(hamiltonian_type), intent(inout) :: h
  real(r8), intent(in) :: delta, emin
  sub_name = 'hamiltonian_span'; call push_sub()
    h%spectral_middle_point = ((M_Pi**2/(2*delta**2))+emin)/2._r8
    h%spectral_half_span    = ((M_Pi**2/(2*delta**2))-emin)/2._r8
  call pop_sub(); return
end subroutine hamiltonian_span

subroutine hamiltonian_output(h, m, dir, outp)
  type(hamiltonian_type), intent(IN) :: h
  type(mesh_type), intent(IN) :: m
  character(len=*), intent(IN) :: dir
  type(output_type), intent(IN) :: outp

  integer :: is
  character(len=80) :: fname  
  real(r8) :: u
  sub_name = 'hamiltonian_output'; call push_sub()

  u = units_out%energy%factor
  if(outp%what(output_potential)) then
    call doutput_function(outp, dir, "v0", m, h%Vpsl, u)

    if(h%classic_pot > 0) then
      call doutput_function(outp, dir, "vc", m, h%Vclassic, u)
    end if

    if(.not.h%ip_app) then
      call doutput_function(outp, dir, "vh", m, h%Vhartree, u)

      do is = 1, min(h%ispin, 2)
        write(fname, '(a,i1)') 'vxc-', is
        call doutput_function(outp, dir, fname, m, h%Vxc(:, is), u)
      end do
    end if
  end if

  call pop_sub(); return
end subroutine hamiltonian_output

#include "h_external_pot.F90"

#include "undef.F90"
#include "real.F90"
#include "h_inc.F90"
#include "h_forces.F90"

#include "undef.F90"
#include "complex.F90"
#include "h_inc.F90"
#include "h_forces.F90"

#if defined(COMPLEX_WFNS)
#include "h_so.F90"
#endif
!#include "h_rel.F90"

end module hamiltonian
