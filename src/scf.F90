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

module scf
use io
use units
use output
use system
use states
use hamiltonian
use eigen_solver
use mix
use lcao

implicit none

type scf_type  ! some variables used for the scf cycle
  integer :: max_iter ! maximum number of scf iterations
  real(r8) :: conv_abs_dens, conv_rel_dens, &
       conv_abs_ener, conv_rel_ener ! several convergence criteria
  
  real(r8) :: abs_dens, rel_dens, abs_ener, rel_ener

  logical :: lcao_restricted

  type(mix_type) :: smix
  type(eigen_solver_type) :: eigens
end type scf_type
  
contains

subroutine scf_init(scf, sys)
  type(scf_type), intent(inout) :: scf
  type(system_type), intent(IN) :: sys

  sub_name = 'systm_scf_init'; call push_sub()

  call oct_parse_int(C_string("MaximumIter"), 200, scf%max_iter)
  call oct_parse_double(C_string("ConvAbsDens"), 1e-5_r8, scf%conv_abs_dens)
  call oct_parse_double(C_string("ConvRelDens"),   0._r8, scf%conv_rel_dens)
  call oct_parse_double(C_string("ConvAbsEnergy"), 0._r8, scf%conv_abs_ener)
  call oct_parse_double(C_string("ConvRelEnergy"), 0._r8, scf%conv_rel_ener)

  if(scf%max_iter <= 0 .and. &
      scf%conv_abs_dens <= 0.0_r8 .and. scf%conv_rel_dens <= 0.0_r8 .and. &
      scf%conv_abs_ener <= 0.0_r8 .and. scf%conv_rel_ener <= 0.0_r8) then
    message(1) = "Input: Not all convergence criteria can be <= 0"
    message(2) = "Please set one of the following:"
    message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEnergy | ConvRelEnergy"
    call write_fatal(3)
  end if

  if(scf%max_iter <= 0) scf%max_iter = huge(scf%max_iter)

  ! Handle mixing now...
  call mix_init(scf%smix, sys%m, sys%st)

  ! now the eigen solver stuff
  call eigen_solver_init(scf%eigens)

  ! Should the calculation be restricted to LCAO subspace?
  call oct_parse_logical(C_string("SCFinLCAO"), .false., scf%lcao_restricted)
  if(scf%lcao_restricted) then
    message(1) = 'Info: SCF restricted to LCAO subspace'
    call write_info(1)
  endif

  call pop_sub()
  return
end subroutine scf_init

subroutine scf_end(scf)
  type(scf_type), intent(inout) :: scf

  call mix_end(scf%smix)

  return
end subroutine scf_end

subroutine scf_run(scf, sys, h)
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  type(scf_type), intent(inout) :: scf

  integer :: iter, iunit, ik, ist, id
  real(r8) :: old_etot
  real(r8), allocatable :: diff(:, :)
  logical :: finish

  sub_name = 'scf_run'; call push_sub()

  if(scf%lcao_restricted) call lcao_init(sys, h)

  allocate(diff(sys%st%nst, sys%st%nik))

  do iter = 1, scf%max_iter
    if(scf%lcao_restricted) then
      call lcao_wf(sys, h)
    else
      call eigen_solver_run(scf%eigens, sys, h, iter, diff)
    endif

    ! compute new density
    call mix_dens(scf%smix, iter, sys%st, sys%m, scf%abs_dens)
    scf%rel_dens = scf%abs_dens / sys%st%qtot

    ! compute new potentials
    call R_FUNC(hamiltonian_setup) (h, sys%m, sys%st, sys)

    ! occupations
    call states_fermi(sys%st, sys%m)

    ! output eigenvalues
    call states_write_eigenvalues(stdout, sys%st%nst, sys%st, diff)

    ! now compute total energy
    old_etot = h%etot
    call hamiltonian_energy(h, sys%st, sys%eii, -1)
    scf%abs_ener = abs(old_etot - h%etot)
    scf%rel_ener = scf%abs_ener / abs(h%etot)

    ! are we finished?
    finish = &
        (scf%conv_abs_dens > 0.0_r8 .and. scf%abs_dens <= scf%conv_abs_dens) .or. &
        (scf%conv_rel_dens > 0.0_r8 .and. scf%rel_dens <= scf%conv_rel_dens) .or. &
        (scf%conv_abs_ener > 0.0_r8 .and. scf%abs_ener <= scf%conv_abs_ener) .or. &
        (scf%conv_rel_ener > 0.0_r8 .and. scf%rel_ener <= scf%conv_rel_ener)

    write(message(1), '(a,i4,a,e14.8,a,e14.8)') &
         'Info: iter = ', iter, ' abs_dens = ', scf%abs_dens, &
         ' abs_ener = ', scf%abs_ener
    write(message(2), '(a)') ''
    call write_info(2)

    ! save restart information
    if(finish.or.(modulo(iter, 3) == 0).or.clean_stop()) &
         call R_FUNC(states_write_restart)("restart.static", sys%m, sys%st)

    if(finish) then
      write(message(1), '(a, i4, a)')'Info: SCF converged in ', iter, ' iterations'
      call write_info(1)
      if(scf%lcao_restricted) call lcao_end
      exit
    end if

    if(clean_stop()) exit
  end do

  if(.not.finish) then
    message(1) = 'SCF *not* converged!'
    call write_warning(1)
  end if

  ! calculate forces
  call R_FUNC(forces)(h, sys)

  ! output final information
  call scf_write_static("static", "info")
  call R_FUNC(states_output) (sys%st, sys%m, "static", sys%outp)
  if(sys%outp%what(output_geometry)) call system_write_xyz("static", "geometry", sys)
  call hamiltonian_output(h, sys%m, "static", sys%outp)

  deallocate(diff)
  call pop_sub()

contains

subroutine scf_write_static(dir, fname)
  character(len=*), intent(in) :: dir, fname

  integer iunit, i, j

  call oct_mkdir(C_string(trim(dir)))
  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(dir)+"/"+trim(fname))

  ! mesh
  write(iunit, '(a,a)') 'System name: ', sys%sysname
  write(iunit, '(1x)')

  write(iunit, '(a)') 'Mesh:'
  call mesh_write_info(sys%m, iunit)
  write(iunit,'(1x)')

  if(.not. h%ip_app) then
    write(iunit, '(a)') 'Exchange and correlation functionals:'
    call xc_write_info(h%xc, iunit)
  else
    write(iunit, '(a)') 'Independent Particles'
  end if
  write(iunit,'(1x)')

  ! scf information
  if(finish) then
    write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
  else
    write(iunit, '(a)') 'SCF *not* converged!'
  end if
  write(iunit, '(1x)')

  call states_write_eigenvalues(iunit, sys%st%nst, sys%st, diff)
  write(iunit, '(1x)')

  write(iunit, '(a)') 'Energy:'
  call hamiltonian_energy(h, sys%st, sys%eii, iunit)
  write(iunit, '(1x)')

  if(sys%st%ispin > 1) then
    call write_magnet(iunit, sys%st)
  end if

  do i = 1, sys%st%nspin
    write(iunit, '(a, i1, 4a)') 'Dipole (', i, ')', &
         ' [', trim(units_out%length%abbrev), ']:'
    do j = 1, conf%dim
      write(iunit, '(6x,a,i1,a,es17.8)') '<x', j, '> = ', &
           sum(sys%m%Lxyz(j,:)*sys%st%rho(:,i))*sys%m%vol_pp*sys%m%h(j) / units_out%length%factor
    end do
    write(iunit, '(1x)')
  end do

  write(iunit, '(a)') 'Convergence:'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'abs_dens = ', scf%abs_dens, &
      ' (', scf%conv_abs_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_dens = ', scf%rel_dens, &
      ' (', scf%conv_rel_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'abs_ener = ', scf%abs_ener, &
      ' (', scf%conv_abs_ener / units_out%energy%factor, ')', &
      ' [',  trim(units_out%energy%abbrev), ']'
  write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'rel_ener = ', scf%rel_ener, &
      ' (', scf%conv_rel_ener / units_out%energy%factor, ')', &
      ' [',  trim(units_out%energy%abbrev), ']'
  write(iunit,'(1x)') 

  write(iunit,'(3a)') 'Forces on the ions [', trim(units_out%force%abbrev), "]"
  write(iunit,'(a,10x,14x,a,14x,a,14x,a)') ' Ion','x','y','z'
  do i = 1,sys%natoms
    write(iunit,'(i4,a10,3f15.6)') i, trim(sys%atom(i)%spec%label), &
         sys%atom(i)%f(:) / units_out%force%factor
  end do

  call io_close(iunit)
end subroutine scf_write_static

subroutine write_magnet(iunit, st)
  integer, intent(in) :: iunit
  type(states_type), intent(IN) :: st
  
  real(r8) :: m(3), sign
  R_TYPE :: c
  integer :: i, ik, ist
  
  write(iunit, '(a)') 'Magnetization:'
  if(st%ispin == 2) then ! collinear spin
    sign = 1._r8
    m(3) = 0._r8
    do ik = 1, st%nik
      do ist = 1, st%nst
        m(3) = m(3) + sign*st%kweights(ik)*st%occ(ist, ik)
      end do
      sign = -sign
    end do
    write(iunit, '(a,f15.6)') ' mz = ', m(3)
    
  else if(st%ispin == 3) then ! non-collinear
    m = 0._r8
    do ik = 1, st%nik
      do ist = 1, st%nst
        do i = 1, sys%m%np
          c = R_CONJ(st%R_FUNC(psi) (i, 1, ist, ik)) * st%R_FUNC(psi) (i, 2, ist, ik)
          m(1) = m(1) + st%kweights(ik)*st%occ(ist, ik)* 2._r8*R_REAL(c)
          m(2) = m(2) + st%kweights(ik)*st%occ(ist, ik)* 2._r8*R_AIMAG(c)
          c = R_ABS(st%R_FUNC(psi) (i, 1, ist, ik))**2 - R_ABS(st%R_FUNC(psi) (i, 2, ist, ik))**2
          m(3) = m(3) + st%kweights(ik)*st%occ(ist, ik)* 2._r8*R_REAL(c)
        end do
      end do
    end do
    m = m*sys%m%vol_pp
    write(iunit, '(a,f15.6)') ' mx = ', m(1)
    write(iunit, '(a,f15.6)') ' my = ', m(2)
    write(iunit, '(a,f15.6)') ' mz = ', m(3)
  end if
  
  write(iunit,'(1x)') 
  
end subroutine write_magnet

end subroutine scf_run

end module scf
