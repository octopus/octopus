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

module timedep
use global
use io
use system
use states
use hamiltonian
use mix
use lasers
#ifndef DISABLE_PES
use PES
#endif

implicit none

type td_type
  complex(r8), pointer :: zpsi(:,:,:) ! the complex wavefunctions

  integer :: max_iter  ! maximum number of iterations to perform
  integer :: iter      ! the actual iteration

  integer :: evolution_method ! which evolution method to use
  integer :: exp_method       ! which method is used to apply the exponential
  real(r8) :: lanczos_tol
  integer  :: exp_order

  real(r8) :: dt            ! time step
  integer  :: move_ions     ! how do we move the ions?

  real(r8) :: delta_strength  ! strength of the delta excitation
  real(r8), pointer :: pol(:) ! direction of the polarization of the efield

  integer :: lmax        ! maximum multipole moment to write

  logical :: occ_analysis ! do we perform occupational analysis?

  logical :: harmonic_spectrum   ! write harmonic spectrum

#ifndef DISABLE_PES
  type(PES_type) :: PESv
#endif

  real(r8), pointer :: v_old(:, :, :)
end type td_type

  ! Parameters.
  integer, parameter :: STATIC_IONS     = 0,    &
                        NORMAL_VERLET   = 3,  &
                        VELOCITY_VERLET = 4

  integer, parameter :: SIMPLE_EXP           = 0, &
                        OLD_REVERSAL         = 1, &
                        REVERSAL             = 2, &
                        APP_REVERSAL         = 3, &
                        EXPONENTIAL_MIDPOINT = 4

  integer, parameter :: FOURTH_ORDER       = 1, &
                        LANCZOS_EXPANSION  = 2, &
                        SPLIT_OPERATOR     = 3, &
                        SUZUKI_TROTTER     = 4, &
                        CHEBYSHEV          = 5

contains

subroutine td_run(td, u_st, sys, h)
  type(td_type), intent(inout) :: td
  type(states_type), intent(IN) :: u_st
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h

  integer :: i, ii, j, idim, ist, ik
  real(r8), allocatable :: dipole(:,:), multipole(:,:,:), x(:,:,:), v(:,:,:), f(:,:,:), &
                           x1(:,:), x2(:,:), f1(:,:), ke(:), pe(:), tacc(:, :)
  real(r8) :: etime
  complex(r8), allocatable :: projections(:,:,:,:)
  character(len=100) :: filename

  sub_name = 'td_run'; call push_sub()
  
  allocate(dipole(sys%st%nspin, sys%outp%iter))
  allocate(multipole((td%lmax + 1)**2, sys%st%nspin, sys%outp%iter))
  if(td%move_ions > 0) then
     allocate(x(sys%outp%iter, sys%natoms, conf%dim), v(sys%outp%iter, sys%natoms, conf%dim), &
              f(sys%outp%iter, sys%natoms, conf%dim))
     allocate(ke(sys%outp%iter), pe(sys%outp%iter))
  endif
  if(td%harmonic_spectrum) allocate(tacc(sys%outp%iter, conf%dim))

  ! occupational analysis stuff
  if(td%occ_analysis) then
    allocate(projections(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik, sys%outp%iter))
  end if

  ! Calculate initial forces and kinetic energy
  if(td%move_ions > 0) then 
    if(td%iter > 0) then
      call td_read_nbo()
      call generate_external_pot(h, sys)
      sys%eii = ion_ion_energy(sys%natoms, sys%atom)
    end if

    call zforces(h, sys, td%iter*td%dt, reduce=.true.)
    sys%kinetic_energy = kinetic_energy(sys%natoms, sys%atom)
    select case(td%move_ions)
      case(NORMAL_VERLET)
        allocate(x1(sys%natoms, conf%dim), x2(sys%natoms, conf%dim))
        do j = 1, sys%natoms
           if(sys%atom(j)%move) then
             x1(j, :) = sys%atom(j)%x(:) - td%dt*sys%atom(j)%v(:) + &
                        0.5_r8 * td%dt**2/sys%atom(j)%spec%weight * &
                        sys%atom(j)%f(:)
           else
             x1(j, :) = sys%atom(j)%x(:)
           endif
        enddo
      case(VELOCITY_VERLET)
        allocate(f1(sys%natoms, conf%dim))
    end select
  endif

  if(td%iter == 0) call td_run_zero_iter(sys%m)
  !call td_check_trotter(td, sys, h)
  td%iter = td%iter + 1

  write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
  call write_info(1)

  ii = 1
  etime = oct_clock()
  do i = td%iter, td%max_iter
    if(clean_stop()) exit
    ! Move the ions.
    if( td%move_ions > 0 ) then
      select case(td%move_ions)
        case(NORMAL_VERLET)
          x2 = x1
          do j = 1, sys%natoms
             if(sys%atom(j)%move) then
                x1(j, :) = sys%atom(j)%x(:)
                sys%atom(j)%x(:) = 2._r8*x1(j, :) - x2(j, :) + &
                   td%dt**2/sys%atom(j)%spec%weight * sys%atom(j)%f(:)
                sys%atom(j)%v(:) = (sys%atom(j)%x(:) - x2(j, :)) / (2._r8*td%dt)
             endif
          enddo
        case(VELOCITY_VERLET)
          do j=1, sys%natoms
             if(sys%atom(j)%move) then
                sys%atom(j)%x(:) = sys%atom(j)%x(:) +  td%dt*sys%atom(j)%v(:) + &
                   0.5_r8*td%dt**2/sys%atom(j)%spec%weight * sys%atom(j)%f(:)
             endif
          enddo
      end select

      call generate_external_pot(h, sys)
      sys%eii = ion_ion_energy(sys%natoms, sys%atom)
    endif

    ! time iterate wavefunctions
    call td_rti(h, sys%m, sys%st, sys, td, i*td%dt)

    ! update density
    call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call zhamiltonian_setup(h, sys%m, sys%st, sys)
    call hamiltonian_energy(h, sys%st, sys%eii, -1, reduce=.true.)

    ! Recalculate forces, update velocities...
    if(td%move_ions > 0) then
      pe(ii) = h%etot
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, sys%natoms
           f1(j, :) = sys%atom(j)%f(:)
        enddo
      endif
      call zforces(h, sys, i*td%dt, reduce=.true.)
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, sys%natoms
           if(sys%atom(j)%move) then
              sys%atom(j)%v(:) = sys%atom(j)%v(:) + &
                                 td%dt/(2._r8*sys%atom(j)%spec%weight) * &
                                 (f1(j, :) + sys%atom(j)%f(:))             
           endif
        enddo
      endif
      sys%kinetic_energy = kinetic_energy(sys%natoms, sys%atom)
      ke(ii) = sys%kinetic_energy

      ! store data
      do j=1, sys%natoms
        x(ii, j, 1:conf%dim) = sys%atom(j)%x(1:conf%dim)
        v(ii, j, 1:conf%dim) = sys%atom(j)%v(1:conf%dim)
        f(ii, j, 1:conf%dim) = sys%atom(j)%f(1:conf%dim)
      end do
    end if

    ! If harmonic spectrum is desired, get the acceleration
    if(td%harmonic_spectrum) call td_calc_tacc(tacc(ii, 1:conf%dim), td%dt*i, reduce = .true.)

    ! measuring
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, &
         dipole(:,ii), multipole(:,:,ii))
    if(td%occ_analysis) then
      call calc_projection(u_st, sys%st, sys%m, projections(:, :, :, ii))
    end if

#ifndef DISABLE_PES
    call PES_doit(td%PESv, sys%m, sys%st, ii, td%dt, h%ab_pot)
#endif

    ! mask function?
    if(h%ab == 2) then
      do ik = 1, sys%st%nik
        do ist = sys%st%st_start, sys%st%st_end
          do idim = 1, sys%st%dim
            sys%st%zpsi(1:sys%m%np, idim, ist, ik) = sys%st%zpsi(1:sys%m%np, idim, ist, ik) * &
                 (1._r8 - h%ab_pot(1:sys%m%np))
          end do
        end do
      end do
    end if

    ! write info
    write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         (oct_clock() - etime)/1e6
    call write_info(1)
    etime = oct_clock()

    ! write down data
    ii = ii + 1
    save: if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output
      if(i == td%max_iter) sys%outp%iter = ii - 1
      ii = 1

      ! first resume file
      write(filename, '(a,i3.3)') "tmp/restart.td.", mpiv%node
      call zstates_write_restart(trim(filename), sys%m, sys%st, &
           iter=i, v1=td%v_old(:, :, 2), v2=td%v_old(:, :, 3))

      call td_write_data()

      ! now write down the rest
      write(filename, '(a,i7.7)') "td.", i  ! name of directory
      call zstates_output(sys%st, sys%m, filename, sys%outp)
      if(sys%outp%what(output_geometry)) call system_write_xyz(filename, "geometry", sys)
      call hamiltonian_output(h, sys%m, filename, sys%outp)
      
    end if save

  end do

  deallocate(dipole, multipole)
  if(td%occ_analysis) then
    deallocate(projections)
  end if
contains

  subroutine td_run_zero_iter(m)
    type(mesh_type), intent(IN) :: m

    integer :: i, iunit
    real(r8) :: x(conf%dim), t_acc(conf%dim)
    complex(r8) :: c
    real(r8), allocatable :: dipole(:), multipole(:,:), pos(:,:), vel(:,:), for(:,:)

    sub_name = 'td_run_zero_iter'; call push_sub()

    ! create general subdir
    call oct_mkdir(C_string("td.general"))

    ! output static dipole (iter 0)
    allocate(dipole(sys%st%nspin), multipole((td%lmax + 1)**2, sys%st%nspin))
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, &
         dipole, multipole)
    call io_assign(iunit)
    open(iunit, status='unknown', file='td.general/multipoles')
    call td_write_multipole(iunit, 0_i4, 0._r8, dipole, multipole, .true.)
    call io_close(iunit)
    deallocate(dipole, multipole)

    ! we now apply the delta(0) impulse to the wf
    if(td%delta_strength .ne. M_ZERO) then
      write(message(1),'(a,f11.6)')  'Info: Applying delta kick: k = ', td%delta_strength
      call write_info(1)
      do i = 1, m%np
        call mesh_xyz(m, i, x)
        c = exp(M_zI * td%delta_strength * sum(x(1:conf%dim)*td%pol(1:conf%dim)))
        sys%st%zpsi(i,:,:,:) = c * sys%st%zpsi(i,:,:,:)
      end do
    end if

    select case(sys%st%ispin)
    case(UNPOLARIZED, SPIN_POLARIZED)
      do i = 1, sys%st%nspin
         td%v_old(:, i, 2) = h%vhartree(:) + h%vxc(:, i)
      enddo
    case(SPINORS)
      td%v_old(:, 1, 2) = h%vhartree(:) + h%vxc(:, 1)
      td%v_old(:, 2, 2) = h%vhartree(:) + h%vxc(:, 2)
      td%v_old(:, 3, 2) = h%vxc(:, 3)
      td%v_old(:, 4, 2) = h%vxc(:, 4)
    end select
    td%v_old(:, :, 3) = td%v_old(:, :, 2)

    ! And same thing for the ions. The formulas are:
    ! If phi(0+) = e^(ikz)phi_0 (kick for the electrons, as done before,
    ! we have applied an electrical field in the form E_0 delta(t),
    ! being E_0 = - k \hbar / e,
    ! and thus the kick for the nuclei is:
    ! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
    ! where M and Z are the ionic mass and charge, respectively.
    if(td%move_ions > 0) then
      do j = 1, sys%natoms
          sys%atom(j)%v(1:conf%dim) = sys%atom(j)%v(1:conf%dim) - &
            td%pol(1:conf%dim)*sys%atom(j)%spec%z_val*td%delta_strength / sys%atom(j)%spec%weight
      enddo
    endif

    ! output laser
    if(h%output_laser) then
      call io_assign(iunit)
      open(unit=iunit, file='td.general/laser', status='unknown')
      call td_write_laser(iunit, 0, 0._r8, .true.)
      call io_close(iunit)
    end if

    ! output occupational analysis data
    if(td%occ_analysis) then
      call io_assign(iunit)
      write(filename, '(a,i3.3)') 'td.general/projections.', mpiv%node
      open(iunit, status='unknown', file=filename)
      write(iunit,'(a,4i5)') '#',sys%st%nik, sys%st%st_start, sys%st%st_end, u_st%nst
      call io_close(iunit)
    end if

    ! output positions, velocities, forces...
    if(td%move_ions > 0) then
       call io_assign(iunit)
       open(iunit, file='td.general/coordinates')
       allocate(pos(sys%natoms, conf%dim), vel(sys%natoms, conf%dim), for(sys%natoms, conf%dim))
       do j=1, sys%natoms
          pos(j, 1:conf%dim) = sys%atom(j)%x(1:conf%dim)
          vel(j, 1:conf%dim) = sys%atom(j)%v(1:conf%dim)
          for(j, 1:conf%dim) = sys%atom(j)%f(1:conf%dim)
       enddo
       call td_write_nbo(iunit, 0, 0._r8, sys%kinetic_energy, h%etot, pos, vel, for, header = .true.)
       call io_close(iunit)
       deallocate(pos, vel, for)
    endif

    ! output harmonic spectrum
    if(td%harmonic_spectrum) then
       call io_assign(iunit)
       open(iunit, file='td.general/acceleration')
       call td_calc_tacc(t_acc, 0.0_r8, reduce = .true.)
       call td_write_acc(iunit, 0, 0.0_r8, t_acc, header=.true.)
       call io_close(iunit)
    endif

    call pop_sub(); return    
  end subroutine td_run_zero_iter

  subroutine td_read_nbo() ! reads the pos and vel from coordinates file
    logical :: found
    integer :: i, iunit

    inquire(file='td.general/coordinates', exist=found)
    if(.not.found) then
      message(1) = "Could not open file 'td.general/coordinates'"
      message(2) = "Starting simulation from initial geometry"
      call write_warning(2)
      return
    end if
    
    call io_assign(iunit)
    open(iunit, file='td.general/coordinates', status='old')

    read(iunit, *); read(iunit, *) ! skip header
    do i = 0, td%iter - 1
      read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
    end do
    read(iunit, '(88x)', advance='no') ! skip unrelevant information

    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%x(1:conf%dim)
      sys%atom(i)%x(:) = sys%atom(i)%x(:) * units_out%length%factor
    end do
    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%v(1:conf%dim)
      sys%atom(i)%v(:) = sys%atom(i)%v(:) * units_out%velocity%factor
    end do
    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%f(1:conf%dim)
      sys%atom(i)%f(:) = sys%atom(i)%f(:) * units_out%force%factor
    end do

    call io_close(iunit)
  end subroutine td_read_nbo

#include "td_calc.F90"
#include "td_write.F90"

end subroutine td_run

!!$subroutine td_check_trotter(td, sys, h)
!!$  type(td_type), intent(inout)          :: td
!!$  type(system_type), intent(in)         :: sys
!!$  type(hamiltonian_type), intent(inout) :: h
!!$
!!$  integer :: unit, order, i, j
!!$  complex(r8), allocatable :: expzpsi(:, :), goodzpsi(:, :), hzpsi(:, :)
!!$  real(r8) :: res(8), dt
!!$  complex(r8) :: r
!!$
!!$  call io_assign(unit)
!!$  open(unit, file = 'td_check')
!!$  close(unit)
!!$
!!$  allocate(expzpsi (0:sys%m%np, sys%st%dim), &
!!$           goodzpsi(0:sys%m%np, sys%st%dim), &
!!$           hzpsi   (  sys%m%np, sys%st%dim))
!!$
!!$ do i = -300, 100, 10
!!$
!!$     ! Calculate the "good" valee
!!$
!!$     dt = exp(real(i, r8)/100._r8 * log(10.0_r8))
!!$     write(*, *) 'Calculating the "good" value...'
!!$     goodzpsi = sys%st%zpsi(:, :, 1, 1)
!!$     td%lanczos_tol = 1.0e-14_r8
!!$     td%exp_order   = 64
!!$     td%exp_method  = LANCZOS_EXPANSION
!!$     call td_dtexp(h, sys, td, 1, goodzpsi, dt, 0.0_r8)
!!$     !goodzpsi(:, :) = exp(-M_zI*dt*sys%st%eigenval(1, 1))*sys%st%zpsi(:, :, 1, 1)
!!$     write(*, *) 'Done.'
!!$
!!$     !do j = 1, 8
!!$       expzpsi = sys%st%zpsi(:, :, 1, 1)
!!$       td%exp_order   = 20
!!$       td%lanczos_tol = 1.0e-4_r8
!!$       td%exp_method  = SUZUKI_TROTTER
!!$       call td_dtexp(h, sys, td, 1, expzpsi,  dt, 0.0_r8)
!!$            expzpsi = expzpsi - goodzpsi
!!$       res(1) = zstates_nrm2(sys%m, sys%st%dim, expzpsi(1:sys%m%np, 1:sys%st%dim))
!!$     !enddo
!!$       open(unit, file='td_check', position='append')
!!$       write(unit, '(2es18.6, i5)') dt, res(1), order
!!$       close(unit)
!!$  end do
!!$
!!$  deallocate(expzpsi, goodzpsi)
!!$  call io_close(unit)
!!$  stop
!!$end subroutine td_check_trotter

#include "td_init.F90"
#include "td_rti.F90"
#include "td_exp.F90"

end module timedep
