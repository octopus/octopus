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
use math
use oct_parser
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

  type(fft_type) :: fft  ! for the split operator method

  !variables controlling the output
  logical :: out_multip  ! multipoles
  logical :: out_coords  ! coordinates
  logical :: out_gsp     ! projection onto the ground state.
  logical :: out_acc     ! electronic acceleration
  logical :: out_laser   ! laser field
  logical :: out_energy  ! several components of the electronic energy
  logical :: out_proj    ! projection onto the GS KS eigenfunctions

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
  integer(POINTER_SIZE) :: out_multip, out_coords, out_gsp, out_acc, &
       out_laser, out_energy, out_proj

  real(r8), allocatable ::  x1(:,:), x2(:,:), f1(:,:) ! stuff for verlet
  real(r8) :: etime
  character(len=100) :: filename

  call push_sub('td_run')
  
  if(mpiv%node==0) then
    write(filename, '(i3.3)') mpiv%node
    if(td%out_multip) &
         call write_iter_init(out_multip, td%iter, td%dt/units_out%time%factor, "td.general/multipoles")
    if(td%out_coords) &
         call write_iter_init(out_coords, td%iter, td%dt/units_out%time%factor, "td.general/coordinates")
    if(td%out_gsp) &
         call write_iter_init(out_gsp,    td%iter, td%dt/units_out%time%factor, "td.general/gs_projection")
    if(td%out_acc) &
         call write_iter_init(out_acc,    td%iter, td%dt/units_out%time%factor, "td.general/acceleration")
    if(td%out_laser) &
         call write_iter_init(out_laser,  td%iter, td%dt/units_out%time%factor, "td.general/laser")
    if(td%out_energy) &
         call write_iter_init(out_energy, td%iter, td%dt/units_out%time%factor, "td.general/el_energy")
    if(td%out_proj) &
         call write_iter_init(out_proj,   td%iter, td%dt/units_out%time%factor, "td.general/projections."//trim(filename))
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

    ! update density
    call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call zhamiltonian_setup(h, sys%m, sys%st, sys)
    call hamiltonian_energy(h, sys%st, sys%eii, -1, reduce=.true.)

    ! Recalculate forces, update velocities...
    if(td%move_ions > 0) then
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, sys%natoms
          f1(j, :) = sys%atom(j)%f(:)
        end do
      end if
      call zforces(h, sys, i*td%dt, reduce=.true.)
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, sys%natoms
          if(sys%atom(j)%move) then
            sys%atom(j)%v(:) = sys%atom(j)%v(:) + &
                 td%dt/(2._r8*sys%atom(j)%spec%weight) * &
                 (f1(j, :) + sys%atom(j)%f(:))             
          end if
        end do
      end if
      sys%kinetic_energy = kinetic_energy(sys%natoms, sys%atom)
    end if

    ! output multipoles
    if(td%out_multip) call td_write_multipole(i)

    ! output positions, vels, etc.
    if(td%out_coords) call td_write_nbo(i, sys%kinetic_energy, h%etot)

    ! If harmonic spectrum is desired, get the acceleration
    if(td%out_acc) call td_write_acc(i)

    ! output laser field
    if(td%out_laser) call td_write_laser(i)

    ! output electronic energy
    if(td%out_energy) call td_write_el_energy(i)

    ! output projections onto the GS KS eigenfunctions
    if(td%out_proj) call td_write_proj(i)
    
#ifndef DISABLE_PES
    call PES_doit(td%PESv, sys%m, sys%st, ii, td%dt, h%ab_pot)
#endif

    ! write info
    write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         (oct_clock() - etime)/1e6
    call write_info(1)
    etime = oct_clock()

    ! write down data
    ii = ii + 1
    if(ii==sys%outp%iter+1 .or. i == td%max_iter) then ! output
      if(i == td%max_iter) sys%outp%iter = ii - 1
      ii = 1

      call td_write_data(i)
    end if

  end do

  ! close all buffers
  if(mpiv%node==0) then
    if(td%out_multip) call write_iter_end(out_multip)
    if(td%out_coords) call write_iter_end(out_coords)
    if(td%out_gsp)    call write_iter_end(out_gsp)
    if(td%out_acc)    call write_iter_end(out_acc)
    if(td%out_laser)  call write_iter_end(out_laser)
    if(td%out_energy) call write_iter_end(out_energy)
    if(td%out_proj)   call write_iter_end(out_proj)
  end if

contains

  subroutine td_run_zero_iter(m)
    type(mesh_type), intent(IN) :: m

    integer :: i, iunit
    real(r8) :: x(conf%dim)
    complex(r8) :: c

    call push_sub('td_run_zero_iter')

    ! create general subdir
    call oct_mkdir("td.general")

    if(td%out_multip) call td_write_multipole(0)
    if(td%out_proj)   call td_write_proj(0)

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

    ! create files for output and output headers
    if(td%out_coords) call td_write_nbo(0, sys%kinetic_energy, h%etot)    
    if(td%out_acc)    call td_write_acc(0)
    if(td%out_laser)  call td_write_laser(0)
    if(td%out_energy) call td_write_el_energy(0)
    call td_write_data(0)

    call pop_sub()
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
!!$  allocate(expzpsi (sys%m%np, sys%st%dim), &
!!$           goodzpsi(sys%m%np, sys%st%dim), &
!!$           hzpsi   (sys%m%np, sys%st%dim))
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
