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

module opt_control
use global
use io
use system
use states
use hamiltonian
use mix
use lasers
use timedep

implicit none

contains

subroutine opt_control_run(td, sys, h)
  type(td_type), intent(inout) :: td
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h

  type(states_type), pointer :: psi_i
  type(states_type)          :: psi_f
  real(r8), pointer :: v_old_i(:,:,:), v_old_f(:,:,:)
  real(r8), pointer :: laser_i(:,:), laser_f(:,:)

  integer :: i, iunit
  real(r8) :: alpha, functional

  call init()

  ! first propagate psi_f to ti
  call prop_psi_f()

  ! we now setup v_old arrays
  do i = 1, 3
    v_old_f(:,:,i) = td%v_old(:,:,1) ! this one comes from prop_psi_f
  end do

  call zcalcdens(psi_i, sys%m%np, psi_i%rho, reduce=.true.)
  call zhamiltonian_setup(h, sys%m, psi_i, sys)
  do i = 1, sys%st%nspin
    v_old_i(:, i, 2) = h%Vhartree(:) + h%Vxc(:, i)
  end do
  v_old_i(:, :, 3) = v_old_i(:, :, 2)

  do
    ! first propagate both states forward
    message(1) = "Info: Propagating forward"
    call write_info(1)

    call read_state(psi_i, "wf.initial")
    call oct_progress_bar(-1, td%max_iter-1)
    h%lasers(1)%dt = td%dt
    functional = M_ZERO
    do i = 1, td%max_iter
      call prop_iter1(i)
      call oct_progress_bar(i-1, td%max_iter-1)
    end do
    write(stdout, '(1x)')

    ! calculate overlap
    call overlap()

    ! and now backward
    message(1) = "Info: Propagating backward"
    call write_info(1)

    call read_state(psi_f, "wf.final")
    call oct_progress_bar(-1, td%max_iter-1)
    td%dt = - td%dt
    h%lasers(1)%dt = td%dt
    functional = M_ZERO
    do i = td%max_iter-1, 0, -1
      call prop_iter2(i)
      call oct_progress_bar(td%max_iter-1-i, td%max_iter-1)
    end do
    td%dt = - td%dt
    write(stdout, '(1x)')

    ! write new field to file
    call io_assign(iunit)
    open(iunit, file='opt-control/laser', status='unknown')
    do i = 0, 2*td%max_iter
      write(iunit, '(4es20.12)') i*td%dt/M_TWO, laser_f(:, i)
    end do
    call io_close(iunit)

  end do
    
  ! clean up
  td%v_old => v_old_i
  nullify(h%lasers(1)%numerical)
  deallocate(v_old_f, laser_i, laser_f)
  call states_end(psi_f)

contains
  subroutine overlap()
    integer :: ik, p, dim, i
    real(r8) :: d1, j

    message(1) = "Overlap between wavefunctions"
    call write_info(1)
    do ik = 1, psi_i%nik
      do p  = psi_i%st_start, psi_i%st_end
        d1 = M_z0
        do dim = 1, psi_i%nik
          do i = 1, sys%m%np
            d1 = d1 + conjg(psi_i%zpsi(i, dim, p, ik))*psi_f%zpsi(i, dim, p, ik)
          end do
        end do
        d1 = d1*sys%m%vol_pp
        write(message(1), '(6x,i3,x,i3,a,2f16.10)') ik, p, " => ", abs(d1)**2, abs(d1)**2-alpha*functional
        call write_info(1)
      end do
    end do
  end subroutine overlap

  subroutine prop_iter1(iter)
    integer, intent(in) :: iter

    ! new electric field
    if(iter < td%max_iter) call update_field(iter, laser_i)

    ! psi_f
    td%v_old => v_old_f
    h%lasers(1)%numerical => laser_f
    call zcalcdens(psi_f, sys%m%np, sys%st%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_f, sys)
    call td_rti(h, sys%m, psi_f, sys, td, abs(iter*td%dt))

    ! psi_i
    td%v_old => v_old_i
    h%lasers(1)%numerical => laser_i
    call zcalcdens(psi_i, sys%m%np, sys%st%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_i, sys)
    call td_rti(h, sys%m, psi_i, sys, td, abs(iter*td%dt))

  end subroutine prop_iter1

  subroutine prop_iter2(iter)
    integer, intent(in) :: iter

    ! new electric field
    if(iter > 0) call update_field(iter, laser_f)

    ! psi_i
    td%v_old => v_old_i
    h%lasers(1)%numerical => laser_i
    call zcalcdens(psi_i, sys%m%np, sys%st%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_i, sys)
    call td_rti(h, sys%m, psi_i, sys, td, abs(iter*td%dt))

    ! psi_f
    td%v_old => v_old_f
    h%lasers(1)%numerical => laser_f
    call zcalcdens(psi_f, sys%m%np, sys%st%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_f, sys)
    call td_rti(h, sys%m, psi_f, sys, td, abs(iter*td%dt))

  end subroutine prop_iter2

  subroutine update_field(iter, l)
    integer, intent(in) :: iter
    real(r8), intent(inout) :: l(0:td%max_iter, conf%dim)

    complex(r8), allocatable :: grad(:,:)
    complex(r8) :: d1, d2(conf%dim), d3(conf%dim)
    integer :: ik, p, dim, i, j

    d1 = M_z0; d2 = M_z0; d3 = M_z0
    allocate(grad(3, sys%m%np))
    do ik = 1, psi_i%nik
      do p  = psi_i%st_start, psi_i%st_end
        do dim = 1, psi_i%nik
          call zmesh_derivatives(sys%m, psi_i%zpsi(:, dim, p, ik), grad=grad)
!          do i = 1, sys%m%np
!            write(70,*) sys%m%Lxyz(1, i)*sys%m%h(1), real(psi_i%zpsi(i, dim, p, ik)), real(grad(1, i))
!          end do
!          stop
          do i = 1, sys%m%np
            d1 = d1 + conjg(psi_i%zpsi(i, dim, p, ik))*psi_f%zpsi(i, dim, p, ik)
            d2(1:conf%dim) = d2(1:conf%dim) + &
                 conjg(psi_f%zpsi(i, dim, p, ik))*sys%m%Lxyz(1:conf%dim, i)*psi_i%zpsi(i, dim, p, ik)
            d3(1:conf%dim) = d3(1:conf%dim) + &
                 conjg(psi_f%zpsi(i, dim, p, ik))*grad(1:conf%dim, i)
          end do
        end do
      end do
    end do

    d1 = d1*sys%m%vol_pp
    d2(1:conf%dim) = -d2(1:conf%dim)*sys%m%h(1:conf%dim)*sys%m%vol_pp
    d3 = -d3*sys%m%vol_pp * M_zI*td%dt/M_TWO

    l(2*iter, 1:conf%dim) = -aimag(d1*d2(1:conf%dim))/alpha
    functional = functional + sum(l(iter, 1:conf%dim)**2)*abs(td%dt)

    ! extrapolate to t+-dt/2
    i = int(sign(M_ONE, td%dt))
    l(  i, 1:conf%dim) = l(0, 1:conf%dim)
    l(2*i, 1:conf%dim) = l(0, 1:conf%dim)
    
    l(2*iter+  i, 1:conf%dim) = M_HALF*(M_THREE*l(2*iter, 1:conf%dim) -       l(2*iter-2*i, 1:conf%dim))
    l(2*iter+2*i, 1:conf%dim) = M_HALF*( M_FOUR*l(2*iter, 1:conf%dim) - M_TWO*l(2*iter-2*i, 1:conf%dim))

    !print *, iter, l(2*iter, 1:conf%dim), l(2*iter+1, 1:conf%dim)

    !l(2*iter+i, 1:conf%dim)   = l(2*iter, 1:conf%dim) - aimag(d1*d3(1:conf%dim))/alpha
    !l(2*iter+2*i, 1:conf%dim) = l(2*iter, 1:conf%dim) - M_TWO*aimag(d1*d3(1:conf%dim))/alpha

  end subroutine update_field

  subroutine prop_psi_f()
    message(1) = "Info: Backward propagating Psi_f"
    call write_info(1)

    ! read final state
    call read_state(psi_f, "wf.final")

    ! setup the hamiltonian
    call zcalcdens(psi_f, sys%m%np, psi_f%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_f, sys)
    
    ! setup start of the propagation
    do i = 1, sys%st%nspin
      td%v_old(:, i, 2) = h%Vhartree(:) + h%Vxc(:, i)
    end do
    td%v_old(:, :, 3) = td%v_old(:, :, 2)

    !do i = 1, sys%m%np
    !  write(70,*) sys%m%Lxyz(1,i)*sys%m%h(1), real(psi_f%zpsi(i,1,1,1)), aimag(psi_f%zpsi(i,1,1,1))
    !end do
    td%dt = -td%dt
    call oct_progress_bar(-1, td%max_iter-1)
    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti(h, sys%m, psi_f, sys, td, abs(i*td%dt))
      ! update
      call zcalcdens(psi_f, sys%m%np, sys%st%rho, reduce=.true.)
      call zhamiltonian_setup(h, sys%m, psi_f, sys)

      call oct_progress_bar(td%max_iter-1-i, td%max_iter-1)
    end do
    td%dt = -td%dt
    message(1) = ""; call write_info(1)
    !do i = 1, sys%m%np
    !  write(71,*) sys%m%Lxyz(1,i)*sys%m%h(1), real(psi_f%zpsi(i,1,1,1)), aimag(psi_f%zpsi(i,1,1,1))
    !end do
  end subroutine prop_psi_f

  subroutine read_state(st, filename)
    type(states_type), intent(out) :: st
    character(len=*), intent(in) :: filename

    if(.not.zstates_load_restart("opt-control/"+trim(filename), sys%m, st)) then
      message(1) = "Could not read file 'opt-control/"+trim(filename)+"'"
      call write_fatal(1)
    end if
  end subroutine read_state

  subroutine init()
    integer :: t

    ! psi_i is initialized in system_init
    psi_i => sys%st
    v_old_i => td%v_old;
    
    ! now we initialize psi_f. This will repeat some stuff
    call states_init(psi_f, sys%m, sys%val_charge)
    if(sys%nlcc) then
      allocate(psi_f%rho_core(sys%m%np))
      psi_f%rho_core(sys%m%np) = psi_i%rho_core(sys%m%np)
    end if
    psi_f%st_start = psi_i%st_start
    psi_f%st_end = psi_i%st_end
    allocate(psi_f%zpsi(0:sys%m%np, psi_f%dim, psi_f%st_start:psi_f%st_end, psi_f%nik))
    allocate(v_old_f(sys%m%np, psi_f%nspin, 3))
    
    ! prepare the initial laser
    if(h%no_lasers.ne.0) then
      message(1) = "Please turn off any lasers for optimum control"
      call write_fatal(1)
    end if
    h%no_lasers = 1
    allocate(h%lasers(1))
    h%lasers(1)%envelope = 99 ! internal type
    h%lasers(1)%dt = td%dt
    allocate(laser_i(conf%dim, 0:2*td%max_iter), laser_f(conf%dim, 0:2*td%max_iter))
    laser_i = M_ZERO; 
    !laser_f = M_ZERO;
    laser_f = 1e-4_r8*exp(-((real(i)-real(td%max_iter)/M_THREE)*8._r8/real(td%max_iter))**2)
    h%lasers(1)%numerical => laser_f

    ! read alpha from input
    call oct_parse_double(C_string("OptControlAlpha"), M_ONE, alpha)

  end subroutine init

end subroutine opt_control_run

end module opt_control
