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

  integer :: i
  real(r8) :: alpha

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
    do i = 1, td%max_iter
      call prop_iter(i)
    end do

    ! calculate overlap
    call overlap()

    ! and now backward
    td%dt = - td%dt
    do i = td%max_iter-1, 0, -1
      call prop_iter(i)
    end do
    td%dt = - td%dt
  end do
    
  ! clean up
  td%v_old => v_old_i
  deallocate(v_old_f)
  call states_end(psi_f)

contains
  subroutine overlap()
    integer :: ik, p, dim, i
    real(r8) :: d1

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
        write(message(1), '(6x,i3,x,i3,a,f10.4)') ik, p, " => ", d1
        call write_info(1)
      end do
    end do
  end subroutine overlap

  subroutine prop_iter(iter)
    integer, intent(in) :: iter

    ! psi_i
    if(iter.ne.0) then
      td%v_old => v_old_i
      call zcalcdens(psi_i, sys%m%np, sys%st%rho, reduce=.true.)
      call zhamiltonian_setup(h, sys%m, psi_i, sys)
      call td_rti(h, sys%m, psi_i, sys, td, iter*td%dt)
    else
      call read_state(psi_i, "wf.initial")
    end if

    ! psi_f
    if(iter.ne.td%max_iter) then
      td%v_old => v_old_f
      call zcalcdens(psi_f, sys%m%np, sys%st%rho, reduce=.true.)
      call zhamiltonian_setup(h, sys%m, psi_f, sys)
      call td_rti(h, sys%m, psi_f, sys, td, iter*td%dt)
    else
      call read_state(psi_f, "wf.final")
    end if

    ! new electric field
    call update_field(iter)

  end subroutine prop_iter

  subroutine update_field(iter)
    integer, intent(in) :: iter

    complex(r8) :: d1, d2(conf%dim)
    integer :: ik, p, dim, i, j

    d1 = M_z0
    d2 = M_z0
    do ik = 1, psi_i%nik
      do p  = psi_i%st_start, psi_i%st_end
        do dim = 1, psi_i%nik
          do i = 1, sys%m%np
            d1    = d1    + conjg(psi_i%zpsi(i, dim, p, ik))*psi_f%zpsi(i, dim, p, ik)
            d2(:) = d2(:) + psi_i%zpsi(i, dim, p, ik)*sys%m%Lxyz(:,i)*conjg(psi_f%zpsi(i, dim, p, ik))
          end do
        end do
      end do

      d1    = d1*sys%m%vol_pp
      d2(:) = d2(:)*sys%m%h(:)*sys%m%vol_pp
      h%lasers(1)%numerical(:, iter) = - aimag(d1*d2(:))/alpha
    end do
  end subroutine update_field

  subroutine init()
    real(r8) :: d

    ! psi_i is initialized in system_init
    psi_i => sys%st
    v_old_i => td%v_old;
    
    ! now we initialize psi_f. This will repeat some stuff
    d = 0._r8
    do i = 1, sys%natoms
      d = d - sys%atom(i)%spec%Z_val
    enddo
    call states_init(psi_f, sys%m, d)
    psi_f%nlcc = psi_i%nlcc
    if(psi_f%nlcc) allocate(psi_f%rho_core(sys%m%np))
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
    allocate(h%lasers(1)%numerical(3, 0:td%max_iter))
    h%lasers(1)%numerical(:,:) = 0._r8

    ! read alpha from input
    call oct_parse_double(C_string("OptControlAlpha"), 1._r8, alpha)
  end subroutine init

  subroutine prop_psi_f()
    ! read initial state
    call read_state(psi_f, "wf.final")

    ! setup the hamiltonian
    call zcalcdens(psi_f, sys%m%np, psi_f%rho, reduce=.true.)
    call zhamiltonian_setup(h, sys%m, psi_f, sys)
    
    ! setup start of the propagation
    do i = 1, sys%st%nspin
      td%v_old(:, i, 2) = h%Vhartree(:) + h%Vxc(:, i)
    end do
    td%v_old(:, :, 3) = td%v_old(:, :, 2)

    td%dt = -td%dt
    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti(h, sys%m, psi_f, sys, td, i*td%dt)

      ! update
      call zcalcdens(psi_f, sys%m%np, sys%st%rho, reduce=.true.)
      call zhamiltonian_setup(h, sys%m, psi_f, sys)
    end do
    td%dt = -td%dt
    
  end subroutine prop_psi_f

  subroutine read_state(st, filename)
    type(states_type), intent(out) :: st
    character(len=*), intent(in) :: filename

    if(.not.dstates_load_restart("opt-control/"+trim(filename), sys%m, st)) then
      message(1) = "Could not read file 'opt-control/"+trim(filename)+"'"
      call write_fatal(1)
    end if
  end subroutine read_state

end subroutine opt_control_run

end module opt_control
