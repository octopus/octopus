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

module opt_control
use global
use lib_oct
use lib_oct_parser
use io
use states
use hamiltonian
use timedep
use td_rti

implicit none

contains

subroutine opt_control_run(td, m, st, val_charge, h, outp)
  type(td_type), intent(inout) :: td
  type(mesh_type), intent(in)  :: m
  type(states_type), intent(inout), target :: st
  FLOAT,                   intent(in) :: val_charge
  type(hamiltonian_type), intent(inout) :: h
  type(output_type), intent(in) :: outp

  type(states_type), pointer :: psi_i
  type(states_type)          :: psi_f
  FLOAT, pointer :: v_old_i(:,:,:), v_old_f(:,:,:)
  FLOAT, pointer :: laser_i(:,:), laser_f(:,:)

  integer :: i, ctr_iter, ctr_iter_max
  FLOAT :: eps, alpha, overlap, functional, old_functional, laser_init
  character(len=80) :: filename

  call init()

  ! first propagate psi_f to ti
  call prop_psi_f()

  old_functional = -CNST(1e10)
  ctr_iter = 1
  ctr_loop: do
    ! first propagate both states forward
    write(message(1), '(a,i3)') 'Info: Optimum control iteration #', ctr_iter
    message(2) = "Info: Propagating forward"
    call write_info(2)

    ! setup forward propagation
    call read_state(psi_i, "wf.initial")
    call zcalcdens(psi_i, m%np, psi_i%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_i, calc_eigenval=.true.)
    do i = 1, st%d%nspin
      v_old_i(:, i, 2) = h%vhxc(:, i)
    end do
    v_old_i(:, :, 3) = v_old_i(:, :, 2)
    
    do i = 2, 3
      v_old_f(:,:,i) = v_old_f(:,:,1) ! this one comes from the previous propagation
    end do

    call loct_progress_bar(-1, td%max_iter-1)
    h%ep%lasers(1)%dt = td%dt
    functional = M_ZERO
    do i = 1, td%max_iter
      call prop_iter1(i)
      call loct_progress_bar(i-1, td%max_iter-1)
    end do
    write(stdout, '(1x)')

    ! write new field to file
    write(filename,'(a,i3.3)') 'opt-control/laser.', ctr_iter
    call write_field(filename, laser_i)

    ! setup backward propagation
    call read_state(psi_f, "wf.final")
    call zcalcdens(psi_f, m%np, psi_f%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_f, calc_eigenval=.true.)
    do i = 1, st%d%nspin
      v_old_f(:, i, 2) = h%vhxc(:, i)
    end do
    v_old_f(:, :, 3) = v_old_f(:, :, 2)
    
    do i = 2, 3
      v_old_i(:,:,i) = v_old_i(:,:,1) ! this one comes the previous propagation
    end do

    ! calculate overlap
    call calc_overlap()
    if((ctr_iter==ctr_iter_max).or.(eps>M_ZERO.and.abs(functional-old_functional) < eps)) &
         exit ctr_loop
    ctr_iter = ctr_iter + 1
    old_functional = functional
    
    ! and now backward
    message(1) = "Info: Propagating backward"
    call write_info(1)

    call loct_progress_bar(-1, td%max_iter-1)
    td%dt = -td%dt
    h%ep%lasers(1)%dt = td%dt
    functional = M_ZERO
    do i = td%max_iter-1, 0, -1
      call prop_iter2(i)
      call loct_progress_bar(td%max_iter-1-i, td%max_iter-1)
    end do
    td%dt = -td%dt
    write(stdout, '(1x)')
  end do ctr_loop

  ! output some useful information
  call output()

  ! clean up
  td%tr%v_old => v_old_i
  nullify(h%ep%lasers(1)%numerical)
  deallocate(v_old_f, laser_i, laser_f)
  call states_end(psi_f)

contains
  subroutine prop_iter1(iter)
    integer, intent(in) :: iter

    ! new electric field
    call update_field(iter-1, laser_i)

    ! psi_f
    td%tr%v_old => v_old_f
    h%ep%lasers(1)%numerical => laser_f
    call zcalcdens(psi_f, m%np, st%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_f, calc_eigenval=.true.)
    call td_rti_dt(h, m, psi_f, td%tr, abs(iter*td%dt), abs(td%dt))

    ! psi_i
    td%tr%v_old => v_old_i
    h%ep%lasers(1)%numerical => laser_i
    call zcalcdens(psi_i, m%np, st%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_i, calc_eigenval=.true.)
    call td_rti_dt(h, m, psi_i, td%tr, abs(iter*td%dt), abs(td%dt))

  end subroutine prop_iter1

  subroutine prop_iter2(iter)
    integer, intent(in) :: iter

    ! new electric field
    call update_field(iter+1, laser_f)

    ! psi_i
    td%tr%v_old => v_old_i
    h%ep%lasers(1)%numerical => laser_i
    call zcalcdens(psi_i, m%np, st%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_i, calc_eigenval=.true.)
    call td_rti_dt(h, m, psi_i, td%tr, abs(iter*td%dt), abs(td%dt))

    ! psi_f
    td%tr%v_old => v_old_f
    h%ep%lasers(1)%numerical => laser_f
    call zcalcdens(psi_f, m%np, st%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_f, calc_eigenval=.true.)
    call td_rti_dt(h, m, psi_f, td%tr, abs(iter*td%dt), abs(td%dt))

  end subroutine prop_iter2

  subroutine update_field(iter, l)
    integer, intent(in) :: iter
    FLOAT, intent(inout) :: l(0:2*td%max_iter, conf%dim)

    CMPLX, allocatable :: grad(:,:)
    CMPLX :: d1, d2(conf%dim)
    integer :: ik, p, dim, i

    d1 = M_z0; d2 = M_z0
    allocate(grad(3, m%np))
    do ik = 1, psi_i%nik
      do p  = psi_i%st_start, psi_i%st_end
        do dim = 1, psi_i%dim
          do i = 1, m%np
            d1 = d1 + conjg(psi_i%zpsi(i, dim, p, ik))*psi_f%zpsi(i, dim, p, ik)
            d2(1:conf%dim) = d2(1:conf%dim) + &
                 conjg(psi_f%zpsi(i, dim, p, ik))*m%Lxyz(1:conf%dim, i)*psi_i%zpsi(i, dim, p, ik)
          end do
        end do
      end do
    end do

    d1 = d1*m%vol_pp
    d2(1:conf%dim) = -d2(1:conf%dim)*m%h(1:conf%dim)*m%vol_pp

    l(2*iter, 1:conf%dim) = -aimag(d1*d2(1:conf%dim))/alpha
    functional = functional + sum(l(2*iter, 1:conf%dim)**2)*abs(td%dt)

    ! extrapolate to t+-dt/2
    i = int(sign(M_ONE, td%dt))
    if(iter==0.or.iter==td%max_iter) then
      l(2*iter+  i, 1:conf%dim) = l(2*iter, 1:conf%dim)
      l(2*iter+2*i, 1:conf%dim) = l(2*iter, 1:conf%dim)
    else
      l(2*iter+  i, 1:conf%dim) = M_HALF*(M_THREE*l(2*iter, 1:conf%dim) -       l(2*iter-2*i, 1:conf%dim))
      l(2*iter+2*i, 1:conf%dim) = M_HALF*( M_FOUR*l(2*iter, 1:conf%dim) - M_TWO*l(2*iter-2*i, 1:conf%dim))
    end if

  end subroutine update_field

  subroutine prop_psi_f()
    message(1) = "Info: Backward propagating Psi_f"
    call write_info(1)

    ! read final state
    call read_state(psi_f, "wf.final")

    ! setup the hamiltonian
    call zcalcdens(psi_f, m%np, psi_f%rho, reduce=.true.)
    call zh_calc_vhxc(h, m, psi_f, calc_eigenval=.true.)
    
    ! setup start of the propagation
    do i = 1, st%d%nspin
      v_old_f(:, i, 2) = h%vhxc(:, i)
    end do
    v_old_f(:, :, 3) = v_old_f(:, :, 2)
    
    h%ep%lasers(1)%numerical => laser_f
    td%tr%v_old => v_old_f

    td%dt = -td%dt
    h%ep%lasers(1)%dt = td%dt
    call loct_progress_bar(-1, td%max_iter-1)
    do i = td%max_iter-1, 0, -1
      ! time iterate wavefunctions
      call td_rti_dt(h, m, psi_f, td%tr, abs(i*td%dt), abs(td%dt))
      ! update
      call zcalcdens(psi_f, m%np, st%rho, reduce=.true.)
      call zh_calc_vhxc(h, m, psi_f, calc_eigenval=.true.)

      call loct_progress_bar(td%max_iter-1-i, td%max_iter-1)
    end do
    td%dt = -td%dt
    message(1) = ""; call write_info(1)
  end subroutine prop_psi_f

  subroutine calc_overlap()
    integer :: ik, p

    message(1) = "Overlap between wavefunctions"
    call write_info(1)
    do ik = 1, psi_i%nik
      do p  = psi_i%st_start, psi_i%st_end
        ! WARNING gives garbage when calculated through zstates_dotp
        overlap = abs(sum(conjg(psi_i%zpsi(1:,:, p, ik))*psi_f%zpsi(1:,:, p, ik))*m%vol_pp)**2
        functional = overlap - alpha*functional
        write(message(1), '(6x,i3,1x,i3,a,f16.10,a,f16.10)') ik, p, " => overlap:", overlap, "  functional: " , functional
        call write_info(1)
      end do
    end do
  end subroutine calc_overlap

  subroutine read_state(st, filename)
    type(states_type), intent(out) :: st
    character(len=*), intent(in) :: filename

    if(.not.zstates_load_restart("opt-control/" // trim(filename), m, st)) then
      message(1) = "Could not read file 'opt-control/" // trim(filename) // "'"
      call write_fatal(1)
    end if
  end subroutine read_state

  subroutine write_field(filename, las)
    character(len=*), intent(in) :: filename
    FLOAT, intent(in) :: las(1:conf%dim,0:2*td%max_iter)

    integer :: i, iunit

    call io_assign(iunit)
    open(iunit, file=trim(filename), status='unknown')
    do i = 0, 2*td%max_iter
      write(iunit, '(4es20.12)') i*td%dt/M_TWO, las(:, i)
    end do
    call io_close(iunit)
  
  end subroutine write_field

  subroutine output()
    integer :: iunit

    call io_assign(iunit)
    open(iunit, file='opt-control/info', status='unknown')
    write(iunit, '(a,i4)')     'Iterations = ', ctr_iter
    write(iunit, '(a,f14.8)') 'Overlap    = ', overlap
    write(iunit, '(a,f14.8)') 'Functional = ', functional
    call io_close(iunit)

    ! should output wavefunctions ;)
    call zstates_output(psi_i, m, "opt-control", outp)
  end subroutine output

  subroutine init()
    ! psi_i is initialized in system_init
    psi_i => st
    v_old_i => td%tr%v_old
    
    ! now we initialize psi_f. This will repeat some stuff
    call states_init(psi_f, m, val_charge)
    if(h%ep%nvnl > 0) then
      allocate(psi_f%rho_core(m%np))
      psi_f%rho_core(m%np) = psi_i%rho_core(m%np)
    end if
    psi_f%st_start = psi_i%st_start
    psi_f%st_end = psi_i%st_end
    allocate(psi_f%zpsi(m%np, psi_f%dim, psi_f%st_start:psi_f%st_end, psi_f%nik))
    allocate(v_old_f(m%np, psi_f%nspin, 3))
    
    ! prepare the initial laser
    if(h%ep%no_lasers.ne.0) then
      message(1) = "Please turn off any lasers for optimum control"
      call write_fatal(1)
    end if
    h%ep%no_lasers = 1
    allocate(h%ep%lasers(1))
    h%ep%lasers(1)%envelope = 99 ! internal type
    h%ep%lasers(1)%dt = td%dt
    allocate(laser_i(conf%dim, 0:2*td%max_iter), laser_f(conf%dim, 0:2*td%max_iter))

    ! read parameters from input
    ! we assume atomic units
    call loct_parse_float("OptControlInitLaser", M_ZERO, laser_init)
    laser_i = laser_init
    laser_f = laser_init
    call loct_parse_float("OptControlAlpha", M_ONE, alpha)
    call loct_parse_float("OptControlEps", CNST(1e-3), eps)
    call loct_parse_int("OptControlMaxIter", 10, ctr_iter_max)
    if(ctr_iter_max < 0.and.eps<M_ZERO) then
      message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
      call write_fatal(1)
    end if

    if(ctr_iter_max < 0) ctr_iter_max = huge(ctr_iter_max)

  end subroutine init

end subroutine opt_control_run

end module opt_control
