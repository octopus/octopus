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
!!
!! $Id$

#include "global.h"

module opt_control_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_m
  use lib_oct_parser_m
  use io_m
  use mesh_m
  use grid_m
  use states_m
  use system_m
  use restart_m
  use v_ks_m
  use hamiltonian_m
  use timedep_m
  use td_rti_m

  implicit none

  private
  public :: opt_control_run

contains

  ! ---------------------------------------------------------
  subroutine opt_control_run(sys, h)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h

    type(td_t)                :: td
    type(grid_t),     pointer :: gr   ! some shortcuts
    type(states_t),   pointer :: st

    type(states_t), pointer :: psi_i
    type(states_t)          :: psi_f
    FLOAT, pointer :: v_old_i(:,:,:), v_old_f(:,:,:)
    FLOAT, pointer :: laser_i(:,:), laser_f(:,:)

    integer :: i, ctr_iter, ctr_iter_max
    FLOAT :: eps, alpha, overlap, functional, old_functional, laser_init
    character(len=80) :: filename

    call init_()

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
      call zstates_calc_dens(psi_i, NP, psi_i%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_i, calc_eigenval=.true.)
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
      call zstates_calc_dens(psi_f, NP, psi_f%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_f, calc_eigenval=.true.)
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
    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine prop_iter1(iter)
      integer, intent(in) :: iter

      ! new electric field
      call update_field(iter-1, laser_i)

      ! psi_f
      td%tr%v_old => v_old_f
      h%ep%lasers(1)%numerical => laser_f
      call zstates_calc_dens(psi_f, NP, st%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_f, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi_f, td%tr, abs(iter*td%dt), abs(td%dt))

      ! psi_i
      td%tr%v_old => v_old_i
      h%ep%lasers(1)%numerical => laser_i
      call zstates_calc_dens(psi_i, NP, st%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_i, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi_i, td%tr, abs(iter*td%dt), abs(td%dt))

    end subroutine prop_iter1


    ! ---------------------------------------------------------
    subroutine prop_iter2(iter)
      integer, intent(in) :: iter

      ! new electric field
      call update_field(iter+1, laser_f)

      ! psi_i
      td%tr%v_old => v_old_i
      h%ep%lasers(1)%numerical => laser_i
      call zstates_calc_dens(psi_i, NP, st%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_i, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi_i, td%tr, abs(iter*td%dt), abs(td%dt))

      ! psi_f
      td%tr%v_old => v_old_f
      h%ep%lasers(1)%numerical => laser_f
      call zstates_calc_dens(psi_f, NP, st%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_f, calc_eigenval=.true.)
      call td_rti_dt(sys%ks, h, gr, psi_f, td%tr, abs(iter*td%dt), abs(td%dt))

    end subroutine prop_iter2


    ! ---------------------------------------------------------
    subroutine update_field(iter, l)
      integer, intent(in)    :: iter
      FLOAT,   intent(inout) :: l(0:2*td%max_iter, NDIM)

      CMPLX :: d1, d2(MAX_DIM)
      integer :: ik, p, dim, i

      d1 = M_z0; d2 = M_z0
      do ik = 1, psi_i%d%nik
        do p  = psi_i%st_start, psi_i%st_end
          do dim = 1, psi_i%d%dim
            do i = 1, NP
              d1 = d1 + conjg(psi_i%zpsi(i, dim, p, ik))*psi_f%zpsi(i, dim, p, ik) * gr%m%vol_pp(i)
              d2(1:NDIM) = d2(1:NDIM) - &
                conjg(psi_f%zpsi(i, dim, p, ik))*gr%m%x(i, 1:NDIM)*psi_i%zpsi(i, dim, p, ik) * gr%m%vol_pp(i)
            end do
          end do
        end do
      end do

      l(2*iter, 1:NDIM) = -aimag(d1*d2(1:NDIM))/alpha
      functional = functional + sum(l(2*iter, 1:NDIM)**2)*abs(td%dt)

      ! extrapolate to t+-dt/2
      i = int(sign(M_ONE, td%dt))
      if(iter==0.or.iter==td%max_iter) then
        l(2*iter+  i, 1:NDIM) = l(2*iter, 1:NDIM)
        l(2*iter+2*i, 1:NDIM) = l(2*iter, 1:NDIM)
      else
        l(2*iter+  i, 1:NDIM) = M_HALF*(M_THREE*l(2*iter, 1:NDIM) -       l(2*iter-2*i, 1:NDIM))
        l(2*iter+2*i, 1:NDIM) = M_HALF*( M_FOUR*l(2*iter, 1:NDIM) - M_TWO*l(2*iter-2*i, 1:NDIM))
      end if

    end subroutine update_field


    ! ---------------------------------------------------------
    subroutine prop_psi_f()
      message(1) = "Info: Backward propagating Psi_f"
      call write_info(1)

      ! read final state
      call read_state(psi_f, "wf.final")

      ! setup the hamiltonian
      call zstates_calc_dens(psi_f, NP, psi_f%rho)
      call zv_ks_calc(gr, sys%ks, h, psi_f, calc_eigenval=.true.)

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
        call td_rti_dt(sys%ks, h, gr, psi_f, td%tr, abs(i*td%dt), abs(td%dt))
        ! update
        call zstates_calc_dens(psi_f, NP, st%rho)
        call zv_ks_calc(gr, sys%ks, h, psi_f, calc_eigenval=.true.)

        call loct_progress_bar(td%max_iter-1-i, td%max_iter-1)
      end do
      td%dt = -td%dt
      message(1) = ""; call write_info(1)
    end subroutine prop_psi_f


    ! ---------------------------------------------------------
    subroutine calc_overlap()
      integer :: ik, p

      message(1) = "Overlap between wavefunctions"
      call write_info(1)
      do ik = 1, psi_i%d%nik
        do p  = psi_i%st_start, psi_i%st_end
          ! WARNING gives garbage when calculated through zstates_dotp
          overlap = zstates_dotp(gr%m, psi_i%d%dim, psi_i%zpsi(:,:, p, ik), psi_f%zpsi(:,:, p, ik))
          functional = overlap - alpha*functional
          write(message(1), '(6x,i3,1x,i3,a,f16.10,a,f16.10)') ik, p, " => overlap:", overlap, "  functional: " , functional
          call write_info(1)
        end do
      end do
    end subroutine calc_overlap


    ! ---------------------------------------------------------
    subroutine read_state(st, filename)
      type(states_t),   intent(out) :: st
      character(len=*), intent(in)  :: filename

      integer :: ierr

      call zrestart_read('opt-control/'//trim(filename), st, gr, ierr)
      if(ierr.ne.0) then
        message(1) = "Unsuccesfull read of states in 'opt-control/" // trim(filename) // "'"
        call write_fatal(1)
      end if
    end subroutine read_state


    ! ---------------------------------------------------------
    subroutine write_field(filename, las)
      character(len=*), intent(in) :: filename
      FLOAT,            intent(in) :: las(1:NDIM,0:2*td%max_iter)

      integer :: i, iunit

      iunit = io_open(filename, action='write')
      do i = 0, 2*td%max_iter
        write(iunit, '(4es20.12)') i*td%dt/M_TWO, las(:, i)
      end do
      call io_close(iunit)

    end subroutine write_field


    ! ---------------------------------------------------------
    subroutine output()
      integer :: iunit

      iunit = io_open('opt-control/info', action='write')
      write(iunit, '(a,i4)')    'Iterations = ', ctr_iter
      write(iunit, '(a,f14.8)') 'Overlap    = ', overlap
      write(iunit, '(a,f14.8)') 'Functional = ', functional
      call io_close(iunit)

      ! should output wavefunctions ;)
      call zstates_output(psi_i, gr, 'opt-control', sys%outp)
    end subroutine output


    ! ---------------------------------------------------------
    subroutine init_()
      ! some shortcuts
      gr  => sys%gr
      st  => sys%st

      call td_init(gr, td, sys%st, sys%outp)

      ! allocate memory
      ALLOCATE(st%zpsi(NP, st%d%dim, st%st_start:st%st_end, st%d%nik), NP*st%d%dim*(st%st_end-st%st_start+1)*st%d%nik)

      ! psi_i is initialized in system_init
      psi_i => st
      v_old_i => td%tr%v_old

      ! now we initialize psi_f. This will repeat some stuff
      call states_init(psi_f, gr)
      if(h%ep%nvnl > 0) then
        ALLOCATE(psi_f%rho_core(NP), NP)
        psi_f%rho_core(NP) = psi_i%rho_core(NP)
      end if
      psi_f%st_start = psi_i%st_start
      psi_f%st_end = psi_i%st_end
      ALLOCATE(psi_f%zpsi(NP, psi_f%d%dim, psi_f%st_start:psi_f%st_end, psi_f%d%nik), NP*psi_f%d%dim*(psi_f%st_end-psi_f%st_start+1)*psi_f%d%nik)
      ALLOCATE(v_old_f(NP, psi_f%d%nspin, 3), NP*psi_f%d%nspin*3)

      ! prepare the initial laser
      if(h%ep%no_lasers.ne.0) then
        message(1) = "Please turn off any lasers for optimum control"
        call write_fatal(1)
      end if
      h%ep%no_lasers = 1
      ALLOCATE(h%ep%lasers(1), 1)
      h%ep%lasers(1)%envelope = 99 ! internal type
      h%ep%lasers(1)%dt = td%dt
      ALLOCATE(laser_i(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))
      ALLOCATE(laser_f(NDIM, 0:2*td%max_iter), NDIM*(2*td%max_iter+1))

      ! read parameters from input
      ! we assume atomic units
      call loct_parse_float(check_inp('OptControlInitLaser'), M_ZERO, laser_init)
      laser_i = laser_init
      laser_f = laser_init
      call loct_parse_float(check_inp('OptControlAlpha'), M_ONE, alpha)
      call loct_parse_float(check_inp('OptControlEps'), CNST(1e-3), eps)
      call loct_parse_int(check_inp('OptControlMaxIter'), 10, ctr_iter_max)
      if(ctr_iter_max < 0.and.eps<M_ZERO) then
        message(1) = "OptControlMaxIter and OptControlEps can not be both <0"
        call write_fatal(1)
      end if

      if(ctr_iter_max < 0) ctr_iter_max = huge(ctr_iter_max)

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      call td_end(td)
    end subroutine end_

  end subroutine opt_control_run

end module opt_control_m
