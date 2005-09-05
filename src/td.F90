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

module timedep
  use global
  use messages
  use lib_oct
  use lib_oct_parser
  use geometry
  use mesh
  use grid
  use mesh_function
  use functions
  use states
  use output
  use restart
  use lasers
  use v_ks
  use hamiltonian
  use external_pot
  use system
  use td_rti
  use td_write
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  use PES
#endif
  use grid

  implicit none

  private
  public :: td_type, &
       td_run, &
       td_init, &
       td_end

  type td_type
     type(td_rti_type) :: tr             ! contains the details of the time evolution
     FLOAT             :: dt             ! time step
     integer           :: max_iter       ! maximum number of iterations to perform
     integer           :: iter           ! the actual iteration
     integer           :: epot_regenerate! Every epot_regenerate, the external potential
                                         ! regenerated *exactly*.
     integer           :: move_ions      ! how do we move the ions?
     FLOAT             :: pol(3)         ! the direction of the polarization of the field

     type(td_write_type) :: write_handler

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
     type(PES_type) :: PESv
#endif

  end type td_type

  ! Parameters.
  integer, parameter :: STATIC_IONS     = 0,    &
       NORMAL_VERLET   = 3,  &
       VELOCITY_VERLET = 4

contains

  integer function td_run(sys, h, fromScratch) result(ierr)
    type(system_type), target, intent(inout) :: sys
    type(hamiltonian_type),    intent(inout) :: h
    logical,                   intent(inout) :: fromScratch

    type(td_type)                :: td
    type(grid_type),     pointer :: gr   ! some shortcuts
    type(states_type),   pointer :: st
    type(geometry_type), pointer :: geo

    logical :: stopping
    integer :: i, ii, j, idim, ist, ik

    FLOAT, allocatable ::  x1(:,:), x2(:,:), f1(:,:) ! stuff for verlet
    FLOAT :: etime

    call init_()
    ierr = init_wfs()
    if(ierr.ne.0) then
       call end_()
       return
    end if
    call td_write_init_iter(td%write_handler, td%iter, td%dt)

    ! Calculate initial forces and kinetic energy
    if(td%move_ions > 0) then
       if(td%iter > 0) then
          call td_read_nbo()
          call epot_generate(h%ep, gr%m, gr%sb, geo, st, h%reltype)
          geo%eii = ion_ion_energy(geo)
          h%eii = geo%eii
       end if

       call zepot_forces(gr, h%ep, st, td%iter*td%dt, reduce_=.true.)
       geo%kinetic_energy = kinetic_energy(geo)
       select case(td%move_ions)
       case(NORMAL_VERLET)
          allocate(x1(geo%natoms, NDIM), x2(geo%natoms, NDIM))
          do j = 1, geo%natoms
             if(geo%atom(j)%move) then
                x1(j, :) = geo%atom(j)%x(:) - td%dt*geo%atom(j)%v(:) + &
                     M_HALF * td%dt**2/geo%atom(j)%spec%weight * &
                     geo%atom(j)%f(:)
             else
                x1(j, :) = geo%atom(j)%x(:)
             endif
          enddo
       case(VELOCITY_VERLET)
          allocate(f1(geo%natoms, NDIM))
       end select
    endif

    
    if(td%iter == 0) then
       call apply_delta_field()
       call td_run_zero_iter()
    endif
    !call td_check_trotter(td, sys, h)
    td%iter = td%iter + 1

    write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
    call write_info(1)

    ii = 1
    stopping = .false.
    etime = loct_clock()
    do i = td%iter, td%max_iter
       if(clean_stop()) stopping = .true.
       ! Move the ions.
       if( td%move_ions > 0 ) then
          select case(td%move_ions)
          case(NORMAL_VERLET)
             x2 = x1
             do j = 1, geo%natoms
                if(geo%atom(j)%move) then
                   x1(j, :) = geo%atom(j)%x(:)
                   geo%atom(j)%x(:) = M_TWO*x1(j, :) - x2(j, :) + &
                        td%dt**2/geo%atom(j)%spec%weight * geo%atom(j)%f(:)
                   geo%atom(j)%v(:) = (geo%atom(j)%x(:) - x2(j, :)) / (M_TWO*td%dt)
                endif
             enddo
          case(VELOCITY_VERLET)
             do j=1, geo%natoms
                if(geo%atom(j)%move) then
                   geo%atom(j)%x(:) = geo%atom(j)%x(:) +  td%dt*geo%atom(j)%v(:) + &
                        M_HALF*td%dt**2/geo%atom(j)%spec%weight * geo%atom(j)%f(:)
                endif
             enddo
          end select

          if(mod(i, td%epot_regenerate) == 0) then
             call epot_generate(h%ep, gr%m, gr%sb, geo, st, h%reltype)
          else
             call epot_generate(h%ep, gr%m, gr%sb, geo, st, h%reltype, fast_generation = .true.)
          endif
          geo%eii = ion_ion_energy(geo)
          h%eii = geo%eii
       endif

       ! time iterate wavefunctions
       call td_rti_dt(sys%ks, h, gr, st, td%tr, i*td%dt, td%dt)

       ! mask function?
       if(h%ab == MASK_ABSORBING) then
          do ik = 1, st%d%nik
             do ist = st%st_start, st%st_end
                do idim = 1, st%d%dim
                   st%zpsi(1:NP, idim, ist, ik) = st%zpsi(1:NP, idim, ist, ik) * &
                        (M_ONE - h%ab_pot(1:NP))
                end do
             end do
          end do
       end if

       ! update density
       call zstates_calc_dens(st, NP, st%rho, reduce=.true.)

       ! update hamiltonian and eigenvalues (fermi is *not* called)
       call zv_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)
       call hamiltonian_energy(h, st, geo%eii, -1, reduce=.true.)

       ! Recalculate forces, update velocities...
       if(td%move_ions > 0) then
          if(td%move_ions == VELOCITY_VERLET) then
             do j = 1, geo%natoms
                f1(j, :) = geo%atom(j)%f(:)
             end do
          end if
          call zepot_forces(gr, h%ep, st, i*td%dt, reduce_=.true.)
          if(td%move_ions == VELOCITY_VERLET) then
             do j = 1, geo%natoms
                if(geo%atom(j)%move) then
                   geo%atom(j)%v(:) = geo%atom(j)%v(:) + &
                        td%dt/(M_TWO*geo%atom(j)%spec%weight) * &
                        (f1(j, :) + geo%atom(j)%f(:))
                end if
             end do
          end if
          geo%kinetic_energy = kinetic_energy(geo)
       end if

       call td_write_iter(td%write_handler, gr, st, h, geo, td%pol, td%dt, i)

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
       call PES_doit(td%PESv, gr%m, st, ii, td%dt, h%ab_pot)
#endif

       ! write info
       write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
            i*td%dt       / units_out%time%factor, &
            (h%etot + geo%kinetic_energy) / units_out%energy%factor, &
            (loct_clock() - etime)/1e6
       call write_info(1)
       etime = loct_clock()

       ! write down data
       ii = ii + 1
       if(ii==sys%outp%iter+1 .or. i == td%max_iter .or. stopping) then ! output
          if(i == td%max_iter) sys%outp%iter = ii - 1
          ii = 1
          call td_save_restart(i)
          call td_write_data(td%write_handler, gr, st, h, sys%outp, geo, td%dt, i)
       end if

       if (stopping) exit
    end do

    call td_write_end_iter(td%write_handler)
    call end_()

  contains

    subroutine init_()
      call push_sub('td.td_run')

      ! some shortcuts
      gr  => sys%gr
      geo => sys%gr%geo
      st  => sys%st

      call states_distribute_nodes(st)

      ! allocate memory
      allocate(st%zpsi(NP, st%d%dim, st%st_start:st%st_end, st%d%nik))

      call td_init(gr, td, st, h, sys%outp)

    end subroutine init_

    subroutine end_()
      ! free memory
      deallocate(st%zpsi)
      call td_end(td)

      call pop_sub()
    end subroutine end_

    integer function init_wfs() result(ierr)
      integer :: i, is, err
      character(len=50) :: filename
      FLOAT :: x

      ierr = 0
      if(.not.fromScratch) then
         call zrestart_read(trim(tmpdir)//'restart_td', st, gr%m, err, td%iter)
         if(err.ne.0) then
            message(1) = "Could not load "//trim(tmpdir)//"restart_td: Starting from scratch"
            call write_warning(1)

            fromScratch = .true.
         else
            ! read potential from previous interactions
            do i = 1, 2
               do is = 1, st%d%nspin
                  write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'restart_td/vprev_', i, is
                  call dinput_function(filename, gr%m, td%tr%v_old(1:NP, is, i), ierr)
               end do
            end do

         end if
      end if

      if(fromScratch) then
         call zrestart_read(trim(tmpdir)//'restart_gs', st, gr%m, ierr)
         if(ierr.ne.0) then
            message(1) = "Could not load '"//trim(tmpdir)//"restart_gs': Starting from scratch"
            call write_warning(1)

            ierr = 1
         end if
      end if

      if(ierr==0) then
         call zstates_calc_dens(st, NP, st%rho, reduce=.true.)
         call zv_ks_calc(gr, sys%ks, h, st, calc_eigenval=.true.)
         x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
         call MPI_BCAST(x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD, i)
#endif
         call hamiltonian_span(h, minval(gr%m%h(1:NDIM)), x)
         call hamiltonian_energy(h, st, geo%eii, -1, reduce=.true.)
      end if

    end function init_wfs

    subroutine td_run_zero_iter()
      call push_sub('td.td_run_zero_iter')

      call io_mkdir('td.general')
      call td_write_iter(td%write_handler, gr, st, h, geo, td%pol, td%dt, 0)
      call td_save_restart(0)
      call td_write_data(td%write_handler, gr, st, h, sys%outp, geo, td%dt, 0)
      call td_rti_run_zero_iter(h, td%tr)

      call pop_sub()
    end subroutine td_run_zero_iter

!!! Applies the delta function electric field E(t) = E_0 delta(t)
!!! where E_0 = - k \hbar / e
    subroutine apply_delta_field()
      integer :: i, mode
      FLOAT   :: k
      CMPLX   :: c(2), kick

      call push_sub('td.apply_delta_field')

!!! units are 1/length
      call loct_parse_float(check_inp('TDDeltaStrength'), M_ZERO, k)
      k = k / units_inp%length%factor
      call loct_parse_int(check_inp('TDDeltaStrengthMode'), 0, mode)
      select case (mode)
      case (0)
      case (1:2)
         if (st%d%ispin == UNPOLARIZED) then
            message(1) = "TDDeltaStrengthMode 1 or 2 can not be used when"
            message(2) = "performing unpolarized calculations"
            call write_fatal(1)
         end if
      case default
         write(message(1), '(a,i2,a)') "Input: '", mode, "' is not a valid TDDeltaStrengthMode"
         message(2) = '(0 <= TDDeltaStrengthMode <= 2 )'
         call write_fatal(2)
      end select

!!! The wave-functions at time delta t read
!!! psi(delta t) = psi(t) exp(i k x)
      if(k .ne. M_ZERO) then
         write(message(1),'(a,f11.6)')  'Info: Applying delta kick: k = ', k
         select case (mode)
         case (0)
            message(2) = "Info: Delta kick mode: Density mode"
         case (1)
            message(2) = "Info: Delta kick mode: Spin mode"
         case (2)
            message(2) = "Info: Delta kick mode: Density + Spin modes"
         case (3)
         end select
         call write_info(2)
         do i = 1, NP
            kick = M_zI * k * sum(gr%m%x(i, 1:NDIM)*td%pol(1:NDIM))

            select case (mode)
            case (0)
               c(1) = exp(kick)
               st%zpsi(i,:,:,:) = c(1) * st%zpsi(i,:,:,:)

            case (1)
               c(1) = exp(kick)
               c(2) = exp(-kick)
               select case (st%d%ispin)
               case (SPIN_POLARIZED)
                  do ik = 1, st%d%nik, 2
                     st%zpsi(i,:,:,ik)   = c(1) * st%zpsi(i,:,:,ik)
                     st%zpsi(i,:,:,ik+1) = c(2) * st%zpsi(i,:,:,ik+1)
                  end do
               case (SPINORS)
                  st%zpsi(i,1,:,:) = c(1) * st%zpsi(i,1,:,:)
                  st%zpsi(i,2,:,:) = c(2) * st%zpsi(i,2,:,:)
               end select

            case (2)
               c(1) = exp(M_TWO*kick)
               select case (st%d%ispin)
               case (SPIN_POLARIZED)
                  do ik = 1, st%d%nik, 2
                     st%zpsi(i,:,:,ik) = c(1) * st%zpsi(i,:,:,ik)
                  end do
               case (SPINORS)
                  st%zpsi(i,1,:,:) = c(1) * st%zpsi(i,1,:,:)
               end select

            end select
         end do
      end if

!!! the nuclei velocity will be changed by
!!! Delta v_z = ( Z*e*E_0 / M) = - ( Z*k*\hbar / M)
!!! where M and Z are the ionic mass and charge, respectively.
      if(td%move_ions > 0) then
         do i = 1, geo%natoms
            geo%atom(i)%v(1:NDIM) = geo%atom(i)%v(1:NDIM) - &
                 k*td%pol(1:NDIM)*geo%atom(i)%spec%z_val / geo%atom(i)%spec%weight
         end do
      end if

      call pop_sub()
    end subroutine apply_delta_field

    subroutine td_read_nbo() ! reads the pos and vel from coordinates file
      integer :: i, iunit, record_length

      record_length = 100 + 3*geo%natoms*3*20
      call io_assign(iunit)
      open(unit = iunit, file = 'td.general/coordinates', &
           action='read', status='old', recl = record_length)
      if(iunit < 0) then
         message(1) = "Could not open file 'td.general/coordinates'"
         message(2) = "Starting simulation from initial geometry"
         call write_warning(2)
         return
      end if

      read(iunit, *); read(iunit, *) ! skip header
      do i = 0, td%iter - 1
         read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
      end do
      read(iunit, '(88x)', advance='no') ! skip unrelevant information

      do i = 1, geo%natoms
         read(iunit, '(3es20.12)', advance='no') geo%atom(i)%x(1:NDIM)
         geo%atom(i)%x(:) = geo%atom(i)%x(:) * units_out%length%factor
      end do
      do i = 1, geo%natoms
         read(iunit, '(3es20.12)', advance='no') geo%atom(i)%v(1:NDIM)
         geo%atom(i)%v(:) = geo%atom(i)%v(:) * units_out%velocity%factor
      end do
      do i = 1, geo%natoms
         read(iunit, '(3es20.12)', advance='no') geo%atom(i)%f(1:NDIM)
         geo%atom(i)%f(:) = geo%atom(i)%f(:) * units_out%force%factor
      end do

      call io_close(iunit)
    end subroutine td_read_nbo

    subroutine td_save_restart(iter)
      integer, intent(in) :: iter

      integer :: i, is, ierr
      character(len=256) :: filename

      call push_sub('td.td_save_restart')

      ! first write resume file
      call zrestart_write(trim(tmpdir)//'restart_td', st, gr, ierr, iter)
      if(ierr.ne.0) then
         message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_td"'
         call write_fatal(1)
      end if

      ! write potential from previous interactions
      if(mpiv%node==0) then
         do i = 1, 2
            do is = 1, st%d%nspin
               write(filename,'(a6,i2.2,i3.3)') 'vprev_', i, is

               call doutput_function(restart_format, trim(tmpdir)//"restart_td", &
                    filename, gr%m, gr%sb, td%tr%v_old(1:NP, is, i), M_ONE, ierr)

               if(ierr.ne.0) then
                  write(message(1), '(3a)') 'Unsuccesfull write of "', trim(filename), '"'
                  call write_fatal(1)
               end if
            end do
         end do
      end if

      call pop_sub()
    end subroutine td_save_restart

  end function td_run

#include "td_init.F90"

end module timedep
