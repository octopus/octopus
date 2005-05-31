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
use lib_oct
use lib_oct_parser
use geometry
use mesh
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
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
use PES
#endif

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
  integer           :: lmax           ! maximum multipole moment to output
  FLOAT             :: lmm_r          ! radius of the sphere used to compute the local magnetic moments

  !variables controlling the output
  logical           :: out_multip     ! multipoles
  logical           :: out_coords     ! coordinates
  logical           :: out_gsp        ! projection onto the ground state.
  logical           :: out_acc        ! electronic acceleration
  logical           :: out_laser      ! laser field
  logical           :: out_energy     ! several components of the electronic energy
  logical           :: out_proj       ! projection onto the GS KS eigenfunctions
  logical           :: out_angular    ! total angular momentum
  logical           :: out_spin       ! total spin
  logical           :: out_magnets    ! local magnetic moments

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
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h
  logical,                intent(inout) :: fromScratch

  type(td_type)                :: td
  type(mesh_type),     pointer :: m   ! some shortcuts
  type(states_type),   pointer :: st
  type(geometry_type), pointer :: geo

  logical :: stopping
  integer :: i, ii, j, idim, ist, ik
  integer(POINTER_SIZE) :: out_multip, out_coords, out_gsp, out_acc, &
       out_laser, out_energy, out_proj, out_angular, out_spin, out_magnets

  FLOAT, allocatable ::  x1(:,:), x2(:,:), f1(:,:) ! stuff for verlet
  FLOAT :: etime

  call init_()
  ierr = init_wfs()
  if(ierr.ne.0) then
    call end_()
    return
  end if
  call init_iter_output()

  ! Calculate initial forces and kinetic energy
  if(td%move_ions > 0) then
    if(td%iter > 0) then
      call td_read_nbo()
      call epot_generate(h%ep, m, st, geo, h%reltype)
      geo%eii = ion_ion_energy(geo)
      h%eii = geo%eii
    end if

    call zepot_forces(h%ep, m, st, geo, td%iter*td%dt, reduce_=.true.)
    geo%kinetic_energy = kinetic_energy(geo)
    select case(td%move_ions)
      case(NORMAL_VERLET)
        allocate(x1(geo%natoms, conf%dim), x2(geo%natoms, conf%dim))
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
        allocate(f1(geo%natoms, conf%dim))
    end select
  endif

  if(td%iter == 0) call td_run_zero_iter()
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
          call epot_generate(h%ep, m, st, geo, h%reltype)
      else
          call epot_generate(h%ep, m, st, geo, h%reltype, fast_generation = .true.)
      endif
      geo%eii = ion_ion_energy(geo)
      h%eii = geo%eii
    endif

    ! time iterate wavefunctions
    call td_rti_dt(sys%ks, h, m, sys%f_der, st, td%tr, i*td%dt, td%dt)

    ! mask function?
    if(h%ab == MASK_ABSORBING) then
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            st%zpsi(1:m%np, idim, ist, ik) = st%zpsi(1:m%np, idim, ist, ik) * &
                 (M_ONE - h%ab_pot(1:m%np))
          end do
        end do
      end do
    end if

    ! update density
    call zstates_calc_dens(st, m%np, st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call zh_calc_vhxc(sys%ks, h, m, sys%f_der, st, calc_eigenval=.true.)
    call hamiltonian_energy(h, st, geo%eii, -1, reduce=.true.)

    ! Recalculate forces, update velocities...
    if(td%move_ions > 0) then
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, geo%natoms
          f1(j, :) = geo%atom(j)%f(:)
        end do
      end if
      call zepot_forces(h%ep, m, st, geo, i*td%dt, reduce_=.true.)
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

    call iter_output()

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    call PES_doit(td%PESv, m, st, ii, td%dt, h%ab_pot)
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
      call td_write_data(i)
    end if

    if (stopping) exit
  end do

  call end_iter_output()
  call end_()

contains

  subroutine init_()
    call push_sub('td_run')

    ! some shortcuts
    m   => sys%m
    st  => sys%st
    geo => sys%geo

    call states_distribute_nodes(st)

    ! allocate memory
    allocate(st%zpsi(m%np, st%d%dim, st%st_start:st%st_end, st%d%nik))

    call td_init(td, m, st, geo, h, sys%outp)

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
      call zrestart_read(trim(tmpdir)//'restart_td', st, m, err, td%iter)
      if(err.ne.0) then
        message(1) = "Could not load "//trim(tmpdir)//"restart_td: Starting from scratch"
        call write_warning(1)
        
        fromScratch = .true.
      else
        ! read potential from previous interactions
        do i = 1, 2
          do is = 1, st%d%nspin
            write(filename,'(a,i2.2,i3.3)') trim(tmpdir)//'restart_td/vprev_', i, is
            call dinput_function(filename, m, td%tr%v_old(1:m%np, is, i), ierr)
          end do
        end do

      end if
    end if

    if(fromScratch) then
      call zrestart_read(trim(tmpdir)//'restart_gs', st, m, ierr)
      if(ierr.ne.0) then
        message(1) = "Could not load '"//trim(tmpdir)//"restart_gs': Starting from scratch"
        call write_warning(1)

        ierr = 1
      end if
    end if

    if(ierr==0) then
      call zstates_calc_dens(st, m%np, st%rho, reduce=.true.)
      call zh_calc_vhxc(sys%ks, h, m, sys%f_der, st, calc_eigenval=.true.)
      x = minval(st%eigenval(st%st_start, :))
#ifdef HAVE_MPI
      call MPI_BCAST(x, 1, MPI_FLOAT, 0, MPI_COMM_WORLD, i)
#endif
      call hamiltonian_span(h, minval(m%h(1:conf%dim)), x)
      call hamiltonian_energy(h, st, geo%eii, -1, reduce=.true.)
    end if

  end function init_wfs

  ! initialize buffers for output
  subroutine init_iter_output()
    character(len=256) :: filename

    integer :: first

    if (td%iter == 0) then
      first = 0
    else
      first = td%iter + 1
    end if

    if(mpiv%node==0) then
      if(td%out_multip)  call write_iter_init(out_multip,  first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/multipoles")))
      if(td%out_angular) call write_iter_init(out_angular, first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/angular")))
      if(td%out_spin)    call write_iter_init(out_spin,    first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/spin")))
      if(td%out_magnets) call write_iter_init(out_magnets, first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/magnetic_moments")))
      if(td%out_coords)  call write_iter_init(out_coords,  first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/coordinates")))
      if(td%out_gsp)     call write_iter_init(out_gsp,     first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/gs_projection")))
      if(td%out_acc)     call write_iter_init(out_acc,     first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/acceleration")))
      if(td%out_laser)   call write_iter_init(out_laser,   first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/laser")))
      if(td%out_energy)  call write_iter_init(out_energy,  first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/el_energy")))
    end if

    if(td%out_proj) then
      write(filename, '(i3.3)') mpiv%node
      call write_iter_init(out_proj, first, td%dt/units_out%time%factor, &
         trim(io_workpath(trim(current_label)//"td.general/projections."//trim(filename))) )
    end if

  end subroutine init_iter_output


  subroutine end_iter_output()
    ! close all buffers
    if(mpiv%node==0) then
      if(td%out_multip)  call write_iter_end(out_multip)
      if(td%out_angular) call write_iter_end(out_angular)
      if(td%out_spin)    call write_iter_end(out_spin)
      if(td%out_magnets) call write_iter_end(out_magnets)
      if(td%out_coords)  call write_iter_end(out_coords)
      if(td%out_gsp)     call write_iter_end(out_gsp)
      if(td%out_acc)     call write_iter_end(out_acc)
      if(td%out_laser)   call write_iter_end(out_laser)
      if(td%out_energy)  call write_iter_end(out_energy)
    end if
    if(td%out_proj)    call write_iter_end(out_proj)
    
  end subroutine end_iter_output


  subroutine iter_output()
    ! output multipoles
    if(td%out_multip) call td_write_multipole(out_multip, m, st, geo, td, i)

    ! output angular momentum
    if(td%out_angular) call td_write_angular(out_angular, m, sys%f_der, st, i)

    ! output spin
    if(td%out_spin) call td_write_spin(out_spin, m, st, i)

    ! output atoms magnetization
    if(td%out_magnets) call td_write_local_magnetic_moments(out_magnets, m, st, geo, td, i)
    
    ! output projections onto the GS KS eigenfunctions
!!$    if(td%out_proj) call td_write_proj(out_proj, m, st, u_st, i)
    
    ! output positions, vels, etc.
    if(td%out_coords) call td_write_nbo(out_coords, geo, td, i, geo%kinetic_energy, h%etot)
    
    ! If harmonic spectrum is desired, get the acceleration
    if(td%out_acc) call td_write_acc(out_acc, m, sys%f_der, st, geo, h, td, i)
    
    ! output laser field
    if(td%out_laser) call td_write_laser(out_laser, h, td, i)
    
    ! output electronic energy
    if(td%out_energy) call td_write_el_energy(out_energy, h, i)
  end subroutine iter_output


  subroutine td_run_zero_iter()

    call push_sub('td_run_zero_iter')

    ! create general subdir
    call io_mkdir(trim(current_label)//'td.general')

    if(td%out_multip)  call td_write_multipole(out_multip, m, st, geo, td, 0)
    if(td%out_angular) call td_write_angular(out_angular, m, sys%f_der, st, 0)
    if(td%out_spin)    call td_write_spin(out_spin, m, st, 0)
    if(td%out_magnets) call td_write_local_magnetic_moments(out_magnets, m, st, geo, td, 0)
!!$    if(td%out_proj)    call td_write_proj(out_proj, m, st, u_st, 0)

    call apply_delta_field()

    ! create files for output and output headers
    if(td%out_coords) call td_write_nbo(out_coords, geo, td, 0, geo%kinetic_energy, h%etot)    
    if(td%out_acc)    call td_write_acc(out_acc, m, sys%f_der, st, geo, h, td, 0)
    if(td%out_laser)  call td_write_laser(out_laser, h, td, 0)
    if(td%out_energy) call td_write_el_energy(out_energy, h, 0)

    call td_write_data(0)

    call td_rti_run_zero_iter(h, td%tr)

    call pop_sub()
  end subroutine td_run_zero_iter

  !!! Applies the delta function electric field E(t) = E_0 delta(t)
  !!! where E_0 = - k \hbar / e
  subroutine apply_delta_field()
    integer :: i, mode
    FLOAT   :: k, x(conf%dim)
    CMPLX   :: c(2), kick

    call push_sub('apply_delta_field')

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
      do i = 1, m%np
        call mesh_xyz(m, i, x)
        kick = M_zI * k * sum(x(1:conf%dim)*td%pol(1:conf%dim))

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
        geo%atom(i)%v(1:conf%dim) = geo%atom(i)%v(1:conf%dim) - &
             k*td%pol(1:conf%dim)*geo%atom(i)%spec%z_val / geo%atom(i)%spec%weight
      end do
    end if
    
    call pop_sub()
  end subroutine apply_delta_field

  subroutine td_read_nbo() ! reads the pos and vel from coordinates file
    integer :: i, iunit, record_length

    record_length = 100 + 3*geo%natoms*3*20
    call io_assign(iunit)
    open(unit = iunit, file = trim(current_label)//'td.general/coordinates', &
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
      read(iunit, '(3es20.12)', advance='no') geo%atom(i)%x(1:conf%dim)
      geo%atom(i)%x(:) = geo%atom(i)%x(:) * units_out%length%factor
    end do
    do i = 1, geo%natoms
      read(iunit, '(3es20.12)', advance='no') geo%atom(i)%v(1:conf%dim)
      geo%atom(i)%v(:) = geo%atom(i)%v(:) * units_out%velocity%factor
    end do
    do i = 1, geo%natoms
      read(iunit, '(3es20.12)', advance='no') geo%atom(i)%f(1:conf%dim)
      geo%atom(i)%f(:) = geo%atom(i)%f(:) * units_out%force%factor
    end do

    call io_close(iunit)
  end subroutine td_read_nbo

  subroutine td_write_data(iter)
    integer, intent(in) :: iter

    integer :: i, is, ierr
    character(len=256) :: filename
  
    call push_sub('td_write_data')
    
    ! first write resume file
    call zrestart_write(trim(tmpdir)//'restart_td', st, m, ierr, iter)
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
               filename, m, td%tr%v_old(1:m%np, is, i), M_ONE, ierr)

          if(ierr.ne.0) then
            write(message(1), '(3a)') 'Unsuccesfull write of "', trim(filename), '"'
            call write_fatal(1)
          end if
        end do
      end do
    end if

    ! calculate projection onto the ground state
    if(td%out_gsp) call td_write_gsp(out_gsp, m, st, td, iter)
    
    if(mpiv%node==0) then
      if(td%out_multip)  call write_iter_flush(out_multip)
      if(td%out_angular) call write_iter_flush(out_angular)
      if(td%out_spin)    call write_iter_flush(out_spin)
      if(td%out_magnets) call write_iter_flush(out_magnets)
      if(td%out_coords)  call write_iter_flush(out_coords)
      if(td%out_gsp)     call write_iter_flush(out_gsp)
      if(td%out_acc)     call write_iter_flush(out_acc)
      if(td%out_laser)   call write_iter_flush(out_laser)
      if(td%out_energy)  call write_iter_flush(out_energy)
    end if
    
    ! now write down the rest
    write(filename, '(a,i7.7)') trim(current_label)//"td.", iter  ! name of directory

    call zstates_output(st, m, sys%f_der, filename, sys%outp)
    if(sys%outp%what(output_geometry)) &
         call atom_write_xyz(filename, "geometry", geo)
    call hamiltonian_output(h, m, filename, sys%outp)
    
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    call PES_output(td%PESv, m, st, iter, sys%outp%iter, td%dt)
#endif
    
    call pop_sub()
  end subroutine td_write_data

end function td_run

#include "td_write.F90"
#include "td_init.F90"

end module timedep
