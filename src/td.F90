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

module timedep
use global
use lib_oct_parser
use mesh
use atom
use states
use hamiltonian
use external_pot
use td_rti
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
use PES
#endif

implicit none

type td_type
  type(td_rti_type) :: tr             ! contains the details of the time evolution
  FLOAT          :: dt             ! time step
  integer           :: max_iter       ! maximum number of iterations to perform
  integer           :: iter           ! the actual iteration

  integer           :: move_ions      ! how do we move the ions?

  FLOAT          :: pol(3)         ! the direction of the polarization of the field
  integer           :: lmax           ! maximum multipole moment to output

  !variables controlling the output
  logical           :: out_multip     ! multipoles
  logical           :: out_coords     ! coordinates
  logical           :: out_gsp        ! projection onto the ground state.
  logical           :: out_acc        ! electronic acceleration
  logical           :: out_laser      ! laser field
  logical           :: out_energy     ! several components of the electronic energy
  logical           :: out_proj       ! projection onto the GS KS eigenfunctions
  logical           :: out_angular    ! total angular momentum.

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  type(PES_type) :: PESv
#endif

end type td_type

  ! Parameters.
  integer, parameter :: STATIC_IONS     = 0,    &
                        NORMAL_VERLET   = 3,  &
                        VELOCITY_VERLET = 4

contains

subroutine td_run(td, u_st, m, st, geo, h, outp)
  type(td_type),          intent(inout) :: td
  type(states_type),      intent(IN)    :: u_st
  type(mesh_type),        intent(IN)    :: m
  type(states_type),      intent(inout) :: st
  type(geometry_type),    intent(inout) :: geo
  type(hamiltonian_type), intent(inout) :: h
  type(output_type),      intent(inout) :: outp

  integer :: i, ii, j, idim, ist, ik
  integer(POINTER_SIZE) :: out_multip, out_coords, out_gsp, out_acc, &
       out_laser, out_energy, out_proj, out_angular

  FLOAT, allocatable ::  x1(:,:), x2(:,:), f1(:,:) ! stuff for verlet
  FLOAT :: etime
  character(len=100) :: filename

  call push_sub('td_run')
  
  if(mpiv%node==0) then
    write(filename, '(i3.3)') mpiv%node
    if(td%out_multip) &
         call write_iter_init(out_multip,  td%iter, td%dt/units_out%time%factor, "td.general/multipoles")
    if(td%out_angular) &
         call write_iter_init(out_angular, td%iter, td%dt/units_out%time%factor, "td.general/angular")
    if(td%out_coords) &
         call write_iter_init(out_coords,  td%iter, td%dt/units_out%time%factor, "td.general/coordinates")
    if(td%out_gsp) &
         call write_iter_init(out_gsp,     td%iter, td%dt/units_out%time%factor, "td.general/gs_projection")
    if(td%out_acc) &
         call write_iter_init(out_acc,     td%iter, td%dt/units_out%time%factor, "td.general/acceleration")
    if(td%out_laser) &
         call write_iter_init(out_laser,   td%iter, td%dt/units_out%time%factor, "td.general/laser")
    if(td%out_energy) &
         call write_iter_init(out_energy,  td%iter, td%dt/units_out%time%factor, "td.general/el_energy")
    if(td%out_proj) &
         call write_iter_init(out_proj,    td%iter, td%dt/units_out%time%factor, "td.general/projections."//trim(filename))
  end if

  ! Calculate initial forces and kinetic energy
  if(td%move_ions > 0) then
    if(td%iter > 0) then
      call td_read_nbo()
      call epot_generate(h%ep, m, st, geo, h%Vpsl, h%reltype)
      geo%eii = ion_ion_energy(geo)
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
  etime = loct_clock()
  do i = td%iter, td%max_iter
    if(clean_stop()) exit
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

      call epot_generate(h%ep, m, st, geo, h%Vpsl, h%reltype)
      geo%eii = ion_ion_energy(geo)
    endif

    ! time iterate wavefunctions
    call td_rti_dt(h, m, st, td%tr, i*td%dt, td%dt)

    ! mask function?
    if(h%ab == MASK_ABSORBING) then
      do ik = 1, st%nik
        do ist = st%st_start, st%st_end
          do idim = 1, st%dim
            st%zpsi(1:m%np, idim, ist, ik) = st%zpsi(1:m%np, idim, ist, ik) * &
                 (M_ONE - h%ab_pot(1:m%np))
          end do
        end do
      end do
    end if

    ! update density
    call zcalcdens(st, m%np, st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call zh_calc_vhxc(h, m, st, calc_eigenval=.true.)
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

    ! output multipoles
    if(td%out_multip) call td_write_multipole(out_multip, m, st, geo, td, i)

    ! output angular momentum
    if(td%out_angular) call td_write_angular(out_angular, m, st, td, i)

    ! output projections onto the GS KS eigenfunctions
    if(td%out_proj) call td_write_proj(out_proj, m, st, u_st, i)
    
    ! output positions, vels, etc.
    if(td%out_coords) call td_write_nbo(out_coords, geo, td, i, geo%kinetic_energy, h%etot)

    ! If harmonic spectrum is desired, get the acceleration
    if(td%out_acc) call td_write_acc(out_acc, m, st, geo, h, td, i)

    ! output laser field
    if(td%out_laser) call td_write_laser(out_laser, h, td, i)

    ! output electronic energy
    if(td%out_energy) call td_write_el_energy(out_energy, h, i)

#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    call PES_doit(td%PESv, m, st, ii, td%dt, h%ab_pot)
#endif

    ! write info
    write(message(1), '(i7,1x,2f14.6,f14.3, i10)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         (loct_clock() - etime)/1e6
    call write_info(1)
    etime = loct_clock()

    ! write down data
    ii = ii + 1
    if(ii==outp%iter+1 .or. i == td%max_iter) then ! output
      if(i == td%max_iter) outp%iter = ii - 1
      ii = 1

      call td_write_data(i)
    end if

  end do

  ! close all buffers
  if(mpiv%node==0) then
    if(td%out_multip)  call write_iter_end(out_multip)
    if(td%out_angular) call write_iter_end(out_angular)
    if(td%out_coords)  call write_iter_end(out_coords)
    if(td%out_gsp)     call write_iter_end(out_gsp)
    if(td%out_acc)     call write_iter_end(out_acc)
    if(td%out_laser)   call write_iter_end(out_laser)
    if(td%out_energy)  call write_iter_end(out_energy)
    if(td%out_proj)    call write_iter_end(out_proj)
  end if

contains

  subroutine td_run_zero_iter()

    call push_sub('td_run_zero_iter')

    ! create general subdir
    call loct_mkdir("td.general")

    if(td%out_multip)  call td_write_multipole(out_multip, m, st, geo, td, 0)
    if(td%out_angular) call td_write_angular(out_angular, m, st, td, 0)
    if(td%out_proj)    call td_write_proj(out_proj, m, st, u_st, 0)

    call apply_delta_field()

    ! create files for output and output headers
    if(td%out_coords) call td_write_nbo(out_coords, geo, td, 0, geo%kinetic_energy, h%etot)    
    if(td%out_acc)    call td_write_acc(out_acc, m, st, geo, h, td, 0)
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

    !!! units are 1/length
    call loct_parse_float("TDDeltaStrength", M_ZERO, k)
    k = k / units_inp%length%factor
    call loct_parse_int("TDDeltaStrengthMode", 0, mode)
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
      call write_info(1)
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
            do ik = 1, st%nik, 2
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
            do ik = 1, st%nik, 2
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
    
  end subroutine apply_delta_field

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

    character(len=50) :: filename
  
    call push_sub('td_write_data')
    
    ! first resume file
    write(filename, '(a,i3.3)') "tmp/restart.td.", mpiv%node
    call zstates_write_restart(trim(filename), m, st, iter=iter, &
                 v1=td%tr%v_old(:, :, 1), v2=td%tr%v_old(:, :, 2))
    
    ! calculate projection onto the ground state
    if(td%out_gsp) call td_write_gsp(out_gsp, m, st, td, iter)
    
    if(mpiv%node==0) then
      if(td%out_multip)  call write_iter_flush(out_multip)
      if(td%out_angular) call write_iter_flush(out_angular)
      if(td%out_coords)  call write_iter_flush(out_coords)
      if(td%out_gsp)     call write_iter_flush(out_gsp)
      if(td%out_acc)     call write_iter_flush(out_acc)
      if(td%out_laser)   call write_iter_flush(out_laser)
      if(td%out_energy)  call write_iter_flush(out_energy)
    end if
    
    ! now write down the rest
    write(filename, '(a,i7.7)') "td.", iter  ! name of directory
    call zstates_output(st, m, filename, outp)
    if(outp%what(output_geometry)) &
         call atom_write_xyz(filename, "geometry", geo)
    call hamiltonian_output(h, m, filename, outp)
    
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    call PES_output(td%PESv, m, st, iter, outp%iter, td%dt)
#endif
    
    call pop_sub()
  end subroutine td_write_data

end subroutine td_run

#include "td_write.F90"
#include "td_init.F90"

end module timedep
