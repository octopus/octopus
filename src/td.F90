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
  integer :: save_iter ! save every save_iter iterations
  integer :: iter      ! the actual iteration
  integer :: evolution_method ! which evolution method to use

  real(r8) :: dt            ! time step
  integer  :: move_ions     ! how do we move the ions?

  real(r8) :: delta_strength ! strength of the delta excitation
  real(r8) :: pol(3)         ! direction of the polarization of the efield

  integer :: lmax        ! maximum multipole moment to write

  integer :: gauge ! in which gauge shall we work in
                   ! 1 = length gauge
                   ! 2 = velocity gauge

  ! lasers stuff
  integer :: no_lasers ! number of laser pulses used
  logical :: output_laser ! write laser field
  type(laser_type), pointer :: lasers(:)

  ! absorbing boundaries
  integer  :: ab         ! do we have absorbing boundaries?
  real(r8) :: ab_width   ! width of the absorbing boundary
  real(r8) :: ab_height  ! height of the absorbing boundary
  real(r8), pointer :: ab_pot(:) ! where we store the ab potential

  ! occupational analysis
  logical :: occ_analysis ! do we perform occupational analysis?

  ! write harmonic spectrum
  logical :: harmonic_spectrum

#ifndef DISABLE_PES
  type(PES_type) :: PESv
#endif

  real(r8), pointer :: v_old1(:,:), v_old2(:,:)
end type td_type

  integer, parameter :: STATIC_IONS = 0,    &
                        NORMAL_VERLET = 3,  &
                        VELOCITY_VERLET = 4

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
  complex(r4), allocatable :: projections(:,:,:,:)
  character(len=100) :: filename

  sub_name = 'td_run'; call push_sub()

  write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Elapsed Time '
  call write_info(1)
  
  allocate(dipole(sys%st%nspin, td%save_iter))
  allocate(multipole((td%lmax + 1)**2, sys%st%nspin, td%save_iter))
  if(td%move_ions > 0) then
     allocate(x(td%save_iter, sys%natoms, 3), v(td%save_iter, sys%natoms, 3), &
              f(td%save_iter, sys%natoms, 3))
     allocate(ke(td%save_iter), pe(td%save_iter))
  endif
  if(td%harmonic_spectrum) allocate(tacc(td%save_iter, 3))

  ! occupational analysis stuff
  if(td%occ_analysis) then
    allocate(projections(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik, td%save_iter))
  end if

  ! Calculate initial forces and kinetic energy
  if(td%move_ions > 0) then 
    if(td%iter > 0) then
      call td_read_nbo()
      call generate_external_pot(h, sys)
      sys%eii = ion_ion_energy(sys%natoms, sys%atom)
    end if

    call zforces(h, sys, td%iter*td%dt, td%no_lasers, td%lasers, reduce=.true.)
    sys%kinetic_energy = kinetic_energy(sys%natoms, sys%atom)
    select case(td%move_ions)
      case(NORMAL_VERLET)
        allocate(x1(sys%natoms, 3), x2(sys%natoms, 3))
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
        allocate(f1(sys%natoms, 3))
    end select
  endif

  if(td%iter == 0) call td_run_zero_iter(sys%m)
  td%iter = td%iter + 1

  ii = 1
  etime = elapsed_time()
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
    call td_rti(sys, h, td, i*td%dt)

    ! update density
    call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call zhamiltonian_setup(h, sys)
    call zhamiltonian_eigenval (h, sys, sys%st%st_start, sys%st%st_end) ! eigenvalues
    call hamiltonian_energy(h, sys, -1, reduce=.true.)

    ! Recalculate forces, update velocities...
    if(td%move_ions > 0) then
      pe(ii) = h%etot
      if(td%move_ions == VELOCITY_VERLET) then
        do j = 1, sys%natoms
           f1(j, :) = sys%atom(j)%f(:)
        enddo
      endif
      call zforces(h, sys, i*td%dt, td%no_lasers, td%lasers, reduce=.true.)
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
        x(ii, j, 1:3) = sys%atom(j)%x(1:3)
        v(ii, j, 1:3) = sys%atom(j)%v(1:3)
        f(ii, j, 1:3) = sys%atom(j)%f(1:3)
      end do
    end if

    ! If harmonic spectrum is desired, get the acceleration
    if(td%harmonic_spectrum) call td_calc_tacc(tacc(ii, 1:3), td%dt*i)

    ! measuring
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, &
         dipole(:,ii), multipole(:,:,ii))
    if(td%occ_analysis) then
      call td_calc_projection(projections(:,:,:,ii))
    end if

#ifndef DISABLE_PES
    call PES_doit(td%PESv, sys%m, sys%st, ii, td%dt, td%ab_pot)
#endif

    ! mask function?
    if(td%ab == 2) then
      do ik = 1, sys%st%nik
        do ist = sys%st%st_start, sys%st%st_end
          do idim = 1, sys%st%dim
            sys%st%zpsi(1:sys%m%np, idim, ist, ik) = sys%st%zpsi(1:sys%m%np, idim, ist, ik) * &
                 (1._r8 - td%ab_pot(1:sys%m%np))
          end do
        end do
      end do
    end if

    ! write info
    write(message(1), '(i7,1x,3f14.6)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         elapsed_time() - etime
    call write_info(1)
    etime = elapsed_time()

    ! write down data
    ii = ii + 1
    save: if(ii==td%save_iter+1 .or. i == td%max_iter) then ! output
      if(i == td%max_iter) td%save_iter = ii - 1
      ii = 1

      ! first resume file
      write(filename, '(a,i3.3)') "restart.td.", mpiv%node
      call zstates_write_restart(trim(filename), sys%m, sys%st, &
           iter=i, v1=td%v_old1, v2=td%v_old2)

      if(mpiv%node == 0) call td_write_data()

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
    real(r8) :: x(3), t_acc(3)
    complex(r8) :: c
    real(r8), allocatable :: dipole(:), multipole(:,:), pos(:,:), vel(:,:), for(:,:)

    sub_name = 'td_run_zero_iter'; call push_sub()

    do i = 1, sys%st%nspin
      td%v_old1(:, i) = h%Vhartree(:) + h%Vxc(:, i)
      td%v_old2(:, i) = td%v_old1(:, i)
    end do

    ! we now apply the delta(0) impulse to the wf
    if(td%delta_strength .ne. 0._r8) then
      do i = 1, m%np
        call mesh_xyz(m, i, x)
        c = exp(M_zI * td%delta_strength * sum(x(:)*td%pol(:)))

        sys%st%zpsi(i,:,:,:) = c * sys%st%zpsi(i,:,:,:)
      end do
    end if

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

    ! output laser
    if(td%output_laser) then
      call io_assign(iunit)
      open(unit=iunit, file='td.general/laser', status='unknown')
      call td_write_laser(iunit, 0, 0._r8, .true.)
      call io_close(iunit)
    end if

    ! output occupational analysis data
    if(td%occ_analysis) then
      call io_assign(iunit)
      write(filename, '(a,i3.3)') 'td.general/projections.', mpiv%node
      open(iunit, form='unformatted', status='unknown', file=filename)
      write(iunit) sys%st%nik, sys%st%st_start, sys%st%st_end, u_st%nst
      call io_close(iunit)
    end if

    ! output positions, velocities, forces...
    if(td%move_ions > 0) then
       call io_assign(iunit)
       open(iunit, file='td.general/coordinates')
       allocate(pos(sys%natoms, 3), vel(sys%natoms, 3), for(sys%natoms, 3))
       do j=1, sys%natoms
          pos(j, 1:3) = sys%atom(j)%x(1:3)
          vel(j, 1:3) = sys%atom(j)%v(1:3)
          for(j, 1:3) = sys%atom(j)%f(1:3)
       enddo
       call td_write_nbo(iunit, 0, 0._r8, sys%kinetic_energy, h%etot, pos, vel, for, header = .true.)
       call io_close(iunit)
       deallocate(pos, vel, for)
    endif

    ! output harmonic spectrum
    if(td%harmonic_spectrum) then
       call io_assign(iunit)
       open(iunit, file='td.general/acceleration')
       call td_calc_tacc(t_acc, 0.0_r8)
       call td_write_acc(iunit, 0, 0.0_r8, t_acc, header=.true.)
       call io_close(iunit)
    endif

    call push_sub(); return    
  end subroutine td_run_zero_iter

  subroutine td_write_data()
    integer :: iunit, j, jj, ist, ik, uist

    ! output multipoles
    call io_assign(iunit)
    open(iunit, position='append', file="td.general/multipoles")
    do j = 1, td%save_iter
      jj = i - td%save_iter + j
      call td_write_multipole(iunit, jj, jj*td%dt, &
           dipole(:, j), multipole(:,:, j), .false.)
    end do
    call io_close(iunit)

    ! output the laser field
    if(td%output_laser) then
      call io_assign(iunit)
      open(iunit, position='append', file='td.general/laser')
      do j = 1, td%save_iter
        jj = i - td%save_iter + j
        call td_write_laser(iunit, jj, jj*td%dt, .false.)
      end do
      call io_close(iunit)
    end if

    ! and now we should output the projections
    if(td%occ_analysis) then
      call io_assign(iunit)
      write(filename, '(a,i3.3)') 'td.general/projections.', mpiv%node
      open(iunit, form='unformatted', position='append', file=filename)
      do j = 1, td%save_iter
        do ik = 1, sys%st%nik
          do ist = 1, sys%st%st_start, sys%st%st_end
            write(iunit) (projections(uist, ist, ik, j), uist = 1, u_st%nst)
          end do
        end do
      end do
      call io_close(iunit)
    end if

    ! output positions, vels...
    if(td%move_ions > 0) then
    call io_assign(iunit)
    open(iunit, position='append', file='td.general/coordinates')
    do j = 1, td%save_iter
       jj = i -td%save_iter + j
       call td_write_nbo(iunit, jj, jj*td%dt, ke(j), pe(j), x(j, :, :), v(j, :, :), f(j, :, :), & 
                         header=.false.)
    enddo
    call io_close(iunit)
    endif

    ! output electron acceleration if desired
    if(td%harmonic_spectrum) then
       call io_assign(iunit)
       open(iunit, position='append', file="td.general/acceleration")
       do j = 1, td%save_iter
          jj = i - td%save_iter + j
          call td_write_acc(iunit, jj, jj*td%dt, tacc(j, 1:3), header=.false.)
       enddo
       call io_close(iunit)
    endif

#ifndef DISABLE_PES
    call PES_output(td%PESv, sys%m, sys%st, i, td%save_iter, td%dt)
#endif

  end subroutine td_write_data

  subroutine td_write_laser(iunit, iter, t, header)
    integer, intent(in) :: iunit, iter
    real(r8), intent(IN) :: t
    logical, intent(in) :: header

    integer :: i
    real(r8) :: l_field(1:3), l_vector_field(1:3)
    character(len=80) :: aux

    ! TODO -> confirm these stupid units, especially for the vector field
    if(header) then
      ! first line
      write(iunit, '(a7,e20.12,3a)') '# dt = ', td%dt/units_out%time%factor, &
           "[", trim(units_out%time%abbrev), "]"

      ! second line
      write(iunit, '(a8,a20)', advance='no') '# Iter  ', str_center('t', 20)
      do i = 1, 3
        write(aux, '(a,i1,a)') 'E(', i, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      do i = 1, 3
        write(aux, '(a,i1,a)') 'A(', i, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      write(iunit, '(1x)', advance='yes')

      ! third line
      write(iunit, '(a8,a20)', advance='no') '#       ', &
           str_center('['//trim(units_out%time%abbrev)//']', 20)
      aux = str_center('['//trim(units_out%energy%abbrev) // ' / ' // &
           trim(units_inp%length%abbrev) // ']', 20)
      write(iunit, '(3a20)', advance='no') aux, aux, aux
      aux = str_center('[1/'// trim(units_inp%length%abbrev) // ']', 20)
      write(iunit, '(3a20)', advance='no') aux, aux, aux
      write(iunit, '(1x)', advance='yes')
    end if

    call laser_field(td%no_lasers, td%lasers, t, l_field)
    call laser_vector_field(td%no_lasers, td%lasers, t, l_vector_field)
    write(iunit,'(i8,7es20.12)') iter, t/units_out%time%factor, &
         l_field(1:3) * units_inp%length%factor / units_inp%energy%factor, &
         l_vector_field(1:3) * units_inp%length%factor
    
  end subroutine td_write_laser

  subroutine td_write_acc(iunit, iter, t, acc, header)
    integer, intent(in)  :: iunit, iter
    real(r8), intent(in) :: t, acc(3)
    logical, intent(in)  :: header

    ! first line: column names
    if(header) then
      ! first line -> column names
      write(iunit, '(a8,4a20)') '# Iter  ', str_center('t', 20), str_center('Acc(1)', 20), &
                                            str_center('Acc(2)', 20), str_center('Acc(3)', 20)
      ! second line -> units
      write(iunit, '(a8,2a20)') '#       ',                            &
           str_center('['//trim(units_out%time%abbrev)//']', 20),         &
           str_center('['//trim(units_out%acceleration%abbrev)//']', 20)
    endif

    write(iunit,'(i8,4es20.12)') iter, t/units_out%time%factor, acc/units_out%acceleration%factor

  end subroutine td_write_acc

  subroutine td_write_nbo(iunit, iter, t, ke, pe, x, v, f, header)
    integer, intent(in)  :: iunit, iter
    real(r8), intent(in) :: t
    logical, intent(in)  :: header
    real(r8)             :: ke, pe, x(:, :), v(:, :), f(:, :)

    integer :: i, j
    character(len=50) :: aux

    ! first line: column names
    if(header) then
      write(iunit, '(a8,4a20)', advance='no') '# Iter  ', str_center('t', 20), &
           str_center('Ekin', 20), str_center('Epot', 20), str_center('Etot', 20)
      do i=1, sys%natoms
        do j=1, 3
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',',j,')'
          write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
        end do
      end do
      do i=1, sys%natoms
        do j=1, 3
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
          write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
        enddo
      end do
      do i=1, sys%natoms
        do j=1, 3
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
          write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
        end do
      end do
      write(iunit,'(1x)', advance='yes')
      
      ! second line: units
      write(iunit,'(9a)') '#       ', &
           'Energy in ',       trim(units_out%energy%abbrev),   &
           ', Positions in ',  trim(units_out%length%abbrev),   &
           ', Velocities in ', trim(units_out%velocity%abbrev), &
           ', Forces in ',     trim(units_out%force%abbrev)
    end if

    write(iunit, '(i8, es20.12)', advance='no') iter, t/units_out%time%factor
    write(iunit, '(3es20.12)', advance='no') &
        ke/units_out%energy%factor, &
        pe/units_out%energy%factor, &
        (ke + pe)/units_out%energy%factor
    do i=1, sys%natoms
       write(iunit, '(3es20.12)', advance='no') &
         x(i, 1:3)/units_out%length%factor
    enddo
    do i=1, sys%natoms
       write(iunit, '(3es20.12)', advance='no') &
         v(i, 1:3)/units_out%velocity%factor
    enddo
    do i=1, sys%natoms
       write(iunit, '(3es20.12)', advance='no') & 
         f(i, 1:3)/units_out%force%factor
    enddo
    write(iunit, '(1x)', advance='yes')
    
  end subroutine td_write_nbo

  subroutine td_read_nbo() ! reads the pos and vel from .nbo file
    logical :: found
    integer :: i, iunit

    inquire(file=trim(sys%sysname)//'.nbo', exist=found)
    if(.not.found) then
      message(1) = "Could not open file '"//trim(sys%sysname)//".nbo'"
      message(2) = "Starting simulation from initial geometry"
      call write_warning(2)
      return
    end if
    
    call io_assign(iunit)
    open(iunit, file=trim(sys%sysname)//'.nbo', status='old')

    read(iunit, *); read(iunit, *) ! skip header
    do i = 0, td%iter - 1
      read(iunit, *) ! skip previous iterations... sorry, but no seek in Fortran
    end do
    read(iunit, '(88x)', advance='no') ! skip unrelevant information

    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%x(1:3)
      sys%atom(i)%x(:) = sys%atom(i)%x(:) * units_out%length%factor
    end do
    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%v(1:3)
      sys%atom(i)%v(:) = sys%atom(i)%v(:) * units_out%velocity%factor
    end do
    do i = 1, sys%natoms
      read(iunit, '(3es20.12)', advance='no') sys%atom(i)%f(1:3)
      sys%atom(i)%f(:) = sys%atom(i)%f(:) * units_out%force%factor
    end do

    call io_close(iunit)
  end subroutine td_read_nbo

  subroutine td_write_multipole(iunit, iter, t, dipole, multipole, header)
    integer, intent(in) :: iunit, iter
    real(r8), intent(IN) :: t, dipole(sys%st%nspin), multipole((td%lmax+1)**2, sys%st%nspin)
    logical, intent(in) :: header

    integer :: is, l, m, add_lm
    character(len=50) :: aux
    
    if(header) then
      ! first line
      write(iunit, '(a10,i2,a8,i2)') '# nspin = ', sys%st%nspin, ' lmax = ', td%lmax

      ! second line -> columns name
      write(iunit, '(a8,a20)', advance='no') '# Iter  ', str_center('t', 20)
      do is = 1, sys%st%nspin
        write(aux, '(a,i1,a)') 'dipole(', is, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
            write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')

      ! third line -> units
      write(iunit, '(a8,a20)', advance='no') '#       ', &
           str_center('['//trim(units_out%time%abbrev)//']', 20)
      do is = 1, sys%st%nspin
        write(iunit, '(a20)', advance='no') &
             str_center('['//trim(units_out%length%abbrev)//']', 20)
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            select case(l)
            case(0)
              write(iunit, '(a20)', advance='no') ' '
            case(1)
              write(iunit, '(a20)', advance='no') &
                   str_center('['//trim(units_out%length%abbrev)//']', 20)
            case default
              write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
              write(iunit, '(a20)', advance='no') str_center('['//trim(aux)//']', 20)
            end select
            add_lm = add_lm + 1
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
      
    end if
    
    write(iunit, '(i8, es20.12)', advance='no') iter, t/units_out%time%factor
    do is = 1, sys%st%nspin
      write(iunit, '(es20.12)', advance='no') dipole(is)/units_out%length%factor
    end do
    
    do is = 1, sys%st%nspin
      add_lm = 1
      do l = 0, td%lmax
        do m = -l, l
          write(iunit, '(es20.12)', advance='no') multipole(add_lm, is)/units_out%length%factor**l
          add_lm = add_lm + 1
        end do
      end do
    end do
    write(iunit, '(1x)', advance='yes')
    
  end subroutine td_write_multipole

  subroutine td_calc_projection(p)
    complex(r4), intent(out) :: p(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik)

    integer :: uist, uik, ist, ik

    do ik = 1, sys%st%nik
      do ist = sys%st%st_start, sys%st%st_end
        do uist = 1, u_st%nst
          p(uist, ist, ik) = cmplx(sum(sys%st%zpsi(1:sys%m%np,:, ist, ik)* &
               u_st%R_FUNC(psi) (1:sys%m%np,:, uist, ik)), kind=r4)*sys%m%vol_pp
        end do
      end do
    end do
  end subroutine td_calc_projection

  subroutine td_calc_tacc(acc, t)
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: acc(3)

    real(r8) :: field(3), x(3), y(3), r, &
                vl, dvl, d, charge
    real(r8), allocatable :: V(:), dV(:,:)
    complex(r8) :: uvpsi, p
    integer  :: j, k, is, i, ik, ist, idim, add_lm, l, m, ii, jj
    type(atom_type), pointer :: atm

    sub_name = 'td_calc_tacc'; call push_sub()

    acc(1:3) = 0.0_r8

    ! This should be the ionic acceleration, if it was implemented...
    if(td%move_ions > 0) then
      message(1) = 'Error. If harmonic spectrum is to be calculated, moves should not move'
      message(2) = '(In present version)'
      call write_fatal(2)
      !do j=1, sys%natoms
      !  acc(1:3) = acc(1:3) + sys%atom(j)%spec%z_val*sys%atom(j)%f(1:3)/sys%atom(j)%spec%weight
      !end do
    end if

    ! Gets the gradient of the external pot
    do j = 1, sys%natoms
      do k = 1, sys%m%np
        call mesh_r(sys%m, k, r, x=x, a=sys%atom(j)%x)
        if(r < r_small) cycle
        
        vl  = splint(sys%atom(j)%spec%ps%vlocal, r)
        dvl = splint(sys%atom(j)%spec%ps%dvlocal, r)
        
        d = sum(sys%st%rho(k, :)) * sys%m%vol_pp* &
               (dvl - (vl - sys%atom(j)%spec%Z_val)/r)/r**2
        acc(:) = acc - d * x(:)
      end do
    end do

   ! And now the non-local part...
   ! this comes first to do the reduce...
    x(1:3) = 0.0_r8
    atm_loop: do i = 1, sys%natoms
      atm => sys%atom(i)

      if(atm%spec%local) cycle

      ik_loop: do ik = 1, sys%st%nik
        st_loop: do ist = sys%st%st_start, sys%st%st_end
          dim_loop: do idim = 1, sys%st%dim
            add_lm = 1
            l_loop: do l = 0, atm%spec%ps%L_max
              if(l == atm%spec%ps%L_loc) then
                add_lm = add_lm + (2*l + 1)
                cycle l_loop
              end if

              m_loop: do m = -l, l
                do ii = 1, atm%spec%ps%kbc
                do jj = 1, atm%spec%ps%kbc
                   uVpsi = sum(atm%uV(:, add_lm, ii) * sys%st%occ(ist, ik)  * &
                           sys%st%zpsi(atm%Jxyz(:), idim, ist, ik)) * &
                           sys%m%vol_pp**2 * atm%uVu(add_lm, ii, jj)

                   do j = 1, 3
                      p = sum(atm%duV(j, :, add_lm, jj) * R_CONJ(sys%st%zpsi (atm%Jxyz(:), idim, ist, ik)))
                      x(j) = x(j) + 2._r8 * R_REAL(uvpsi*p)
                   end do
                end do
                end do

                add_lm = add_lm + 1
              end do m_loop
            end do l_loop

          end do dim_loop
        end do st_loop
      end do ik_loop

#if defined(HAVE_MPI) && defined(MPI_TD)
      if(present(reduce) .and. reduce) then
        call MPI_ALLREDUCE(x(1), y(1), 3, &
             MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        x = y
      end if
#endif

    end do atm_loop

    acc(1:3) = acc(1:3) - x(1:3)

    ! Adds the laser contribution
    if(td%no_lasers > 0) then
      call laser_field(td%no_lasers, td%lasers, t, field)
      charge = sum(sys%st%rho(:,:))*sys%m%vol_pp
      acc(1:3) = acc(1:3) - field(1:3)*charge
    end if

    call pop_sub(); return 
  end subroutine td_calc_tacc

end subroutine td_run

#include "td_init.F90"
#include "td_rti.F90"

end module timedep
