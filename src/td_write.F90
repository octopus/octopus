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

module td_write
  use global
  use messages
  use io
  use lib_oct
  use geometry
  use mesh
  use grid
  use mesh_function
  use functions
  use states
  use output
  use restart
  use lasers
  use hamiltonian
  use external_pot
  use grid

  implicit none

  private
  public :: td_write_type, &
            td_write_init, &
            td_write_end, &
            td_write_iter, &
            td_write_data


  type td_write_type
    integer(POINTER_SIZE) :: out_multip, &
                             out_coords, &
                             out_gsp,    &
                             out_acc,    &
                             out_laser,  &
                             out_energy, &
                             out_proj,   &
                             out_angular,&
                             out_spin,   &
                             out_magnets

    integer           :: lmax           ! maximum multipole moment to output
    FLOAT             :: lmm_r          ! radius of the sphere used to compute the local magnetic moments
    type(states_type) :: gs_st          ! The states_type where the ground state is stored, in order to
                                        ! calculate the projections(s) onto it.
  end type td_write_type


contains

subroutine td_write_init(w, gr, st, geo, ions_move, there_are_lasers, iter, dt)
  type(td_write_type)             :: w
  type(grid_type),     intent(in) :: gr
  type(states_type),   intent(in) :: st
  type(geometry_type), intent(in) :: geo
  logical,             intent(in) :: ions_move, there_are_lasers
  integer,             intent(in) :: iter
  FLOAT,               intent(in) :: dt
     

  FLOAT :: rmin
  integer :: ierr, nus, first
  logical :: log
  character(len=256) :: filename

  call push_sub('td_write.td_write_handler')

  call loct_parse_logical(check_inp('TDOutputMultipoles'), .true., log)
       w%out_multip = 0; if(log) w%out_multip = 1
  call loct_parse_logical(check_inp('TDOutputCoordinates'), .true., log)
       if(.not.(ions_move))  log = .false.
       w%out_coords = 0; if(log) w%out_coords = 1
  call loct_parse_logical(check_inp('TDOutputAngularMomentum'), .false., log)
       w%out_angular = 0; if(log) w%out_angular = 1
  call loct_parse_logical(check_inp('TDOutputSpin'), .false., log)
       w%out_spin = 0; if(log) w%out_spin = 1
  call loct_parse_logical(check_inp('TDOutputGSProjection'), .false., log)
       w%out_gsp = 0; if(log) w%out_gsp = 1
  call loct_parse_logical(check_inp('TDOutputAcceleration'), .false., log)
       w%out_acc = 0; if(log) w%out_acc = 1
  call loct_parse_logical(check_inp('TDOutputLaser'), there_are_lasers, log)
       w%out_laser = 0; if(log) w%out_laser = 1
  call loct_parse_logical(check_inp('TDOutputElEnergy'), .false., log)
       w%out_energy = 0; if(log) w%out_energy = 1
  call loct_parse_logical(check_inp('TDOutputOccAnalysis'), .false., log)
       w%out_proj = 0; if(log) w%out_proj = 1
  call loct_parse_logical(check_inp('TDOutputLocalMagneticMoments'), .false., log)
       w%out_magnets = 0; if(log) w%out_magnets = 1

  call loct_parse_int(check_inp('TDDipoleLmax'), 1, w%lmax)
  if (w%lmax < 0 .or. w%lmax > 4) then
    write(message(1), '(a,i6,a)') "Input: '", w%lmax, "' is not a valid TDDipoleLmax"
    message(2) = '(0 <= TDDipoleLmax <= 4 )'
    call write_fatal(2)
  end if

  ! Compatibility test
  if( (w%out_acc.ne.0) .and. ions_move ) then
       message(1) = 'Error. If harmonic spectrum is to be calculated'
       message(2) = 'Atoms should not be allowed to move'
       call write_fatal(2)
  endif

  call geometry_min_distance(geo, rmin)
  call loct_parse_float(check_inp('LocalMagneticMomentsSphereRadius'), rmin*M_HALF/units_inp%length%factor, w%lmm_r)
  w%lmm_r = w%lmm_r * units_inp%length%factor

  if( (w%out_proj.ne.0)  .or.  (w%out_gsp.ne.0) ) then
     call states_copy(w%gs_st, st)
     ! Now include the unoccupied states
     deallocate(w%gs_st%zpsi)
     call loct_parse_int(check_inp('NumberUnoccStates'), 5, nus)
     if(nus < 0) then
       message(1) = "Input: NumberUnoccStates must be >= 0"
       call write_fatal(1)
     end if
     ! fix states: THIS IS NOT OK
     w%gs_st%nst    = w%gs_st%nst + nus
     call states_distribute_nodes(w%gs_st)
     ! allocate memory
     allocate(w%gs_st%zpsi(NP, w%gs_st%d%dim, w%gs_st%st_start:w%gs_st%st_end, w%gs_st%d%nik))
     call zrestart_read(trim(tmpdir)//'restart_gs', w%gs_st, gr%m, ierr)
     if(ierr.ne.0) then
        message(1) = "Could not load "//trim(tmpdir)//"restart_gs"
        call write_fatal(1)
     endif
  endif

  if (iter == 0) then
      first = 0
  else
      first = iter + 1
  end if

  if(mpiv%node==0) then
    if(w%out_multip.ne.0)  call write_iter_init(w%out_multip,  first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/multipoles")))
    if(w%out_angular.ne.0) call write_iter_init(w%out_angular, first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/angular")))
    if(w%out_spin.ne.0)    call write_iter_init(w%out_spin,    first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/spin")))
    if(w%out_magnets.ne.0) call write_iter_init(w%out_magnets, first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/magnetic_moments")))
    if(w%out_coords.ne.0)  call write_iter_init(w%out_coords,  first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/coordinates")))
    if(w%out_gsp.ne.0)     call write_iter_init(w%out_gsp,     first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/gs_projection")))
    if(w%out_acc.ne.0)     call write_iter_init(w%out_acc,     first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/acceleration")))
    if(w%out_laser.ne.0)   call write_iter_init(w%out_laser,   first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/laser")))
    if(w%out_energy.ne.0)  call write_iter_init(w%out_energy,  first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/el_energy")))
    if(w%out_proj.ne.0)    call write_iter_init(w%out_proj,    first, dt/units_out%time%factor, &
              trim(io_workpath("td.general/projections")))
  end if

  call pop_sub()
end subroutine td_write_init

subroutine td_write_end(w)
  type(td_write_type)       :: w
  call push_sub('td_write.td_write_end')

  if(mpiv%node==0) then
    if(w%out_multip.ne.0)  call write_iter_end(w%out_multip)
    if(w%out_angular.ne.0) call write_iter_end(w%out_angular)
    if(w%out_spin.ne.0)    call write_iter_end(w%out_spin)
    if(w%out_magnets.ne.0) call write_iter_end(w%out_magnets)
    if(w%out_coords.ne.0)  call write_iter_end(w%out_coords)
    if(w%out_gsp.ne.0)     call write_iter_end(w%out_gsp)
    if(w%out_acc.ne.0)     call write_iter_end(w%out_acc)
    if(w%out_laser.ne.0)   call write_iter_end(w%out_laser)
    if(w%out_energy.ne.0)  call write_iter_end(w%out_energy)
    if(w%out_proj.ne.0)    call write_iter_end(w%out_proj)
  end if
  call states_end(w%gs_st)

  call pop_sub()
end subroutine td_write_end

subroutine td_write_iter(w, gr, st, h, geo, pol, dt, i)
  type(td_write_type),    intent(in) :: w
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(in)    :: st
  type(hamiltonian_type), intent(in)    :: h
  type(geometry_type),    intent(in)    :: geo
  FLOAT,                  intent(in) :: pol(3)
  FLOAT,                  intent(in) :: dt
  integer,                intent(in) :: i

  call push_sub('td_write.td_write_iter')

  if(w%out_multip.ne.0)   call td_write_multipole(w%out_multip, gr, st, w%lmax, pol, i)
  if(w%out_angular.ne.0)  call td_write_angular(w%out_angular, gr, st, i)
  if(w%out_spin.ne.0)     call td_write_spin(w%out_spin, gr%m, st, i)
  if(w%out_magnets.ne.0)  call td_write_local_magnetic_moments(w%out_magnets, gr%m, st, geo, w%lmm_r, i)
  if(w%out_proj.ne.0)     call td_write_proj(w%out_proj, gr, st, w%gs_st, i)
  if(w%out_coords.ne.0)   call td_write_nbo(w%out_coords, gr, i, geo%kinetic_energy, h%etot)
  if(w%out_acc.ne.0)      call td_write_acc(w%out_acc, gr, st, h, dt, i)
  if(w%out_laser.ne.0)    call td_write_laser(w%out_laser, gr, h, dt, i)
  if(w%out_energy.ne.0)   call td_write_el_energy(w%out_energy, h, i)

  call pop_sub()
end subroutine td_write_iter

subroutine td_write_data(w, gr, st, h, outp, geo, dt, iter)
  type(td_write_type), intent(in) :: w
  type(grid_type),       intent(inout) :: gr
  type(states_type),     intent(in)    :: st
  type(hamiltonian_type), intent(in)    :: h
  type(output_type),     intent(in)    :: outp
  type(geometry_type),   intent(in)    :: geo
  FLOAT, intent(in) :: dt
  integer, intent(in) :: iter

 
  character(len=256) :: filename

  call push_sub('td.td_write_data')

  ! calculate projection onto the ground state
  if(w%out_gsp.ne.0) call td_write_gsp(w%out_gsp, gr%m, st, w%gs_st, dt, iter)

  if(mpiv%node==0) then
      if(w%out_multip.ne.0)  call write_iter_flush(w%out_multip)
      if(w%out_angular.ne.0) call write_iter_flush(w%out_angular)
      if(w%out_spin.ne.0)    call write_iter_flush(w%out_spin)
      if(w%out_magnets.ne.0) call write_iter_flush(w%out_magnets)
      if(w%out_coords.ne.0)  call write_iter_flush(w%out_coords)
      if(w%out_gsp.ne.0)     call write_iter_flush(w%out_gsp)
      if(w%out_acc.ne.0)     call write_iter_flush(w%out_acc)
      if(w%out_laser.ne.0)   call write_iter_flush(w%out_laser)
      if(w%out_energy.ne.0)  call write_iter_flush(w%out_energy)
  end if

  if(w%out_proj.ne.0) call write_iter_flush(w%out_proj)

  ! now write down the rest
  write(filename, '(a,i7.7)') "td.", iter  ! name of directory

  call zstates_output(st, gr, filename, outp)
  if(outp%what(output_geometry)) &
      call atom_write_xyz(filename, "geometry", geo)
  call hamiltonian_output(h, gr%m, gr%sb, filename, outp)

!!$#if !defined(DISABLE_PES) && defined(HAVE_FFT)
!!$  call PES_output(td%PESv, gr%m, st, iter, outp%iter, dt)
!!$#endif

  call pop_sub()
end subroutine td_write_data

subroutine td_write_spin(out_spin, m, st, iter)
  integer(POINTER_SIZE), intent(in) :: out_spin
  type(mesh_type),       intent(in) :: m
  type(states_type),     intent(in) :: st
  integer,               intent(in) :: iter

  character(len=130) :: aux
  FLOAT :: spin(3)

  call push_sub('td_write.td_write_spin')

  if(mpiv%node == 0) then ! only first node outputs

     ! The expectation value of the spin operator is half the total magnetic moment
     call states_magnetic_moment(m, st, st%rho, spin)
     spin = M_HALF*spin

     if(iter ==0) then
        !empty file
        call write_iter_clear(out_spin)

        !fist line -> now unused.
        write(aux, '(a)') '#'
        call write_iter_string(out_spin, aux)
        call write_iter_nl(out_spin)

        !second line -> columns name
        call write_iter_header_start(out_spin)
        if (st%d%ispin == SPINORS) then
           write(aux, '(a2,18x)') 'Sx'
           call write_iter_header(out_spin, aux)
           write(aux, '(a2,18x)') 'Sy'
           call write_iter_header(out_spin, aux)
        end if
        write(aux, '(a2,18x)') 'Sz'
        call write_iter_header(out_spin, aux)
        call write_iter_nl(out_spin)
     endif

     call write_iter_start(out_spin)
     select case (st%d%ispin)
     case (SPIN_POLARIZED)
        call write_iter_double(out_spin, spin(3), 1)
     case (SPINORS)
        call write_iter_double(out_spin, spin(1:3), 3)
     end select
     call write_iter_nl(out_spin)

  end if

  call pop_sub()
end subroutine td_write_spin

subroutine td_write_local_magnetic_moments(out_magnets, m, st, geo, lmm_r, iter)
  integer(POINTER_SIZE), intent(in) :: out_magnets
  type(mesh_type),       intent(in) :: m
  type(states_type),     intent(in) :: st
  type(geometry_type),   intent(in) :: geo
  FLOAT,                 intent(in) :: lmm_r
  integer,               intent(in) :: iter

  integer :: ia
  character(len=50) :: aux
  FLOAT, allocatable :: lmm(:,:)

  call push_sub('td_write.td_write_local_magnetic_moments')

  if(mpiv%node == 0) then ! only first node outputs

    !get the atoms magnetization
    allocate(lmm(3, geo%natoms))
    call states_local_magnetic_moments(m, st, geo, st%rho, lmm_r, lmm)

    if(iter ==0) then
      !empty file
      call write_iter_clear(out_magnets)

      !fist line ->  now unused.
      write(aux, '(a)') '#'
      call write_iter_string(out_magnets, aux)
      call write_iter_nl(out_magnets)

      !second line -> columns name
      call write_iter_header_start(out_magnets)
      do ia = 1, geo%natoms
        if (st%d%ispin == SPINORS) then
          write(aux, '(a2,i2.2,16x)') 'mx', ia
          call write_iter_header(out_magnets, aux)
          write(aux, '(a2,i2.2,16x)') 'my', ia
          call write_iter_header(out_magnets, aux)
        end if
        write(aux, '(a2,i2.2,16x)') 'mz', ia
        call write_iter_header(out_magnets, aux)
      end do
      call write_iter_nl(out_magnets)

    endif

    call write_iter_start(out_magnets)
    do ia = 1, geo%natoms
      select case (st%d%ispin)
      case (SPIN_POLARIZED)
        call write_iter_double(out_magnets, lmm(3, ia), 1)
      case (SPINORS)
        call write_iter_double(out_magnets, lmm(1:3, ia), 3)
      end select
    end do
    call write_iter_nl(out_magnets)
    deallocate(lmm)
  end if

  call pop_sub()
end subroutine td_write_local_magnetic_moments

subroutine td_write_angular(out_angular, gr, st, iter)
  integer(POINTER_SIZE), intent(in)    :: out_angular
  type(grid_type),       intent(inout) :: gr
  type(states_type),     intent(in)    :: st
  integer,               intent(in)    :: iter

  character(len=130) :: aux
  FLOAT :: angular(3)

  call push_sub('td_write.td_write_angular')

  ! The angular momentum has to be calculated by all nodes...
  call zstates_calc_angular(gr, st, angular)

  if(mpiv%node == 0) then ! Only first node outputs

    if(iter ==0) then
      !empty file
      call write_iter_clear(out_angular)

      !fist line -> now unused.
      write(aux, '(a)') '#'
      call write_iter_string(out_angular, aux)
      call write_iter_nl(out_angular)

      !second line -> columns name
      call write_iter_header_start(out_angular)
      write(aux, '(a2,18x)') 'Lx'
      call write_iter_header(out_angular, aux)
      write(aux, '(a2,18x)') 'Ly'
      call write_iter_header(out_angular, aux)
      write(aux, '(a2,18x)') 'Lz'
      call write_iter_header(out_angular, aux)
      call write_iter_nl(out_angular)

      !third line -> should hold the units. Now unused (assumes atomic units)
      call write_iter_string(out_angular, '##########')
      call write_iter_nl(out_angular)
    endif

    call write_iter_start(out_angular)
    call write_iter_double(out_angular, angular(1:3), 3)
    call write_iter_nl(out_angular)

  endif

  call pop_sub()
end subroutine td_write_angular

subroutine td_write_multipole(out_multip, gr, st, lmax, pol, iter)
  integer(POINTER_SIZE), intent(in) :: out_multip
  type(grid_type),       intent(in) :: gr
  type(states_type),     intent(in) :: st
  integer,               intent(in) :: lmax   
  FLOAT,                 intent(in) :: pol(3)
  integer,               intent(in) :: iter

  integer :: is, j, l, m, add_lm
  character(len=50) :: aux
  FLOAT, allocatable :: dipole(:), nuclear_dipole(:), multipole(:,:)

  if(mpiv%node.ne.0) return ! only first node outputs

  call push_sub('td_write.td_write_multipole')

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_multip)

    ! first line
    write(aux, '(a10,i2,a8,i2)') '# nspin = ', st%d%nspin, ' lmax = ', lmax
    call write_iter_string(out_multip, aux)
    call write_iter_nl(out_multip)

    ! second line -> columns name
    call write_iter_header_start(out_multip)

    do is = 1, st%d%nspin
      write(aux, '(a,i1,a)') 'dipole(', is, ')'
      call write_iter_header(out_multip, aux)
    end do
    do j = 1, NDIM
      write(aux,'(a,i1,a)')  'n.dip.(', j,  ')'
      call write_iter_header(out_multip, aux)
    enddo
    do is = 1, st%d%nspin
      do l = 0, lmax
        do m = -l, l
          write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
          call write_iter_header(out_multip, aux)
        end do
      end do
    end do
    call write_iter_nl(out_multip)

    ! third line -> units
    call write_iter_string(out_multip, '##########')
    call write_iter_header(out_multip, '[' // trim(units_out%time%abbrev) // ']')

    do is = 1, st%d%nspin
      call write_iter_header(out_multip, '[' // trim(units_out%length%abbrev) // ']')
    end do
    do j = 1, NDIM
      call write_iter_header(out_multip, '[' // trim(units_out%length%abbrev) // ']')
    enddo

    do is = 1, st%d%nspin
      do l = 0, lmax
        do m = -l, l
          select case(l)
          case(0)
            call write_iter_header(out_multip, '-')
          case(1)
            call write_iter_header(out_multip, '[' // trim(units_out%length%abbrev) // ']')
          case default
            write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
            call write_iter_header(out_multip, '[' // trim(aux) // ']')
          end select
        end do
      end do
    end do
    call write_iter_nl(out_multip)
  end if

  allocate(dipole(st%d%nspin), nuclear_dipole(NDIM), multipole((lmax + 1)**2, st%d%nspin))
  call states_calculate_multipoles(gr, st, pol, dipole, lmax, multipole)
  call geometry_dipole(gr%geo, nuclear_dipole)

  call write_iter_start(out_multip)
  do is = 1, st%d%nspin
    call write_iter_double(out_multip, dipole(is)/units_out%length%factor, 1)
  end do
  do j = 1, NDIM
    call write_iter_double(out_multip, nuclear_dipole(j)/units_out%length%factor, 1)
  enddo
  do is = 1, st%d%nspin
    add_lm = 1
    do l = 0, lmax
      do m = -l, l
        call write_iter_double(out_multip, multipole(add_lm, is)/units_out%length%factor**l, 1)
        add_lm = add_lm + 1
      end do
    end do
  end do
  call write_iter_nl(out_multip)

  deallocate(dipole, nuclear_dipole, multipole)

  call pop_sub()
end subroutine td_write_multipole


subroutine td_write_nbo(out_coords, gr, iter, ke, pe)
  integer(POINTER_SIZE), intent(in) :: out_coords
  type(grid_type),       intent(in) :: gr
  integer,               intent(in) :: iter
  FLOAT,                 intent(in) :: ke, pe

  integer :: i, j
  character(len=50) :: aux

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_coords)

    ! first line: column names
    call write_iter_header_start(out_coords)
    call write_iter_header(out_coords, 'Ekin')
    call write_iter_header(out_coords, 'Epot')
    call write_iter_header(out_coords, 'Etot')

    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
        call write_iter_header(out_coords, aux)
      end do
    end do
    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
        call write_iter_header(out_coords, aux)
      end do
    end do
    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
        call write_iter_header(out_coords, aux)
      end do
    end do
    call write_iter_nl(out_coords)

    ! second line: units
    call write_iter_string(out_coords, '##########')
    call write_iter_header(out_coords, '[' // trim(units_out%time%abbrev) // ']')
    call write_iter_string(out_coords, &
         'Energy in '      // trim(units_out%energy%abbrev)   //   &
         ', Positions in ' // trim(units_out%length%abbrev)   //   &
         ', Velocities in '// trim(units_out%velocity%abbrev) //   &
         ', Forces in '    // trim(units_out%force%abbrev))
    call write_iter_nl(out_coords)
  end if

  call write_iter_start(out_coords)
  call write_iter_double(out_coords,     ke /units_out%energy%factor, 1)
  call write_iter_double(out_coords,     pe /units_out%energy%factor, 1)
  call write_iter_double(out_coords, (ke+pe)/units_out%energy%factor, 1)

  do i = 1, gr%geo%natoms
    call write_iter_double(out_coords, gr%geo%atom(i)%x(1:NDIM)/units_out%length%factor,   NDIM)
  end do
  do i = 1, gr%geo%natoms
    call write_iter_double(out_coords, gr%geo%atom(i)%v(1:NDIM)/units_out%velocity%factor, NDIM)
  end do
  do i = 1, gr%geo%natoms
    call write_iter_double(out_coords, gr%geo%atom(i)%f(1:NDIM)/units_out%force%factor,    NDIM)
  end do
  call write_iter_nl(out_coords)

end subroutine td_write_nbo

subroutine td_write_gsp(out_gsp, m, st, gs_st, dt, iter)
  integer(POINTER_SIZE), intent(in) :: out_gsp
  type(mesh_type),       intent(in) :: m
  type(states_type),     intent(in) :: st
  type(states_type),     intent(in) :: gs_st
  FLOAT,                 intent(in) :: dt
  integer,               intent(in) :: iter

  CMPLX :: gsp
  integer :: ierr

  call push_sub('td_write.td_write_gsp')

  ! all processors calculate the projection
  gsp = zstates_mpdotp(m, 1, gs_st, st)

  if(mpiv%node .eq. 0) then
    if(iter == 0) then
      ! empty file
      call write_iter_clear(out_gsp)

      ! first line -> column names
      call write_iter_header_start(out_gsp)
      call write_iter_header(out_gsp, 'Re <Phi_gs|Phi(t)>')
      call write_iter_header(out_gsp, 'Im <Phi_gs|Phi(t)>')
      call write_iter_nl(out_gsp)

      ! second line -> units
      call write_iter_string(out_gsp, '##########')
      call write_iter_header(out_gsp, '[' // trim(units_out%time%abbrev) // ']')
      call write_iter_nl(out_gsp)
    end if

    ! can not call write_iter_start, for the step is not 1
    call write_iter_int(out_gsp, iter, 1)
    call write_iter_double(out_gsp, iter*dt/units_out%time%factor,  1)
    call write_iter_double(out_gsp, real(gsp),  1)
    call write_iter_double(out_gsp, aimag(gsp), 1)
    call write_iter_nl(out_gsp)
  endif

  call pop_sub()
end subroutine td_write_gsp

subroutine td_write_acc(out_acc, gr, st, h, dt, iter)
  integer(POINTER_SIZE),  intent(in)    :: out_acc
  type(grid_type),        intent(inout) :: gr
  type(states_type),      intent(in)    :: st
  type(hamiltonian_type), intent(in)    :: h
  FLOAT,                  intent(in)    :: dt
  integer,                intent(in)    :: iter

  integer :: i
  character(len=7) :: aux
  FLOAT :: acc(3)

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_acc)

    ! first line -> column names
    call write_iter_header_start(out_acc)
    do i = 1, NDIM
      write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
      call write_iter_header(out_acc, aux)
    end do
    call write_iter_nl(out_acc)

    ! second line: units
    call write_iter_string(out_acc, '##########')
    call write_iter_header(out_acc, '[' // trim(units_out%time%abbrev) // ']')
    do i = 1, NDIM
      call write_iter_header(out_acc, '[' // trim(units_out%acceleration%abbrev) // ']')
    end do
    call write_iter_nl(out_acc)
  endif

  call td_calc_tacc(gr, st, h, acc, dt*i, reduce = .true.)

  call write_iter_start(out_acc)
  call write_iter_double(out_acc, acc/units_out%acceleration%factor, NDIM)
  call write_iter_nl(out_acc)

end subroutine td_write_acc


subroutine td_write_laser(out_laser, gr, h, dt, iter)
  integer(POINTER_SIZE),  intent(in) :: out_laser
  type(grid_type),        intent(in) :: gr
  type(hamiltonian_type), intent(in) :: h
  FLOAT,                  intent(in) :: dt
  integer,                intent(in) :: iter

  integer :: i
  FLOAT :: field(3)
  character(len=80) :: aux

  if(mpiv%node.ne.0) return ! only first node outputs

  ! TODO -> confirm these stupid units, especially for the vector field
  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_laser)

    ! first line
    write(aux, '(a7,e20.12,3a)') '# dt = ', dt/units_out%time%factor, &
         " [", trim(units_out%time%abbrev), "]"
    call write_iter_string(out_laser, aux)
    call write_iter_nl(out_laser)

    ! second line -> column names
    call write_iter_header_start(out_laser)
    do i = 1, NDIM
      write(aux, '(a,i1,a)') 'E(', i, ')'
      call write_iter_header(out_laser, aux)
    end do
    do i = 1, NDIM
      write(aux, '(a,i1,a)') 'A(', i, ')'
      call write_iter_header(out_laser, aux)
    end do
    call write_iter_nl(out_laser)

    ! third line -> units
    call write_iter_string(out_laser, '##########')
    call write_iter_header(out_laser, '[' // trim(units_out%time%abbrev) // ']')

    aux = '[' // trim(units_out%energy%abbrev) // ' / ' // trim(units_inp%length%abbrev) // ']'
    do i = 1, NDIM
      call write_iter_header(out_laser, aux)
    end do

    aux = '[1/' // trim(units_inp%length%abbrev) // ']'
    do i = 1, NDIM
      call write_iter_header(out_laser, aux)
    end do
    call write_iter_nl(out_laser)
  end if

  call write_iter_start(out_laser)

  field = M_ZERO

  call laser_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
  field = field * units_inp%length%factor / units_inp%energy%factor
  call write_iter_double(out_laser, field, NDIM)

  call laser_vector_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
  field = field  * units_inp%length%factor
  call write_iter_double(out_laser, field, NDIM)

  call write_iter_nl(out_laser)

end subroutine td_write_laser


subroutine td_write_el_energy(out_energy, h, iter)
  integer(POINTER_SIZE),  intent(in) :: out_energy
  type(hamiltonian_type), intent(in) :: h
  integer,                intent(in) :: iter

  integer :: i

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_energy)

    ! first line -> column names
    call write_iter_header_start(out_energy)
    call write_iter_header(out_energy, 'Total')
    call write_iter_header(out_energy, 'Ion-Ion')
    call write_iter_header(out_energy, 'Exchange')
    call write_iter_header(out_energy, 'Correlation')
    call write_iter_header(out_energy, 'Potentials')
    call write_iter_nl(out_energy)

    ! second line: units
    call write_iter_string(out_energy, '##########')
    call write_iter_header(out_energy, '[' // trim(units_out%time%abbrev) // ']')
    do i = 1, 5
      call write_iter_header(out_energy, '[' // trim(units_out%energy%abbrev) // ']')
    end do
    call write_iter_nl(out_energy)
  endif

  call write_iter_start(out_energy)
  call write_iter_double(out_energy, h%etot/units_out%energy%factor, 1)
  call write_iter_double(out_energy, h%eii /units_out%energy%factor, 1)
  call write_iter_double(out_energy, h%ex  /units_out%energy%factor, 1)
  call write_iter_double(out_energy, h%ec  /units_out%energy%factor, 1)
  call write_iter_double(out_energy, h%epot/units_out%energy%factor, 1)
  call write_iter_nl(out_energy)


end subroutine td_write_el_energy

subroutine td_write_proj(out_proj, gr, st, u_st, iter)
  integer(POINTER_SIZE), intent(in) :: out_proj
  type(grid_type),       intent(in) :: gr
  type(states_type),     intent(in) :: st
  type(states_type),     intent(in) :: u_st
  integer,               intent(in) :: iter

  CMPLX, allocatable :: projections(:,:,:)
  character(len=20) :: aux
  integer :: ik, ist, uist

  call push_sub('td_write.td_write_proj')

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out_proj)

    ! first line -> column names
    call write_iter_header_start(out_proj)
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do uist = 1, u_st%nst
          write(aux, '(i3,a,i3)') ist, ' -> ', uist
          call write_iter_header(out_proj, 'Re {'//trim(aux)//'}')
          call write_iter_header(out_proj, 'Im {'//trim(aux)//'}')
        end do
      end do
    end do
    call write_iter_nl(out_proj)
  endif

  allocate(projections(u_st%nst, st%st_start:st%st_end, st%d%nik))
  call calc_projection(gr%m, u_st, st, projections)

  call write_iter_start(out_proj)
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do uist = 1, u_st%nst
        call write_iter_double(out_proj,  real(projections(uist, ist, ik)), 1)
        call write_iter_double(out_proj, aimag(projections(uist, ist, ik)), 1)
      end do
    end do
  end do
  call write_iter_nl(out_proj)

  deallocate(projections)
  call pop_sub()
end subroutine td_write_proj

#include "td_calc.F90"

end module td_write
