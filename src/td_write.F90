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
  public :: td_write_spin, &
            td_write_local_magnetic_moments, &
            td_write_angular, &
            td_write_multipole, &
            td_write_nbo, &
            td_write_el_energy, &
            td_write_laser, &
            td_write_acc, &
            td_write_gsp, &
            td_write_proj


contains

subroutine td_write_spin(out, m, st, iter)
  integer(POINTER_SIZE), intent(in) :: out
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
        call write_iter_clear(out)

        !fist line -> now unused.
        write(aux, '(a)') '#'
        call write_iter_string(out, aux)
        call write_iter_nl(out)

        !second line -> columns name
        call write_iter_header_start(out)
        if (st%d%ispin == SPINORS) then
           write(aux, '(a2,18x)') 'Sx'
           call write_iter_header(out, aux)
           write(aux, '(a2,18x)') 'Sy'
           call write_iter_header(out, aux)
        end if
        write(aux, '(a2,18x)') 'Sz'
        call write_iter_header(out, aux)
        call write_iter_nl(out)
     endif

     call write_iter_start(out)
     select case (st%d%ispin)
     case (SPIN_POLARIZED)
        call write_iter_double(out, spin(3), 1)
     case (SPINORS)
        call write_iter_double(out, spin(1:3), 3)
     end select
     call write_iter_nl(out)

  end if

  call pop_sub()
end subroutine td_write_spin

subroutine td_write_local_magnetic_moments(out, m, st, geo, lmm_r, iter)
  integer(POINTER_SIZE), intent(in) :: out
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
      call write_iter_clear(out)

      !fist line ->  now unused.
      write(aux, '(a)') '#'
      call write_iter_string(out, aux)
      call write_iter_nl(out)

      !second line -> columns name
      call write_iter_header_start(out)
      do ia = 1, geo%natoms
        if (st%d%ispin == SPINORS) then
          write(aux, '(a2,i2.2,16x)') 'mx', ia
          call write_iter_header(out, aux)
          write(aux, '(a2,i2.2,16x)') 'my', ia
          call write_iter_header(out, aux)
        end if
        write(aux, '(a2,i2.2,16x)') 'mz', ia
        call write_iter_header(out, aux)
      end do
      call write_iter_nl(out)

    endif

    call write_iter_start(out)
    do ia = 1, geo%natoms
      select case (st%d%ispin)
      case (SPIN_POLARIZED)
        call write_iter_double(out, lmm(3, ia), 1)
      case (SPINORS)
        call write_iter_double(out, lmm(1:3, ia), 3)
      end select
    end do
    call write_iter_nl(out)
    deallocate(lmm)
  end if

  call pop_sub()
end subroutine td_write_local_magnetic_moments

subroutine td_write_angular(out, gr, st, iter)
  integer(POINTER_SIZE), intent(in)    :: out
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
      call write_iter_clear(out)

      !fist line -> now unused.
      write(aux, '(a)') '#'
      call write_iter_string(out, aux)
      call write_iter_nl(out)

      !second line -> columns name
      call write_iter_header_start(out)
      write(aux, '(a2,18x)') 'Lx'
      call write_iter_header(out, aux)
      write(aux, '(a2,18x)') 'Ly'
      call write_iter_header(out, aux)
      write(aux, '(a2,18x)') 'Lz'
      call write_iter_header(out, aux)
      call write_iter_nl(out)

      !third line -> should hold the units. Now unused (assumes atomic units)
      call write_iter_string(out, '##########')
      call write_iter_nl(out)
    endif

    call write_iter_start(out)
    call write_iter_double(out, angular(1:3), 3)
    call write_iter_nl(out)

  endif

  call pop_sub()
end subroutine td_write_angular

subroutine td_write_multipole(gr, out, st, lmax, pol, iter)
  type(grid_type),       intent(in) :: gr
  integer(POINTER_SIZE), intent(in) :: out
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
    call write_iter_clear(out)

    ! first line
    write(aux, '(a10,i2,a8,i2)') '# nspin = ', st%d%nspin, ' lmax = ', lmax
    call write_iter_string(out, aux)
    call write_iter_nl(out)

    ! second line -> columns name
    call write_iter_header_start(out)

    do is = 1, st%d%nspin
      write(aux, '(a,i1,a)') 'dipole(', is, ')'
      call write_iter_header(out, aux)
    end do
    do j = 1, NDIM
      write(aux,'(a,i1,a)')  'n.dip.(', j,  ')'
      call write_iter_header(out, aux)
    enddo
    do is = 1, st%d%nspin
      do l = 0, lmax
        do m = -l, l
          write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
          call write_iter_header(out, aux)
        end do
      end do
    end do
    call write_iter_nl(out)

    ! third line -> units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')

    do is = 1, st%d%nspin
      call write_iter_header(out, '[' // trim(units_out%length%abbrev) // ']')
    end do
    do j = 1, NDIM
      call write_iter_header(out, '[' // trim(units_out%length%abbrev) // ']')
    enddo

    do is = 1, st%d%nspin
      do l = 0, lmax
        do m = -l, l
          select case(l)
          case(0)
            call write_iter_header(out, '-')
          case(1)
            call write_iter_header(out, '[' // trim(units_out%length%abbrev) // ']')
          case default
            write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
            call write_iter_header(out, '[' // trim(aux) // ']')
          end select
        end do
      end do
    end do
    call write_iter_nl(out)
  end if

  allocate(dipole(st%d%nspin), nuclear_dipole(NDIM), multipole((lmax + 1)**2, st%d%nspin))
  call states_calculate_multipoles(gr, st, pol, dipole, lmax, multipole)
  call geometry_dipole(gr%geo, nuclear_dipole)

  call write_iter_start(out)
  do is = 1, st%d%nspin
    call write_iter_double(out, dipole(is)/units_out%length%factor, 1)
  end do
  do j = 1, NDIM
    call write_iter_double(out, nuclear_dipole(j)/units_out%length%factor, 1)
  enddo
  do is = 1, st%d%nspin
    add_lm = 1
    do l = 0, lmax
      do m = -l, l
        call write_iter_double(out, multipole(add_lm, is)/units_out%length%factor**l, 1)
        add_lm = add_lm + 1
      end do
    end do
  end do
  call write_iter_nl(out)

  deallocate(dipole, nuclear_dipole, multipole)

  call pop_sub()
end subroutine td_write_multipole


subroutine td_write_nbo(gr, out, iter, ke, pe)
  type(grid_type),       intent(in) :: gr
  integer(POINTER_SIZE), intent(in) :: out
  integer,               intent(in) :: iter
  FLOAT,                 intent(in) :: ke, pe

  integer :: i, j
  character(len=50) :: aux

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)

    ! first line: column names
    call write_iter_header_start(out)
    call write_iter_header(out, 'Ekin')
    call write_iter_header(out, 'Epot')
    call write_iter_header(out, 'Etot')

    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, gr%geo%natoms
      do j = 1, NDIM
        write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
        call write_iter_header(out, aux)
      end do
    end do
    call write_iter_nl(out)

    ! second line: units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
    call write_iter_string(out, &
         'Energy in '      // trim(units_out%energy%abbrev)   //   &
         ', Positions in ' // trim(units_out%length%abbrev)   //   &
         ', Velocities in '// trim(units_out%velocity%abbrev) //   &
         ', Forces in '    // trim(units_out%force%abbrev))
    call write_iter_nl(out)
  end if

  call write_iter_start(out)
  call write_iter_double(out,     ke /units_out%energy%factor, 1)
  call write_iter_double(out,     pe /units_out%energy%factor, 1)
  call write_iter_double(out, (ke+pe)/units_out%energy%factor, 1)

  do i = 1, gr%geo%natoms
    call write_iter_double(out, gr%geo%atom(i)%x(1:NDIM)/units_out%length%factor,   NDIM)
  end do
  do i = 1, gr%geo%natoms
    call write_iter_double(out, gr%geo%atom(i)%v(1:NDIM)/units_out%velocity%factor, NDIM)
  end do
  do i = 1, gr%geo%natoms
    call write_iter_double(out, gr%geo%atom(i)%f(1:NDIM)/units_out%force%factor,    NDIM)
  end do
  call write_iter_nl(out)

end subroutine td_write_nbo

subroutine td_write_gsp(out, m, st, dt, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(mesh_type),       intent(in) :: m
  type(states_type),     intent(in) :: st
  FLOAT,                 intent(in) :: dt
  integer,               intent(in) :: iter

  CMPLX :: gsp
  type(states_type) :: stgs
  integer :: ierr

  call push_sub('td_write.td_write_gsp')

  ! all processors calculate the projection
  call states_copy(stgs, st)
  call zrestart_read (trim(tmpdir)//'restart_gs', stgs, m, ierr)
  if(ierr.ne.0) then
    message(1) = 'Error loading GS in zstates_project_gs'
    call write_fatal(1)
  endif
  gsp = zstates_mpdotp(m, 1, stgs, st)
  call states_end(stgs)

  if(mpiv%node .eq. 0) then
    if(iter == 0) then
      ! empty file
      call write_iter_clear(out)

      ! first line -> column names
      call write_iter_header_start(out)
      call write_iter_header(out, 'Re <Phi_gs|Phi(t)>')
      call write_iter_header(out, 'Im <Phi_gs|Phi(t)>')
      call write_iter_nl(out)

      ! second line -> units
      call write_iter_string(out, '##########')
      call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
      call write_iter_nl(out)
    end if

    ! can not call write_iter_start, for the step is not 1
    call write_iter_int(out, iter, 1)
    call write_iter_double(out, iter*dt/units_out%time%factor,  1)
    call write_iter_double(out, real(gsp),  1)
    call write_iter_double(out, aimag(gsp), 1)
    call write_iter_nl(out)
  endif

  call pop_sub()
end subroutine td_write_gsp

subroutine td_write_acc(gr, out, st, h, dt, iter)
  type(grid_type),        intent(inout) :: gr
  integer(POINTER_SIZE),  intent(in)    :: out
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
    call write_iter_clear(out)

    ! first line -> column names
    call write_iter_header_start(out)
    do i = 1, NDIM
      write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)

    ! second line: units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
    do i = 1, NDIM
      call write_iter_header(out, '[' // trim(units_out%acceleration%abbrev) // ']')
    end do
    call write_iter_nl(out)
  endif

  call td_calc_tacc(gr, st, h, acc, dt*i, reduce = .true.)

  call write_iter_start(out)
  call write_iter_double(out, acc/units_out%acceleration%factor, NDIM)
  call write_iter_nl(out)

end subroutine td_write_acc


subroutine td_write_laser(gr, out, h, dt, iter)
  type(grid_type),        intent(in) :: gr
  integer(POINTER_SIZE),  intent(in) :: out
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
    call write_iter_clear(out)

    ! first line
    write(aux, '(a7,e20.12,3a)') '# dt = ', dt/units_out%time%factor, &
         " [", trim(units_out%time%abbrev), "]"
    call write_iter_string(out, aux)
    call write_iter_nl(out)

    ! second line -> column names
    call write_iter_header_start(out)
    do i = 1, NDIM
      write(aux, '(a,i1,a)') 'E(', i, ')'
      call write_iter_header(out, aux)
    end do
    do i = 1, NDIM
      write(aux, '(a,i1,a)') 'A(', i, ')'
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)

    ! third line -> units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')

    aux = '[' // trim(units_out%energy%abbrev) // ' / ' // trim(units_inp%length%abbrev) // ']'
    do i = 1, NDIM
      call write_iter_header(out, aux)
    end do

    aux = '[1/' // trim(units_inp%length%abbrev) // ']'
    do i = 1, NDIM
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)
  end if

  call write_iter_start(out)

  field = M_ZERO

  call laser_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
  field = field * units_inp%length%factor / units_inp%energy%factor
  call write_iter_double(out, field, NDIM)

  call laser_vector_field(gr%sb, h%ep%no_lasers, h%ep%lasers, iter*dt, field)
  field = field  * units_inp%length%factor
  call write_iter_double(out, field, NDIM)

  call write_iter_nl(out)

end subroutine td_write_laser


subroutine td_write_el_energy(out, h, iter)
  integer(POINTER_SIZE),  intent(in) :: out
  type(hamiltonian_type), intent(in) :: h
  integer,                intent(in) :: iter

  integer :: i

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)

    ! first line -> column names
    call write_iter_header_start(out)
    call write_iter_header(out, 'Total')
    call write_iter_header(out, 'Ion-Ion')
    call write_iter_header(out, 'Exchange')
    call write_iter_header(out, 'Correlation')
    call write_iter_header(out, 'Potentials')
    call write_iter_nl(out)

    ! second line: units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
    do i = 1, 5
      call write_iter_header(out, '[' // trim(units_out%energy%abbrev) // ']')
    end do
    call write_iter_nl(out)
  endif

  call write_iter_start(out)
  call write_iter_double(out, h%etot/units_out%energy%factor, 1)
  call write_iter_double(out, h%eii /units_out%energy%factor, 1)
  call write_iter_double(out, h%ex  /units_out%energy%factor, 1)
  call write_iter_double(out, h%ec  /units_out%energy%factor, 1)
  call write_iter_double(out, h%epot/units_out%energy%factor, 1)
  call write_iter_nl(out)

end subroutine td_write_el_energy

subroutine td_write_proj(gr, out, st, u_st, iter)
  type(grid_type),       intent(in) :: gr
  integer(POINTER_SIZE), intent(in) :: out
  type(states_type),     intent(in) :: st
  type(states_type),     intent(in) :: u_st
  integer,               intent(in) :: iter

  CMPLX, allocatable :: projections(:,:,:)
  character(len=20) :: aux
  integer :: ik, ist, uist

  call push_sub('td_write.td_write_proj')

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)

    ! first line -> column names
    call write_iter_header_start(out)
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do uist = 1, u_st%nst
          write(aux, '(i3,a,i3)') ist, ' -> ', uist
          call write_iter_header(out, 'Re {'//trim(aux)//'}')
          call write_iter_header(out, 'Im {'//trim(aux)//'}')
        end do
      end do
    end do
    call write_iter_nl(out)
  endif

  allocate(projections(u_st%nst, st%st_start:st%st_end, st%d%nik))
  call calc_projection(gr%m, u_st, st, projections)

  call write_iter_start(out)
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do uist = 1, u_st%nst
        call write_iter_double(out,  real(projections(uist, ist, ik)), 1)
        call write_iter_double(out, aimag(projections(uist, ist, ik)), 1)
      end do
    end do
  end do
  call write_iter_nl(out)

  deallocate(projections)
  call pop_sub()
end subroutine td_write_proj

#include "td_calc.F90"

end module td_write
