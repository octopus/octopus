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

subroutine td_write_spin(out, m, st, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(mesh_type),       intent(in) :: m
  type(states_type),     intent(in) :: st
  type(td_type),         intent(in) :: td
  integer,               intent(in) :: iter

  character(len=130) :: aux
  FLOAT :: spin(3)

  call push_sub('td_write_spin')

  ! The spin has to be calculated by all nodes...
  ! The expectation value of the spin operator is half the magnetization value
  call zstates_calculate_magnetization(m, st, spin)
  spin = M_HALF*spin

  if(mpiv%node.ne.0) then ! only first node outputs
    call pop_sub(); return
  end if

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

  call pop_sub()
end subroutine td_write_spin

subroutine td_write_angular(out, m, f_der, st, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(mesh_type),       intent(IN) :: m
  type(f_der_type),      intent(inout) :: f_der
  type(states_type),     intent(IN) :: st
  type(td_type),         intent(IN) :: td
  integer,               intent(in) :: iter

  integer :: ierr
  character(len=130) :: aux
  FLOAT :: angular(3)

  call push_sub('td_write_angular')

  ! The angular momentum has to be calculated by all nodes...
  call zstates_calculate_angular(m, f_der, st, angular)

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

subroutine td_write_multipole(out, mesh, st, geo, td, iter)
  integer(POINTER_SIZE), intent(IN) :: out
  type(mesh_type),       intent(IN) :: mesh
  type(states_type),     intent(IN) :: st
  type(geometry_type),   intent(IN) :: geo
  type(td_type),         intent(IN) :: td
  integer,               intent(in) :: iter
  
  integer :: is, j, l, m, add_lm
  character(len=50) :: aux
  FLOAT, allocatable :: dipole(:), nuclear_dipole(:), multipole(:,:)

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)
    
    ! first line
    write(aux, '(a10,i2,a8,i2)') '# nspin = ', st%d%nspin, ' lmax = ', td%lmax
    call write_iter_string(out, aux)
    call write_iter_nl(out)
    
    ! second line -> columns name
    call write_iter_header_start(out)

    do is = 1, st%d%nspin
      write(aux, '(a,i1,a)') 'dipole(', is, ')'
      call write_iter_header(out, aux)
    end do
    do j = 1, conf%dim
      write(aux,'(a,i1,a)')  'n.dip.(', j,  ')'
      call write_iter_header(out, aux)
    enddo
    do is = 1, st%d%nspin
      do l = 0, td%lmax
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
    do j = 1, conf%dim
      call write_iter_header(out, '[' // trim(units_out%length%abbrev) // ']')
    enddo

    do is = 1, st%d%nspin
      do l = 0, td%lmax
        do m = -l, l
          select case(l)
          case(0)
            call write_iter_header(out, ' ')
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
    
  allocate(dipole(st%d%nspin), nuclear_dipole(conf%dim), multipole((td%lmax + 1)**2, st%d%nspin))
  call states_calculate_multipoles(mesh, st, td%pol, dipole, td%lmax, multipole)
  call geometry_dipole(geo, nuclear_dipole)
  
  call write_iter_start(out)
  do is = 1, st%d%nspin
    call write_iter_double(out, dipole(is)/units_out%length%factor, 1)
  end do
  do j = 1, conf%dim
    call write_iter_double(out, nuclear_dipole(j)/units_out%length%factor, 1)
  enddo
  do is = 1, st%d%nspin
    add_lm = 1
    do l = 0, td%lmax
      do m = -l, l
        call write_iter_double(out, multipole(add_lm, is)/units_out%length%factor**l, 1)
        add_lm = add_lm + 1
      end do
    end do
  end do
  call write_iter_nl(out)
  
  deallocate(dipole, nuclear_dipole, multipole)

end subroutine td_write_multipole

subroutine td_write_nbo(out, geo, td, iter, ke, pe)
  integer(POINTER_SIZE), intent(IN) :: out
  type(geometry_type),   intent(IN) :: geo
  type(td_type),         intent(IN) :: td
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
    
    do i = 1, geo%natoms
      do j = 1, conf%dim
        write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, geo%natoms
      do j = 1, conf%dim
        write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, geo%natoms
      do j = 1, conf%dim
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

  do i = 1, geo%natoms
    call write_iter_double(out, geo%atom(i)%x(1:conf%dim)/units_out%length%factor,   conf%dim)
  end do
  do i = 1, geo%natoms
    call write_iter_double(out, geo%atom(i)%v(1:conf%dim)/units_out%velocity%factor, conf%dim)
  end do
  do i = 1, geo%natoms
    call write_iter_double(out, geo%atom(i)%f(1:conf%dim)/units_out%force%factor,    conf%dim)
  end do
  call write_iter_nl(out)
  
end subroutine td_write_nbo

subroutine td_write_gsp(out, m, st, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(mesh_type),       intent(IN) :: m
  type(states_type),     intent(IN) :: st
  type(td_type),         intent(IN) :: td
  integer,               intent(in) :: iter

  CMPLX :: gsp
  type(states_type) :: stgs
  integer :: ierr

  call push_sub('td_write_gsp')

  ! all processors calculate the projection
  stgs = st
  call X(restart_read)("tmp/restart_gs", stgs, m, ierr)
  if(ierr.ne.0) then
    message(1) = 'Error loading GS in zstates_project_gs'
    call write_fatal(1)
  endif
  gsp = zstates_mpdotp(m, 1, stgs, st)
  call states_end(stgs)

  ! but only first node outputs
  if(mpiv%node.ne.0) then
    call pop_sub; return
  end if
  
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
  call write_iter_double(out, iter*td%dt/units_out%time%factor,  1)
  call write_iter_double(out, real(gsp),  1)
  call write_iter_double(out, aimag(gsp), 1)
  call write_iter_nl(out)
  
  call pop_sub()
end subroutine td_write_gsp

subroutine td_write_acc(out, mesh, f_der, st, geo, h, td, iter)
  integer(POINTER_SIZE),  intent(IN)    :: out
  type(mesh_type),        intent(IN)    :: mesh
  type(f_der_type),       intent(inout) :: f_der
  type(states_type),      intent(inout) :: st
  type(geometry_type),    intent(inout) :: geo
  type(hamiltonian_type), intent(IN)    :: h
  type(td_type),          intent(IN)    :: td
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
    do i = 1, conf%dim
      write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)
    
    ! second line: units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
    do i = 1, conf%dim
      call write_iter_header(out, '[' // trim(units_out%acceleration%abbrev) // ']')
    end do
    call write_iter_nl(out)
  endif
  
  call td_calc_tacc(acc, td%dt*i, reduce = .true.)
  
  call write_iter_start(out)
  call write_iter_double(out, acc/units_out%acceleration%factor, conf%dim)
  call write_iter_nl(out)
  
contains

#include "td_calc.F90"

end subroutine td_write_acc

subroutine td_write_laser(out, h, td, iter)
  integer(POINTER_SIZE),  intent(in) :: out
  type(hamiltonian_type), intent(IN) :: h
  type(td_type),          intent(IN) :: td
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
    write(aux, '(a7,e20.12,3a)') '# dt = ', td%dt/units_out%time%factor, &
         " [", trim(units_out%time%abbrev), "]"
    call write_iter_string(out, aux)
    call write_iter_nl(out)
    
    ! second line -> column names
    call write_iter_header_start(out)
    do i = 1, conf%dim
      write(aux, '(a,i1,a)') 'E(', i, ')'
      call write_iter_header(out, aux)
    end do
    do i = 1, conf%dim
      write(aux, '(a,i1,a)') 'A(', i, ')'
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)
    
    ! third line -> units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '[' // trim(units_out%time%abbrev) // ']')
    
    aux = '[' // trim(units_out%energy%abbrev) // ' / ' // trim(units_inp%length%abbrev) // ']'
    do i = 1, conf%dim
      call write_iter_header(out, aux)
    end do
    
    aux = '[1/' // trim(units_inp%length%abbrev) // ']'
    do i = 1, conf%dim
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)
  end if
  
  call write_iter_start(out)
  
  field = M_ZERO
  call epot_laser_field(h%ep, iter*td%dt, field)
  field = field * units_inp%length%factor / units_inp%energy%factor
  call write_iter_double(out, field, conf%dim)
  
  call epot_laser_vector_field(h%ep, iter*td%dt, field)
  field = field  * units_inp%length%factor
  call write_iter_double(out, field, conf%dim)
  
  call write_iter_nl(out)
  
end subroutine td_write_laser
    
subroutine td_write_el_energy(out, h, iter)
  integer(POINTER_SIZE),  intent(in) :: out
  type(hamiltonian_type), intent(IN) :: h
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

subroutine td_write_proj(out, m, st, u_st, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(mesh_type),       intent(IN) :: m
  type(states_type),     intent(IN) :: st
  type(states_type),     intent(IN) :: u_st
  integer,               intent(in) :: iter

  CMPLX, allocatable :: projections(:,:,:)
  character(len=20) :: aux
  integer :: ik, ist, uist

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)
    
    ! first line -> column names
    call write_iter_header_start(out)
    do ik = 1, st%nik
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

  allocate(projections(u_st%nst, st%st_start:st%st_end, st%nik))
  call calc_projection(u_st, st, m, projections)

  call write_iter_start(out)
  do ik = 1, st%nik
    do ist = 1, st%nst
      do uist = 1, u_st%nst
        call write_iter_double(out,  real(projections(uist, ist, ik)), 1)
        call write_iter_double(out, aimag(projections(uist, ist, ik)), 1)
      end do
    end do
  end do
  call write_iter_nl(out)

  deallocate(projections)
end subroutine td_write_proj
