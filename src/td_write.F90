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

subroutine td_write_angular(out, sys, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(system_type), intent(in) :: sys
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter

  integer :: is, i
  character(len=130) :: aux
  real(r8) :: angular(3)

  if(mpiv%node.ne.0) return ! only first node outputs

  if(iter ==0) then
    !empty file
    call write_iter_clear(out)

    !first line -> columns name
    call write_iter_header_start(out)
    write(aux, '(a2,18x)') 'Lx'
    call write_iter_header(out, aux)
    write(aux, '(a2,18x)') 'Ly'
    call write_iter_header(out, aux)
    write(aux, '(a2,18x)') 'Lz'
    call write_iter_header(out, aux)
    call write_iter_nl(out)
  endif

  call zstates_calculate_angular(sys%m, sys%st, angular)

  call write_iter_start(out)
  call write_iter_double(out, angular(1:3), 3)
  call write_iter_nl(out)
  
end subroutine td_write_angular

subroutine td_write_multipole(out, sys, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(system_type), intent(in) :: sys
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter
  
  integer :: is, l, m, add_lm
  character(len=50) :: aux
  real(r8), allocatable :: dipole(:), multipole(:,:)
  
  if(mpiv%node.ne.0) return ! only first node outputs
  
  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)
    
    ! first line
    write(aux, '(a10,i2,a8,i2)') '# nspin = ', sys%st%nspin, ' lmax = ', td%lmax
    call write_iter_string(out, aux)
    call write_iter_nl(out)
    
    ! second line -> columns name
    call write_iter_header_start(out)

    do is = 1, sys%st%nspin
      write(aux, '(a,i1,a)') 'dipole(', is, ')'
      call write_iter_header(out, aux)
    end do
    do is = 1, sys%st%nspin
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
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
    
    do is = 1, sys%st%nspin
      call write_iter_header(out, '['+trim(units_out%length%abbrev)+']')
    end do
    
    do is = 1, sys%st%nspin
      do l = 0, td%lmax
        do m = -l, l
          select case(l)
          case(0)
            call write_iter_header(out, ' ')
          case(1)
            call write_iter_header(out, '['+trim(units_out%length%abbrev)+']')
          case default
            write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
            call write_iter_header(out, '['+trim(aux)+']')
          end select
        end do
      end do
    end do
    call write_iter_nl(out)
  end if
    
  allocate(dipole(sys%st%nspin), multipole((td%lmax + 1)**2, sys%st%nspin))
  call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, dipole, multipole)
  
  call write_iter_start(out)
  do is = 1, sys%st%nspin
    call write_iter_double(out, dipole(is)/units_out%length%factor, 1)
  end do
  do is = 1, sys%st%nspin
    add_lm = 1
    do l = 0, td%lmax
      do m = -l, l
        call write_iter_double(out, multipole(add_lm, is)/units_out%length%factor**l, 1)
        add_lm = add_lm + 1
      end do
    end do
  end do
  call write_iter_nl(out)
  
  deallocate(dipole, multipole)
end subroutine td_write_multipole

subroutine td_write_nbo(out, sys, td, iter, ke, pe)
  integer(POINTER_SIZE), intent(in) :: out
  type(system_type), intent(in) :: sys
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter
  real(r8),          intent(in) :: ke, pe

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
    
    do i = 1, sys%natoms
      do j = 1, conf%dim
        write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, sys%natoms
      do j = 1, conf%dim
        write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
        call write_iter_header(out, aux)
      end do
    end do
    do i = 1, sys%natoms
      do j = 1, conf%dim
        write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
        call write_iter_header(out, aux)
      end do
    end do
    call write_iter_nl(out)
    
    ! second line: units
    call write_iter_string(out, '##########')
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
    call write_iter_string(out, &
         'Energy in '      + trim(units_out%energy%abbrev)   +   &
         ', Positions in ' + trim(units_out%length%abbrev)   +   &
         ', Velocities in '+ trim(units_out%velocity%abbrev) + &
         ', Forces in '    + trim(units_out%force%abbrev))
    call write_iter_nl(out)
  end if
  
  call write_iter_start(out)
  call write_iter_double(out,     ke /units_out%energy%factor, 1)
  call write_iter_double(out,     pe /units_out%energy%factor, 1)
  call write_iter_double(out, (ke+pe)/units_out%energy%factor, 1)

  do i = 1, sys%natoms
    call write_iter_double(out, sys%atom(i)%x(1:conf%dim)/units_out%length%factor,   conf%dim)
  end do
  do i = 1, sys%natoms
    call write_iter_double(out, sys%atom(i)%v(1:conf%dim)/units_out%velocity%factor, conf%dim)
  end do
  do i = 1, sys%natoms
    call write_iter_double(out, sys%atom(i)%f(1:conf%dim)/units_out%force%factor,    conf%dim)
  end do
  call write_iter_nl(out)
  
end subroutine td_write_nbo

subroutine td_write_gsp(out, sys, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(system_type), intent(in) :: sys
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter

  complex(r8) :: gsp
  
  call push_sub('td_write_gsp')
  
  ! all processors calculate the projection
  call zstates_project_gs(sys%st, sys%m, gsp)
  
  ! but only first node outputs
  if(mpiv%node.ne.0) return
  
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
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
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

subroutine td_write_acc(out, sys, h, td, iter)
  integer(POINTER_SIZE),  intent(in) :: out
  type(system_type),      intent(inout) :: sys
  type(hamiltonian_type), intent(in) :: h
  type(td_type),          intent(in) :: td
  integer,                intent(in) :: iter

  integer :: i
  character(len=7) :: aux
  real(r8) :: acc(3)
  
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
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
    do i = 1, conf%dim
      call write_iter_header(out, '['+trim(units_out%acceleration%abbrev)+']')
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
  integer(POINTER_SIZE), intent(in) :: out
  type(hamiltonian_type), intent(in) :: h
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter

  integer :: i
  real(r8) :: field(3)
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
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
    
    aux = '['+trim(units_out%energy%abbrev) + ' / ' + trim(units_inp%length%abbrev) + ']'
    do i = 1, conf%dim
      call write_iter_header(out, aux)
    end do
    
    aux = '[1/'+ trim(units_inp%length%abbrev) + ']'
    do i = 1, conf%dim
      call write_iter_header(out, aux)
    end do
    call write_iter_nl(out)
  end if
  
  call write_iter_start(out)
  
  field = M_ZERO
  call laser_field(h%no_lasers, h%lasers, iter*td%dt, field)
  field = field * units_inp%length%factor / units_inp%energy%factor
  call write_iter_double(out, field, conf%dim)
  
  call laser_vector_field(h%no_lasers, h%lasers, iter*td%dt, field)
  field = field  * units_inp%length%factor
  call write_iter_double(out, field, conf%dim)
  
  call write_iter_nl(out)
  
end subroutine td_write_laser
    
subroutine td_write_el_energy(out, h, td, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(hamiltonian_type), intent(in) :: h
  type(td_type),     intent(in) :: td
  integer,           intent(in) :: iter

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
    call write_iter_header(out, '['+trim(units_out%time%abbrev)+']')
    do i = 1, 5
      call write_iter_header(out, '['+trim(units_out%energy%abbrev)+']')
    end do
    call write_iter_nl(out)
  endif
  
  call write_iter_start(out)
  call write_iter_double(out, h%etot/units_out%acceleration%factor, 1)
  call write_iter_double(out, h%eii /units_out%acceleration%factor, 1)
  call write_iter_double(out, h%ex  /units_out%acceleration%factor, 1)
  call write_iter_double(out, h%ec  /units_out%acceleration%factor, 1)
  call write_iter_double(out, h%epot/units_out%acceleration%factor, 1)
  call write_iter_nl(out)
  
end subroutine td_write_el_energy

subroutine td_write_proj(out, sys, u_st, iter)
  integer(POINTER_SIZE), intent(in) :: out
  type(system_type), intent(in) :: sys
  type(states_type), intent(in) :: u_st
  integer,           intent(in) :: iter

  complex(r8), allocatable :: projections(:,:,:)
  character(len=20) :: aux
  integer :: ik, ist, uist

  if(iter == 0) then
    ! empty file
    call write_iter_clear(out)
    
    ! first line -> column names
    call write_iter_header_start(out)
    do ik = 1, sys%st%nik
      do ist = 1, sys%st%nst
        do uist = 1, u_st%nst
          write(aux, '(i3,a,i3)') ist, ' -> ', uist
          call write_iter_header(out, 'Re {'//trim(aux)//'}')
          call write_iter_header(out, 'Im {'//trim(aux)//'}')
        end do
      end do
    end do
    call write_iter_nl(out)
  endif

  allocate(projections(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik))
  call calc_projection(u_st, sys%st, sys%m, projections)

  call write_iter_start(out)
  do ik = 1, sys%st%nik
    do ist = 1, sys%st%nst
      do uist = 1, u_st%nst
        call write_iter_double(out,  real(projections(uist, ist, ik)), 1)
        call write_iter_double(out, aimag(projections(uist, ist, ik)), 1)
      end do
    end do
  end do
  call write_iter_nl(out)

  deallocate(projections)
end subroutine td_write_proj
