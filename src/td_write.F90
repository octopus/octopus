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

  subroutine td_write_data(iter)
    integer, intent(in) :: iter
    character(len=50) :: filename
 
    call push_sub('td_write_data')

    ! first resume file
    write(filename, '(a,i3.3)') "tmp/restart.td.", mpiv%node
    call zstates_write_restart(trim(filename), sys%m, sys%st, &
         iter=iter, v1=td%v_old(:, :, 1), v2=td%v_old(:, :, 2))
    
    ! calculate projection onto the ground state
    if(td%out_gsp) call td_write_gsp(iter)
    
    if(mpiv%node==0) then
      if(td%out_multip) call write_iter_flush(out_multip)
      if(td%out_coords) call write_iter_flush(out_coords)
      if(td%out_gsp)    call write_iter_flush(out_gsp)
      if(td%out_acc)    call write_iter_flush(out_acc)
      if(td%out_laser)  call write_iter_flush(out_laser)
      if(td%out_energy) call write_iter_flush(out_energy)
    end if
    
    ! now write down the rest
    write(filename, '(a,i7.7)') "td.", iter  ! name of directory
    call zstates_output(sys%st, sys%m, filename, sys%outp)
    if(sys%outp%what(output_geometry)) call system_write_xyz(filename, "geometry", sys)
    call hamiltonian_output(h, sys%m, filename, sys%outp)
    
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
    call PES_output(td%PESv, sys%m, sys%st, iter, sys%outp%iter, td%dt)
#endif

    call pop_sub()
  end subroutine td_write_data

  subroutine td_write_multipole(iter)
    integer, intent(in) :: iter

    integer :: is, l, m, add_lm
    character(len=50) :: aux
    real(r8), allocatable :: dipole(:), multipole(:,:)
    
    if(mpiv%node.ne.0) return ! only first node outputs

    if(iter == 0) then
      ! empty file
      call write_iter_clear(out_multip)

      ! first line
      write(aux, '(a10,i2,a8,i2)') '# nspin = ', sys%st%nspin, ' lmax = ', td%lmax
      call write_iter_string(out_multip, aux)
      call write_iter_nl(out_multip)

      ! second line -> columns name
      call write_iter_header_start(out_multip)

      do is = 1, sys%st%nspin
        write(aux, '(a,i1,a)') 'dipole(', is, ')'
        call write_iter_header(out_multip, aux)
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
            call write_iter_header(out_multip, aux)
          end do
        end do
      end do
      call write_iter_nl(out_multip)

      ! third line -> units
      call write_iter_string(out_multip, '##########')
      call write_iter_header(out_multip, '['+trim(units_out%time%abbrev)+']')

      do is = 1, sys%st%nspin
        call write_iter_header(out_multip, '['+trim(units_out%length%abbrev)+']')
      end do

      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            select case(l)
            case(0)
              call write_iter_header(out_multip, ' ')
            case(1)
              call write_iter_header(out_multip, '['+trim(units_out%length%abbrev)+']')
            case default
              write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
              call write_iter_header(out_multip, '['+trim(aux)+']')
            end select
          end do
        end do
      end do
      call write_iter_nl(out_multip)
    end if
    
    allocate(dipole(sys%st%nspin), multipole((td%lmax + 1)**2, sys%st%nspin))
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, dipole, multipole)
   
    call write_iter_start(out_multip)
    do is = 1, sys%st%nspin
      call write_iter_double(out_multip, dipole(is)/units_out%length%factor, 1)
    end do
    do is = 1, sys%st%nspin
      add_lm = 1
      do l = 0, td%lmax
        do m = -l, l
          call write_iter_double(out_multip, multipole(add_lm, is)/units_out%length%factor**l, 1)
          add_lm = add_lm + 1
        end do
      end do
    end do
    call write_iter_nl(out_multip)
    
    deallocate(dipole, multipole)
  end subroutine td_write_multipole

  subroutine td_write_nbo(iter, ke, pe)
    integer, intent(in) :: iter
    real(r8),intent(in) :: ke, pe

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

      do i = 1, sys%natoms
        do j = 1, conf%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'x(', i, ',', j, ')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, sys%natoms
        do j = 1, conf%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'v(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      do i = 1, sys%natoms
        do j = 1, conf%dim
          write(aux, '(a2,i3,a1,i3,a1)') 'f(', i, ',',j,')'
          call write_iter_header(out_coords, aux)
        end do
      end do
      call write_iter_nl(out_coords)
      
      ! second line: units
      call write_iter_string(out_coords, '##########')
      call write_iter_header(out_coords, '['+trim(units_out%time%abbrev)+']')
      call write_iter_string(out_coords, &
           'Energy in '      + trim(units_out%energy%abbrev)   +   &
           ', Positions in ' + trim(units_out%length%abbrev)   +   &
           ', Velocities in '+ trim(units_out%velocity%abbrev) + &
           ', Forces in '    + trim(units_out%force%abbrev))
      call write_iter_nl(out_coords)
    end if

    call write_iter_start(out_coords)
    call write_iter_double(out_coords,     ke /units_out%energy%factor, 1)
    call write_iter_double(out_coords,     pe /units_out%energy%factor, 1)
    call write_iter_double(out_coords, (ke+pe)/units_out%energy%factor, 1)

    do i = 1, sys%natoms
      call write_iter_double(out_coords, sys%atom(i)%x(1:conf%dim)/units_out%length%factor,   conf%dim)
    end do
    do i = 1, sys%natoms
      call write_iter_double(out_coords, sys%atom(i)%v(1:conf%dim)/units_out%velocity%factor, conf%dim)
    end do
    do i = 1, sys%natoms
      call write_iter_double(out_coords, sys%atom(i)%f(1:conf%dim)/units_out%force%factor,    conf%dim)
    end do
    call write_iter_nl(out_coords)
    
  end subroutine td_write_nbo

  subroutine td_write_gsp(iter)
    integer, intent(in)  :: iter
    complex(r8) :: gsp

    call push_sub('td_write_gsp')

    ! all processors calculate the projection
    call zstates_project_gs(sys%st, sys%m, gsp)
    print *, oct_getmem()

    ! but only first node outputs
    if(mpiv%node.ne.0) return

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
      call write_iter_header(out_gsp, '['+trim(units_out%time%abbrev)+']')
      call write_iter_nl(out_gsp)
    end if
    
    ! can not call write_iter_start, for the step is not 1
    call write_iter_int(out_gsp, iter, 1)
    call write_iter_double(out_gsp, iter*td%dt/units_out%time%factor,  1)
    call write_iter_double(out_gsp, real(gsp),  1)
    call write_iter_double(out_gsp, aimag(gsp), 1)
    call write_iter_nl(out_gsp)

    call pop_sub()
  end subroutine td_write_gsp

  subroutine td_write_acc(iter)
    integer, intent(in)  :: iter

    integer :: i
    character(len=7) :: aux
    real(r8) :: acc(3)

    if(mpiv%node.ne.0) return ! only first node outputs

    if(iter == 0) then
      ! empty file
      call write_iter_clear(out_acc)

      ! first line -> column names
      call write_iter_header_start(out_acc)
      do i = 1, conf%dim
        write(aux, '(a4,i1,a1)') 'Acc(', i, ')'
        call write_iter_header(out_acc, aux)
      end do
      call write_iter_nl(out_acc)

      ! second line: units
      call write_iter_string(out_acc, '##########')
      call write_iter_header(out_acc, '['+trim(units_out%time%abbrev)+']')
      do i = 1, conf%dim
        call write_iter_header(out_acc, '['+trim(units_out%acceleration%abbrev)+']')
      end do
      call write_iter_nl(out_acc)
    endif

    call td_calc_tacc(acc, td%dt*i, reduce = .true.)

    call write_iter_start(out_acc)
    call write_iter_double(out_acc, acc/units_out%acceleration%factor, conf%dim)
    call write_iter_nl(out_acc)

  end subroutine td_write_acc

  subroutine td_write_laser(iter)
    integer, intent(in) :: iter

    integer :: i
    real(r8) :: field(3)
    character(len=80) :: aux

    if(mpiv%node.ne.0) return ! only first node outputs

    ! TODO -> confirm these stupid units, especially for the vector field
    if(iter == 0) then
      ! empty file
      call write_iter_clear(out_laser)

      ! first line
      write(aux, '(a7,e20.12,3a)') '# dt = ', td%dt/units_out%time%factor, &
           " [", trim(units_out%time%abbrev), "]"
      call write_iter_string(out_laser, aux)
      call write_iter_nl(out_laser)

      ! second line -> column names
      call write_iter_header_start(out_laser)
      do i = 1, conf%dim
        write(aux, '(a,i1,a)') 'E(', i, ')'
        call write_iter_header(out_laser, aux)
      end do
      do i = 1, conf%dim
        write(aux, '(a,i1,a)') 'A(', i, ')'
        call write_iter_header(out_laser, aux)
      end do
      call write_iter_nl(out_laser)

      ! third line -> units
      call write_iter_string(out_laser, '##########')
      call write_iter_header(out_laser, '['+trim(units_out%time%abbrev)+']')

      aux = '['+trim(units_out%energy%abbrev) + ' / ' + trim(units_inp%length%abbrev) + ']'
      do i = 1, conf%dim
        call write_iter_header(out_laser, aux)
      end do

      aux = '[1/'+ trim(units_inp%length%abbrev) + ']'
      do i = 1, conf%dim
        call write_iter_header(out_laser, aux)
      end do
      call write_iter_nl(out_laser)
    end if

    call write_iter_start(out_laser)

    field = M_ZERO
    call laser_field(h%no_lasers, h%lasers, iter*td%dt, field)
    field = field * units_inp%length%factor / units_inp%energy%factor
    call write_iter_double(out_laser, field, conf%dim)
    
    call laser_vector_field(h%no_lasers, h%lasers, iter*td%dt, field)
    field = field  * units_inp%length%factor
    call write_iter_double(out_laser, field, conf%dim)

    call write_iter_nl(out_laser)

  end subroutine td_write_laser
    
  subroutine td_write_el_energy(iter)
    integer, intent(in)  :: iter

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
      call write_iter_header(out_energy, '['+trim(units_out%time%abbrev)+']')
      do i = 1, 5
        call write_iter_header(out_energy, '['+trim(units_out%energy%abbrev)+']')
      end do
      call write_iter_nl(out_energy)
    endif

    call write_iter_start(out_energy)
    call write_iter_double(out_energy, h%etot/units_out%acceleration%factor, 1)
    call write_iter_double(out_energy, h%eii /units_out%acceleration%factor, 1)
    call write_iter_double(out_energy, h%ex  /units_out%acceleration%factor, 1)
    call write_iter_double(out_energy, h%ec  /units_out%acceleration%factor, 1)
    call write_iter_double(out_energy, h%epot/units_out%acceleration%factor, 1)
    call write_iter_nl(out_energy)

  end subroutine td_write_el_energy

  subroutine td_write_proj(iter)
    integer, intent(in)  :: iter

    complex(r8), allocatable :: projections(:,:,:)
    character(len=20) :: aux
    integer :: ik, ist, uist

    if(iter == 0) then
      ! empty file
      call write_iter_clear(out_proj)

      ! first line -> column names
      call write_iter_header_start(out_proj)
      do ik = 1, sys%st%nik
        do ist = 1, sys%st%nst
          do uist = 1, u_st%nst
            write(aux, '(i3,a,i3)') ist, ' -> ', uist
            call write_iter_header(out_proj, 'Re {'//trim(aux)//'}')
            call write_iter_header(out_proj, 'Im {'//trim(aux)//'}')
          end do
        end do
      end do
      call write_iter_nl(out_proj)
    endif

    allocate(projections(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik))
    call calc_projection(u_st, sys%st, sys%m, projections)

    call write_iter_start(out_proj)
    do ik = 1, sys%st%nik
      do ist = 1, sys%st%nst
        do uist = 1, u_st%nst
          call write_iter_double(out_proj,  real(projections(uist, ist, ik)), 1)
          call write_iter_double(out_proj, aimag(projections(uist, ist, ik)), 1)
        end do
      end do
    end do
    call write_iter_nl(out_proj)

    deallocate(projections)
  end subroutine td_write_proj
