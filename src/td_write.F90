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

  subroutine td_write_data()
    integer :: iunit, j, jj, ist, ik, uist
    character(len=50) :: fmt
 
    ! output multipoles
    if(mpiv%node==0) then
    call io_assign(iunit)
    open(iunit, position='append', file="td.general/multipoles")
    do j = 1, td%save_iter
      jj = i - td%save_iter + j
      call td_write_multipole(iunit, jj, jj*td%dt, &
           dipole(:, j), multipole(:,:, j), .false.)
    end do
    call io_close(iunit)
    endif

    ! output the laser field
    if(h%output_laser .and. mpiv%node==0) then
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
      write(fmt,'(i5)') 2*u_st%nst+1
      write(fmt,'(a)') '('//trim(adjustl(fmt))//'f14.8)'
      call io_assign(iunit)
      write(filename, '(a,i3.3)') 'td.general/projections.', mpiv%node
      open(iunit, position='append', file=filename)
      do j = 1, td%save_iter
        jj = i - td%save_iter + j
        do ik = 1, sys%st%nik
          do ist = 1, sys%st%nst
            write(iunit, fmt) jj*td%dt, (projections(uist, ist, ik, j), uist = 1, u_st%nst)
          end do
        end do
      end do
      call io_close(iunit)
    end if

    ! output positions, vels...
    if(td%move_ions > 0 .and. mpiv%node==0) then
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
    if(td%harmonic_spectrum .and. mpiv%node==0) then
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
    real(r8) :: l_field(3), l_vector_field(3)
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

    call laser_field(h%no_lasers, h%lasers, t, l_field)
    call laser_vector_field(h%no_lasers, h%lasers, t, l_vector_field)
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
    !write(iunit,'(i8,4es20.12)') iter, t/units_out%time%factor, acc/units_out%velocity%factor

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
