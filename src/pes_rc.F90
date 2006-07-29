!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine PES_rc_init(v, m, st, save_iter)
  type(mesh_t),   intent(in) :: m
  type(states_t), intent(in) :: st
  integer,        intent(in) :: save_iter
  type(PES_rc_t), intent(out) :: v

  integer(POINTER_SIZE) :: blk
  integer  :: i
  FLOAT ::  x(MAX_DIM)

  message(1) = 'Info: Calculating PES using rc technique'
  call write_info(1)

  !%Variable PES_rc_points
  !%Type block
  !%Section Time Dependent::PES
  !%Description
  !% List of points where to calculate the photo-electron spectrum a la Suraud.
  !% The exact syntax is:
  !%
  !% <tt>%TDPES_rc_points
  !% <br>&nbsp;&nbsp;x1 | y1 | z1
  !% <br>%
  !% </tt>
  !%End
  if (loct_parse_block(check_inp('PES_rc_points'), blk) < 0) call input_error('PES_rc_points')

  v%npoints = loct_parse_block_n(blk)

  ! setup filenames and read points
  ALLOCATE(v%filenames(v%npoints), v%npoints)
  ALLOCATE(v%points   (v%npoints), v%npoints)
  do i = 1, v%npoints
    write(v%filenames(i), '(a,i2.2,a)') 'PES_rc.', i, '.out'

    call loct_parse_block_float(blk, i-1, 0, x(1))
    call loct_parse_block_float(blk, i-1, 1, x(2))
    call loct_parse_block_float(blk, i-1, 2, x(3))

    ! adjust units
    x = x*units_inp%length%factor

    v%points(i) = m%Lxyz_inv(int(x(1)/m%h(1)), int(x(2)/m%h(2)), int(x(3)/m%h(3)))
  end do

  call loct_parse_block_end(blk)

  i = v%npoints*st%d%dim*(st%st_end-st%st_start+1)*st%d%nik*save_iter
  ALLOCATE(v%wf(v%npoints, st%d%dim, st%st_start:st%st_end, st%d%nik, save_iter), i)
end subroutine PES_rc_init


! ---------------------------------------------------------
subroutine PES_rc_end(v)
  type(PES_rc_t), intent(inout) :: v

  if(associated(v%filenames)) then
    deallocate(v%filenames, v%points, v%wf)
    nullify   (v%filenames, v%points, v%wf)
  end if
end subroutine PES_rc_end


! ---------------------------------------------------------
subroutine PES_rc_doit(v, st, ii)
  type(PES_rc_t), intent(inout) :: v
  type(states_t), intent(in) :: st
  integer,        intent(in) :: ii

  integer :: ix, ik, p, idim

  do ix = 1, v%npoints
    do ik = 1, st%d%nik
      do p = st%st_start, st%st_end
        do idim = 1, st%d%dim
          v%wf(ix, idim, p, ik, ii) = st%occ(p, ik)*st%zpsi(v%points(ix), idim, p, ik)
        end do
      end do
    end do
  end do

end subroutine PES_rc_doit


! ---------------------------------------------------------
subroutine PES_rc_output(v, st, iter, save_iter, dt)
  type(PES_rc_t), intent(in) :: v
  type(states_t), intent(in) :: st
  integer,        intent(in) :: iter, save_iter
  FLOAT,          intent(in) :: dt

  integer :: ix, iunit, j, jj, ik, p, idim

  if(iter == 1) then
    do ix = 1, v%npoints
      iunit = io_open(v%filenames(ix), action='write')
      write(iunit, '(a7,f17.6,3a)') &
        '# dt = ', dt/units_inp%time%factor, ' [', trim(units_inp%time%abbrev), ']'
      write(iunit, '(a3,14x)', advance='no') '# t'
      do ik = 1, st%d%nik
        do p = st%st_start, st%st_end
          do idim = 1, st%d%dim
            write(iunit, '(3x,a8,i3,a7,i3,a8,i3,3x)', advance='no') &
              "ik = ", ik, " ist = ", p, " idim = ", idim
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
      call io_close(iunit)
    end do
  end if

  do ix = 1, v%npoints
    call io_assign(iunit)
    iunit= io_open(v%filenames(ix), action='write', position='append')
    do j = 1, save_iter
      jj = iter - save_iter + j
      write(iunit, '(e17.10)', advance='no') jj*dt/units_inp%time%factor
      do ik = 1, st%d%nik
        do p = st%st_start, st%st_end
          do idim = 1, st%d%dim
            write(iunit, '(1x,a1,e17.10,a1,e17.10,a1)', advance='no') &
              '(', real(v%wf(ix, idim, p, ik, j))*units_out%length%factor**1.5, &
              ',', aimag(v%wf(ix, idim, p, ik, j))*units_out%length%factor**1.5, ')'
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
    end do
    call io_close(iunit)
  end do

end subroutine PES_rc_output
