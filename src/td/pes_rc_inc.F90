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
!! $Id$

! ---------------------------------------------------------
subroutine PES_rc_init(v, m, st, save_iter)
  type(mesh_t),   intent(in) :: m
  type(states_t), intent(in) :: st
  integer,        intent(in) :: save_iter
  type(PES_rc_t), intent(out) :: v

  type(block_t) :: blk
  integer  :: i
  FLOAT ::  x(MAX_DIM)

  message(1) = 'Info: Calculating PES using rc technique'
  call write_info(1)

  !%Variable PES_rc_points
  !%Type block
  !%Section Time-Dependent::PES
  !%Description
  !% List of points where to calculate the photoelectron spectrum a la Suraud.
  !% The exact syntax is:
  !%
  !% <tt>%TDPES_rc_points
  !% <br>&nbsp;&nbsp;x1 | y1 | z1
  !% <br>%
  !% </tt>
  !%End
  if (parse_block(datasets_check('PES_rc_points'), blk) < 0) call input_error('PES_rc_points')

  v%npoints = parse_block_n(blk)

  ! setup filenames and read points
  SAFE_ALLOCATE(v%filenames(1:v%npoints))
  SAFE_ALLOCATE(v%points   (1:v%npoints))
  do i = 1, v%npoints
    write(v%filenames(i), '(a,i2.2,a)') 'PES_rc.', i, '.out'

    call parse_block_float(blk, i-1, 0, x(1), units_inp%length)
    call parse_block_float(blk, i-1, 1, x(2), units_inp%length)
    call parse_block_float(blk, i-1, 2, x(3), units_inp%length)

    v%points(i) = m%idx%Lxyz_inv(int(x(1)/m%spacing(1)), int(x(2)/m%spacing(2)), int(x(3)/m%spacing(3)))
  end do

  call parse_block_end(blk)

  SAFE_ALLOCATE(v%wf(1:v%npoints, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, 1:save_iter))
end subroutine PES_rc_init


! ---------------------------------------------------------
subroutine PES_rc_end(v)
  type(PES_rc_t), intent(inout) :: v

  if(associated(v%filenames)) then
    SAFE_DEALLOCATE_P(v%filenames)
    SAFE_DEALLOCATE_P(v%points)
    SAFE_DEALLOCATE_P(v%wf)
  end if
end subroutine PES_rc_end


! ---------------------------------------------------------
subroutine PES_rc_doit(v, st, ii)
  type(PES_rc_t), intent(inout) :: v
  type(states_t), intent(in) :: st
  integer,        intent(in) :: ii

  integer :: ix, ik, p, idim

  do ix = 1, v%npoints
    do ik = st%d%kpt%start, st%d%kpt%end
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
  CMPLX :: vfu

  if(iter == 1) then
    do ix = 1, v%npoints
      iunit = io_open(v%filenames(ix), action='write')
      write(iunit, '(a7,f17.6,3a)') &
        '# dt = ', units_from_atomic(units_inp%time, dt), ' [', trim(units_abbrev(units_inp%time)), ']'
      write(iunit, '(a3,14x)', advance='no') '# t'
      do ik = st%d%kpt%start, st%d%kpt%end
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
      write(iunit, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj*dt)
      do ik = 1, st%d%nik
        do p = st%st_start, st%st_end
          do idim = 1, st%d%dim
            vfu = v%wf(ix, idim, p, ik, j)
            vfu = units_from_atomic(sqrt(units_out%length**3), vfu)
            write(iunit, '(1x,a1,e17.10,a1,e17.10,a1)', advance='no') &
              '(', real(vfu), ',', aimag(vfu), ')'
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
    end do
    call io_close(iunit)
  end do

end subroutine PES_rc_output

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
