!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
!! $Id: states.F90 5022 2009-03-03 17:47:58Z nitsche $

#include "global.h"

subroutine modelmb_wf_write(dir, gr, mm, wf)

  !use datasets_m
  use global_m
  use grid_m
  !use io_function_m
  use index_m
  use io_m
  !use lalg_adv_m
  !use loct_m
  !use loct_parser_m
  use messages_m
  !use modelmb_particles_m
  !use mpi_m
  !use mpi_lib_m
  use profiling_m

  implicit none

! args
  character(len=*),   intent(in) :: dir
  type(grid_t),       intent(in) :: gr
  integer,            intent(in) :: mm
  CMPLX,              intent(in) :: wf(1:gr%mesh%np_part_global)

!local
  integer :: ip, idir, iunit
  integer, allocatable :: ix(:)
  character(len=200) :: filename

  call push_sub('states.modelmb_wf_write')

  SAFE_ALLOCATE(ix(1:MAX_DIM))

  write(filename,'(a,i4.4)') trim(dir)//'/wf_imb', mm
  iunit = io_open(filename, action='write')
  write(iunit, '(a)', ADVANCE='no') '#'
  do idir = 1, gr%sb%dim
     write(iunit, '(a)', ADVANCE='no') '      position          '
  end do
   write(iunit, '(a)') '         Re(wf)                    Im(wf)'

  ! for each point in whole box, without boundary points
  do ip = 1, gr%mesh%np_global

    ! get coordinates
    call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)

    ! print out wavefunction
    do idir = 1, gr%sb%dim
      write(iunit,'(es24.15)', ADVANCE='no') ix(idir)*gr%mesh%h(idir)+gr%sb%box_offset(idir)
    end do ! idir
    write(iunit,'(es24.15,es24.15)') real(wf(ip)), aimag(wf(ip))

  end do ! ip

  call io_close(iunit)

  SAFE_DEALLOCATE_A(ix)

  call pop_sub()

end subroutine modelmb_wf_write


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
