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

program make_st
  use global
  use messages
  use syslabels
  use io
  use units
  use lib_oct_parser
  use mesh
  use grid
  use states
  use restart
  use system

  implicit none

  type(system_type) :: sys
  integer :: i, n, ik, ist, idim, type, err
  integer(POINTER_SIZE) :: blk

  ! Initialize stuff
  call global_init()
  call parser_init()
  call io_init()
  call syslabels_init(1)
  if(in_debug_mode) then
    call io_mkdir('debug')
  end if
  call units_init()
  call system_init(sys)

#ifdef T_REAL
  deallocate(sys%st%dpsi)
#endif

  allocate(sys%st%zpsi (sys%gr%m%np_part, sys%st%d%dim, sys%st%nst, sys%st%d%nik), &
    sys%st%eigenval(sys%st%nst, sys%st%d%nik))

  call X(restart_read)(trim(tmpdir)//'restart_gs', sys%st, sys%gr%m, err)
  if(err < 0) then
    message(1) = "Error opening 'restart.static' file"
    call write_fatal(1)
  end if

  if(loct_parse_block(check_inp('MakeStates'), blk).ne.0) then
    message(1) = "Block '%MakeStates' must be defined"
    call write_fatal(1)
  end if

  n = loct_parse_block_n(blk)
  do i = 1, n
    call loct_parse_block_int(blk, i-1, 0, ik)
    call loct_parse_block_int(blk, i-1, 1, ist)
    call loct_parse_block_int(blk, i-1, 2, idim)
    call loct_parse_block_int(blk, i-1, 3, type)
    select case(type)
    case(1)
      call wf_gaussian(sys%gr, sys%st, i-1)
    end select
  end do
  call loct_parse_block_end(blk)

  ! renormalize functions
  call wf_renormalize(sys%gr%m, sys%st)

  ! save wfs in a new static file
  call X(restart_write) (trim(tmpdir)//'restart_gs_new', sys%st, sys%gr, err)
  if(err.ne.0) then
    message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_gs_new"'
    call write_fatal(1)
  end if

  call system_end(sys)
  call syslabels_end()
  call io_end()
  call parser_end()
  call global_end()

contains

  ! ---------------------------------------------------------
  subroutine wf_gaussian(gr, st, line)
    type(grid_type),   intent(in)    :: gr
    type(states_type), intent(inout) :: st
    integer,           intent(in)    :: line

    FLOAT :: x1(3), x(3), s, k(3)
    integer :: i, j

    ! read gaussian parameters
    x1 = M_ZERO; k = M_ZERO
    call loct_parse_block_float(blk, line, 4,  s)
    s = s * units_inp%length%factor

    j = 5
    do i = 1, NDIM
      call loct_parse_block_float(blk, line, j,  x1(i))
      x1(i) = x1(i) * units_inp%length%factor
      j = j + 1
    end do

    do i = 1, NDIM
      call loct_parse_block_float(blk, line, j,  k(i))
      k(i) = k(i) / units_inp%length%factor
      j = j + 1
    end do

    ! build a gaussian
    do i = 1, NP
      x = gr%m%x(i,:)
      st%zpsi(i, idim, ist, ik) = exp(-sum((x(1:NDIM)-x1(1:NDIM))**2)/(2*s*s) + &
        M_zI*sum(k(1:NDIM)*(x(1:NDIM)-x1(1:NDIM))))
    end do

  end subroutine wf_gaussian


  ! ---------------------------------------------------------
  subroutine wf_renormalize(m, st)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(inout) :: st

    integer :: ik, ist
    FLOAT :: nrm2

    do ik = 1, st%d%nik
      do ist = 1, st%nst
        nrm2 = zstates_nrm2 (m, st%d%dim, st%zpsi(:,:, ist, ik))
        st%zpsi(:,:, ist, ik) = st%zpsi(:,:, ist, ik)/sqrt(nrm2)
      end do
    end do
  end subroutine wf_renormalize

end program make_st
