!! Copyright (C) 2009 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"
#include "io_csv.h"

module io_csv_m
  use global_m
  use messages_m

  implicit none 

  private

  public ::        &
    dread_csv,     &
    io_csv_get_info

contains

  subroutine dread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer(8),          intent(in)  :: np
    real(8),             intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    PUSH_SUB(dread_csv)

    call read_csv(np, ff(1), TYPE_DOUBLE, ierr, trim(fname))

    POP_SUB(dread_csv)
  end subroutine dread_csv
  
  subroutine io_csv_get_info(fname, dims, ierr)
    character(len=*),    intent(in)    :: fname
    integer(8),          intent(inout) :: dims(:)
    integer,             intent(out)   :: ierr

    PUSH_SUB(io_csv_get_info)

    call get_info_csv(dims, ierr, trim(fname))

    POP_SUB(io_csv_get_info)
  end subroutine io_csv_get_info

end module io_csv_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
