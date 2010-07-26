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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: global.F90 3587 2007-11-22 16:43:00Z xavier $

#include "global.h"
#include "io_csv.h"

module io_csv_m

  use messages_m

  implicit none 

  private

  public ::             &
    io_csv_read,     &
    io_csv_get_info

  interface io_csv_read
    module procedure sread_csv, dread_csv,  cread_csv, zread_csv, iread_csv, lread_csv
  end interface

contains

  subroutine sread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer(8),             intent(in)   :: np
    real(4),             intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr

    integer, parameter :: type = TYPE_FLOAT

    call push_sub('io_csv_f.sread_csv')

    ierr = 0
    call read_csv(np, ff(1), type, ierr, trim(fname))

    call pop_sub('io_csv_f.sread_csv')
  end subroutine sread_csv

  subroutine dread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)  :: fname
    integer(8),             intent(in)  :: np
    real(8),             intent(out) :: ff(:)
    integer,             intent(out) :: ierr

    integer, parameter :: type = TYPE_DOUBLE

    call push_sub('io_csv_f.dread_csv')

    ierr = 0
    call read_csv(np, ff(1), type, ierr, trim(fname))

    call pop_sub('io_csv_f.dread_csv')
  end subroutine dread_csv
  
  subroutine cread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer(8),             intent(in)   :: np
    complex(4),          intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    
    ierr = 1
  end subroutine cread_csv
  
! some ugly stub functions just to satisfy the linker

  subroutine zread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer(8),             intent(in)   :: np
    complex(8),          intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    
    ierr = 1
  end subroutine zread_csv
  
  subroutine iread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer(8),             intent(in)   :: np
    integer(4),          intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    
    ierr = 1
  end subroutine iread_csv
  
  subroutine lread_csv(fname, np, ff, ierr)
    character(len=*),    intent(in)   :: fname
    integer(8),             intent(in)   :: np
    integer(8),          intent(out)  :: ff(:)
    integer,             intent(out)  :: ierr
    
    ierr = 1
  end subroutine lread_csv

! some ugly stub functions end
  
  subroutine io_csv_get_info(fname, dims, ierr)
    character(len=*),    intent(in)    :: fname
    integer(8),          intent(inout) :: dims(:)
    integer,             intent(out)   :: ierr

    integer :: type
    
    call push_sub('io_csv_f.io_csv_get_info')

    type = 0
    ierr = 0
    call get_info_csv(dims, ierr, trim(fname))

    call pop_sub('io_csv_f.io_csv_get_info')
  end subroutine io_csv_get_info

end module io_csv_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
