!! Copyright (C) 2010 X. Andrade
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
!! $Id: partitioner.F90 6396 2010-03-26 08:51:58Z mjv500 $

#include "global.h"

module partition_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::                &
       partition_init


contains

  ! ----------------------------------------------------------------------
  subroutine partition_init()
    PUSH_SUB(partition_init)

    POP_SUB(partition_init)
  end subroutine partition_init

end module partition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
