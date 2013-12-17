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

module papi_m

  implicit none
  private

  public ::                             &
    papi_init,                          &
    papi_end

  ! functions implemented in papi.c
  interface
    subroutine papi_init()
      implicit none
    end subroutine papi_init

    subroutine papi_end()
      implicit none
    end subroutine papi_end
    
    subroutine papi_get_count_and_reset(fp_ops)
      implicit none
      real(8), intent(out) :: fp_ops
    end subroutine papi_get_count_and_reset
  end interface

end module papi_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
