!! Copyright (C) 2019 X. Andrade
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

#include "global.h"

module cuda_oct_m
  use global_oct_m
  use io_oct_m
  use loct_oct_m
  use messages_oct_m
  use profiling_oct_m
  
  implicit none

  private
  public ::                    &
    cuda_mem_alloc
  
  interface
    
    ! -------------------------------------------------
    
    subroutine cuda_mem_alloc(cuda_ptr, size)
      use iso_c_binding
      implicit none
      
      type(c_ptr), intent(inout) :: cuda_ptr
      integer(8),  intent(in)    :: size
    end subroutine cuda_mem_alloc
    
  end interface
  
end module cuda_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
