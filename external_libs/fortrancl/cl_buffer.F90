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
  !! $Id: cl.F90 3587 2007-11-22 16:43:00Z xavier $

#include "config_F90.h"
 
module cl_buffer_m
  use cl_types_m

  implicit none 

  private

  public ::                          &
    flCreateBuffer,                  &
    flReleaseMemObject

  interface

    subroutine flReleaseMemObject(memobj, status)
      use cl_types_m

      implicit none

      type(cl_mem),           intent(inout) :: memobj
      integer,                intent(out)   :: status
    end subroutine flReleaseMemObject

  end interface

  contains 

    type(cl_mem) function flCreateBuffer(context, flags, size, errcode_ret) result(buffer)
      type(cl_context), intent(in)    :: context
      integer,          intent(in)    :: flags
      integer(8),       intent(in)    :: size
      integer,          intent(out)   :: errcode_ret
      
      interface
        
        subroutine flCreateBuffer_low(context, flags, size, errcode_ret, buffer)
          use cl_types_m
          
          implicit none
          
          type(cl_context),        intent(in)    :: context
          integer,                 intent(in)    :: flags
          integer(SIZEOF_SIZE_T),  intent(in)    :: size
          integer,                 intent(out)   :: errcode_ret
          type(cl_mem),            intent(out)   :: buffer
        end subroutine flCreateBuffer_low
        
      end interface
  
      call flCreateBuffer_low(context, flags, int(size, SIZEOF_SIZE_T), errcode_ret, buffer)
      
    end function flCreateBuffer

end module cl_buffer_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
