!! Copyright (C) 2019 X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the Lesser GNU General Public License as published by
!! the Free Software Foundation; either version 3, or (at your option)
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

  
module alloc_cache_oct_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use global_oct_m
  use iso_c_binding

  implicit none

  private

  public ::                                  &
    alloc_cache_t,                           &
    alloc_cache_init,                        &
    alloc_cache_end,                         &
    alloc_cache_put,                         &
    alloc_cache_get
  
  type alloc_cache_t
    private
    integer(8) :: dummy
  end type alloc_cache_t

  integer(8), public, parameter ::              &
    ALLOC_CACHE_ANY_SIZE = -1_8
  
  interface

    subroutine alloc_cache_init(alloc_cache, max_size)
      import :: alloc_cache_t
      implicit none

      type(alloc_cache_t),   intent(out) :: alloc_cache
      integer(8),            intent(in)  :: max_size
    end subroutine alloc_cache_init

    ! -------------------------------------------------
    
    subroutine alloc_cache_end(alloc_cache, hits, misses, vol_hits, vol_misses)
      use iso_c_binding
      import :: alloc_cache_t
      implicit none

      type(alloc_cache_t),   intent(inout) :: alloc_cache
      integer(8),            intent(out)   :: hits
      integer(8),            intent(out)   :: misses
      real(c_double),        intent(out)   :: vol_hits
      real(c_double),        intent(out)   :: vol_misses
    end subroutine alloc_cache_end

  end interface

contains

  subroutine alloc_cache_put(alloc_cache, size, loc, put)
    type(alloc_cache_t),   intent(inout) :: alloc_cache
    integer(8),            intent(in)    :: size
#ifdef HAVE_OPENCL
    type(cl_mem),          intent(in)    :: loc
#else
    type(c_ptr),           intent(in)    :: loc
#endif
    logical,               intent(out)   :: put

    interface
      subroutine alloc_cache_put_low(alloc_cache, size, loc, put)
#ifdef HAVE_OPENCL
        use cl
#endif
        use iso_c_binding
        import :: alloc_cache_t
        implicit none

        type(alloc_cache_t),   intent(inout) :: alloc_cache
        integer(8),            intent(in)    :: size
#ifdef HAVE_OPENCL
        type(cl_mem),          intent(in)    :: loc
#else
        type(c_ptr),           intent(in)    :: loc
#endif
        integer,               intent(out)   :: put
      end subroutine alloc_cache_put_low
    end interface

    integer :: iput
    
    call alloc_cache_put_low(alloc_cache, size, loc, iput)

    put = (iput /= 0)
    
  end subroutine alloc_cache_put

  ! -------------------------------------------------
  
  subroutine alloc_cache_get(alloc_cache, size, found, loc)
      type(alloc_cache_t),   intent(inout) :: alloc_cache
      integer(8),            intent(in)    :: size
      logical,               intent(out)   :: found
#ifdef HAVE_OPENCL
      type(cl_mem),          intent(out)   :: loc
#else
      type(c_ptr),           intent(out)   :: loc
#endif

    interface
      subroutine alloc_cache_get_low(alloc_cache, size, found, loc)
#ifdef HAVE_OPENCL
        use cl
#endif
        use iso_c_binding
        import :: alloc_cache_t
        implicit none
        
        type(alloc_cache_t),   intent(inout) :: alloc_cache
        integer(8),            intent(in)    :: size
        integer,               intent(out)   :: found
#ifdef HAVE_OPENCL
        type(cl_mem),          intent(out)   :: loc
#else
        type(c_ptr),           intent(out)   :: loc
#endif
      end subroutine alloc_cache_get_low
    end interface

    integer :: ifound

    call alloc_cache_get_low(alloc_cache, size, ifound, loc)

    found = (ifound /= 0)
    
  end subroutine alloc_cache_get
  

end module alloc_cache_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
