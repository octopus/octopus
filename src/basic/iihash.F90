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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module implements a simple hash table for non-negative integer keys
!! and integer values.
!!
!! The collision resolution used is separate chaining (cf. Knuth, 1973, vol. 3)
!! using associative lists. The reason to use separate chaining is that we
!! look up a lot of keys not in the table and, e.g., open addressing is
!! very slow in that case.
!! The hash function is simply (key mod size) but size is taken to be a prime
!! number, i.e. the table is usually slightly larger than the user requests.

module iihash_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private
  public ::        &
    iihash_t,      &
    iihash_init,   &
    iihash_end,    &
    iihash_insert, &
    iihash_lookup, &
    lihash_t,      &
    lihash_init,   &
    lihash_end,    &
    lihash_insert, &
    lihash_lookup

  type iihash_t
    private
    type(c_ptr) :: map
  end type iihash_t

  type lihash_t
    private
    type(c_ptr) :: map
  end type lihash_t

contains

  ! ---------------------------------------------------------
  !> Initialize a hash table h
  subroutine iihash_init(h)
    type(iihash_t), intent(out) :: h

    interface
      subroutine iihash_map_init(map)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
      end subroutine iihash_map_init
    end interface


    PUSH_SUB(iihash_init)

    call iihash_map_init(h%map)

    POP_SUB(iihash_init)
  end subroutine iihash_init


  ! ---------------------------------------------------------
  !> Free a hash table.
  subroutine iihash_end(h)
    type(iihash_t), intent(inout) :: h

    interface
      subroutine iihash_map_end(map)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
      end subroutine iihash_map_end
    end interface

    PUSH_SUB(iihash_end)

    call iihash_map_end(h%map)

    POP_SUB(iihash_end)
  end subroutine iihash_end


  ! ---------------------------------------------------------
  !> Insert a (key, val) pair into the hash table h.
  subroutine iihash_insert(h, key, val)
    type(iihash_t),    intent(inout) :: h
    integer,           intent(in)    :: key
    integer,           intent(in)    :: val

    interface
      subroutine iihash_map_insert(map, key, val)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
        integer,     intent(in)    :: key
        integer,     intent(in)    :: val
      end subroutine iihash_map_insert
    end interface

    call iihash_map_insert(h%map, key, val)
  end subroutine iihash_insert


  ! ---------------------------------------------------------
  !> Look up a value in the hash table h. If found is present, it
  !! indicates if key could be found in the table. If found = .false.,
  !! the return value of iihash_lookup is meaningless (and essentially
  !! undefined).
  integer function iihash_lookup(h, key, found)
    type(iihash_t),    intent(in)  :: h
    integer,           intent(in)  :: key
    logical, optional, intent(out) :: found

    interface
      subroutine iihash_map_lookup(map, key, ifound, val)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(in)    :: map
        integer,     intent(in)    :: key
        integer,     intent(out)   :: ifound
        integer,     intent(out)   :: val
      end subroutine iihash_map_lookup
    end interface

    integer :: ifound, val

    call iihash_map_lookup(h%map, key, ifound, val)

    found = (ifound == 1)

    iihash_lookup = -1
    if(found) iihash_lookup = val

  end function iihash_lookup

  ! ---------------------------------------------------------
  !> Initialize a hash table h
  subroutine lihash_init(h)
    type(lihash_t), intent(out) :: h

    interface
      subroutine lihash_map_init(map)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
      end subroutine lihash_map_init
    end interface


    PUSH_SUB(lihash_init)

    call lihash_map_init(h%map)

    POP_SUB(lihash_init)
  end subroutine lihash_init


  ! ---------------------------------------------------------
  !> Free a hash table.
  subroutine lihash_end(h)
    type(lihash_t), intent(inout) :: h

    interface
      subroutine lihash_map_end(map)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
      end subroutine lihash_map_end
    end interface

    PUSH_SUB(lihash_end)

    call lihash_map_end(h%map)

    POP_SUB(lihash_end)
  end subroutine lihash_end


  ! ---------------------------------------------------------
  !> Insert a (key, val) pair into the hash table h.
  subroutine lihash_insert(h, key, val)
    type(lihash_t),    intent(inout) :: h
    integer(8),        intent(in)    :: key
    integer,           intent(in)    :: val

    interface
      subroutine lihash_map_insert(map, key, val)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(inout) :: map
        integer(8),  intent(in)    :: key
        integer,     intent(in)    :: val
      end subroutine lihash_map_insert
    end interface

    call lihash_map_insert(h%map, key, val)
  end subroutine lihash_insert


  ! ---------------------------------------------------------
  !> Look up a value in the hash table h. If found is present, it
  !! indicates if key could be found in the table. If found = .false.,
  !! the return value of lihash_lookup is meaningless (and essentially
  !! undefined).
  integer function lihash_lookup(h, key, found)
    type(lihash_t),    intent(in)  :: h
    integer(8),        intent(in)  :: key
    logical, optional, intent(out) :: found

    interface
      subroutine lihash_map_lookup(map, key, ifound, val)
        use iso_c_binding
        implicit none

        type(c_ptr), intent(in)    :: map
        integer(8),  intent(in)    :: key
        integer,     intent(out)   :: ifound
        integer,     intent(out)   :: val
      end subroutine lihash_map_lookup
    end interface

    integer :: ifound, val

    call lihash_map_lookup(h%map, key, ifound, val)

    found = (ifound == 1)

    lihash_lookup = -1
    if(found) lihash_lookup = val

  end function lihash_lookup

end module iihash_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
