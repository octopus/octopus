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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: io.F90 3613 2007-11-29 16:47:41Z xavier $

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

module iihash_m
  use ialist_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::        &
    iihash_t,      &
    iihash_init,   &
    iihash_end,    &
    iihash_insert, &
    iihash_lookup, &
    get_next_prime

  type iihash_t
    integer                 :: size
    type(ialist_t), pointer :: keyval(:)
  end type iihash_t
  
contains

  ! ---------------------------------------------------------
  !> Initialize a hash table h with size entries. Since we use separate
  !! chaining, the number of entries in the hash table is, in
  !! principle, unlimited. We take the smallest prime number as table
  !! size that is greater or equal than the requested size to reduce
  !! collisions.
  subroutine iihash_init(h, size)
    type(iihash_t), intent(out) :: h
    integer,        intent(in)  :: size

    integer :: prime_size, i, min_size

    PUSH_SUB(iihash_init)

    if(size.lt.2) then
      min_size = 3
    else
      min_size = size
    end if

    prime_size = get_next_prime(min_size)
    SAFE_ALLOCATE(h%keyval(0:prime_size-1))
    do i = 0, prime_size-1
      call ialist_init(h%keyval(i))
    end do
    h%size = prime_size

    POP_SUB(iihash_init)
  end subroutine iihash_init


  ! ---------------------------------------------------------
  !> Free a hash table.
  subroutine iihash_end(h)
    type(iihash_t), intent(inout) :: h

    integer :: i

    PUSH_SUB(iihash_end)

    do i = 0, h%size-1
      call ialist_end(h%keyval(i))
    end do
    SAFE_DEALLOCATE_P(h%keyval)

    POP_SUB(iihash_end)
  end subroutine iihash_end


  ! ---------------------------------------------------------
  !> Insert a (key, val) pair into the hash table h.
  subroutine iihash_insert(h, key, val)
    type(iihash_t),    intent(inout) :: h
    integer,           intent(in)    :: key
    integer,           intent(in)    :: val

    integer :: k

    ASSERT(key.ge.0)

    k = hash(h, key)
    call ialist_insert(key, val, h%keyval(k))
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

    integer :: k

    k = hash(h, key)

    iihash_lookup = ialist_lookup(key, h%keyval(k), found)
  end function iihash_lookup


  ! ---------------------------------------------------------
  !> The hash function.
  integer function hash(h, key)
    type(iihash_t), intent(in) :: h
    integer,        intent(in) :: key

    hash = mod(key, h%size)
  end function hash


  ! ---------------------------------------------------------
  !> Returns the smallest prime number that is greater than k
  !! using the Sieve of Eratosthenes.
  integer function get_next_prime(k)
    integer, intent(in) :: k

    integer          :: size, i, j
    logical, pointer :: primes(:), primes1(:)
    logical          :: found, searching

    found     = .false.
    searching = .true.
    size      = k+k/4

    SAFE_ALLOCATE(primes(1:size))

    primes    = .true.
    primes(1) = .false.

    j = 1

    do while(searching)
      call sieve()

      ! Search for prime greater or equal k.
      do i = k, size
        if(primes(i)) then
          found = .true.
          exit
        end if
      end do

      ! If no suitable prime could be found, enlarge
      ! the array and continue sieving.
      if(found) then
        searching = .false.
        get_next_prime = i
      else
        SAFE_ALLOCATE(primes1(1:size+size/4))
        primes1(1:size)             = primes
        primes1(size+1:size+size/4) = .true.
        size                        = size+size/4
        SAFE_DEALLOCATE_P(primes)
        primes => primes1
      end if
    end do
    SAFE_DEALLOCATE_P(primes)

  contains

    ! This routine implements the Sieve of Eratosthenes.
    subroutine sieve()
      do while(j.lt.size)
        do while(.not.primes(j).and.j.lt.size)
          j = j + 1
        end do
        i = 2*j
        do while(i.le.size)
          primes(i) = .false.
          i = i + j
        end do
        j = j + 1
      end do
    end subroutine sieve
  end function get_next_prime
end module iihash_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
