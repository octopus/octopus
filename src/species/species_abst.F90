!! Copyright (C) 2020 N. Tancogne-Dejean
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
!! along with st program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module species_abst_oct_m
  use global_oct_m
  use loct_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use namespace_oct_m

  implicit none

  private

  public ::                           &
    species_abst_t,                   &
    species_iterator_t

  integer, public, parameter :: LABEL_LEN=15

  type, abstract :: species_abst_t

  character(len=LABEL_LEN) :: label !< Identifier for the species

  contains
    procedure(species_abst_debug),      deferred :: debug
    procedure get_label => species_abst_label
  end type species_abst_t
 
  abstract interface
    subroutine species_abst_debug(this, namespace, dir)
      import species_abst_t
      import namespace_t
      class(species_abst_t),   intent(inout) :: this
      type(namespace_t),       intent(in)    :: namespace
      character(len=*),        intent(in)    :: dir
    end subroutine species_abst_debug
  end interface

 
  !> This class extends the list iterator and adds one method to get the
  !! interaction as a pointer of type class(species_abst_t).
  type, extends(list_iterator_t) :: species_iterator_t
    private
  contains
    procedure :: get_next_interaction => species_iterator_get_next
  end type species_iterator_t

  contains

    ! ---------------------------------------------------------
    function species_iterator_get_next(this) result(value)
      class(species_iterator_t), intent(inout) :: this
      class(species_abst_t),     pointer       :: value

      class(*), pointer :: ptr

      PUSH_SUB(species_iterator_get_next)

      ptr => this%get_next()
      select type (ptr)
      class is (species_abst_t)
        value => ptr
      class default
        ASSERT(.false.)
      end select

      POP_SUB(species_iterator_get_next)
    end function species_iterator_get_next

   ! ---------------------------------------------------------
   character(len=LABEL_LEN) pure function species_abst_label(this)
     class(species_abst_t), intent(in) :: this

     species_abst_label = trim(this%label)
   end function species_abst_label


end module species_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

