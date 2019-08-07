!! Copyright (C) 2019 N. Tancogne-Dejean
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

module states_abst_oct_m
  use global_oct_m
  use loct_oct_m
  use messages_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    states_abst_t,                    &
    states_are_complex,               &
    states_are_real,                  &
    states_set_complex

  type, abstract :: states_abst_t
    private
    type(type_t), public  :: wfs_type         !< real (TYPE_FLOAT) or complex (TYPE_CMPLX) wavefunctions
    integer, public  :: nst                   !< Number of states in each irreducible subspace
    logical, public  :: packed

  contains

    procedure(nullify),    deferred :: nullify
    procedure(pack),       deferred :: pack
    procedure(unpack),     deferred :: unpack  
    procedure(write_info), deferred :: write_info
    procedure(set_zero),   deferred :: set_zero
    procedure, non_overridable      :: are_packed
    procedure, non_overridable      :: get_type
  end type states_abst_t

  abstract interface
    subroutine nullify(st)
      import states_abst_t
      class(states_abst_t), intent(inout) :: st
    end subroutine nullify

    subroutine set_zero(st)
      import states_abst_t
      class(states_abst_t), intent(inout) :: st
    end subroutine set_zero

    subroutine write_info(st)
      import states_abst_t
      class(states_abst_t), intent(in) :: st
    end subroutine write_info

    subroutine pack(st, copy)
      import states_abst_t
      class(states_abst_t), intent(inout) :: st
      logical, optional,    intent(in)    :: copy 
    end subroutine pack

    subroutine unpack(st, copy)
      import states_abst_t
      class(states_abst_t), intent(inout) :: st
      logical, optional,    intent(in)    :: copy
    end subroutine unpack
  end interface

contains

  ! ---------------------------------------------------------
  subroutine states_set_complex(st)
    class(states_abst_t),    intent(inout) :: st

    PUSH_SUB(states_set_complex)

    st%wfs_type = TYPE_CMPLX

    POP_SUB(states_set_complex)
  end subroutine states_set_complex

  ! ---------------------------------------------------------
  pure logical function states_are_complex(st) result (wac)
    class(states_abst_t),    intent(in) :: st

    wac = (st%wfs_type == TYPE_CMPLX)

  end function states_are_complex


  ! ---------------------------------------------------------
  pure logical function states_are_real(st) result (war)
    class(states_abst_t),    intent(in) :: st

    war = (st%wfs_type == TYPE_FLOAT)

  end function states_are_real

  ! -----------------------------------------------------------

  logical pure function are_packed(st) result(packed)
    class(states_abst_t),    intent(in) :: st

    packed = st%packed
  end function are_packed

  pure type(type_t) function get_type(st) result(res)
    class(states_abst_t),    intent(in) :: st

    res = st%wfs_type
  end function get_type 


end module states_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
