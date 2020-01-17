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
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module propagation_ops_abst_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                  &
    propagation_ops_abst_t

  type, abstract :: propagation_ops_abst_t
    private

  contains

    !Below are the list of operations that needs to be implemented by the propagation_opss
    procedure(propagation_ops_init), deferred    :: init
    procedure(propagation_ops_end),  deferred    :: end
  end type propagation_ops_abst_t


  abstract interface
    subroutine propagation_ops_init(wo)
      import propagation_ops_abst_t
      class(propagation_ops_abst_t), intent(inout) :: wo
    end subroutine propagation_ops_init

    subroutine propagation_ops_end(wo)
      import propagation_ops_abst_t
      class(propagation_ops_abst_t), intent(inout) :: wo
    end subroutine propagation_ops_end

  end interface

contains

end module propagation_ops_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
