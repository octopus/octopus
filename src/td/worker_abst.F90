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

module worker_abst_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                  &
    worker_abst_t

  type, abstract :: worker_abst_t
    private

  contains

    !Below are the list of operations that needs to be implemented by the workers
    procedure(worker_init), deferred    :: init
    procedure(worker_end),  deferred    :: end
  end type worker_abst_t


  abstract interface
    subroutine worker_init(wo)
      import worker_abst_t
      class(worker_abst_t), intent(inout) :: wo
    end subroutine worker_init

    subroutine worker_end(wo)
      import worker_abst_t
      class(worker_abst_t), intent(inout) :: wo
    end subroutine worker_end

  end interface 

contains

end module worker_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
