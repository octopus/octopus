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

module mesh_cube_map_oct_m
  use accel_oct_m
  use global_oct_m
  use index_oct_m
  use messages_oct_m
  use profiling_oct_m
  use types_oct_m

  implicit none

  private
  public ::                               &
    mesh_cube_map_t,                      &
    mesh_cube_map_init,                   &
    mesh_cube_map_end

  type mesh_cube_map_t
    ! Components are public by default
    integer            :: nmap      !< The number of maps
    integer, pointer   :: map(:, :)
    type(accel_mem_t) :: map_buffer
  end type mesh_cube_map_t

  integer, public, parameter :: MCM_POINT = 4, MCM_COUNT = 5

contains

  ! ---------------------------------------------------------

  subroutine mesh_cube_map_init(this, idx, np_global)
    type(mesh_cube_map_t),         intent(out) :: this
    type(index_t),                 intent(in)  :: idx
    integer,                       intent(in)  :: np_global

    integer :: i1(1:3), i2(1:3)
    integer :: step, ip

    PUSH_SUB(mesh_cube_map_init)

    if (idx%dim <= 3) then
      do step = 1, 2
        if(step == 2) then
          SAFE_ALLOCATE(this%map(1:5, 1:this%nmap))
        end if

        this%nmap = 0
        i2 = 0
        do ip = 1, np_global

          i1 = 0
          call index_to_coords(idx, ip, i1(1:3))

          if(any(i1(1:2) /= i2(1:2)) .or. i1(3) /= i2(3) + 1) then
            INCR(this%nmap, 1)
            if(step == 2) then
              call index_to_coords(idx, ip, this%map(1:, this%nmap))
              this%map(idx%dim + 1:3, this%nmap) = 0
              this%map(MCM_POINT, this%nmap) = ip
              this%map(MCM_COUNT, this%nmap) = 1
            end if
          else
            if(step == 2) INCR(this%map(MCM_COUNT, this%nmap), 1)
          end if
          i2 = i1
        end do
      end do

      if(accel_is_enabled()) then
        call accel_create_buffer(this%map_buffer, ACCEL_MEM_READ_ONLY, TYPE_INTEGER, this%nmap*5)
        call accel_write_buffer(this%map_buffer, this%nmap*5, this%map)
      end if

    else
      nullify(this%map)
    end if

    POP_SUB(mesh_cube_map_init)
  end subroutine mesh_cube_map_init

  ! ---------------------------------------------------------

  subroutine mesh_cube_map_end(this)
    type(mesh_cube_map_t), intent(inout) :: this

    PUSH_SUB(mesh_cube_map_end)

    if(associated(this%map)) then

      SAFE_DEALLOCATE_P(this%map)
      
      if(accel_is_enabled()) call accel_release_buffer(this%map_buffer)

    end if

    POP_SUB(mesh_cube_map_end)
  end subroutine mesh_cube_map_end

  ! ---------------------------------------------------------
  
end module mesh_cube_map_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
