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
!! $Id: species.F90 13307 2015-03-09 06:55:02Z xavier $

#include "global.h"

module volume_oct_m
  use iso_c_binding
  use parser_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m

  implicit none

  private
  public :: &
       volume_t,                &
       volume_init,             &
       volume_end,              &
       volume_read_from_block,  &
       volume_in_volume

  type volume_t
    integer          :: n_elements
    logical, pointer :: join(:)      => NULL()         ! Add or subtract the volume
    integer, pointer :: type(:)      => NULL()         ! sphere, slab, etc.
    FLOAT,   pointer :: params(:,:)  => NULL()         ! parameters of the elements
  end type volume_t

contains

  subroutine volume_init(vol)
    type(volume_t),   intent(out) :: vol

    vol%n_elements = 0
  end subroutine volume_init

  subroutine volume_end(vol)
    type(volume_t),   intent(inout) :: vol

    SAFE_DEALLOCATE_P(vol%join)
    SAFE_DEALLOCATE_P(vol%type)
    SAFE_DEALLOCATE_P(vol%params)
  end subroutine volume_end

  subroutine volume_read_from_block(vol, block_name)
    type(volume_t),   intent(inout) :: vol
    character(len=*), intent(in)    :: block_name

    type(block_t) :: blk
    integer :: i, j, n_par
    character(len=100) :: str

    !%Variable Volume
    !%Type block
    !%Default 
    !%Section Utilities::
    !%Description
    !% Describes a volume in space defined through the addition and substraction of
    !% spheres. The first field is always "+" (include points inside the volume) or "-"
    !% (exclude points inside the volume)
    !%Option vol_sphere 10001
    !%
    !% <tt>%Volume
    !% <br>&nbsp;&nbsp; "+"/"-" | vol_sphere | center_x | center_y | center_z | radius 
    !% <br>%</tt>
    !%Option vol_slab 10002
    !%
    !% <tt>%Volume
    !% <br>&nbsp;&nbsp; "+"/"-" | vol_slab | thickness
    !% <br>%</tt>
    !%
    !%End

    if(parse_block(block_name, blk, check_varinfo_=.false.) == 0) then
      vol%n_elements = parse_block_n(blk)

      SAFE_ALLOCATE(vol%join(1:vol%n_elements))
      SAFE_ALLOCATE(vol%type(1:vol%n_elements))
      SAFE_ALLOCATE(vol%params(1:8, 1:vol%n_elements))

      vol%params = M_ZERO

      do i = 1, vol%n_elements
        call parse_block_string(blk, i-1, 0, str)
        if(str == '+') then
          vol%join(i) = .true.
        else
          vol%join(i) = .false.
        end if

        call parse_block_integer(blk, i-1, 1, vol%type(i))
        select case(vol%type(i))
        case(OPTION__VOLUME__VOL_SPHERE)
          n_par = 4 ! center point + radius
        case(OPTION__VOLUME__VOL_SLAB)
          n_par = 1 ! thickness of the slab
        case default
          call messages_input_error('Species', "Unknown type for volume")
        end select

        do j = 1, n_par
          call parse_block_float(blk, i-1, i+j, vol%params(j, i))
        end do

      end do
    else
      call messages_input_error('Volume')
    end if
  end subroutine volume_read_from_block


  logical function volume_in_volume(sb, vol, xx, rr) result(in_vol)
    type(simul_box_t), intent(in) :: sb
    type(volume_t),    intent(in) :: vol
    FLOAT,             intent(in) :: xx(:), rr

    logical :: in_partial_volume
    integer :: i
    FLOAT   :: r, x(MAX_DIM)

    in_vol = .false.
    do i = 1, vol%n_elements
      select case(vol%type(i))
      case(OPTION__VOLUME__VOL_SPHERE)
        x(1:sb%dim) = xx(1:sb%dim) - vol%params(1:sb%dim, i)
        r = sqrt(dot_product(x(1:sb%dim), x(1:sb%dim)))
        in_partial_volume = (r <= vol%params(4, i))

      case(OPTION__VOLUME__VOL_SLAB)
        r = abs(xx(3))
        in_partial_volume = (r <= vol%params(1, i))
      end select

      if(vol%join(i)) then
        in_vol = in_vol .or. in_partial_volume
      else
        in_vol = in_vol .and. .not. in_partial_volume
      end if
    end do

  end function volume_in_volume

end module volume_oct_m