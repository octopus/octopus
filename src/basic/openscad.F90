!! Copyright (C) 2013 X. Andrade
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
!! $Id: openscad.F90 9479 2012-10-04 22:11:39Z dstrubbe $

#include "global.h"

module openscad_m
  use global_m
  use messages_m
  use io_m

  implicit none

  private

  public ::                & 
    openscad_file_t,       &
    openscad_file_init,    &
    openscad_file_end,     &
    openscad_file_sphere,  &
    openscad_file_bond

  type openscad_file_t
    private
    
    integer :: iunit
    FLOAT   :: scale
    integer :: curve_res
  end type openscad_file_t

contains

  subroutine openscad_file_init(this, filename, scale)
    type(openscad_file_t), intent(out) :: this
    character(len=*),      intent(in)  :: filename
    FLOAT,                 intent(in)  :: scale    
    
    PUSH_SUB(openscad_file_init)

    this%iunit = io_open(filename, action = 'write')
    this%scale = scale
    this%curve_res = 50

    write(this%iunit, '(a,i10,a)') '$fn = ', this%curve_res, ';'
    write(this%iunit, '(a,f12.6,a)') 'scale(', this%scale, '){'

    POP_SUB(openscad_file_init)
  end subroutine openscad_file_init

  !-------------------------------------------------------

  subroutine openscad_file_sphere(this, position, radius)
    type(openscad_file_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: position(:)
    FLOAT,                 intent(in)    :: radius

    PUSH_SUB(openscad_file_sphere)

    call write_translate(this, position)
    write(this%iunit, '(a,f12.6,a)') '  sphere(', radius, ');'

    POP_SUB(openscad_file_sphere)
  end subroutine openscad_file_sphere
  
  !-------------------------------------------------------

  subroutine openscad_file_bond(this, pos1, pos2, radius)
    type(openscad_file_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: pos1(:)
    FLOAT,                 intent(in)    :: pos2(:)
    FLOAT,                 intent(in)    :: radius

    FLOAT :: length, vec(1:3), angles(1:3)

    PUSH_SUB(openscad_file_bond)

    vec(1:3) = pos2(1:3) - pos1(1:3)
    length = sqrt(sum(vec(1:3)**2))
    
    angles(1) = -acos(vec(3)/length)
    angles(2) = CNST(0.0)
    angles(3) = -atan2(vec(1), vec(2))

    angles(1:3) = angles(1:3)*CNST(180.0)/M_PI
    
    call write_translate(this, CNST(0.5)*(pos1(1:3) + pos2(1:3)))
    call write_rotate(this, angles)
    write(this%iunit, '(a,f12.6,a,f12.6,a)') '  cylinder(h=', length, ', r =', radius, ', center = true);'

    POP_SUB(openscad_file_bond)
  end subroutine openscad_file_bond
  
  !-------------------------------------------------------

  subroutine openscad_file_end(this)
    type(openscad_file_t), intent(inout) :: this
    
    PUSH_SUB(openscad_file_end)

    write(this%iunit, '(a)') '}'

    call io_close(this%iunit)

    POP_SUB(openscad_file_end)
  end subroutine openscad_file_end
  
  !-------------------------------------------------------

  subroutine write_translate(this, position)
    type(openscad_file_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: position(:)

    write(this%iunit, '(a,f12.6,a,f12.6,a,f12.6,a)') 'translate([', position(1), ',',  position(2), ',',  position(3), '])'
    
  end subroutine write_translate

  !-------------------------------------------------------

  subroutine write_rotate(this, angles)
    type(openscad_file_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: angles(:)

    write(this%iunit, '(a,f12.6,a,f12.6,a,f12.6,a)') 'rotate([', angles(1), ',',  angles(2), ',',  angles(3), '])'
    
  end subroutine write_rotate
  
end module openscad_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
