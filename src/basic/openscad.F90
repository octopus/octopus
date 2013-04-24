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

  public ::                         &
    openscad_file_t,                &
    openscad_file_init,             &
    openscad_file_end,              &
    openscad_file_sphere,           &
    openscad_file_bond,             &
    openscad_file_cube,             &
    openscad_file_define_variable,  &
    openscad_file_triangle

  type openscad_file_t
    private
    
    integer :: iunit
    integer :: curve_res
  end type openscad_file_t

contains

  subroutine openscad_file_init(this, filename)
    type(openscad_file_t), intent(out) :: this
    character(len=*),      intent(in)  :: filename
    
    PUSH_SUB(openscad_file_init)

    this%iunit = io_open(filename, action = 'write')
    this%curve_res = 50

    write(this%iunit, '(a,i10,a)') '$fn = ', this%curve_res, ';'

    POP_SUB(openscad_file_init)
  end subroutine openscad_file_init

  !-------------------------------------------------------

  subroutine openscad_file_define_variable(this, variable_name, variable_value)
    type(openscad_file_t), intent(inout) :: this
    character(len=*),      intent(in)    :: variable_name
    FLOAT,                 intent(in)    :: variable_value

    PUSH_SUB(openscad_file_define_variable)

    write(this%iunit, '(a,f12.6,a)') variable_name//' = ', variable_value, ';'

    POP_SUB(openscad_file_define_variable)
  end subroutine openscad_file_define_variable
  
  !-------------------------------------------------------

  subroutine openscad_file_sphere(this, position, radius, radius_variable)
    type(openscad_file_t),      intent(inout) :: this
    FLOAT,                      intent(in)    :: position(:)
    FLOAT,            optional, intent(in)    :: radius
    character(len=*), optional, intent(in)    :: radius_variable

    PUSH_SUB(openscad_file_sphere)

    ASSERT(.not. present(radius) .eqv. present(radius_variable))

    call write_translate(this, position)

    if(present(radius)) then
      write(this%iunit, '(a,f12.6,a)') '  sphere(', radius, ');'
    else
      write(this%iunit, '(a,a,a)') '  sphere(', trim(radius_variable), ');'
    end if

    POP_SUB(openscad_file_sphere)
  end subroutine openscad_file_sphere

  !-------------------------------------------------------

  subroutine openscad_file_cube(this, position, sizes)
    type(openscad_file_t),      intent(inout) :: this
    FLOAT,                      intent(in)    :: position(:)
    FLOAT,            optional, intent(in)    :: sizes(:)

    PUSH_SUB(openscad_file_cube)
    
    call write_translate(this, position)
    write(this%iunit, '(a,f12.6,a,f12.6,a,f12.6,a)') '  cube([', sizes(1), ',', sizes(2), ',', sizes(3), '], center = true);'

    POP_SUB(openscad_file_cube)
  end subroutine openscad_file_cube

  !-------------------------------------------------------

  subroutine openscad_file_bond(this, pos1, pos2, radius, radius_variable)
    type(openscad_file_t),      intent(inout) :: this
    FLOAT,                      intent(in)    :: pos1(:)
    FLOAT,                      intent(in)    :: pos2(:)
    FLOAT,            optional, intent(in)    :: radius
    character(len=*), optional, intent(in)    :: radius_variable

    FLOAT :: length, vec(1:3), angles(1:3)

    PUSH_SUB(openscad_file_bond)

    ASSERT(.not. present(radius) .eqv. present(radius_variable))

    vec(1:3) = pos2(1:3) - pos1(1:3)
    length = sqrt(sum(vec(1:3)**2))
    
    angles(1) = -acos(vec(3)/length)
    angles(2) = CNST(0.0)
    angles(3) = -atan2(vec(1), vec(2))

    angles(1:3) = angles(1:3)*CNST(180.0)/M_PI
    
    call write_translate(this, CNST(0.5)*(pos1(1:3) + pos2(1:3)))
    call write_rotate(this, angles)

    if(present(radius)) then
      write(this%iunit, '(a,f12.6,a,f12.6,a)') '  cylinder(h=', length, ', r =', radius, ', center = true);'
    else
      write(this%iunit, '(a,f12.6,a,a,a)') '  cylinder(h=', length, ', r =', trim(radius_variable), ', center = true);'
    end if

    POP_SUB(openscad_file_bond)
  end subroutine openscad_file_bond
  
  !-------------------------------------------------------

  subroutine openscad_file_end(this)
    type(openscad_file_t), intent(inout) :: this
    
    PUSH_SUB(openscad_file_end)

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
  
  !------------------------------------------------------

  subroutine openscad_file_triangle(this, pos1, pos2, pos3)
    type(openscad_file_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: pos1(:)
    FLOAT,                 intent(in)    :: pos2(:)
    FLOAT,                 intent(in)    :: pos3(:)

    write(this%iunit, '(a,a,f12.6,a,f12.6,a,f12.6,a,a,f12.6,a,f12.6,a,f12.6,a,a,f12.6,a,f12.6,a,f12.6,a,a)') &
      ' polyhedron( points = [ ', &
      '[', pos1(1), ',', pos1(2), ',', pos1(3), '],',       &
      '[', pos2(1), ',', pos2(2), ',', pos2(3), '],',       &
      '[', pos3(1), ',', pos3(2), ',', pos3(3), '],',       &
      '], triangles=[ [2,1,0] ]);'
    
  end subroutine openscad_file_triangle

end module openscad_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
