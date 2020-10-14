!! Copyright (C) 2014 Vladimir Fuka
!! Copyright (C) 2020 Heiko Appel
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, version 3
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


module cgal_polyhedra_oct_m
  use iso_c_binding

  implicit none

  private

  public cgal_polyhedron_read, &
         cgal_polyhedron_point_inside, &
         cgal_polyhedron_finalize

  type, bind(C) :: d3
    real(c_double) :: x, y, z
  end type

  interface

    subroutine polyhedron_from_file(ptree, fname, verbose, ierr) bind(C,name="polyhedron_from_file")
      import
      type(c_ptr), intent(out) :: ptree
      character(kind=c_char), intent(in) :: fname(*)
      integer(c_int), value :: verbose                 ! avoid bool in C++, not equal to c_bool
      integer(c_int), intent(out) :: ierr
    end subroutine

    function polyhedron_point_inside(poly, vec) result(res) bind(C,name="polyhedron_point_inside")
      import
      logical(c_bool) :: res
      type(c_ptr), value :: poly
      type(d3),intent(in) :: vec
    end function

    subroutine polyhedron_finalize(ptree) bind(C,name="polyhedron_finalize")
      import
      type(c_ptr), intent(inout) :: ptree
    end subroutine

  end interface

  interface cgal_polyhedron_point_inside
    module procedure cgal_polyhedron_point_inside_s
    module procedure cgal_polyhedron_point_inside_d
  end interface


contains

  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_read(ptree, fname)
    type(c_ptr), intent(out) :: ptree
    character(*), intent(in) :: fname
    integer(c_int) :: ierr

    call polyhedron_from_file(ptree, fname//c_null_char, 1, ierr)

    if (ierr==1) then
      write (*,*) "Error reading file "//fname//", it appears to be empty."
      stop
    else if (ierr==2) then
      write (*,*) "Error reading file "//fname//"."
      stop
    end if
  end subroutine

  ! ---------------------------------------------------------
  function cgal_polyhedron_point_inside_s(poly, xq,yq,zq) result(res)
    logical :: res
    type(c_ptr), intent(in) :: poly
    real(c_float), intent(in) :: xq,yq,zq
    type(d3) :: query

    query = d3(real(xq,c_double),real(yq,c_double),real(zq,c_double))
    res = polyhedron_point_inside(poly, query)
  end function

  ! ---------------------------------------------------------
  function cgal_polyhedron_point_inside_d(poly, xq,yq,zq) result(res)
    logical :: res
    type(c_ptr), intent(in) :: poly
    real(c_double), intent(in) :: xq,yq,zq
    type(d3) :: query

    query = d3(xq,yq,zq)
    res = polyhedron_point_inside(poly, query)
  end function

  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_finalize(ptree)
    type(c_ptr), intent(inout) :: ptree

    call polyhedron_finalize(ptree)
  end subroutine

end module cgal_polyhedra_oct_m
