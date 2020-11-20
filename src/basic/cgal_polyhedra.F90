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

#include "global.h"

module cgal_polyhedra_oct_m
  use global_oct_m
  use iso_c_binding
  use messages_oct_m

  implicit none

  private

  public ::                          &
    cgal_polyhedra_t,                &
    cgal_polyhedron_init,            &
    cgal_polyhedron_read,            &
    cgal_polyhedron_build_AABB_tree, &
    cgal_polyhedron_point_inside,    &
    cgal_polyhedron_end

  type, bind(C) :: d3
    real(c_double) :: x, y, z
  end type

  type cgal_polyhedra_t
    private
    type(c_ptr) :: AABB_tree
    type(c_ptr) :: polyhedron
  end type cgal_polyhedra_t

#ifdef HAVE_CGAL
  interface

    subroutine polyhedron_from_file(ptree, fname, verbose, ierr) bind(C,name="polyhedron_from_file")
      import
      type(c_ptr), intent(out)            :: ptree
      character(kind=c_char), intent(in)  :: fname(*)
      integer(c_int), value               :: verbose    ! avoid bool in C++, not equal to c_bool
      integer(c_int), intent(out)         :: ierr
    end subroutine

    subroutine polyhedron_build_AABB_tree(tree, poly) bind(C,name="polyhedron_build_AABB_tree")
      import
      type(c_ptr), intent(out) :: tree
      type(c_ptr), intent(in)  :: poly
    end subroutine

    function polyhedron_point_inside(tree, vec) result(res) bind(C,name="polyhedron_point_inside")
      import
      logical(c_bool)         :: res
      type(c_ptr), intent(in) :: tree
      type(d3),    intent(in) :: vec
    end function

    subroutine polyhedron_finalize_AABB_tree(tree) bind(C,name="polyhedron_finalize_AABB_tree")
      import
      type(c_ptr), intent(inout) :: tree
    end subroutine

    subroutine polyhedron_finalize_polyhedron(polyhedron) bind(C,name="polyhedron_finalize_polyhedron")
      import
      type(c_ptr), intent(inout) :: polyhedron
    end subroutine

  end interface
#endif


contains

  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_init(cgal_poly, fname, verbose)
    type(cgal_polyhedra_t), intent(inout) :: cgal_poly
    character(*),              intent(in) :: fname
    logical,                   intent(in) :: verbose

    PUSH_SUB(cgal_polyhedron_init)

    call cgal_polyhedron_read(cgal_poly, fname, verbose)
    call cgal_polyhedron_build_AABB_tree(cgal_poly)

    POP_SUB(cgal_polyhedron_init)
  end subroutine cgal_polyhedron_init


  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_read(cgal_poly, fname, verbose)
    type(cgal_polyhedra_t), intent(inout) :: cgal_poly
    character(*),              intent(in) :: fname
    logical,                   intent(in) :: verbose

    integer(c_int) :: verb, ierr

    PUSH_SUB(cgal_polyhedron_read)

    verb = 0
    ierr = 0

    if (verbose) verb = 1
#ifdef HAVE_CGAL
    call polyhedron_from_file(cgal_poly%polyhedron, fname//c_null_char, verb, ierr)
#endif
#ifndef HAVE_CGAL
    ierr = 3
#endif

    select case(ierr)
    case(0)
      message(1) = "Info: finished reading polyhedron from file " // fname
      call messages_info(1)
    case(1)
      message(1) = "Error reading file " // fname // ", it appears to be empty."
      call messages_fatal(1)
    case(2)
      message(1) = "Error reading file " // fname // "."
      call messages_fatal(1)
    case(3)
      message(1) = "You are trying to read polyhedron data from " // fname
      message(2) = "For this feature Octopus has to be linked with the CGAL library."
      call messages_fatal(2)
    case default
      message(1) = "Error: Error status not implemented in CGAL Fortran interface."
      call messages_fatal(1)
    end select

    POP_SUB(cgal_polyhedron_read)
  end subroutine cgal_polyhedron_read


  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_build_AABB_tree(cgal_poly)
    type(cgal_polyhedra_t), intent(inout) :: cgal_poly

    PUSH_SUB(cgal_polyhedron_build_AABB_tree)

#ifdef HAVE_CGAL
    call polyhedron_build_AABB_tree(cgal_poly%AABB_tree, cgal_poly%polyhedron)
#endif

    POP_SUB(cgal_polyhedron_build_AABB_tree)
  end subroutine cgal_polyhedron_build_AABB_tree


  ! ---------------------------------------------------------
  function cgal_polyhedron_point_inside(cgal_poly, xq, yq, zq) result(res)
    logical                            :: res
    type(cgal_polyhedra_t), intent(in) :: cgal_poly
    real(c_double),         intent(in) :: xq, yq, zq

    type(d3) :: query

    ! no push/pop here since this is called too frequently
    res = .false.
    query = d3(xq, yq, zq)
#ifdef HAVE_CGAL
    res = polyhedron_point_inside(cgal_poly%AABB_tree, query)
#endif

  end function cgal_polyhedron_point_inside

  ! ---------------------------------------------------------
  subroutine cgal_polyhedron_end(cgal_poly)
    type(cgal_polyhedra_t), intent(inout) :: cgal_poly

    PUSH_SUB(cgal_polyhedron_end)

#ifdef HAVE_CGAL
    call polyhedron_finalize_AABB_tree(cgal_poly%AABB_tree)
    call polyhedron_finalize_polyhedron(cgal_poly%polyhedron)
#endif

    POP_SUB(cgal_polyhedron_end)
  end subroutine cgal_polyhedron_end

end module cgal_polyhedra_oct_m
