!! Copyright (C) 2015 X. Andrade
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
!! $Id$

#include "global.h"

module sphericalmesh_oct_m
  use global_oct_m
  use io_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_interpolation_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  
  implicit none
  
  private
  public ::                        &
    sphericalmesh_t,               &
    sphericalmesh_init,            &
    sphericalmesh_end,             &
    sphericalmesh_integrate,       &
    sphericalmesh_from_mesh
  
  type sphericalmesh_t
    integer            :: np
    FLOAT, allocatable :: x(:, :)
    FLOAT, allocatable :: weight(:)
    integer            :: nangular
    integer            :: nradial
  end type sphericalmesh_t
    
contains

  subroutine sphericalmesh_init(this, center, radius)
    type(sphericalmesh_t), intent(out)  :: this
    FLOAT,              intent(in)   :: center(:)
    FLOAT,              intent(in)   :: radius

    integer, parameter :: nlebedev = 32
    integer :: lebedev(1:2, 1:nlebedev)
    integer :: precision
    integer :: ii, iunit, ir, ia
    FLOAT, allocatable :: angular(:, :), radial(:, :)
    FLOAT :: rr, theta, phi
    character(len=MAX_PATH_LEN) :: filename
    
    PUSH_SUB(sphericalmesh_init)

    precision = 17
    
    ! the information of the radial quadrature, first number is
    ! precision, second number is number of points
    lebedev = reshape((/&
      3,   6,    &
      5,   14,   &
      7,   26,   &
      9,   38,   &
      11,  50,   &
      13,  74,   &
      15,  86,   &
      17,  110,  &
      19,  146,  &
      21,  170,  &
      23,  194,  &
      25,  230,  &
      27,  266,  &
      29,  302,  &
      31,  350,  &
      35,  434,  &
      41,  590,  &
      47,  770,  &
      53,  974,  &
      59,  1202, &
      65,  1454, &
      71,  1730, &
      77,  2030, &
      83,  2354, &
      89,  2702, &
      95,  3074, &
      101, 3470, &
      107, 3890, &
      113, 4334, &
      119, 4802, &
      125, 5294, &
      131, 5810  &
      /), shape(lebedev))

    do ii = 1, nlebedev
      if(lebedev(1, ii) >= precision) then
        precision = lebedev(1, ii)
        this%nangular = lebedev(2, ii)
        exit
      end if
    end do

    SAFE_ALLOCATE(angular(1:3, this%nangular))

    write(filename, '(2a,i3.3,a)') trim(conf%share), '/quadrature/lebedev_', precision, '.txt'
    
    iunit = io_open(trim(filename), action = 'read', status = 'old', die = .true.)

    do ii = 1, this%nangular
      read(iunit, *) angular(1, ii), angular(2, ii), angular(3, ii)
      angular(1:2, ii) = angular(1:2, ii)*2*M_PI/CNST(360.0)
    end do

    call io_close(iunit)

    this%nradial = 20
    
    SAFE_ALLOCATE(radial(1:2, this%nradial))

    call radial_init()


    this%np = this%nradial*this%nangular

    print*, "NP", this%np
    
    SAFE_ALLOCATE(this%x(1:3, 1:this%np))
    SAFE_ALLOCATE(this%weight(1:this%np))

    ii = 1
    do ia = 1, this%nangular
      phi = angular(1, ia)
      theta = angular(2, ia)
      do ir = 1, this%nradial
        rr = radial(1, ir)
        this%x(1:3, ii) = center(1:3) + (/rr*sin(theta)*cos(phi), rr*sin(theta)*sin(phi), rr*cos(theta)/)
        this%weight(ii) = CNST(4.0)*M_PI*radial(2, ir)*angular(3, ia)
        ii = ii + 1
      end do
    end do

    do ii = 1, this%np
      write(12, *) this%x(1:3, ii)
    end do
    
    SAFE_DEALLOCATE_A(angular)
    SAFE_DEALLOCATE_A(radial)
    
    POP_SUB(sphericalmesh_init)

  contains

    subroutine radial_init()
      integer :: ir
      FLOAT   :: factor, xi, dr

#if 0
      
      ! this is the radial grid described in:
      ! Becke and Dickson JCP 89 2993-2997 (1988) http://dx.doi.org/10.1063/1.455005
      
      factor = M_PI/(this%nradial + CNST(1.0))

      do ir = 1, this%nradial
        xi = cos(factor*ir)
        radial(1, ir) = radius*(CNST(1.0) + xi)/(CNST(1.0) - xi)
        radial(2, ir) = radius*CNST(2.0)*factor*sin(factor*ir)/(xi - CNST(1.0))**2*radial(1, ir)**2
        !        print*, radial(:, ir)
      end do

#else
      ! this is a uniform grid
      
      dr = radius/this%nradial
      do ir = 1, this%nradial
        radial(1, ir) = dr*ir
        radial(2, ir) = dr !*radial(1, ir)**2
        print*, radial(:, ir)
      end do

      radial(2, 1) = CNST(1.5)*radial(2, 1)
      radial(2, this%nradial) = CNST(0.5)*radial(2, this%nradial)

      print*, "RADIAL SUM", sum(radial(2, :))
      
      do ir = 1, this%nradial
        print*, radial(:, ir)
      end do
#endif
      
    end subroutine radial_init
    
  end subroutine sphericalmesh_init

  ! -------------------------------------------------------------
  
  subroutine sphericalmesh_end(this)
    type(sphericalmesh_t), intent(inout) :: this
    
    PUSH_SUB(sphericalmesh_end)

    SAFE_DEALLOCATE_A(this%x)
    SAFE_DEALLOCATE_A(this%weight)
    
    POP_SUB(sphericalmesh_end)
  end subroutine sphericalmesh_end

  ! -------------------------------------------------------------

  FLOAT function sphericalmesh_integrate(this, ff) result(integral)
    type(sphericalmesh_t), intent(inout) :: this
    FLOAT,                 intent(in)    :: ff(:)

    PUSH_SUB(sphericalmesh_integrate)
    
    integral = sum(this%weight(1:this%np)*ff(1:this%np))

    POP_SUB(sphericalmesh_integrate)
  end function sphericalmesh_integrate

  ! -------------------------------------------------------------

  subroutine sphericalmesh_from_mesh(this, mesh, meshff, sphff)
    type(sphericalmesh_t), intent(in)    :: this
    type(mesh_t),          intent(in)    :: mesh
    FLOAT,                 intent(in)    :: meshff(:)
    FLOAT,                 intent(inout) :: sphff(:)

    type(profile_t), save :: prof
    type(mesh_interpolation_t) :: interp
    
    PUSH_SUB(sphericalmesh_from_mesh)

    call profiling_in(prof, 'TO_SPHERICAL')
    
    call mesh_interpolation_init(interp, mesh)

    call mesh_interpolation_evaluate(interp, this%np, meshff, this%x, sphff)
    
    call mesh_interpolation_end(interp)

    call profiling_out(prof)
    
    POP_SUB(sphericalmesh_from_mesh)
  end subroutine sphericalmesh_from_mesh
  
  
end module sphericalmesh_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
