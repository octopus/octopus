!! Copyright (C) 2005-2006 Heiko Appel
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module scalar_mesh_m
  use global_m
  use messages_m
  use datasets_m
  use math_m
  use io_m
  use lib_oct_parser_m

  implicit none

  private
  public ::                                   &
    scalar_mesh_t,                            &
    scalar_mesh_init,                         &
    scalar_mesh_create,                       &
    scalar_mesh_end,                          &
    scalar_mesh_integrate,                    &
    scalar_mesh_write

  integer, parameter ::                       &
    MESH_LINEAR         = 1,                  &
    MESH_LOG            = 2,                  &
    MESH_DOUBLE_LOG     = 3,                  &
    MESH_SINH           = 4,                  &
    MESH_GAUSS_LEGENDRE = 5,                  &
    MESH_MINVAL         = MESH_LINEAR,        &
    MESH_MAXVAL         = MESH_GAUSS_LEGENDRE

  type scalar_mesh_t
    integer :: mtype                ! mesh type (see MESH_* variables above)
    integer :: np                   ! number of points in the mesh
    FLOAT   :: min, max             ! lower and upper boundary for scalar grid
    FLOAT   :: center               ! highest density of gridpoints is around center
    FLOAT   :: alpha1 = CNST(0.3)   ! parameters for logarithmic grid below(1) and above(2) center
    FLOAT   :: alpha2 = CNST(0.3)
    FLOAT   :: dx                   ! grid spacing for linear grid
    FLOAT, pointer    :: mesh(:)    ! the mesh points
    FLOAT, pointer    :: weights(:) ! weights for integrations
    character(len=32) :: label      ! label for namespacing
  end type scalar_mesh_t

contains


  ! ---------------------------------------------------------
  subroutine scalar_mesh_init(sm, label)
    type(scalar_mesh_t), intent(inout) :: sm
    character(len=*),    intent(in)    :: label

    FLOAT :: etmp

    call push_sub('scalar_mesh.scalar_mesh_init')
    !%Variable ScalarMeshType
    !%Type integer
    !%Default mesh_sinh
    !%Section Math::General
    !%Description
    !% Specifies what kind of scalar mesh will be used
    !%Option mesh_linear 1
    !% Linear mesh
    !%Option mesh_double_log 2
    !% Logarithmic mesh
    !%Option mesh_log 3
    !% Double logarithmic mesh
    !%Option mesh_sinh 4
    !% Sinh mesh
    !%Option gauss_legendre 5
    !% Gauss-Legendre mesh
    !%End
    call loct_parse_int  (check_inp(trim(label)//'MeshType'),             4, sm%mtype)
    if( sm%mtype.lt.MESH_MINVAL.or.sm%mtype.gt.MESH_MAXVAL ) then
      call input_error(check_inp(trim(label)//'MeshType'))
    end if
    call loct_parse_int  (check_inp(trim(label)//'MeshNPoints'),         20, sm%np)
    call loct_parse_float(check_inp(trim(label)//'MeshMin'),      CNST(0.0), sm%min)
    call loct_parse_float(check_inp(trim(label)//'MeshMax'),      CNST(4.0), sm%max)
    call loct_parse_float(check_inp(trim(label)//'MeshCenter'),   CNST(2.0), sm%center)
    call loct_parse_float(check_inp(trim(label)//'MeshAlpha1'),   CNST(0.3), sm%alpha1)
    call loct_parse_float(check_inp(trim(label)//'MeshAlpha2'),   CNST(0.3), sm%alpha2)

    ! a few sanity checks
    if(sm%min.ge.sm%max) then
      message(1) = 'Warning: scalar_mesh: '//trim(label)//'MeshMin is equal to or larger'
      message(2) = 'than '//trim(label)//'MeshMax.'
      message(3) = 'Reverting order.'
      call write_warning(3)
      etmp   = sm%min
      sm%min = sm%max
      sm%max = etmp
    end if

    if(sm%center.lt.sm%min.or.sm%center.gt.sm%max) then
      message(1) = 'Warning: scalar_mesh: '//trim(label)//'MeshCenter is not in interval'
      message(2) = '['//trim(label)//'MeshMin,'//trim(label)//'MeshMax].'
      message(3) = 'Will use the average of '//trim(label)//'MeshMin and '//trim(label)//'MeshMax.'
      call write_warning(3)
      sm%center = M_HALF*(sm%min+sm%max)
    end if

    sm%label = label

    call scalar_mesh_create(sm, label)

    call pop_sub()
  end subroutine scalar_mesh_init


  ! ---------------------------------------------------------
  subroutine scalar_mesh_create(sm, label)
    type(scalar_mesh_t), intent(inout) :: sm
    character(len=*),    intent(in)    :: label

    FLOAT :: xmin1, xmax1, xmin2, xmax2
    FLOAT, allocatable :: gl(:), gw(:)
    character(len=128) :: filename
    integer :: i

    call push_sub('scalar_mesh.scalar_mesh_create')

    xmin1    = M_ZERO
    xmax1    = sm%center - sm%min
    xmin2    = M_ZERO
    xmax2    = sm%max - sm%center
    sm%label = label
    sm%dx    = (sm%max-sm%min)/sm%np

    ! all meshes, except the linear one, need an odd number of mesh 
    ! points. Stop the code, if the user is giving wrong input.
    if(sm%mtype.ne.MESH_LINEAR) then
      if(modulo(sm%np, 2) .eq. 0) then
        message(1) = 'scalar_mesh: Need an odd number of points for "'//trim(sm%label)//'" mesh'
        message(2) = 'Please correct your input file.'
        call write_fatal(2)
      end if
    end if

    ALLOCATE(sm%mesh   (sm%np), sm%np)
    ALLOCATE(sm%weights(sm%np), sm%np)

    ! if we have less that 11 points we change the requested mesh type
    ! to a linear mesh and issue a warning
    if (sm%np.le.10.and.sm%mtype.ne.MESH_LINEAR) then
      sm%mtype = MESH_LINEAR
      message(1) = 'Less than 11 points in scalar mesh "'//trim(sm%label)//'".'
      message(2) = 'Switching to linear scalar mesh.'
      call write_warning(2)
    end if

    select case(sm%mtype)
    case(MESH_LINEAR)
      message(1) = 'Info: scalar_mesh: Using linear scalar mesh for '//trim(sm%label)
      call write_info(1)

      ! setup linear mesh
      if (sm%np.eq.1) then
        sm%mesh = M_HALF*(sm%max+sm%min)
      else
        do i = 1, sm%np
          sm%mesh(i) = sm%min + (sm%max-sm%min)/(sm%np-1)*(i-1)
        end do
      end if
      ! setup integration weights
      if (sm%np.lt.10) then
        ! manual integration
        sm%weights = sm%dx
      else
        ! Euler MacLaurin integration
        sm%weights = (sm%max-sm%min)/(sm%np-1)
        do i = 1, 5
          sm%weights(i)           = EMcLCoeff(i)*(sm%max-sm%min)/(sm%np-1)
          sm%weights(sm%np+1 - i) = EMcLCoeff(i)*(sm%max-sm%min)/(sm%np-1)
        end do
      end if


    case(MESH_LOG)
      message(1) = 'Info: scalar_mesh: Using logarithmic scalar mesh for '//trim(sm%label)
      call write_info(1)

      do i = 1, sm%np
        sm%mesh(i)      =   sm%min   + (sm%max - sm%min)*                       &
          ( exp(-sm%alpha1*(i-1)) - M_ONE )/( exp(-sm%alpha1*(sm%np-1))-M_ONE )

        sm%weights(i)   = - sm%alpha1 * (sm%max - sm%min)*                      &
          ( exp(-sm%alpha1*(i-1))         )/( exp(-sm%alpha1*(sm%np-1))-M_ONE )
      end do
      do i = 1, 5
        sm%weights(i)           = EMcLCoeff(i)*sm%weights(i)
        sm%weights(sm%np+1 - i) = EMcLCoeff(i)*sm%weights(sm%np+1 - i)
      end do


    case(MESH_DOUBLE_LOG)
      message(1) = 'Info: scalar_mesh: Using double logarithmic scalar mesh for '//trim(sm%label)
      call write_info(1)

      ! setup double logarithmic mesh
      sm%mesh = sm%center
      do i = 1, (sm%np-1)/2
        sm%mesh((sm%np-1)/2-i+1) = - ( sm%center + xmin1 + (xmax1-xmin1) * &
          ( exp(sm%alpha1*(i)) - M_ONE ) /                                 &
          ( exp(sm%alpha1*((sm%np-1)/2)) - M_ONE ) ) + 2*sm%center

        sm%mesh((sm%np-1)/2+i+1) =     sm%center + xmin2 + (xmax2-xmin2) * &
          ( exp(sm%alpha2*(i)) - M_ONE ) /                                 &
          ( exp(sm%alpha2*((sm%np-1)/2)) - M_ONE )
      end do
      ! setup integration weights
      do i = 1, (sm%np-1)/2
        sm%weights(i)             = sm%alpha1 * (xmax1-xmin1) /            &
          ( exp(sm%alpha1*((sm%np-1)/2))     - M_ONE ) *                   &
          ( exp(sm%alpha1*((sm%np-1)/2-i+1))         )
        sm%weights(sm%np+1 - i)   = sm%alpha2 * (xmax2-xmin2) /            &
          ( exp(sm%alpha2*((sm%np-1)/2))     - M_ONE ) *                   &
          ( exp(sm%alpha2*((sm%np-1)/2-i+1))         )
      end do
      sm%weights((sm%np-1)/2+1)        = EMcLCoeff(1)*                           &
        ( sm%alpha1 * (xmax1-xmin1) / ( exp(sm%alpha1*((sm%np-1)/2)) - M_ONE )   &
        + sm%alpha2 * (xmax2-xmin2) / ( exp(sm%alpha2*((sm%np-1)/2)) - M_ONE ) )
      sm%weights(1)                    = EMcLCoeff(1)*sm%weights(1)
      sm%weights(sm%np)                = EMcLCoeff(1)*sm%weights(sm%np)
      do i = 2, 5
        sm%weights(i)                  = EMcLCoeff(i)*sm%weights(i)
        sm%weights((sm%np-1)/2+2 - i)  = EMcLCoeff(i)*sm%weights((sm%np-1)/2+2 - i)
        sm%weights((sm%np-1)/2+i)      = EMcLCoeff(i)*sm%weights((sm%np-1)/2+i)
        sm%weights(sm%np+1 - i)        = EMcLCoeff(i)*sm%weights(sm%np+1 - i)
      end do


    case(MESH_SINH)
      message(1) = 'Info: scalar_mesh: Using sinh scalar mesh for '//trim(sm%label)
      call write_info(1)

      ! setup sinh mesh
      do i = 1, sm%np
        sm%mesh(i) = M_HALF*(sm%max-sm%min) + M_HALF*(sm%max-sm%min) *         &
          sinh(sm%alpha1*(i-(sm%np-1)/2-1))/sinh(sm%alpha1*((sm%np-1)/2))
      end do
      do i = 1, sm%np
        sm%weights(i)   =  sm%alpha1 * M_HALF*(sm%max-sm%min) /   &
          sinh(sm%alpha1*((sm%np-1)/2)) * cosh(sm%alpha1*((sm%np-1)/2-i+1))
      end do
      do i = 1, 5
        sm%weights(i)           = EMcLCoeff(i)*sm%weights(i)
        sm%weights(sm%np+1 - i) = EMcLCoeff(i)*sm%weights(sm%np+1 - i)
      end do


    case(MESH_GAUSS_LEGENDRE)
      message(1) = 'Info: scalar_mesh: Using Gauss-Legendre scalar mesh for '//trim(sm%label)
      call write_info(1)

      ALLOCATE(gl(sm%np), sm%np)
      ALLOCATE(gw(sm%np), sm%np)

      ! TODO: need to compute roots of Legendre Polynomials and weights for integration
      ! a call of the corresponding routine in math.F90 should be placed here

      !do i=1, 2*sm%np+1
      !   sm%mesh(i)    = M_HALF*(sm%max-sm%min) * gl(i) + M_HALF*(sm%max+sm%min)
      !   sm%weights(i) = gw(i)*M_HALF*(sm%max-sm%min)
      !end do
      deallocate(gl, gw)

      message(1) = 'scalar_mesh: Gauss-Legendre scalar mesh not fully implemented'
      call write_fatal(1)


    case default
      write(message(1), '(a,i4,a)') "Input: '", sm%mtype, &
        "' is not a valid scalar mesh type"
      message(2) = '( ScalarMeshType =  Mesh_Linear | Mesh_Log | Mesh_Double_Log | Mesh_Sinh )'
      call write_fatal(2)
    end select

    filename = trim(sm%label)//'_mesh.dat'
    call scalar_mesh_write(sm, filename)

    call pop_sub()
  end subroutine scalar_mesh_create


  ! ---------------------------------------------------------
  subroutine scalar_mesh_write(sm, filename)
    type(scalar_mesh_t), intent(in) :: sm
    character(len=*)                :: filename

    integer :: i, iunit

    call push_sub('scalar_mesh.scalar_mesh_write')
    ! output scalar mesh
    iunit = io_open('static/'//trim(filename), action='write')
    write(iunit,'(a)') '# index    gridpoint             integration weight'
    do i = 1, sm%np
      write(iunit,'(1x,i4,2x,2F22.16)') i, sm%mesh(i), sm%weights(i)
    end do
    call io_close(iunit)

    call pop_sub()
  end subroutine scalar_mesh_write


  ! ---------------------------------------------------------
  FLOAT function scalar_mesh_integrate(sm, func) result(res)
    type(scalar_mesh_t), intent(in) :: sm
    FLOAT, intent(in) :: func(:)

    call push_sub('scalar_mesh.scalar_mesh_integrate')

    res = sum(sm%weights(:)*func(:))

    call pop_sub()
  end function scalar_mesh_integrate


  ! ---------------------------------------------------------
  subroutine scalar_mesh_end(sm)
    type(scalar_mesh_t), intent(inout) :: sm

    call push_sub('scalar_mesh.scalar_mesh_end')

    deallocate(sm%mesh, sm%weights)

    call pop_sub()
  end subroutine scalar_mesh_end

end module scalar_mesh_m
