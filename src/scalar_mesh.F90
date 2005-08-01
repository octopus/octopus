!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

#include "global.h"

module scalar_mesh
  use global
  use messages
  use syslabels
  use math
  use io

  implicit none

  private
  public :: scalar_mesh_type
  public :: scalar_mesh_init, scalar_mesh_create, scalar_mesh_end
  public :: scalar_mesh_integrate, scalar_mesh_write

  integer, parameter ::         &
       MESH_LINEAR         = 1, &
       MESH_LOG            = 2, &
       MESH_DOUBLE_LOG     = 3, &
       MESH_SINH           = 4, &
       MESH_GAUSS_LEGENDRE = 5

  type scalar_mesh_type
     integer :: mtype            ! mesh type (see MESH_* variables above)
     integer :: np               ! number of points in the mesh
     FLOAT   :: min, max         ! lower and upper boundary for scalar grid
     FLOAT   :: center           ! highest density of gridpoints is around center
     FLOAT   :: alpha1, alpha2   ! parameters for logarithmic grid below(1) and above(2) center
     FLOAT   :: dx               ! grid spacing for linear grid
     FLOAT, pointer :: mesh(:)   ! the mesh points
     FLOAT, pointer :: w(:)      ! weights for integrations
  end type scalar_mesh_type

contains

  ! ---------------------------------------------------------
  subroutine scalar_mesh_init(sm)
    type(scalar_mesh_type), intent(inout) :: sm

    FLOAT :: etmp

    call push_sub('scalar_mesh.scalar_mesh_init')
    !%Variable ScalarMeshType
    !%Type integer
    !%Section 1 Generalities
    !%Description
    !% Specifies what kind of scalar mesh will be used
    !%Option mesh_linear 1
    !% Linear mesh
    !%Option mesh_double_log 2
    !% Logarithmic mesh
    !%Option mesh_double_log 3
    !% Double logarithmic mesh
    !%Option mesh_sinh 4
    !% Sinh mesh
    !%Option gauss_legendre 5
    !% Gauss-Legendre mesh
    !%End
    call loct_parse_int  (check_inp('ScalarMeshType'),             4, sm%mtype)
    call loct_parse_int  (check_inp('ScalarMeshNPoints'),         20, sm%np)
    call loct_parse_float(check_inp('ScalarMeshMin'),      CNST(0.0), sm%min)
    call loct_parse_float(check_inp('ScalarMeshMax'),      CNST(4.0), sm%max)
    call loct_parse_float(check_inp('ScalarMeshCenter'),   CNST(2.0), sm%center)
    call loct_parse_float(check_inp('ScalarMeshAlpha1'),   CNST(0.3), sm%alpha1)
    call loct_parse_float(check_inp('ScalarMeshAlpha2'),   CNST(0.3), sm%alpha2)

    ! a few sanity checks
    if(sm%min.ge.sm%max) then
       message(1) = 'Warning: scalar_mesh: ScalarMeshMin is equal to or larger than ScalarMeshMax.'
       message(2) = 'Reverting order.'
       call write_warning(2)
       etmp    = sm%min
       sm%min = sm%max
       sm%max = etmp
    endif

    if(sm%center.lt.sm%min.or.sm%center.gt.sm%max) then
       message(1) = 'Warning: scalar_mesh: ScalarMeshCenter is not in interval [ScalarMeshMin,ScalarMeshMax].'
       message(2) = 'Will use the average of ScalarMeshMin and ScalarMeshMax.'
       call write_warning(2)
       sm%center = M_HALF*(sm%min+sm%max)
    endif

    call scalar_mesh_create(sm)

    call pop_sub()
  end subroutine scalar_mesh_init


  ! ---------------------------------------------------------
  subroutine scalar_mesh_create(sm)
    type(scalar_mesh_type), intent(out) :: sm

    FLOAT :: xmin1, xmax1, xmin2, xmax2, offset
    FLOAT, allocatable :: gl(:), gw(:)
    integer :: i, iunit

    call push_sub('scalar_mesh.scalar_mesh_create')

    allocate(sm%mesh(2*sm%np+1),sm%w(2*sm%np+1))

    xmin1  = M_ZERO
    xmax1  = sm%center - sm%min
    xmin2  = M_ZERO
    xmax2  = sm%max - sm%center
    offset = sm%center
    sm%dx  = (sm%max-sm%min)/(2*sm%np)

    select case(sm%mtype)
    case(MESH_LINEAR)
       message(1) = 'Info: scalar_mesh: Using linear scalar mesh.'
       call write_info(1)

       ! setup linear mesh
       if (sm%np.eq.0) then
          sm%mesh = sm%center
       else
          do i = 1, 2*sm%np+1
             sm%mesh(i) = sm%min + (sm%max-sm%min)/(2*sm%np)*(i-1)
          enddo
       endif
       ! setup integration weights
       sm%w = sm%dx
       if (2*sm%np+1.lt.10) then
          sm%w = sm%dx
       else
          do i = 1, 5
             sm%w(i)             = EMcLCoeff(i)*sm%dx
             sm%w(2*sm%np+2 - i) = EMcLCoeff(i)*sm%dx
          enddo
       endif


    case(MESH_LOG)
       message(1) = 'Info: scalar_mesh: Using logarithmic scalar mesh.'
       call write_info(1)

       do i = 1, 2*sm%np+1
          sm%mesh(i) =   sm%min   + (sm%max - sm%min)* &
               ( exp(-sm%alpha1*(i-1)) - M_ONE )/(exp(-sm%alpha1*(2*sm%np))-M_ONE)

          sm%w(i)   = - sm%alpha1 * (sm%max - sm%min)* &
               ( exp(-sm%alpha1*(i-1))         )/(exp(-sm%alpha1*(2*sm%np))-M_ONE)
       enddo
       do i = 1, 5
          sm%w(i)             = EMcLCoeff(i)*sm%w(i)
          sm%w(2*sm%np+2 - i) = EMcLCoeff(i)*sm%w(2*sm%np+2 - i)
       enddo


    case(MESH_DOUBLE_LOG)
       message(1) = 'Info: scalar_mesh: Using double logarithmic scalar mesh.'
       call write_info(1)

       ! setup double logarithmic mesh
       sm%mesh = offset
       do i = 1, sm%np
          sm%mesh(sm%np-i+1) = - ( offset + xmin1 + (xmax1-xmin1) * &
               ( exp(sm%alpha1*(i))     - M_ONE ) /                 &
               ( exp(sm%alpha1*(sm%np)) - M_ONE ) ) + 2*offset

          sm%mesh(sm%np+i+1) =     offset + xmin2 + (xmax2-xmin2) * &
               ( exp(sm%alpha2*(i))     - M_ONE ) /                 &
               ( exp(sm%alpha2*(sm%np)) - M_ONE )
       enddo
       ! setup integration weights
       do i = 1, sm%np
          sm%w(i)             = sm%alpha1 * (xmax1-xmin1) /        &
               ( exp(sm%alpha1*(sm%np))     - M_ONE ) *             &
               ( exp(sm%alpha1*(sm%np-i+1))         )
          sm%w(2*sm%np+2 - i) = sm%alpha2 * (xmax2-xmin2) /        &
               ( exp(sm%alpha2*(sm%np))     - M_ONE ) *             &
               ( exp(sm%alpha2*(sm%np-i+1))         )
       enddo
       sm%w(sm%np+1)   = EMcLCoeff(1)*                                           &
            ( sm%alpha1 * (xmax1-xmin1) / ( exp(sm%alpha1*(sm%np)) - M_ONE ) &
            + sm%alpha2 * (xmax2-xmin2) / ( exp(sm%alpha2*(sm%np)) - M_ONE ) )
       sm%w(1)         = EMcLCoeff(1)*sm%w(1)
       sm%w(2*sm%np+1) = EMcLCoeff(1)*sm%w(2*sm%np+1)
       do i = 2, 5
          sm%w(i)             = EMcLCoeff(i)*sm%w(i)
          sm%w(sm%np+2 - i)   = EMcLCoeff(i)*sm%w(sm%np+2 - i)
          sm%w(sm%np+i)       = EMcLCoeff(i)*sm%w(sm%np+i)
          sm%w(2*sm%np+2 - i) = EMcLCoeff(i)*sm%w(2*sm%np+2 - i)
       enddo


    case(MESH_SINH)
       message(1) = 'Info: scalar_mesh: Using sinh scalar mesh.'
       call write_info(1)

       ! setup sinh mesh
       do i = 1, 2*sm%np+1
          sm%mesh(i) = sm%center + M_HALF*(sm%max-sm%min) *        &
               sinh(sm%alpha1*(i-sm%np-1))/sinh(sm%alpha1*(sm%np))
       enddo
       do i = 1, 2*sm%np+1
          sm%w(i)   =  sm%alpha1 * M_HALF*(sm%max-sm%min) /        &
               sinh(sm%alpha1*(sm%np)) * cosh(sm%alpha1*(sm%np-i+1))
       enddo
       do i = 1, 5
          sm%w(i)             = EMcLCoeff(i)*sm%w(i)
          sm%w(2*sm%np+2 - i) = EMcLCoeff(i)*sm%w(2*sm%np+2 - i)
       enddo

    case(MESH_GAUSS_LEGENDRE)
       message(1) = 'Info: scalar_mesh: Using Gauss-Legendre scalar mesh.'
       call write_info(1)

       allocate(gl(2*sm%np+1), gw(2*sm%np+1))

       ! TODO: need to compute roots of Legendre Polynomials and weights for integration
       ! a call of the corresponding routine in math.F90 should be placed here

       do i=1, 2*sm%np+1
          sm%mesh(i) = M_HALF*(sm%max-sm%min) * gl(i) + M_HALF*(sm%max+sm%min)
          sm%w(i)   = gw(i)*M_HALF*(sm%max-sm%min)
       enddo

       deallocate(gl, gw)

    case default
       write(message(1), '(a,i4,a)') "Input: '", sm%mtype, &
            "' is not a valid scalar mesh type"
       message(2) = '( ScalarMeshType =  Mesh_Linear | Mesh_Log | Mesh_Double_Log | Mesh_Sinh )'
       call write_fatal(2)
    end select

    call scalar_mesh_write(sm, 'scalar_mesh.dat')

    call pop_sub()
  end subroutine scalar_mesh_create


  ! ---------------------------------------------------------
  subroutine scalar_mesh_write(sm, filename)
    type(scalar_mesh_type), intent(in) :: sm
    character(len=*)                   :: filename

    integer i, iunit

    call push_sub('scalar_mesh.scalar_mesh_write')
    ! output scalar mesh
    iunit = io_open('static/'//trim(filename), action='write')
    do i = 1, 2*sm%np+1
       write(iunit,'(1x,i,2x,2F22.16)') i,sm%mesh(i),sm%w(i)
    enddo
    call io_close(iunit)
    call pop_sub()
  end subroutine scalar_mesh_write


  ! ---------------------------------------------------------
  FLOAT function scalar_mesh_integrate(sm, func) result(res)
    type(scalar_mesh_type), intent(in) :: sm
    FLOAT, intent(in) :: func(:)

    integer :: i

    call push_sub('scalar_mesh.scalar_mesh_end')

    res = sum(sm%w(:)*func(:))

    call pop_sub()
  end function scalar_mesh_integrate


  ! ---------------------------------------------------------
  subroutine scalar_mesh_end(sm)
    type(scalar_mesh_type), intent(inout) :: sm

    call push_sub('scalar_mesh.scalar_mesh_end')

    deallocate(sm%mesh, sm%w)

    call pop_sub()
  end subroutine scalar_mesh_end

end module scalar_mesh
