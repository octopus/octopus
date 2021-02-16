!! Copyright (C) 2021 N. Tancogne-Dejean
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

module lattice_vectors_oct_m
  use global_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m

  implicit none

  private

  public ::                   &
    lattice_vectors_t,           &
    build_metric_from_angles, &
    reciprocal_lattice

  type lattice_vectors_t
    ! Components are public by default
    FLOAT :: rlattice_primitive(MAX_DIM,MAX_DIM)   !< lattice primitive vectors
    FLOAT :: rlattice          (MAX_DIM,MAX_DIM)   !< lattice vectors
    FLOAT :: klattice_primitive(MAX_DIM,MAX_DIM)   !< reciprocal-lattice primitive vectors
    FLOAT :: klattice          (MAX_DIM,MAX_DIM)   !< reciprocal-lattice vectors
    logical :: nonorthogonal

  end type lattice_vectors_t

contains

  subroutine build_metric_from_angles(this, angles)
    type(lattice_vectors_t),    intent(inout) :: this
    FLOAT,                   intent(in)    :: angles(3)

    FLOAT :: cosang, a2, aa, cc
    FLOAT, parameter :: tol_angle = CNST(1.0e-6)

    PUSH_SUB(build_metric_from_angles)
   
    !Converting the angles to LatticeVectors
    !See 57_iovars/ingeo.F90 in Abinit for details
    if( abs(angles(1)-angles(2))< tol_angle .and. abs(angles(2)-angles(3))< tol_angle .and.  &
      (abs(angles(1)-CNST(90.0))+abs(angles(2)-CNST(90.0))+abs(angles(3)-CNST(90.0)))> tol_angle ) then

      cosang=cos(M_PI*angles(1)/CNST(180.0));
      a2=M_TWO/M_THREE*(M_ONE-cosang);
      aa=sqrt(a2);
      cc=sqrt(M_ONE-a2);
      this%rlattice_primitive(1,1) = aa
      this%rlattice_primitive(2,1) = M_ZERO
      this%rlattice_primitive(3,1) = cc
      this%rlattice_primitive(1,2) =-M_HALF*aa
      this%rlattice_primitive(2,2) = M_HALF*sqrt(M_THREE)*aa
      this%rlattice_primitive(3,2) = cc
      this%rlattice_primitive(1,3) =-M_HALF*aa
      this%rlattice_primitive(2,3) =-M_HALF*sqrt(M_THREE)*aa
      this%rlattice_primitive(3,3) = cc
    else
      this%rlattice_primitive(1,1) = M_ONE
      this%rlattice_primitive(2,1) = M_ZERO
      this%rlattice_primitive(3,1) = M_ZERO
      this%rlattice_primitive(1,2) = cos(M_PI*angles(3)/CNST(180.0))
      this%rlattice_primitive(2,2) = sin(M_PI*angles(3)/CNST(180.0))
      this%rlattice_primitive(3,2) = M_ZERO
      this%rlattice_primitive(1,3) = cos(M_PI*angles(2)/CNST(180.0))
      this%rlattice_primitive(2,3) = (cos(M_PI*angles(1)/CNST(180.0))-this%rlattice_primitive(1,2)* this%rlattice_primitive(1,3))&
                                       /this%rlattice_primitive(2,2)
      this%rlattice_primitive(3,3) = sqrt(M_ONE-this%rlattice_primitive(1,3)**2-this%rlattice_primitive(2,3)**2)
    end if

    if(any(abs(angles-CNST(90.0)) > M_EPSILON )) then
      this%nonorthogonal = .true.
    else
      this%nonorthogonal = .false.
    end if

    POP_SUB(build_metric_from_angles)
  end subroutine build_metric_from_angles


    !--------------------------------------------------------------
  subroutine reciprocal_lattice(rv, kv, volume, dim, namespace)
    FLOAT,             intent(in)  :: rv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: kv(:,:) !< (1:MAX_DIM, 1:MAX_DIM)
    FLOAT,             intent(out) :: volume
    integer,           intent(in)  :: dim
    type(namespace_t), intent(in)  :: namespace

    integer :: ii
    FLOAT :: cross(1:3), rv3(1:3, 1:3)

    PUSH_SUB(reciprocal_lattice)

    kv(:,:) = M_ZERO

    select case(dim)
    case(3)
      cross(1:3) = dcross_product(rv(1:3, 2), rv(1:3, 3)) 
      volume = dot_product(rv(1:3, 1), cross(1:3))

      kv(1:3, 1) = dcross_product(rv(:, 2), rv(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv(:, 3), rv(:, 1))/volume
      kv(1:3, 3) = dcross_product(rv(:, 1), rv(:, 2))/volume    
    case(2)
      rv3(1:3, 1:3) = M_ZERO
      rv3(1:2, 1:2) = rv(1:2, 1:2)
      rv3(3, 3) = M_ONE
      cross(1:3) = dcross_product(rv3(1:3, 1), rv3(1:3, 2)) 
      volume = dot_product(rv3(1:3, 3), cross(1:3))

      kv(1:3, 1) = dcross_product(rv3(:, 2), rv3(:, 3))/volume
      kv(1:3, 2) = dcross_product(rv3(:, 3), rv3(:, 1))/volume
    case(1)
      volume = rv(1, 1)
      kv(1, 1) = M_ONE / rv(1, 1)
    case default ! dim > 3
      message(1) = "Reciprocal lattice for dim > 3 assumes no periodicity."
      call messages_warning(1, namespace=namespace)
      volume = M_ONE
      do ii = 1, dim
        kv(ii, ii) = M_ONE/rv(ii,ii)
        !  At least initialize the thing
        volume = volume * sqrt(sum(rv(:, ii)**2))
      end do
    end select

    if ( volume < M_ZERO ) then 
      message(1) = "Your lattice vectors form a left-handed system."
      call messages_fatal(1, namespace=namespace)
    end if

    POP_SUB(reciprocal_lattice)
  end subroutine reciprocal_lattice


end module lattice_vectors_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

