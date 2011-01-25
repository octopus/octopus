!! Copyright (C) 2006 Hyllios
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
!! $Id: em_resp.F90 2647 2007-01-09 18:02:46Z lorenzen $

#include "global.h"

module gauss_legendre_m
  use global_m

  implicit none

  private
  public :: &
    gauss_legendre_points
  

  ! 2 point
  FLOAT, target :: GL_points_2(1) = (/                                   &
    CNST(-0.57735026919) /)
  FLOAT, target :: GL_weights_2(1) = (/                                  &
    CNST(1.0) /)

  ! 3 point
  FLOAT, target :: GL_points_3(2) = (/                                   &
    CNST(-0.774596669241), CNST(0.0) /)
  FLOAT, target :: GL_weights_3(2) = (/                                  &
    CNST(0.555555555556), CNST(0.888888888889) /)

  ! 4 point
  FLOAT, target :: GL_points_4(2) = (/                                   &
    CNST(-0.861136311594), CNST(-0.339981043585) /)
  FLOAT, target :: GL_weights_4(2) = (/                                  &
    CNST(0.347854845137), CNST(0.652145154863) /)

  ! 5 point
  FLOAT, target :: GL_points_5(3) = (/                                   &
    CNST(-0.906179845939), CNST(-0.538469310106), CNST(0.0) /)
  FLOAT, target :: GL_weights_5(3) = (/                                  &
    CNST(0.236926885056), CNST(0.478628670499), CNST(0.568888888889) /)

  ! 6 point
  FLOAT, target :: GL_points_6(3) = (/                                   &
    CNST(-0.932469514203), CNST(-0.661209386466), CNST(-0.238619186083) /)
  FLOAT, target :: GL_weights_6(3) = (/                                  &
    CNST(0.171324492379), CNST(0.360761573048), CNST(0.467913934573) /)
  
  ! 7 point
  FLOAT, target :: GL_points_7(4) = (/                                   &
    CNST(-0.949107912343), CNST(-0.741531185599), CNST(-0.405845151377), &
    CNST(0.0) /)
  FLOAT, target :: GL_weights_7(4) = (/                                  &
    CNST(0.129484966169), CNST(0.279705391489), CNST(0.381830050505),    &
    CNST(0.417959183673) /)

  ! 8 point
  FLOAT, target :: GL_points_8(4) = (/                                   &
    CNST(-0.960289856498), CNST(-0.796666477414), CNST(-0.525532409916), &
    CNST(-0.183434642496) /)
  FLOAT, target :: GL_weights_8(4) = (/                                  &
    CNST(0.10122853629), CNST(0.222381034453), CNST(0.313706645878),     &
    CNST(0.362683783378) /)

  ! 9 point
  FLOAT, target :: GL_points_9(5) = (/                                   &
    CNST(-0.968160239508), CNST(-0.836031107327), CNST(-0.613371432701), &
    CNST(-0.324253423404), CNST(0.0) /)
  FLOAT, target :: GL_weights_9(5) = (/                                  &
    CNST(0.0812743883616), CNST(0.180648160695), CNST(0.260610696403),   &
    CNST(0.31234707704), CNST(0.330239355001) /)

  ! 10 point
  FLOAT, target :: GL_points_10(5) = (/                                  &
    CNST(-0.973906528517), CNST(-0.865063366689), CNST(-0.679409568299), &
    CNST(-0.433395394129), CNST(-0.148874338982) /)
  FLOAT, target :: GL_weights_10(5) = (/                                 &
    CNST(0.0666713443087), CNST(0.149451349151), CNST(0.219086362516),   &
    CNST(0.26926671931),   CNST(0.295524224715) /)

  ! 11 point
  FLOAT, target :: GL_points_11(6) = (/                                  &
    CNST(-0.978228658146), CNST(-0.887062599768), CNST(-0.730152005574), &
    CNST(-0.519096129207), CNST(-0.269543155952), CNST(0.0) /)
  FLOAT, target :: GL_weights_11(6) = (/                                 &
    CNST(0.0556685671162), CNST(0.125580369465), CNST(0.186290210928),   &
    CNST(0.233193764592),  CNST(0.26280454451),  CNST(0.272925086778) /)

  ! 12 point
  FLOAT, target :: GL_points_12(6) = (/                                  &
    CNST(-0.981560634247), CNST(-0.90411725637),  CNST(-0.769902674194), &
    CNST(-0.587317954287), CNST(-0.367831498998), CNST(-0.125233408511) /)
  FLOAT, target :: GL_weights_12(6) = (/                                 &
    CNST(0.0471753363866), CNST(0.106939325995), CNST(0.160078328543),   &
    CNST(0.203167426723),  CNST(0.233492536538), CNST(0.249147045813) /)

contains

  ! --------------------------------------------------------------------
  subroutine gauss_legendre_points(n, points, weights)
    integer, intent(in)  :: n
    FLOAT,   intent(out) :: points(:)
    FLOAT,   intent(out) :: weights(:)

    integer :: i

    FLOAT, pointer :: points_ref(:), weights_ref(:)

    nullify(points_ref)
    nullify(weights_ref)

    ASSERT(n>=2.and.n<=12)
    select case(n)
    case(2);  points_ref  => GL_points_2;  weights_ref => GL_weights_2
    case(3);  points_ref  => GL_points_3;  weights_ref => GL_weights_3
    case(4);  points_ref  => GL_points_4;  weights_ref => GL_weights_4
    case(5);  points_ref  => GL_points_5;  weights_ref => GL_weights_5
    case(6);  points_ref  => GL_points_6;  weights_ref => GL_weights_6
    case(7);  points_ref  => GL_points_7;  weights_ref => GL_weights_7
    case(8);  points_ref  => GL_points_8;  weights_ref => GL_weights_8
    case(9);  points_ref  => GL_points_9;  weights_ref => GL_weights_9
    case(10); points_ref  => GL_points_10; weights_ref => GL_weights_10
    case(11); points_ref  => GL_points_11; weights_ref => GL_weights_11
    case(12); points_ref  => GL_points_12; weights_ref => GL_weights_12
    end select
    
    do i = 1, n/2
      points (i)     =  points_ref(i)
      points (n-i+1) = -points_ref(i)
      weights(i)     =  weights_ref(i)
      weights(n-i+1) =  weights_ref(i)
    end do
    if(mod(n, 2) == 1) then
      points (n/2 + 1) = points_ref (n/2 + 1)
      weights(n/2 + 1) = weights_ref(n/2 + 1)
    end if

  end subroutine gauss_legendre_points

end module gauss_legendre_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
