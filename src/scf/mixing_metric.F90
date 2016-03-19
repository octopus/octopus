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
!! $Id: mix.F90 14715 2015-10-30 05:54:23Z xavier $

#include "global.h"

module mixing_metric_oct_m
  use derivatives_oct_m
  use global_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use nl_operator_oct_m
  use profiling_oct_m
  use stencil_cube_oct_m

  implicit none

  private
  public ::                               &
    mixing_metric_init,                   &
    mixing_metric_end,                    &
    mixing_metric_dotp

  ! Note, since this object has to export mixing_metric_dotp without
  ! an extra argument, we need to make this a global information
  ! object.

  type(derivatives_t), pointer, save :: der_ptr
  type(nl_operator_t), save          :: op
  
contains

  ! ---------------------------------------------------------
  subroutine mixing_metric_init(der, weight)
    type(derivatives_t), target, intent(in) :: der
    FLOAT,                       intent(in) :: weight

    integer :: ip, maxp, is, inb

    PUSH_SUB(mixing_metric_init)

    ASSERT(.not. der%mesh%use_curvilinear)

    der_ptr => der

    ! This is the metric operator used by GPAW:
    !
    ! https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html
    !
    ! The actual coefficients were copied from the file gpaw/mixing.py
    ! from the gpaw source code.

    call nl_operator_init(op, "Mix metric")
    call stencil_cube_get_lapl(op%stencil, der%mesh%sb%dim, 1)
    
    call nl_operator_build(der%mesh, op, der%mesh%np, const_w = .not. der%mesh%use_curvilinear)

    if (op%const_w) then
      maxp = 1
    else
      maxp = der%mesh%np
    end if
    
    do ip = 1, maxp
      
      do is = 1, op%stencil%size

        inb = sum(abs(op%stencil%points(1:der%mesh%sb%dim, is)))

        select case(inb)
        case(0)
          !          op%w_re(is, ip) = CNST(1.0) + weight/CNST(8.0)
          op%w_re(is, ip) = CNST(0.125)*(weight + CNST(7.0))
        case(1)
          !          op%w_re(is, ip) = weight/CNST(16.0)
          op%w_re(is, ip) = CNST(0.0625)*(weight - CNST(1.0))
        case(2)
          !          op%w_re(is, ip) = weight/CNST(32.0)
          op%w_re(is, ip) = CNST(0.03125)*(weight - CNST(1.0))
        case(3)
          !          op%w_re(is, ip) = weight/CNST(64.0)
          op%w_re(is, ip) = CNST(0.015625)*(weight - CNST(1.0))
        case default
          ASSERT(.false.)
        end select
        
      end do

      print*, "COSA", ip, sum(op%w_re(1:op%stencil%size, ip))

    end do
    
    call nl_operator_update_weights(op)

    POP_SUB(mixing_metric_init)
  end subroutine mixing_metric_init

  ! ---------------------------------------------------------
  subroutine mixing_metric_end()

    PUSH_SUB(mixing_metric_end)

    call nl_operator_end(op)
    nullify(der_ptr)

    POP_SUB(mixing_metric_end)
  end subroutine mixing_metric_end

  ! --------------------------------------------------------

  FLOAT function mixing_metric_dotp(xx, yy) result(res)
    FLOAT, intent(in) :: xx(:)
    FLOAT, intent(in) :: yy(:)

    FLOAT, allocatable :: yyc(:), opyy(:)

    PUSH_SUB(mixing_metric_dotp)

    SAFE_ALLOCATE(yyc(1:der_ptr%mesh%np_part))
    SAFE_ALLOCATE(opyy(1:der_ptr%mesh%np))

    yyc(1:der_ptr%mesh%np) = yy(1:der_ptr%mesh%np)

    call dderivatives_perform(op, der_ptr, yyc, opyy)

    res = dmf_dotp(der_ptr%mesh, xx, opyy)

    print*, res, dmf_dotp(der_ptr%mesh, xx, yy)
    
    SAFE_DEALLOCATE_A(yyc)
    SAFE_DEALLOCATE_A(opyy)

    POP_SUB(mixing_metric_dotp)

  end function mixing_metric_dotp

end module mixing_metric_oct_m
  
