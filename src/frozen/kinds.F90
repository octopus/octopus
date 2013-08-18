#include "global.h"

module kinds_m

  implicit none

  private
  public ::               &
    operator(.equal.),    &
    operator(.notequal.), &
    operator(.less.),     &
    operator(.more.)

  FLOAT, parameter :: kind_parm=CNST(1.0)

  integer, public, parameter :: wp=kind(kind_parm)
  integer, public, parameter :: sp=kind(1.0E+0)
  integer, public, parameter :: dp=kind(1.0D+0)

  integer, parameter :: delta_step = 2

  interface operator(.equal.)
    module procedure real_sp_equal
    module procedure real_dp_equal
  end interface operator(.equal.)

  interface operator(.notequal.)
    module procedure real_sp_not_equal
    module procedure real_dp_not_equal
  end interface operator(.notequal.)

  interface operator(.less.)
    module procedure real_sp_less
    module procedure real_dp_less
  end interface operator(.less.)

  interface operator(.more.)
    module procedure real_sp_more
    module procedure real_dp_more
  end interface operator(.more.)

contains

  ! -----------------------------------------------------
  elemental function real_sp_delta(this, that) result(delta)
    real(kind=sp), intent(in) :: this
    real(kind=sp), intent(in) :: that
    !
    real(kind=sp) :: delta
    !
    delta=real(delta_step, kind=sp)*spacing(max(this, that))
    return
  end function real_sp_delta
  
  ! -----------------------------------------------------
  elemental function real_sp_equal(this, that) result(eqv)
    real(kind=sp), intent(in) :: this
    real(kind=sp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(abs(this-that)<real_sp_delta(this, that))
    return
  end function real_sp_equal
  
  ! -----------------------------------------------------
  elemental function real_sp_not_equal(this, that) result(neqv)
    real(kind=sp), intent(in) :: this
    real(kind=sp), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=.not.real_sp_equal(this, that)
    return
  end function real_sp_not_equal

  ! -----------------------------------------------------
  elemental function real_sp_less(this, that) result(eqv)
    real(kind=sp), intent(in) :: this
    real(kind=sp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(this<that-real_sp_delta(this, that))
    return
  end function real_sp_less
  
  ! -----------------------------------------------------
  elemental function real_sp_more(this, that) result(eqv)
    real(kind=sp), intent(in) :: this
    real(kind=sp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(this>that+real_sp_delta(this, that))
    return
  end function real_sp_more
  
  ! -----------------------------------------------------
  elemental function real_dp_delta(this, that) result(delta)
    real(kind=dp), intent(in) :: this
    real(kind=dp), intent(in) :: that
    !
    real(kind=dp) :: delta
    !
    delta=real(delta_step, kind=dp)*spacing(max(this, that))
    return
  end function real_dp_delta
  
  ! -----------------------------------------------------
  elemental function real_dp_equal(this, that) result(eqv)
    real(kind=dp), intent(in) :: this
    real(kind=dp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(abs(this-that)<real_dp_delta(this, that))
    return
  end function real_dp_equal
  
  ! -----------------------------------------------------
  elemental function real_dp_not_equal(this, that) result(neqv)
    real(kind=dp), intent(in) :: this
    real(kind=dp), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=.not.real_dp_equal(this, that)
    return
  end function real_dp_not_equal

  ! -----------------------------------------------------
  elemental function real_dp_less(this, that) result(eqv)
    real(kind=dp), intent(in) :: this
    real(kind=dp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(this<that-real_dp_delta(this, that))
    return
  end function real_dp_less
  
  ! -----------------------------------------------------
  elemental function real_dp_more(this, that) result(eqv)
    real(kind=dp), intent(in) :: this
    real(kind=dp), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(this>that+real_dp_delta(this, that))
    return
  end function real_dp_more
  
end module kinds_m

!! Local Variables:
!! mode: f90
!! End:
