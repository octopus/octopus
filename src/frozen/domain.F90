#include "global.h"

module domain_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m,     only: wp
  use geometry_m,  only: geometry_t
  use simul_box_m, only: simul_box_t, simul_box_in_box

  implicit none

  private
  public ::           &
    domain_init,      &
    domain_in_domain, &
    domain_copy,      &
    domain_end

  type, public :: domain_t
    private
    type(simul_box_t), pointer :: sb  =>null()
    type(geometry_t),  pointer :: geo =>null()
  end type domain_t

contains
  
  ! ---------------------------------------------------------
  subroutine domain_init(this, sb, geo)
    type(domain_t),            intent(out) :: this
    type(simul_box_t), target, intent(in)  :: sb
    type(geometry_t),  target, intent(in)  :: geo
    !
    PUSH_SUB(domain_init)
    this%sb=>sb
    this%geo=>geo
    POP_SUB(domain_init)
    return
  end subroutine domain_init

  ! ---------------------------------------------------------
  function domain_in_domain(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x
    !
    logical :: in
    !
    PUSH_SUB(domain_in_domain)
    in=simul_box_in_box(this%sb, this%geo, x)
    POP_SUB(domain_in_domain)
    return
  end function domain_in_domain

  ! ---------------------------------------------------------
  subroutine domain_copy(this_out, this_in)
    type(domain_t), intent(out) :: this_out
    type(domain_t), intent(in)  :: this_in
    !
    PUSH_SUB(domain_copy)
    this_out%sb=>this_in%sb
    this_out%geo=>this_in%geo
    POP_SUB(domain_copy)
    return
  end subroutine domain_copy

  ! ---------------------------------------------------------
  elemental subroutine domain_end(this)
    type(domain_t), intent(inout) :: this
    !
    nullify(this%sb, this%geo)
    return
  end subroutine domain_end

end module domain_m

!! Local Variables:
!! mode: f90
!! End:
