#include "global.h"

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_KEY_TYPE_MODULE_NAME
#undef HASH_KEY_FUNCTION_NAME
#undef HASH_KEY_FUNCTION_MODULE_NAME
#undef HASH_VAL_TEMPLATE_NAME
#undef HASH_VAL_TYPE_NAME
#undef HASH_VAL_TYPE_MODULE_NAME
#undef HASH_INCLUDE_PREFIX
#undef HASH_INCLUDE_HEADER
#undef HASH_INCLUDE_BODY

#define HASH_TEMPLATE_NAME domain
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME domain

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module domain_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_object_t, json_hash
  use kinds_m,  only: wp

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

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: domain_t
    private
    type(simul_box_t), pointer :: sb  =>null()
    type(geometry_t),  pointer :: geo =>null()
    type(domain_hash_t)        :: hash
  end type domain_t

contains
  
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine domain_init(this, sb, geo)
    type(domain_t),            intent(out) :: this
    type(simul_box_t), target, intent(in)  :: sb
    type(geometry_t),  target, intent(in)  :: geo
    !
    PUSH_SUB(domain_init)
    this%sb=>sb
    this%geo=>geo
    call domain_hash_init(this%hash)
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
    in=.false.
    in=simul_box_in_box(this%sb, this%geo, x)
    POP_SUB(domain_in_domain)
    return
  end function domain_in_domain

  ! ---------------------------------------------------------
  subroutine domain_copy(this, that)
    type(domain_t), intent(out) :: this
    type(domain_t), intent(in)  :: that
    !
    PUSH_SUB(domain_copy)
    this%sb=>that%sb
    this%geo=>that%geo
    call domain_hash_copy(this%hash, that%hash)
    POP_SUB(domain_copy)
    return
  end subroutine domain_copy

  ! ---------------------------------------------------------
  subroutine domain_end(this)
    type(domain_t), intent(inout) :: this
    !
    nullify(this%sb, this%geo)
    call domain_hash_end(this%hash)
    return
  end subroutine domain_end

end module domain_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
