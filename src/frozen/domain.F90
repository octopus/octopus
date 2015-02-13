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
  public ::        &
    domain__add__

  public ::           &
    domain_init,      &
    domain_start,     &
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

  type, public :: domain_iterator_t
    private
    type(domain_t),      pointer :: self =>null()
    type(domain_hash_iterator_t) :: iter
  end type domain_iterator_t

  interface domain_init
    module procedure domain_init_domain
    module procedure domain_iterator_init
  end interface domain_init

  interface domain_next
    module procedure domain_iterator_next_config_domain
    module procedure domain_iterator_next_config
    module procedure domain_iterator_next_domain
  end interface domain_next

  interface domain_copy
    module procedure domain_copy_domain
    module procedure domain_iterator_copy
  end interface domain_copy

  interface domain_end
    module procedure domain_end_domain
    module procedure domain_iterator_end
  end interface domain_end

  integer, parameter :: DOMAIN_OK          = DOMAIN_HASH_OK
  integer, parameter :: DOMAIN_KEY_ERROR   = DOMAIN_HASH_KEY_ERROR
  integer, parameter :: DOMAIN_EMPTY_ERROR = DOMAIN_HASH_EMPTY_ERROR

contains
  
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine domain_init_domain(this)
    type(domain_t), intent(out) :: this
    !
    PUSH_SUB(domain_init_domain)
    call domain_hash_init(this%hash)
    POP_SUB(domain_init_domain)
    return
  end subroutine domain_init_domain

  ! ---------------------------------------------------------
  subroutine domain_start(this, sb, geo)
    type(domain_t),            intent(out) :: this
    type(simul_box_t), target, intent(in)  :: sb
    type(geometry_t),  target, intent(in)  :: geo
    !
    PUSH_SUB(domain_start)
    ASSERT(.not.associated(this%sb))
    ASSERT(.not.associated(this%geo))
    this%sb=>sb
    this%geo=>geo
    POP_SUB(domain_start)
    return
  end subroutine domain_start

  ! ---------------------------------------------------------
  subroutine domain__add__(this, that, config)
    type(domain_t),      intent(inout) :: this
    type(domain_t),      intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(domain__add__)
    call domain_hash_set(this%hash, config, that)
    POP_SUB(domain__add__)
    return
  end subroutine domain__add__

  ! ---------------------------------------------------------
  function domain_in_domain_aux(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x
    !
    logical :: in
    !
    PUSH_SUB(domain_in_domain_aux)
    ASSERT(associated(this%sb))
    ASSERT(associated(this%geo))
    in=simul_box_in_box(this%sb, this%geo, x)
    POP_SUB(domain_in_domain_aux)
    return
  end function domain_in_domain_aux

  ! ---------------------------------------------------------
  function domain_in_domain(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x
    !
    logical :: in
    !
    type(domain_hash_iterator_t) :: iter
    type(domain_t),      pointer :: domain
    integer                      :: ierr
    !
    PUSH_SUB(domain_in_domain)
    nullify(domain)
    in=domain_in_domain_aux(this, x)
    if(.not.in)then
      call domain_hash_init(iter, this%hash)
      do
        nullify(domain)
        call domain_hash_next(iter, domain, ierr)
        if(ierr/=DOMAIN_HASH_OK)exit
        in=domain_in_domain_aux(domain, x)
        if(in)exit
      end do
      call domain_hash_end(iter)
    end if
    POP_SUB(domain_in_domain)
    return
  end function domain_in_domain

  ! ---------------------------------------------------------
  subroutine domain_copy_domain(this, that)
    type(domain_t), intent(out) :: this
    type(domain_t), intent(in)  :: that
    !
    PUSH_SUB(domain_copy_domain)
    this%sb=>that%sb
    this%geo=>that%geo
    call domain_hash_copy(this%hash, that%hash)
    POP_SUB(domain_copy_domain)
    return
  end subroutine domain_copy_domain

  ! ---------------------------------------------------------
  subroutine domain_end_domain(this)
    type(domain_t), intent(inout) :: this
    !
    PUSH_SUB(domain_end_domain)
    nullify(this%sb, this%geo)
    call domain_hash_end(this%hash)
    POP_SUB(domain_end_domain)
    return
  end subroutine domain_end_domain

  ! ---------------------------------------------------------
  subroutine domain_iterator_init(this, that)
    type(domain_iterator_t), intent(out) :: this
    type(domain_t),  target, intent(in)  :: that
    !
    PUSH_SUB(domain_iterator_init)
    this%self=>that
    call domain_hash_init(this%iter, that%hash)
    POP_SUB(domain_iterator_init)
    return
  end subroutine domain_iterator_init

  ! ---------------------------------------------------------
  subroutine domain_iterator_next_config_domain(this, config, sim, ierr)
    type(domain_iterator_t), intent(inout) :: this
    type(json_object_t),    pointer        :: config
    type(domain_t),         pointer        :: sim
    integer,       optional, intent(out)   :: ierr
    !
    PUSH_SUB(domain_iterator_next_config_domain)
    call domain_hash_next(this%iter, config, sim, ierr)
    POP_SUB(domain_iterator_next_config_domain)
    return
  end subroutine domain_iterator_next_config_domain

  ! ---------------------------------------------------------
  subroutine domain_iterator_next_config(this, that, ierr)
    type(domain_iterator_t), intent(inout) :: this
    type(json_object_t),    pointer        :: that
    integer,       optional, intent(out)   :: ierr
    !
    PUSH_SUB(domain_iterator_next_config)
    call domain_hash_next(this%iter, that, ierr)
    POP_SUB(domain_iterator_next_config)
    return
  end subroutine domain_iterator_next_config

  ! ---------------------------------------------------------
  subroutine domain_iterator_next_domain(this, that, ierr)
    type(domain_iterator_t), intent(inout) :: this
    type(domain_t),         pointer        :: that
    integer,       optional, intent(out)   :: ierr
    !
    PUSH_SUB(domain_iterator_next_domain)
    call domain_hash_next(this%iter, that, ierr)
    POP_SUB(domain_iterator_next_domain)
    return
  end subroutine domain_iterator_next_domain

  ! ---------------------------------------------------------
  subroutine domain_iterator_copy(this, that)
    type(domain_iterator_t), intent(inout) :: this
    type(domain_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(domain_iterator_copy)
    this%self=>that%self
    call domain_hash_copy(this%iter, that%iter)
    POP_SUB(domain_iterator_copy)
    return
  end subroutine domain_iterator_copy

  ! ---------------------------------------------------------
  subroutine domain_iterator_end(this)
    type(domain_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(domain_iterator_end)
    nullify(this%self)
    call domain_hash_end(this%iter)
    POP_SUB(domain_iterator_end)
    return
  end subroutine domain_iterator_end

end module domain_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
