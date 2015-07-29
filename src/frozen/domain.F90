#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

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

module domain_m

  use global_m
  use messages_m
  use profiling_m

#define LIST_TEMPLATE_NAME domain
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

  use kinds_m
  use geometry_m
  use simul_box_m
  use json_m
  use config_dict_m

#define TEMPLATE_PREFIX domain
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private
  public ::          &
    domain__init__,  &
    domain__start__, &
    domain__add__,   &
    domain__copy__,  &
    domain__end__

  public ::           &
    domain_new,       &
    domain_del,       &
    domain_init,      &
    domain_start,     &
    domain_in_domain, &
    domain_get,       &
    domain_copy,      &
    domain_end

#define LIST_TEMPLATE_NAME domain
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: domain_t
    private
    type(simul_box_t), pointer :: sb   =>null()
    type(geometry_t),  pointer :: geo  =>null()
    type(domain_t),    pointer :: prnt =>null()
    type(config_dict_t)        :: dict
    type(domain_hash_t)        :: hash
    type(domain_list_t)        :: list
  end type domain_t

  interface domain_init
    module procedure domain_init_domain
    module procedure domain_init_copy
  end interface domain_init

  interface domain_get
    module procedure domain_get_domain_by_config
    module procedure domain_get_domain_by_name
  end interface domain_get

  interface domain_copy
    module procedure domain_copy_domain
  end interface domain_copy

  interface domain_end
    module procedure domain_end_domain
  end interface domain_end

  integer, public, parameter :: DOMAIN_OK          = DOMAIN_HASH_OK
  integer, public, parameter :: DOMAIN_KEY_ERROR   = DOMAIN_HASH_KEY_ERROR
  integer, public, parameter :: DOMAIN_EMPTY_ERROR = DOMAIN_HASH_EMPTY_ERROR

#define TEMPLATE_PREFIX domain
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains
  
#define LIST_TEMPLATE_NAME domain
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine domain_new(this, that)
    type(domain_t),  target, intent(inout) :: this
    type(domain_t), pointer                :: that
    !
    PUSH_SUB(domain_new)
    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt=>this
    call domain_list_push(this%list, that)
    POP_SUB(domain_new)
    return
  end subroutine domain_new

  ! ---------------------------------------------------------
  subroutine domain__idel__(this)
    type(domain_t), pointer :: this
    !
    PUSH_SUB(domain__idel__)
    SAFE_DEALLOCATE_P(this)
    nullify(this)
    POP_SUB(domain__idel__)
    return
  end subroutine domain__idel__

  ! ---------------------------------------------------------
  subroutine domain_del(this)
    type(domain_t), pointer :: this
    !
    PUSH_SUB(domain_del)
    if(associated(this))then
      if(associated(this%prnt))then
        call domain_list_del(this%prnt%list, this)
        call domain_end(this)
        call domain__idel__(this)
      end if
    end if
    POP_SUB(domain_del)
    return
  end subroutine domain_del

  ! ---------------------------------------------------------
  subroutine domain__inull__(this)
    type(domain_t), intent(inout) :: this
    !
    PUSH_SUB(domain__inull__)
    nullify(this%sb, this%geo, this%prnt)
    POP_SUB(domain__inull__)
    return
  end subroutine domain__inull__

  ! ---------------------------------------------------------
  subroutine domain__init__(this)
    type(domain_t), intent(out) :: this
    !
    PUSH_SUB(domain__init__)
    call domain__inull__(this)
    call config_dict_init(this%dict)
    call domain_hash_init(this%hash)
    call domain_list_init(this%list)
    POP_SUB(domain__init__)
    return
  end subroutine domain__init__

  ! ---------------------------------------------------------
  subroutine domain_init_domain(this)
    type(domain_t), intent(out) :: this
    !
    PUSH_SUB(domain_init_domain)
    call domain__init__(this)
    POP_SUB(domain_init_domain)
    return
  end subroutine domain_init_domain

  ! ---------------------------------------------------------
  recursive subroutine domain_init_copy(this, that)
    type(domain_t), intent(out) :: this
    type(domain_t), intent(in)  :: that
    !
    type(domain_iterator_t)      :: iter
    type(domain_t),      pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(domain_init_copy)
    nullify(cnfg, osub, isub)
    call domain__init__(this)
    call domain_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call domain_next(iter, cnfg, isub, ierr)
      if(ierr/=DOMAIN_OK)exit
      call domain_new(this, osub)
      call domain_init(osub, isub)
      call domain__add__(this, osub, cnfg)
    end do
    call domain_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(domain_init_copy)
    return
  end subroutine domain_init_copy

  ! ---------------------------------------------------------
  subroutine domain__istart__(this, sb, geo)
    type(domain_t),            intent(out) :: this
    type(simul_box_t), target, intent(in)  :: sb
    type(geometry_t),  target, intent(in)  :: geo
    !
    PUSH_SUB(domain__istart__)
    ASSERT(.not.associated(this%sb))
    ASSERT(.not.associated(this%geo))
    this%sb=>sb
    this%geo=>geo
    POP_SUB(domain__istart__)
    return
  end subroutine domain__istart__

  ! ---------------------------------------------------------
  subroutine domain__start__(this, sb, geo)
    type(domain_t),              intent(out) :: this
    type(simul_box_t), optional, intent(in)  :: sb
    type(geometry_t),  optional, intent(in)  :: geo
    !
    PUSH_SUB(domain__start__)
    if(present(sb).and.present(geo))then
      call domain__istart__(this, sb, geo)
    else
      if((.not.associated(this%sb)).and.(.not.associated(this%geo)))then
        ASSERT(associated(this%prnt))
        ASSERT(associated(this%prnt%sb))
        ASSERT(associated(this%prnt%geo))
        call domain__istart__(this, this%prnt%sb, this%prnt%geo)
      end if
    end if
    POP_SUB(domain__start__)
    return
  end subroutine domain__start__

  ! ---------------------------------------------------------
  subroutine domain_start(this, sb, geo)
    type(domain_t),              intent(out) :: this
    type(simul_box_t), optional, intent(in)  :: sb
    type(geometry_t),  optional, intent(in)  :: geo
    !
    PUSH_SUB(domain_start)
    call domain__start__(this, sb, geo)
    POP_SUB(domain_start)
    return
  end subroutine domain_start

  ! ---------------------------------------------------------
  subroutine domain__add__(this, that, config)
    type(domain_t),      intent(inout) :: this
    type(domain_t),      intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(domain__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call domain_hash_set(this%hash, config, that)
    POP_SUB(domain__add__)
    return
  end subroutine domain__add__

  ! ---------------------------------------------------------
  subroutine domain_get_domain_by_config(this, config, that)
    type(domain_t),      intent(in) :: this
    type(json_object_t), intent(in) :: config
    type(domain_t),     pointer     :: that
    !
    integer :: ierr
    !
    PUSH_SUB(domain_get_domain_by_config)
    nullify(that)
    call domain_hash_get(this%hash, config, that, ierr)
    if(ierr/=DOMAIN_OK)nullify(that)
    POP_SUB(domain_get_domain_by_config)
    return
  end subroutine domain_get_domain_by_config

  ! ---------------------------------------------------------
  subroutine domain_get_domain_by_name(this, name, that)
    type(domain_t),   intent(in) :: this
    character(len=*), intent(in) :: name
    type(domain_t),  pointer     :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(domain_get_domain_by_name)
    nullify(config, that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)&
      call domain_get(this, config, that)
    POP_SUB(domain_get_domain_by_name)
    return
  end subroutine domain_get_domain_by_name

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
  subroutine domain__copy__(this, that)
    type(domain_t), intent(inout) :: this
    type(domain_t), intent(in)    :: that
    !
    PUSH_SUB(domain__copy__)
    call domain__end__(this)
    call domain__init__(this)
    if(associated(that%sb).and.associated(that%geo))&
      call domain__start__(this, that%sb, that%geo)
    POP_SUB(domain__copy__)
    return
  end subroutine domain__copy__

  ! ---------------------------------------------------------
  recursive subroutine domain_copy_domain(this, that)
    type(domain_t), intent(inout) :: this
    type(domain_t), intent(in)    :: that
    !
    type(domain_iterator_t)      :: iter
    type(domain_t),      pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(domain_copy_domain)
    nullify(cnfg, osub, isub)
    call domain_end(this)
    call domain__copy__(this, that)
    call domain_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call domain_next(iter, cnfg, isub, ierr)
      if(ierr/=DOMAIN_OK)exit
      call domain_new(this, osub)
      call domain_copy(osub, isub)
      call domain__add__(this, osub, cnfg)
    end do
    call domain_end(iter)
    nullify(cnfg, osub, isub)
    POP_SUB(domain_copy_domain)
    return
  end subroutine domain_copy_domain

  ! ---------------------------------------------------------
  subroutine domain__end__(this)
    type(domain_t), intent(inout) :: this
    !
    PUSH_SUB(domain__end__)
    call domain__inull__(this)
    call config_dict_end(this%dict)
    call domain_hash_end(this%hash)
    call domain_list_end(this%list)
    POP_SUB(domain__end__)
    return
  end subroutine domain__end__

  ! ---------------------------------------------------------
  recursive subroutine domain_end_domain(this)
    type(domain_t), intent(inout) :: this
    !
    type(domain_t), pointer :: subs
    !
    PUSH_SUB(domain_end_domain)
    do
      nullify(subs)
      call domain_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call domain_end(subs)
      call domain__idel__(subs)
    end do
    nullify(subs)
    call domain__end__(this)
    POP_SUB(domain_end_domain)
    return
  end subroutine domain_end_domain

#define TEMPLATE_PREFIX domain
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module domain_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
