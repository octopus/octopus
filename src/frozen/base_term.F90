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

#define HASH_TEMPLATE_NAME base_term
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_term

module base_term_m

  use base_system_m
  use config_dict_m
  use global_m
  use json_m
  use kinds_m
  use messages_m
  use profiling_m

#define LIST_TEMPLATE_NAME base_term
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_PREFIX
#include "thash_inc.F90"
#undef HASH_INCLUDE_PREFIX

#define TEMPLATE_PREFIX base_term
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                &
    BASE_TERM_OK,          &
    BASE_TERM_KEY_ERROR,   &
    BASE_TERM_EMPTY_ERROR

  public ::      &
    base_term_t

  public ::              &
    base_term__init__,   &
    base_term__update__, &
    base_term__reset__,  &
    base_term__acc__,    &
    base_term__sub__,    &
    base_term__copy__,   &
    base_term__end__

  public ::           &
    base_term_new,    &
    base_term_del,    &
    base_term_init,   &
    base_term_update, &
    base_term_sets,   &
    base_term_gets,   &
    base_term_set,    &
    base_term_get,    &
    base_term_copy,   &
    base_term_end

#define LIST_TEMPLATE_NAME base_term
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_HEADER
#include "thash_inc.F90"
#undef HASH_INCLUDE_HEADER

  integer, parameter :: BASE_TERM_OK          = BASE_TERM_HASH_OK
  integer, parameter :: BASE_TERM_KEY_ERROR   = BASE_TERM_HASH_KEY_ERROR
  integer, parameter :: BASE_TERM_EMPTY_ERROR = BASE_TERM_HASH_EMPTY_ERROR

  type :: base_term_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_system_t), pointer :: sys    =>null()
    type(base_term_t),   pointer :: prnt   =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(config_dict_t)          :: dict
    type(base_term_hash_t)       :: hash
    type(base_term_list_t)       :: list
  end type base_term_t

  interface base_term__init__
    module procedure base_term__init__term
    module procedure base_term__init__copy
  end interface base_term__init__

  interface base_term_init
    module procedure base_term_init_term
    module procedure base_term_init_copy
  end interface base_term_init

  interface base_term_set
    module procedure base_term_set_info
  end interface base_term_set

  interface base_term_gets
    module procedure base_term_gets_config
    module procedure base_term_gets_name
  end interface base_term_gets

  interface base_term_get
    module procedure base_term_get_info
    module procedure base_term_get_config
    module procedure base_term_get_system
  end interface base_term_get

  interface base_term_copy
    module procedure base_term_copy_term
  end interface base_term_copy

  interface base_term_end
    module procedure base_term_end_term
  end interface base_term_end

#define TEMPLATE_PREFIX base_term
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_term
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define HASH_INCLUDE_BODY
#include "thash_inc.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_term_new(this, that)
    type(base_term_t),  target, intent(inout) :: this
    type(base_term_t), pointer                :: that

    PUSH_SUB(base_term_new)

    nullify(that)
    SAFE_ALLOCATE(that)
    that%prnt => this
    call base_term_list_push(this%list, that)

    POP_SUB(base_term_new)
  end subroutine base_term_new

  ! ---------------------------------------------------------
  subroutine base_term__idel__(this)
    type(base_term_t), pointer :: this

    PUSH_SUB(base_term__idel__)

    SAFE_DEALLOCATE_P(this)
    nullify(this)

    POP_SUB(base_term__idel__)
  end subroutine base_term__idel__

  ! ---------------------------------------------------------
  subroutine base_term_del(this)
    type(base_term_t), pointer :: this

    PUSH_SUB(base_term_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_term_list_del(this%prnt%list, this)
        call base_term_end(this)
        call base_term__idel__(this)
      end if
    end if

    POP_SUB(base_term_del)
  end subroutine base_term_del

  ! ---------------------------------------------------------
  subroutine base_term__inull__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__inull__)

    nullify(this%config, this%sys, this%prnt)
    this%energy = 0.0_wp

    POP_SUB(base_term__inull__)
  end subroutine base_term__inull__

  ! ---------------------------------------------------------
  subroutine base_term__init__term(this, sys, config)
    type(base_term_t),           intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_term__init__term)

    call base_term__inull__(this)
    this%config => config
    this%sys => sys
    call config_dict_init(this%dict)
    call base_term_hash_init(this%hash)
    call base_term_list_init(this%list)

    POP_SUB(base_term__init__term)
  end subroutine base_term__init__term

  ! ---------------------------------------------------------
  subroutine base_term__init__copy(this, that)
    type(base_term_t), intent(out) :: this
    type(base_term_t), intent(in)  :: that

    PUSH_SUB(base_term__init__copy)

    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    call base_term__init__(this, that%sys, that%config)

    POP_SUB(base_term__init__copy)
  end subroutine base_term__init__copy

  ! ---------------------------------------------------------
  subroutine base_term_init_term(this, sys, config)
    type(base_term_t),   intent(out) :: this
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_term_init_term)

    call base_term__init__(this, sys, config)

    POP_SUB(base_term_init_term)
  end subroutine base_term_init_term

  ! ---------------------------------------------------------
  recursive subroutine base_term_init_copy(this, that)
    type(base_term_t), intent(out) :: this
    type(base_term_t), intent(in)  :: that

    type(base_term_iterator_t)   :: iter
    type(base_term_t),   pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_term_init_copy)

    nullify(cnfg, osub, isub)
    call base_term__init__(this, that)
    call base_term_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_term_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_TERM_OK)exit
      call base_term_new(this, osub)
      call base_term_init(osub, isub)
      call base_term_sets(this, osub, cnfg)
    end do
    call base_term_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_term_init_copy)
  end subroutine base_term_init_copy

  ! ---------------------------------------------------------
  subroutine base_term__update__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__update__)

    ASSERT(associated(this%config))

    POP_SUB(base_term__update__)
  end subroutine base_term__update__

  ! ---------------------------------------------------------
  recursive subroutine base_term_update(this)
    type(base_term_t), intent(inout) :: this

    type(base_term_iterator_t) :: iter
    type(base_term_t), pointer :: subs
    integer                    :: ierr

    PUSH_SUB(base_term_update)

    nullify(subs)
    call base_term_init(iter, this)
    do
      nullify(subs)
      call base_term_next(iter, subs, ierr)
      if(ierr/=BASE_TERM_OK)exit
      call base_term_update(subs)
    end do
    call base_term_end(iter)
    nullify(subs)
    call base_term__update__(this)

    POP_SUB(base_term_update)
  end subroutine base_term_update

  ! ---------------------------------------------------------
  subroutine base_term__reset__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__reset__)

    ASSERT(associated(this%config))
    this%energy = 0.0_wp

    POP_SUB(base_term__reset__)
  end subroutine base_term__reset__

  ! ---------------------------------------------------------
  subroutine base_term__acc__(this, that)
    type(base_term_t), intent(inout) :: this
    type(base_term_t), intent(in)    :: that

    PUSH_SUB(base_term__acc__)

    ASSERT(associated(this%config))
    this%energy = this%energy + that%energy

    POP_SUB(base_term__acc__)
  end subroutine base_term__acc__

  ! ---------------------------------------------------------
  subroutine base_term__sub__(this, that)
    type(base_term_t), intent(inout) :: this
    type(base_term_t), intent(in)    :: that

    PUSH_SUB(base_term__sub__)

    ASSERT(associated(this%config))
    this%energy = this%energy - that%energy

    POP_SUB(base_term__sub__)
  end subroutine base_term__sub__

  ! ---------------------------------------------------------
  subroutine base_term_sets(this, that, config)
    type(base_term_t),   intent(inout) :: this
    type(base_term_t),   intent(in)    :: that
    type(json_object_t), intent(in)    :: config

    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr

    PUSH_SUB(base_term_init_sets)

    ASSERT(associated(this%config))
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_term_hash_set(this%hash, config, that)

    POP_SUB(base_term_init_sets)
  end subroutine base_term_sets

  ! ---------------------------------------------------------
  subroutine base_term_gets_config(this, config, that)
    type(base_term_t),   intent(in) :: this
    type(json_object_t), intent(in) :: config
    type(base_term_t),  pointer     :: that

    integer :: ierr

    PUSH_SUB(base_term_gets_config)

    nullify(that)
    ASSERT(associated(this%config))
    call base_term_hash_get(this%hash, config, that, ierr)
    if(ierr/=BASE_TERM_OK)nullify(that)

    POP_SUB(base_term_gets_config)
  end subroutine base_term_gets_config

  ! ---------------------------------------------------------
  subroutine base_term_gets_name(this, name, that)
    type(base_term_t),  intent(in) :: this
    character(len=*),   intent(in) :: name
    type(base_term_t), pointer     :: that

    type(json_object_t), pointer :: config
    integer                      :: ierr

    PUSH_SUB(base_term_gets_name)

    nullify(that)
    ASSERT(associated(this%config))
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK) call base_term_gets(this, config, that)

    POP_SUB(base_term_gets_name)
  end subroutine base_term_gets_name

  ! ---------------------------------------------------------
  subroutine base_term_set_info(this, energy)
    type(base_term_t),       intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_term_set_info)

    if(present(energy)) this%energy = energy

    POP_SUB(base_term_set_info)
  end subroutine base_term_set_info

  ! ---------------------------------------------------------
  subroutine base_term_get_info(this, energy)
    type(base_term_t),       intent(in)  :: this
    real(kind=wp), optional, intent(out) :: energy

    PUSH_SUB(base_term_get_info)

    if(present(energy)) energy = this%energy

    POP_SUB(base_term_get_info)
  end subroutine base_term_get_info

  ! ---------------------------------------------------------
  subroutine base_term_get_config(this, that)
    type(base_term_t),    target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_term_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_term_get_config)
  end subroutine base_term_get_config

  ! ---------------------------------------------------------
  subroutine base_term_get_system(this, that)
    type(base_term_t),    target, intent(in) :: this
    type(base_system_t), pointer             :: that

    PUSH_SUB(base_term_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_term_get_system)
  end subroutine base_term_get_system

  ! ---------------------------------------------------------
  subroutine base_term__copy__(this, that)
    type(base_term_t), intent(inout) :: this
    type(base_term_t), intent(in)    :: that

    type(base_term_t), pointer :: prnt

    PUSH_SUB(base_term__copy__)

    prnt => this%prnt
    call  base_term__end__(this)
    if(associated(that%config).and.associated(that%sys))&
      call base_term__init__(this, that%sys, that%config)
    this%prnt => prnt

    POP_SUB(base_term__copy__)
  end subroutine base_term__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_term_copy_term(this, that)
    type(base_term_t), intent(inout) :: this
    type(base_term_t), intent(in)    :: that

    type(base_term_iterator_t)   :: iter
    type(base_term_t),   pointer :: osub, isub
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_term_copy_term)

    nullify(cnfg, osub, isub)
    call base_term_end(this)
    call base_term__copy__(this, that)
    call base_term_init(iter, that)
    do
      nullify(cnfg, osub, isub)
      call base_term_next(iter, cnfg, isub, ierr)
      if(ierr/=BASE_TERM_OK)exit
      call base_term_new(this, osub)
      call base_term_copy(osub, isub)
      call base_term_sets(this, osub, cnfg)
    end do
    call base_term_end(iter)
    nullify(cnfg, osub, isub)

    POP_SUB(base_term_copy_term)
  end subroutine base_term_copy_term

  ! ---------------------------------------------------------
  subroutine base_term__end__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__end__)

    call base_term__inull__(this)
    call config_dict_end(this%dict)
    call base_term_hash_end(this%hash)
    call base_term_list_end(this%list)

    POP_SUB(base_term__end__)
  end subroutine base_term__end__

  ! ---------------------------------------------------------
  recursive subroutine base_term_end_term(this)
    type(base_term_t), intent(inout) :: this

    type(base_term_t), pointer :: subs

    PUSH_SUB(base_term_end_term)

    do
      nullify(subs)
      call base_term_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_term_end(subs)
      call base_term__idel__(subs)
    end do
    nullify(subs)
    call base_term__end__(this)

    POP_SUB(base_term_end_term)
  end subroutine base_term_end_term

#define TEMPLATE_PREFIX base_term
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_term_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
