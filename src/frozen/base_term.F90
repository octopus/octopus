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

#define HASH_TEMPLATE_NAME base_term
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_term

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_term_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

  use config_dict_m, only: &
    CONFIG_DICT_OK,        &
    CONFIG_DICT_NAME_LEN

  use config_dict_m, only: &
    config_dict_t,         &
    config_dict_init,      &
    config_dict_set,       &
    config_dict_get,       &
    config_dict_copy,      &
    config_dict_end

  use base_system_m, only: &
    base_system_t

  implicit none

  private
  public ::              &
    base_term__update__, &
    base_term__add__

  public ::               &
    base_term_init,       &
    base_term_update,     &
    base_term_next,       &
    base_term_set,        &
    base_term_get,        &
    base_term_get_energy, &
    base_term_set_energy, &
    base_term_copy,       &
    base_term_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_term_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_system_t), pointer :: sys    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(config_dict_t)          :: dict
    type(base_term_hash_t)       :: hash
  end type base_term_t

  type, public :: base_term_iterator_t
    private
    type(base_term_t),      pointer :: self =>null()
    type(base_term_hash_iterator_t) :: iter
  end type base_term_iterator_t

  interface base_term_init
    module procedure base_term_init_term
    module procedure base_term_iterator_init
  end interface base_term_init

  interface base_term_next
    module procedure base_term_iterator_next_config_term
    module procedure base_term_iterator_next_config
    module procedure base_term_iterator_next_term
  end interface base_term_next

  interface base_term_set
    module procedure base_term_set_info
  end interface base_term_set

  interface base_term_get
    module procedure base_term_get_info
    module procedure base_term_get_config
    module procedure base_term_get_system
  end interface base_term_get

  interface base_term_copy
    module procedure base_term_copy_term
    module procedure base_term_iterator_copy
  end interface base_term_copy

  interface base_term_end
    module procedure base_term_end_term
    module procedure base_term_iterator_end
  end interface base_term_end

  integer, public, parameter :: BASE_TERM_OK          = BASE_TERM_HASH_OK
  integer, public, parameter :: BASE_TERM_KEY_ERROR   = BASE_TERM_HASH_KEY_ERROR
  integer, public, parameter :: BASE_TERM_EMPTY_ERROR = BASE_TERM_HASH_EMPTY_ERROR

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_term_init_term(this, sys, config)
    type(base_term_t),           intent(out) :: this
    type(base_system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(base_term_init_term)
    ASSERT(.not.associated(this%config))
    this%config=>config
    this%sys=>sys
    this%energy=0.0_wp
    call config_dict_init(this%dict)
    call base_term_hash_init(this%hash)
    POP_SUB(base_term_init_term)
    return
  end subroutine base_term_init_term

  ! ---------------------------------------------------------
  subroutine base_term__update__(this)
    type(base_term_t), intent(inout) :: this
    !
    PUSH_SUB(base_term__update__)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sys))
    POP_SUB(base_term__update__)
    return
  end subroutine base_term__update__

  ! ---------------------------------------------------------
  recursive subroutine base_term_update(this)
    type(base_term_t), intent(inout) :: this
    !
    type(base_term_iterator_t) :: iter
    type(base_term_t), pointer :: subs
    integer                    :: ierr
    !
    PUSH_SUB(base_term_update)
    nullify(subs)
    this%energy=0.0_wp
    call base_term_init(iter, this)
    do
      nullify(subs)
      call base_term_next(iter, subs, ierr)
      if(ierr/=BASE_TERM_OK)exit
      call base_term_update(subs)
      this%energy=this%energy+subs%energy
    end do
    call base_term_end(iter)
    nullify(subs)
    call base_term__update__(this)
    POP_SUB(base_term_update)
    return
  end subroutine base_term_update

  ! ---------------------------------------------------------
  subroutine base_term__add__(this, that, config)
    type(base_term_t),   intent(inout) :: this
    type(base_term_t),   intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    character(len=CONFIG_DICT_NAME_LEN) :: name
    integer                             :: ierr
    !
    PUSH_SUB(base_term_init__add__)
    call json_get(config, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call config_dict_set(this%dict, trim(adjustl(name)), config)
    call base_term_hash_set(this%hash, config, that)
    POP_SUB(base_term_init__add__)
    return
  end subroutine base_term__add__

  ! ---------------------------------------------------------
  subroutine base_term__get__(this, name, that)
    type(base_term_t),  intent(inout) :: this
    character(len=*),   intent(in)    :: name
    type(base_term_t), pointer        :: that
    !
    type(json_object_t), pointer :: config
    integer                      :: ierr
    !
    PUSH_SUB(base_term__get__)
    nullify(that)
    call config_dict_get(this%dict, trim(adjustl(name)), config, ierr)
    if(ierr==CONFIG_DICT_OK)then
      call base_term_hash_get(this%hash, config, that, ierr)
      if(ierr/=BASE_TERM_OK)nullify(that)
    end if
    POP_SUB(base_term__get__)
    return
  end subroutine base_term__get__

  ! ---------------------------------------------------------
  elemental subroutine base_term_set_energy(this, that)
    type(base_term_t), intent(inout) :: this
    real(kind=wp),     intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine base_term_set_energy

  ! ---------------------------------------------------------
  subroutine base_term_set_info(this, energy)
    type(base_term_t),       intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy
    !
    PUSH_SUB(base_term_set_info)
    if(present(energy))this%energy=energy
    POP_SUB(base_term_set_info)
    return
  end subroutine base_term_set_info

  ! ---------------------------------------------------------
  elemental function base_term_get_energy(this) result(that)
    type(base_term_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function base_term_get_energy

  ! ---------------------------------------------------------
  subroutine base_term_get_info(this, energy)
    type(base_term_t),       intent(in)  :: this
    real(kind=wp), optional, intent(out) :: energy
    !
    PUSH_SUB(base_term_get_info)
    if(present(energy))energy=this%energy
    POP_SUB(base_term_get_info)
    return
  end subroutine base_term_get_info

  ! ---------------------------------------------------------
  subroutine base_term_get_config(this, that)
    type(base_term_t),    target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_term_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_term_get_config)
    return
  end subroutine base_term_get_config

  ! ---------------------------------------------------------
  subroutine base_term_get_system(this, that)
    type(base_term_t),    target, intent(in) :: this
    type(base_system_t), pointer             :: that
    !
    PUSH_SUB(base_term_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(base_term_get_system)
    return
  end subroutine base_term_get_system

  ! ---------------------------------------------------------
  subroutine base_term_copy_term(this, that)
    type(base_term_t), intent(out) :: this
    type(base_term_t), intent(in)  :: that
    !
    PUSH_SUB(base_term_copy_term)
    this%config=>that%config
    this%sys=>that%sys
    this%energy=that%energy
    call config_dict_copy(this%dict, that%dict)
    call base_term_hash_copy(this%hash, that%hash)
    POP_SUB(base_term_copy_term)
    return
  end subroutine base_term_copy_term

  ! ---------------------------------------------------------
  subroutine base_term_end_term(this)
    type(base_term_t), intent(inout) :: this
    !
    PUSH_SUB(base_term_end_term)
    nullify(this%config, this%sys)
    this%energy=0.0_wp
    call config_dict_end(this%dict)
    call base_term_hash_end(this%hash)
    POP_SUB(base_term_end_term)
    return
  end subroutine base_term_end_term

  ! ---------------------------------------------------------
  subroutine base_term_iterator_init(this, that)
    type(base_term_iterator_t), intent(out) :: this
    type(base_term_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_term_iterator_init)
    This%self=>that
    call base_term_hash_init(this%iter, that%hash)
    POP_SUB(base_term_iterator_init)
    return
  end subroutine base_term_iterator_init

  ! ---------------------------------------------------------
  subroutine base_term_iterator_next_config_term(this, config, that, ierr)
    type(base_term_iterator_t), intent(inout) :: this
    type(json_object_t),       pointer        :: config
    type(base_term_t),         pointer        :: that
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_term_iterator_next_config_term)
    call base_term_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_term_iterator_next_config_term)
    return
  end subroutine base_term_iterator_next_config_term

  ! ---------------------------------------------------------
  subroutine base_term_iterator_next_config(this, that, ierr)
    type(base_term_iterator_t), intent(inout) :: this
    type(json_object_t),       pointer        :: that
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_term_iterator_next_config)
    call base_term_hash_next(this%iter, that, ierr)
    POP_SUB(base_term_iterator_next_config)
    return
  end subroutine base_term_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_term_iterator_next_term(this, that, ierr)
    type(base_term_iterator_t), intent(inout) :: this
    type(base_term_t),         pointer        :: that
    integer,          optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_term_iterator_next_term)
    call base_term_hash_next(this%iter, that, ierr)
    POP_SUB(base_term_iterator_next_term)
    return
  end subroutine base_term_iterator_next_term

  ! ---------------------------------------------------------
  subroutine base_term_iterator_copy(this, that)
    type(base_term_iterator_t), intent(inout) :: this
    type(base_term_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_term_iterator_copy)
    this%self=>that%self
    call base_term_hash_copy(this%iter, that%iter)
    POP_SUB(base_term_iterator_copy)
    return
  end subroutine base_term_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_term_iterator_end(this)
    type(base_term_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_term_iterator_end)
    nullify(this%self)
    call base_term_hash_end(this%iter)
    POP_SUB(base_term_iterator_end)
    return
  end subroutine base_term_iterator_end

end module base_term_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
