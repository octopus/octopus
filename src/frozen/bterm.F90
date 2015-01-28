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

#define HASH_TEMPLATE_NAME bterm
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bterm

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bterm_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

  use bsyst_m, only:     &
    system_t => bsyst_t

  implicit none

  private
  public ::           &
    bterm_init,       &
    bterm_get,        &
    bterm_get_energy, &
    bterm_set_energy, &
    bterm_copy,       &
    bterm_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bterm_t
    private
    type(json_object_t), pointer :: config =>null()
    type(system_t),      pointer :: sys    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(bterm_hash_t)           :: hash
  end type bterm_t

  type, public :: bterm_iterator_t
    private
    type(bterm_t),      pointer :: self =>null()
    type(bterm_hash_iterator_t) :: iter
  end type bterm_iterator_t

  interface bterm_init
    module procedure bterm_init_start
    module procedure bterm_init_build
    module procedure bterm_iterator_init
  end interface bterm_init

  interface bterm_next
    module procedure bterm_iterator_next_config_bterm
    module procedure bterm_iterator_next_config
    module procedure bterm_iterator_next_bterm
  end interface bterm_next

  interface bterm_get
    module procedure bterm_get_config
    module procedure bterm_get_system
  end interface bterm_get

  interface bterm_copy
    module procedure bterm_copy_bterm
    module procedure bterm_iterator_copy
  end interface bterm_copy

  interface bterm_end
    module procedure bterm_end_bterm
    module procedure bterm_iterator_end
  end interface bterm_end

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bterm_init_start(this, sys, config)
    type(bterm_t),               intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(bterm_init_start)
    ASSERT(.not.associated(this%config))
    this%config=>config
    this%sys=>sys
    this%energy=0.0_wp
    call bterm_hash_init(this%hash)
    POP_SUB(bterm_init_start)
    return
  end subroutine bterm_init_start

  ! ---------------------------------------------------------
  subroutine bterm_init_build(this, that, config)
    type(bterm_t),       intent(out) :: this
    type(bterm_t),       intent(in)  :: that
    type(json_object_t), intent(in)  :: config
    !
    PUSH_SUB(bterm_init_build)
    call bterm_hash_set(this%hash, config, that)
    POP_SUB(bterm_init_build)
    return
  end subroutine bterm_init_build

 ! ---------------------------------------------------------
  elemental function bterm_get_energy(this) result(that)
    type(bterm_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function bterm_get_energy

  ! ---------------------------------------------------------
  subroutine bterm_get_config(this, that)
    type(bterm_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bterm_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bterm_get_config)
    return
  end subroutine bterm_get_config

  ! ---------------------------------------------------------
  subroutine bterm_get_system(this, that)
    type(bterm_t),   target, intent(in) :: this
    type(system_t), pointer             :: that
    !
    PUSH_SUB(bterm_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(bterm_get_system)
    return
  end subroutine bterm_get_system

  ! ---------------------------------------------------------
  elemental subroutine bterm_set_energy(this, that)
    type(bterm_t), intent(inout) :: this
    real(kind=wp), intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine bterm_set_energy

  ! ---------------------------------------------------------
  subroutine bterm_copy_bterm(this, that)
    type(bterm_t), intent(out) :: this
    type(bterm_t), intent(in)  :: that
    !
    PUSH_SUB(bterm_copy_bterm)
    this%config=>that%config
    this%sys=>that%sys
    this%energy=that%energy
    call bterm_hash_copy(this%hash, that%hash)
    POP_SUB(bterm_copy_bterm)
    return
  end subroutine bterm_copy_bterm

  ! ---------------------------------------------------------
  subroutine bterm_end_bterm(this)
    type(bterm_t), intent(inout) :: this
    !
    PUSH_SUB(bterm_end_bterm)
    nullify(this%config, this%sys)
    this%energy=0.0_wp
    call bterm_hash_end(this%hash)
    POP_SUB(bterm_end_bterm)
    return
  end subroutine bterm_end_bterm

  ! ---------------------------------------------------------
  subroutine bterm_iterator_init(this, that)
    type(bterm_iterator_t), intent(out) :: this
    type(bterm_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bterm_iterator_init)
    This%self=>that
    call bterm_hash_init(this%iter, that%hash)
    POP_SUB(bterm_iterator_init)
    return
  end subroutine bterm_iterator_init

  ! ---------------------------------------------------------
  subroutine bterm_iterator_next_config_bterm(this, config, term, ierr)
    type(bterm_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bterm_t),         pointer        :: term
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bterm_iterator_next_config_bterm)
    call bterm_hash_next(this%iter, config, term, ierr)
    POP_SUB(bterm_iterator_next_config_bterm)
    return
  end subroutine bterm_iterator_next_config_bterm

  ! ---------------------------------------------------------
  subroutine bterm_iterator_next_config(this, that, ierr)
    type(bterm_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bterm_iterator_next_config)
    call bterm_hash_next(this%iter, that, ierr)
    POP_SUB(bterm_iterator_next_config)
    return
  end subroutine bterm_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bterm_iterator_next_bterm(this, that, ierr)
    type(bterm_iterator_t), intent(inout) :: this
    type(bterm_t),         pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bterm_iterator_next_bterm)
    call bterm_hash_next(this%iter, that, ierr)
    POP_SUB(bterm_iterator_next_bterm)
    return
  end subroutine bterm_iterator_next_bterm

  ! ---------------------------------------------------------
  subroutine bterm_iterator_copy(this, that)
    type(bterm_iterator_t), intent(inout) :: this
    type(bterm_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bterm_iterator_copy)
    this%self=>that%self
    call bterm_hash_copy(this%iter, that%iter)
    POP_SUB(bterm_iterator_copy)
    return
  end subroutine bterm_iterator_copy

  ! ---------------------------------------------------------
  subroutine bterm_iterator_end(this)
    type(bterm_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bterm_iterator_end)
    nullify(this%self)
    call bterm_hash_end(this%iter)
    POP_SUB(bterm_iterator_end)
    return
  end subroutine bterm_iterator_end

end module bterm_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
