#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

#define BASE_LEAF_TYPE

module base_term_oct_m

  use base_system_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use message_oct_m
  use messages_oct_m
  use msgbus_oct_m
  use profiling_oct_m

#define BASE_TEMPLATE_NAME base_term
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::      &
    base_term_t

  public ::              &
    base_term__init__,   &
    base_term__update__, &
    base_term__reset__,  &
    base_term__copy__,   &
    base_term__end__

  public ::           &
    base_term_new,    &
    base_term_del,    &
    base_term_init,   &
    base_term_acc,    &
    base_term_calc,   &
    base_term_update, &
    base_term_set,    &
    base_term_get,    &
    base_term_copy,   &
    base_term_end

#define BASE_TEMPLATE_NAME base_term
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_term_t
    private
    type(json_object_t), pointer :: config =>null()
    type(base_system_t), pointer :: sys    =>null()
    type(refcount_t),    pointer :: rcnt   =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(msgbus_t)               :: msgb
    type(base_term_dict_t)       :: dict
  end type base_term_t

  interface base_term__init__
    module procedure base_term__init__type
    module procedure base_term__init__copy
  end interface base_term__init__

  interface base_term_new
    module procedure base_term_new_type
    module procedure base_term_new_pass
  end interface base_term_new

  interface base_term_init
    module procedure base_term_init_type
  end interface base_term_init

  interface base_term_set
    module procedure base_term_set_info
  end interface base_term_set

  interface base_term_get
    module procedure base_term_get_info
    module procedure base_term_get_config
    module procedure base_term_get_system
    module procedure base_term_get_msgbus
  end interface base_term_get

contains

#define BASE_TEMPLATE_NAME base_term
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_term_new_type(sys, config) result(this)
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    type(base_term_t), pointer :: this

    PUSH_SUB(base_term_new_type)

    this => base_term_new(sys, config, base_term_init_type)

    POP_SUB(base_term_new_type)
  end function base_term_new_type

  ! ---------------------------------------------------------
  function base_term_new_pass(sys, config, init) result(this)
    type(base_system_t), intent(in) :: sys
    type(json_object_t), intent(in) :: config

    type(base_term_t), pointer :: this

    interface
      subroutine init(this, sys, config)
        use json_oct_m
        use base_system_oct_m
        import :: base_term_t
        type(base_term_t),   intent(out) :: this
        type(base_system_t), intent(in)  :: sys
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface

    PUSH_SUB(base_term_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, sys, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_term_new_pass)
  end function base_term_new_pass

  ! ---------------------------------------------------------
  subroutine base_term__init__type(this, sys, config)
    type(base_term_t),           intent(out) :: this
    type(base_system_t), target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_term__init__type)

    this%config => config
    this%sys => sys
    this%rcnt => refcount_new()
    call msgbus_init(this%msgb, number=2)
    call base_term_dict_init(this%dict)

    POP_SUB(base_term__init__type)
  end subroutine base_term__init__type

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
  subroutine base_term_init_type(this, sys, config)
    type(base_term_t),   intent(out) :: this
    type(base_system_t), intent(in)  :: sys
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_term_init_type)

    call base_term__init__(this, sys, config)

    POP_SUB(base_term_init_type)
  end subroutine base_term_init_type

  ! ---------------------------------------------------------
  subroutine base_term__acc__(this, that, config)
    type(base_term_t),             intent(inout) :: this
    type(base_term_t),             intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_term__acc__)

    ASSERT(associated(this%config))
    if(present(config)) continue
    this%energy = this%energy + that%energy

    POP_SUB(base_term__acc__)
  end subroutine base_term__acc__

  ! ---------------------------------------------------------
  subroutine base_term__sub__(this, that, config)
    type(base_term_t),             intent(inout) :: this
    type(base_term_t),             intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    PUSH_SUB(base_term__sub__)

    ASSERT(associated(this%config))
    if(present(config)) continue
    this%energy = this%energy - that%energy

    POP_SUB(base_term__sub__)
  end subroutine base_term__sub__

  ! ---------------------------------------------------------
  subroutine base_term__update__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__update__)

    ASSERT(associated(this%config))

    POP_SUB(base_term__update__)
  end subroutine base_term__update__

  ! ---------------------------------------------------------
  subroutine base_term__reset__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__reset__)

    ASSERT(associated(this%config))
    this%energy = 0.0_wp

    POP_SUB(base_term__reset__)
  end subroutine base_term__reset__

  ! ---------------------------------------------------------
  subroutine base_term_acc(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term_acc)

    ASSERT(base_term_dict_len(this%dict)>0)
    call base_term__reset__(this)
    call base_term__reduce__(this, base_term__acc__)
    call base_term__update__(this)
    
    POP_SUB(base_term_acc)
  end subroutine base_term_acc

  ! ---------------------------------------------------------
  subroutine base_term_calc(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term_calc)

    call base_term_acc(this)
    
    POP_SUB(base_term_calc)
  end subroutine base_term_calc

  ! ---------------------------------------------------------
  subroutine base_term_update(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term_update)

    call base_term__apply__(this, base_term__update__)

    POP_SUB(base_term_update)
  end subroutine base_term_update

  ! ---------------------------------------------------------
  subroutine base_term_reset(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term_reset)

    call base_term__apply__(this, base_term__reset__)

    POP_SUB(base_term_reset)
  end subroutine base_term_reset

  ! ---------------------------------------------------------
  subroutine base_term_set_info(this, energy)
    type(base_term_t),       intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: energy

    PUSH_SUB(base_term_set_info)

    ASSERT(associated(this%config))
    if(present(energy)) this%energy = energy

    POP_SUB(base_term_set_info)
  end subroutine base_term_set_info

  ! ---------------------------------------------------------
  subroutine base_term_get_info(this, energy)
    type(base_term_t),       intent(in)  :: this
    real(kind=wp), optional, intent(out) :: energy

    PUSH_SUB(base_term_get_info)

    ASSERT(associated(this%config))
    if(present(energy)) energy = this%energy

    POP_SUB(base_term_get_info)
  end subroutine base_term_get_info

  ! ---------------------------------------------------------
  subroutine base_term_get_config(this, that)
    type(base_term_t),    target, intent(in)  :: this
    type(json_object_t), pointer, intent(out) :: that

    PUSH_SUB(base_term_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_term_get_config)
  end subroutine base_term_get_config

  ! ---------------------------------------------------------
  subroutine base_term_get_system(this, that)
    type(base_term_t),    target, intent(in)  :: this
    type(base_system_t), pointer, intent(out) :: that

    PUSH_SUB(base_term_get_system)

    nullify(that)
    if(associated(this%sys)) that => this%sys

    POP_SUB(base_term_get_system)
  end subroutine base_term_get_system

  ! ---------------------------------------------------------
  subroutine base_term_get_msgbus(this, that)
    type(base_term_t), target, intent(in)  :: this
    type(msgbus_t),   pointer, intent(out) :: that

    PUSH_SUB(base_term_get_msgbus)

    nullify(that)
    if(associated(this%config)) that => this%msgb

    POP_SUB(base_term_get_msgbus)
  end subroutine base_term_get_msgbus

  ! ---------------------------------------------------------
  subroutine base_term__copy__(this, that)
    type(base_term_t), intent(inout) :: this
    type(base_term_t), intent(in)    :: that

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_term__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call  base_term__end__(this)
    if(associated(that%config).and.associated(that%sys))then
      call base_term__init__(this, that)
      call refcount_del(this%rcnt)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_term__copy__)
  end subroutine base_term__copy__

  ! ---------------------------------------------------------
  subroutine base_term__end__(this)
    type(base_term_t), intent(inout) :: this

    PUSH_SUB(base_term__end__)

    nullify(this%config, this%sys)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    this%energy = 0.0_wp
    call msgbus_end(this%msgb)
    ASSERT(base_term_dict_len(this%dict)==0)
    call base_term_dict_end(this%dict)

    POP_SUB(base_term__end__)
  end subroutine base_term__end__

end module base_term_oct_m

#undef BASE_LEAF_TYPE

!! Local Variables:
!! mode: f90
!! End:
