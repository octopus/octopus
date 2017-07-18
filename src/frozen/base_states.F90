#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#undef DICT_TEMPLATE_NAME
#undef DICT_TYPE_NAME
#undef DICT_TYPE_MODULE_NAME
#undef DICT_INCLUDE_PREFIX
#undef DICT_INCLUDE_HEADER
#undef DICT_INCLUDE_BODY

module base_states_oct_m

  use base_density_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m

#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_states
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

#define TEMPLATE_PREFIX base_states
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                  &
    BASE_STATES_OK,          &
    BASE_STATES_KEY_ERROR,   &
    BASE_STATES_EMPTY_ERROR

  public ::        &
    base_states_t

  public ::                &
    base_states__init__,   &
    base_states__build__,  &
    base_states__start__,  &
    base_states__acc__,    &
    base_states__update__, &
    base_states__reset__,  &
    base_states__stop__,   &
    base_states__copy__,   &
    base_states__end__

  public ::             &
    base_states_new,    &
    base_states_del,    &
    base_states_init,   &
    base_states_start,  &
    base_states_acc,    &
    base_states_update, &
    base_states_reset,  &
    base_states_stop,   &
    base_states_set,    &
    base_states_get,    &
    base_states_copy,   &
    base_states_end

#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_states
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

  integer, parameter :: BASE_STATES_OK          = BASE_STATES_HASH_OK
  integer, parameter :: BASE_STATES_KEY_ERROR   = BASE_STATES_HASH_KEY_ERROR
  integer, parameter :: BASE_STATES_EMPTY_ERROR = BASE_STATES_HASH_EMPTY_ERROR

  type :: base_states_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_states_t), pointer :: prnt   =>null()
    type(base_density_t)         :: density
    type(base_states_dict_t)     :: dict
    type(base_states_list_t)     :: list
  end type base_states_t

  interface base_states__init__
    module procedure base_states__init__type
    module procedure base_states__init__copy
  end interface base_states__init__

  interface base_states_init
    module procedure base_states_init_type
    module procedure base_states_init_copy
  end interface base_states_init

  interface base_states_set
    module procedure base_states_set_sub
  end interface base_states_set

  interface base_states_get
    module procedure base_states_get_info
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_charge
    module procedure base_states_get_density
    module procedure base_states_get_density_1d
    module procedure base_states_get_density_2d
    module procedure base_states_get_sub
    module procedure base_states_get_sub_density
    module procedure base_states_get_sub_density_1d
    module procedure base_states_get_sub_density_2d
  end interface base_states_get

  interface base_states_copy
    module procedure base_states_copy_type
  end interface base_states_copy

  interface base_states_end
    module procedure base_states_end_type
  end interface base_states_end

#define TEMPLATE_PREFIX base_states
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains
    
#define LIST_TEMPLATE_NAME base_states
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_states
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_states__new__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_states__new__)
  end subroutine base_states__new__

  ! ---------------------------------------------------------
  subroutine base_states__del__(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_states__del__)
  end subroutine base_states__del__

  ! ---------------------------------------------------------
  subroutine base_states_new(this, that)
    type(base_states_t),  target, intent(inout) :: this
    type(base_states_t), pointer                :: that

    PUSH_SUB(base_states_new)

    nullify(that)
    call base_states__new__(that)
    that%prnt => this
    call base_states_list_push(this%list, that)

    POP_SUB(base_states_new)
  end subroutine base_states_new

  ! ---------------------------------------------------------
  subroutine base_states_del(this)
    type(base_states_t), pointer :: this

    PUSH_SUB(base_states_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_states_list_del(this%prnt%list, this)
        call base_states_end(this)
        call base_states__del__(this)
      end if
    end if
    nullify(this)

    POP_SUB(base_states_del)
  end subroutine base_states_del

  ! ---------------------------------------------------------
  subroutine base_states__iinit__(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    PUSH_SUB(base_states__iinit__)

    this%config => config
    call base_states_dict_init(this%dict)
    call base_states_list_init(this%list)

    POP_SUB(base_states__iinit__)
  end subroutine base_states__iinit__
    
  ! ---------------------------------------------------------
  subroutine base_states__init__type(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_states__init__type)

    nullify(cnfg)
    call base_states__iinit__(this, config)
    call json_get(this%config, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_density__init__(this%density, cnfg)
    nullify(cnfg)

    POP_SUB(base_states__init__type)
  end subroutine base_states__init__type
    
  ! ---------------------------------------------------------
  subroutine base_states__init__copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    PUSH_SUB(base_states__init__copy)

    ASSERT(associated(that%config))
    call base_states__iinit__(this, that%config)
    call base_density__init__(this%density, that%density)
    if(associated(that%sim)) call base_states__start__(this, that%sim)

    POP_SUB(base_states__init__copy)
  end subroutine base_states__init__copy

  ! ---------------------------------------------------------
  subroutine base_states__build__(this, build)
    type(base_states_t), intent(inout) :: this

    interface
      subroutine build(this, that, name)
        import :: base_states_t
        type(base_states_t), intent(inout) :: this
        type(base_states_t), intent(in)    :: that
        character(len=*),    intent(in)    :: name
      end subroutine build
    end interface

    type(base_states_iterator_t)        :: iter
    character(len=BASE_STATES_NAME_LEN) :: name
    type(base_states_t),        pointer :: subs
    integer                             :: ierr

    PUSH_SUB(base_states__build__)

    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, name, subs, ierr)
      if(ierr/=BASE_STATES_OK)exit
      call build(this, subs, name)
    end do
    call base_states_end(iter)
    nullify(subs)

    POP_SUB(base_states__build__)
  end subroutine base_states__build__

  ! ---------------------------------------------------------
  subroutine base_states_init_type(this, config)
    type(base_states_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_states_init_type)

    call base_states__init__(this, config)

    POP_SUB(base_states_init_type)
  end subroutine base_states_init_type
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_init_copy(this, that)
    type(base_states_t), intent(out) :: this
    type(base_states_t), intent(in)  :: that

    PUSH_SUB(base_states_init_copy)

    call base_states__init__(this, that)
    call base_states__build__(this, init)
    
    POP_SUB(base_states_init_copy)

  contains

    recursive subroutine init(this, isub, name)
      type(base_states_t), intent(inout) :: this
      type(base_states_t), intent(in)    :: isub
      character(len=*),    intent(in)    :: name

      type(base_states_t), pointer :: osub
      
      POP_SUB(base_states_init_copy.init)

      nullify(osub)
      if(base_states_list_index(that%list, isub)>0)then
        call base_states_new(this, osub)
        call base_states_init(osub, isub)
        call base_states_set(this, name, osub)
        nullify(osub)
      else
        call base_states_set(this, name, isub)
      end if

      PUSH_SUB(base_states_init_copy.init)
    end subroutine init

  end subroutine base_states_init_copy

  ! ---------------------------------------------------------
  subroutine base_states__start__(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_states__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call base_density__start__(this%density, sim)

    POP_SUB(base_states__start__)
  end subroutine base_states__start__
    
  ! ---------------------------------------------------------
  subroutine base_states__acc__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sim))
    call base_density__acc__(this%density, that%density)

    POP_SUB(base_states__acc__)
  end subroutine base_states__acc__

  ! ---------------------------------------------------------
  subroutine base_states__update__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__update__(this%density)

    POP_SUB(base_states__update__)
  end subroutine base_states__update__

  ! ---------------------------------------------------------
  subroutine base_states__reset__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_density__reset__(this%density)

    POP_SUB(base_states__reset__)
  end subroutine base_states__reset__
    
  ! ---------------------------------------------------------
  subroutine base_states__stop__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_density__stop__(this%density)

    POP_SUB(base_states__stop__)
  end subroutine base_states__stop__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states__apply__(this, operation)
    type(base_states_t), intent(inout) :: this

    interface
      subroutine operation(this)
        import :: base_states_t
        type(base_states_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(base_states_list_iterator_t) :: iter
    type(base_states_t),      pointer :: subs

    PUSH_SUB(base_states__apply__)

    call base_states_list_init(iter, this%list)
    do
      nullify(subs)
      call base_states_list_next(iter, subs)
      if(.not.associated(subs))exit
      call base_states__apply__(subs, operation)
    end do
    call base_states_list_end(iter)
    nullify(subs)
    call operation(this)

    POP_SUB(base_states__apply__)
  end subroutine base_states__apply__

  ! ---------------------------------------------------------
  subroutine base_states_start(this, sim)
    type(base_states_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim

    PUSH_SUB(base_states_start)

    call base_states__apply__(this, start)
    
    POP_SUB(base_states_start)
    
  contains

    subroutine start(this)
      type(base_states_t), intent(inout) :: this

      PUSH_SUB(base_states_start.start)
      
      call base_states__start__(this, sim)

      POP_SUB(base_states_start.start)
    end subroutine start

  end subroutine base_states_start

  ! ---------------------------------------------------------
  recursive subroutine base_states_acc(this)
    type(base_states_t), intent(inout) :: this

    character(len=BASE_STATES_NAME_LEN) :: name

    type(base_states_iterator_t) :: iter
    type(base_states_t), pointer :: subs
    logical                      :: accu

    PUSH_SUB(base_states_acc)

    ASSERT(associated(this%config))
    call base_states_get(this, reduce=accu)
    if(accu) call base_states__reset__(this)
    call base_states_init(iter, this)
    do
      nullify(subs)
      call base_states_next(iter, name, subs)
      if(.not.associated(subs))exit
      call base_states_acc(subs)
      if(accu) call base_states__acc__(this, subs)
    end do
    call base_states_end(iter)
    nullify(subs)
    if(accu) call base_states__update__(this)
    
    POP_SUB(base_states_acc)
  end subroutine base_states_acc

  ! ---------------------------------------------------------
  subroutine base_states_update(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_update)

    call base_states__apply__(this, base_states__update__)
    
    POP_SUB(base_states_update)
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states_reset(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_reset)

    call base_states__apply__(this, base_states__reset__)
    
    POP_SUB(base_states_reset)
  end subroutine base_states_reset

  ! ---------------------------------------------------------
  subroutine base_states_stop(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states_stop)

    call base_states__apply__(this, base_states__stop__)
    
    POP_SUB(base_states_stop)
  end subroutine base_states_stop

  ! ---------------------------------------------------------
  subroutine base_states__set__(this, name, that)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__set__)

    call base_density_set(this%density, name, that%density)

    POP_SUB(base_states__set__)
  end subroutine base_states__set__
    
  ! ---------------------------------------------------------
  subroutine base_states_set_sub(this, name, that)
    type(base_states_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states_set_sub)

    ASSERT(associated(this%config))
    call base_states_dict_set(this%dict, trim(adjustl(name)), that)
    call base_states__set__(this, trim(adjustl(name)), that)

    POP_SUB(base_states_set_sub)
  end subroutine base_states_set_sub
    
  ! ---------------------------------------------------------
  subroutine base_states_get_sub(this, name, that)
    type(base_states_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_states_t), pointer     :: that

    PUSH_SUB(base_states_get_sub)

    ASSERT(associated(this%config))
    nullify(that)
    call base_states_dict_get(this%dict, trim(adjustl(name)), that)

    POP_SUB(base_states_get_sub)
  end subroutine base_states_get_sub

  ! ---------------------------------------------------------
  subroutine base_states_get_sub_density(this, name, that)
    type(base_states_t),   intent(in) :: this
    character(len=*),      intent(in) :: name
    type(base_density_t), pointer     :: that

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_get_sub_density)

    nullify(that, subs)
    call base_states_get(this, name, subs)
    if(associated(subs)) call base_states_get(subs, that)

    POP_SUB(base_states_get_sub_density)
  end subroutine base_states_get_sub_density
    
  ! ---------------------------------------------------------
  subroutine base_states_get_sub_density_1d(this, name, that, spin, total)
    type(base_states_t),          intent(in) :: this
    character(len=*),             intent(in) :: name
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_get_sub_density_1d)

    nullify(that, dnst)
    call base_states_get(this, name, dnst)
    if(associated(dnst)) call base_density_get(dnst, that, spin, total)

    POP_SUB(base_states_get_sub_density_1d)
  end subroutine base_states_get_sub_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_get_sub_density_2d(this, name, that, total)
    type(base_states_t),            intent(in) :: this
    character(len=*),               intent(in) :: name
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    type(base_density_t), pointer :: dnst

    PUSH_SUB(base_states_get_sub_density_2d)

    nullify(that, dnst)
    call base_states_get(this, name, dnst)
    if(associated(dnst)) call base_density_get(dnst, that, total)

    POP_SUB(base_states_get_sub_density_2d)
  end subroutine base_states_get_sub_density_2d

  ! ---------------------------------------------------------
  subroutine base_states_get_info(this, nspin, reduce)
    type(base_states_t), intent(in)  :: this
    integer,   optional, intent(out) :: nspin
    logical,   optional, intent(out) :: reduce

    PUSH_SUB(base_states_get_info)

    call base_density_get(this%density, nspin=nspin, reduce=reduce)

    POP_SUB(base_states_get_info)
  end subroutine base_states_get_info
    
  ! ---------------------------------------------------------
  subroutine base_states_get_charge(this, charge, spin, total)
    type(base_states_t), intent(in)  :: this
    real(kind=wp),       intent(out) :: charge
    integer,   optional, intent(out) :: spin
    logical,   optional, intent(out) :: total

    PUSH_SUB(base_states_get_charge)

    call base_density_get(this%density, charge, spin=spin, total=total)

    POP_SUB(base_states_get_charge)
  end subroutine base_states_get_charge
    
  ! ---------------------------------------------------------
  subroutine base_states_get_config(this, that)
    type(base_states_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_states_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_states_get_config)
  end subroutine base_states_get_config
    
  ! ---------------------------------------------------------
  subroutine base_states_get_simulation(this, that)
    type(base_states_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_states_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_states_get_simulation)
  end subroutine base_states_get_simulation
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density(this, that)
    type(base_states_t),   target, intent(in) :: this
    type(base_density_t), pointer             :: that

    PUSH_SUB(base_states_get_density)

    that => this%density

    POP_SUB(base_states_get_density)
  end subroutine base_states_get_density
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density_1d(this, that, spin, total)
    type(base_states_t),          intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    integer,            optional, intent(in) :: spin
    logical,            optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_1d)

    nullify(that)
    call base_density_get(this%density, that, spin, total)

    POP_SUB(base_states_get_density_1d)
  end subroutine base_states_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_states_get_density_2d(this, that, total)
    type(base_states_t),            intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    logical,              optional, intent(in) :: total

    PUSH_SUB(base_states_get_density_2d)

    nullify(that)
    call base_density_get(this%density, that, total)

    POP_SUB(base_states_get_density_2d)
  end subroutine base_states_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_states__copy__(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states__copy__)

    call base_states__end__(this)
    if(associated(that%config))then
      call base_states__iinit__(this, that%config)
      this%sim => that%sim
      call base_density__copy__(this%density, that%density)
    end if

    POP_SUB(base_states__copy__)
  end subroutine base_states__copy__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_copy_type(this, that)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that

    PUSH_SUB(base_states_copy_type)

    call base_states_end(this)
    call base_states__copy__(this, that)
    call base_states__build__(this, copy)

    POP_SUB(base_states_copy_type)
    
  contains

    recursive subroutine copy(this, isub, name)
      type(base_states_t), intent(inout) :: this
      type(base_states_t), intent(in)    :: isub
      character(len=*),    intent(in)    :: name

      type(base_states_t), pointer :: osub
      
      POP_SUB(base_states_copy_type.copy)

      nullify(osub)
      if(base_states_list_index(that%list, isub)>0)then
        call base_states_new(this, osub)
        call base_states_copy(osub, isub)
        call base_states_set(this, name, osub)
        nullify(osub)
      else
        call base_states_set(this, name, isub)
      end if

      PUSH_SUB(base_states_copy_type.copy)
    end subroutine copy

  end subroutine base_states_copy_type

  ! ---------------------------------------------------------
  subroutine base_states__end__(this)
    type(base_states_t), intent(inout) :: this

    PUSH_SUB(base_states__end__)

    nullify(this%config, this%sim, this%prnt)
    call base_density__end__(this%density)
    call base_states_dict_end(this%dict)
    call base_states_list_end(this%list)

    POP_SUB(base_states__end__)
  end subroutine base_states__end__
    
  ! ---------------------------------------------------------
  recursive subroutine base_states_end_type(this)
    type(base_states_t), intent(inout) :: this

    type(base_states_t), pointer :: subs

    PUSH_SUB(base_states_end_type)

    do
      nullify(subs)
      call base_states_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_states_end(subs)
      call base_states__del__(subs)
    end do
    nullify(subs)
    call base_states__end__(this)

    POP_SUB(base_states_end_type)
  end subroutine base_states_end_type

#define TEMPLATE_PREFIX base_states
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_states_oct_m

!! Local Variables:
!! mode: f90
!! End:
