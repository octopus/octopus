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

module base_system_oct_m

  use base_density_oct_m
  use base_geometry_oct_m
  use base_states_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_system
#define DICT_INCLUDE_PREFIX
#include "tdict_inc.F90"
#undef DICT_INCLUDE_PREFIX
#undef DICT_TEMPLATE_NAME

#define TEMPLATE_PREFIX base_system
#define INCLUDE_PREFIX
#include "iterator_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_PREFIX

  implicit none

  private

  public ::                  &
    BASE_SYSTEM_OK,          &
    BASE_SYSTEM_KEY_ERROR,   &
    BASE_SYSTEM_EMPTY_ERROR

  public ::        &
    base_system_t

  public ::                &
    base_system__init__,   &
    base_system__build__,  &
    base_system__start__,  &
    base_system__acc__,    &
    base_system__update__, &
    base_system__reset__,  &
    base_system__stop__,   &
    base_system__copy__,   &
    base_system__end__

  public ::             &
    base_system_new,    &
    base_system_del,    &
    base_system_init,   &
    base_system_start,  &
    base_system_acc,    &
    base_system_update, &
    base_system_reset,  &
    base_system_stop,   &
    base_system_set,    &
    base_system_get,    &
    base_system_copy,   &
    base_system_end

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_system
#define DICT_INCLUDE_HEADER
#include "tdict_inc.F90"
#undef DICT_INCLUDE_HEADER
#undef DICT_TEMPLATE_NAME

  integer, parameter :: BASE_SYSTEM_OK          = BASE_SYSTEM_HASH_OK
  integer, parameter :: BASE_SYSTEM_KEY_ERROR   = BASE_SYSTEM_HASH_KEY_ERROR
  integer, parameter :: BASE_SYSTEM_EMPTY_ERROR = BASE_SYSTEM_HASH_EMPTY_ERROR

  type :: base_system_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_system_t), pointer :: prnt   =>null()
    type(space_t)                :: space
    type(base_geometry_t)        :: geom
    type(base_states_t)          :: st
    type(base_system_dict_t)     :: dict
    type(base_system_list_t)     :: list
  end type base_system_t

  interface base_system__init__
    module procedure base_system__init__type
    module procedure base_system__init__copy
  end interface base_system__init__

  interface base_system_init
    module procedure base_system_init_type
    module procedure base_system_init_copy
  end interface base_system_init

  interface base_system_set
    module procedure base_system_set_simulation
    module procedure base_system_set_sub
  end interface base_system_set

  interface base_system_get
    module procedure base_system_get_info
    module procedure base_system_get_config
    module procedure base_system_get_simulation
    module procedure base_system_get_space
    module procedure base_system_get_geom
    module procedure base_system_get_geometry
    module procedure base_system_get_charge
    module procedure base_system_get_states
    module procedure base_system_get_density
    module procedure base_system_get_sub
  end interface base_system_get

  interface base_system_copy
    module procedure base_system_copy_type
  end interface base_system_copy

  interface base_system_end
    module procedure base_system_end_type
  end interface base_system_end

#define TEMPLATE_PREFIX base_system
#define INCLUDE_HEADER
#include "iterator_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_PREFIX

contains

#define LIST_TEMPLATE_NAME base_system
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY
#undef LIST_TEMPLATE_NAME

#define DICT_TEMPLATE_NAME base_system
#define DICT_INCLUDE_BODY
#include "tdict_inc.F90"
#undef DICT_INCLUDE_BODY
#undef DICT_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_system__new__(this)
    type(base_system_t), pointer :: this

    PUSH_SUB(base_system__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_system__new__)
  end subroutine base_system__new__

  ! ---------------------------------------------------------
  subroutine base_system__del__(this)
    type(base_system_t), pointer :: this

    PUSH_SUB(base_system__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_system__del__)
  end subroutine base_system__del__

  ! ---------------------------------------------------------
  subroutine base_system_new(this, that)
    type(base_system_t),  target, intent(inout) :: this
    type(base_system_t), pointer                :: that

    PUSH_SUB(base_system_new)

    nullify(that)
    call base_system__new__(that)
    that%prnt => this
    call base_system_list_push(this%list, that)

    POP_SUB(base_system_new)
  end subroutine base_system_new

  ! ---------------------------------------------------------
  subroutine base_system_del(this)
    type(base_system_t), pointer :: this

    PUSH_SUB(base_system_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_system_list_del(this%prnt%list, this)
        call base_system_end(this)
        call base_system__del__(this)
      end if
    end if

    POP_SUB(base_system_del)
  end subroutine base_system_del

  ! ---------------------------------------------------------
  subroutine base_system__iinit__(this, config)
    type(base_system_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system__iinit__)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call space_init(this%space, cnfg)
    nullify(cnfg)
    call base_system_dict_init(this%dict)
    call base_system_list_init(this%list)

    POP_SUB(base_system__iinit__)
  end subroutine base_system__iinit__

  ! ---------------------------------------------------------
  subroutine base_system__init__type(this, config)
    type(base_system_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system__init__type)

    nullify(cnfg)
    call base_system__iinit__(this, config)
    call json_get(this%config, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_geometry__init__(this%geom, this%space, cnfg)
    nullify(cnfg)
    call json_get(this%config, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_states__init__(this%st, cnfg)
    nullify(cnfg)

    POP_SUB(base_system__init__type)
  end subroutine base_system__init__type

  ! ---------------------------------------------------------
  subroutine base_system__init__copy(this, that)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that

    PUSH_SUB(base_system__init__copy)

    ASSERT(associated(that%config))
    call base_system__iinit__(this, that%config)
    call base_geometry__init__(this%geom, that%geom)
    call base_states__init__(this%st, that%st)
    if(associated(that%sim)) call base_system__start__(this, that%sim)

    POP_SUB(base_system__init__copy)
  end subroutine base_system__init__copy

  ! ---------------------------------------------------------
  subroutine base_system__build__(this, build)
    type(base_system_t), intent(inout) :: this

    interface
      subroutine build(this, that, name)
        import :: base_system_t
        type(base_system_t), intent(inout) :: this
        type(base_system_t), intent(in)    :: that
        character(len=*),    intent(in)    :: name
      end subroutine build
    end interface

    type(base_system_iterator_t)        :: iter
    character(len=BASE_SYSTEM_NAME_LEN) :: name
    type(base_system_t),        pointer :: subs
    integer                             :: ierr

    PUSH_SUB(base_system__build__)

    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, name, subs, ierr)
      if(ierr/=BASE_SYSTEM_OK)exit
      call build(this, subs, name)
    end do
    call base_system_end(iter)
    nullify(subs)

    POP_SUB(base_system__build__)
  end subroutine base_system__build__

  ! ---------------------------------------------------------
  subroutine base_system_init_type(this, config)
    type(base_system_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_system_init_type)

    call base_system__init__(this, config)

    POP_SUB(base_system_init_type)
  end subroutine base_system_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_system_init_copy(this, that)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that

    PUSH_SUB(base_system_init_copy)

    call base_system__init__(this, that)
    call base_system__build__(this, init)

    POP_SUB(base_system_init_copy)

  contains

    recursive subroutine init(this, isub, name)
      type(base_system_t), intent(inout) :: this
      type(base_system_t), intent(in)    :: isub
      character(len=*),    intent(in)    :: name

      type(base_system_t), pointer :: osub
      
      POP_SUB(base_system_init_copy.init)

      nullify(osub)
      if(base_system_list_index(that%list, isub)>0)then
        call base_system_new(this, osub)
        call base_system_init(osub, isub)
        call base_system_set(this, name, osub)
        nullify(osub)
      else
        call base_system_set(this, name, isub)
      end if

      PUSH_SUB(base_system_init_copy.init)
    end subroutine init

  end subroutine base_system_init_copy

  ! ---------------------------------------------------------
  subroutine base_system__start__(this, sim)
    type(base_system_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_system__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call base_states__start__(this%st, sim)

    POP_SUB(base_system__start__)
  end subroutine base_system__start__

  ! ---------------------------------------------------------
  subroutine base_system__acc__(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__acc__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%config))
    ASSERT(associated(that%sim))
    call base_states__acc__(this%st, that%st)

    POP_SUB(base_system__acc__)
  end subroutine base_system__acc__

  ! ---------------------------------------------------------
  subroutine base_system__update__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__update__(this%st)

    POP_SUB(base_system__update__)
  end subroutine base_system__update__

  ! ---------------------------------------------------------
  subroutine base_system__reset__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_states__reset__(this%st)

    POP_SUB(base_system__reset__)
  end subroutine base_system__reset__

  ! ---------------------------------------------------------
  subroutine base_system__stop__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_states__stop__(this%st)

    POP_SUB(base_system__stop__)
  end subroutine base_system__stop__

  ! ---------------------------------------------------------
  recursive subroutine base_system__apply__(this, operation)
    type(base_system_t), intent(inout) :: this

    interface
      subroutine operation(this)
        import :: base_system_t
        type(base_system_t), intent(inout) :: this
      end subroutine operation
    end interface

    type(base_system_list_iterator_t) :: iter
    type(base_system_t),      pointer :: subs

    PUSH_SUB(base_system__apply__)

    call base_system_list_init(iter, this%list)
    do
      nullify(subs)
      call base_system_list_next(iter, subs)
      if(.not.associated(subs))exit
      call base_system__apply__(subs, operation)
    end do
    call base_system_list_end(iter)
    nullify(subs)
    call operation(this)

    POP_SUB(base_system__apply__)
  end subroutine base_system__apply__

  ! ---------------------------------------------------------
  subroutine base_system_start(this, sim)
    type(base_system_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim

    PUSH_SUB(base_system_start)

    call base_system__apply__(this, start)

    POP_SUB(base_system_start)
    
  contains

    subroutine start(this)
      type(base_system_t), intent(inout) :: this

      PUSH_SUB(base_system_start.start)
      
      call base_system__start__(this, sim)

      POP_SUB(base_system_start.start)
    end subroutine start

  end subroutine base_system_start

  ! ---------------------------------------------------------
  recursive subroutine base_system_acc(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_iterator_t) :: iter
    type(base_system_t), pointer :: subs
    logical                      :: accu

    PUSH_SUB(base_system_acc)

    ASSERT(associated(this%config))
    call base_system_get(this, reduce=accu)
    if(accu) call base_system__reset__(this)
    call base_system_init(iter, this)
    do
      nullify(subs)
      call base_system_next(iter, subs)
      if(.not.associated(subs))exit
      call base_system_acc(subs)
      if(accu) call base_system__acc__(this, subs)
    end do
    call base_system_end(iter)
    nullify(subs)
    if(accu) call base_system__update__(this)
    
    POP_SUB(base_system_acc)
  end subroutine base_system_acc

  ! ---------------------------------------------------------
  subroutine base_system_update(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system_update)

    call base_system__apply__(this, base_system__update__)

    POP_SUB(base_system_update)
  end subroutine base_system_update

  ! ---------------------------------------------------------
  subroutine base_system_reset(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system_reset)

    call base_system__apply__(this, base_system__reset__)
    
    POP_SUB(base_system_reset)
  end subroutine base_system_reset

  ! ---------------------------------------------------------
  subroutine base_system_stop(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system_stop)

    call base_system__apply__(this, base_system__stop__)

    POP_SUB(base_system_stop)
  end subroutine base_system_stop

  ! ---------------------------------------------------------
  subroutine base_system__set__(this, name, that)
    type(base_system_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__set__)

    ASSERT(this%space==that%space)
    call base_geometry_set(this%geom, name, that%geom)
    call base_states_set(this%st, name, that%st)

    POP_SUB(base_system__set__)
  end subroutine base_system__set__

  ! ---------------------------------------------------------
  subroutine base_system_set_sub(this, name, that)
    type(base_system_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system_set_sub)

    ASSERT(associated(this%config))
    call base_system_dict_set(this%dict, trim(adjustl(name)), that)
    call base_system__set__(this, trim(adjustl(name)), that)

    POP_SUB(base_system_set_sub)
  end subroutine base_system_set_sub

  ! ---------------------------------------------------------
  subroutine base_system_get_sub(this, name, that)
    type(base_system_t),  intent(in) :: this
    character(len=*),     intent(in) :: name
    type(base_system_t), pointer     :: that

    PUSH_SUB(base_system_get_sub)

    nullify(that)
    ASSERT(associated(this%config))
    call base_system_dict_get(this%dict, trim(adjustl(name)), that)

    POP_SUB(base_system_get_sub)
  end subroutine base_system_get_sub

  ! ---------------------------------------------------------
  subroutine base_system_set_simulation(this, that)
    type(base_system_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that

    PUSH_SUB(base_system_set_simulation)

    ASSERT(.not.associated(this%sim))
    this%sim => that

    POP_SUB(base_system_set_simulation)
  end subroutine base_system_set_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_info(this, nspin, reduce)
    type(base_system_t), intent(in)  :: this
    integer,   optional, intent(out) :: nspin
    logical,   optional, intent(out) :: reduce

    PUSH_SUB(base_system_get_info)

    call base_states_get(this%st, nspin=nspin, reduce=reduce)

    POP_SUB(base_system_get_info)
  end subroutine base_system_get_info
 
  ! ---------------------------------------------------------
  subroutine base_system_get_charge(this, charge, spin)
    type(base_system_t), intent(in)  :: this
    real(kind=wp),       intent(out) :: charge
    integer,   optional, intent(out) :: spin

    PUSH_SUB(base_system_get_charge)

    call base_states_get(this%st, charge, spin=spin)

    POP_SUB(base_system_get_charge)
  end subroutine base_system_get_charge
 
  ! ---------------------------------------------------------
  subroutine base_system_get_config(this, that)
    type(base_system_t),  target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_system_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_system_get_config)
  end subroutine base_system_get_config

  ! ---------------------------------------------------------
  subroutine base_system_get_simulation(this, that)
    type(base_system_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_system_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_system_get_simulation)
  end subroutine base_system_get_simulation

  ! ---------------------------------------------------------
  subroutine base_system_get_space(this, that)
    type(base_system_t), target, intent(in) :: this
    type(space_t),      pointer             :: that

    PUSH_SUB(base_system_get_space)

    that => this%space

    POP_SUB(base_system_get_space)
  end subroutine base_system_get_space

  ! ---------------------------------------------------------
  subroutine base_system_get_geom(this, that)
    type(base_system_t),     target, intent(in) :: this
    type(base_geometry_t),  pointer             :: that

    PUSH_SUB(base_system_get_geom)

    that => this%geom

    POP_SUB(base_system_get_geom)
  end subroutine base_system_get_geom

  ! ---------------------------------------------------------
  subroutine base_system_get_geometry(this, that)
    type(base_system_t), intent(in) :: this
    type(geometry_t),   pointer     :: that

    PUSH_SUB(base_system_get_geometry)

    call base_geometry_get(this%geom, that)

    POP_SUB(base_system_get_geometry)
  end subroutine base_system_get_geometry

  ! ---------------------------------------------------------
  subroutine base_system_get_states(this, that)
    type(base_system_t),  target, intent(in) :: this
    type(base_states_t), pointer             :: that

    PUSH_SUB(base_system_get_states)

    that => this%st

    POP_SUB(base_system_get_states)
  end subroutine base_system_get_states

  ! ---------------------------------------------------------
  subroutine base_system_get_density(this, that)
    type(base_system_t),   intent(in) :: this
    type(base_density_t), pointer     :: that

    PUSH_SUB(base_system_get_density)

    call base_states_get(this%st, that)

    POP_SUB(base_system_get_density)
  end subroutine base_system_get_density

  ! ---------------------------------------------------------
  subroutine base_system__copy__(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__copy__)

    call base_system__end__(this)
    if(associated(that%config))then
      call base_system__iinit__(this, that%config)
      this%sim => that%sim
      call base_geometry__copy__(this%geom, that%geom)
      call base_states__copy__(this%st, that%st)
    end if

    POP_SUB(base_system__copy__)
  end subroutine base_system__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_system_copy_type(this, that)
    type(base_system_t), intent(inout) :: this
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system_copy_type)

    call base_system_end(this)
    call base_system__copy__(this, that)
    call base_system__build__(this, copy)

    POP_SUB(base_system_copy_type)

  contains

    recursive subroutine copy(this, isub, name)
      type(base_system_t), intent(inout) :: this
      type(base_system_t), intent(in)    :: isub
      character(len=*),    intent(in)    :: name

      type(base_system_t), pointer :: osub
      
      POP_SUB(base_system_copy_type.copy)

      nullify(osub)
      if(base_system_list_index(that%list, isub)>0)then
        call base_system_new(this, osub)
        call base_system_copy(osub, isub)
        call base_system_set(this, name, osub)
        nullify(osub)
      else
        call base_system_set(this, name, isub)
      end if

      PUSH_SUB(base_system_copy_type.copy)
    end subroutine copy

  end subroutine base_system_copy_type

  ! ---------------------------------------------------------
  subroutine base_system__end__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__end__)

    nullify(this%config, this%sim, this%prnt)
    call space_end(this%space)
    call base_geometry__end__(this%geom)
    call base_states__end__(this%st)
    call base_system_dict_end(this%dict)
    call base_system_list_end(this%list)

    POP_SUB(base_system__end__)
  end subroutine base_system__end__

  ! ---------------------------------------------------------
  recursive subroutine base_system_end_type(this)
    type(base_system_t), intent(inout) :: this

    type(base_system_t), pointer :: subs

    PUSH_SUB(base_system_end_type)

    do
      nullify(subs)
      call base_system_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_system_end(subs)
      call base_system__del__(subs)
    end do
    nullify(subs)
    call base_system__end__(this)

    POP_SUB(base_system_end_type)
  end subroutine base_system_end_type

#define TEMPLATE_PREFIX base_system
#define INCLUDE_BODY
#include "iterator_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_PREFIX

end module base_system_oct_m

!! Local Variables:
!! mode: f90
!! End:
