#include "global.h"
 
#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

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

#define BASE_TEMPLATE_NAME base_system
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::        &
    base_system_t

  public ::                &
    base_system__init__,   &
    base_system__start__,  &
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
    base_system_update, &
    base_system_reset,  &
    base_system_stop,   &
    base_system_set,    &
    base_system_get,    &
    base_system_copy,   &
    base_system_end

#define BASE_TEMPLATE_NAME base_system
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_system_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(refcount_t),    pointer :: rcnt   =>null()
    type(space_t)                :: space
    type(base_geometry_t)        :: geom
    type(base_states_t)          :: st
    type(base_system_dict_t)     :: dict
  end type base_system_t

  interface base_system__init__
    module procedure base_system__init__type
    module procedure base_system__init__copy
  end interface base_system__init__

  interface base_system__sets__
    module procedure base_system__sets__info
    module procedure base_system__sets__type
  end interface base_system__sets__

  interface base_system_new
    module procedure base_system_new_type
    module procedure base_system_new_pass
  end interface base_system_new

  interface base_system_init
    module procedure base_system_init_type
  end interface base_system_init

  interface base_system_set
    module procedure base_system_set_simulation
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
  end interface base_system_get

contains

#define BASE_TEMPLATE_NAME base_system
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  function base_system_new_type(config) result(this)
    type(json_object_t), intent(in)  :: config

    type(base_system_t), pointer :: this

    PUSH_SUB(base_system_new_type)

    this => base_system_new(config, base_system_init_type)

    POP_SUB(base_system_new_type)
  end function base_system_new_type

  ! ---------------------------------------------------------
  function base_system_new_pass(config, init) result(this)
    type(json_object_t), intent(in) :: config

    type(base_system_t), pointer :: this
    
    interface
      subroutine init(this, config)
        use json_oct_m
        import :: base_system_t
        type(base_system_t), intent(out) :: this
        type(json_object_t), intent(in)  :: config
      end subroutine init
    end interface
    
    PUSH_SUB(base_system_new_pass)

    nullify(this)
    SAFE_ALLOCATE(this)
    call init(this, config)
    ASSERT(associated(this%rcnt))
    call refcount_set(this%rcnt, dynamic=.true.)

    POP_SUB(base_system_new_pass)
  end function base_system_new_pass

  ! ---------------------------------------------------------
  subroutine base_system__iinit__(this, config)
    type(base_system_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_system__iinit__)

    nullify(cnfg)
    this%config => config
    this%rcnt => refcount_new()
    call json_get(this%config, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call space_init(this%space, cnfg)
    nullify(cnfg)
    call base_system_dict_init(this%dict)

    POP_SUB(base_system__iinit__)
  end subroutine base_system__iinit__

  ! ---------------------------------------------------------
  subroutine base_system__init__type(this, config)
    type(base_system_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

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
  subroutine base_system__init__copy(this, that, start)
    type(base_system_t), intent(out) :: this
    type(base_system_t), intent(in)  :: that
    logical,   optional, intent(in)  :: start

    logical :: istr

    PUSH_SUB(base_system__init__copy)

    ASSERT(associated(that%config))
    call base_system__iinit__(this, that%config)
    call base_geometry__init__(this%geom, that%geom)
    call base_states__init__(this%st, that%st)
    istr = .true.
    if(present(start)) istr = start
    if(istr)then
      if(present(start))then
        ASSERT(associated(that%sim))
      end if
      if(associated(that%sim)) call base_system__start__(this, that%sim)
    end if

    POP_SUB(base_system__init__copy)
  end subroutine base_system__init__copy

  ! ---------------------------------------------------------
  subroutine base_system_init_type(this, config)
    type(base_system_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_system_init_type)

    call base_system__init__(this, config)

    POP_SUB(base_system_init_type)
  end subroutine base_system_init_type

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
  subroutine base_system__sets__info(this, name, lock, active)
    type(base_system_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    logical,   optional, intent(in)    :: lock
    logical,   optional, intent(in)    :: active

    PUSH_SUB(base_system__sets__info)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    call base_geometry_sets(this%geom, trim(adjustl(name)), lock=lock, active=active)
    call base_states_sets(this%st, trim(adjustl(name)), lock=lock, active=active)

    POP_SUB(base_system__sets__info)
  end subroutine base_system__sets__info

  ! ---------------------------------------------------------
  subroutine base_system__sets__type(this, name, that, config, lock, active)
    type(base_system_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_system_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    logical,   optional, intent(in)    :: lock
    logical,   optional, intent(in)    :: active

    PUSH_SUB(base_system__sets__type)

    ASSERT(associated(this%config))
    ASSERT(len_trim(adjustl(name))>0)
    ASSERT(associated(that%config))
    ASSERT(this%space==that%space)
    call base_geometry_sets(this%geom, trim(adjustl(name)), that%geom, config, lock=lock, active=active)
    call base_states_sets(this%st, trim(adjustl(name)), that%st, config, lock=lock, active=active)

    POP_SUB(base_system__sets__type)
  end subroutine base_system__sets__type

  ! ---------------------------------------------------------
  subroutine base_system__dels__(this, name, that)
    type(base_system_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_system_t), intent(in)    :: that

    PUSH_SUB(base_system__dels__)

    ASSERT(associated(this%config))
    ASSERT(len_trim(name)>0)
    ASSERT(associated(that%config))
    ASSERT(this%space==that%space)
    call base_geometry_dels(this%geom, trim(adjustl(name)), that%geom)
    call base_states_dels(this%st, trim(adjustl(name)), that%st)

    POP_SUB(base_system__dels__)
  end subroutine base_system__dels__

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
  subroutine base_system_get_info(this, nspin)
    type(base_system_t), intent(in)  :: this
    integer,   optional, intent(out) :: nspin

    PUSH_SUB(base_system_get_info)

    call base_states_get(this%st, nspin=nspin)

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

    type(refcount_t), pointer :: rcnt

    PUSH_SUB(base_system__copy__)

    rcnt => this%rcnt
    nullify(this%rcnt)
    call base_system__end__(this)
    if(associated(that%config))then
      call base_system__iinit__(this, that%config)
      call refcount_del(this%rcnt)
      this%sim => that%sim
      call base_geometry__copy__(this%geom, that%geom)
      call base_states__copy__(this%st, that%st)
    end if
    this%rcnt => rcnt
    nullify(rcnt)

    POP_SUB(base_system__copy__)
  end subroutine base_system__copy__

  ! ---------------------------------------------------------
  subroutine base_system__end__(this)
    type(base_system_t), intent(inout) :: this

    PUSH_SUB(base_system__end__)

    nullify(this%config, this%sim)
    if(associated(this%rcnt)) call refcount_del(this%rcnt)
    call space_end(this%space)
    call base_geometry__end__(this%geom)
    call base_states__end__(this%st)
    ASSERT(base_system_dict_len(this%dict)==0)
    call base_system_dict_end(this%dict)

    POP_SUB(base_system__end__)
  end subroutine base_system__end__

end module base_system_oct_m

!! Local Variables:
!! mode: f90
!! End:
