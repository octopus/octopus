#include "global.h"

#undef BASE_TEMPLATE_NAME
#undef BASE_TYPE_NAME
#undef BASE_TYPE_MODULE_NAME
#undef BASE_INCLUDE_PREFIX
#undef BASE_INCLUDE_HEADER
#undef BASE_INCLUDE_BODY

module base_model_oct_m

  use base_density_oct_m
  use base_geometry_oct_m
  use base_hamiltonian_oct_m
  use base_states_oct_m
  use base_system_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

#define BASE_TEMPLATE_NAME base_model
#define BASE_INCLUDE_PREFIX
#include "tbase_inc.F90"
#undef BASE_INCLUDE_PREFIX
#undef BASE_TEMPLATE_NAME

  implicit none

  private

  public ::       &
    base_model_t

  public ::               &
    base_model__init__,   &
    base_model__start__,  &
    base_model__update__, &
    base_model__reset__,  &
    base_model__stop__,   &
    base_model__copy__,   &
    base_model__end__

  public ::            &
    base_model_new,    &
    base_model_del,    &
    base_model_init,   &
    base_model_start,  &
    base_model_update, &
    base_model_reset,  &
    base_model_stop,   &
    base_model_get,    &
    base_model_copy,   &
    base_model_end

#define BASE_TEMPLATE_NAME base_model
#define BASE_INCLUDE_HEADER
#include "tbase_inc.F90"
#undef BASE_INCLUDE_HEADER
#undef BASE_TEMPLATE_NAME

  type :: base_model_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(base_model_t),  pointer :: prnt   =>null()
    type(base_system_t)          :: sys
    type(base_hamiltonian_t)     :: hm
    type(base_model_dict_t)      :: dict
    type(base_model_list_t)      :: list
  end type base_model_t

  interface base_model__init__
    module procedure base_model__init__type
    module procedure base_model__init__copy
  end interface base_model__init__

  interface base_model_init
    module procedure base_model_init_type
    module procedure base_model_init_copy
  end interface base_model_init

  interface base_model_get
    module procedure base_model_get_config
    module procedure base_model_get_simulation
    module procedure base_model_get_space
    module procedure base_model_get_geometry
    module procedure base_model_get_geom
    module procedure base_model_get_density
    module procedure base_model_get_states
    module procedure base_model_get_system
    module procedure base_model_get_hamiltonian
  end interface base_model_get

  interface base_model_copy
    module procedure base_model_copy_type
  end interface base_model_copy

  interface base_model_end
    module procedure base_model_end_type
  end interface base_model_end

contains

#define BASE_TEMPLATE_NAME base_model
#define BASE_INCLUDE_BODY
#include "tbase_inc.F90"
#undef BASE_INCLUDE_BODY
#undef BASE_TEMPLATE_NAME

  ! ---------------------------------------------------------
  subroutine base_model__new__(this)
    type(base_model_t), pointer :: this

    PUSH_SUB(base_model__new__)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(base_model__new__)
  end subroutine base_model__new__

  ! ---------------------------------------------------------
  subroutine base_model__del__(this)
    type(base_model_t), pointer :: this

    PUSH_SUB(base_model__del__)

    if(associated(this))then
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(base_model__del__)
  end subroutine base_model__del__

  ! ---------------------------------------------------------
  subroutine base_model_new(this, that)
    type(base_model_t),  target, intent(inout) :: this
    type(base_model_t), pointer                :: that

    PUSH_SUB(base_model_new)

    nullify(that)
    call base_model__new__(that)
    that%prnt => this
    call base_model_list_push(this%list, that)

    POP_SUB(base_model_new)
  end subroutine base_model_new

  ! ---------------------------------------------------------
  subroutine base_model_del(this)
    type(base_model_t), pointer :: this

    PUSH_SUB(base_model_del)

    if(associated(this))then
      if(associated(this%prnt))then
        call base_model_list_del(this%prnt%list, this)
        call base_model_end(this)
        call base_model__del__(this)
      end if
    end if

    POP_SUB(base_model_del)
  end subroutine base_model_del

  ! ---------------------------------------------------------
  subroutine base_model__init__type(this, config)
    type(base_model_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_model__init__type)

    nullify(cnfg)
    this%config => config
    call json_get(this%config, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_system__init__(this%sys, cnfg)
    nullify(cnfg)
    call json_get(this%config, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call base_hamiltonian__init__(this%hm, this%sys, cnfg)
    nullify(cnfg)
    call base_model_dict_init(this%dict)
    call base_model_list_init(this%list)

    POP_SUB(base_model__init__type)
  end subroutine base_model__init__type

  ! ---------------------------------------------------------
  subroutine base_model__init__copy(this, that, start)
    type(base_model_t), intent(out) :: this
    type(base_model_t), intent(in)  :: that
    logical,  optional, intent(in)  :: start

    logical :: istr

    PUSH_SUB(base_model__init__copy)

    ASSERT(associated(that%config))
    call base_model__init__(this, that%config)
    istr = .true.
    if(present(start)) istr = start
    if(istr)then
      if(present(start))then
        ASSERT(associated(that%sim))
      end if
      if(associated(that%sim)) call base_model__start__(this, that%sim)
    end if

    POP_SUB(base_model__init__copy)
  end subroutine base_model__init__copy

  ! ---------------------------------------------------------
  subroutine base_model_init_type(this, config)
    type(base_model_t),  intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(base_model_init_type)

    call base_model__init__(this, config)

    POP_SUB(base_model_init_type)
  end subroutine base_model_init_type

  ! ---------------------------------------------------------
  recursive subroutine base_model_init_copy(this, that)
    type(base_model_t), intent(out) :: this
    type(base_model_t), intent(in)  :: that

    type(base_model_iterator_t)        :: iter
    character(len=BASE_MODEL_NAME_LEN) :: name
    type(base_model_t),        pointer :: osub, isub
    integer                            :: ierr

    PUSH_SUB(base_model_init_copy)

    call base_model__init__(this, that)
    call base_model_init(iter, that)
    do
      nullify(osub, isub)
      call base_model_next(iter, name, isub, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_new(this, osub)
      call base_model_init(osub, isub)
      call base_model_sets(this, name, osub)
    end do
    call base_model_end(iter)
    nullify(osub, isub)

    POP_SUB(base_model_init_copy)
  end subroutine base_model_init_copy

  ! ---------------------------------------------------------
  subroutine base_model__start__(this, sim)
    type(base_model_t),         intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    PUSH_SUB(base_model__start__)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call base_system__start__(this%sys, this%sim)
    call base_hamiltonian__start__(this%hm, this%sim)

    POP_SUB(base_model__start__)
  end subroutine base_model__start__

  ! ---------------------------------------------------------
  subroutine base_model__update__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__update__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_system__update__(this%sys)
    call base_hamiltonian__update__(this%hm)

    POP_SUB(base_model__update__)
  end subroutine base_model__update__

  ! ---------------------------------------------------------
  subroutine base_model__reset__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__reset__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call base_system__reset__(this%sys)
    call base_hamiltonian__reset__(this%hm)

    POP_SUB(base_model__reset__)
  end subroutine base_model__reset__

  ! ---------------------------------------------------------
  subroutine base_model__stop__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__stop__)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call base_system__stop__(this%sys)
    call base_hamiltonian__stop__(this%hm)

    POP_SUB(base_model__stop__)
  end subroutine base_model__stop__

  ! ---------------------------------------------------------
  subroutine base_model_start(this, sim)
    type(base_model_t), intent(inout) :: this
    type(simulation_t), intent(in)    :: sim

    PUSH_SUB(base_model_start)

    call base_model__apply__(this, start)

    POP_SUB(base_model_start)
    
  contains

    subroutine start(this)
      type(base_model_t), intent(inout) :: this

      PUSH_SUB(base_model_start.start)
      
      call base_model__start__(this, sim)

      POP_SUB(base_model_start.start)
    end subroutine start

  end subroutine base_model_start

  ! ---------------------------------------------------------
  subroutine base_model_update(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model_update)

    call base_model__apply__(this, base_model__update__)

    POP_SUB(base_model_update)
  end subroutine base_model_update

  ! ---------------------------------------------------------
  subroutine base_model_reset(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model_reset)

    call base_model__apply__(this, base_model__reset__)
    
    POP_SUB(base_model_reset)
  end subroutine base_model_reset

  ! ---------------------------------------------------------
  subroutine base_model_stop(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model_stop)

    call base_model__apply__(this, base_model__stop__)

    POP_SUB(base_model_stop)
  end subroutine base_model_stop

  ! ---------------------------------------------------------
  subroutine base_model__sets__(this, name, that)
    type(base_model_t),  intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(base_model_t),  intent(in)    :: that

    PUSH_SUB(base_model__sets__)

    call base_system_sets(this%sys, trim(adjustl(name)), that%sys)
    call base_hamiltonian_sets(this%hm, trim(adjustl(name)), that%hm)

    POP_SUB(base_model__sets__)
  end subroutine base_model__sets__

  ! ---------------------------------------------------------
  subroutine base_model__dels__(this, name, ierr)
    type(base_model_t),  intent(inout) :: this
    character(len=*),    intent(in)    :: name
    integer,             intent(out)   :: ierr

    PUSH_SUB(base_model__dels__)

    call base_hamiltonian_dels(this%hm, trim(adjustl(name)), ierr)
    call base_system_dels(this%sys, trim(adjustl(name)), ierr)

    POP_SUB(base_model__dels__)
  end subroutine base_model__dels__

  ! ---------------------------------------------------------
  subroutine base_model_get_config(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(base_model_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(base_model_get_config)
  end subroutine base_model_get_config

  ! ---------------------------------------------------------
  subroutine base_model_get_simulation(this, that)
    type(base_model_t),  target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(base_model_get_simulation)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(base_model_get_simulation)
  end subroutine base_model_get_simulation

  ! ---------------------------------------------------------
  subroutine base_model_get_space(this, that)
    type(base_model_t), intent(in) :: this
    type(space_t),     pointer     :: that

    PUSH_SUB(base_model_get_space)

    nullify(that)
    if(associated(this%config)) call base_system_get(this%sys, that)

    POP_SUB(base_model_get_space)
  end subroutine base_model_get_space

  ! ---------------------------------------------------------
  subroutine base_model_get_geometry(this, that)
    type(base_model_t), intent(in) :: this
    type(geometry_t),  pointer     :: that

    PUSH_SUB(base_model_get_geometry)

    nullify(that)
    if(associated(this%config)) call base_system_get(this%sys, that)

    POP_SUB(base_model_get_geometry)
  end subroutine base_model_get_geometry

  ! ---------------------------------------------------------
  subroutine base_model_get_geom(this, that)
    type(base_model_t),     intent(in) :: this
    type(base_geometry_t), pointer     :: that

    PUSH_SUB(base_model_get_geom)

    nullify(that)
    if(associated(this%config)) call base_system_get(this%sys, that)

    POP_SUB(base_model_get_geom)
  end subroutine base_model_get_geom

  ! ---------------------------------------------------------
  subroutine base_model_get_density(this, that)
    type(base_model_t) ,   intent(in) :: this
    type(base_density_t), pointer     :: that

    PUSH_SUB(base_model_get_density)

    nullify(that)
    if(associated(this%config)) call base_system_get(this%sys, that)

    POP_SUB(base_model_get_density)
  end subroutine base_model_get_density

  ! ---------------------------------------------------------
  subroutine base_model_get_states(this, that)
    type(base_model_t) ,  intent(in) :: this
    type(base_states_t), pointer     :: that

    PUSH_SUB(base_model_get_states)

    nullify(that)
    if(associated(this%config)) call base_system_get(this%sys, that)

    POP_SUB(base_model_get_states)
  end subroutine base_model_get_states

  ! ---------------------------------------------------------
  subroutine base_model_get_system(this, that)
    type(base_model_t),   target, intent(in) :: this
    type(base_system_t), pointer             :: that

    PUSH_SUB(base_model_get_system)

    nullify(that)
    if(associated(this%config)) that => this%sys

    POP_SUB(base_model_get_system)
  end subroutine base_model_get_system

  ! ---------------------------------------------------------
  subroutine base_model_get_hamiltonian(this, that)
    type(base_model_t),        target, intent(in) :: this
    type(base_hamiltonian_t), pointer             :: that

    PUSH_SUB(base_model_get_hamiltonian)

    nullify(that)
    if(associated(this%config)) that => this%hm

    POP_SUB(base_model_get_hamiltonian)
  end subroutine base_model_get_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_model__copy__(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    PUSH_SUB(base_model__copy__)

    call base_model__end__(this)
    if(associated(that%config))then
      call base_model__init__(this, that)
      if(associated(that%sim))then
        call base_system__copy__(this%sys, that%sys)
        call base_hamiltonian__copy__(this%hm, that%hm)
      end if
    end if

    POP_SUB(base_model__copy__)
  end subroutine base_model__copy__

  ! ---------------------------------------------------------
  recursive subroutine base_model_copy_type(this, that)
    type(base_model_t), intent(inout) :: this
    type(base_model_t), intent(in)    :: that

    type(base_model_iterator_t)        :: iter
    character(len=BASE_MODEL_NAME_LEN) :: name
    type(base_model_t),        pointer :: osub, isub
    integer                            :: ierr

    PUSH_SUB(base_model_copy_type)

    nullify(osub, isub)
    call base_model_end(this)
    call base_model__copy__(this, that)
    call base_model_init(iter, that)
    do
      nullify(osub, isub)
      call base_model_next(iter, name, isub, ierr)
      if(ierr/=BASE_MODEL_OK)exit
      call base_model_new(this, osub)
      call base_model_copy(osub, isub)
      call base_model_sets(this, name, osub)
    end do
    call base_model_end(iter)
    nullify(osub, isub)

    POP_SUB(base_model_copy_type)
  end subroutine base_model_copy_type

  ! ---------------------------------------------------------
  subroutine base_model__end__(this)
    type(base_model_t), intent(inout) :: this

    PUSH_SUB(base_model__end__)

    nullify(this%config, this%sim, this%prnt)
    call base_hamiltonian__end__(this%hm)
    call base_system__end__(this%sys)
    call base_model_dict_end(this%dict)
    call base_model_list_end(this%list)

    POP_SUB(base_model__end__)
  end subroutine base_model__end__

  ! ---------------------------------------------------------
  recursive subroutine base_model_end_type(this)
    type(base_model_t), intent(inout) :: this

    type(base_model_t), pointer :: subs

    PUSH_SUB(base_model_end_type)

    do
      nullify(subs)
      call base_model_list_pop(this%list, subs)
      if(.not.associated(subs))exit
      call base_model_end(subs)
      call base_model__del__(subs)
    end do
    nullify(subs)
    call base_model__end__(this)

    POP_SUB(base_model_end_type)
  end subroutine base_model_end_type

end module base_model_oct_m

!! Local Variables:
!! mode: f90
!! End:
