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

#define HASH_TEMPLATE_NAME base_states
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_states

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_states_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

  use simulation_m, only: &
    simulation_t

  use base_density_m, only:                &
    density_t      => base_density_t,      &
    density_init   => base_density_init,   &
    density_start  => base_density_start,  &
    density_update => base_density_update, &
    density_stop   => base_density_stop,   &
    density_get    => base_density_get,    &
    density_copy   => base_density_copy,   &
    density_end    => base_density_end

  implicit none

  private
  public ::                 &
    base_states_init,       &
    base_states_start,      &
    base_states_update,     &
    base_states_stop,       &
    base_states_set,        &
    base_states_get,        &
    base_states_get_charge, &
    base_states_set_charge, &
    base_states_copy,       &
    base_states_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_states_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    real(kind=wp)                :: charge = 0.0_wp
    type(density_t)              :: density
    type(base_states_hash_t)     :: hash
  end type base_states_t

  type, public :: base_states_iterator_t
    private
    type(base_states_t),      pointer :: self =>null()
    type(base_states_hash_iterator_t) :: iter
  end type base_states_iterator_t

  interface base_states_init
    module procedure base_states_init_begin
    module procedure base_states_init_build
  end interface base_states_init

  interface base_states_set
    module procedure base_states_set_simulation
  end interface base_states_set

  interface base_states_get
    module procedure base_states_get_config
    module procedure base_states_get_simulation
    module procedure base_states_get_density
  end interface base_states_get

  interface base_states_copy
    module procedure base_states_copy_states
    module procedure base_states_iterator_copy
  end interface base_states_copy

  interface base_states_end
    module procedure base_states_end_states
    module procedure base_states_iterator_end
  end interface base_states_end

contains
    
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_states_init_begin(this, config)
    type(base_states_t),         intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(base_states_init_begin)
    this%config=>config
    nullify(cnfg, this%sim)
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK)this%charge=0.0_wp
    call json_get(this%config, "density", cnfg, ierr)
    if(ierr==JSON_OK)call density_init(this%density, cnfg)
    nullify(cnfg)
    call base_states_hash_init(this%hash)
    POP_SUB(base_states_init_begin)
    return
  end subroutine base_states_init_begin
    
  ! ---------------------------------------------------------
  subroutine base_states_init_build(this, that, config)
    type(base_states_t), intent(inout) :: this
    type(base_states_t), intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(base_states_init_build)
    this%charge=this%charge+base_states_get_charge(that)
    call density_init(this%density, that%density, config)
    call base_states_hash_set(this%hash, config, that)
    POP_SUB(base_states_init_build)
    return
  end subroutine base_states_init_build
    
  ! ---------------------------------------------------------
  subroutine base_states_start(this, sim)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_states_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call density_start(this%density, sim)
    POP_SUB(base_states_start)
    return
  end subroutine base_states_start
    
  ! ---------------------------------------------------------
  subroutine base_states_update(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call density_update(this%density)
    POP_SUB(base_states_update)
    return
  end subroutine base_states_update

  ! ---------------------------------------------------------
  subroutine base_states_stop(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call density_stop(this%density)
    POP_SUB(base_states_stop)
    return
  end subroutine base_states_stop
    
  ! ---------------------------------------------------------
  elemental subroutine base_states_set_charge(this, that)
    type(base_states_t), intent(inout) :: this
    real(kind=wp),       intent(in)    :: that
    !
    this%charge=that
    return
  end subroutine base_states_set_charge
    
  ! ---------------------------------------------------------
  subroutine base_states_set_simulation(this, that)
    type(base_states_t),        intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_states_set_simulation)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>that
    POP_SUB(base_states_set_simulation)
    return
  end subroutine base_states_set_simulation
    
  ! ---------------------------------------------------------
  elemental function base_states_get_charge(this) result(that)
    type(base_states_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%charge
    return
  end function base_states_get_charge
    
  ! ---------------------------------------------------------
  subroutine base_states_get_config(this, that)
    type(base_states_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_states_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_states_get_config)
    return
  end subroutine base_states_get_config
    
  ! ---------------------------------------------------------
  subroutine base_states_get_simulation(this, that)
    type(base_states_t), target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(base_states_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_states_get_simulation)
    return
  end subroutine base_states_get_simulation
    
  ! ---------------------------------------------------------
  subroutine base_states_get_density(this, that)
    type(base_states_t), target, intent(in) :: this
    type(density_t),    pointer             :: that
    !
    PUSH_SUB(base_states_get_density)
    that=>this%density
    POP_SUB(base_states_get_density)
    return
  end subroutine base_states_get_density
    
  ! ---------------------------------------------------------
  subroutine base_states_copy_states(this, that)
    type(base_states_t),         intent(out) :: this
    type(base_states_t), target, intent(in)  :: that
    !
    PUSH_SUB(base_states_copy_states)
    this%config=>that%config
    this%sim=>that%sim
    this%charge=that%charge
    call density_copy(this%density, that%density)
    call base_states_hash_copy(this%hash, that%hash)
    POP_SUB(base_states_copy_states)
    return
  end subroutine base_states_copy_states
    
  ! ---------------------------------------------------------
  subroutine base_states_end_states(this)
    type(base_states_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_end_states)
    nullify(this%config, this%sim)
    this%charge=0.0_wp
    call density_end(this%density)
    call base_states_hash_end(this%hash)
    POP_SUB(base_states_end_states)
    return
  end subroutine base_states_end_states
    
  ! ---------------------------------------------------------
  subroutine base_states_iterator_init(this, that)
    type(base_states_iterator_t), intent(out) :: this
    type(base_states_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_states_iterator_init)
    this%self=>that
    call base_states_hash_init(this%iter, that%hash)
    POP_SUB(base_states_iterator_init)
    return
  end subroutine base_states_iterator_init

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_config_states(this, config, system, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: config
    type(base_states_t),         pointer        :: system
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_config_states)
    call base_states_hash_next(this%iter, config, system, ierr)
    POP_SUB(base_states_iterator_next_config_states)
    return
  end subroutine base_states_iterator_next_config_states

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_config(this, that, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(json_object_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_config)
    call base_states_hash_next(this%iter, that, ierr)
    POP_SUB(base_states_iterator_next_config)
    return
  end subroutine base_states_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_states_iterator_next_system(this, that, ierr)
    type(base_states_iterator_t), intent(inout) :: this
    type(base_states_t),         pointer        :: that
    integer,            optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_states_iterator_next_system)
    call base_states_hash_next(this%iter, that, ierr)
    POP_SUB(base_states_iterator_next_system)
    return
  end subroutine base_states_iterator_next_system

  ! ---------------------------------------------------------
  subroutine base_states_iterator_copy(this, that)
    type(base_states_iterator_t), intent(inout) :: this
    type(base_states_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_states_iterator_copy)
    this%self=>that%self
    call base_states_hash_copy(this%iter, that%iter)
    POP_SUB(base_states_iterator_copy)
    return
  end subroutine base_states_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_states_iterator_end(this)
    type(base_states_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_states_iterator_end)
    nullify(this%self)
    call base_states_hash_end(this%iter)
    POP_SUB(base_states_iterator_end)
    return
  end subroutine base_states_iterator_end

end module base_states_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
