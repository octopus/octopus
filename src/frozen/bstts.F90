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

#define HASH_TEMPLATE_NAME bstts
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bstts

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bstts_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,  only: JSON_OK, json_object_t, json_get

  use simulation_m, only: &
    simulation_t

  use bdnst_m, only:                      &
    density_t         => bdnst_t,         &
    density_init      => bdnst_init,      &
    density_start     => bdnst_start,     &
    density_update    => bdnst_update,    &
    density_get       => bdnst_get,       &
    density_get_size  => bdnst_get_size,  &
    density_get_nspin => bdnst_get_nspin, &
    density_copy      => bdnst_copy,      &
    density_end       => bdnst_end

  implicit none

  private
  public ::            &
    bstts_init,       &
    bstts_start,      &
    bstts_update,     &
    bstts_set,        &
    bstts_get,        &
    bstts_get_charge, &
    bstts_set_charge, &
    bstts_copy,       &
    bstts_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bstts_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    real(kind=wp)                :: charge = 0.0_wp
    type(density_t)              :: density
    type(bstts_hash_t)           :: hash
  end type bstts_t

  interface bstts_init
    module procedure bstts_init_bstts
    module procedure bstts_init_build
  end interface bstts_init

  interface bstts_set
    module procedure bstts_set_simulation
  end interface bstts_set

  interface bstts_get
    module procedure bstts_get_config
    module procedure bstts_get_simulation
    module procedure bstts_get_density
  end interface bstts_get

contains
    
#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bstts_init_bstts(this, config)
    type(bstts_t),               intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(bstts_init_bstts)
    this%config=>config
    nullify(cnfg, this%sim)
    call json_get(this%config, "charge", this%charge, ierr)
    if(ierr/=JSON_OK)this%charge=0.0_wp
    call json_get(this%config, "density", cnfg, ierr)
    if(ierr==JSON_OK)call density_init(this%density, cnfg)
    nullify(cnfg)
    call bstts_hash_init(this%hash)
    POP_SUB(bstts_init_bstts)
    return
  end subroutine bstts_init_bstts
    
  ! ---------------------------------------------------------
  subroutine bstts_init_build(this, that, config)
    type(bstts_t),       intent(inout) :: this
    type(bstts_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bstts_init_build)
    this%charge=this%charge+bstts_get_charge(that)
    call density_init(this%density, that%density, config)
    call bstts_hash_set(this%hash, config, that)
    POP_SUB(bstts_init_build)
    return
  end subroutine bstts_init_build
    
  ! ---------------------------------------------------------
  subroutine bstts_start(this, sim)
    type(bstts_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(bstts_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call density_start(this%density, sim)
    POP_SUB(bstts_start)
    return
  end subroutine bstts_start
    
  ! ---------------------------------------------------------
  subroutine bstts_update(this)
    type(bstts_t), intent(inout) :: this
    !
    PUSH_SUB(bstts_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call density_update(this%density)
    POP_SUB(bstts_update)
    return
  end subroutine bstts_update

  ! ---------------------------------------------------------
  elemental subroutine bstts_set_charge(this, that)
    type(bstts_t), intent(inout) :: this
    real(kind=wp), intent(in)    :: that
    !
    this%charge=that
    return
  end subroutine bstts_set_charge
    
  ! ---------------------------------------------------------
  subroutine bstts_set_simulation(this, that)
    type(bstts_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: that
    !
    PUSH_SUB(bstts_set_simulation)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>that
    POP_SUB(bstts_set_simulation)
    return
  end subroutine bstts_set_simulation
    
  ! ---------------------------------------------------------
  elemental function bstts_get_charge(this) result(that)
    type(bstts_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%charge
    return
  end function bstts_get_charge
    
  ! ---------------------------------------------------------
  subroutine bstts_get_config(this, that)
    type(bstts_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bstts_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bstts_get_config)
    return
  end subroutine bstts_get_config
    
  ! ---------------------------------------------------------
  subroutine bstts_get_simulation(this, that)
    type(bstts_t),       target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(bstts_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(bstts_get_simulation)
    return
  end subroutine bstts_get_simulation
    
  ! ---------------------------------------------------------
  subroutine bstts_get_density(this, that)
    type(bstts_t),    target, intent(in) :: this
    type(density_t), pointer             :: that
    !
    PUSH_SUB(bstts_get_density)
    that=>this%density
    POP_SUB(bstts_get_density)
    return
  end subroutine bstts_get_density
    
  ! ---------------------------------------------------------
  subroutine bstts_copy(this, that)
    type(bstts_t),         intent(out) :: this
    type(bstts_t), target, intent(in)  :: that
    !
    PUSH_SUB(bstts_copy)
    this%config=>that%config
    this%sim=>that%sim
    this%charge=that%charge
    call density_copy(this%density, that%density)
    call bstts_hash_copy(this%hash, that%hash)
    POP_SUB(bstts_copy)
    return
  end subroutine bstts_copy
    
  ! ---------------------------------------------------------
  subroutine bstts_end(this)
    type(bstts_t), intent(inout) :: this
    !
    PUSH_SUB(bstts_end)
    nullify(this%config, this%sim)
    this%charge=0.0_wp
    call density_end(this%density)
    call bstts_hash_end(this%hash)
    POP_SUB(bstts_end)
    return
  end subroutine bstts_end
    
end module bstts_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
