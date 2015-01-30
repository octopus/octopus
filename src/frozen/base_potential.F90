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

#define HASH_TEMPLATE_NAME base_potential
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_potential

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_potential_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,   only: JSON_OK, json_object_t, json_get

  use storage_m, only:     &
    storage_t,             &
    storage_init,          &
    storage_start,         &
    storage_update,        &
    storage_stop,          &
    storage_eval,          &
    storage_get,           &
    storage_get_size,      &
    storage_get_dimension, &
    storage_copy,          &
    storage_end

  use storage_m, only: &
    storage_intrpl_t

  use storage_m, only:                             &
    BASE_POTENTIAL_INTRPL_OK => STORAGE_INTRPL_OK, &
    BASE_POTENTIAL_INTRPL_OD => STORAGE_INTRPL_OD, &
    BASE_POTENTIAL_INTRPL_NI => STORAGE_INTRPL_NI

  use simulation_m, only: &
    simulation_t

  use base_system_m, only:     &
    system_t => base_system_t

  implicit none

  private
  public ::                &
    base_potential_init,   &
    base_potential_start,  &
    base_potential_update, &
    base_potential_stop,   &
    base_potential_eval,   &
    base_potential_set,    &
    base_potential_get,    &
    base_potential_copy,   &
    base_potential_end

  public ::                   &
    BASE_POTENTIAL_INTRPL_OK, &
    BASE_POTENTIAL_INTRPL_OD, &
    BASE_POTENTIAL_INTRPL_NI

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_potential_t
    private
    type(json_object_t), pointer :: config =>null()
    type(system_t),      pointer :: sys    =>null()
    type(simulation_t),  pointer :: sim    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(storage_t)              :: data
    type(base_potential_hash_t)  :: hash
  end type base_potential_t

  type, public :: base_potential_iterator_t
    private
    type(base_potential_t),      pointer :: self =>null()
    type(base_potential_hash_iterator_t) :: iter
  end type base_potential_iterator_t

  type, public :: base_potential_intrpl_t
    private
    type(base_potential_t), pointer :: self =>null()
    type(storage_intrpl_t)          :: intrp
  end type base_potential_intrpl_t

  interface base_potential_init
    module procedure base_potential_init_begin
    module procedure base_potential_init_copy
    module procedure base_potential_init_build
    module procedure base_potential_iterator_init
    module procedure base_potential_intrpl_init
  end interface base_potential_init

  interface base_potential_set
    module procedure base_potential_set_energy
  end interface base_potential_set

  interface base_potential_get
    module procedure base_potential_get_info
    module procedure base_potential_get_config
    module procedure base_potential_get_system
    module procedure base_potential_get_simulation
    module procedure base_potential_get_potential_1d
    module procedure base_potential_get_potential_md
    module procedure base_potential_intrpl_get
  end interface base_potential_get

  interface base_potential_next
    module procedure base_potential_iterator_next_config_potential
    module procedure base_potential_iterator_next_config
    module procedure base_potential_iterator_next_potential
  end interface base_potential_next

  interface base_potential_eval
    module procedure base_potential_intrpl_eval_1d
    module procedure base_potential_intrpl_eval_md
  end interface base_potential_eval

  interface base_potential_copy
    module procedure base_potential_copy_potential
    module procedure base_potential_iterator_copy
    module procedure base_potential_intrpl_copy
  end interface base_potential_copy

  interface base_potential_end
    module procedure base_potential_end_potential
    module procedure base_potential_iterator_end
    module procedure base_potential_intrpl_end
  end interface base_potential_end

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_potential_init_begin(this, sys, config)
    type(base_potential_t),      intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: nspin, ierr
    !
    PUSH_SUB(base_potential_init_begin)
    ASSERT(.not.associated(this%config))
    this%config=>config
    this%sys=>sys
    nullify(this%sim)
    this%energy=0.0_wp
    call json_get(this%config, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=1
    call storage_init(this%data, nspin)
    call base_potential_hash_init(this%hash)
    POP_SUB(base_potential_init_begin)
    return
  end subroutine base_potential_init_begin

  ! ---------------------------------------------------------
  subroutine base_potential_init_copy(this, that)
    type(base_potential_t),         intent(out) :: this
    type(base_potential_t), target, intent(in)  :: that
    !
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(system_t),      pointer :: sys    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(storage_t)              :: data

    PUSH_SUB(base_potential_init_copy)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    ASSERT(associated(that%sim))
    call base_potential_init_begin(this, that%sys, that%config)
    call base_potential_start(this, that%sim)
    POP_SUB(base_potential_init_copy)
    return
  end subroutine base_potential_init_copy

  ! ---------------------------------------------------------
  subroutine base_potential_init_build(this, that, config)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config
    !
    PUSH_SUB(base_potential_init_build)
    call base_potential_hash_set(this%hash, config, that)
    POP_SUB(base_potential_init_build)
    return
  end subroutine base_potential_init_build

  ! ---------------------------------------------------------
  subroutine base_potential_start(this, sim)
    type(base_potential_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_potential_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call storage_start(this%data, sim)
    POP_SUB(base_potential_start)
    return
  end subroutine base_potential_start

  ! ---------------------------------------------------------
  subroutine base_potential_update(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_update_finish)
    ASSERT(associated(this%sim))
    call storage_update(this%data)
    POP_SUB(base_potential_update_finish)
    return
  end subroutine base_potential_update

  ! ---------------------------------------------------------
  subroutine base_potential_stop(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)
    POP_SUB(base_potential_stop)
    return
  end subroutine base_potential_stop

  ! ---------------------------------------------------------
  elemental subroutine base_potential_set_energy(this, that)
    type(base_potential_t), intent(inout) :: this
    real(kind=wp),          intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine base_potential_set_energy

  ! ---------------------------------------------------------
  elemental function base_potential_get_size(this) result(np)
    type(base_potential_t), intent(in) :: this
    !
    integer :: np
    !
    np=storage_get_size(this%data)
    return
  end function base_potential_get_size

  ! ---------------------------------------------------------
  elemental function base_potential_get_nspin(this) result(that)
    type(base_potential_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%data)
    return
  end function base_potential_get_nspin

  ! ---------------------------------------------------------
  elemental function base_potential_get_energy(this) result(that)
    type(base_potential_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function base_potential_get_energy

  ! ---------------------------------------------------------
  elemental subroutine base_potential_get_info(this, size, nspin, energy)
    type(base_potential_t),  intent(in)  :: this
    integer,       optional, intent(out) :: size
    integer,       optional, intent(out) :: nspin
    real(kind=wp), optional, intent(out) :: energy
    !
    if(present(size))&
      size=base_potential_get_size(this)
    if(present(nspin))&
      nspin=base_potential_get_nspin(this)
    if(present(energy))&
      energy=base_potential_get_energy(this)
    return
  end subroutine base_potential_get_info

  ! ---------------------------------------------------------
  subroutine base_potential_get_config(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(json_object_t),   pointer             :: that
    !
    PUSH_SUB(base_potential_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_potential_get_config)
    return
  end subroutine base_potential_get_config

  ! ---------------------------------------------------------
  subroutine base_potential_get_system(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(system_t),        pointer             :: that
    !
    PUSH_SUB(base_potential_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(base_potential_get_system)
    return
  end subroutine base_potential_get_system

  ! ---------------------------------------------------------
  subroutine base_potential_get_simulation(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(simulation_t),    pointer             :: that
    !
    PUSH_SUB(base_potential_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_potential_get_simulation)
    return
  end subroutine base_potential_get_simulation

  ! ---------------------------------------------------------
  subroutine base_potential_get_storage(this, that)
    type(base_potential_t), target, intent(in) :: this
    type(storage_t),       pointer             :: that
    !
    PUSH_SUB(base_potential_get_storage)
    that=>this%data
    POP_SUB(base_potential_get_storage)
    return
  end subroutine base_potential_get_storage

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_1d(this, that)
    type(base_potential_t),       intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(base_potential_get_potential_1d)
    call storage_get(this%data, that)
    POP_SUB(base_potential_get_potential_1d)
    return
  end subroutine base_potential_get_potential_1d

  ! ---------------------------------------------------------
  subroutine base_potential_get_potential_md(this, that)
    type(base_potential_t),         intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(base_potential_get_potential_md)
    call storage_get(this%data, that)
    POP_SUB(base_potential_get_potential_md)
    return
  end subroutine base_potential_get_potential_md

  ! ---------------------------------------------------------
  subroutine base_potential_copy_potential(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    !
    PUSH_SUB(base_potential_copy_potential)
    this%config=>that%config
    this%sys=>that%sys
    this%sim=>that%sim
    this%energy=that%energy
    call storage_copy(this%data, that%data)
    call base_potential_hash_copy(this%hash, that%hash)
    POP_SUB(base_potential_copy_potential)
    return
  end subroutine base_potential_copy_potential

  ! ---------------------------------------------------------
  subroutine base_potential_end_potential(this)
    type(base_potential_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_end_potential)
    nullify(this%config, this%sys, this%sim)
    this%energy=0.0_wp
    call storage_end(this%data)
    call base_potential_hash_end(this%hash)
    POP_SUB(base_potential_end_potential)
    return
  end subroutine base_potential_end_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_init(this, that)
    type(base_potential_iterator_t), intent(out) :: this
    type(base_potential_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_potential_iterator_init)
    This%self=>that
    call base_potential_hash_init(this%iter, that%hash)
    POP_SUB(base_potential_iterator_init)
    return
  end subroutine base_potential_iterator_init

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_config_potential(this, config, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(json_object_t),            pointer        :: config
    type(base_potential_t),         pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_config_potential)
    call base_potential_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_potential_iterator_next_config_potential)
    return
  end subroutine base_potential_iterator_next_config_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_config(this, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(json_object_t),            pointer        :: that
    integer,               optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_config)
    call base_potential_hash_next(this%iter, that, ierr)
    POP_SUB(base_potential_iterator_next_config)
    return
  end subroutine base_potential_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_next_potential(this, that, ierr)
    type(base_potential_iterator_t), intent(inout) :: this
    type(base_potential_t),         pointer        :: that
    integer,               optional, intent(out)    :: ierr
    !
    PUSH_SUB(base_potential_iterator_next_potential)
    call base_potential_hash_next(this%iter, that, ierr)
    POP_SUB(base_potential_iterator_next_potential)
    return
  end subroutine base_potential_iterator_next_potential

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_copy(this, that)
    type(base_potential_iterator_t), intent(inout) :: this
    type(base_potential_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_potential_iterator_copy)
    this%self=>that%self
    call base_potential_hash_copy(this%iter, that%iter)
    POP_SUB(base_potential_iterator_copy)
    return
  end subroutine base_potential_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_potential_iterator_end(this)
    type(base_potential_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_iterator_end)
    nullify(this%self)
    call base_potential_hash_end(this%iter)
    POP_SUB(base_potential_iterator_end)
    return
  end subroutine base_potential_iterator_end

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_init(this, that, type)
    type(base_potential_intrpl_t),  intent(out) :: this
    type(base_potential_t), target, intent(in)  :: that
    integer,              optional, intent(in)  :: type
    !
    PUSH_SUB(base_potential_intrpl_init)
    this%self=>that
    call storage_init(this%intrp, that%data, type)
    POP_SUB(base_potential_intrpl_init)
    return
  end subroutine base_potential_intrpl_init

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_eval_1d(this, x, v, ierr)
    type(base_potential_intrpl_t), intent(in)  :: this
    real(kind=wp),   dimension(:), intent(in)  :: x
    real(kind=wp),                 intent(out) :: v
    integer,                       intent(out) :: ierr
    !
    PUSH_SUB(base_potential_intrpl_eval_1d)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(base_potential_intrpl_eval_1d)
    return
  end subroutine base_potential_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_eval_md(this, x, v, ierr)
    type(base_potential_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: v
    integer,                     intent(out) :: ierr
    !
    PUSH_SUB(base_potential_intrpl_eval_md)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(base_potential_intrpl_eval_md)
    return
  end subroutine base_potential_intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_get(this, that)
    type(base_potential_intrpl_t), intent(in) :: this
    type(base_potential_t),       pointer     :: that
    !
    PUSH_SUB(base_potential_intrpl_get)
    nullify(that)
    if(associated(this%self))&
      that=>this%self
    POP_SUB(base_potential_intrpl_get)
    return
  end subroutine base_potential_intrpl_get

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_copy(this, that)
    type(base_potential_intrpl_t), intent(out) :: this
    type(base_potential_intrpl_t), intent(in)  :: that
    !
    PUSH_SUB(base_potential_intrpl_copy)
    this%self=>that%self
    call storage_copy(this%intrp, that%intrp)
    POP_SUB(base_potential_intrpl_copy)
    return
  end subroutine base_potential_intrpl_copy

  ! ---------------------------------------------------------
  subroutine base_potential_intrpl_end(this)
    type(base_potential_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(base_potential_intrpl_end)
    nullify(this%self)
    call storage_end(this%intrp)
    POP_SUB(base_potential_intrpl_end)
    return
  end subroutine base_potential_intrpl_end

end module base_potential_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
