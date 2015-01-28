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

#define HASH_TEMPLATE_NAME bpotn
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bpotn

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bpotn_m

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
    storage_get_storage,   &
    storage_copy,          &
    storage_end

  use storage_m, only: &
    storage_intrpl_t

  use storage_m, only:                    &
    BPOTN_INTRPL_OK => STORAGE_INTRPL_OK, &
    BPOTN_INTRPL_OD => STORAGE_INTRPL_OD, &
    BPOTN_INTRPL_NI => STORAGE_INTRPL_NI

  use simulation_m, only:  &
    simulation_t

  use bsyst_m, only:     &
    system_t => bsyst_t

  implicit none

  private
  public ::           &
    bpotn_init,       &
    bpotn_start,      &
    bpotn_update,     &
    bpotn_stop,       &
    bpotn_eval,       &
    bpotn_get_size,   &
    bpotn_get_nspin,  &
    bpotn_get,        &
    bpotn_set_energy, &
    bpotn_get_energy, &
    bpotn_copy,       &
    bpotn_end

  public ::          &
    BPOTN_INTRPL_OK, &
    BPOTN_INTRPL_OD, &
    BPOTN_INTRPL_NI

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bpotn_t
    private
    type(json_object_t), pointer :: config =>null()
    type(system_t),      pointer :: sys    =>null()
    type(simulation_t),  pointer :: sim    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(storage_t)              :: data
    type(bpotn_hash_t)           :: hash
  end type bpotn_t

  type, public :: bpotn_iterator_t
    private
    type(bpotn_t),      pointer :: self =>null()
    type(bpotn_hash_iterator_t) :: iter
  end type bpotn_iterator_t

  type, public :: bpotn_intrpl_t
    private
    type(bpotn_t), pointer :: self =>null()
    type(storage_intrpl_t) :: intrp
  end type bpotn_intrpl_t

  interface bpotn_init
    module procedure bpotn_init_start
    module procedure bpotn_init_copy
    module procedure bpotn_init_build
    module procedure bpotn_iterator_init
    module procedure bpotn_intrpl_init
  end interface bpotn_init

  interface bpotn_get
    module procedure bpotn_get_config
    module procedure bpotn_get_system
    module procedure bpotn_get_simulation
    module procedure bpotn_get_potential_1d
    module procedure bpotn_get_potential_md
    module procedure bpotn_intrpl_get
  end interface bpotn_get

  interface bpotn_next
    module procedure bpotn_iterator_next_config_bpotn
    module procedure bpotn_iterator_next_config
    module procedure bpotn_iterator_next_bpotn
  end interface bpotn_next

  interface bpotn_eval
    module procedure bpotn_intrpl_eval_1d
    module procedure bpotn_intrpl_eval_md
  end interface bpotn_eval

  interface bpotn_copy
    module procedure bpotn_copy_bpotn
    module procedure bpotn_iterator_copy
    module procedure bpotn_intrpl_copy
  end interface bpotn_copy

  interface bpotn_end
    module procedure bpotn_end_bpotn
    module procedure bpotn_iterator_end
    module procedure bpotn_intrpl_end
  end interface bpotn_end

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bpotn_init_start(this, sys, config)
    type(bpotn_t),               intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: nspin, ierr
    !
    PUSH_SUB(bpotn_init_start)
    ASSERT(.not.associated(this%config))
    this%config=>config
    this%sys=>sys
    nullify(this%sim)
    this%energy=0.0_wp
    call json_get(this%config, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=1
    call storage_init(this%data, nspin)
    call bpotn_hash_init(this%hash)
    POP_SUB(bpotn_init_start)
    return
  end subroutine bpotn_init_start

  ! ---------------------------------------------------------
  subroutine bpotn_init_copy(this, that)
    type(bpotn_t),         intent(out) :: this
    type(bpotn_t), target, intent(in)  :: that
    !
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(system_t),      pointer :: sys    =>null()
    real(kind=wp)                :: energy = 0.0_wp
    type(storage_t)              :: data

    PUSH_SUB(bpotn_init_copy)
    ASSERT(associated(that%config))
    ASSERT(associated(that%sys))
    ASSERT(associated(that%sim))
    call bpotn_init_start(this, that%sys, that%config)
    call bpotn_start(this, that%sim)
    POP_SUB(bpotn_init_copy)
    return
  end subroutine bpotn_init_copy

  ! ---------------------------------------------------------
  subroutine bpotn_init_build(this, that, config)
    type(bpotn_t),       intent(inout) :: this
    type(bpotn_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bpotn_init_build)
    call bpotn_hash_set(this%hash, config, that)
    POP_SUB(bpotn_init_build)
    return
  end subroutine bpotn_init_build

  ! ---------------------------------------------------------
  subroutine bpotn_start(this, sim)
    type(bpotn_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(bpotn_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call storage_start(this%data, sim)
    POP_SUB(bpotn_start)
    return
  end subroutine bpotn_start

  ! ---------------------------------------------------------
  subroutine bpotn_update(this)
    type(bpotn_t), intent(inout) :: this
    !
    PUSH_SUB(bpotn_update_finish)
    ASSERT(associated(this%sim))
    call storage_update(this%data)
    POP_SUB(bpotn_update_finish)
    return
  end subroutine bpotn_update

  ! ---------------------------------------------------------
  subroutine bpotn_stop(this)
    type(bpotn_t), intent(inout) :: this
    !
    PUSH_SUB(bpotn_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call storage_stop(this%data)
    POP_SUB(bpotn_stop)
    return
  end subroutine bpotn_stop

  ! ---------------------------------------------------------
  elemental function bpotn_get_size(this) result(np)
    type(bpotn_t), intent(in) :: this
    !
    integer :: np
    !
    np=storage_get_size(this%data)
    return
  end function bpotn_get_size

  ! ---------------------------------------------------------
  elemental function bpotn_get_nspin(this) result(that)
    type(bpotn_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%data)
    return
  end function bpotn_get_nspin

  ! ---------------------------------------------------------
  elemental function bpotn_get_energy(this) result(that)
    type(bpotn_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function bpotn_get_energy

  ! ---------------------------------------------------------
  subroutine bpotn_get_config(this, that)
    type(bpotn_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bpotn_get_config)
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bpotn_get_config)
    return
  end subroutine bpotn_get_config

  ! ---------------------------------------------------------
  subroutine bpotn_get_system(this, that)
    type(bpotn_t),   target, intent(in) :: this
    type(system_t), pointer             :: that
    !
    PUSH_SUB(bpotn_get_system)
    nullify(that)
    if(associated(this%sys))&
      that=>this%sys
    POP_SUB(bpotn_get_system)
    return
  end subroutine bpotn_get_system

  ! ---------------------------------------------------------
  subroutine bpotn_get_simulation(this, that)
    type(bpotn_t),       target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(bpotn_get_simulation)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(bpotn_get_simulation)
    return
  end subroutine bpotn_get_simulation

  ! ---------------------------------------------------------
  subroutine bpotn_get_storage(this, that)
    type(bpotn_t),    target, intent(in) :: this
    type(storage_t), pointer             :: that
    !
    PUSH_SUB(bpotn_get_storage)
    that=>this%data
    POP_SUB(bpotn_get_storage)
    return
  end subroutine bpotn_get_storage

  ! ---------------------------------------------------------
  subroutine bpotn_get_potential_1d(this, that)
    type(bpotn_t),                intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(bpotn_get_potential_1d)
    call storage_get_storage(this%data, that)
    POP_SUB(bpotn_get_potential_1d)
    return
  end subroutine bpotn_get_potential_1d

  ! ---------------------------------------------------------
  subroutine bpotn_get_potential_md(this, that)
    type(bpotn_t),                  intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(bpotn_get_potential_md)
    call storage_get_storage(this%data, that)
    POP_SUB(bpotn_get_potential_md)
    return
  end subroutine bpotn_get_potential_md

  ! ---------------------------------------------------------
  elemental subroutine bpotn_set_energy(this, that)
    type(bpotn_t), intent(inout) :: this
    real(kind=wp),     intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine bpotn_set_energy

  ! ---------------------------------------------------------
  subroutine bpotn_copy_bpotn(this, that)
    type(bpotn_t), intent(inout) :: this
    type(bpotn_t), intent(in)    :: that
    !
    PUSH_SUB(bpotn_copy_bpotn)
    this%config=>that%config
    this%sys=>that%sys
    this%sim=>that%sim
    this%energy=that%energy
    call storage_copy(this%data, that%data)
    call bpotn_hash_copy(this%hash, that%hash)
    POP_SUB(bpotn_copy_bpotn)
    return
  end subroutine bpotn_copy_bpotn

  ! ---------------------------------------------------------
  subroutine bpotn_end_bpotn(this)
    type(bpotn_t), intent(inout) :: this
    !
    PUSH_SUB(bpotn_end_bpotn)
    nullify(this%config, this%sys, this%sim)
    this%energy=0.0_wp
    call storage_end(this%data)
    call bpotn_hash_end(this%hash)
    POP_SUB(bpotn_end_bpotn)
    return
  end subroutine bpotn_end_bpotn

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_init(this, that)
    type(bpotn_iterator_t), intent(out) :: this
    type(bpotn_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bpotn_iterator_init)
    This%self=>that
    call bpotn_hash_init(this%iter, that%hash)
    POP_SUB(bpotn_iterator_init)
    return
  end subroutine bpotn_iterator_init

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_next_config_bpotn(this, config, bpotn, ierr)
    type(bpotn_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bpotn_t),         pointer        :: bpotn
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bpotn_iterator_next_config_bpotn)
    call bpotn_hash_next(this%iter, config, bpotn, ierr)
    POP_SUB(bpotn_iterator_next_config_bpotn)
    return
  end subroutine bpotn_iterator_next_config_bpotn

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_next_config(this, that, ierr)
    type(bpotn_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bpotn_iterator_next_config)
    call bpotn_hash_next(this%iter, that, ierr)
    POP_SUB(bpotn_iterator_next_config)
    return
  end subroutine bpotn_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_next_bpotn(this, that, ierr)
    type(bpotn_iterator_t), intent(inout) :: this
    type(bpotn_t),         pointer        :: that
    integer,     optional, intent(out)    :: ierr
    !
    PUSH_SUB(bpotn_iterator_next_bpotn)
    call bpotn_hash_next(this%iter, that, ierr)
    POP_SUB(bpotn_iterator_next_bpotn)
    return
  end subroutine bpotn_iterator_next_bpotn

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_copy(this, that)
    type(bpotn_iterator_t), intent(inout) :: this
    type(bpotn_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bpotn_iterator_copy)
    this%self=>that%self
    call bpotn_hash_copy(this%iter, that%iter)
    POP_SUB(bpotn_iterator_copy)
    return
  end subroutine bpotn_iterator_copy

  ! ---------------------------------------------------------
  subroutine bpotn_iterator_end(this)
    type(bpotn_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bpotn_iterator_end)
    nullify(this%self)
    call bpotn_hash_end(this%iter)
    POP_SUB(bpotn_iterator_end)
    return
  end subroutine bpotn_iterator_end

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_init(this, that)
    type(bpotn_intrpl_t),  intent(out) :: this
    type(bpotn_t), target, intent(in)  :: that
    !
    PUSH_SUB(bpotn_intrpl_init)
    this%self=>that
    call storage_init(this%intrp, that%data)
    POP_SUB(bpotn_intrpl_init)
    return
  end subroutine bpotn_intrpl_init

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_eval_1d(this, x, v, ierr)
    type(bpotn_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    integer,                     intent(out) :: ierr
    !
    PUSH_SUB(bpotn_intrpl_eval_1d)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(bpotn_intrpl_eval_1d)
    return
  end subroutine bpotn_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_eval_md(this, x, v, ierr)
    type(bpotn_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: v
    integer,                     intent(out) :: ierr
    !
    PUSH_SUB(bpotn_intrpl_eval_md)
    call storage_eval(this%intrp, x, v, ierr)
    POP_SUB(bpotn_intrpl_eval_md)
    return
  end subroutine bpotn_intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_get(this, that)
    type(bpotn_intrpl_t), intent(in) :: this
    type(bpotn_t),       pointer     :: that
    !
    PUSH_SUB(bpotn_intrpl_get)
    nullify(that)
    if(associated(this%self))&
      that=>this%self
    POP_SUB(bpotn_intrpl_get)
    return
  end subroutine bpotn_intrpl_get

!!$  ! ---------------------------------------------------------
!!$  subroutine bpotn_intrpl_set(this, that)
!!$    type(bpotn_intrpl_t), intent(inout) :: this
!!$    type(basis_t),                             intent(in)    :: that
!!$    !
!!$    P1USH_SUB(bpotn_intrpl_set)
!!$    call storage_intrpl_set(this%intrp, that)
!!$    P1OP_SUB(bpotn_intrpl_set)
!!$    return
!!$  end subroutine bpotn_intrpl_set

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_copy(this, that)
    type(bpotn_intrpl_t), intent(out) :: this
    type(bpotn_intrpl_t), intent(in)  :: that
    !
    PUSH_SUB(bpotn_intrpl_copy)
    this%self=>that%self
    call storage_copy(this%intrp, that%intrp)
    POP_SUB(bpotn_intrpl_copy)
    return
  end subroutine bpotn_intrpl_copy

  ! ---------------------------------------------------------
  subroutine bpotn_intrpl_end(this)
    type(bpotn_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(bpotn_intrpl_end)
    nullify(this%self)
    call storage_end(this%intrp)
    POP_SUB(bpotn_intrpl_end)
    return
  end subroutine bpotn_intrpl_end

end module bpotn_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
