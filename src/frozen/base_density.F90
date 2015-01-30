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

#define HASH_TEMPLATE_NAME base_density
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME base_density

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module base_density_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,   only: operator(==), json_hash
  use kinds_m,  only: wp

  use json_m,   only: JSON_OK, json_object_t, json_get

  use simulation_m, only: &
    simulation_t

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
    STORAGE_INTRPL_OK, &
    storage_intrpl_t

  implicit none

  private
  public ::              &
    base_density_init,   &
    base_density_start,  &
    base_density_update, &
    base_density_stop,   &
    base_density_eval,   &
    base_density_get,    &
    base_density_copy,   &
    base_density_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: base_density_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(storage_t),     pointer :: total  =>null()
    type(storage_t)              :: data
    type(base_density_hash_t)    :: hash
  end type base_density_t

  type, public :: base_density_iterator_t
    private
    type(base_density_t),      pointer :: self =>null()
    type(base_density_hash_iterator_t) :: iter
  end type base_density_iterator_t

  type, public :: base_density_intrpl_t
    private
    type(base_density_t), pointer :: self =>null()
    type(storage_intrpl_t)        :: intrp
  end type base_density_intrpl_t
 
  interface base_density_init
    module procedure base_density_init_begin
    module procedure base_density_init_build
    module procedure base_density_iterator_init
    module procedure base_density_intrpl_init
  end interface base_density_init

  interface base_density_get
    module procedure base_density_get_info
    module procedure base_density_get_config
    module procedure base_density_get_simulation
    module procedure base_density_get_density_1d
    module procedure base_density_get_density_2d
  end interface base_density_get

  interface base_density_next
    module procedure base_density_iterator_next_config_density
    module procedure base_density_iterator_next_config
    module procedure base_density_iterator_next_density
  end interface base_density_next

  interface base_density_eval
    module procedure base_density_intrpl_eval_1d
    module procedure base_density_intrpl_eval_2d
  end interface base_density_eval

  interface base_density_copy
    module procedure base_density_copy_density
    module procedure base_density_iterator_copy
    module procedure base_density_intrpl_copy
  end interface base_density_copy

  interface base_density_end
    module procedure base_density_end_density
    module procedure base_density_iterator_end
    module procedure base_density_intrpl_end
  end interface base_density_end

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine base_density_init_begin(this, config)
    type(base_density_t),       target, intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: nspin, ierr
    !
    PUSH_SUB(base_density_init_begin)
    this%config=>config
    nullify(this%sim, this%total)
    call json_get(this%config, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=1
    if(nspin>1)then
      SAFE_ALLOCATE(this%total)
      call storage_init(this%total)
    else
      this%total=>this%data
    end if
    call storage_init(this%data, nspin)
    call base_density_hash_init(this%hash)
    POP_SUB(base_density_init_begin)
    return
  end subroutine base_density_init_begin

  ! ---------------------------------------------------------
  subroutine base_density_init_build(this, that, config)
    type(base_density_t),       intent(inout) :: this
    type(base_density_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(base_density_init_build)
    call base_density_hash_set(this%hash, config, that)
    POP_SUB(base_density_init_build)
    return
  end subroutine base_density_init_build

  ! ---------------------------------------------------------
  subroutine base_density_start(this, sim)
    type(base_density_t),       intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(base_density_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    if(base_density_get_nspin(this)>1)&
      call storage_start(this%total, sim)
    call storage_start(this%data, sim)
    POP_SUB(base_density_start)
    return
  end subroutine base_density_start

  ! ---------------------------------------------------------
  subroutine base_density_update(this)
    type(base_density_t), intent(inout) :: this
    !
    real(kind=wp), dimension(:,:), pointer :: prho
    real(kind=wp), dimension(:),   pointer :: trho
    !
    PUSH_SUB(base_density_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(prho, trho)
    call storage_update(this%data)
    if(base_density_get_nspin(this)>1)then
      call storage_get(this%data, prho)
      call storage_get(this%total, trho)
      trho=sum(prho, dim=2)
      call storage_update(this%total)
      nullify(prho, trho)
    end if
    POP_SUB(base_density_update)
    return
  end subroutine base_density_update

  ! ---------------------------------------------------------
  subroutine base_density_stop(this)
    type(base_density_t), intent(inout) :: this
    !
    PUSH_SUB(base_density_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(base_density_get_nspin(this)>1)&
      call storage_stop(this%total)
    call storage_stop(this%data)
    POP_SUB(base_density_stop)
    return
  end subroutine base_density_stop

  ! ---------------------------------------------------------
  elemental function base_density_get_size(this) result(that)
    type(base_density_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_size(this%data)
    return
  end function base_density_get_size

  ! ---------------------------------------------------------
  elemental function base_density_get_nspin(this) result(that)
    type(base_density_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%data)
    return
  end function base_density_get_nspin

  ! ---------------------------------------------------------
  elemental subroutine base_density_get_info(this, size, nspin)
    type(base_density_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin
    !
    if(present(size))&
      size=base_density_get_size(this)
    if(present(nspin))&
      nspin=base_density_get_nspin(this)
    return
  end subroutine base_density_get_info

  ! ---------------------------------------------------------
  subroutine base_density_get_config(this, that)
    type(base_density_t), target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(base_density_get_config)
    that=>null()
    if(associated(this%config))&
      that=>this%config
    POP_SUB(base_density_get_config)
    return
  end subroutine base_density_get_config

  ! ---------------------------------------------------------
  subroutine base_density_get_simulation(this, that)
    type(base_density_t), target, intent(in) :: this
    type(simulation_t),  pointer             :: that
    !
    PUSH_SUB(base_density_get_simulation)
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(base_density_get_simulation)
    return
  end subroutine base_density_get_simulation

  ! ---------------------------------------------------------
  subroutine base_density_get_storage(this, that)
    type(base_density_t), target, intent(in) :: this
    type(storage_t),     pointer             :: that
    !
    PUSH_SUB(base_density_get_storage)
    that=>this%data
    POP_SUB(base_density_get_storage)
    return
  end subroutine base_density_get_storage

  ! ---------------------------------------------------------
  subroutine base_density_get_density_1d(this, that, total)
    type(base_density_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total
    !
    PUSH_SUB(base_density_get_density_1d)
    call storage_get(this%data, that)
    if(present(total))&
      call base_density_get_total_density(this, total)
    POP_SUB(base_density_get_density_1d)
    return
  end subroutine base_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine base_density_get_density_2d(this, that)
    type(base_density_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(density_get_base_density_2d)
    call storage_get(this%data, that)
    POP_SUB(base_density_get_density_2d)
    return
  end subroutine base_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine base_density_get_total_density(this, that)
    type(base_density_t),         intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(base_density_get_total_density)
    call storage_get(this%total, that)
    POP_SUB(base_density_get_total_density)
    return
  end subroutine base_density_get_total_density

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_1_n(this, that)
    real(kind=wp),               intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(base_density_adjust_spin_1_n)
    select case(size(that))
    case(1)
      this=that(1)
    case(2)
      this=sum(that)
    case(4)
      ASSERT(.false.)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(base_density_adjust_spin_1_n)
    return
  end subroutine base_density_adjust_spin_1_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_2_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(base_density_adjust_spin_2_n)
    select case(size(that))
    case(1)
      this=0.5_wp*that(1)
    case(2)
      this=that
    case(4)
      ASSERT(.false.)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(base_density_adjust_spin_2_n)
    return
  end subroutine base_density_adjust_spin_2_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin_4_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(base_density_adjust_spin_4_n)
    select case(size(that))
    case(1)
      ASSERT(.false.)
    case(2)
      ASSERT(.false.)
    case(4)
      this=that
    case default
      ASSERT(.false.)
    end select
    POP_SUB(base_density_adjust_spin_4_n)
    return
  end subroutine base_density_adjust_spin_4_n

  ! ---------------------------------------------------------
  !pure 
  subroutine base_density_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(base_density_adjust_spin)
    select case(size(this))
    case(1)
      call base_density_adjust_spin_1_n(this(1), that)
    case(2)
      call base_density_adjust_spin_2_n(this, that)
    case(4)
      call base_density_adjust_spin_4_n(this, that)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(base_density_adjust_spin)
    return
  end subroutine base_density_adjust_spin

  ! ---------------------------------------------------------
  subroutine base_density_copy_density(this, that)
    type(base_density_t), target, intent(inout) :: this
    type(base_density_t), target, intent(in)    :: that
    !
    PUSH_SUB(base_density_copy_density)
    if(associated(that%config))then
      this%config=>that%config
      if(associated(that%sim))then
        this%sim=>that%sim
        call storage_copy(this%data, that%data)
        if(base_density_get_nspin(this)>1)then
          call storage_copy(this%total, that%total)
        else
          this%total=>this%data
        end if
        call base_density_hash_copy(this%hash, that%hash)
      end if
    end if
    POP_SUB(base_density_copy_density)
    return
  end subroutine base_density_copy_density

  ! ---------------------------------------------------------
  subroutine base_density_end_density(this)
    type(base_density_t), intent(inout) :: this
    !
    PUSH_SUB(base_density_end_density)
    if(associated(this%config))then
      call base_density_hash_end(this%hash)
      if(base_density_get_nspin(this)>1)then
        call storage_end(this%total)
        SAFE_DEALLOCATE_P(this%total)
      end if
      call storage_end(this%data)
    end if
    nullify(this%total, this%sim, this%config)
    POP_SUB(base_density_end_density)
    return
  end subroutine base_density_end_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_init(this, that)
    type(base_density_iterator_t), intent(out) :: this
    type(base_density_t),  target, intent(in)  :: that
    !
    PUSH_SUB(base_density_iterator_init)
    This%self=>that
    call base_density_hash_init(this%iter, that%hash)
    POP_SUB(base_density_iterator_init)
    return
  end subroutine base_density_iterator_init

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_config_density(this, config, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(json_object_t),          pointer        :: config
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_density_iterator_next_config_density)
    call base_density_hash_next(this%iter, config, that, ierr)
    POP_SUB(base_density_iterator_next_config_density)
    return
  end subroutine base_density_iterator_next_config_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_config(this, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(json_object_t),          pointer        :: that
    integer,             optional, intent(out)   :: ierr
    !
    PUSH_SUB(base_density_iterator_next_config)
    call base_density_hash_next(this%iter, that, ierr)
    POP_SUB(base_density_iterator_next_config)
    return
  end subroutine base_density_iterator_next_config

  ! ---------------------------------------------------------
  subroutine base_density_iterator_next_density(this, that, ierr)
    type(base_density_iterator_t), intent(inout) :: this
    type(base_density_t),         pointer        :: that
    integer,             optional, intent(out)    :: ierr
    !
    PUSH_SUB(base_density_iterator_next_density)
    call base_density_hash_next(this%iter, that, ierr)
    POP_SUB(base_density_iterator_next_density)
    return
  end subroutine base_density_iterator_next_density

  ! ---------------------------------------------------------
  subroutine base_density_iterator_copy(this, that)
    type(base_density_iterator_t), intent(inout) :: this
    type(base_density_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(base_density_iterator_copy)
    this%self=>that%self
    call base_density_hash_copy(this%iter, that%iter)
    POP_SUB(base_density_iterator_copy)
    return
  end subroutine base_density_iterator_copy

  ! ---------------------------------------------------------
  subroutine base_density_iterator_end(this)
    type(base_density_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(base_density_iterator_end)
    nullify(this%self)
    call base_density_hash_end(this%iter)
    POP_SUB(base_density_iterator_end)
    return
  end subroutine base_density_iterator_end

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_init(this, that, type)
    type(base_density_intrpl_t),  intent(out) :: this
    type(base_density_t), target, intent(in)  :: that
    integer,            optional, intent(in)  :: type
    !
    PUSH_SUB(base_density_intrpl_init)
    this%self=>that
    call storage_init(this%intrp, that%data, type)
    POP_SUB(base_density_intrpl_init)
    return
  end subroutine base_density_intrpl_init

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_eval_1d(this, x, val)
    type(base_density_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val
    !
    real(kind=wp), dimension(base_density_get_nspin(this%self)) :: tvl
    real(kind=wp), dimension(1)                                 :: tv1
    integer                                                     :: ierr
    !
    PUSH_SUB(base_density_intrpl_eval_1d)
    if(base_density_get_nspin(this%self)==1)then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call base_density_adjust_spin(tv1, tvl)
      val=tv1(1)
    end if
    POP_SUB(base_density_intrpl_eval_1d)
    return
  end subroutine base_density_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_eval_2d(this, x, val)
    type(base_density_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val
    !
    real(kind=wp), dimension(base_density_get_nspin(this%self)) :: tvl
    integer                                                     :: ierr
    !
    PUSH_SUB(base_density_intrpl_eval_2d)
    if(base_density_get_nspin(this%self)==size(val))then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call base_density_adjust_spin(val, tvl)
    end if
    POP_SUB(base_density_intrpl_eval_2d)
    return
  end subroutine base_density_intrpl_eval_2d

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_copy(this, that)
    type(base_density_intrpl_t), intent(out) :: this
    type(base_density_intrpl_t), intent(in)  :: that
    !
    PUSH_SUB(base_density_intrpl_copy)
    this%self=>that%self
    call storage_copy(this%intrp, that%intrp)
    POP_SUB(base_density_intrpl_copy)
    return
  end subroutine base_density_intrpl_copy

  ! ---------------------------------------------------------
  subroutine base_density_intrpl_end(this)
    type(base_density_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(base_density_intrpl_end)
    nullify(this%self)
    call storage_end(this%intrp)
    POP_SUB(base_density_intrpl_end)
    return
  end subroutine base_density_intrpl_end

end module base_density_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

