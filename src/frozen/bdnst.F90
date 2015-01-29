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

#define HASH_TEMPLATE_NAME bdnst
#define HASH_KEY_TEMPLATE_NAME json
#define HASH_KEY_TYPE_NAME json_object_t
#define HASH_VAL_TEMPLATE_NAME bdnst

#define HASH_INCLUDE_PREFIX
#include "thash.F90"
#undef HASH_INCLUDE_PREFIX

module bdnst_m

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
  public ::       &
    bdnst_init,   &
    bdnst_start,  &
    bdnst_update, &
    bdnst_stop,   &
    bdnst_eval,   &
    bdnst_get,    &
    bdnst_copy,   &
    bdnst_end

#define HASH_INCLUDE_HEADER
#include "thash.F90"
#undef HASH_INCLUDE_HEADER

  type, public :: bdnst_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(storage_t),     pointer :: total  =>null()
    type(storage_t)              :: data
    type(bdnst_hash_t)           :: hash
  end type bdnst_t

  type, public :: bdnst_iterator_t
    private
    type(bdnst_t),      pointer :: self =>null()
    type(bdnst_hash_iterator_t) :: iter
  end type bdnst_iterator_t

  type, public :: bdnst_intrpl_t
    private
    type(bdnst_t), pointer :: self =>null()
    type(storage_intrpl_t) :: intrp
  end type bdnst_intrpl_t
 
  interface bdnst_init
    module procedure bdnst_init_bdnst
    module procedure bdnst_init_build
    module procedure bdnst_iterator_init
    module procedure bdnst_intrpl_init
  end interface bdnst_init

  interface bdnst_get
    module procedure bdnst_get_info
    module procedure bdnst_get_config
    module procedure bdnst_get_simulation
    module procedure bdnst_get_density_1d
    module procedure bdnst_get_density_2d
  end interface bdnst_get

  interface bdnst_next
    module procedure bdnst_iterator_next_config_bdnst
    module procedure bdnst_iterator_next_config
    module procedure bdnst_iterator_next_bdnst
  end interface bdnst_next

  interface bdnst_eval
    module procedure bdnst_intrpl_eval_1d
    module procedure bdnst_intrpl_eval_2d
  end interface bdnst_eval

  interface bdnst_copy
    module procedure bdnst_copy_bdnst
    module procedure bdnst_iterator_copy
    module procedure bdnst_intrpl_copy
  end interface bdnst_copy

  interface bdnst_end
    module procedure bdnst_end_bdnst
    module procedure bdnst_iterator_end
    module procedure bdnst_intrpl_end
  end interface bdnst_end

contains

#define HASH_INCLUDE_BODY
#include "thash.F90"
#undef HASH_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine bdnst_init_bdnst(this, config)
    type(bdnst_t),       target, intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: nspin, ierr
    !
    PUSH_SUB(bdnst_init_bdnst)
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
    call bdnst_hash_init(this%hash)
    POP_SUB(bdnst_init_bdnst)
    return
  end subroutine bdnst_init_bdnst

  ! ---------------------------------------------------------
  subroutine bdnst_init_build(this, that, config)
    type(bdnst_t),       intent(inout) :: this
    type(bdnst_t),       intent(in)    :: that
    type(json_object_t), intent(in)    :: config
    !
    PUSH_SUB(bdnst_init_build)
    call bdnst_hash_set(this%hash, config, that)
    POP_SUB(bdnst_init_build)
    return
  end subroutine bdnst_init_build

  ! ---------------------------------------------------------
  subroutine bdnst_start(this, sim)
    type(bdnst_t),              intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(bdnst_start)
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    if(bdnst_get_nspin(this)>1)&
      call storage_start(this%total, sim)
    call storage_start(this%data, sim)
    POP_SUB(bdnst_start)
    return
  end subroutine bdnst_start

  ! ---------------------------------------------------------
  subroutine bdnst_update(this)
    type(bdnst_t), intent(inout) :: this
    !
    real(kind=wp), dimension(:,:), pointer :: prho
    real(kind=wp), dimension(:),   pointer :: trho
    !
    PUSH_SUB(bdnst_update)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(prho, trho)
    call storage_update(this%data)
    if(bdnst_get_nspin(this)>1)then
      call storage_get(this%data, prho)
      call storage_get(this%total, trho)
      trho=sum(prho, dim=2)
      call storage_update(this%total)
      nullify(prho, trho)
    end if
    POP_SUB(bdnst_update)
    return
  end subroutine bdnst_update

  ! ---------------------------------------------------------
  subroutine bdnst_stop(this)
    type(bdnst_t), intent(inout) :: this
    !
    PUSH_SUB(bdnst_stop)
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(bdnst_get_nspin(this)>1)&
      call storage_stop(this%total)
    call storage_stop(this%data)
    POP_SUB(bdnst_stop)
    return
  end subroutine bdnst_stop

  ! ---------------------------------------------------------
  elemental function bdnst_get_size(this) result(that)
    type(bdnst_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_size(this%data)
    return
  end function bdnst_get_size

  ! ---------------------------------------------------------
  elemental function bdnst_get_nspin(this) result(that)
    type(bdnst_t), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_dimension(this%data)
    return
  end function bdnst_get_nspin

  ! ---------------------------------------------------------
  elemental subroutine bdnst_get_info(this, size, nspin)
    type(bdnst_t),     intent(in)  :: this
    integer, optional, intent(out) :: size
    integer, optional, intent(out) :: nspin
    !
    if(present(size))&
      size=bdnst_get_size(this)
    if(present(nspin))&
      nspin=bdnst_get_nspin(this)
    return
  end subroutine bdnst_get_info

  ! ---------------------------------------------------------
  subroutine bdnst_get_config(this, that)
    type(bdnst_t),        target, intent(in) :: this
    type(json_object_t), pointer             :: that
    !
    PUSH_SUB(bdnst_get_config)
    that=>null()
    if(associated(this%config))&
      that=>this%config
    POP_SUB(bdnst_get_config)
    return
  end subroutine bdnst_get_config

  ! ---------------------------------------------------------
  subroutine bdnst_get_simulation(this, that)
    type(bdnst_t),       target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(bdnst_get_simulation)
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(bdnst_get_simulation)
    return
  end subroutine bdnst_get_simulation

  ! ---------------------------------------------------------
  subroutine bdnst_get_storage(this, that)
    type(bdnst_t),    target, intent(in) :: this
    type(storage_t), pointer             :: that
    !
    PUSH_SUB(bdnst_get_storage)
    that=>this%data
    POP_SUB(bdnst_get_storage)
    return
  end subroutine bdnst_get_storage

  ! ---------------------------------------------------------
  subroutine bdnst_get_density_1d(this, that, total)
    type(bdnst_t),                          intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total
    !
    PUSH_SUB(bdnst_get_density_1d)
    call storage_get(this%data, that)
    if(present(total))&
      call bdnst_get_total_density(this, total)
    POP_SUB(bdnst_get_density_1d)
    return
  end subroutine bdnst_get_density_1d

  ! ---------------------------------------------------------
  subroutine bdnst_get_density_2d(this, that)
    type(bdnst_t),                  intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    PUSH_SUB(density_get_bdnst_2d)
    call storage_get(this%data, that)
    POP_SUB(bdnst_get_density_2d)
    return
  end subroutine bdnst_get_density_2d

  ! ---------------------------------------------------------
  subroutine bdnst_get_total_density(this, that)
    type(bdnst_t),                intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    !
    PUSH_SUB(bdnst_get_total_density)
    call storage_get(this%total, that)
    POP_SUB(bdnst_get_total_density)
    return
  end subroutine bdnst_get_total_density

  ! ---------------------------------------------------------
  !pure 
  subroutine bdnst_adjust_spin_1_n(this, that)
    real(kind=wp),               intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(bdnst_adjust_spin_1_n)
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
    POP_SUB(bdnst_adjust_spin_1_n)
    return
  end subroutine bdnst_adjust_spin_1_n

  ! ---------------------------------------------------------
  !pure 
  subroutine bdnst_adjust_spin_2_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(bdnst_adjust_spin_2_n)
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
    POP_SUB(bdnst_adjust_spin_2_n)
    return
  end subroutine bdnst_adjust_spin_2_n

  ! ---------------------------------------------------------
  !pure 
  subroutine bdnst_adjust_spin_4_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(bdnst_adjust_spin_4_n)
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
    POP_SUB(bdnst_adjust_spin_4_n)
    return
  end subroutine bdnst_adjust_spin_4_n

  ! ---------------------------------------------------------
  !pure 
  subroutine bdnst_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    PUSH_SUB(bdnst_adjust_spin)
    select case(size(this))
    case(1)
      call bdnst_adjust_spin_1_n(this(1), that)
    case(2)
      call bdnst_adjust_spin_2_n(this, that)
    case(4)
      call bdnst_adjust_spin_4_n(this, that)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(bdnst_adjust_spin)
    return
  end subroutine bdnst_adjust_spin

  ! ---------------------------------------------------------
  subroutine bdnst_copy_bdnst(this, that)
    type(bdnst_t), target, intent(inout) :: this
    type(bdnst_t), target, intent(in)    :: that
    !
    PUSH_SUB(bdnst_copy_bdnst)
    if(associated(that%config))then
      this%config=>that%config
      if(associated(that%sim))then
        this%sim=>that%sim
        call storage_copy(this%data, that%data)
        if(bdnst_get_nspin(this)>1)then
          call storage_copy(this%total, that%total)
        else
          this%total=>this%data
        end if
        call bdnst_hash_copy(this%hash, that%hash)
      end if
    end if
    POP_SUB(bdnst_copy_bdnst)
    return
  end subroutine bdnst_copy_bdnst

  ! ---------------------------------------------------------
  subroutine bdnst_end_bdnst(this)
    type(bdnst_t), intent(inout) :: this
    !
    PUSH_SUB(bdnst_end_bdnst)
    if(associated(this%config))then
      call bdnst_hash_end(this%hash)
      if(bdnst_get_nspin(this)>1)then
        call storage_end(this%total)
        SAFE_DEALLOCATE_P(this%total)
      end if
      call storage_end(this%data)
    end if
    nullify(this%total, this%sim, this%config)
    POP_SUB(bdnst_end_bdnst)
    return
  end subroutine bdnst_end_bdnst

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_init(this, that)
    type(bdnst_iterator_t), intent(out) :: this
    type(bdnst_t),  target, intent(in)  :: that
    !
    PUSH_SUB(bdnst_iterator_init)
    This%self=>that
    call bdnst_hash_init(this%iter, that%hash)
    POP_SUB(bdnst_iterator_init)
    return
  end subroutine bdnst_iterator_init

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_next_config_bdnst(this, config, bdnst, ierr)
    type(bdnst_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: config
    type(bdnst_t),         pointer        :: bdnst
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bdnst_iterator_next_config_bdnst)
    call bdnst_hash_next(this%iter, config, bdnst, ierr)
    POP_SUB(bdnst_iterator_next_config_bdnst)
    return
  end subroutine bdnst_iterator_next_config_bdnst

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_next_config(this, that, ierr)
    type(bdnst_iterator_t), intent(inout) :: this
    type(json_object_t),   pointer        :: that
    integer,      optional, intent(out)   :: ierr
    !
    PUSH_SUB(bdnst_iterator_next_config)
    call bdnst_hash_next(this%iter, that, ierr)
    POP_SUB(bdnst_iterator_next_config)
    return
  end subroutine bdnst_iterator_next_config

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_next_bdnst(this, that, ierr)
    type(bdnst_iterator_t), intent(inout) :: this
    type(bdnst_t),         pointer        :: that
    integer,      optional, intent(out)    :: ierr
    !
    PUSH_SUB(bdnst_iterator_next_bdnst)
    call bdnst_hash_next(this%iter, that, ierr)
    POP_SUB(bdnst_iterator_next_bdnst)
    return
  end subroutine bdnst_iterator_next_bdnst

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_copy(this, that)
    type(bdnst_iterator_t), intent(inout) :: this
    type(bdnst_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(bdnst_iterator_copy)
    this%self=>that%self
    call bdnst_hash_copy(this%iter, that%iter)
    POP_SUB(bdnst_iterator_copy)
    return
  end subroutine bdnst_iterator_copy

  ! ---------------------------------------------------------
  subroutine bdnst_iterator_end(this)
    type(bdnst_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(bdnst_iterator_end)
    nullify(this%self)
    call bdnst_hash_end(this%iter)
    POP_SUB(bdnst_iterator_end)
    return
  end subroutine bdnst_iterator_end

  ! ---------------------------------------------------------
  subroutine bdnst_intrpl_init(this, that, type)
    type(bdnst_intrpl_t),  intent(out) :: this
    type(bdnst_t), target, intent(in)  :: that
    integer,     optional, intent(in)  :: type
    !
    PUSH_SUB(bdnst_intrpl_init)
    this%self=>that
    call storage_init(this%intrp, that%data, type)
    POP_SUB(bdnst_intrpl_init)
    return
  end subroutine bdnst_intrpl_init

  ! ---------------------------------------------------------
  subroutine bdnst_intrpl_eval_1d(this, x, val)
    type(bdnst_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val
    !
    real(kind=wp), dimension(bdnst_get_nspin(this%self)) :: tvl
    real(kind=wp), dimension(1)                          :: tv1
    integer                                              :: ierr
    !
    PUSH_SUB(bdnst_intrpl_eval_1d)
    if(bdnst_get_nspin(this%self)==1)then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call bdnst_adjust_spin(tv1, tvl)
      val=tv1(1)
    end if
    POP_SUB(bdnst_intrpl_eval_1d)
    return
  end subroutine bdnst_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine bdnst_intrpl_eval_2d(this, x, val)
    type(bdnst_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val
    !
    real(kind=wp), dimension(bdnst_get_nspin(this%self)) :: tvl
    integer                                              :: ierr
    !
    PUSH_SUB(bdnst_intrpl_eval_2d)
    if(bdnst_get_nspin(this%self)==size(val))then
      call storage_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRPL_OK)val=0.0_wp
    else
      call storage_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRPL_OK)tvl=0.0_wp
      call bdnst_adjust_spin(val, tvl)
    end if
    POP_SUB(bdnst_intrpl_eval_2d)
    return
  end subroutine bdnst_intrpl_eval_2d

  ! ---------------------------------------------------------
  subroutine bdnst_intrpl_copy(this, that)
    type(bdnst_intrpl_t), intent(out) :: this
    type(bdnst_intrpl_t), intent(in)  :: that
    !
    PUSH_SUB(bdnst_intrpl_copy)
    this%self=>that%self
    call storage_copy(this%intrp, that%intrp)
    POP_SUB(bdnst_intrpl_copy)
    return
  end subroutine bdnst_intrpl_copy

  ! ---------------------------------------------------------
  subroutine bdnst_intrpl_end(this)
    type(bdnst_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(bdnst_intrpl_end)
    nullify(this%self)
    call storage_end(this%intrp)
    POP_SUB(bdnst_intrpl_end)
    return
  end subroutine bdnst_intrpl_end

end module bdnst_m

#undef HASH_TEMPLATE_NAME
#undef HASH_KEY_TEMPLATE_NAME
#undef HASH_KEY_TYPE_NAME
#undef HASH_VAL_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:

