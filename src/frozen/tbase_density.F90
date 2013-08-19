#include "global.h"
#include "template.h"

module TEMPLATE(base_density_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use TEMPLATE(simulation_m), only:         &
    simulation_t !=> TEMPLATE(simulation_t)

  use storage_m, only: &
    operator(+),       &
    operator(-)

  use storage_m, only:     &
    storage_t,             &
    storage_init,          &
    storage_start,         &
    storage_update,        &
    storage_get_size,      &
    storage_get_dimension, &
    storage_get_storage,   &
    storage_copy,          &
    storage_end

  use storage_m, only:          &
    STORAGE_INTRP_OK,           &
    storage_interpolation_t,    &
    storage_interpolation_init, &
    storage_interpolation_eval, &
    storage_interpolation_copy, &
    storage_interpolation_end

#ifdef SUBTEMPLATE_NAME
  use SUBTEMPLATE(m), only:      &
    sub_t   => SUBTEMPLATE(t),   &
    sub_get => SUBTEMPLATE(get)

  use SUBTEMPLATE(list_m), only:         &
    list_t    => SUBTEMPLATE(list_t),    &
    list_init => SUBTEMPLATE(list_init), &
    list_push => SUBTEMPLATE(list_push), &
    list_copy => SUBTEMPLATE(list_copy), &
    list_end  => SUBTEMPLATE(list_end)
#endif

  implicit none

  private
  public ::      &
    operator(+), &
    operator(-)

  public ::                                   &
    TEMPLATE(base_density_init),              &
    TEMPLATE(base_density_start),             &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(base_density_extend),            &
#endif
    TEMPLATE(base_density_update),            &
    TEMPLATE(base_density_get),               &
    TEMPLATE(base_density_get_size),          &
    TEMPLATE(base_density_get_nspin),         &
    TEMPLATE(base_density_get_density),       &
    TEMPLATE(base_density_get_total_density), &
    TEMPLATE(base_density_copy),              &
    TEMPLATE(base_density_end)

  public ::                                    &
    TEMPLATE(base_density_interpolation_init), &
    TEMPLATE(base_density_interpolation_eval), &
    TEMPLATE(base_density_interpolation_copy), &
    TEMPLATE(base_density_interpolation_end)

  type, public :: TEMPLATE(base_density_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(storage_t),     pointer :: total  =>null()
    integer                      :: nspin
    type(storage_t)              :: density
#ifdef SUBTEMPLATE_NAME
    logical                      :: block  = .false.
    type(list_t)                 :: list
#endif
  end type TEMPLATE(base_density_t)

  type, public :: TEMPLATE(base_density_interpolation_t)
    private
    type(TEMPLATE(base_density_t)), pointer :: self =>null()
    type(storage_interpolation_t)           :: intrp
  end type TEMPLATE(base_density_interpolation_t)
 
  interface operator(+)
    module procedure TEMPLATE(base_density_add)
  end interface operator(+)

  interface operator(-)
    module procedure TEMPLATE(base_density_sub)
  end interface operator(-)

  interface TEMPLATE(base_density_get)
    module procedure TEMPLATE(base_density_get_config)
    module procedure TEMPLATE(base_density_get_simulation)
  end interface TEMPLATE(base_density_get)

  interface TEMPLATE(base_density_get_density)
    module procedure TEMPLATE(base_density_get_density_1d)
    module procedure TEMPLATE(base_density_get_density_2d)
  end interface TEMPLATE(base_density_get_density)

  interface TEMPLATE(base_density_interpolation_eval)
    module procedure TEMPLATE(base_density_interpolation_eval_1d)
    module procedure TEMPLATE(base_density_interpolation_eval_2d)
  end interface TEMPLATE(base_density_interpolation_eval)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_init)(this, config)
    type(TEMPLATE(base_density_t)), target, intent(out) :: this
    type(json_object_t),            target, intent(in)  :: config
    !
    integer :: ierr
    !
    this%config=>config
    nullify(this%sim, this%total)
    call json_get(this%config, "SpinComponents", this%nspin, ierr)
    if(ierr/=JSON_OK)this%nspin=1
    if(this%nspin>1)then
      SAFE_ALLOCATE(this%total)
      call storage_init(this%total, 1, 0.0_wp)
    else
      this%total=>this%density
    end if
    call storage_init(this%density, this%nspin, 0.0_wp)
#ifdef SUBTEMPLATE_NAME
    this%block=.false.
    call list_init(this%list)
#endif
    return
  end subroutine TEMPLATE(base_density_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_start)(this, sim)
    type(TEMPLATE(base_density_t)), target, intent(inout) :: this
    type(simulation_t),             target, intent(in)    :: sim
    !
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    if(this%nspin>1)&
      call storage_start(this%total, this%sim, .true.)
    call storage_start(this%density, this%sim, .true.)
    return
  end subroutine TEMPLATE(base_density_start)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_extend)(this, that)
    type(TEMPLATE(base_density_t)), intent(inout) :: this
    type(sub_t),          optional, intent(in)    :: that
    !
    type(json_object_t), pointer :: cnfg
    logical                      :: temp
    integer                      :: ierr
    !
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    nullify(cnfg)
    if(present(that))then
      ASSERT(.not.this%block)
      call sub_get(that, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "temporary", temp, ierr)
      if(ierr/=JSON_OK)temp=.false.
      if(.not.temp)call list_push(this%list, that)
      nullify(cnfg)
    else
      this%block=.true.
    end if
    return
  end subroutine TEMPLATE(base_density_extend)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_update)(this)
    type(TEMPLATE(base_density_t)), intent(inout) :: this
    !
    real(kind=wp), dimension(:,:), pointer :: prho
    real(kind=wp), dimension(:),   pointer :: trho
    !
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    nullify(prho, trho)
    call storage_update(this%density)
    if(this%nspin>1)then
      call storage_get_storage(this%density, prho)
      call storage_get_storage(this%total, trho)
      trho=sum(prho, dim=2)
      call storage_update(this%total)
      nullify(prho, trho)
    end if
    return
  end subroutine TEMPLATE(base_density_update)

  ! ---------------------------------------------------------
  function TEMPLATE(base_density_add)(this, that) result(resl)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    type(TEMPLATE(base_density_t)), intent(in) :: that
    !
    type(TEMPLATE(base_density_t)) :: resl
    !
    resl%density=this%density+that%density
    call TEMPLATE(base_density_update)(resl)
    return
  end function TEMPLATE(base_density_add)

  ! ---------------------------------------------------------
  function TEMPLATE(base_density_sub)(this, that) result(resl)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    type(TEMPLATE(base_density_t)), intent(in) :: that
    !
    type(TEMPLATE(base_density_t)) :: resl
    !
    resl%density=this%density-that%density
    call TEMPLATE(base_density_update)(resl)
    return
  end function TEMPLATE(base_density_sub)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(base_density_get_size)(this) result(that)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    !
    integer :: that
    !
    that=storage_get_size(this%density)
    return
  end function TEMPLATE(base_density_get_size)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(base_density_get_nspin)(this) result(that)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    !
    integer :: that
    !
    that=this%nspin
    return
  end function TEMPLATE(base_density_get_nspin)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_get_config)(this, that)
    type(TEMPLATE(base_density_t)), target, intent(in) :: this
    type(json_object_t),           pointer             :: that
    !
    that=>null()
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(base_density_get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_get_simulation)(this, that)
    type(TEMPLATE(base_density_t)), target, intent(in) :: this
    type(simulation_t),            pointer             :: that
    !
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    return
  end subroutine TEMPLATE(base_density_get_simulation)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_get_density_1d)(this, that)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    real(kind=wp),   dimension(:), pointer     :: that
    !
    ASSERT(this%nspin==1)
    call storage_get_storage(this%density, that)
    return
  end subroutine TEMPLATE(base_density_get_density_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_get_density_2d)(this, that)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    that=>null()
    if(this%nspin>0)&
      call storage_get_storage(this%density, that)
    return
  end subroutine TEMPLATE(base_density_get_density_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_get_total_density)(this, that)
    type(TEMPLATE(base_density_t)), intent(in) :: this
    real(kind=wp),   dimension(:), pointer     :: that
    !
    that=>null()
    if(this%nspin>0)&
      call storage_get_storage(this%total, that)
    return
  end subroutine TEMPLATE(base_density_get_total_density)

  ! ---------------------------------------------------------
  !pure 
  subroutine TEMPLATE(base_density_adjust_spin)(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that
    !
    select case(size(this))
    case(1)
      select case(size(that))
      case(1)
        this=that
      case(2)
        this=sum(that)
      !case(4)
      case default
        ASSERT(.false.)
      end select
    case(2)
      select case(size(that))
      case(1)
        this=0.5_wp*that
      case(2)
        this=that
      !case(4)
      case default
        ASSERT(.false.)
      end select
    case(4)
      select case(size(that))
      !case(1)
      !case(2)
      case(4)
        this=that
      case default
        ASSERT(.false.)
      end select
    case default
      !
    end select
    return
  end subroutine TEMPLATE(base_density_adjust_spin)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_copy)(this, that)
    type(TEMPLATE(base_density_t)), target, intent(out) :: this
    type(TEMPLATE(base_density_t)), target, intent(in)  :: that
    !
    call TEMPLATE(base_density_end)(this)
    if(that%nspin>0)then
      this%config=>that%config
      this%sim=>that%sim
      this%total=>null()
      this%nspin=that%nspin
      if(this%nspin>1)then
        SAFE_ALLOCATE(this%total)
        call storage_copy(this%total, that%total)
      else
        this%total=>this%density
      end if
      call storage_copy(this%density, that%density)
#ifdef SUBTEMPLATE_NAME
      this%block=that%block
      call list_copy(this%list, that%list)
#endif
    end if
    return
  end subroutine TEMPLATE(base_density_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_end)(this)
    type(TEMPLATE(base_density_t)), intent(inout) :: this
    !
#ifdef SUBTEMPLATE_NAME
    call list_end(this%list)
    this%block=.false.
#endif
    if(this%nspin>0)then
      call storage_end(this%density)
      if(this%nspin>1)then
        call storage_end(this%total)
        SAFE_DEALLOCATE_P(this%total)
      end if
    end if
    this%nspin=0
    nullify(this%total, this%sim, this%config)
    return
  end subroutine TEMPLATE(base_density_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_interpolation_init)(this, that)
    type(TEMPLATE(base_density_interpolation_t)), intent(out) :: this
    type(TEMPLATE(base_density_t)),       target, intent(in)  :: that
    !
    integer  :: type, ierr
    !
    this%self=>that
    call json_get(that%config, "interpolation", type, ierr)
    if(ierr==JSON_OK)then
      call storage_interpolation_init(this%intrp, that%density, type)
    else
      call storage_interpolation_init(this%intrp, that%density)
    end if
    return
  end subroutine TEMPLATE(base_density_interpolation_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_interpolation_eval_1d)(this, x, val)
    type(TEMPLATE(base_density_interpolation_t)), intent(in)  :: this
    real(kind=wp),                  dimension(:), intent(in)  :: x
    real(kind=wp),                                intent(out) :: val
    !
    real(kind=wp), dimension(this%self%nspin) :: tvl
    real(kind=wp), dimension(1)               :: tv1
    integer                                   :: ierr
    !
    if(this%self%nspin==1)then
      call storage_interpolation_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRP_OK)val=0.0_wp
    else
      call storage_interpolation_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRP_OK)tvl=0.0_wp
      call TEMPLATE(base_density_adjust_spin)(tv1, tvl)
      val=tv1(1)
    end if
    return
  end subroutine TEMPLATE(base_density_interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_interpolation_eval_2d)(this, x, val)
    type(TEMPLATE(base_density_interpolation_t)), intent(in)  :: this
    real(kind=wp),                  dimension(:), intent(in)  :: x
    real(kind=wp),                  dimension(:), intent(out) :: val
    !
    real(kind=wp), dimension(this%self%nspin) :: tvl
    integer                                   :: ierr
    !
    if(this%self%nspin==size(val))then
      call storage_interpolation_eval(this%intrp, x, val, ierr)
      if(ierr/=STORAGE_INTRP_OK)val=0.0_wp
    else
      call storage_interpolation_eval(this%intrp, x, tvl, ierr)
      if(ierr/=STORAGE_INTRP_OK)tvl=0.0_wp
      call TEMPLATE(base_density_adjust_spin)(val, tvl)
    end if
    return
  end subroutine TEMPLATE(base_density_interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_interpolation_copy)(this, that)
    type(TEMPLATE(base_density_interpolation_t)), intent(out) :: this
    type(TEMPLATE(base_density_interpolation_t)), intent(in)  :: that
    !
    this%self=>that%self
    call storage_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(base_density_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(base_density_interpolation_end)(this)
    type(TEMPLATE(base_density_interpolation_t)), intent(inout) :: this
    !
    call storage_interpolation_end(this%intrp)
    this%self=>null()
    return
  end subroutine TEMPLATE(base_density_interpolation_end)

end module TEMPLATE(base_density_m)

!! Local Variables:
!! mode: f90
!! End:

