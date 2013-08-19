#include "global.h"
#include "template.h"

module TEMPLATE(states_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use TEMPLATE(simulation_m), only:         &
    simulation_t !=> TEMPLATE(simulation_t)

  use TEMPLATE(density_m), only:                              &
    density_t         => TEMPLATE(density_t),                 &
    density_init      => TEMPLATE(density_init),              &
    density_start     => TEMPLATE(density_start),             &
#ifdef SUBTEMPLATE_NAME
    density_extend    => TEMPLATE(density_extend),            &
#endif
    density_update    => TEMPLATE(density_update),            &
    density_get       => TEMPLATE(density_get),               &
    density_get_size  => TEMPLATE(density_get_size),          &
    density_get_nspin => TEMPLATE(density_get_nspin),         &
    density_copy      => TEMPLATE(density_copy),              &
    density_end       => TEMPLATE(density_end)

  use TEMPLATE(density_m), only:                                        &
    density_interpolation_t    => TEMPLATE(density_interpolation_t),    &
    density_interpolation_init => TEMPLATE(density_interpolation_init), &
    density_interpolation_eval => TEMPLATE(density_interpolation_eval), &
    density_interpolation_copy => TEMPLATE(density_interpolation_copy), &
    density_interpolation_end  => TEMPLATE(density_interpolation_end)

#ifdef SUBTEMPLATE_NAME
  use SUBTEMPLATE(m), only:  &
    sub_t => SUBTEMPLATE(t)
#endif

  implicit none

  private
  public ::                     &
    TEMPLATE(states_init),      &
    TEMPLATE(states_start),     &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(states_extend),    &
#endif
    TEMPLATE(states_update),    &
    TEMPLATE(states_get),       &
    TEMPLATE(states_get_size),  &
    TEMPLATE(states_get_nspin), &
    TEMPLATE(states_copy),      &
    TEMPLATE(states_end)

  public ::                              &
    TEMPLATE(states_interpolation_init), &
    TEMPLATE(states_interpolation_eval), &
    TEMPLATE(states_interpolation_copy), &
    TEMPLATE(states_interpolation_end)

  type, public :: TEMPLATE(states_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(density_t)              :: density
  end type TEMPLATE(states_t)

  type, public :: TEMPLATE(states_interpolation_t)
    private
    type(TEMPLATE(states_t)), pointer :: self =>null()
    type(density_interpolation_t)     :: intrp
  end type TEMPLATE(states_interpolation_t)

  interface TEMPLATE(states_get)
    module procedure TEMPLATE(states_get_config)
    module procedure TEMPLATE(states_get_simulation)
    module procedure TEMPLATE(states_get_density)
  end interface TEMPLATE(states_get)

  interface TEMPLATE(states_update)
#ifdef SUBTEMPLATE_NAME
    module procedure TEMPLATE(states_update_build)
#endif
    module procedure TEMPLATE(states_update_finish)
  end interface TEMPLATE(states_update)

  interface TEMPLATE(states_interpolation_eval)
    module procedure TEMPLATE(states_interpolation_eval_1d)
    module procedure TEMPLATE(states_interpolation_eval_2d)
  end interface TEMPLATE(states_interpolation_eval)

contains
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_init)(this, config)
    type(TEMPLATE(states_t)),    intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer :: ierr
    !
    this%config=>config
    nullify(cnfg, this%sim)
    call json_get(this%config, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call density_init(this%density, cnfg)
    nullify(cnfg)
    return
  end subroutine TEMPLATE(states_init)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_start)(this, sim)
    type(TEMPLATE(states_t)),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call density_start(this%density, sim)
    return
  end subroutine TEMPLATE(states_start)
    
#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_extend)(this, that)
    type(TEMPLATE(states_t)), intent(inout) :: this
    type(sub_t),    optional, intent(in)    :: that
    !
    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    call density_extend(this%density, that)
    return
  end subroutine TEMPLATE(states_extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_update_build)(this, that)
    type(TEMPLATE(states_t)), intent(inout) :: this
    type(sub_t),              intent(in)    :: that
    !
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call density_update(this%density, that)
    return
  end subroutine TEMPLATE(states_update_build)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_update_finish)(this)
    type(TEMPLATE(states_t)), intent(inout) :: this
    !
    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call density_update(this%density)
    return
  end subroutine TEMPLATE(states_update_finish)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(states_get_size)(this) result(that)
    type(TEMPLATE(states_t)), intent(in) :: this
    !
    integer :: that
    !
    that=density_get_size(this%density)
    return
  end function TEMPLATE(states_get_size)
    
  ! ---------------------------------------------------------
  elemental function TEMPLATE(states_get_nspin)(this) result(that)
    type(TEMPLATE(states_t)), intent(in) :: this
    !
    integer :: that
    !
    that=density_get_nspin(this%density)
    return
  end function TEMPLATE(states_get_nspin)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_get_config)(this, that)
    type(TEMPLATE(states_t)), target, intent(in) :: this
    type(json_object_t),     pointer             :: that
    !
    that=>null()
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(states_get_config)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_get_simulation)(this, that)
    type(TEMPLATE(states_t)), target, intent(in) :: this
    type(simulation_t),      pointer             :: that
    !
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    return
  end subroutine TEMPLATE(states_get_simulation)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_get_density)(this, that)
    type(TEMPLATE(states_t)), target, intent(in) :: this
    type(density_t),         pointer             :: that
    !
    that=>this%density
    return
  end subroutine TEMPLATE(states_get_density)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_copy)(this, that)
    type(TEMPLATE(states_t)),         intent(out) :: this
    type(TEMPLATE(states_t)), target, intent(in)  :: that
    !
    this%config=>that%config
    this%sim=>that%sim
    call density_copy(this%density, that%density)
    return
  end subroutine TEMPLATE(states_copy)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_end)(this)
    type(TEMPLATE(states_t)), intent(inout) :: this
    !
    call density_end(this%density)
    nullify(this%sim, this%config)
    return
  end subroutine TEMPLATE(states_end)
    
  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_interpolation_init)(this, that, what)
    type(TEMPLATE(states_interpolation_t)), intent(out) :: this
    type(TEMPLATE(states_t)),       target, intent(in)  :: that
    type(density_t),                        intent(in)  :: what
    !
    this%self=>that
    call density_interpolation_init(this%intrp, what)
    return
  end subroutine TEMPLATE(states_interpolation_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_interpolation_eval_1d)(this, x, val)
    type(TEMPLATE(states_interpolation_t)), intent(in)  :: this
    real(kind=wp),            dimension(:), intent(in)  :: x
    real(kind=wp),                          intent(out) :: val
    !
    call density_interpolation_eval(this%intrp, x, val)
    return
  end subroutine TEMPLATE(states_interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_interpolation_eval_2d)(this, x, val)
    type(TEMPLATE(states_interpolation_t)), intent(in)  :: this
    real(kind=wp),            dimension(:), intent(in)  :: x
    real(kind=wp),            dimension(:), intent(out) :: val
    !
    call density_interpolation_eval(this%intrp, x, val)
    return
  end subroutine TEMPLATE(states_interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_interpolation_copy)(this, that)
    type(TEMPLATE(states_interpolation_t)), intent(out) :: this
    type(TEMPLATE(states_interpolation_t)), intent(in)  :: that
    !
    this%self=>that%self
    call density_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(states_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(states_interpolation_end)(this)
    type(TEMPLATE(states_interpolation_t)), intent(inout) :: this
    !
    call density_interpolation_end(this%intrp)
    this%self=>null()
    return
  end subroutine TEMPLATE(states_interpolation_end)

end module TEMPLATE(states_m)

!! Local Variables:
!! mode: f90
!! End:
