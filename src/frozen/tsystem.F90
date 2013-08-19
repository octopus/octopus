#include "global.h"
#include "template.h"

module TEMPLATE(system_m)

  use global_m
  use messages_m
  use profiling_m

  use space_m, only: space_t, operator(==), space_init, space_copy, space_end
  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use TEMPLATE(density_m), only:      &
    density_t => TEMPLATE(density_t)

  use TEMPLATE(geo_m), only:        &
    geo_t    => TEMPLATE(geo_t),    &
    geo_init => TEMPLATE(geo_init), &
    geo_get  => TEMPLATE(geo_get),  &
    geo_copy => TEMPLATE(geo_copy), &
    geo_end  => TEMPLATE(geo_end)

  use TEMPLATE(geo_m), only:                          &
    geo_iterator_init => TEMPLATE(geo_iterator_init)

  use TEMPLATE(geometry_m), only:       &
    geometry_t !=> TEMPLATE(geometry_t)

  use TEMPLATE(simulation_m), only:         &
    simulation_t !=> TEMPLATE(simulation_t)

  use TEMPLATE(states_m), only:               &
    states_t      => TEMPLATE(states_t),      &
    states_init   => TEMPLATE(states_init),   &
    states_start  => TEMPLATE(states_start),  &
    states_update => TEMPLATE(states_update), &
    states_get    => TEMPLATE(states_get),    &
    states_copy   => TEMPLATE(states_copy),   &
    states_end    => TEMPLATE(states_end)

  use TEMPLATE(geo_m), only:                                       &
    TEMPLATE(SYSTEM_ITER_OK)       => TEMPLATE(GEO_ITER_OK),       &
    TEMPLATE(SYSTEM_ITER_END)      => TEMPLATE(GEO_ITER_END),      & 
    TEMPLATE(system_iterator_t)    => TEMPLATE(geo_iterator_t),    &
    TEMPLATE(system_iterator_next) => TEMPLATE(geo_iterator_next), &
    TEMPLATE(system_iterator_copy) => TEMPLATE(geo_iterator_copy), &
    TEMPLATE(system_iterator_end)  => TEMPLATE(geo_iterator_end)

  use TEMPLATE(states_m), only:                                       &
    states_interpolation_t    => TEMPLATE(states_interpolation_t),    &
    states_interpolation_init => TEMPLATE(states_interpolation_init), &
    states_interpolation_eval => TEMPLATE(states_interpolation_eval), &
    states_interpolation_copy => TEMPLATE(states_interpolation_copy), &
    states_interpolation_end  => TEMPLATE(states_interpolation_end)

#ifdef SUBTEMPLATE_NAME
  use TEMPLATE(geo_m), only:            &
    geo_extend => TEMPLATE(geo_extend)

  use TEMPLATE(states_m), only:               &
    states_extend => TEMPLATE(states_extend)

  use SUBTEMPLATE(m), only:      &
    sub_t   => SUBTEMPLATE(t),   &
    sub_get => SUBTEMPLATE(get)
#endif

  implicit none

  private
  public ::                  &
    TEMPLATE(system_init),   &
    TEMPLATE(system_start),  &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(system_extend), &
#endif
    TEMPLATE(system_update), &
    TEMPLATE(system_get),    &
    TEMPLATE(system_copy),   &
    TEMPLATE(system_end)

  public ::                         &
    TEMPLATE(SYSTEM_ITER_OK),       &
    TEMPLATE(SYSTEM_ITER_END),      & 
    TEMPLATE(system_iterator_t),    &
    TEMPLATE(system_iterator_init), &
    TEMPLATE(system_iterator_next), &
    TEMPLATE(system_iterator_copy), &
    TEMPLATE(system_iterator_end)

  public ::                              &
    TEMPLATE(system_interpolation_init), &
    TEMPLATE(system_interpolation_eval), &
    TEMPLATE(system_interpolation_copy), &
    TEMPLATE(system_interpolation_end)

  type, public :: TEMPLATE(system_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(geometry_t),    pointer :: geo    =>null()
    type(space_t),       pointer :: space  =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(geo_t)                  :: gmt
    type(states_t)               :: st
  end type TEMPLATE(system_t)

  type, public :: TEMPLATE(system_interpolation_t)
    private
    type(TEMPLATE(system_t)), pointer :: self =>null()
    type(states_interpolation_t)      :: intrp
  end type TEMPLATE(system_interpolation_t)

  interface TEMPLATE(system_get)
    module procedure TEMPLATE(system_get_config)
    module procedure TEMPLATE(system_get_space)
    module procedure TEMPLATE(system_get_geometry)
    module procedure TEMPLATE(system_get_states)
    module procedure TEMPLATE(system_get_density)
  end interface TEMPLATE(system_get)

  interface TEMPLATE(system_update)
#ifdef SUBTEMPLATE_NAME
    module procedure TEMPLATE(system_update_build)
#endif
    module procedure TEMPLATE(system_update_finish)
  end interface TEMPLATE(system_update)

  interface TEMPLATE(system_interpolation_eval)
    module procedure TEMPLATE(system_interpolation_eval_1d)
    module procedure TEMPLATE(system_interpolation_eval_2d)
  end interface TEMPLATE(system_interpolation_eval)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_init_common)(this, space, config)
    type(TEMPLATE(system_t)),    intent(inout) :: this
    type(space_t),       target, intent(in)    :: space
    type(json_object_t), target, intent(in)    :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    this%config=>config
    call geo_get(this%gmt, this%geo)
    ASSERT(associated(this%geo))
    this%space=>space
    nullify(cnfg)
    this%sim=>null()
    call json_get(this%config, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call states_init(this%st, cnfg)
    nullify(cnfg)
    return
  end subroutine TEMPLATE(system_init_common)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_init)(this, space, config)
    type(TEMPLATE(system_t)), intent(out) :: this
    type(space_t),            intent(in)  :: space
    type(json_object_t),      intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(config, "geometry", cnfg, ierr)
    if(ierr==JSON_OK)then
      call geo_init(this%gmt, space, cnfg)
    else
      call geo_init(this%gmt, space)
    end if
    nullify(cnfg)
    call TEMPLATE(system_init_common)(this, space, config)
    return
  end subroutine TEMPLATE(system_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_extend)(this, that)
    type(TEMPLATE(system_t)), intent(inout) :: this
    type(sub_t),    optional, intent(in)    :: that
    !
    ASSERT(.not.associated(this%sim))
    call geo_extend(this%gmt, that)
    call states_extend(this%st, that)
    return
  end subroutine TEMPLATE(system_extend)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_update_build)(this, that)
    type(TEMPLATE(system_t)), intent(inout) :: this
    type(sub_t),              intent(in)    :: that
    !
    ASSERT(associated(this%sim))
    print *, "***: ", "TEMPLATE_NAME" , ": states_update_build"
    call states_update(this%st, that)
    return
  end subroutine TEMPLATE(system_update_build)
#else
  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_init)(this, space, config)
    type(TEMPLATE(system_t)), target, intent(out) :: this
    type(space_t),                    intent(in)  :: space
    type(json_object_t),      target, intent(in)  :: config
    !
    type(json_object_t), pointer :: cnfg
    type(space_t)                :: spc
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(config, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call space_init(spc, cnfg)
    nullify(cnfg)
    ASSERT(space==spc)
    call space_end(spc)
    call json_get(config, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call geo_init(this%gmt, this%space, cnfg)
    nullify(cnfg)
    call TEMPLATE(system_init_common)(this, space, config)
    return
  end subroutine TEMPLATE(system_init)
#endif

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_start)(this, sim)
    type(TEMPLATE(system_t)),   intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    print *, "***: states_start"
    call states_start(this%st, sim)
    return
  end subroutine TEMPLATE(system_start)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_update_finish)(this)
    type(TEMPLATE(system_t)), intent(inout) :: this
    !
    ASSERT(associated(this%sim))
    print *, "***: ", "TEMPLATE_NAME" , ": states_update_finish"
    call states_update(this%st)
    return
  end subroutine TEMPLATE(system_update_finish)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_get_config)(this, that)
    type(TEMPLATE(system_t)), target, intent(in) :: this
    type(json_object_t),     pointer             :: that
    !
    nullify(that)
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(system_get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_get_space)(this, that)
    type(TEMPLATE(system_t)), target, intent(in) :: this
    type(space_t),           pointer             :: that
    !
    that=>this%space
    return
  end subroutine TEMPLATE(system_get_space)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_get_geometry)(this, that)
    type(TEMPLATE(system_t)), target, intent(in) :: this
    type(geometry_t),        pointer             :: that
    !
    that=>this%geo
    return
  end subroutine TEMPLATE(system_get_geometry)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_get_states)(this, that)
    type(TEMPLATE(system_t)),  target, intent(in) :: this
    type(states_t),           pointer             :: that
    !
    that=>this%st
    return
  end subroutine TEMPLATE(system_get_states)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_get_density)(this, that)
    type(TEMPLATE(system_t)),  target, intent(in) :: this
    type(density_t),           pointer             :: that
    !
    call states_get(this%st, that)
    return
  end subroutine TEMPLATE(system_get_density)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_copy)(this, that)
    type(TEMPLATE(system_t)), intent(out) :: this
    type(TEMPLATE(system_t)), intent(in)  :: that
    !
    this%config=>that%config
    this%geo=>that%geo
    this%space=>that%space
    this%sim=>that%sim
    call geo_copy(this%gmt, that%gmt)
    call states_copy(this%st, that%st)
    return
  end subroutine TEMPLATE(system_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_end)(this)
    type(TEMPLATE(system_t)), intent(inout) :: this
    !
    call states_end(this%st)
    call geo_end(this%gmt)
    nullify(this%sim, this%space, this%geo, this%config)
    return
  end subroutine TEMPLATE(system_end)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_iterator_init)(this, that)
    type(TEMPLATE(system_iterator_t)), intent(out) :: this
    type(TEMPLATE(system_t)),          intent(in)  :: that
    !
    call geo_iterator_init(this, that%gmt)
    return
  end subroutine TEMPLATE(system_iterator_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_interpolation_init)(this, that, what)
    type(TEMPLATE(system_interpolation_t)), intent(out) :: this
    type(TEMPLATE(system_t)),       target, intent(in)  :: that
    type(density_t),                        intent(in)  :: what
    !
    this%self=>that
    call states_interpolation_init(this%intrp, that%st, what)
    return
  end subroutine TEMPLATE(system_interpolation_init)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_interpolation_eval_1d)(this, x, val)
    type(TEMPLATE(system_interpolation_t)), intent(in)  :: this
    real(kind=wp),            dimension(:), intent(in)  :: x
    real(kind=wp),                          intent(out) :: val
    !
    call states_interpolation_eval(this%intrp, x, val)
    return
  end subroutine TEMPLATE(system_interpolation_eval_1d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_interpolation_eval_2d)(this, x, val)
    type(TEMPLATE(system_interpolation_t)), intent(in)  :: this
    real(kind=wp),            dimension(:), intent(in)  :: x
    real(kind=wp),            dimension(:), intent(out) :: val
    !
    call states_interpolation_eval(this%intrp, x, val)
    return
  end subroutine TEMPLATE(system_interpolation_eval_2d)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_interpolation_copy)(this, that)
    type(TEMPLATE(system_interpolation_t)), intent(out) :: this
    type(TEMPLATE(system_interpolation_t)), intent(in)  :: that
    !
    this%self=>that%self
    call states_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(system_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(system_interpolation_end)(this)
    type(TEMPLATE(system_interpolation_t)), intent(inout) :: this
    !
    call states_interpolation_end(this%intrp)
    this%self=>null()
    return
  end subroutine TEMPLATE(system_interpolation_end)

end module TEMPLATE(system_m)

!! Local Variables:
!! mode: f90
!! End:
