#include "global.h"
#include "template.h"

module TEMPLATE(term_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use TEMPLATE(density_m), only:      &
    density_t => TEMPLATE(density_t)

  use TEMPLATE(geometry_m), only:       &
    geometry_t !=> TEMPLATE(geometry_t)

  use TEMPLATE(states_m), only:     &
    states_t => TEMPLATE(states_t)

  use TEMPLATE(system_m), only:         &
    system_t   => TEMPLATE(system_t),   &
    system_get => TEMPLATE(system_get)

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
  public ::                    &
    TEMPLATE(term_init),       &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(term_extend),     &
#endif
    TEMPLATE(term_get),        &
    TEMPLATE(term_set_energy), &
    TEMPLATE(term_get_energy), &
    TEMPLATE(term_copy),       &
    TEMPLATE(term_end)

  type, public :: TEMPLATE(term_t)
    private
    type(json_object_t), pointer :: config =>null()
    type(system_t),      pointer :: sys    =>null()
    real(kind=wp)                :: energy = 0.0_wp
#ifdef SUBTEMPLATE_NAME
    logical                      :: temp   = .false.
    logical                      :: block  = .false.
    type(list_t)                 :: list
#endif
  end type TEMPLATE(term_t)

  interface TEMPLATE(term_get)
    module procedure TEMPLATE(term_get_config)
    module procedure TEMPLATE(term_get_system)
    module procedure TEMPLATE(term_get_geometry)
    module procedure TEMPLATE(term_get_states)
    module procedure TEMPLATE(term_get_density)
  end interface TEMPLATE(term_get)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_init)(this, sys, config)
    type(TEMPLATE(term_t)),      intent(out) :: this
    type(system_t),      target, intent(in)  :: sys
    type(json_object_t), target, intent(in)  :: config
    !
    this%config=>config
    this%sys=>sys
    this%energy=0.0_wp
    return
  end subroutine TEMPLATE(term_init)

#ifdef SUBTEMPLATE_NAME
  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_extend)(this, that)
    type(TEMPLATE(term_t)), intent(inout) :: this
    type(sub_t),  optional, intent(in)    :: that
    !
    type(json_object_t), pointer :: cnfg
    logical                      :: temp
    integer                      :: ierr
    !
    ASSERT(associated(this%config))
    if(present(that))then
      ASSERT(.not.this%block)
      nullify(cnfg)
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
  end subroutine TEMPLATE(term_extend)
#endif

  ! ---------------------------------------------------------
  elemental subroutine TEMPLATE(term_set_energy)(this, that)
    type(TEMPLATE(term_t)), intent(inout) :: this
    real(kind=wp),          intent(in)    :: that
    !
    this%energy=that
    return
  end subroutine TEMPLATE(term_set_energy)

  ! ---------------------------------------------------------
  elemental function TEMPLATE(term_get_energy)(this) result(that)
    type(TEMPLATE(term_t)), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=this%energy
    return
  end function TEMPLATE(term_get_energy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_get_config)(this, that)
    type(TEMPLATE(term_t)), target, intent(in) :: this
    type(json_object_t),   pointer             :: that
    !
    that=>null()
    if(associated(this%config))&
      that=>this%config
    return
  end subroutine TEMPLATE(term_get_config)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_get_system)(this, that)
    type(TEMPLATE(term_t)), target, intent(in) :: this
    type(system_t),        pointer             :: that
    !
    that=>null()
    if(associated(this%sys))&
      that=>this%sys
    return
  end subroutine TEMPLATE(term_get_system)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_get_geometry)(this, that)
    type(TEMPLATE(term_t)), intent(in) :: this
    type(geometry_t),      pointer     :: that
    !
    that=>null()
    if(associated(this%sys))&
      call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(term_get_geometry)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_get_states)(this, that)
    type(TEMPLATE(term_t)), intent(in) :: this
    type(states_t),        pointer     :: that
    !
    that=>null()
    if(associated(this%sys))&
      call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(term_get_states)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_get_density)(this, that)
    type(TEMPLATE(term_t)), intent(in) :: this
    type(density_t),       pointer     :: that
    !
    that=>null()
    if(associated(this%sys))&
      call system_get(this%sys, that)
    return
  end subroutine TEMPLATE(term_get_density)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_copy)(this, that)
    type(TEMPLATE(term_t)), intent(out) :: this
    type(TEMPLATE(term_t)), intent(in)  :: that
    !
    this%config=>that%config
    this%sys=>that%sys
    this%energy=that%energy
#ifdef SUBTEMPLATE_NAME
    this%block=that%block
    call list_copy(this%list, that%list)
#endif
    return
  end subroutine TEMPLATE(term_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(term_end)(this)
    type(TEMPLATE(term_t)), intent(inout) :: this
    !
#ifdef SUBTEMPLATE_NAME
    call list_end(this%list)
    this%block=.false.
#endif
    this%energy=0.0_wp
    nullify(this%sys, this%config)
    return
  end subroutine TEMPLATE(term_end)

end module TEMPLATE(term_m)

!! Local Variables:
!! mode: f90
!! End:
