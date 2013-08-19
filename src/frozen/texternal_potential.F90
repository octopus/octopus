#include "global.h"
#include "template.h"

module TEMPLATE(external_potential_m)

  use global_m
  use messages_m
  use profiling_m

  use json_m,          only: json_object_t
  use kinds_m,         only: wp
  use species_m,       only: species_zval

  use TEMPLATE(geometry_m), only:       &
    geometry_t !=> TEMPLATE(geometry_t)

  use TEMPLATE(system_m), only:         &
    system_t   => TEMPLATE(system_t),   &
    system_get => TEMPLATE(system_get)

  use TEMPLATE(potential_m), only:                                &
    potential_get_potential => TEMPLATE(potential_get_potential)

  use TEMPLATE(potential_m), only:                                             &
    TEMPLATE(external_potential_t)          => TEMPLATE(potential_t),          &
    TEMPLATE(external_potential_init)       => TEMPLATE(potential_init),       &
    TEMPLATE(external_potential_start)      => TEMPLATE(potential_start),      &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(external_potential_extend)     => TEMPLATE(potential_extend),     &
#endif
    TEMPLATE(external_potential_update)     => TEMPLATE(potential_update),     &
    TEMPLATE(external_potential_set_energy) => TEMPLATE(potential_set_energy), &
    TEMPLATE(external_potential_get)        => TEMPLATE(potential_get),        &
    TEMPLATE(external_potential_get_size)   => TEMPLATE(potential_get_size),   &
    TEMPLATE(external_potential_get_energy) => TEMPLATE(potential_get_energy), &
    TEMPLATE(external_potential_copy)       => TEMPLATE(potential_copy),       &
    TEMPLATE(external_potential_end)        => TEMPLATE(potential_end)

  use TEMPLATE(potential_m), only:                                          &
    potential_interpolation_t    => TEMPLATE(potential_interpolation_t),    &
    potential_interpolation_init => TEMPLATE(potential_interpolation_init), &
    potential_interpolation_eval => TEMPLATE(potential_interpolation_eval), &
    potential_interpolation_copy => TEMPLATE(potential_interpolation_copy), &
    potential_interpolation_end  => TEMPLATE(potential_interpolation_end)

  use TEMPLATE(potential_m), only:                                         &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_OK) => TEMPLATE(POTENTIAL_INTRP_OK), &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_OD) => TEMPLATE(POTENTIAL_INTRP_OD), &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_NI) => TEMPLATE(POTENTIAL_INTRP_NI)

  implicit none

  private
  public ::                                     &
    TEMPLATE(external_potential_t),             &
    TEMPLATE(external_potential_init),          &
    TEMPLATE(external_potential_start),         &
#ifdef SUBTEMPLATE_NAME
    TEMPLATE(external_potential_extend),        &
#endif
    TEMPLATE(external_potential_update),        &
    TEMPLATE(external_potential_set_energy),    &
    TEMPLATE(external_potential_get),           &
    TEMPLATE(external_potential_get_size),      &
    TEMPLATE(external_potential_get_energy),    &
    TEMPLATE(external_potential_get_potential), &
    TEMPLATE(external_potential_copy),          &
    TEMPLATE(external_potential_end)

  public ::                                          &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_OK),           &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_OD),           &
    TEMPLATE(EXTERNAL_POTENTIAL_INTRP_NI),           &
    TEMPLATE(external_potential_interpolation_init), &
    TEMPLATE(external_potential_interpolation_eval), &
    TEMPLATE(external_potential_interpolation_copy), &
    TEMPLATE(external_potential_interpolation_end)

  type, public :: TEMPLATE(external_potential_interpolation_t)
    private
    type(TEMPLATE(external_potential_t)), pointer :: self =>null()
    type(potential_interpolation_t)               :: intrp
  end type TEMPLATE(external_potential_interpolation_t)

contains

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_get_potential)(this, that)
    type(TEMPLATE(external_potential_t)),  intent(in) :: this
    real(kind=wp),          dimension(:), pointer     :: that
    !
    call potential_get_potential(this, that)
    return
  end subroutine TEMPLATE(external_potential_get_potential)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_interpolation_init)(this, that)
    type(TEMPLATE(external_potential_interpolation_t)), intent(out) :: this
    type(TEMPLATE(external_potential_t)),       target, intent(in)  :: that
    !
    this%self=>that
    call potential_interpolation_init(this%intrp, that)
    return
  end subroutine TEMPLATE(external_potential_interpolation_init)

  ! ---------------------------------------------------------
  pure function TEMPLATE(external_potential_calc)(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c
    !
    real(kind=wp) :: v
    !
    real(kind=wp) :: r
    !
    r=sqrt(sum((x-y)**2))
    if(r<r_small) r=r_small
    v=-c/r
    return
  end function TEMPLATE(external_potential_calc)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_classical)(this, x, v)
    type(TEMPLATE(external_potential_t)), intent(in)  :: this
    real(kind=wp),          dimension(:), intent(in)  :: x
    real(kind=wp),                        intent(out) :: v
    !
    type(geometry_t), pointer :: geo
    integer                   :: i
    !
    v=0.0_wp
    geo=>null()
    call TEMPLATE(external_potential_get)(this, geo)
    do i = 1, geo%natoms
      v=v+TEMPLATE(external_potential_calc)(x, geo%atom(i)%x, species_zval(geo%atom(i)%spec))
    end do
    do i = 1, geo%ncatoms
      v=v+TEMPLATE(external_potential_calc)(x, geo%catom(i)%x, geo%catom(i)%charge)
    end do
    return
  end subroutine TEMPLATE(external_potential_classical)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_interpolation_eval)(this, x, v)
    type(TEMPLATE(external_potential_interpolation_t)), intent(in)  :: this
    real(kind=wp),                        dimension(:), intent(in)  :: x
    real(kind=wp),                                      intent(out) :: v
    !
    integer :: type, ierr
    !
    call potential_interpolation_eval(this%intrp, x, v, ierr)
    if(ierr/=TEMPLATE(EXTERNAL_POTENTIAL_INTRP_OK))&
      call TEMPLATE(external_potential_classical)(this%self, x, v)
    return
  end subroutine TEMPLATE(external_potential_interpolation_eval)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_interpolation_copy)(this, that)
    type(TEMPLATE(external_potential_interpolation_t)), intent(out) :: this
    type(TEMPLATE(external_potential_interpolation_t)), intent(in)  :: that
    !
    this%self=>that%self
    call potential_interpolation_copy(this%intrp, that%intrp)
    return
  end subroutine TEMPLATE(external_potential_interpolation_copy)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(external_potential_interpolation_end)(this)
    type(TEMPLATE(external_potential_interpolation_t)), intent(inout) :: this
    !
    call potential_interpolation_end(this%intrp)
    nullify(this%self)
    return
  end subroutine TEMPLATE(external_potential_interpolation_end)

end module TEMPLATE(external_potential_m)

!! Local Variables:
!! mode: f90
!! End:
