#include "global.h"

module fio_geometry_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,    only: json_object_t
  use space_m,   only: space_t
  use species_m, only: SPEC_FROZEN, species_set_type

  use igeometry_m, only: &
    geometry_t,          &
    geometry_init

  use igeometry_m, only:                &
    fio_geometry_copy => geometry_copy, &
    fio_geometry_end  => geometry_end

  implicit none

  private
  public ::     &
    geometry_t

  public ::            &
    fio_geometry_init, &
    fio_geometry_copy, &
    fio_geometry_end

contains

  ! ---------------------------------------------------------
  subroutine fio_geometry_init(this, space, json)
    type(geometry_t),      intent(out) :: this
    type(space_t), target, intent(in)  :: space
    type(json_object_t),   intent(in)  :: json
    !
    integer :: i
    !
    PUSH_SUB(fio_geometry_init)
    call geometry_init(this, space, json)
    do i= 1, this%nspecies
      call species_set_type(this%species(i), SPEC_FROZEN)
    end do
    POP_SUB(fio_geometry_init)
    return
  end subroutine fio_geometry_init

end module fio_geometry_m

!! Local Variables:
!! mode: f90
!! End:
