#include "global.h"

module fio_curvilinear_oct_m

  use curv_briggs_oct_m
  use curv_gygi_oct_m
  use curv_modine_oct_m
  use curvilinear_oct_m
  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m

  implicit none

  private

  public ::               &
    fio_curvilinear_init, &
    fio_curvilinear_copy, &
    fio_curvilinear_end

contains

  ! ---------------------------------------------------------
  subroutine fio_curvilinear_init(this, sb, geo, config)
    type(curvilinear_t), intent(out) :: this
    type(simul_box_t),   intent(in)  :: sb
    type(geometry_t),    intent(in)  :: geo
    type(json_object_t), intent(in)  :: config

    real(kind=wp), dimension(MAX_DIM) :: spcng
    integer                           :: ierr

    PUSH_SUB(fio_curvilinear_init)

    call json_get(config, "method", this%method, ierr)
    if(ierr/=JSON_OK) this%method = CURV_METHOD_UNIFORM
    select case(this%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_init(this%gygi, sb, geo)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_init(this%briggs, sb)
    case(CURV_METHOD_MODINE)
      call json_get(config, "spacing", spcng(1:sb%dim), ierr)
      ASSERT(ierr==JSON_OK)
      if(size(spcng)>sb%dim) spcng(sb%dim+1:) = 0.0_wp
      call curv_modine_init(this%modine, sb, geo, spcng)
    end select

    POP_SUB(fio_curvilinear_init)
  end subroutine fio_curvilinear_init

  ! ---------------------------------------------------------
  subroutine fio_curvilinear_copy(this, that)
    type(curvilinear_t), intent(inout) :: this
    type(curvilinear_t), intent(in)    :: that

    PUSH_SUB(fio_curvilinear_copy)

    call curvilinear_copy(this, that)

    POP_SUB(fio_curvilinear_copy)
  end subroutine fio_curvilinear_copy

  ! ---------------------------------------------------------
  subroutine fio_curvilinear_end(this)
    type(curvilinear_t), intent(inout) :: this

    PUSH_SUB(fio_curvilinear_end)

    call curvilinear_end(this)

    POP_SUB(fio_curvilinear_end)
  end subroutine fio_curvilinear_end

end module fio_curvilinear_oct_m

!! Local Variables:
!! mode: f90
!! End:
