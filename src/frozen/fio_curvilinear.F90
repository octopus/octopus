#include "global.h"

module fio_curvilinear_m

  use global_m
  use messages_m
  use profiling_m

  use curv_briggs_m, only: curv_briggs_init
  use curv_gygi_m,   only: curv_gygi_init
  use curv_modine_m, only: curv_modine_init
  use geometry_m,    only: geometry_t
  use json_m,        only: JSON_OK, json_object_t, json_get
  use kinds_m,       only: wp

  use curvilinear_m, only: &
    CURV_METHOD_UNIFORM,   &
    CURV_METHOD_GYGI,      &
    CURV_METHOD_BRIGGS,    &
    CURV_METHOD_MODINE

  use curvilinear_m, only:                    &
    fio_curvilinear_t    => curvilinear_t,    &
    fio_curvilinear_copy => curvilinear_copy, &
    fio_curvilinear_end  => curvilinear_end

  use base_geom_m, only:           &
    fio_geom_t   => base_geom_t,   &
    fio_geom_get => base_geom_get

  use fio_simul_box_m, only: &
    fio_simul_box_t

  implicit none

  private
  public ::               &
    fio_curvilinear_t,    &
    fio_curvilinear_init, &
    fio_curvilinear_copy, &
    fio_curvilinear_end

contains

  ! ---------------------------------------------------------
  subroutine fio_curvilinear_init(this, sb, geom, config)
    type(fio_curvilinear_t), intent(out) :: this
    type(fio_simul_box_t),   intent(in)  :: sb
    type(fio_geom_t),        intent(in)  :: geom
    type(json_object_t),     intent(in)  :: config
    !
    type(geometry_t),         pointer :: geo
    real(kind=wp), dimension(MAX_DIM) :: spcng
    integer                           :: ierr
    !
    PUSH_SUB(fio_curvilinear_init)
    nullify(geo)
    call json_get(config, "method", this%method, ierr)
    if(ierr/=JSON_OK)this%method=CURV_METHOD_UNIFORM
    call fio_geom_get(geom, geo)
    ASSERT(associated(geo))
    select case(this%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_init(this%gygi, sb, geo)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_init(this%briggs, sb)
    case(CURV_METHOD_MODINE)
      call json_get(config, "spacing", spcng, ierr)
      ASSERT(ierr==JSON_OK)
      if(sb%dim<MAX_DIM)spcng(sb%dim+1:)=0.0_wp
      call curv_modine_init(this%modine, sb, geo, spcng)
    end select
    nullify(geo)
    POP_SUB(fio_curvilinear_init)
    return
  end subroutine fio_curvilinear_init

end module fio_curvilinear_m

!! Local Variables:
!! mode: f90
!! End:
