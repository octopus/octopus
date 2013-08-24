#include "global.h"

module lda_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,       only: json_object_t
  use kinds_m,      only: wp
  use simulation_m, only: simulation_t

  use interface_xc_m, only:   &
    interface_xc_t,           &
    interface_xc_init,        &
    interface_xc_start,       &
    interface_xc_get_kind,    &
    interface_xc_lda_exc,     &
    interface_xc_lda_vxc,     &
    interface_xc_lda_exc_vxc, &
    interface_xc_copy,        &
    interface_xc_end

  implicit none

  private
  public ::              &
    lda_init,            &
    lda_start,           &
    lda_get_kind,        &
    lda_get_exc,         &
    lda_get_exc_and_vxc, &
    lda_get_vxc,         &
    lda_copy,            &
    lda_end

  type, public :: lda_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(interface_xc_t)         :: funct
  end type lda_t

contains

  ! ---------------------------------------------------------
  subroutine lda_init(this, config)
    type(lda_t),                 intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(lda_init)
    this%config=>config
    this%sim=>null()
    call interface_xc_init(this%funct, config)
    POP_SUB(lda_init)
    return
  end subroutine lda_init

  ! ---------------------------------------------------------
  subroutine lda_start(this, sim)
    type(lda_t),                intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    PUSH_SUB(lda_start)
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    call interface_xc_start(this%funct, sim)
    POP_SUB(lda_start)
    return
  end subroutine lda_start

  ! ---------------------------------------------------------
  elemental function lda_get_kind(this) result(kind)
    type(lda_t), intent(in) :: this
    !
    integer :: kind
    !
    kind=interface_xc_get_kind(this%funct)
    return
  end function lda_get_kind

  ! ---------------------------------------------------------
  subroutine lda_get_exc(this, density, exc)
    type(lda_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    !
    PUSH_SUB(lda_get_exc)
    call interface_xc_lda_exc(this%funct, density, exc)
    POP_SUB(lda_get_exc)
    return
  end subroutine lda_get_exc

  ! ---------------------------------------------------------
  subroutine lda_get_exc_and_vxc(this, density, exc, vxc)
    type(lda_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    PUSH_SUB(lda_get_exc_and_vxc)
    call interface_xc_lda_exc_vxc(this%funct, density, exc, vxc)
    POP_SUB(lda_get_exc_and_vxc)
    return
  end subroutine lda_get_exc_and_vxc

  ! ---------------------------------------------------------
  subroutine lda_get_vxc(this, density, vxc)
    type(lda_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    PUSH_SUB(lda_get_vxc)
    call interface_xc_lda_vxc(this%funct, density, vxc)
    POP_SUB(lda_get_vxc)
    return
  end subroutine lda_get_vxc

  ! ---------------------------------------------------------
  subroutine lda_copy(this, that)
    type(lda_t), intent(out) :: this
    type(lda_t), intent(in)  :: that
    !
    PUSH_SUB(lda_copy)
    this%config=>that%config
    this%sim=>that%sim
    call interface_xc_copy(this%funct, that%funct)
    POP_SUB(lda_copy)
    return
  end subroutine lda_copy

  ! ---------------------------------------------------------
  subroutine lda_end(this)
    type(lda_t), intent(inout) :: this
    !
    PUSH_SUB(lda_end)
    call interface_xc_end(this%funct)
    nullify(this%sim, this%config)
    POP_SUB(lda_end)
    return
  end subroutine lda_end

end module lda_m

!! Local Variables:
!! mode: f90
!! End:
