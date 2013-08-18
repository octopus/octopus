module lda_m

  use interface_xc_m, only: interface_xc_t, interface_xc_init, interface_xc_end, &
    interface_xc_set_parameter, interface_xc_get_kind, &
    interface_xc_lda_exc, interface_xc_lda_vxc, interface_xc_lda_exc_vxc
  use json_m,         only: json_object_t
  use kinds_m,        only: wp

  implicit none

  private
  public ::              &
    lda_init,            &
    lda_set_parameter,   &
    lda_get_kind,        &
    lda_get_exc,         &
    lda_get_exc_and_vxc, &
    lda_get_vxc,         &
    lda_end

  type, public :: lda_t
    private
    type(interface_xc_t)  :: funct
  end type lda_t

contains

  ! ---------------------------------------------------------
  subroutine lda_init(this, id, ndim, nspin)
    type(lda_t), intent(out) :: this
    integer,     intent(in)  :: id
    integer,     intent(in)  :: ndim
    integer,     intent(in)  :: nspin
    !
    call interface_xc_init(this%funct, id, ndim, nspin)
    return
  end subroutine lda_init

  ! ---------------------------------------------------------
  subroutine lda_set_parameter(this, config)
    type(lda_t),         intent(inout) :: this
    type(json_object_t), intent(in)    :: config
    !
    call interface_xc_set_parameter(this%funct, config)
    return
  end subroutine lda_set_parameter

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
    call interface_xc_lda_exc(this%funct, density, exc)
    return
  end subroutine lda_get_exc

  ! ---------------------------------------------------------
  subroutine lda_get_exc_and_vxc(this, density, exc, vxc)
    type(lda_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    call interface_xc_lda_exc_vxc(this%funct, density, exc, vxc)
    return
  end subroutine lda_get_exc_and_vxc

  ! ---------------------------------------------------------
  subroutine lda_get_vxc(this, density, vxc)
    type(lda_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    call interface_xc_lda_vxc(this%funct, density, vxc)
    return
  end subroutine lda_get_vxc

  ! ---------------------------------------------------------
  subroutine lda_end(this)
    type(lda_t), intent(inout) :: this
    !
    call interface_xc_end(this%funct)
    return
  end subroutine lda_end

end module lda_m

!! Local Variables:
!! mode: f90
!! End:
