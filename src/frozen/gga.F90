#include "global.h"

module gga_m

  use global_m
  use messages_m
  use profiling_m

  use derivatives_m, only: derivatives_t, dderivatives_grad, dderivatives_div
  use grid_m,        only: grid_t
  use json_m,        only: json_object_t
  use kinds_m,       only: wp
  use simulation_m,  only: simulation_t, simulation_get

  use interface_xc_m, only: &
    XC_NONE,                &
    XC_POLARIZED,           &
    XC_UNPOLARIZED

  use interface_xc_m, only:    &
    interface_xc_t,            &
    interface_xc_init,         &
    interface_xc_start,        &
    interface_xc_is_polarized, &
    interface_xc_get_nspin,    &
    interface_xc_get_kind,     &
    interface_xc_density,      &
    interface_xc_gga_exc,      &
    interface_xc_gga_vxc,      &
    interface_xc_gga_exc_vxc,  &
    interface_xc_copy,         &
    interface_xc_end

  implicit none

  private
  public ::              &
    gga_init,            &
    gga_start,           &
    gga_get_kind,        &
    gga_get_exc,         &
    gga_get_exc_and_vxc, &
    gga_get_vxc,         &
    gga_copy,            &
    gga_end

  type, public :: gga_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(derivatives_t), pointer :: der    =>null()
    integer                      :: spin   = XC_NONE ! XC_UNPOLARIZED | XC_POLARIZED
    integer                      :: nspin  = 0
    integer                      :: ndim   = 0
    type(interface_xc_t)         :: funct
  end type gga_t

contains

  ! ---------------------------------------------------------
  subroutine gga_init(this, config)
    type(gga_t),                 intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    PUSH_SUB(gga_init)
    this%config=>config
    nullify(this%sim, this%der)
    call interface_xc_init(this%funct, config)
    this%nspin=interface_xc_get_nspin(this%funct)
    this%spin=XC_UNPOLARIZED
    if(interface_xc_is_polarized(this%funct))&
      this%spin=XC_POLARIZED
    POP_SUB(gga_init)
    return
  end subroutine gga_init

  ! ---------------------------------------------------------
  subroutine gga_start(this, sim)
    type(gga_t),                intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    type(grid_t), pointer :: grid
    !
    PUSH_SUB(gga_start)
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    nullify(grid)
    call simulation_get(this%sim, grid)
    ASSERT(associated(grid))
    this%der=>grid%fine%der
    this%ndim=this%der%mesh%sb%dim
    call interface_xc_start(this%funct, sim)
    POP_SUB(gga_start)
    return
  end subroutine gga_start

  ! ---------------------------------------------------------
  elemental function gga_get_kind(this) result(kind)
    type(gga_t), intent(in) :: this
    !
    integer :: kind
    !
    kind=interface_xc_get_kind(this%funct)
    return
  end function gga_get_kind

  ! ---------------------------------------------------------
  subroutine gga_potential_correction(this, potential, dedg)
    type(gga_t),                     intent(in)    :: this
    real(kind=wp), dimension(:,:),   intent(inout) :: potential
    real(kind=wp), dimension(:,:,:), intent(inout) :: dedg !Should be intent in only
    !
    real(kind=wp), dimension(size(potential,dim=1)) :: div
    integer                                         :: i
    !
    PUSH_SUB(gga_potential_correction)
    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    do i = 1, this%spin
      call dderivatives_div(this%der, dedg(:,:,i), div)
      potential(:,i)=potential(:,i)-div
    end do
    POP_SUB(gga_potential_correction)
    return
  end subroutine gga_potential_correction

  ! ---------------------------------------------------------
  subroutine gga_get_exc(this, density, exc)
    type(gga_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    !
    real(kind=wp), dimension(size(density,dim=1),this%spin) :: dens
    real(kind=wp), dimension(size(exc),this%ndim,this%spin) :: grad
    integer                                                 :: isp
    !
    PUSH_SUB(gga_get_exc)
    call interface_xc_density(this%funct, dens, density)
    do isp = 1, this%spin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_gga_exc(this%funct, dens, grad, exc)
    POP_SUB(gga_get_exc)
    return
  end subroutine gga_get_exc

  ! ---------------------------------------------------------
  subroutine gga_get_exc_and_vxc(this, density, exc, vxc)
    type(gga_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    real(kind=wp), dimension(size(density,dim=1),          this%spin ) :: dens
    real(kind=wp), dimension(size(exc),          this%ndim,this%nspin) :: grad
    real(kind=wp), dimension(size(density,dim=1),this%ndim,this%spin ) :: dedg
    integer                                                            :: isp
    !
    PUSH_SUB(gga_get_exc_and_vxc)
    call interface_xc_density(this%funct, dens, density)
    do isp = 1, this%nspin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_gga_exc_vxc(this%funct, dens, grad, exc, vxc, dedg)
    dedg(this%der%mesh%np+1:,:,:)=0.0_wp
    call gga_potential_correction(this, vxc, dedg)
    POP_SUB(gga_get_exc_and_vxc)
    return
  end subroutine gga_get_exc_and_vxc

  ! ---------------------------------------------------------
  subroutine gga_get_vxc(this, density, potential)
    type(gga_t),                   intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: potential
    !
    real(kind=wp), dimension(this%der%mesh%np_part,this%spin)                      :: dens
    real(kind=wp), dimension(this%der%mesh%np,this%der%mesh%sb%dim,this%spin)      :: grad
    real(kind=wp), dimension(this%der%mesh%np_part,this%der%mesh%sb%dim,this%spin) :: dedg
    integer                                                                        :: isp
    !
    PUSH_SUB(gga_get_vxc)
    call interface_xc_density(this%funct, dens, density)
    do isp = 1, this%spin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_gga_vxc(this%funct, dens(1:this%der%mesh%np,:), grad, &
      potential, dedg(1:this%der%mesh%np,:,:))
    dedg(this%der%mesh%np+1:,:,:)=0.0_wp
    call gga_potential_correction(this, potential, dedg)
    POP_SUB(gga_get_vxc)
    return
  end subroutine gga_get_vxc

  ! ---------------------------------------------------------
  subroutine gga_copy(this, that)
    type(gga_t), intent(out) :: this
    type(gga_t), intent(in)  :: that
    !
    PUSH_SUB(gga_copy)
    this%config=>that%config
    this%sim=>that%sim
    this%der=>that%der
    this%spin=that%spin
    this%nspin=that%nspin
    this%ndim=that%ndim
    call interface_xc_copy(this%funct, that%funct)
    POP_SUB(gga_copy)
    return
  end subroutine gga_copy

  ! ---------------------------------------------------------
  subroutine gga_end(this)
    type(gga_t), intent(inout) :: this
    !
    PUSH_SUB(gga_end)
    call interface_xc_end(this%funct)
    this%ndim=0
    this%nspin=0
    this%spin=XC_NONE
    nullify(this%der, this%sim, this%config)
    POP_SUB(gga_end)
    return
  end subroutine gga_end

end module gga_m

!! Local Variables:
!! mode: f90
!! End:
