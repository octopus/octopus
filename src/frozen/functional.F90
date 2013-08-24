#include "global.h"

module functional_m

  use global_m
  use messages_m
  use profiling_m

  use grid_m,          only: grid_t
  use json_m,          only: JSON_OK, json_object_t, json_get
  use kinds_m,         only: wp
  use mesh_m,          only: mesh_t
  use mesh_function_m, only: dmf_integrate
  use simulation_m,    only: simulation_t, simulation_get

  use interface_xc_m, only: &
    XC_UNKNOWN,             &   
    XC_NONE

  use interface_xc_m, only:  &
    XC_EXCHANGE,             &
    XC_CORRELATION,          &
    XC_EXCHANGE_CORRELATION, &
    XC_KINETIC              

  use interface_xc_m, only:  &
    XC_FAMILY_LDA,           &
    XC_FAMILY_GGA,           &
    XC_FAMILY_MGGA,          &
    XC_FAMILY_LCA,           &
    XC_FAMILY_OEP,           &
    XC_FAMILY_HYB_GGA

  use gga_m, only:       &
    gga_t,               &
    gga_init,            &
    gga_start,           &
    gga_get_kind,        &
    gga_get_exc,         &
    gga_get_exc_and_vxc, &
    gga_get_vxc,         &
    gga_copy,            &
    gga_end

  use lda_m, only:       &
    lda_t,               &
    lda_init,            &
    lda_start,           &
    lda_get_kind,        &
    lda_get_exc,         &
    lda_get_exc_and_vxc, &
    lda_get_vxc,         &
    lda_copy,            &
    lda_end

  implicit none

  private
  public ::     &
    XC_UNKNOWN, &
    XC_NONE

  public ::                  &
    XC_EXCHANGE,             &
    XC_CORRELATION,          &
    XC_EXCHANGE_CORRELATION, &
    XC_KINETIC              

  public ::            &
    XC_FAMILY_LDA,     &
    XC_FAMILY_GGA,     &
    XC_FAMILY_MGGA,    &
    XC_FAMILY_LCA,     &
    XC_FAMILY_OEP,     &
    XC_FAMILY_HYB_GGA

  public ::                              &
    functional_init,                     &
    functional_start,                    &
    functional_get_kind,                 &
    functional_get_energy,               &
    functional_get_energy_and_potential, &
    functional_get_potential,            &
    functional_copy,                     &
    functional_end

  type, public :: functional_t
    private
    type(json_object_t), pointer :: config =>null()
    type(simulation_t),  pointer :: sim    =>null()
    type(mesh_t),        pointer :: mesh   =>null()
    type(lda_t),         pointer :: lda    =>null()
    type(gga_t),         pointer :: gga    =>null()
    integer                      :: family = XC_NONE ! LDA, GGA, etc.
  end type functional_t

  interface functional_get_energy
    module procedure functional_get_energy_from_exc
    module procedure functional_get_energy_from_density
  end interface functional_get_energy

contains

  ! ---------------------------------------------------------
  subroutine functional_init(this, config)
    type(functional_t),          intent(out) :: this
    type(json_object_t), target, intent(in)  :: config
    !
    integer :: ierr
    !
    PUSH_SUB(functional_init)
    this%config=>config
    nullify(this%sim, this%mesh)
    call json_get(this%config, "family", this%family, ierr)
    if(ierr/=JSON_OK)this%family=XC_NONE
    print *, "***: functional_init: ", this%family
    ASSERT(this%family>XC_UNKNOWN)
    if(this%family>XC_NONE)then
      select case(this%family)
      case(XC_FAMILY_LDA)
        SAFE_ALLOCATE(this%lda)
        call lda_init(this%lda, config)
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        SAFE_ALLOCATE(this%gga)
        call gga_init(this%gga, config)
      !case(XC_FAMILY_MGGA)
      !case(XC_FAMILY_LCA)
      !case(XC_FAMILY_OEP)
      case default
        call functional_end(this)
        ASSERT(.false.)
      end select
    end if
    POP_SUB(functional_init)
    return
  end subroutine functional_init

  ! ---------------------------------------------------------
  subroutine functional_start(this, sim)
    type(functional_t),         intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    !
    type(grid_t), pointer :: grid
    !
    PUSH_SUB(functional_start)
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    nullify(grid)
    call simulation_get(this%sim, grid)
    ASSERT(associated(grid))
    this%mesh=>grid%fine%mesh
    if(this%family>XC_NONE)then
      select case(this%family)
      case(XC_FAMILY_LDA)
        call lda_start(this%lda, sim)
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call gga_start(this%gga, sim)
      !case(XC_FAMILY_MGGA)
      !case(XC_FAMILY_LCA)
      !case(XC_FAMILY_OEP)
      case default
        call functional_end(this)
        ASSERT(.false.)
      end select
    end if
    POP_SUB(functional_start)
    return
  end subroutine functional_start

  ! ---------------------------------------------------------
  elemental function functional_get_kind(this) result(kind)
    type(functional_t), intent(in) :: this
    !
    integer :: kind
    !
    select case(this%family)
    case(XC_FAMILY_LDA)
      kind=lda_get_kind(this%lda)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      kind=gga_get_kind(this%gga)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      kind=XC_UNKNOWN
    end select
    return
  end function functional_get_kind

  ! ---------------------------------------------------------
  function functional_get_energy_from_exc(this, density, exc) result(energy)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(in)  :: exc
    !
    real(kind=wp) :: energy
    !
    real(kind=wp), dimension(size(exc)) :: enrg
    !
    PUSH_SUB(functional_get_energy_from_exc)
    enrg=sum(density, dim=2)*exc
    energy=dmf_integrate(this%mesh, enrg)
    POP_SUB(functional_get_energy_from_exc)
    return
  end function functional_get_energy_from_exc

  ! ---------------------------------------------------------
  function functional_get_energy_from_density(this, density) result(energy)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    !
    real(kind=wp) :: energy
    !
    real(kind=wp), dimension(size(density,dim=1)) :: exc
    !
    PUSH_SUB(functional_get_energy_from_density)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc(this%lda, density, exc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc(this%gga, density, exc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      ASSERT(.false.)
    end select
    energy=functional_get_energy_from_exc(this, density, exc)
    POP_SUB(functional_get_energy_from_density)
    return
  end function functional_get_energy_from_density

  ! ---------------------------------------------------------
  subroutine functional_get_exc(this, density, exc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    !
    PUSH_SUB(functional_get_exc)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc(this%lda, density, exc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc(this%gga, density, exc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(functional_get_exc)
    return
  end subroutine functional_get_exc

  ! ---------------------------------------------------------
  subroutine functional_get_exc_and_vxc(this, density, exc, vxc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    PUSH_SUB(functional_get_exc_and_vxc)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc_and_vxc(this%lda, density, exc, vxc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc_and_vxc(this%gga, density, exc, vxc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(functional_get_exc_and_vxc)
    return
  end subroutine functional_get_exc_and_vxc

  ! ---------------------------------------------------------
  subroutine functional_get_energy_and_potential(this, density, energy, potential)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp),                 intent(out) :: energy
    real(kind=wp), dimension(:,:), intent(out) :: potential
    !
    real(kind=wp), dimension(size(density,dim=1)) :: exc
    !
    PUSH_SUB(functional_get_energy_and_potential)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc_and_vxc(this%lda, density, exc, potential)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc_and_vxc(this%gga, density, exc, potential)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      ASSERT(.false.)
    end select
    energy=functional_get_energy_from_exc(this, density, exc)
    POP_SUB(functional_get_energy_and_potential)
    return
  end subroutine functional_get_energy_and_potential

  ! ---------------------------------------------------------
  subroutine functional_get_potential(this, density, potential)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: potential
    !
    PUSH_SUB(functional_get_potential)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_vxc(this%lda, density, potential)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_vxc(this%gga, density, potential)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      ASSERT(.false.)
    end select
    POP_SUB(functional_get_potential)
    return
  end subroutine functional_get_potential

  ! ---------------------------------------------------------
  subroutine functional_copy(this, that)
    type(functional_t), intent(out) :: this
    type(functional_t), intent(in)  :: that
    !
    PUSH_SUB(functional_copy)
    this%config=>that%config
    this%sim=>that%sim
    this%mesh=>that%mesh
    this%family=that%family
    select case(this%family)
    case(XC_FAMILY_LDA)
      SAFE_ALLOCATE(this%lda)
      call lda_copy(this%lda, that%lda)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      SAFE_ALLOCATE(this%gga)
      call gga_copy(this%gga, that%gga)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      nullify(this%lda, this%lda)
      ASSERT(this%family==XC_NONE)
    end select
    POP_SUB(functional_copy)
    return
  end subroutine functional_copy

  ! ---------------------------------------------------------
  subroutine functional_end(this)
    type(functional_t), intent(inout) :: this
    !
    PUSH_SUB(functional_end)
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_end(this%lda)
      SAFE_DEALLOCATE_P(this%lda)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_end(this%gga)
      SAFE_DEALLOCATE_P(this%gga)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    end select
    this%family=XC_NONE
    nullify(this%gga, this%lda, this%mesh, this%sim, this%config)
    POP_SUB(functional_end)
    return
  end subroutine functional_end

end module functional_m

!! Local Variables:
!! mode: f90
!! End:
