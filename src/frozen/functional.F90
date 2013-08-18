#include "global.h"

module functional_m

  use datasets_m, only: datasets_check
  use derivatives_m, only: derivatives_t
  use global_m
  use gga_m, only: gga_t, gga_init, gga_end, gga_set_parameter, gga_get_kind, &
    gga_get_exc, gga_get_exc_and_vxc, gga_get_vxc
  use interface_xc_m, only: XC_UNKNOWN, XC_NONE, INT_SOFT_COULOMB
  use json_m, only: json_object_t, json_init, json_set
  use kinds_m, only: wp
  use lda_m, only: lda_t, lda_init, lda_end, lda_set_parameter, lda_get_kind, &
    lda_get_exc, lda_get_exc_and_vxc, lda_get_vxc
  use mesh_function_m, only: dmf_integrate
  use messages_m
  use multigrid_m, only: INJECTION, dmultigrid_fine2coarse
  use parser_m, only: parse_integer, parse_float
  use profiling_m
  use XC_F90(lib_m), only: XC_EXCHANGE, XC_CORRELATION, XC_EXCHANGE_CORRELATION, XC_KINETIC, &
    XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_MGGA, XC_FAMILY_LCA, XC_FAMILY_OEP, XC_FAMILY_HYB_GGA, &
    XC_F90(family_from_id)

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
    functional_parse_config,             &
    functional_init,                     &
    functional_set_parameter,            &
    functional_get_kind,                 &
    functional_get_energy,               &
    functional_get_energy_and_potential, &
    functional_get_potential,            &
    functional_copy,                     &
    functional_end

  type, public :: functional_t
    private
    integer                      :: family = XC_NONE ! LDA, GGA, etc.
    integer                      :: id     = XC_NONE ! identifier
    integer                      :: nspin  = 0
    type(derivatives_t), pointer :: der    =>null()
    type(lda_t),         pointer :: lda    =>null()
    type(gga_t),         pointer :: gga    =>null()
  end type functional_t

  interface functional_get_energy
    module procedure functional_get_energy_from_exc
    module procedure functional_get_energy_from_density
  end interface functional_get_energy

contains

  ! ---------------------------------------------------------
  subroutine functional_init(this, der, id, nspin)
    type(functional_t),          intent(out) :: this
    type(derivatives_t), target, intent(in)  :: der
    integer,                     intent(in)  :: id
    integer,                     intent(in)  :: nspin
    !
    integer :: family
    !
    this%family=XC_NONE
    ASSERT(id>XC_UNKNOWN)
    if(id>XC_NONE)then
      this%der=>der
      family=XC_F90(family_from_id)(id)
      print *, "***: functional_init: ", id, family
      ASSERT(family>XC_UNKNOWN)
      if(family>XC_NONE)then
        this%id=id
        this%family=family
        select case(this%family)
        case(XC_FAMILY_LDA)
          SAFE_ALLOCATE(this%lda)
          call lda_init(this%lda, id, this%der%mesh%sb%dim, nspin)
        case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          SAFE_ALLOCATE(this%gga)
          call gga_init(this%gga, this%der, id, nspin)
        !case(XC_FAMILY_MGGA)
        !case(XC_FAMILY_LCA)
        !case(XC_FAMILY_OEP)
        case default
          ASSERT(.false.)
        end select
      end if
    end if
    return
  end subroutine functional_init

  ! ---------------------------------------------------------
  subroutine functional_parse_config(this, nel)
    type(json_object_t),     intent(out) :: this
    real(kind=wp), optional, intent(in)  :: nel
    !
    real(kind=wp) :: rtmp
    integer       :: itmp
    !
    call json_init(this)
    if(present(nel))call json_set(this, "nel", nel)
    call parse_float(datasets_check('Xalpha'), 1.0_wp, rtmp)
    call json_set(this, "Xalpha", rtmp)
    call parse_integer(datasets_check('Interaction1D'), INT_SOFT_COULOMB, itmp)
    call json_set(this, "Interaction1D", itmp)
    call parse_float(datasets_check("Interaction1DScreening"), 1.0_wp, rtmp)
    call json_set(this, "Interaction1DScreening", rtmp)
    call parse_integer(datasets_check("LB94_modified"), 0, itmp)
    call json_set(this, "LB94_modified", itmp)
    call parse_float(datasets_check('LB94_threshold'), 1.0e-6_wp, rtmp)
    call json_set(this, "LB94_threshold", rtmp)
    return
  end subroutine functional_parse_config
  

  ! ---------------------------------------------------------
  subroutine functional_set_parameter(this, config)
    type(functional_t),  intent(inout) :: this
    type(json_object_t), intent(in)    :: config
    !
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_set_parameter(this%lda, config)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_set_parameter(this%gga, config)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    end select
    return
  end subroutine functional_set_parameter

  ! ---------------------------------------------------------
  elemental function functional_get_kind(this) result(kind)
    type(functional_t), intent(in) :: this
    !
    integer :: kind
    !
    kind=XC_UNKNOWN
    select case(this%family)
    case(XC_FAMILY_LDA)
      kind=lda_get_kind(this%lda)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      kind=gga_get_kind(this%gga)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
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
    enrg=sum(density, dim=2)*exc
    energy=dmf_integrate(this%der%mesh, enrg)
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
    logical                                       :: calc
    !
    calc=.true.
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc(this%lda, density, exc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc(this%gga, density, exc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      calc=.false.
      energy=0.0_wp
    end select
    if(calc)&
      energy=functional_get_energy_from_exc(this, density, exc)
    return
  end function functional_get_energy_from_density

  ! ---------------------------------------------------------
  subroutine functional_get_exc(this, density, exc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    !
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc(this%lda, density, exc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc(this%gga, density, exc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      exc=0.0_wp
    end select
    return
  end subroutine functional_get_exc

  ! ---------------------------------------------------------
  subroutine functional_get_exc_and_vxc(this, density, exc, vxc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc_and_vxc(this%lda, density, exc, vxc)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc_and_vxc(this%gga, density, exc, vxc)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      exc=0.0_wp
      vxc=0.0_wp
    end select
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
    logical                                       :: calc
    !
    calc=.true.
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_exc_and_vxc(this%lda, density, exc, potential)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_exc_and_vxc(this%gga, density, exc, potential)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      calc=.false.
      energy=0.0_wp
      potential=0.0_wp
    end select
    if(calc)&
      energy=functional_get_energy_from_exc(this, density, exc)
    return
  end subroutine functional_get_energy_and_potential

  ! ---------------------------------------------------------
  subroutine functional_get_potential(this, density, potential)
    type(functional_t),                    intent(in)  :: this
    real(kind=wp), dimension(:,:),         intent(in)  :: density
    real(kind=wp), dimension(:,:), target, intent(out) :: potential
    !
    select case(this%family)
    case(XC_FAMILY_LDA)
      call lda_get_vxc(this%lda, density, potential)
    case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
      call gga_get_vxc(this%gga, density, potential)
    !case(XC_FAMILY_MGGA)
    !case(XC_FAMILY_LCA)
    !case(XC_FAMILY_OEP)
    case default
      potential=0.0_wp
    end select
    return
  end subroutine functional_get_potential

  ! ---------------------------------------------------------
  subroutine functional_copy(this, that)
    type(functional_t), intent(out) :: this
    type(functional_t), intent(in)  :: that
    !
    ASSERT(.false.)
!!$    integer                      :: family = XC_NONE ! LDA, GGA, etc.
!!$    integer                      :: id     = XC_NONE ! identifier
!!$    integer                      :: nspin  = 0
!!$    type(derivatives_t), pointer :: der    =>null()
!!$    type(lda_t),         pointer :: lda    =>null()
!!$    type(gga_t),         pointer :: gga    =>null()
    return
  end subroutine functional_copy

  ! ---------------------------------------------------------
  subroutine functional_end(this)
    type(functional_t), intent(inout) :: this
    !
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
    nullify(this%gga, this%lda, this%der)
    this%nspin=0
    this%id=XC_NONE
    this%family=XC_NONE
    return
  end subroutine functional_end

end module functional_m

!! Local Variables:
!! mode: f90
!! End:
