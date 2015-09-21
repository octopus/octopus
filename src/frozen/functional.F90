#include "global.h"

module functional_m

  use derivatives_m
  use global_m
  use grid_m
  use interface_xc_m
  use json_m
  use kinds_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use simulation_m
  use XC_F90(lib_m)

  implicit none

  private

  public ::           &
    FUNCT_XC_NONE

  public ::       &
    functional_t

  public ::           &
    functional_init,  &
    functional_start, &
    functional_calc,  &
    functional_get,   &
    functional_copy,  &
    functional_end

  integer, parameter :: FUNCT_XC_NONE = XC_NONE

  type :: functional_t
    private
    type(simulation_t), pointer :: sim  =>null()
    type(mesh_t),       pointer :: mesh =>null()
    type(interface_xc_t)        :: funct
  end type functional_t

  interface functional_get
    module procedure functional_get_info
    module procedure functional_get_sim
    module procedure functional_get_mesh
  end interface functional_get

  interface functional_calc
    module procedure functional_calc_energy
    module procedure functional_calc_potential
    module procedure functional_calc_energy_and_potential
  end interface functional_calc

contains

  ! ---------------------------------------------------------
  subroutine functional_init(this, id, nspin)
    type(functional_t), intent(out) :: this
    integer,            intent(in)  :: id
    integer,            intent(in)  :: nspin

    PUSH_SUB(functional_init)

    nullify(this%sim,this%mesh)
    call interface_xc_init(this%funct, id, nspin)

    POP_SUB(functional_init)
  end subroutine functional_init

  ! ---------------------------------------------------------
  subroutine functional_start(this, sim, fine)
    type(functional_t),         intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    logical,          optional, intent(in)    :: fine

    type(grid_t),        pointer :: grid
    type(derivatives_t), pointer :: der
    logical                      :: fn

    PUSH_SUB(functional_start)

    ASSERT(.not.associated(this%sim))
    ASSERT(.not.associated(this%mesh))
    this%sim=>sim
    call simulation_get(this%sim, grid)
    ASSERT(associated(grid))
    fn=.false.
    if(grid%have_fine_mesh)then
      if(present(fine))fn=fine
    end if
    call simulation_get(this%sim, this%mesh, fine=fn)
    ASSERT(associated(this%mesh))
    call simulation_get(this%sim, der, fine=fn)
    ASSERT(associated(der))
    call interface_xc_start(this%funct, this%mesh, der)

    POP_SUB(functional_start)
  end subroutine functional_start

  ! ---------------------------------------------------------
  subroutine functional_get_info(this, id, family, kind, nspin, polarized, ndim)
    type(functional_t), intent(in)  :: this
    integer,  optional, intent(out) :: id
    integer,  optional, intent(out) :: family
    integer,  optional, intent(out) :: kind
    integer,  optional, intent(out) :: nspin
    logical,  optional, intent(out) :: polarized
    integer,  optional, intent(out) :: ndim

    PUSH_SUB(functional_get_info)

    call interface_xc_get(this%funct, id, family, kind, nspin, polarized, ndim)

    POP_SUB(functional_get_info)
  end subroutine functional_get_info

  ! ---------------------------------------------------------
  subroutine functional_get_sim(this, that)
    type(functional_t),  target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(functional_get_sim)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(functional_get_sim)
  end subroutine functional_get_sim

  ! ---------------------------------------------------------
  subroutine functional_get_mesh(this, that)
    type(functional_t), target, intent(in) :: this
    type(mesh_t),      pointer             :: that

    PUSH_SUB(functional_get_mesh)

    nullify(that)
    if(associated(this%mesh)) that => this%mesh

    POP_SUB(functional_get_mesh)
  end subroutine functional_get_mesh

  ! ---------------------------------------------------------
  subroutine functional_exc(this, density, exc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc

    integer :: family

    PUSH_SUB(functional_exc)

    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_exc(this%funct, density, exc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_exc(this%funct, density, exc)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(functional_exc)
  end subroutine functional_exc

  ! ---------------------------------------------------------
  subroutine functional_vxc(this, density, vxc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    integer :: family

    PUSH_SUB(functional_vxc)

    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_vxc(this%funct, density, vxc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_vxc(this%funct, density, vxc)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(functional_vxc)
  end subroutine functional_vxc

  ! ---------------------------------------------------------
  subroutine functional_exc_vxc(this, density, exc, vxc)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    integer :: family

    PUSH_SUB(functional_exc_vxc)

    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_exc_vxc(this%funct, density, exc, vxc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_exc_vxc(this%funct, density, exc, vxc)
    case default
      ASSERT(.false.)
    end select

    POP_SUB(functional_exc_vxc)
  end subroutine functional_exc_vxc

  ! ---------------------------------------------------------
  function functional_energy(this, density, exc) result(energy)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(in)  :: exc

    real(kind=wp) :: energy

    real(kind=wp), dimension(this%mesh%np) :: enrg
    integer                                :: ip

    PUSH_SUB(functional_energy)

    do ip = 1, this%mesh%np
      enrg(ip)=sum(density(ip,:))*exc(ip)
    end do
    energy=dmf_integrate(this%mesh, enrg)

    POP_SUB(functional_energy)
  end function functional_energy

  ! ---------------------------------------------------------
  subroutine functional_calc_energy(this, density, energy)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp),                 intent(out) :: energy

    real(kind=wp), dimension(this%mesh%np) :: exc

    PUSH_SUB(functional_calc_energy)

    call functional_exc(this, density, exc)
    energy=functional_energy(this, density, exc)

    POP_SUB(functional_calc_energy)
  end subroutine functional_calc_energy

  ! ---------------------------------------------------------
  subroutine functional_calc_potential(this, density, potential)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: potential

    PUSH_SUB(functional_calc_potential)

    call functional_vxc(this, density, potential)

    POP_SUB(functional_calc_potential)
  end subroutine functional_calc_potential

  ! ---------------------------------------------------------
  subroutine functional_calc_energy_and_potential(this, density, energy, potential)
    type(functional_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp),                 intent(out) :: energy
    real(kind=wp), dimension(:,:), intent(out) :: potential

    real(kind=wp), dimension(this%mesh%np) :: exc

    PUSH_SUB(functional_calc_energy_and_potential)

    call functional_exc_vxc(this, density, exc, potential)
    energy=functional_energy(this, density, exc)

    POP_SUB(functional_calc_energy_and_potential)
  end subroutine functional_calc_energy_and_potential

  ! ---------------------------------------------------------
  subroutine functional_copy(this, that)
    type(functional_t), intent(inout) :: this
    type(functional_t), intent(in)    :: that

    PUSH_SUB(functional_copy)

    call functional_end(this)
    this%sim => that%sim
    this%mesh => that%mesh
    call interface_xc_copy(this%funct, that%funct)

    POP_SUB(functional_copy)
  end subroutine functional_copy

  ! ---------------------------------------------------------
  subroutine functional_end(this)
    type(functional_t), intent(inout) :: this

    PUSH_SUB(functional_end)

    nullify(this%sim,this%mesh)
    call interface_xc_end(this%funct)

    POP_SUB(functional_end)
  end subroutine functional_end

end module functional_m

!! Local Variables:
!! mode: f90
!! End:
