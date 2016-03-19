#include "global.h"

module functional_oct_m

  use derivatives_oct_m
  use global_oct_m
  use interface_xc_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m
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
    functional_stop,  &
    functional_calc,  &
    functional_get,   &
    functional_copy,  &
    functional_end

  integer, parameter :: FUNCT_XC_NONE = XC_NONE

  type :: functional_t
    private
    type(json_object_t), pointer :: config =>null()
    type(storage_t),     pointer :: dnst   =>null()
    type(simulation_t),  pointer :: sim    =>null()
    logical                      :: fine   = .false.
    type(interface_xc_t)         :: funct
  end type functional_t

  interface functional_init
    module procedure functional_init_config
    module procedure functional_init_type
  end interface functional_init

  interface functional_get
    module procedure functional_get_info
    module procedure functional_get_config
    module procedure functional_get_sim
  end interface functional_get

  interface functional_calc
    module procedure functional_calc_energy
    module procedure functional_calc_potential
    module procedure functional_calc_energy_and_potential
  end interface functional_calc

contains

  ! ---------------------------------------------------------
  subroutine functional_init_config(this, id, polarized, blocksize)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: id
    logical,   optional, intent(in)  :: polarized
    integer,   optional, intent(in)  :: blocksize

    PUSH_SUB(functional_init_config)

    call json_init(this)
    if(present(id)) call json_set(this, "functional", id)
    if(present(polarized)) call json_set(this, "polarized", polarized)
    if(present(blocksize)) call json_set(this, "blocksize", blocksize)

    POP_SUB(functional_init_config)
  end subroutine functional_init_config

  ! ---------------------------------------------------------
  subroutine functional_init_type(this, density, config)
    type(functional_t),          intent(out) :: this
    type(storage_t),     target, intent(in)  :: density
    type(json_object_t), target, intent(in)  :: config

    integer :: id, nblck, ierr
    logical :: plrz

    PUSH_SUB(functional_init_type)

    this%config => config
    this%dnst => density
    call json_get(this%config, "functional", id, ierr)
    if(ierr/=JSON_OK) id = FUNCT_XC_NONE
    call json_get(this%config, "polarized", plrz, ierr)
    if(ierr/=JSON_OK) plrz = .true.
    call json_get(this%config, "blocksize", nblck, ierr)
    if(ierr/=JSON_OK) nblck = 0
    if(nblck>0)then
      call interface_xc_init(this%funct, id, plrz, nblck)
    else
      call interface_xc_init(this%funct, id, plrz)
    end if

    POP_SUB(functional_init_type)
  end subroutine functional_init_type

  ! ---------------------------------------------------------
  subroutine functional_start(this, sim)
    type(functional_t),         intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim

    type(mesh_t),        pointer :: mesh
    type(derivatives_t), pointer :: der

    PUSH_SUB(functional_start)

    ASSERT(associated(this%config))
    ASSERT(associated(this%dnst))
    ASSERT(.not.associated(this%sim))
    nullify(mesh, der)
    this%sim => sim
    call storage_get(this%dnst, fine=this%fine)
    call simulation_get(this%sim, mesh, fine=this%fine)
    ASSERT(associated(mesh))
    call simulation_get(this%sim, der, fine=this%fine)
    ASSERT(associated(der))
    call interface_xc_start(this%funct, mesh, der)
    nullify(mesh, der)

    POP_SUB(functional_start)
  end subroutine functional_start

  ! ---------------------------------------------------------
  subroutine functional_stop(this)
    type(functional_t), intent(inout) :: this

    PUSH_SUB(functional_stop)

    ASSERT(associated(this%config))
    ASSERT(associated(this%dnst))
    ASSERT(associated(this%sim))
    nullify(this%sim)
    call interface_xc_stop(this%funct)

    POP_SUB(functional_stop)
  end subroutine functional_stop

  ! ---------------------------------------------------------
  subroutine functional_get_info(this, id, family, kind, polarized, ndim, use)
    type(functional_t), intent(in)  :: this
    integer,  optional, intent(out) :: id
    integer,  optional, intent(out) :: family
    integer,  optional, intent(out) :: kind
    logical,  optional, intent(out) :: polarized
    integer,  optional, intent(out) :: ndim
    logical,  optional, intent(out) :: use

    PUSH_SUB(functional_get_info)

    call interface_xc_get(this%funct, id, family, kind, polarized, ndim, use)

    POP_SUB(functional_get_info)
  end subroutine functional_get_info

  ! ---------------------------------------------------------
  subroutine functional_get_config(this, that)
    type(functional_t),   target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(functional_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(functional_get_config)
  end subroutine functional_get_config

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
  subroutine functional_exc(this, enxc)
    type(functional_t), intent(in)    :: this
    type(storage_t),    intent(inout) :: enxc

    real(kind=wp), dimension(:,:), pointer :: rho
    real(kind=wp), dimension(:),   pointer :: exc
    integer                                :: family

    PUSH_SUB(functional_exc)

    nullify(rho, exc)
    call storage_get(this%dnst, rho)
    ASSERT(associated(rho))
    call storage_get(enxc, exc)
    ASSERT(associated(exc))
    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_exc(this%funct, rho, exc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_exc(this%funct, rho, exc)
    case default
      ASSERT(.false.)
    end select
    nullify(rho, exc)

    POP_SUB(functional_exc)
  end subroutine functional_exc

  ! ---------------------------------------------------------
  subroutine functional_vxc(this, ptxc)
    type(functional_t), intent(in)    :: this
    type(storage_t),    intent(inout) :: ptxc

    real(kind=wp), dimension(:,:), pointer :: rho
    real(kind=wp), dimension(:,:), pointer :: vxc
    integer                                :: family

    PUSH_SUB(functional_vxc)

    nullify(rho, vxc)
    call storage_get(this%dnst, rho)
    ASSERT(associated(rho))
    call storage_get(ptxc, vxc)
    ASSERT(associated(vxc))
    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_vxc(this%funct, rho, vxc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_vxc(this%funct, rho, vxc)
    case default
      ASSERT(.false.)
    end select
    nullify(rho, vxc)

    POP_SUB(functional_vxc)
  end subroutine functional_vxc

  ! ---------------------------------------------------------
  subroutine functional_exc_vxc(this, enxc, ptxc)
    type(functional_t), intent(in)    :: this
    type(storage_t),    intent(inout) :: enxc
    type(storage_t),    intent(inout) :: ptxc

    real(kind=wp), dimension(:,:), pointer :: rho
    real(kind=wp), dimension(:),   pointer :: exc
    real(kind=wp), dimension(:,:), pointer :: vxc
    integer                                :: family

    PUSH_SUB(functional_exc_vxc)

    nullify(rho, exc, vxc)
    call storage_get(this%dnst, rho)
    ASSERT(associated(rho))
    call storage_get(enxc, exc)
    ASSERT(associated(exc))
    call storage_get(ptxc, vxc)
    ASSERT(associated(vxc))
    call functional_get(this, family=family)
    select case(family)
    case(XC_FAMILY_LDA)
      call interface_xc_lda_exc_vxc(this%funct, rho, exc, vxc)
    case(XC_FAMILY_GGA)
      call interface_xc_gga_exc_vxc(this%funct, rho, exc, vxc)
    case default
      ASSERT(.false.)
    end select
    nullify(rho, exc, vxc)

    POP_SUB(functional_exc_vxc)
  end subroutine functional_exc_vxc

  ! ---------------------------------------------------------
  subroutine functional_calc_energy(this, energy)
    type(functional_t), intent(in)  :: this
    real(kind=wp),      intent(out) :: energy

    type(storage_t) :: exc

    PUSH_SUB(functional_calc_energy)

    ASSERT(associated(this%config))
    ASSERT(associated(this%dnst))
    ASSERT(associated(this%sim))
    call storage_init(exc, this%dnst, ndim=1)
    call functional_exc(this, exc)
    call storage_integrate(this%dnst, exc, energy)
    call storage_end(exc)

    POP_SUB(functional_calc_energy)
  end subroutine functional_calc_energy

  ! ---------------------------------------------------------
  subroutine functional_calc_potential(this, potential)
    type(functional_t), intent(in)    :: this
    type(storage_t),    intent(inout) :: potential

    type(storage_t) :: pot
    logical         :: fine

    PUSH_SUB(functional_calc_potential)

    ASSERT(associated(this%config))
    ASSERT(associated(this%dnst))
    ASSERT(associated(this%sim))
    call storage_get(potential, fine=fine)
    if(this%fine.eqv.fine)then
      call functional_vxc(this, potential)
    else
      call storage_init(pot, potential, fine=this%fine)
      call functional_vxc(this, pot)
      call storage_transfer(potential, pot)
      call storage_end(pot)
    end if

    POP_SUB(functional_calc_potential)
  end subroutine functional_calc_potential

  ! ---------------------------------------------------------
  subroutine functional_calc_energy_and_potential(this, energy, potential)
    type(functional_t), intent(in)    :: this
    real(kind=wp),      intent(out)   :: energy
    type(storage_t),    intent(inout) :: potential

    type(storage_t) :: exc, pot
    logical         :: fine

    PUSH_SUB(functional_calc_energy_and_potential)

    ASSERT(associated(this%config))
    ASSERT(associated(this%dnst))
    ASSERT(associated(this%sim))
    call storage_init(exc, this%dnst, ndim=1)
    call storage_get(potential, fine=fine)
    if(this%fine.eqv.fine)then
      call functional_exc_vxc(this, exc, potential)
    else
      call storage_init(pot, potential, fine=this%fine)
      call functional_exc_vxc(this, exc, pot)
      call storage_transfer(potential, pot)
      call storage_end(pot)
    end if
    call storage_integrate(this%dnst, exc, energy)
    call storage_end(exc)

    POP_SUB(functional_calc_energy_and_potential)
  end subroutine functional_calc_energy_and_potential

  ! ---------------------------------------------------------
  subroutine functional_copy(this, that)
    type(functional_t), intent(inout) :: this
    type(functional_t), intent(in)    :: that

    PUSH_SUB(functional_copy)

    this%config => that%config
    this%dnst => that%dnst
    this%sim => that%sim
    this%fine = that%fine
    call interface_xc_copy(this%funct, that%funct)

    POP_SUB(functional_copy)
  end subroutine functional_copy

  ! ---------------------------------------------------------
  subroutine functional_end(this)
    type(functional_t), intent(inout) :: this

    PUSH_SUB(functional_end)

    nullify(this%config, this%dnst, this%sim)
    this%fine = .false.
    call interface_xc_end(this%funct)

    POP_SUB(functional_end)
  end subroutine functional_end

end module functional_oct_m

!! Local Variables:
!! mode: f90
!! End:
