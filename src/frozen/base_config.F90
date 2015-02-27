#include "global.h"

module base_config_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m
  use intrpl_m
  use json_m
  
  use base_hamiltonian_m, only: &
    HMLT_TYPE_HMLT

  use base_handle_m, only: &
    HNDL_TYPE_NONE

  implicit none

  private
  public ::            &
    base_config_parse

  integer, parameter :: default_ndim  = 3
  integer, parameter :: default_nspin = 1

contains

  ! ---------------------------------------------------------
  subroutine base_config_parse_simul_box(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    integer :: idim

    idim = default_ndim
    if(present(ndim))idim = ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)

  end subroutine base_config_parse_simul_box

  ! ---------------------------------------------------------
  subroutine base_config_parse_curvilinear(this, method)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: method

    integer :: mthd

    mthd = CURV_METHOD_UNIFORM
    if(present(method))mthd = method
    call json_init(this)
    call json_set(this, "method", mthd)

  end subroutine base_config_parse_curvilinear

  ! ---------------------------------------------------------
  subroutine base_config_parse_mesh(this)
    type(json_object_t), intent(out) :: this
    
    type(json_array_t), pointer :: list

    nullify(list)
    call json_init(this)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "spacing", list)
    nullify(list)

  end subroutine base_config_parse_mesh

  ! ---------------------------------------------------------
  subroutine base_config_parse_grid(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_simul_box(cnfg, ndim)
    call json_set(this, "simul_box", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_curvilinear(cnfg)
    call json_set(this, "curvilinear", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_mesh(cnfg)
    call json_set(this, "mesh", cnfg)
    nullify(cnfg)

  end subroutine base_config_parse_grid

  ! ---------------------------------------------------------
  subroutine base_config_parse_simulation(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_grid(cnfg, ndim)
    call json_set(this, "grid", cnfg)
    nullify(cnfg)

  end subroutine base_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine base_config_parse_space(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    integer :: idim

    idim = default_ndim
    if(present(ndim))idim = ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)

  end subroutine base_config_parse_space

  ! ---------------------------------------------------------
  subroutine base_config_parse_geometry(this)
    type(json_object_t), intent(out) :: this

    type(json_array_t), pointer :: list

    call json_init(this)
    call json_set(this, "nspecies", 0)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "species", list)
    nullify(list)
    call json_set(this, "natoms", 0)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "atom", list)
    nullify(list)
    call json_set(this, "ncatoms", 0)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "catom", list)
    nullify(list)

  end subroutine base_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine base_config_parse_density(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    integer :: ispin

    ispin = default_nspin
    if(present(nspin)) ispin = nspin
    ASSERT(ispin>0)
    call json_init(this)
    call json_set(this, "nspin", ispin)

  end subroutine base_config_parse_density

  ! ---------------------------------------------------------
  subroutine base_config_parse_states(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_density(cnfg, nspin)
    call json_set(this, "density", cnfg)
    nullify(cnfg)

  end subroutine base_config_parse_states

  ! ---------------------------------------------------------
  subroutine base_config_parse_system(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_space(cnfg, ndim)
    call json_set(this, "space", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_geometry(cnfg)
    call json_set(this, "geometry", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_states(cnfg, nspin)
    call json_set(this, "states", cnfg)
    nullify(cnfg)

  end subroutine base_config_parse_system

  ! ---------------------------------------------------------
  subroutine base_config_parse_hamiltonian(this)
    type(json_object_t), intent(out) :: this

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_HMLT)

  end subroutine base_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_config_parse_model(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_model)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_simulation(cnfg, ndim)
    call json_set(this, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_system(cnfg, ndim, nspin)
    call json_set(this, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_hamiltonian(cnfg)
    call json_set(this, "hamiltonian", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_model)
  end subroutine base_config_parse_model

  ! ---------------------------------------------------------
  subroutine base_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    
    PUSH_SUB(base_config_parse)

    nullify(cnfg, list)
    call json_init(this)
    call json_set(this, "type", HNDL_TYPE_NONE)
    call json_set(this, "name", "base")
    call json_set(this, "interpolation", NEAREST)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_model(cnfg, ndim, nspin)
    call json_set(this, "model", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "positions", list)
    nullify(list)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "systems", list)
    nullify(list)

    POP_SUB(base_config_parse)
  end subroutine base_config_parse

end module base_config_m

!! Local Variables:
!! mode: f90
!! End:

