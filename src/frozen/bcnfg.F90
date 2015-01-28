#include "global.h"

module bcnfg_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m, only: CURV_METHOD_UNIFORM
  use bhmlt_m,       only: HMLT_TYPE
  use intrpl_m,      only: NEAREST
  use json_m,        only: JSON_OK
  use json_m,        only: json_object_t, json_array_t
  use json_m,        only: json_isdef, json_init, json_set, json_get

  implicit none

  private
  public ::       &
    bcnfg_parse

  integer, parameter, public :: NONE_TYPE = 0
  integer, parameter, public :: FNIO_TYPE = 1
  integer, parameter, public :: FRZN_TYPE = 2
  integer, parameter, public :: MAIN_TYPE = 3

  integer, parameter, public :: NAME_LEN = 63

  integer, parameter :: default_ndim  = 3
  integer, parameter :: default_nspin = 1

contains

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_simul_box(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    !
    integer :: idim
    !
    idim=default_ndim
    if(present(ndim))idim=ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)
    return
  end subroutine bcnfg_parse_simul_box

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_curvilinear(this, method)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: method
    !
    integer :: mthd
    !
    mthd=CURV_METHOD_UNIFORM
    if(present(method))mthd=method
    call json_init(this)
    call json_set(this, "method", mthd)
    return
  end subroutine bcnfg_parse_curvilinear

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_mesh(this)
    type(json_object_t), intent(out) :: this
    !
    type(json_array_t), pointer :: list
    !
    nullify(list)
    call json_init(this)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "spacing", list)
    nullify(list)
    return
  end subroutine bcnfg_parse_mesh

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_grid(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_simul_box(cnfg, ndim)
    call json_set(this, "simul_box", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_curvilinear(cnfg)
    call json_set(this, "curvilinear", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_mesh(cnfg)
    call json_set(this, "mesh", cnfg)
    nullify(cnfg)
    return
  end subroutine bcnfg_parse_grid

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_simulation(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_grid(cnfg, ndim)
    call json_set(this, "grid", cnfg)
    nullify(cnfg)
    return
  end subroutine bcnfg_parse_simulation

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_space(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    !
    integer :: idim
    !
    idim=default_ndim
    if(present(ndim))idim=ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)
    return
  end subroutine bcnfg_parse_space

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_geometry(this)
    type(json_object_t), intent(out) :: this
    !
    type(json_array_t), pointer :: list
    !
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
    return
  end subroutine bcnfg_parse_geometry

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_density(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    !
    integer :: ispin
    !
    ispin=default_nspin
    if(present(nspin))ispin=nspin
    ASSERT(ispin>0)
    call json_init(this)
    call json_set(this, "nspin", ispin)
    return
  end subroutine bcnfg_parse_density

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_states(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_density(cnfg, nspin)
    call json_set(this, "density", cnfg)
    nullify(cnfg)
    return
  end subroutine bcnfg_parse_states

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_system(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_space(cnfg, ndim)
    call json_set(this, "space", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_geometry(cnfg)
    call json_set(this, "geometry", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_states(cnfg, nspin)
    call json_set(this, "states", cnfg)
    nullify(cnfg)
    return
  end subroutine bcnfg_parse_system

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_hamiltonian(this)
    type(json_object_t), intent(out) :: this
    !
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE)
    return
  end subroutine bcnfg_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine bcnfg_parse_model(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    !
    PUSH_SUB(bcnfg_parse_model)
    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_simulation(cnfg, ndim)
    call json_set(this, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_system(cnfg, ndim, nspin)
    call json_set(this, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_hamiltonian(cnfg)
    call json_set(this, "hamiltonian", cnfg)
    nullify(cnfg)
    POP_SUB(bcnfg_parse_model)
    return
  end subroutine bcnfg_parse_model

  ! ---------------------------------------------------------
  subroutine bcnfg_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    !
    PUSH_SUB(bcnfg_parse)
    nullify(cnfg, list)
    call json_init(this)
    call json_set(this, "type", NONE_TYPE)
    call json_set(this, "name", "base")
    call json_set(this, "interpolation", NEAREST)
    SAFE_ALLOCATE(cnfg)
    call bcnfg_parse_model(cnfg, ndim, nspin)
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
    POP_SUB(bcnfg_parse)
    return
  end subroutine bcnfg_parse

end module bcnfg_m

!! Local Variables:
!! mode: f90
!! End:

