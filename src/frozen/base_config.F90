#include "global.h"

module base_config_oct_m

  use base_hamiltonian_oct_m
  use base_handle_oct_m
  use global_oct_m
  use intrpl_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m

  implicit none

  private

  public ::               &
    BASE_CONFIG_NAME_LEN
  
  public ::            &
    base_config_use,   &
    base_config_add,   &
    base_config_parse

  integer, parameter :: BASE_CONFIG_NAME_LEN = BASE_HANDLE_NAME_LEN
  
  integer, parameter :: default_ndim  = 3
  integer, parameter :: default_nspin = 1

contains

  ! ---------------------------------------------------------
  function base_config_use(this) result(that)
    type(json_object_t), intent(in) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    logical :: that

    PUSH_SUB(base_config_use)

    nullify(cnfg)
    call json_get(this, "subsystems", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    that = (json_len(cnfg)>0)
    nullify(cnfg)

    POP_SUB(base_config_use)
  end function base_config_use

  ! ---------------------------------------------------------
  subroutine base_config_add(this, name, that)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(json_object_t), intent(in)    :: that

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_config_add)

    nullify(cnfg)
    call json_get(this, "subsystems", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, name, that)
    nullify(cnfg)

    POP_SUB(base_config_add)
  end subroutine base_config_add

  ! ---------------------------------------------------------
  subroutine base_config_parse_space(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    integer :: idim

    PUSH_SUB(base_config_parse_space)

    idim = default_ndim
    if(present(ndim)) idim = ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)

    POP_SUB(base_config_parse_space)
  end subroutine base_config_parse_space

  ! ---------------------------------------------------------
  subroutine base_config_parse_molecule(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(base_config_parse_molecule)

    call json_init(this)

    POP_SUB(base_config_parse_molecule)
  end subroutine base_config_parse_molecule

  ! ---------------------------------------------------------
  subroutine base_config_parse_geometry(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_geometry)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "reduce", .false.)
    call json_set(this, "external", .false.)
    call json_set(this, "default", .true.)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_molecule(cnfg)
    call json_set(this, "molecule", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_geometry)
  end subroutine base_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine base_config_parse_density(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    type(json_object_t),             pointer :: cnfg
    real(kind=wp), dimension(:), allocatable :: chrg
    integer                                  :: ispin

    PUSH_SUB(base_config_parse_density)

    nullify(cnfg)
    ispin = default_nspin
    if(present(nspin)) ispin = nspin
    ASSERT(ispin>0)
    SAFE_ALLOCATE(chrg(1:ispin))
    chrg = 0.0_wp
    call json_init(this)
    call json_set(this, "reduce", .false.)
    call json_set(this, "external", .false.)
    call json_set(this, "default", .true.)
    call json_set(this, "nspin", ispin)
    call json_set(this, "charge", chrg)
    SAFE_DEALLOCATE_A(chrg)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, ndim=ispin, fine=.true.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_density)
  end subroutine base_config_parse_density

  ! ---------------------------------------------------------
  subroutine base_config_parse_states(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_states)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_density(cnfg, nspin)
    call json_set(this, "density", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_states)
  end subroutine base_config_parse_states

  ! ---------------------------------------------------------
  subroutine base_config_parse_system(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_system)

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

    POP_SUB(base_config_parse_system)
  end subroutine base_config_parse_system

  ! ---------------------------------------------------------
  subroutine base_config_parse_hamiltonian(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_hamiltonian)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_HMLT)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false., allocate=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_hamiltonian)
  end subroutine base_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_config_parse_model(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_model)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_system(cnfg, nspin, ndim)
    call json_set(this, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_hamiltonian(cnfg)
    call json_set(this, "hamiltonian", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_model)
  end subroutine base_config_parse_model

  ! ---------------------------------------------------------
  subroutine base_config_parse(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    
    PUSH_SUB(base_config_parse)

    nullify(cnfg, list)
    call json_init(this)
    call json_set(this, "type", HNDL_TYPE_NONE)
    call json_set(this, "name", "base")
    call json_set(this, "interpolation", NEAREST)
    SAFE_ALLOCATE(cnfg)
    call simulation_init(cnfg)
    call json_set(this, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_model(cnfg, nspin, ndim)
    call json_set(this, "model", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "positions", list)
    nullify(list)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(this, "subsystems", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse)
  end subroutine base_config_parse

end module base_config_oct_m

!! Local Variables:
!! mode: f90
!! End:

