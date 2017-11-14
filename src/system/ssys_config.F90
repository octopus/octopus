#include "global.h"

module ssys_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  ! use base_handle_oct_m
  ! use fio_config_oct_m
  ! use fio_handle_oct_m
  use frozen_config_oct_m
  ! use frozen_handle_oct_m
  use functional_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use live_config_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ssys_handle_oct_m
  use space_oct_m
  use states_oct_m
  use storage_oct_m

  implicit none

  private

  public ::            &
    ssys_config_parse

contains

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_geometry(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(ssys_config_parse_geometry)

    call json_set(this, "reduce", .true.)

    POP_SUB(ssys_config_parse_geometry)
  end subroutine ssys_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_density(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(ssys_config_parse_density)

    call json_set(this, "reduce", .true.)

    POP_SUB(ssys_config_parse_density)
  end subroutine ssys_config_parse_density

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_states(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse_states)

    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_density(cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_states)
  end subroutine ssys_config_parse_states

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_system(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse_system)

    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_geometry(cnfg)
    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_states(cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_system)
  end subroutine ssys_config_parse_system

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_external)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", TERM_TYPE_POTN)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_external)
  end subroutine ssys_config_parse_external

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_ionic(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(ssys_config_parse_ionic)

    call json_init(this)
    call json_set(this, "type", TERM_TYPE_TERM)

    POP_SUB(ssys_config_parse_ionic)
  end subroutine ssys_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_tnadd(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
 
    type(json_object_t), pointer :: cnfg
    logical                      :: plrz

    PUSH_SUB(ssys_config_parse_tnadd)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", TERM_TYPE_HMLT)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false., allocate=.true.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)
    call parse_variable('TnaddPolarized', (nspin>1), plrz)
    call json_set(this, "spin", plrz)

    POP_SUB(ssys_config_parse_tnadd)
  end subroutine ssys_config_parse_tnadd

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_hamiltonian(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_tnadd(cnfg, nspin)
    call json_set(this, "tnadd", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_kinetic(cnfg, nspin)
    call json_set(this, "kinetic", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_external(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_ionic(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_hamiltonian)
  end subroutine ssys_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_model(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse_model)

    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_system(cnfg)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_hamiltonian(cnfg, nspin)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_model)
  end subroutine ssys_config_parse_model

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_subsystems(this, st, space)
    type(json_object_t), intent(out) :: this
    type(states_t),      intent(in)  :: st
    type(space_t),       intent(in)  :: space
    
    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_subsystems)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call frozen_config_parse(cnfg, st%d%nspin, space%dim)
    if(json_len(cnfg)>0) call json_set(this, "frozen", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call live_config_parse(cnfg, st, space)
    call json_set(this, "live", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_subsystems)
  end subroutine ssys_config_parse_subsystems

  ! ---------------------------------------------------------
  subroutine ssys_config_parse(this, st, space)
    type(json_object_t), intent(out) :: this
    type(states_t),      intent(in)  :: st
    type(space_t),       intent(in)  :: space

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse)

    nullify(cnfg)
    call base_config_parse(this, st%d%nspin, space%dim)
    call json_set(this, "type", HNDL_TYPE_SSYS)
    call json_set(this, "name", "sigma")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_model(cnfg, st%d%nspin)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_subsystems(cnfg, st, space)
    ASSERT(json_len(cnfg)>0)
    call json_set(this, "subsystems", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse)
  end subroutine ssys_config_parse

end module ssys_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
