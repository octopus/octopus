#include "global.h"

module live_config_m

  use base_config_m
  use base_hamiltonian_m
  use global_m
  use json_m
  use live_handle_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::           &
    live_config_parse

contains

  ! ---------------------------------------------------------
  subroutine live_config_parse_density(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    PUSH_SUB(live_config_parse_density)

    POP_SUB(live_config_parse_density)
  end subroutine live_config_parse_density

  ! ---------------------------------------------------------
  subroutine live_config_parse_states(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_states)

    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_density(cnfg, nspin)
    nullify(cnfg)

    POP_SUB(live_config_parse_states)
  end subroutine live_config_parse_states

  ! ---------------------------------------------------------
  subroutine live_config_parse_system(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_system)

    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_states(cnfg, nspin)
    nullify(cnfg)

    POP_SUB(live_config_parse_system)
  end subroutine live_config_parse_system

  ! ---------------------------------------------------------
  subroutine live_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(live_config_parse_external)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)

    POP_SUB(live_config_parse_external)
  end subroutine live_config_parse_external

  ! ---------------------------------------------------------
  subroutine live_config_parse_ionic(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(live_config_parse_ionic)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(live_config_parse_ionic)
  end subroutine live_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine live_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(live_config_parse_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call live_config_parse_external(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call live_config_parse_ionic(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_hamiltonian)
  end subroutine live_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine live_config_parse_model(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_model)

    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_system(cnfg, nspin)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_hamiltonian(cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_model)
  end subroutine live_config_parse_model

  ! ---------------------------------------------------------
  subroutine live_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse)

    nullify(cnfg)
    call base_config_parse(this, ndim, nspin)
    call json_set(this, "type", HNDL_TYPE_LIVE)
    call json_set(this, "name", "live")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_model(cnfg, nspin)
    nullify(cnfg)

    POP_SUB(live_config_parse)
  end subroutine live_config_parse

end module live_config_m

!! Local Variables:
!! mode: f90
!! End:
