#include "global.h"

module frozen_config_m

  use base_config_m
  use base_hamiltonian_m
  use frozen_handle_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private

  public ::              &
    frozen_config_parse

contains

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(frozen_config_parse_external)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)

    POP_SUB(frozen_config_parse_external)
  end subroutine frozen_config_parse_external

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_ionic(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(frozen_config_parse_ionic)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(frozen_config_parse_ionic)
  end subroutine frozen_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(frozen_config_parse_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_config_parse_external(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call frozen_config_parse_ionic(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse_hamiltonian)
  end subroutine frozen_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_model(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse_model)

    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_hamiltonian(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse_model)
  end subroutine frozen_config_parse_model

  ! ---------------------------------------------------------
  subroutine frozen_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse)

    nullify(cnfg)
    call base_config_parse(this, ndim, nspin)
    call json_set(this, "type", HNDL_TYPE_FRZN)
    call json_set(this, "name", "frozen")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_model(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse)
  end subroutine frozen_config_parse

end module frozen_config_m

!! Local Variables:
!! mode: f90
!! End:
