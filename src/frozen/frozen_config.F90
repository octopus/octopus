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
  subroutine frozen_config_parse_simulation(this)
    type(json_object_t), intent(inout) :: this

    integer :: ierr

    PUSH_SUB(frozen_config_parse_simulation)

    call json_del(this, "grid", ierr)
    ASSERT(ierr==JSON_OK)

    POP_SUB(frozen_config_parse_simulation)
  end subroutine frozen_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_external(this)
    type(json_object_t), intent(inout) :: this

    integer :: type, nspin, ierr

    PUSH_SUB(frozen_config_parse_external)

    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)
    call json_get(this, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)call json_set(this, "nspin", 1)

    POP_SUB(frozen_config_parse_external)
  end subroutine frozen_config_parse_external

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_ionic(this)
    type(json_object_t), intent(inout) :: this

    integer :: type, ierr

    PUSH_SUB(frozen_config_parse_ionic)

    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(frozen_config_parse_ionic)
  end subroutine frozen_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr

    PUSH_SUB(frozen_config_parse_hamiltonian)

    nullify(cnfg)
    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_HMLT)
    call json_get(this, "external", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "external", cnfg)
    end if
    call frozen_config_parse_external(cnfg)
    nullify(cnfg)
    call json_get(this, "ionic", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "ionic", cnfg)
    end if
    call frozen_config_parse_ionic(cnfg)
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
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_simulation(cnfg)
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
