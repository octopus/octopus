#include "global.h"

module frozen_config_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t
  use json_m, only: json_init, json_set, json_get, json_del

  use base_hamiltonian_m, only: &
    HMLT_TYPE_POTN,             &
    HMLT_TYPE_HMLT

  use base_config_m, only: &
    base_config_parse

  use frozen_handle_m, only: &
    HNDL_TYPE_FRZN

  implicit none

  private
  public ::              &
    frozen_config_parse

contains

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_simulation(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: ierr
    !
    call json_del(this, "grid", ierr)
    ASSERT(ierr==JSON_OK)
    return
  end subroutine frozen_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_external(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: type, ierr
    !
    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)
    return
  end subroutine frozen_config_parse_external

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr
    !
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
    return
  end subroutine frozen_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_model(this)
    type(json_object_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_simulation(cnfg)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_hamiltonian(cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_config_parse_model

  ! ---------------------------------------------------------
  subroutine frozen_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call base_config_parse(this, ndim, nspin)
    call json_set(this, "type", HNDL_TYPE_FRZN)
    call json_set(this, "name", "frozen")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_model(cnfg)
    nullify(cnfg)
    return
  end subroutine frozen_config_parse

end module frozen_config_m

!! Local Variables:
!! mode: f90
!! End:
