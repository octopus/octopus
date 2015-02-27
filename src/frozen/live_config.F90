#include "global.h"

module live_config_m

  use global_m
  use messages_m
  use profiling_m

  use geometry_m, only: geometry_t, geometry_create_data_object
  use intrpl_m,   only: NEAREST
  use json_m,     only: JSON_OK, json_object_t
  use json_m,     only: json_init, json_set, json_get, json_del
  use kinds_m,    only: wp

  use base_hamiltonian_m, only: &
    HMLT_TYPE_POTN,             &
    HMLT_TYPE_HMLT

  use base_config_m, only: &
    base_config_parse

  use live_handle_m, only: &
    HNDL_TYPE_LIVE

  implicit none

  private
  public ::           &
    live_config_parse

contains

  ! ---------------------------------------------------------
  subroutine live_config_parse_simulation(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: ierr
    !
    call json_del(this, "grid", ierr)
    ASSERT(ierr==JSON_OK)
    return
  end subroutine live_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine live_config_parse_geometry(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    !
    call geometry_create_data_object(geo, this)
    return
  end subroutine live_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine live_config_parse_system(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
   call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_geometry(cnfg, geo)
    nullify(cnfg)
    return
  end subroutine live_config_parse_system

  ! ---------------------------------------------------------
  subroutine live_config_parse_external(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: type, ierr
    !
    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)
    return
  end subroutine live_config_parse_external

  ! ---------------------------------------------------------
  subroutine live_config_parse_hamiltonian(this)
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
    call live_config_parse_external(cnfg)
    nullify(cnfg  )
    return
  end subroutine live_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine live_config_parse_model(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(live_config_parse_model)
    nullify(cnfg)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_simulation(cnfg)
    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_system(cnfg, geo)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_hamiltonian(cnfg)
    nullify(cnfg)
    POP_SUB(live_config_parse_model)
    return
  end subroutine live_config_parse_model

  ! ---------------------------------------------------------
  subroutine live_config_parse(this, geo, ndim, nspin)
    type(json_object_t), intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(live_config_parse)
    nullify(cnfg)
    call base_config_parse(this, ndim, nspin)
    call json_set(this, "name", "live")
    call json_set(this, "type", HNDL_TYPE_LIVE)
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_model(cnfg, geo)
    nullify(cnfg)
    POP_SUB(live_config_parse)
    return
  end subroutine live_config_parse

end module live_config_m

!! Local Variables:
!! mode: f90
!! End:
