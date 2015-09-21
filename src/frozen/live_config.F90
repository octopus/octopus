#include "global.h"

module live_config_m

  use base_config_m
  use base_hamiltonian_m
  use functional_m
  use geometry_m
  use global_m
  use json_m
  use kinds_m
  use live_handle_m
  use messages_m
  use parser_m
  use profiling_m

  implicit none

  private

  public ::           &
    live_config_parse

contains

  ! ---------------------------------------------------------
  subroutine live_config_parse_simulation(this)
    type(json_object_t), intent(inout) :: this

    integer :: ierr

    PUSH_SUB(live_config_parse_simulation)

    call json_del(this, "grid", ierr)
    ASSERT(ierr==JSON_OK)

    POP_SUB(live_config_parse_simulation)
  end subroutine live_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine live_config_parse_geometry(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo

    PUSH_SUB(live_config_parse_geometry)

    call geometry_create_data_object(geo, this)

    POP_SUB(live_config_parse_geometry)
  end subroutine live_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine live_config_parse_system(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_system)

    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_geometry(cnfg, geo)
    nullify(cnfg)

    POP_SUB(live_config_parse_system)
  end subroutine live_config_parse_system

  ! ---------------------------------------------------------
  subroutine live_config_parse_external(this)
    type(json_object_t), intent(inout) :: this

    integer :: type, ierr

    PUSH_SUB(live_config_parse_external)

    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)

    POP_SUB(live_config_parse_external)
  end subroutine live_config_parse_external

  ! ---------------------------------------------------------
  subroutine live_config_parse_ionic(this)
    type(json_object_t), intent(inout) :: this

    integer :: type, ierr

    PUSH_SUB(live_config_parse_ionic)

    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(live_config_parse_ionic)
  end subroutine live_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine live_config_parse_tnadd(this)
    type(json_object_t), intent(inout) :: this

    real(kind=wp) :: factor
    integer       :: type, id, ierr

    PUSH_SUB(live_config_parse_tnadd)

    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_FNCT)
    call parse_variable("TnaddFunctional", FUNCT_XC_NONE, id)
    call json_set(this, "functional", id)
    if(id>FUNCT_XC_NONE)then
      call parse_variable("TnaddFactor", 1.0_wp, factor)
      call json_set(this, "factor", factor)
    end if

    POP_SUB(live_config_parse_tnadd)
  end subroutine live_config_parse_tnadd

  ! ---------------------------------------------------------
  subroutine live_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr

    PUSH_SUB(live_config_parse_hamiltonian)

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
    nullify(cnfg)
    call json_get(this, "ionic", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "ionic", cnfg)
    end if
    call live_config_parse_ionic(cnfg)
    nullify(cnfg)
    call json_get(this, "tnadd", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "tnadd", cnfg)
    end if
    call live_config_parse_tnadd(cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_hamiltonian)
  end subroutine live_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine live_config_parse_model(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

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
  end subroutine live_config_parse_model

  ! ---------------------------------------------------------
  subroutine live_config_parse(this, geo, ndim, nspin)
    type(json_object_t), intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

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
  end subroutine live_config_parse

end module live_config_m

!! Local Variables:
!! mode: f90
!! End:
