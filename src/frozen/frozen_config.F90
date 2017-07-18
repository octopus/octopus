#include "global.h"

module frozen_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  use fio_handle_oct_m
  use frozen_handle_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use storage_oct_m

  implicit none

  private

  public ::              &
    frozen_config_parse

contains

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_geometry(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(frozen_config_parse_geometry)

    call json_set(this, "default", .false.)

    POP_SUB(frozen_config_parse_geometry)
  end subroutine frozen_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_density(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(frozen_config_parse_density)

    call json_set(this, "default", .true.)
    call json_set(this, "external", .true.)

    POP_SUB(frozen_config_parse_density)
  end subroutine frozen_config_parse_density

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_states(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse_states)

    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_density(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse_states)
  end subroutine frozen_config_parse_states

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_system(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse_system)

    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_geometry(cnfg)
    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_states(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse_system)
  end subroutine frozen_config_parse_system

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(frozen_config_parse_external)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

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
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_system(cnfg)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_hamiltonian(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse_model)
  end subroutine frozen_config_parse_model

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_positions(this)
    type(json_object_t), intent(in) :: this

    type(json_array_t),         pointer :: list
    character(len=BASE_CONFIG_NAME_LEN) :: name
    integer                             :: ierr
    
    PUSH_SUB(frozen_config_parse_positions)

    nullify(list)
    call json_get(this, "name", name, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    if(json_len(list)<1)then
      message(1) = "Frozen subsystem '"//trim(adjustl(name))//"' position(s) were not specified in the 'SubSystemCoordinates' block."
      call messages_fatal(1)
    end if
    nullify(list)
    
    POP_SUB(frozen_config_parse_positions)
  end subroutine frozen_config_parse_positions

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_subsystems(this, dict)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: dict

    type(json_object_iterator_t)        :: iter
    character(len=BASE_CONFIG_NAME_LEN) :: name
    type(json_object_t),        pointer :: ocfg, icfg
    integer                             :: type, ierr
    
    PUSH_SUB(frozen_config_parse_subsystems)

    call json_init(iter, dict)
    do
      nullify(ocfg, icfg)
      call json_next(iter, name, icfg, ierr)
      if(ierr/=JSON_OK)exit
      call json_get(icfg, "type", type, ierr)
      ASSERT(ierr==JSON_OK)
      if(type/=HNDL_TYPE_FNIO)cycle
      SAFE_ALLOCATE(ocfg)
      call json_copy(ocfg, icfg)
      call json_set(this, trim(adjustl(name)), ocfg)
      call frozen_config_parse_positions(ocfg)
    end do
    call json_end(iter)
    nullify(ocfg, icfg)

    POP_SUB(frozen_config_parse_subsystems)
  end subroutine frozen_config_parse_subsystems

  ! ---------------------------------------------------------
  subroutine frozen_config_parse(this, dict, nspin, ndim)
    type(json_object_t), intent(out) :: this
    type(json_object_t), intent(in)  :: dict
    integer,             intent(in)  :: nspin
    integer,             intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse)

    nullify(cnfg)
    call base_config_parse(this, nspin, ndim)
    call json_set(this, "type", HNDL_TYPE_FRZN)
    call json_set(this, "name", "frozen")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_model(cnfg)
    nullify(cnfg)
    call json_get(this, "subsystems", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_subsystems(cnfg, dict)
    nullify(cnfg)

    POP_SUB(frozen_config_parse)
  end subroutine frozen_config_parse

end module frozen_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
