#include "global.h"

module frozen_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  use base_handle_oct_m
  use base_potential_oct_m
  use fio_config_oct_m
  use fio_handle_oct_m
  use frozen_handle_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use parser_oct_m
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

    call json_set(this, "default", .true.)
    call json_set(this, "reduce", .true.)

    POP_SUB(frozen_config_parse_geometry)
  end subroutine frozen_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_density(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(frozen_config_parse_density)

    call json_set(this, "default", .true.)

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
    call json_set(this, "type", TERM_TYPE_POTN)
    call json_set(this, "kind", POTN_KIND_XTRN)
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
    call json_set(this, "type", TERM_TYPE_TERM)

    POP_SUB(frozen_config_parse_ionic)
  end subroutine frozen_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine frozen_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: trms, cnfg, defn
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse_hamiltonian)

    nullify(trms, cnfg, defn)
    call json_get(this, "terms", trms, ierr)
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(trms, "external", cnfg)
    SAFE_ALLOCATE(defn)
    call frozen_config_parse_external(defn)
    call json_set(cnfg, "definition", defn)
    nullify(cnfg, defn)
    SAFE_ALLOCATE(cnfg)
    call json_init(cnfg)
    call json_set(trms, "ionic", cnfg)
    SAFE_ALLOCATE(defn)
    call frozen_config_parse_ionic(defn)
    call json_set(cnfg, "definition", defn)
    nullify(trms, cnfg, defn)

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
  subroutine frozen_config_parse(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    integer,             intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_config_parse)

    nullify(cnfg)
    call base_config_parse(this, nspin, ndim, parse)
    call json_set(this, "type", HNDL_TYPE_FRZN)
    call json_set(this, "name", "frozen")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_config_parse_model(cnfg)
    nullify(cnfg)

    POP_SUB(frozen_config_parse)

  contains
    
    subroutine parse(this, type, block, line, ierr)
      type(json_object_t), intent(out) :: this
      integer,             intent(in)  :: type
      type(block_t),       intent(in)  :: block
      integer,             intent(in)  :: line
      integer,             intent(out) :: ierr

      PUSH_SUB(frozen_config_parse.parse)

      ierr = PARSE_FAIL
      if(type==HNDL_TYPE_FRZN)then
        call fio_config_parse(this, block, line)
        ierr = PARSE_OK
      end if
      
      POP_SUB(frozen_config_parse.parse)
    end subroutine parse
      
  end subroutine frozen_config_parse

end module frozen_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
