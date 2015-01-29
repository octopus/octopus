#include "global.h"

module fio_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_get

  use bhmlt_m, only:                &
    hamiltonian_setn => bhmlt_setn, &
    hamiltonian_getn => bhmlt_getn, &
    hamiltonian_deln => bhmlt_deln

  use bhmlt_m, only:                        &
    fio_hamiltonian_t      => bhmlt_t,      &
    fio_hamiltonian_start  => bhmlt_start,  &
    fio_hamiltonian_stop   => bhmlt_stop,   &
    fio_hamiltonian_get    => bhmlt_get

  use fio_external_m, only: &
    fio_external_t,         &
    fio_external_init,      &
    fio_external_update,    &
    fio_external_copy,      &
    fio_external_end

  use fio_simulation_m, only: &
    fio_simulation_t

  use fio_system_m, only: &
    fio_system_t

  implicit none

  private
  public ::                 &
    fio_hamiltonian_t,      &
    fio_hamiltonian_init,   &
    fio_hamiltonian_start,  &
    fio_hamiltonian_update, &
    fio_hamiltonian_stop,   &
    fio_hamiltonian_get,    &
    fio_hamiltonian_copy,   &
    fio_hamiltonian_end
  
  interface fio_hamiltonian_get
    module procedure fio_hamiltonian_get_external
  end interface !fio_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_init(this, sys, config)
    type(fio_hamiltonian_t), intent(inout) :: this
    type(fio_system_t),      intent(in)    :: sys
    type(json_object_t),     intent(in)    :: config
    !
    type(json_object_t),  pointer :: cnfg
    type(fio_external_t), pointer :: epot
    integer                       :: ierr
    !
    PUSH_SUB(fio_hamiltonian_init)
    nullify(cnfg, epot)
    call json_get(config, "external", cnfg, ierr)
    if(ierr==JSON_OK)then
      SAFE_ALLOCATE(epot)
      call fio_external_init(epot, sys, cnfg)
      call hamiltonian_setn(this, "external", epot)
    end if
    POP_SUB(fio_hamiltonian_init)
    return
  end subroutine fio_hamiltonian_init

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_update(this)
    type(fio_hamiltonian_t), intent(inout) :: this
    !
    type(fio_external_t), pointer :: epot
    !
    PUSH_SUB(fio_hamiltonian_update)
    nullify(epot)
    call fio_hamiltonian_get_external(this, epot)
    ASSERT(associated(epot))
    call fio_external_update(epot)
    POP_SUB(fio_hamiltonian_update)
    return
  end subroutine fio_hamiltonian_update

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_external(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_external_t),   pointer     :: that
    !
    PUSH_SUB(fio_hamiltonian_get_external)
    call hamiltonian_getn(this, "external", that)
    POP_SUB(fio_hamiltonian_get_external)
    return
  end subroutine fio_hamiltonian_get_external

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_copy(this, that)
    type(fio_hamiltonian_t), intent(inout) :: this
    type(fio_hamiltonian_t), intent(in)    :: that
    !
    type(fio_external_t), pointer :: oept, iept
    !
    PUSH_SUB(fio_hamiltonian_copy)
    nullify(oept, iept)
    call fio_hamiltonian_end(this)
    call fio_hamiltonian_get_external(that, iept)
    if(associated(iept))then
      SAFE_ALLOCATE(oept)
      call fio_external_copy(oept, iept)
      call hamiltonian_setn(this, "external", oept)
    end if
    POP_SUB(fio_hamiltonian_copy)
    return
  end subroutine fio_hamiltonian_copy

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_end(this)
    type(fio_hamiltonian_t), intent(inout) :: this
    !
    type(fio_external_t), pointer :: epot
    !
    PUSH_SUB(fio_hamiltonian_end)
    nullify(epot)
    call hamiltonian_deln(this, "external", epot)
    if(associated(epot))then
      call fio_external_end(epot)
      SAFE_DEALLOCATE_P(epot)
      nullify(epot)
    end if
    POP_SUB(fio_hamiltonian_end)
    return
  end subroutine fio_hamiltonian_end

end module fio_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
