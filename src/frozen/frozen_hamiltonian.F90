#include "global.h"

module frozen_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_get

  use fio_external_m, only: &
    fio_external_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t,         &
    fio_hamiltonian_get

  use frozen_system_m, only: &
    frozen_system_t

  use frozen_external_m, only:     &
    frozen_external_t,             &
    frozen_external_init,          &
    frozen_external_start,         &
    frozen_external_update,        &
    frozen_external_copy,          &
    frozen_external_end

  use base_hamiltonian_m, only:                &
    hamiltonian_setn => base_hamiltonian_setn, &
    hamiltonian_getn => base_hamiltonian_getn, &
    hamiltonian_deln => base_hamiltonian_deln

  use base_hamiltonian_m, only:                         &
    frozen_hamiltonian_t     => base_hamiltonian_t,     &
    frozen_hamiltonian_start => base_hamiltonian_start, &
    frozen_hamiltonian_stop  => base_hamiltonian_stop,  &
    frozen_hamiltonian_get   => base_hamiltonian_get

  implicit none

  private
  public ::                    &
    frozen_hamiltonian_t,      &
    frozen_hamiltonian_init,   &
    frozen_hamiltonian_start,  &
    frozen_hamiltonian_update, &
    frozen_hamiltonian_stop,   &
    frozen_hamiltonian_get,    &
    frozen_hamiltonian_copy,   &
    frozen_hamiltonian_end
  
  interface frozen_hamiltonian_get
    module procedure frozen_hamiltonian_get_external
  end interface !frozen_hamiltonian_get

  !interface frozen_hamiltonian_get_energy
  !  module procedure frozen_external_get_energy
  !end interface !frozen_hamiltonian_get_energy

contains

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_init(this, sys, config)
    type(frozen_hamiltonian_t), intent(inout) :: this
    type(frozen_system_t),      intent(in)    :: sys
    type(json_object_t),        intent(in)    :: config
    !
    type(json_object_t),     pointer :: cnfg
    type(frozen_external_t), pointer :: epot
    integer                          :: ierr
    !
    PUSH_SUB(frozen_hamiltonian_init)
    nullify(cnfg, epot)
    call json_get(config, "external", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(epot)
    call frozen_external_init(epot, sys, cnfg)
    call hamiltonian_setn(this, "external", epot)
    POP_SUB(frozen_hamiltonian_init)
    return
  end subroutine frozen_hamiltonian_init

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_update(this, that, config)
    type(frozen_hamiltonian_t), intent(inout) :: this
    type(fio_hamiltonian_t),    intent(in)    :: that
    type(json_object_t),        intent(in)    :: config
    !
    type(frozen_external_t), pointer :: mept
    type(fio_external_t),    pointer :: sept
    !
    PUSH_SUB(frozen_hamiltonian_update)
    nullify(mept, sept)
    call frozen_hamiltonian_get(this, mept)
    ASSERT(associated(mept))
    call fio_hamiltonian_get(that, sept)
    ASSERT(associated(sept))
    call frozen_external_update(mept, sept, config)
    nullify(mept, sept)
    POP_SUB(frozen_hamiltonian_update)
    return
  end subroutine frozen_hamiltonian_update

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_external(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(frozen_external_t),   pointer     :: that
    !
    PUSH_SUB(frozen_hamiltonian_get_external)
    call hamiltonian_getn(this, "external", that)
    POP_SUB(frozen_hamiltonian_get_external)
    return
  end subroutine frozen_hamiltonian_get_external

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_copy(this, that)
    type(frozen_hamiltonian_t), intent(inout) :: this
    type(frozen_hamiltonian_t), intent(in)    :: that
    !
    type(frozen_external_t), pointer :: oept, iept
    !
    PUSH_SUB(frozen_hamiltonian_copy)
    nullify(oept, iept)
    call frozen_hamiltonian_end(this)
    call frozen_hamiltonian_get_external(that, iept)
    if(associated(iept))then
      SAFE_ALLOCATE(oept)
      call frozen_external_copy(oept, iept)
      call hamiltonian_setn(this, "external", oept)
    end if
    POP_SUB(frozen_hamiltonian_copy)
    return
  end subroutine frozen_hamiltonian_copy

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_end(this)
    type(frozen_hamiltonian_t), intent(inout) :: this
    !
    type(frozen_external_t), pointer :: epot
    !
    PUSH_SUB(frozen_hamiltonian_end)
    nullify(epot)
    call hamiltonian_deln(this, "external", epot)
    if(associated(epot))then
      call frozen_external_end(epot)
      SAFE_DEALLOCATE_P(epot)
      nullify(epot)
    end if
    POP_SUB(frozen_hamiltonian_end)
    return
  end subroutine frozen_hamiltonian_end

end module frozen_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
