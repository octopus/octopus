#include "global.h"

module frozen_model_m

  use global_m
  use messages_m
  use profiling_m

  use grid_m,  only: grid_t
  use json_m,  only: JSON_OK, json_object_t, json_get
  use space_m, only: space_t

  use fio_system_m, only: &
    fio_system_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t,         &
    fio_hamiltonian_init,      &
    fio_hamiltonian_copy,      &
    fio_hamiltonian_end

  use fio_model_m, only: &
    fio_model_t,         &
    fio_model_get

  use frozen_system_m, only: &
    frozen_system_t,         &
    frozen_system_update

  use frozen_hamiltonian_m, only: &
    frozen_hamiltonian_t,         &
    frozen_hamiltonian_init,      &
    frozen_hamiltonian_update,    &
    frozen_hamiltonian_copy,      &
    frozen_hamiltonian_end

  use base_model_m, only:                   &
    frozen_model_t     => base_model_t,     &
    frozen_model_start => base_model_start, &
    frozen_model_get   => base_model_get

  implicit none

  private
  public ::              &
    frozen_model_t,      &
    frozen_model_init,   &
    frozen_model_start,  &
    frozen_model_update, &
    frozen_model_get,    &
    frozen_model_copy,   &
    frozen_model_end

contains

  ! ---------------------------------------------------------
  subroutine frozen_model_init(this, config)
    type(frozen_model_t), intent(inout) :: this
    type(json_object_t),  intent(in)    :: config
    !
    type(json_object_t),        pointer :: cnfg
    type(frozen_hamiltonian_t), pointer :: hml
    type(frozen_system_t),      pointer :: sys
    integer                             :: ierr
    !
    PUSH_SUB(frozen_model_init)
    nullify(cnfg, hml, sys)
    call frozen_model_get(this, hml)
    ASSERT(associated(hml))
    call frozen_model_get(this, sys)
    ASSERT(associated(sys))
    call json_get(config, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_hamiltonian_init(hml, sys, cnfg)
    nullify(cnfg, hml, sys)
    POP_SUB(frozen_model_init)
    return
  end subroutine frozen_model_init

  ! ---------------------------------------------------------
  subroutine frozen_model_update(this, that, config)
    type(frozen_model_t), intent(inout) :: this
    type(fio_model_t),    intent(in)    :: that
    type(json_object_t),  intent(in)    :: config
    !
    type(frozen_system_t),      pointer :: msys
    type(frozen_hamiltonian_t), pointer :: mhml
    type(fio_system_t),         pointer :: ssys
    type(fio_hamiltonian_t),    pointer :: shml
    !
    PUSH_SUB(frozen_model_update)
    nullify(msys, mhml, ssys, shml)
    call frozen_model_get(this, msys)
    ASSERT(associated(msys))
    call fio_model_get(that, ssys)
    ASSERT(associated(ssys))
    call frozen_system_update(msys, ssys, config)
    nullify(msys, ssys)
    call frozen_model_get(this, mhml)
    ASSERT(associated(mhml))
    call fio_model_get(that, shml)
    ASSERT(associated(shml))
    call frozen_hamiltonian_update(mhml, shml, config)
    nullify(mhml, shml)
    POP_SUB(frozen_model_update)
    return
  end subroutine frozen_model_update

  ! ---------------------------------------------------------
  subroutine frozen_model_copy(this, that)
    type(frozen_model_t), intent(inout) :: this
    type(frozen_model_t), intent(in)    :: that
    !
    type(frozen_hamiltonian_t), pointer :: ohml, ihml
    !
    PUSH_SUB(frozen_model_copy)
    nullify(ohml, ihml)
    call frozen_model_get(this, ohml)
    ASSERT(associated(ohml))
    call frozen_model_get(that, ihml)
    ASSERT(associated(ihml))
    call frozen_hamiltonian_copy(ohml, ihml)
    nullify(ohml, ihml)
    POP_SUB(frozen_model_copy)
    return
  end subroutine frozen_model_copy

  ! ---------------------------------------------------------
  subroutine frozen_model_end(this)
    type(frozen_model_t), intent(inout) :: this
    !
    type(frozen_hamiltonian_t), pointer :: hml
    !
    PUSH_SUB(frozen_model_end)
    nullify(hml)
    call frozen_model_get(this, hml)
    if(associated(hml))then
      call frozen_hamiltonian_end(hml)
      nullify(hml)
    end if
    POP_SUB(frozen_model_end)
    return
  end subroutine frozen_model_end

end module frozen_model_m

!! Local Variables:
!! mode: f90
!! End:
