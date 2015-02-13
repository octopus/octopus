#include "global.h"

module frozen_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use simulation_m, only: &
    simulation_t

  use fio_external_m, only: &
    fio_external_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t,         &
    fio_hamiltonian_get

  use frozen_system_m, only: &
    frozen_system_t

  use frozen_external_m, only: &
    frozen_external__update__

  use frozen_external_m, only: &
    frozen_external_t

  use base_hamiltonian_m, only: &
    base_hamiltonian_get

  use base_hamiltonian_m, only:                              &
    frozen_hamiltonian_start  => base_hamiltonian__start__,  &
    frozen_hamiltonian_update => base_hamiltonian__update__, &
    frozen_hamiltonian_stop   => base_hamiltonian__stop__

  use base_hamiltonian_m, only:                       &
    frozen_hamiltonian_t    => base_hamiltonian_t,    &
    frozen_hamiltonian_init => base_hamiltonian_init, &
    frozen_hamiltonian_copy => base_hamiltonian_copy, &
    frozen_hamiltonian_end  => base_hamiltonian_end

  implicit none

  private
  public ::                       &
    frozen_hamiltonian__update__

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
    module procedure frozen_hamiltonian_get_config
    module procedure frozen_hamiltonian_get_simulation
    module procedure frozen_hamiltonian_get_system
    module procedure frozen_hamiltonian_get_external
  end interface frozen_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian__update__(this, that, config)
    type(frozen_hamiltonian_t), intent(inout) :: this
    type(fio_hamiltonian_t),    intent(in)    :: that
    type(json_object_t),        intent(in)    :: config
    !
    type(frozen_external_t), pointer :: mept
    type(fio_external_t),    pointer :: sept
    !
    PUSH_SUB(frozen_hamiltonian__update__)
    nullify(mept, sept)
    call frozen_hamiltonian_get(this, mept)
    ASSERT(associated(mept))
    call fio_hamiltonian_get(that, sept)
    ASSERT(associated(sept))
    call frozen_external__update__(mept, sept, config)
    nullify(mept, sept)
    POP_SUB(frozen_hamiltonian__update__)
    return
  end subroutine frozen_hamiltonian__update__

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_config(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(json_object_t),    pointer     :: that
    !
    PUSH_SUB(frozen_hamiltonian_get_config)
    call base_hamiltonian_get(this, that)
    POP_SUB(frozen_hamiltonian_get_config)
    return
  end subroutine frozen_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_system(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(frozen_system_t),     pointer     :: that
    !
    PUSH_SUB(frozen_hamiltonian_get_system)
    call base_hamiltonian_get(this, that)
    POP_SUB(frozen_hamiltonian_get_system)
    return
  end subroutine frozen_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_simulation(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(simulation_t),     pointer     :: that
    !
    PUSH_SUB(frozen_hamiltonian_get_simulation)
    call base_hamiltonian_get(this, that)
    POP_SUB(frozen_hamiltonian_get_simulation)
    return
  end subroutine frozen_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_external(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(frozen_external_t),   pointer     :: that
    !
    PUSH_SUB(frozen_hamiltonian_get_external)
    call base_hamiltonian_get(this, "external", that)
    POP_SUB(frozen_hamiltonian_get_external)
    return
  end subroutine frozen_hamiltonian_get_external

end module frozen_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
