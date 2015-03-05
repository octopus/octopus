#include "global.h"

module live_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use simulation_m, only: &
    simulation_t

  use live_system_m, only: &
    live_system_t

  use live_external_m, only: &
    live_external_t

  use base_hamiltonian_m, only: &
    base_hamiltonian_get

  use base_hamiltonian_m, only:                            &
    live_hamiltonian_init   => base_hamiltonian__init__,   &
    live_hamiltonian_start  => base_hamiltonian__start__,  &
    live_hamiltonian_update => base_hamiltonian__update__, &
    live_hamiltonian_stop   => base_hamiltonian__stop__,   &
    live_hamiltonian_copy   => base_hamiltonian__copy__,   &
    live_hamiltonian_end    => base_hamiltonian__end__

  use base_hamiltonian_m, only:               &
    live_hamiltonian_t => base_hamiltonian_t

  implicit none

  private
  public ::                  &
    live_hamiltonian_t,      &
    live_hamiltonian_init,   &
    live_hamiltonian_start,  &
    live_hamiltonian_update, &
    live_hamiltonian_stop,   &
    live_hamiltonian_get,    &
    live_hamiltonian_copy,   &
    live_hamiltonian_end

  interface live_hamiltonian_get
    module procedure live_hamiltonian_get_config
    module procedure live_hamiltonian_get_simulation
    module procedure live_hamiltonian_get_system
    module procedure live_hamiltonian_get_external
  end interface live_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_config(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(json_object_t),     pointer     :: that
    !
    PUSH_SUB(live_hamiltonian_get_config)
    call base_hamiltonian_get(this, that)
    POP_SUB(live_hamiltonian_get_config)
    return
  end subroutine live_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_system(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(live_system_t),     pointer     :: that
    !
    PUSH_SUB(live_hamiltonian_get_system)
    call base_hamiltonian_get(this, that)
    POP_SUB(live_hamiltonian_get_system)
    return
  end subroutine live_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_simulation(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that
    !
    PUSH_SUB(live_hamiltonian_get_simulation)
    call base_hamiltonian_get(this, that)
    POP_SUB(live_hamiltonian_get_simulation)
    return
  end subroutine live_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_external(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(live_external_t),   pointer     :: that
    !
    PUSH_SUB(live_hamiltonian_get_external)
    call base_hamiltonian_get(this, "external", that)
    POP_SUB(live_hamiltonian_get_external)
    return
  end subroutine live_hamiltonian_get_external

end module live_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
