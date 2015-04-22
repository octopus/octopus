#include "global.h"

module live_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t

  use simulation_m, only: &
    simulation_t

  use base_system_m, only: &
    base_system_t

  use base_potential_m, only: &
    base_potential_t

  use base_hamiltonian_m, only:               &
    live_hamiltonian_t => base_hamiltonian_t

  use base_hamiltonian_m, only: &
    base_hamiltonian__get__

  use base_hamiltonian_m, only: &
    base_hamiltonian_get

  implicit none

  private
  public ::             &
    live_hamiltonian_t

  public ::               &
    live_hamiltonian_get

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

    PUSH_SUB(live_hamiltonian_get_config)

    call base_hamiltonian_get(this, that)

    POP_SUB(live_hamiltonian_get_config)
  end subroutine live_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_system(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(base_system_t),     pointer     :: that

    PUSH_SUB(live_hamiltonian_get_system)

    call base_hamiltonian_get(this, that)

    POP_SUB(live_hamiltonian_get_system)
  end subroutine live_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_simulation(this, that)
    type(live_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that

    PUSH_SUB(live_hamiltonian_get_simulation)

    call base_hamiltonian_get(this, that)

    POP_SUB(live_hamiltonian_get_simulation)
  end subroutine live_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine live_hamiltonian_get_external(this, that)
    type(live_hamiltonian_t),  intent(in) :: this
    type(base_potential_t),   pointer     :: that

    PUSH_SUB(live_hamiltonian_get_external)

    call base_hamiltonian__get__(this, "external", that)

    POP_SUB(live_hamiltonian_get_external)
  end subroutine live_hamiltonian_get_external

end module live_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
