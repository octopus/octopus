#include "global.h"

module fio_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m
  
  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  use fio_simulation_m, only: &
    fio_simulation_t

  use fio_system_m, only: &
    fio_system_t

  use fio_external_m, only: &
    fio_external_t

  use fio_external_m, only: &
    fio_external__load__

  use base_hamiltonian_m, only:              &
    fio_hamiltonian_t => base_hamiltonian_t

  use base_hamiltonian_m, only: &
    base_hamiltonian__start__,  &
    base_hamiltonian__get__

  use base_hamiltonian_m, only: &
    base_hamiltonian_set,       &
    base_hamiltonian_get

  implicit none

  private
  public ::            &
    fio_hamiltonian_t

  public ::                  &
    fio_hamiltonian__load__

  public ::              &
    fio_hamiltonian_get
  
  interface fio_hamiltonian_get
    module procedure fio_hamiltonian_get_config
    module procedure fio_hamiltonian_get_simulation
    module procedure fio_hamiltonian_get_system
    module procedure fio_hamiltonian_get_external
  end interface fio_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian__load__(this)
    type(fio_hamiltonian_t), intent(inout) :: this

    type(fio_external_t), pointer :: epot

    PUSH_SUB(fio_hamiltonian__load__)

    nullify(epot)
    call fio_hamiltonian_get_external(this, epot)
    ASSERT(associated(epot))
    call fio_external__load__(epot)
    nullify(epot)

    POP_SUB(fio_hamiltonian__load__)
  end subroutine fio_hamiltonian__load__

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_config(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(json_object_t),    pointer     :: that

    PUSH_SUB(fio_hamiltonian_get_config)

    call base_hamiltonian_get(this, that)

    POP_SUB(fio_hamiltonian_get_config)
  end subroutine fio_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_system(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_system_t),     pointer     :: that

    PUSH_SUB(fio_hamiltonian_get_system)

    call base_hamiltonian_get(this, that)

    POP_SUB(fio_hamiltonian_get_system)
  end subroutine fio_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_simulation(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_simulation_t), pointer     :: that

    PUSH_SUB(fio_hamiltonian_get_simulation)

    call base_hamiltonian_get(this, that)

    POP_SUB(fio_hamiltonian_get_simulation)
  end subroutine fio_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_external(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_external_t),   pointer     :: that

    PUSH_SUB(fio_hamiltonian_get_external)

    call base_hamiltonian__get__(this, "external", that)

    POP_SUB(fio_hamiltonian_get_external)
  end subroutine fio_hamiltonian_get_external

end module fio_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
