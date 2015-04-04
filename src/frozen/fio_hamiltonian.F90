#include "global.h"

module fio_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  
  use json_m, only: json_object_t

  use fio_simulation_m, only: &
    fio_simulation_t

  use fio_system_m, only: &
    fio_system_t

  use fio_external_m, only: &
    fio_external_t,         &
    fio_external_update

  use root_hamiltonian_m, only: &
    root_hamiltonian__get__

  use root_hamiltonian_m, only:                         &
    fio_hamiltonian_init  => root_hamiltonian__init__,  &
    fio_hamiltonian_start => root_hamiltonian__start__, &
    fio_hamiltonian_stop  => root_hamiltonian__stop__,  &
    fio_hamiltonian_copy  => root_hamiltonian__copy__,  &
    fio_hamiltonian_end   => root_hamiltonian__end__

  use base_hamiltonian_m, only: &
    base_hamiltonian_get

  use base_hamiltonian_m, only:              &
    fio_hamiltonian_t => base_hamiltonian_t

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
    module procedure fio_hamiltonian_get_config
    module procedure fio_hamiltonian_get_simulation
    module procedure fio_hamiltonian_get_system
    module procedure fio_hamiltonian_get_external
  end interface fio_hamiltonian_get

contains

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
  subroutine fio_hamiltonian_get_config(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(json_object_t),    pointer     :: that
    !
    PUSH_SUB(fio_hamiltonian_get_config)
    call base_hamiltonian_get(this, that)
    POP_SUB(fio_hamiltonian_get_config)
    return
  end subroutine fio_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_system(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_system_t),     pointer     :: that
    !
    PUSH_SUB(fio_hamiltonian_get_system)
    call base_hamiltonian_get(this, that)
    POP_SUB(fio_hamiltonian_get_system)
    return
  end subroutine fio_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_simulation(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_simulation_t), pointer     :: that
    !
    PUSH_SUB(fio_hamiltonian_get_simulation)
    call base_hamiltonian_get(this, that)
    POP_SUB(fio_hamiltonian_get_simulation)
    return
  end subroutine fio_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_get_external(this, that)
    type(fio_hamiltonian_t), intent(in) :: this
    type(fio_external_t),   pointer     :: that
    !
    PUSH_SUB(fio_hamiltonian_get_external)
    call root_hamiltonian__get__(this, "external", that)
    POP_SUB(fio_hamiltonian_get_external)
    return
  end subroutine fio_hamiltonian_get_external

end module fio_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
