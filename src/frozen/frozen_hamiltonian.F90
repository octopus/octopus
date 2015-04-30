#include "global.h"

module frozen_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: json_object_t
  use kinds_m, only: wp

  use simulation_m, only: &
    simulation_t

  use fio_external_m, only: &
    fio_external_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_t

  use fio_hamiltonian_m, only: &
    fio_hamiltonian_get

  use frozen_system_m, only: &
    frozen_system_t

  use frozen_external_m, only: &
    frozen_external__acc__

  use frozen_external_m, only: &
    frozen_external_t

  use base_hamiltonian_m, only:                 &
    frozen_hamiltonian_t => base_hamiltonian_t

  use base_hamiltonian_m, only: &
    base_hamiltonian__get__

  use base_hamiltonian_m, only: &
    base_hamiltonian_set,       &
    base_hamiltonian_get

  implicit none

  private
  public ::               &
    frozen_hamiltonian_t

  public ::                    &
    frozen_hamiltonian__acc__

  public ::                 &
    frozen_hamiltonian_get
  
  interface frozen_hamiltonian_get
    module procedure frozen_hamiltonian_get_info
    module procedure frozen_hamiltonian_get_config
    module procedure frozen_hamiltonian_get_simulation
    module procedure frozen_hamiltonian_get_system
    module procedure frozen_hamiltonian_get_external
  end interface frozen_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian__acc__(this, that, config)
    type(frozen_hamiltonian_t), intent(inout) :: this
    type(fio_hamiltonian_t),    intent(in)    :: that
    type(json_object_t),        intent(in)    :: config

    type(frozen_external_t), pointer :: mept
    type(fio_external_t),    pointer :: sept
    real(kind=wp)                    :: energy, enrg

    PUSH_SUB(frozen_hamiltonian__acc__)

    call frozen_hamiltonian_get(this, energy=energy)
    call frozen_hamiltonian_get(that, energy=enrg)
    call base_hamiltonian_set(this, energy=(energy+enrg))
    nullify(mept, sept)
    call frozen_hamiltonian_get(this, mept)
    ASSERT(associated(mept))
    call fio_hamiltonian_get(that, sept)
    ASSERT(associated(sept))
    call frozen_external__acc__(mept, sept, config)
    nullify(mept, sept)

    POP_SUB(frozen_hamiltonian__acc__)
  end subroutine frozen_hamiltonian__acc__

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_info(this, size, nspin, energy)
    type(frozen_hamiltonian_t), intent(in)  :: this
    integer,          optional, intent(out) :: size
    integer,          optional, intent(out) :: nspin
    real(kind=wp),    optional, intent(out) :: energy

    PUSH_SUB(frozen_hamiltonian_get_info)

    call base_hamiltonian_get(this, size=size, nspin=nspin, energy=energy)

    POP_SUB(base_hamiltonian_get_info)
  end subroutine frozen_hamiltonian_get_info

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_config(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(json_object_t),    pointer     :: that

    PUSH_SUB(frozen_hamiltonian_get_config)

    call base_hamiltonian_get(this, that)

    POP_SUB(frozen_hamiltonian_get_config)
  end subroutine frozen_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_system(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(frozen_system_t),     pointer     :: that

    PUSH_SUB(frozen_hamiltonian_get_system)

    call base_hamiltonian_get(this, that)

    POP_SUB(frozen_hamiltonian_get_system)
  end subroutine frozen_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_simulation(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(simulation_t),        pointer     :: that

    PUSH_SUB(frozen_hamiltonian_get_simulation)

    call base_hamiltonian_get(this, that)

    POP_SUB(frozen_hamiltonian_get_simulation)
  end subroutine frozen_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine frozen_hamiltonian_get_external(this, that)
    type(frozen_hamiltonian_t), intent(in) :: this
    type(frozen_external_t),   pointer     :: that

    PUSH_SUB(frozen_hamiltonian_get_external)

    call base_hamiltonian__get__(this, "external", that)

    POP_SUB(frozen_hamiltonian_get_external)
  end subroutine frozen_hamiltonian_get_external

end module frozen_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
