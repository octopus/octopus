#include "global.h"

module ssys_tnadd_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  use base_hamiltonian_m, only: &
    base_hamiltonian__update__, &
    base_hamiltonian__acc__,    &
    base_hamiltonian__get__

  use base_hamiltonian_m, only:                    &
    ssys_tnadd_start => base_hamiltonian__start__, &
    ssys_tnadd_stop  => base_hamiltonian__stop__

  !use base_hamiltonian_m, only: &
  !  base_hamiltonian_set

  use base_hamiltonian_m, only:               &
    ssys_tnadd_t    => base_hamiltonian_t,    &
    ssys_tnadd_init => base_hamiltonian_init, &
    ssys_tnadd_next => base_hamiltonian_next, &
    ssys_tnadd_get  => base_hamiltonian_get,  &
    ssys_tnadd_copy => base_hamiltonian_copy, &
    ssys_tnadd_end  => base_hamiltonian_end

  use base_hamiltonian_m, only:                           &
    ssys_tnadd_iterator_t => base_hamiltonian_iterator_t

  use base_hamiltonian_m, only:                             &
    SSYS_TNADD_OK          => BASE_HAMILTONIAN_OK,          &
    SSYS_TNADD_KEY_ERROR   => BASE_HAMILTONIAN_KEY_ERROR,   &
    SSYS_TNADD_EMPTY_ERROR => BASE_HAMILTONIAN_EMPTY_ERROR

  use base_functional_m, only: &
    base_functional_t,         &
    base_functional_get

  use ssys_functional_m, only: &
    ssys_functional_t,         &
    ssys_functional_get

  implicit none

  private
  public ::            &
    ssys_tnadd_t,      &
    ssys_tnadd_init,   &
    ssys_tnadd_start,  &
    ssys_tnadd_update, &
    ssys_tnadd_stop,   &
    ssys_tnadd_next,   &
    ssys_tnadd_get,    &
    ssys_tnadd_copy,   &
    ssys_tnadd_end

  public ::                &
    ssys_tnadd_iterator_t

  public ::                 &
    SSYS_TNADD_OK,          &
    SSYS_TNADD_KEY_ERROR,   &
    SSYS_TNADD_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_tnadd_update(this)
    type(ssys_tnadd_t), intent(inout) :: this

    real(kind=wp), dimension(:,:), pointer :: lpot, spot, tpot
    type(base_functional_t),       pointer :: live
    type(ssys_functional_t),       pointer :: ssys
    real(kind=wp)                          :: lenr, senr

    PUSH_SUB(ssys_tnadd_update)

    nullify(live, ssys)
    call base_hamiltonian__update__(this)
    !call base_hamiltonian__get__(this, "live", live)
    ASSERT(associated(live))
    call base_functional_get(live, energy=lenr)
    call base_functional_get(live, lpot)
    ASSERT(associated(lpot))
    nullify(live)
    !call base_hamiltonian__get__(this, "total", ssys)
    ASSERT(associated(ssys))
    call ssys_functional_get(ssys, energy=senr)
    call ssys_functional_get(ssys, spot)
    ASSERT(associated(spot))
    nullify(ssys)
    !call base_hamiltonian_set(this, energy=(senr-lenr))
    !call ssys_tnadd_get(this, tpot)
    tpot=spot-lpot
    nullify(lpot, spot, tpot)

    POP_SUB(ssys_tnadd_update)
  end subroutine ssys_tnadd_update

end module ssys_tnadd_m

!! Local Variables:
!! mode: f90
!! End:
