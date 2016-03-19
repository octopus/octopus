#include "global.h"

module ssys_hamiltonian_oct_m

  use base_hamiltonian_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::               &
    ssys_hamiltonian_acc

contains

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_acc(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_hamiltonian_iterator_t)  :: iter
    type(base_hamiltonian_t),  pointer :: subs
    integer                            :: ierr

    PUSH_SUB(ssys_hamiltonian_acc)

    nullify(subs)
    call base_hamiltonian__reset__(this)
    call base_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call base_hamiltonian_next(iter, subs, ierr)
      if(ierr/=BASE_HAMILTONIAN_OK)exit
      call base_hamiltonian__acc__(this, subs)
    end do
    call base_hamiltonian_end(iter)
    nullify(subs)
    call base_hamiltonian__update__(this)

    POP_SUB(ssys_hamiltonian_acc)
  end subroutine ssys_hamiltonian_acc

end module ssys_hamiltonian_oct_m

!! Local Variables:
!! mode: f90
!! End:
