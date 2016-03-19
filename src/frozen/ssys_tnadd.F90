#include "global.h"

module ssys_tnadd_oct_m

  use base_functional_oct_m
  use base_hamiltonian_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::          &
    ssys_tnadd_calc

contains

  ! ---------------------------------------------------------
  subroutine ssys_tnadd_calc(this)
    type(base_hamiltonian_t), intent(inout) :: this

    type(base_functional_t), pointer :: func

    PUSH_SUB(ssys_tnadd_calc)

    nullify(func)
    call base_hamiltonian__reset__(this)
    call base_hamiltonian_get(this, "total", func)
    ASSERT(associated(func))
    call base_functional_calc(func)
    call base_hamiltonian__acc__(this, func)
    nullify(func)
    call base_hamiltonian_get(this, "live", func)
    ASSERT(associated(func))
    call base_functional_calc(func)
    call base_hamiltonian__sub__(this, func)
    nullify(func)

    POP_SUB(ssys_tnadd_calc)
  end subroutine ssys_tnadd_calc

end module ssys_tnadd_oct_m

!! Local Variables:
!! mode: f90
!! End:
