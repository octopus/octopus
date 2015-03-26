#include "global.h"

module ssys_hamiltonian_m
  use global_m
  use messages_m
  use profiling_m
  use ssys_external_m

  implicit none

  private
  public ::               &
    ssys_hamiltonian_t,   &
    ssys_hamiltonian_get

  type :: ssys_hamiltonian_t
    private
  end type ssys_hamiltonian_t

  interface ssys_hamiltonian_get
    module procedure ssys_hamiltonian_get_external
  end interface ssys_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_external(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(ssys_external_t),   pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_external)

    ASSERT(.false.)

    POP_SUB(ssys_hamiltonian_get_external)
  end subroutine ssys_hamiltonian_get_external

end module ssys_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
