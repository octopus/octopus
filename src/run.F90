#include "config.h"

module run_prog
use hamiltonian

implicit none

type(system_type) :: sys
type(hamiltonian_type) :: h
  
contains

subroutine run()
  use global
  use units
  implicit none

  sub_name = 'run'; call push_sub()

  ! initialize some stuff
  call units_init()
  call system_init(sys)
  call hamiltonian_init(h, sys%st%ispin, sys%m%np)

  call system_end(sys)
  call pop_sub()
end subroutine run

end module run_prog
