#include "config.h"

module run_prog
use system

implicit none

type(system_type) :: sys
  
contains

subroutine run()
  use global
  use units
  implicit none

  sub_name = 'run'; call push_sub()

  ! initialize some stuff
  call units_init()
  call system_init(sys)

  call pop_sub()
end subroutine run

end module run_prog
