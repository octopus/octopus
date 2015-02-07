#include "global.h"

module frozen_simulation_m

  use global_m
  use messages_m
  use profiling_m
  use basis_m
  use fio_m
  use simulation_m
  use simulation_m, only:                        &
    frozen_simulation_init  => simulation_init,  &
    frozen_simulation_start => simulation_start, &
    frozen_simulation_copy  => simulation_copy,  &
    frozen_simulation_end   => simulation_end

  implicit none

  private

  public ::                   &
    frozen_simulation_init,   &
    frozen_simulation_start,  &
    frozen_simulation_extend, &
    frozen_simulation_copy,   &
    frozen_simulation_end
  
contains

  ! ---------------------------------------------------------
  subroutine frozen_simulation_extend(this, that)
    type(simulation_t), intent(inout) :: this
    type(fio_t),     optional, intent(in)    :: that
    !
    type(simulation_t), pointer :: sim
    type(basis_t),      pointer :: basis
    !
    nullify(sim, basis)
    if(present(that))then
      call fio_get(that, sim)
      call fio_get(that, basis)
      call simulation_extend(this, sim, basis)
      nullify(sim, basis)
    else
      call simulation_extend(this)
    end if
    return
  end subroutine frozen_simulation_extend

end module frozen_simulation_m

!! Local Variables:
!! mode: f90
!! End:
