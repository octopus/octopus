#include "global.h"

module fio_states_oct_m

  use base_density_oct_m
  use base_states_oct_m
  use fio_density_oct_m
  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::             &
    fio_states__init__, &
    fio_states__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_states__init__(this)
    type(base_states_t), intent(inout) :: this

    type(base_density_t), pointer :: dnst

    PUSH_SUB(fio_states__init__)

    nullify(dnst)
    call base_states_get(this, dnst)
    ASSERT(associated(dnst))
    call fio_density__init__(dnst)
    nullify(dnst)

    POP_SUB(fio_states__init__)
  end subroutine fio_states__init__
    
  ! ---------------------------------------------------------
  subroutine fio_states__load__(this)
    type(base_states_t), intent(inout) :: this

    type(base_density_t), pointer :: density

    PUSH_SUB(fio_states__load__)

    nullify(density)
    call base_states_get(this, density)
    ASSERT(associated(density))
    call fio_density__load__(density)

    POP_SUB(fio_states__load__)
  end subroutine fio_states__load__

end module fio_states_oct_m

!! Local Variables:
!! mode: f90
!! End:
