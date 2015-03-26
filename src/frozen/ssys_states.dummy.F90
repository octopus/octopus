#include "global.h"

module ssys_states_m
  use global_m
  use messages_m
  use profiling_m
  use ssys_density_m

  implicit none

  private
  public ::            &
    ssys_states_t,     &
    ssys_states_new,   &
    ssys_states_del,   &
    ssys_states_init,  &
    ssys_states_start, &
    ssys_states_get,   &
    ssys_states_copy

  type :: ssys_states_t
    private
  end type ssys_states_t

  interface ssys_states_init
    module procedure ssys_states_init_copy
  end interface ssys_states_init

  interface ssys_states_get
    module procedure ssys_states_get_density
  end interface ssys_states_get

  interface ssys_states_copy
    module procedure ssys_states_copy_states
  end interface ssys_states_copy

contains

  ! ---------------------------------------------------------
  subroutine ssys_states_new(this, that)
    type(ssys_states_t),  intent(inout) :: this
    type(ssys_states_t), pointer        :: that

    PUSH_SUB(ssys_states_new)

    ASSERT(.false.)

    POP_SUB(ssys_states_new)
  end subroutine ssys_states_new

  ! ---------------------------------------------------------
  subroutine ssys_states_del(this)
    type(ssys_states_t), pointer :: this

    PUSH_SUB(ssys_states_del)

    ASSERT(.false.)

    POP_SUB(ssys_states_del)
  end subroutine ssys_states_del

  ! ---------------------------------------------------------
  recursive subroutine ssys_states_init_copy(this, that)
    type(ssys_states_t), intent(out) :: this
    type(ssys_states_t), intent(in)  :: that

    PUSH_SUB(ssys_states_init_copy)

    ASSERT(.false.)

    POP_SUB(ssys_states_init_copy)
  end subroutine ssys_states_init_copy

  ! ---------------------------------------------------------
  subroutine ssys_states_start(this)
    type(ssys_states_t), intent(inout) :: this

    PUSH_SUB(ssys_states_start)

    ASSERT(.false.)

    POP_SUB(ssys_states_start)
  end subroutine ssys_states_start

  ! ---------------------------------------------------------
  subroutine ssys_states_get_density(this, that)
    type(ssys_states_t),   intent(in) :: this
    type(ssys_density_t), pointer     :: that

    PUSH_SUB(ssys_states_get_density)

    ASSERT(.false.)

    POP_SUB(ssys_states_get_density)
  end subroutine ssys_states_get_density
    
  ! ---------------------------------------------------------
  subroutine ssys_states_copy_states(this, that)
    type(ssys_states_t), intent(inout) :: this
    type(ssys_states_t), intent(in)    :: that

    PUSH_SUB(ssys_states_copy_states)

    ASSERT(.false.)

    POP_SUB(ssys_states_copy_states)
  end subroutine ssys_states_copy_states

end module ssys_states_m

!! Local Variables:
!! mode: f90
!! End:
