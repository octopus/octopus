#include "global.h"

module base_config_m
  use global_m
  use json_m
  use messages_m
  use profiling_m

  implicit none

  private
  public ::            &
    base_config_parse

contains

  ! ---------------------------------------------------------
  subroutine base_config_parse(this, ndim, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim
    integer,   optional, intent(in)  :: nspin

    PUSH_SUB(base_config_parse)

    message(1)="Octopus was compiled without frozen support."
    message(2)="Please --enable-frozen on configure."
    call messages_fatal(2)

    POP_SUB(base_config_parse)
  end subroutine base_config_parse

end module base_config_m

!! Local Variables:
!! mode: f90
!! End:

