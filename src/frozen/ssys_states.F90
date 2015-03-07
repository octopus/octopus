#include "global.h"

module ssys_states_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  use ssys_density_m, only: &
    ssys_density_t,         &
    ssys_density_update

  use base_states_m, only: &
    base_states_t,         &
    base_states_set,       &
    base_states_get

  use base_states_m, only:                     &
    ssys_states_start => base_states__start__, &
    ssys_states_stop  => base_states__stop__

  use base_states_m, only:                &
    ssys_states_t    => base_states_t,    &
    ssys_states_new  => base_states_new,  &
    ssys_states_del  => base_states_del,  &
    ssys_states_init => base_states_init, &
    ssys_states_next => base_states_next, &
    ssys_states_get  => base_states_get,  &
    ssys_states_copy => base_states_copy, &
    ssys_states_end  => base_states_end

  use base_states_m, only:                            &
    ssys_states_iterator_t => base_states_iterator_t

  use base_states_m, only:                              &
    SSYS_STATES_OK          => BASE_STATES_OK,          &
    SSYS_STATES_KEY_ERROR   => BASE_STATES_KEY_ERROR,   &
    SSYS_STATES_EMPTY_ERROR => BASE_STATES_EMPTY_ERROR

  implicit none

  private
  public ::             &
    ssys_states_t,      &
    ssys_states_new,    &
    ssys_states_del,    &
    ssys_states_init,   &
    ssys_states_start,  &
    ssys_states_update, &
    ssys_states_stop,   &
    ssys_states_next,   &
    ssys_states_get,    &
    ssys_states_copy,   &
    ssys_states_end

  public ::                 &
    ssys_states_iterator_t

  public ::                  &
    SSYS_STATES_OK,          &
    SSYS_STATES_KEY_ERROR,   &
    SSYS_STATES_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_states_update(this)
    type(ssys_states_t), intent(inout) :: this
    !
    type(ssys_states_iterator_t)  :: iter
    type(ssys_density_t), pointer :: dnst
    type(base_states_t),  pointer :: subs
    real(kind=wp)                 :: ochrg, ichrg
    integer                       :: ierr
    !
    PUSH_SUB(ssys_states_update)
    nullify(dnst, subs)
    ochrg=0.0_wp
    ichrg=0.0_wp
    call ssys_states_init(iter, this)
    do
      nullify(subs)
      call ssys_states_next(iter, subs, ierr)
      if(ierr/=SSYS_STATES_OK)exit
      call base_states_get(subs, ichrg)
      ochrg=ochrg+ichrg
    end do
    call ssys_states_end(iter)
    nullify(subs)
    call base_states_set(this, ochrg)
    call ssys_states_get(this, dnst)
    ASSERT(associated(dnst))
    call ssys_density_update(dnst)
    POP_SUB(ssys_states_update)
    return
  end subroutine ssys_states_update

end module ssys_states_m

!! Local Variables:
!! mode: f90
!! End:
