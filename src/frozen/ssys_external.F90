#include "global.h"

module ssys_external_m

  use global_m
  use messages_m
  use profiling_m

  use kinds_m, only: wp

  use storage_m, only: &
    storage_t

  use storage_m, only: &
    storage_integrate

  use ssys_density_m, only: &
    ssys_density_t

  use base_potential_m, only:            &
    ssys_external_t => base_potential_t

  use ssys_density_m, only: &
    ssys_density_get

  use base_potential_m, only: &
    base_potential__update__, &
    base_potential__reset__,  &
    base_potential__acc__

  use base_potential_m, only: &
    base_potential_set

  use base_potential_m, only:                      &
    ssys_external_new    => base_potential_new,    &
    ssys_external_del    => base_potential_del,    &
    ssys_external_init   => base_potential_init,   &
    ssys_external_start  => base_potential_start,  &
    ssys_external_update => base_potential_update, &
    ssys_external_stop   => base_potential_stop,   &
    ssys_external_next   => base_potential_next,   &
    ssys_external_get    => base_potential_get,    &
    ssys_external_copy   => base_potential_copy,   &
    ssys_external_end    => base_potential_end

  use base_potential_m, only:                             &
    ssys_external_iterator_t => base_potential_iterator_t

  use base_potential_m, only:                          &
    SSYS_EXTERNAL_NAME_LEN => BASE_POTENTIAL_NAME_LEN

  use base_potential_m, only:                                &
    SSYS_EXTERNAL_OK          => BASE_POTENTIAL_OK,          &
    SSYS_EXTERNAL_KEY_ERROR   => BASE_POTENTIAL_KEY_ERROR,   &
    SSYS_EXTERNAL_EMPTY_ERROR => BASE_POTENTIAL_EMPTY_ERROR

  implicit none

  private
  public ::          &
    ssys_external_t

  public ::             &
    ssys_external_calc

  public ::               &
    ssys_external_new,    &
    ssys_external_del,    &
    ssys_external_init,   &
    ssys_external_start,  &
    ssys_external_update, &
    ssys_external_stop,   &
    ssys_external_next,   &
    ssys_external_get,    &
    ssys_external_copy,   &
    ssys_external_end

  public ::                   &
    ssys_external_iterator_t

  public ::                 &
    SSYS_EXTERNAL_NAME_LEN

  public ::                    &
    SSYS_EXTERNAL_OK,          &
    SSYS_EXTERNAL_KEY_ERROR,   &
    SSYS_EXTERNAL_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine ssys_external__acc__energy(this, that)
    type(ssys_external_t), intent(inout) :: this
    type(ssys_external_t), intent(in)    :: that

    real(kind=wp) :: oenrg, ienrg

    PUSH_SUB(ssys_external__acc__energy)

    call ssys_external_get(this, energy=oenrg)
    call ssys_external_get(that, energy=ienrg)
    call base_potential_set(this, energy=(oenrg+ienrg))

    POP_SUB(ssys_external__acc__energy)
  end subroutine ssys_external__acc__energy

  ! ---------------------------------------------------------
  subroutine ssys_external__sub__energy(this, that)
    type(ssys_external_t), intent(inout) :: this
    type(ssys_external_t), intent(in)    :: that

    real(kind=wp) :: oenrg, ienrg

    PUSH_SUB(ssys_external__sub__energy)

    call ssys_external_get(this, energy=oenrg)
    call ssys_external_get(that, energy=ienrg)
    call base_potential_set(this, energy=(oenrg-ienrg))

    POP_SUB(ssys_external__sub__energy)
  end subroutine ssys_external__sub__energy

  ! ---------------------------------------------------------
  subroutine ssys_external__calc__energy(this, that)
    type(ssys_external_t), intent(inout) :: this
    type(storage_t),       intent(in)    :: that

    type(ssys_density_t), pointer :: density
    type(storage_t),      pointer :: data
    real(kind=wp)                 :: energy

    PUSH_SUB(ssys_external__calc__energy)

    nullify(density, data)
    call ssys_external_get(this, density)
    ASSERT(associated(density))
    call ssys_density_get(density, data, total=.true.)
    ASSERT(associated(data))
    nullify(density)
    call storage_integrate(that, data, energy)
    nullify(data)
    call base_potential_set(this, energy=energy)

    POP_SUB(ssys_external__calc__energy)
  end subroutine ssys_external__calc__energy

  ! ---------------------------------------------------------
  subroutine ssys_external_calc_energy(this)
    type(ssys_external_t), intent(inout) :: this

    type(ssys_external_iterator_t) :: iter
    type(ssys_external_t), pointer :: subs
    type(storage_t),       pointer :: data
    integer                        :: ierr

    PUSH_SUB(ssys_external_calc_energy)

    nullify(subs, data)
    call base_potential_set(this, energy=0.0_wp)
    call ssys_external_get(this, data)
    ASSERT(associated(data))
    call ssys_external_init(iter, this)
    do
      nullify(subs)
      call ssys_external_next(iter, subs, ierr)
      if(ierr/=SSYS_EXTERNAL_OK)exit
      call ssys_external__calc__energy(subs, data)
      call ssys_external__acc__energy(this, subs)
    end do
    call ssys_external_end(iter)
    nullify(subs, data)

    POP_SUB(ssys_external_calc_energy)
  end subroutine ssys_external_calc_energy

  ! ---------------------------------------------------------
  subroutine ssys_external_calc_potential(this)
    type(ssys_external_t), intent(inout) :: this

    type(ssys_external_iterator_t) :: iter
    type(ssys_external_t), pointer :: subs
    integer                        :: ierr

    PUSH_SUB(ssys_external_calc_potential)

    nullify(subs)
    call base_potential__reset__(this)
    call ssys_external_init(iter, this)
    do
      nullify(subs)
      call ssys_external_next(iter, subs, ierr)
      if(ierr/=SSYS_EXTERNAL_OK)exit
      call base_potential__acc__(this, subs)
    end do
    call ssys_external_end(iter)
    nullify(subs)
    call base_potential__update__(this)

    POP_SUB(ssys_external_calc_potential)
  end subroutine ssys_external_calc_potential

  ! ---------------------------------------------------------
  subroutine ssys_external_calc(this)
    type(ssys_external_t), intent(inout) :: this

    type(ssys_external_t), pointer :: live

    PUSH_SUB(ssys_external_calc)

    nullify(live)
    call ssys_external_calc_potential(this)
    call ssys_external_calc_energy(this)
    call ssys_external_get(this, "live", live)
    ASSERT(associated(live))
    call ssys_external__sub__energy(this, live)
    nullify(live)

    POP_SUB(ssys_external_calc)
  end subroutine ssys_external_calc

end module ssys_external_m

!! Local Variables:
!! mode: f90
!! End:
