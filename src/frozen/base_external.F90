#include "global.h"

module base_external_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t
  use kinds_m,   only: wp
  use species_m, only: species_zval

  use base_geom_m, only: &
    base_geom_t,         &
    base_geom_init,      &
    base_geom_next,      &
    base_geom_end

  use base_geom_m, only: &
    base_geom_iterator_t

  use base_system_m, only: &
    base_system_t,         &
    base_system_get

  use base_potential_m, only: &
    BASE_POTENTIAL_INTRPL_OK

  use base_potential_m, only: &
    base_potential_eval

  use base_potential_m, only:                            &
    base_external__init__   => base_potential__init__,   &
    base_external__start__  => base_potential__start__,  &
    base_external__update__ => base_potential__update__, &
    base_external__stop__   => base_potential__stop__,   &
    base_external__reset__  => base_potential__reset__,  &
    base_external__acc__    => base_potential__acc__,    &
    base_external__add__    => base_potential__add__,    &
    base_external__copy__   => base_potential__copy__,   &
    base_external__end__    => base_potential__end__

  use base_potential_m, only:                      &
    base_external_t      => base_potential_t,      &
    base_external_new    => base_potential_new,    &
    base_external_del    => base_potential_del,    &
    base_external_init   => base_potential_init,   &
    base_external_start  => base_potential_start,  &
    base_external_update => base_potential_update, &
    base_external_stop   => base_potential_stop,   &
    base_external_next   => base_potential_next,   &
    base_external_set    => base_potential_set,    &
    base_external_get    => base_potential_get,    &
    base_external_copy   => base_potential_copy,   &
    base_external_end    => base_potential_end

  use base_potential_m, only:                              &
    base_external_iterator_t => base_potential_iterator_t

  use base_potential_m, only:                          &
    base_external_intrpl_t => base_potential_intrpl_t

  use base_potential_m, only:                          &
    BASE_EXTERNAL_NAME_LEN => BASE_POTENTIAL_NAME_LEN

  use base_potential_m, only:                                &
    BASE_EXTERNAL_OK          => BASE_POTENTIAL_OK,          &
    BASE_EXTERNAL_KEY_ERROR   => BASE_POTENTIAL_KEY_ERROR,   &
    BASE_EXTERNAL_EMPTY_ERROR => BASE_POTENTIAL_EMPTY_ERROR

  implicit none

  private
  public ::                  &
    base_external__init__,   &
    base_external__start__,  &
    base_external__update__, &
    base_external__stop__,   &
    base_external__reset__,  &
    base_external__acc__,    &
    base_external__add__,    &
    base_external__copy__,   &
    base_external__end__

  public ::               &
    base_external_t,      &
    base_external_new,    &
    base_external_del,    &
    base_external_init,   &
    base_external_start,  &
    base_external_update, &
    base_external_stop,   &
    base_external_next,   &
    base_external_eval,   &
    base_external_set,    &
    base_external_get,    &
    base_external_copy,   &
    base_external_end

  public ::                   &
    base_external_iterator_t

  public ::                 &
    base_external_intrpl_t

  public ::                 &
    BASE_EXTERNAL_NAME_LEN

  public ::                    &
    BASE_EXTERNAL_OK,          &
    BASE_EXTERNAL_KEY_ERROR,   &
    BASE_EXTERNAL_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  pure function base_external_calc(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c
    !
    real(kind=wp) :: v
    !
    real(kind=wp) :: r
    !
    r=sqrt(sum((x-y)**2))
    if(r<r_small) r=r_small
    v=-c/r
    return
  end function base_external_calc

  ! ---------------------------------------------------------
  subroutine base_external_classical(this, x, v)
    type(base_external_t),       intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    !
    type(base_system_t),    pointer :: sys
    type(base_geom_t),      pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(base_geom_iterator_t)      :: iter
    !
    PUSH_SUB(base_external_classical)
    v=0.0_wp
    nullify(sys, geom)
    call base_external_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iter, geom)
    do
      nullify(atom)
      call base_geom_next(iter, atom)
      if(.not.associated(atom))exit
      v=v+base_external_calc(x, atom%x, species_zval(atom%species))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call base_geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+base_external_calc(x, catom%x, catom%charge)
    !end do
    call base_geom_end(iter)
    !nullify(catom)
    nullify(sys, geom)
    POP_SUB(base_external_classical)
    return
  end subroutine base_external_classical

  ! ---------------------------------------------------------
  subroutine base_external_eval(this, x, v)
    type(base_external_intrpl_t), intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: v
    !
    type(base_external_t), pointer :: epot
    integer                        :: ierr
    !
    PUSH_SUB(base_external_eval)
    nullify(epot)
    call base_potential_eval(this, x, v, ierr)
    if(ierr/=BASE_POTENTIAL_INTRPL_OK)then
      call base_external_get(this, epot)
      ASSERT(associated(epot))
      call base_external_classical(epot, x, v)
      nullify(epot)
    end if
    POP_SUB(base_external_intrpl_eval)
    return
  end subroutine base_external_eval

end module base_external_m

!! Local Variables:
!! mode: f90
!! End:
