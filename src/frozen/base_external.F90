#include "global.h"

module base_external_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t, atom_classical_t
  use kinds_m,   only: wp
  use species_m, only: species_zval

  use base_geom_m, only:         &
    geom_t    => base_geom_t,    &
    geom_init => base_geom_init, &
    geom_next => base_geom_next, &
    geom_end  => base_geom_end

  use base_geom_m, only:                     &
    geom_iterator_t => base_geom_iterator_t

  use base_system_m, only:         &
    system_t   => base_system_t,   &
    system_get => base_system_get

  use base_potential_m, only:                        &
    POTENTIAL_INTRPL_OK => BASE_POTENTIAL_INTRPL_OK

  use base_potential_m, only:              &
    potential_eval => base_potential_eval

  use base_potential_m, only:                      &
    base_external_t      => base_potential_t,      &
    base_external_init   => base_potential_init,   &
    base_external_start  => base_potential_start,  &
    base_external_update => base_potential_update, &
    base_external_stop   => base_potential_stop,   &
    base_external_set    => base_potential_set,    &
    base_external_get    => base_potential_get,    &
    base_external_copy   => base_potential_copy,   &
    base_external_end    => base_potential_end

  use base_potential_m, only:                          &
    base_external_intrpl_t => base_potential_intrpl_t

  implicit none

  private
  public ::               &
    base_external_t,      &
    base_external_init,   &
    base_external_start,  &
    base_external_update, &
    base_external_stop,   &
    base_external_eval,   &
    base_external_set,    &
    base_external_get,    &
    base_external_copy,   &
    base_external_end

  public ::                 &
    base_external_intrpl_t

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
    type(system_t),         pointer :: sys
    type(geom_t),           pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(geom_iterator_t)           :: iter
    !
    PUSH_SUB(base_external_classical)
    v=0.0_wp
    nullify(sys, geom)
    call base_external_get(this, sys)
    ASSERT(associated(sys))
    call system_get(sys, geom)
    ASSERT(associated(geom))
    call geom_init(iter, geom)
    do
      nullify(atom)
      call geom_next(iter, atom)
      if(.not.associated(atom))exit
      v=v+base_external_calc(x, atom%x, species_zval(atom%spec))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+base_external_calc(x, catom%x, catom%charge)
    !end do
    call geom_end(iter)
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
    call potential_eval(this, x, v, ierr)
    if(ierr/=POTENTIAL_INTRPL_OK)then
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
