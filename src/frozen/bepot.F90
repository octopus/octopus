#include "global.h"

module bepot_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t, atom_classical_t
  use json_m,    only: json_object_t
  use kinds_m,   only: wp
  use species_m, only: species_zval

  use bgeom_m, only:             &
    geom_t    => bgeom_t,    &
    geom_init => bgeom_init, &
    geom_next => bgeom_next, &
    geom_end  => bgeom_end

  use bgeom_m, only:                         &
    geom_iterator_t => bgeom_iterator_t

  use bsyst_m, only:         &
    system_t   => bsyst_t,   &
    system_get => bsyst_get

  use bpotn_m, only:                        &
    POTENTIAL_INTRPL_OK => BPOTN_INTRPL_OK, &
    potential_eval      => bpotn_eval,      &
    potential_get       => bpotn_get

  use bpotn_m, only:                      &
    bepot_t          => bpotn_t,          &
    bepot_init       => bpotn_init,       &
    bepot_start      => bpotn_start,      &
    bepot_update     => bpotn_update,     &
    bepot_stop       => bpotn_stop,       &
    bepot_set_energy => bpotn_set_energy, &
    bepot_get        => bpotn_get,        &
    bepot_get_size   => bpotn_get_size,   &
    bepot_get_energy => bpotn_get_energy, &
    bepot_copy       => bpotn_copy,       &
    bepot_end        => bpotn_end

  use bpotn_m, only:                  &
    bepot_intrpl_t => bpotn_intrpl_t

  implicit none

  private
  public ::              &
    bepot_t,             &
    bepot_init,          &
    bepot_start,         &
    bepot_update,        &
    bepot_stop,          &
    bepot_set_energy,    &
    bepot_get,           &
    bepot_get_size,      &
    bepot_get_energy,    &
    bepot_copy,          &
    bepot_end

  public ::                      &
    bepot_intrpl_t

contains

  ! ---------------------------------------------------------
  pure function bepot_calc(x, y, c) result(v)
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
  end function bepot_calc

  ! ---------------------------------------------------------
  subroutine bepot_classical(this, x, v)
    type(bepot_t),               intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    !
    type(system_t),         pointer :: sys
    type(geom_t),           pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(geom_iterator_t)           :: iter
    !
    PUSH_SUB(bepot_classical)
    v=0.0_wp
    nullify(sys, geom)
    call bepot_get(this, sys)
    ASSERT(associated(sys))
    call system_get(sys, geom)
    ASSERT(associated(geom))
    call geom_init(iter, geom)
    do
      nullify(atom)
      call geom_next(iter, atom)
      if(.not.associated(atom))exit
      v=v+bepot_calc(x, atom%x, species_zval(atom%spec))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+bepot_calc(x, catom%x, catom%charge)
    !end do
    call geom_end(iter)
    !nullify(catom)
    nullify(sys, geom)
    POP_SUB(bepot_classical)
    return
  end subroutine bepot_classical

  ! ---------------------------------------------------------
  subroutine bepot_eval(this, x, v)
    type(bepot_intrpl_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v
    !
    type(bepot_t), pointer :: ep
    integer                :: ierr
    !
    PUSH_SUB(bepot_eval)
    nullify(ep)
    call potential_eval(this, x, v, ierr)
    if(ierr/=POTENTIAL_INTRPL_OK)then
      call potential_get(this, ep)
      ASSERT(associated(ep))
      call bepot_classical(ep, x, v)
      nullify(ep)
    end if
    POP_SUB(bepot_intrpl_eval)
    return
  end subroutine bepot_eval

end module bepot_m

!! Local Variables:
!! mode: f90
!! End:
