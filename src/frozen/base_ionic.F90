#include "global.h"

module base_ionic_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t, atom_distance
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

  use base_term_m, only:  &
    base_term__update__

  use base_term_m, only: &
    base_term_set

  use base_term_m, only:                     &
    base_ionic__init__ => base_term__init__, &
    base_ionic__add__  => base_term__add__,  &
    base_ionic__copy__ => base_term__copy__, &
    base_ionic__end__  => base_term__end__

  use base_term_m, only:               &
    base_ionic_t    => base_term_t,    &
    base_ionic_new  => base_term_new,  &
    base_ionic_del  => base_term_del,  &
    base_ionic_init => base_term_init, &
    base_ionic_next => base_term_next, &
    base_ionic_get  => base_term_get,  &
    base_ionic_copy => base_term_copy, &
    base_ionic_end  => base_term_end

  use base_term_m, only:                           &
    base_ionic_iterator_t => base_term_iterator_t

  use base_term_m, only:                       &
    BASE_IONIC_NAME_LEN => BASE_TERM_NAME_LEN

  use base_term_m, only:                             &
    BASE_IONIC_OK          => BASE_TERM_OK,          &
    BASE_IONIC_KEY_ERROR   => BASE_TERM_KEY_ERROR,   &
    BASE_IONIC_EMPTY_ERROR => BASE_TERM_EMPTY_ERROR

  implicit none

  private
  public ::               &
    base_ionic__init__,   &
    base_ionic__update__, &
    base_ionic__add__,    &
    base_ionic__copy__,   &
    base_ionic__end__

  public ::            &
    base_ionic_t,      &
    base_ionic_new,    &
    base_ionic_del,    &
    base_ionic_init,   &
    base_ionic_update, &
    base_ionic_calc,   &
    base_ionic_next,   &
    base_ionic_get,    &
    base_ionic_copy,   &
    base_ionic_end

  public ::                &
    base_ionic_iterator_t

  public ::               &
    BASE_IONIC_NAME_LEN

  public ::                 &
    BASE_IONIC_OK,          &
    BASE_IONIC_KEY_ERROR,   &
    BASE_IONIC_EMPTY_ERROR

contains

  ! ---------------------------------------------------------
  subroutine base_ionic__update__(this)
    type(base_ionic_t), intent(inout) :: this
    !
    type(base_system_t), pointer :: sys
    type(base_geom_t),   pointer :: geom
    type(atom_t),        pointer :: iatom, jatom
    type(base_geom_iterator_t)   :: iiter, jiter
    real(kind=wp)                :: rr, zi, zj, energy
    !
    PUSH_SUB(base_ionic__update__)
    call base_term__update__(this)
    energy=0.0_wp
    nullify(sys, geom, iatom, jatom)
    call base_ionic_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iiter, geom)
    do
      nullify(iatom)
      call base_geom_next(iiter, iatom)
      if(.not.associated(iatom))exit
      zi=species_zval(iatom%species)
      call base_geom_init(jiter, iiter)
      do
        nullify(jatom)
        call base_geom_next(jiter, jatom)
        if(.not.associated(jatom))exit
        rr=atom_distance(iatom, jatom)
        ASSERT(rr>1.0e-16_wp)
        zj=species_zval(jatom%species)
        energy=energy+zi*zj/rr
      end do
      call base_geom_end(jiter)
    end do
    call base_geom_end(iiter)
    nullify(sys, geom, iatom, jatom)
    call base_term_set(this, energy=energy)
    POP_SUB(base_ionic__update__)
    return
  end subroutine base_ionic__update__

  ! ---------------------------------------------------------
  recursive subroutine base_ionic_update(this)
    type(base_ionic_t), intent(inout) :: this
    !
    type(base_ionic_iterator_t) :: iter
    type(base_ionic_t), pointer :: subs
    real(kind=wp)               :: tnrg, enrg
    integer                     :: ierr
    !
    PUSH_SUB(base_ionic_update)
    nullify(subs)
    tnrg=0.0_wp
    call base_term_set(this, energy=0.0_wp)
    call base_ionic_init(iter, this)
    do
      nullify(subs)
      call base_ionic_next(iter, subs, ierr)
      if(ierr/=BASE_IONIC_OK)exit
      call base_ionic_update(subs)
      call base_ionic_get(subs, energy=enrg)
      tnrg=tnrg+enrg
    end do
    call base_ionic_end(iter)
    nullify(subs)
    call base_ionic__update__(this)
    call base_ionic_get(this, energy=enrg)
    tnrg=tnrg+enrg
    call base_term_set(this, energy=enrg)
    POP_SUB(base_ionic_update)
    return
  end subroutine base_ionic_update

  ! ---------------------------------------------------------
  subroutine base_ionic_calc(this, atom, energy, force)
    type(base_ionic_t),          intent(in)  :: this
    type(atom_t),                intent(in)  :: atom
    real(kind=wp),               intent(out) :: energy
    real(kind=wp), dimension(:), intent(out) :: force
    !
    type(base_system_t),          pointer :: sys
    type(base_geom_t),            pointer :: geom
    type(atom_t),                 pointer :: iatom
    type(base_geom_iterator_t)            :: iter
    real(kind=wp), dimension(size(force)) :: r
    real(kind=wp)                         :: rr, za, zi, dd
    !
    PUSH_SUB(base_ionic_calc)
    energy=0.0_wp
    force=0.0_wp
    nullify(sys, geom, iatom)
    call base_ionic_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    za=species_zval(atom%species)
    call base_geom_init(iter, geom)
    do
      nullify(iatom)
      call base_geom_next(iter, iatom)
      if(.not.associated(iatom))exit
      zi=species_zval(iatom%species)
      r=atom%x-iatom%x
      rr=sqrt(sum(r**2))
      dd=za*zi/rr
      force=force+(dd/rr**2)*r
      energy=energy+dd
    end do
    call base_geom_end(iter)
    nullify(sys, geom, iatom)
    POP_SUB(base_ionic_calc)
    return
  end subroutine base_ionic_calc

end module base_ionic_m

!! Local Variables:
!! mode: f90
!! End:
