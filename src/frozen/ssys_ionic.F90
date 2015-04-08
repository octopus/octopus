#include "global.h"

module ssys_ionic_m

  use global_m
  use messages_m
  use profiling_m
  
  use atom_m,  only: atom_t
  use kinds_m, only: wp

  use base_term_m, only:                   &
    ssys_ionic_t      => base_term_t,      &
    ssys_ionic_new    => base_term_new,    &
    ssys_ionic_del    => base_term_del,    &
    ssys_ionic_init   => base_term_init,   &
    ssys_ionic_next   => base_term_next,   &
    ssys_ionic_get    => base_term_get,    &
    ssys_ionic_copy   => base_term_copy,   &
    ssys_ionic_end    => base_term_end

  use base_term_m, only:                           &
    ssys_ionic_iterator_t => base_term_iterator_t

  use base_term_m, only:                       &
    SSYS_IONIC_NAME_LEN => BASE_TERM_NAME_LEN

  use base_term_m, only:                             &
    SSYS_IONIC_OK          => BASE_TERM_OK,          &
    SSYS_IONIC_KEY_ERROR   => BASE_TERM_KEY_ERROR,   &
    SSYS_IONIC_EMPTY_ERROR => BASE_TERM_EMPTY_ERROR

  implicit none

  private
  public ::            &
    ssys_ionic_t,      &
    ssys_ionic_new,    &
    ssys_ionic_del,    &
    ssys_ionic_init,   &
    ssys_ionic_update, &
    ssys_ionic_calc,   &
    ssys_ionic_next,   &
    ssys_ionic_get,    &
    ssys_ionic_copy,   &
    ssys_ionic_end

  public ::                &
    ssys_ionic_iterator_t

  public ::               &
    SSYS_IONIC_NAME_LEN

  public ::                 &
    SSYS_IONIC_OK,          &
    SSYS_IONIC_KEY_ERROR,   &
    SSYS_IONIC_EMPTY_ERROR

contains

#if 0

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

#endif

  ! ---------------------------------------------------------
  recursive subroutine ssys_ionic_update(this)
    type(ssys_ionic_t), intent(inout) :: this
    !
    !type(base_ionic_iterator_t) :: iter
    !type(base_ionic_t), pointer :: subs
    !real(kind=wp)               :: tnrg, enrg
    !integer                     :: ierr
    !
    PUSH_SUB(ssys_ionic_update)
#if 0
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
#endif
    POP_SUB(ssys_ionic_update)
    return
  end subroutine ssys_ionic_update

  ! ---------------------------------------------------------
  subroutine ssys_ionic_calc(this, atom, energy, force)
    type(ssys_ionic_t),          intent(in)  :: this
    type(atom_t),                intent(in)  :: atom
    real(kind=wp),               intent(out) :: energy
    real(kind=wp), dimension(:), intent(out) :: force
    !
    !type(base_system_t),          pointer :: sys
    !type(base_geom_t),            pointer :: geom
    !type(atom_t),                 pointer :: iatom
    !type(base_geom_iterator_t)            :: iter
    !real(kind=wp), dimension(size(force)) :: r
    !real(kind=wp)                         :: rr, za, zi, dd
    !
    PUSH_SUB(ssys_ionic_calc)
#if 0
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
#endif
    POP_SUB(ssys_ionic_calc)
    return
  end subroutine ssys_ionic_calc

end module ssys_ionic_m

!! Local Variables:
!! mode: f90
!! End:
