#include "global.h"

module ssys_ionic_m

  use global_m
  use messages_m
  use profiling_m

  use atom_m,    only: atom_t, atom_distance
  use kinds_m,   only: wp
  use species_m, only: species_zval

  use base_geom_m, only: &
    base_geom_t

  use base_geom_m, only:  &
    base_geom_iterator_t

  use base_geom_m, only: &
    base_geom_init,      &
    base_geom_next,      &
    base_geom_end

  use ssys_system_m, only: &
    ssys_system_t

  use ssys_system_m, only: &
    ssys_system_get

  use base_term_m, only:         &
    ssys_ionic_t => base_term_t

  use base_term_m, only: &
    base_term__update__, &
    base_term__reset__

  use base_term_m, only: &
    base_term_set

  use base_term_m, only:                   &
    ssys_ionic_new    => base_term_new,    &
    ssys_ionic_del    => base_term_del,    &
    ssys_ionic_init   => base_term_init,   &
    ssys_ionic_update => base_term_update, &
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
  public ::       &
    ssys_ionic_t

  public ::                 &
    ssys_ionic_new,         &
    ssys_ionic_del,         &
    ssys_ionic_init,        &
    ssys_ionic_update,      &
    ssys_ionic_calc,        &
    ssys_ionic_interaction, &
    ssys_ionic_next,        &
    ssys_ionic_get,         &
    ssys_ionic_copy,        &
    ssys_ionic_end

  public ::                &
    ssys_ionic_iterator_t

  public ::               &
    SSYS_IONIC_NAME_LEN

  public ::                 &
    SSYS_IONIC_OK,          &
    SSYS_IONIC_KEY_ERROR,   &
    SSYS_IONIC_EMPTY_ERROR

  interface ssys_ionic__atom_interaction__
    module procedure ssys_ionic__atom_interaction__energy
    module procedure ssys_ionic__atom_interaction__force
  end interface ssys_ionic__atom_interaction__

  interface ssys_ionic__interaction__
    module procedure ssys_ionic__interaction__energy
    module procedure ssys_ionic__interaction__force
  end interface ssys_ionic__interaction__

  interface ssys_ionic_interaction
    module procedure ssys_ionic_interaction_energy
    module procedure ssys_ionic_interaction_force
  end interface ssys_ionic_interaction

contains

  ! ---------------------------------------------------------
  pure function ssys_ionic__in__(name, list) result(in)
    character(len=*),                         intent(in) :: name
    character(len=*), dimension(:), optional, intent(in) :: list

    logical :: in

    integer :: indx

    in=.false.
    if(present(list))then
      do indx = 1, size(list)
        in=(trim(adjustl(name))==trim(adjustl(list(indx))))
        if(in)exit
      end do
    end if

  end function ssys_ionic__in__

  ! ---------------------------------------------------------
  subroutine ssys_ionic__atom_interaction__energy(iatom, jatom, energy)
    type(atom_t),  intent(in)  :: iatom
    type(atom_t),  intent(in)  :: jatom
    real(kind=wp), intent(out) :: energy

    real(kind=wp) :: rr, zi, zj

    PUSH_SUB(ssys_ionic__atom_interaction__energy)

    zi=species_zval(iatom%species)
    zj=species_zval(jatom%species)
    rr=atom_distance(iatom, jatom)
    ASSERT(rr>1.0e-16_wp)
    energy=zi*zj/rr

    POP_SUB(ssys_ionic__atom_interaction__energy)
  end subroutine ssys_ionic__atom_interaction__energy

  ! ---------------------------------------------------------
  subroutine ssys_ionic__atom_interaction__force(iatom, jatom, force)
    type(atom_t),                intent(in)  :: iatom
    type(atom_t),                intent(in)  :: jatom
    real(kind=wp), dimension(:), intent(out) :: force

    real(kind=wp), dimension(size(force)) :: r
    real(kind=wp)                         :: rr, zi, zj, dd

    PUSH_SUB(ssys_ionic__atom_interaction__force)

    force=0.0_wp
    zi=species_zval(iatom%species)
    zj=species_zval(jatom%species)
    r=jatom%x-iatom%x
    rr=sqrt(sum(r**2))
    ASSERT(rr>1.0e-16_wp)
    dd=zi*zj/rr
    force=(dd/rr**2)*r

    POP_SUB(ssys_ionic__atom_interaction__force)
  end subroutine ssys_ionic__atom_interaction__force

  ! ---------------------------------------------------------
  subroutine ssys_ionic__interaction__energy(this, atom, energy)
    type(ssys_ionic_t), intent(in)  :: this
    type(atom_t),       intent(in)  :: atom
    real(kind=wp),      intent(out) :: energy

    type(ssys_system_t), pointer :: sys
    type(base_geom_t),   pointer :: geom
    type(atom_t),        pointer :: iatom
    type(base_geom_iterator_t)   :: iter
    real(kind=wp)                :: enrg

    PUSH_SUB(ssys_ionic__interaction__energy)

    energy=0.0_wp
    nullify(sys, geom, iatom)
    call ssys_ionic_get(this, sys)
    ASSERT(associated(sys))
    call ssys_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iter, geom)
    do
      nullify(iatom)
      call base_geom_next(iter, iatom)
      if(.not.associated(iatom))exit
      call ssys_ionic__atom_interaction__(iatom, atom, enrg)
      energy=energy+enrg
    end do
    call base_geom_end(iter)
    nullify(sys, geom, iatom)

    POP_SUB(ssys_ionic__interaction__energy)
  end subroutine ssys_ionic__interaction__energy

  ! ---------------------------------------------------------
  subroutine ssys_ionic__interaction__force(this, atom, force)
    type(ssys_ionic_t),          intent(in)  :: this
    type(atom_t),                intent(in)  :: atom
    real(kind=wp), dimension(:), intent(out) :: force

    type(ssys_system_t),          pointer :: sys
    type(base_geom_t),            pointer :: geom
    type(atom_t),                 pointer :: iatom
    type(base_geom_iterator_t)            :: iter
    real(kind=wp), dimension(size(force)) :: frc

    PUSH_SUB(ssys_ionic__interaction__force)

    force=0.0_wp
    nullify(sys, geom, iatom)
    call ssys_ionic_get(this, sys)
    ASSERT(associated(sys))
    call ssys_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iter, geom)
    do
      nullify(iatom)
      call base_geom_next(iter, iatom)
      if(.not.associated(iatom))exit
      call ssys_ionic__atom_interaction__(iatom, atom, frc)
      force=force+frc
    end do
    call base_geom_end(iter)
    nullify(sys, geom, iatom)

    POP_SUB(ssys_ionic__interaction__force)
  end subroutine ssys_ionic__interaction__force

  ! ---------------------------------------------------------
  subroutine ssys_ionic__calc__(this)
    type(ssys_ionic_t), intent(inout) :: this

    type(ssys_system_t), pointer :: sys
    type(base_geom_t),   pointer :: geom
    type(atom_t),        pointer :: iatom, jatom
    type(base_geom_iterator_t)   :: iiter, jiter
    real(kind=wp)                :: energy, enrg

    PUSH_SUB(ssys_ionic__calc__)

    call base_term__reset__(this)
    energy=0.0_wp
    nullify(sys, geom, iatom, jatom)
    call ssys_ionic_get(this, sys)
    ASSERT(associated(sys))
    call ssys_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iiter, geom)
    do
      nullify(iatom)
      call base_geom_next(iiter, iatom)
      if(.not.associated(iatom))exit
      call base_geom_init(jiter, iiter)
      do
        nullify(jatom)
        call base_geom_next(jiter, jatom)
        if(.not.associated(jatom))exit
        call ssys_ionic__atom_interaction__(iatom, jatom, enrg)
        energy=energy+enrg
      end do
      call base_geom_end(jiter)
    end do
    call base_geom_end(iiter)
    nullify(sys, geom, iatom, jatom)
    call base_term_set(this, energy=energy)
    call base_term__update__(this)

    POP_SUB(ssys_ionic__calc__)
  end subroutine ssys_ionic__calc__

  ! ---------------------------------------------------------
  recursive subroutine ssys_ionic_calc(this)
    type(ssys_ionic_t), intent(inout) :: this

    type(ssys_ionic_iterator_t) :: iter
    type(ssys_ionic_t), pointer :: subs
    integer                     :: ierr

    PUSH_SUB(ssys_ionic_calc)

    nullify(subs)
    call ssys_ionic_init(iter, this)
    do
      nullify(subs)
      call ssys_ionic_next(iter, subs, ierr)
      if(ierr/=SSYS_IONIC_OK)exit
      call ssys_ionic_calc(subs)
    end do
    call ssys_ionic_end(iter)
    nullify(subs)
    call ssys_ionic__calc__(this)

    POP_SUB(ssys_ionic_calc)
  end subroutine ssys_ionic_calc

  ! ---------------------------------------------------------
  subroutine ssys_ionic_interaction_energy(this, energy, except)
    type(ssys_ionic_t),                       intent(in)  :: this
    real(kind=wp),                            intent(out) :: energy
    character(len=*), dimension(:), optional, intent(in)  :: except

    type(ssys_ionic_t), pointer :: subs
    real(kind=wp)               :: enrg
    integer                     :: indx

    PUSH_SUB(ssys_ionic_interaction_energy)

    nullify(subs)
    call ssys_ionic_get(this, energy=energy)
    if(present(except))then
      do indx = 1, size(except)
        nullify(subs)
        call ssys_ionic_get(this, trim(adjustl(except(indx))), subs)
        ASSERT(associated(subs))
        call ssys_ionic_get(subs, energy=enrg)
        energy=energy-enrg
      end do
      nullify(subs)
    end if

    POP_SUB(ssys_ionic_interaction_energy)
  end subroutine ssys_ionic_interaction_energy

  ! ---------------------------------------------------------
  subroutine ssys_ionic_interaction_force(this, atom, force, except)
    type(ssys_ionic_t),                       intent(in)  :: this
    type(atom_t),                             intent(in)  :: atom
    real(kind=wp),    dimension(:),           intent(out) :: force
    character(len=*), dimension(:), optional, intent(in)  :: except

    type(ssys_ionic_iterator_t)           :: iter
    type(ssys_ionic_t),           pointer :: subs
    character(len=SSYS_IONIC_NAME_LEN)    :: name
    real(kind=wp), dimension(size(force)) :: frce
    integer                               :: ierr

    PUSH_SUB(ssys_ionic_interaction_force)

    force=0.0_wp
    nullify(subs)
    call ssys_ionic_init(iter, this)
    do
      nullify(subs)
      call ssys_ionic_next(iter, name, subs, ierr)
      if(ierr/=SSYS_IONIC_OK)exit
      if(ssys_ionic__in__(name, except))cycle
      call ssys_ionic__interaction__(subs, atom, frce)
      force=force+frce
    end do
    call ssys_ionic_end(iter)
    nullify(subs)

    POP_SUB(ssys_ionic_interaction_force)
  end subroutine ssys_ionic_interaction_force

end module ssys_ionic_m

!! Local Variables:
!! mode: f90
!! End:
