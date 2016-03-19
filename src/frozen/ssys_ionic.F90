#include "global.h"

module ssys_ionic_oct_m

  use base_geometry_oct_m
  use base_system_oct_m
  use base_term_oct_m
  use atom_oct_m
  use global_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use species_oct_m

  implicit none

  private

  public ::                 &
    ssys_ionic_calc,        &
    ssys_ionic_interaction, &
    ssys_ionic_get

  interface ssys_ionic__atom_interaction__
    module procedure ssys_ionic__atom_interaction__energy
    module procedure ssys_ionic__atom_interaction__force
  end interface ssys_ionic__atom_interaction__

  interface ssys_ionic_get
    module procedure ssys_ionic_get_info
  end interface ssys_ionic_get

contains

  ! ---------------------------------------------------------
  pure function ssys_ionic__in__(name, list) result(in)
    character(len=*),                         intent(in) :: name
    character(len=*), dimension(:), optional, intent(in) :: list

    logical :: in

    integer :: indx

    in = .false.
    if(present(list))then
      do indx = 1, size(list)
        in = (trim(adjustl(name))==trim(adjustl(list(indx))))
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

    zi = species_zval(iatom%species)
    zj = species_zval(jatom%species)
    rr = atom_distance(iatom, jatom)
    ASSERT(rr>1.0e-16_wp)
    energy = zi * zj / rr

    POP_SUB(ssys_ionic__atom_interaction__energy)
  end subroutine ssys_ionic__atom_interaction__energy

  ! ---------------------------------------------------------
  subroutine ssys_ionic__atom_interaction__force(iatom, jatom, force)
    type(atom_t),                intent(in)  :: iatom
    type(atom_t),                intent(in)  :: jatom
    real(kind=wp), dimension(:), intent(out) :: force

    real(kind=wp), dimension(MAX_DIM) :: r
    real(kind=wp)                     :: rr, zi, zj, dd
    integer                           :: nd

    PUSH_SUB(ssys_ionic__atom_interaction__force)

    force = 0.0_wp
    zi = species_zval(iatom%species)
    zj = species_zval(jatom%species)
    r = jatom%x - iatom%x
    rr = sqrt(sum(r**2))
    ASSERT(rr>1.0e-16_wp)
    dd = zi * zj / rr
    nd = size(force)
    ASSERT(nd<=MAX_DIM)
    force = ( dd / rr**2 ) * r(1:nd)

    POP_SUB(ssys_ionic__atom_interaction__force)
  end subroutine ssys_ionic__atom_interaction__force

  ! ---------------------------------------------------------
  subroutine ssys_ionic__calc__(this)
    type(base_term_t), intent(inout) :: this

    type(base_system_t),   pointer :: sys
    type(base_geometry_t), pointer :: geom
    type(atom_t),          pointer :: iatom, jatom
    type(base_geometry_iterator_t) :: iiter, jiter
    real(kind=wp)                  :: energy, enrg

    PUSH_SUB(ssys_ionic__calc__)

    energy = 0.0_wp
    nullify(sys, geom, iatom, jatom)
    call base_term__reset__(this)
    call base_term_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geometry_init(iiter, geom)
    do
      nullify(iatom)
      call base_geometry_next(iiter, iatom)
      if(.not.associated(iatom))exit
      call base_geometry_init(jiter, iiter)
      do
        nullify(jatom)
        call base_geometry_next(jiter, jatom)
        if(.not.associated(jatom))exit
        call ssys_ionic__atom_interaction__(iatom, jatom, enrg)
        energy=energy+enrg
      end do
      call base_geometry_end(jiter)
    end do
    call base_geometry_end(iiter)
    nullify(sys, geom, iatom, jatom)
    call base_term_set(this, energy=energy)
    call base_term__update__(this)

    POP_SUB(ssys_ionic__calc__)
  end subroutine ssys_ionic__calc__

  ! ---------------------------------------------------------
  subroutine ssys_ionic_calc(this)
    type(base_term_t), intent(inout) :: this

    type(base_term_iterator_t) :: iter
    type(base_term_t), pointer :: subs
    integer                    :: ierr

    PUSH_SUB(ssys_ionic_calc)

    nullify(subs)
    call base_term_init(iter, this)
    do
      nullify(subs)
      call base_term_next(iter, subs, ierr)
      if(ierr/=BASE_TERM_OK)exit
      call ssys_ionic__calc__(subs)
    end do
    call base_term_end(iter)
    nullify(subs)
    call ssys_ionic__calc__(this)

    POP_SUB(ssys_ionic_calc)
  end subroutine ssys_ionic_calc

  ! ---------------------------------------------------------
  subroutine ssys_ionic__interaction__(this, atom, force)
    type(base_term_t),           intent(in)  :: this
    type(atom_t),                intent(in)  :: atom
    real(kind=wp), dimension(:), intent(out) :: force

    type(base_system_t),          pointer :: sys
    type(base_geometry_t),        pointer :: geom
    type(atom_t),                 pointer :: iatom
    type(base_geometry_iterator_t)        :: iter
    real(kind=wp), dimension(size(force)) :: frc

    PUSH_SUB(ssys_ionic__interaction__)

    force = 0.0_wp
    nullify(sys, geom, iatom)
    call base_term_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geometry_init(iter, geom)
    do
      nullify(iatom)
      call base_geometry_next(iter, iatom)
      if(.not.associated(iatom))exit
      call ssys_ionic__atom_interaction__(iatom, atom, frc)
      force = force + frc
    end do
    call base_geometry_end(iter)
    nullify(sys, geom, iatom)

    POP_SUB(ssys_ionic__interaction__)
  end subroutine ssys_ionic__interaction__

  ! ---------------------------------------------------------
  subroutine ssys_ionic_interaction(this, atom, force, except)
    type(base_term_t),                        intent(in)  :: this
    type(atom_t),                             intent(in)  :: atom
    real(kind=wp),    dimension(:),           intent(out) :: force
    character(len=*), dimension(:), optional, intent(in)  :: except

    type(base_term_iterator_t)           :: iter
    type(base_term_t),           pointer :: subs
    character(len=BASE_TERM_NAME_LEN)    :: name
    real(kind=wp), dimension(size(force)) :: frce
    integer                               :: ierr

    PUSH_SUB(ssys_ionic_interaction)

    force = 0.0_wp
    nullify(subs)
    call base_term_init(iter, this)
    do
      nullify(subs)
      call base_term_next(iter, name, subs, ierr)
      if(ierr/=BASE_TERM_OK)exit
      if(ssys_ionic__in__(name, except))cycle
      call ssys_ionic__interaction__(subs, atom, frce)
      force = force + frce
    end do
    call base_term_end(iter)
    nullify(subs)

    POP_SUB(ssys_ionic_interaction)
  end subroutine ssys_ionic_interaction

  ! ---------------------------------------------------------
  subroutine ssys_ionic_get_term_by_name(this, name, that)
    type(base_term_t),  intent(in) :: this
    character(len=*),    intent(in) :: name
    type(base_term_t), pointer     :: that

    PUSH_SUB(ssys_ionic_get_term_by_name)

    call base_term_gets(this, name, that)

    POP_SUB(ssys_ionic_get_term_by_name)
  end subroutine ssys_ionic_get_term_by_name

  ! ---------------------------------------------------------
  subroutine ssys_ionic_get_info(this, energy, except)
    type(base_term_t),                        intent(in)  :: this
    real(kind=wp),                            intent(out) :: energy
    character(len=*), dimension(:), optional, intent(in)  :: except

    type(base_term_t), pointer :: subs
    real(kind=wp)              :: enrg
    integer                    :: indx
    
    PUSH_SUB(ssys_ionic_get_info)

    nullify(subs)
    call base_term_get(this, energy=energy)
    if(present(except))then
      do indx = 1, size(except)
        nullify(subs)
        call base_term_gets(this, trim(adjustl(except(indx))), subs)
        ASSERT(associated(subs))
        call base_term_get(subs, energy=enrg)
        energy = energy - enrg
      end do
      nullify(subs)
    end if
    
    POP_SUB(ssys_ionic_get_info)
  end subroutine ssys_ionic_get_info

end module ssys_ionic_oct_m

!! Local Variables:
!! mode: f90
!! End:
