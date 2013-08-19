#include "global.h"
#include "template.h"

module TEMPLATE(ionic_term_m)

  use atom_m,    only: atom_t, atom_distance
  use kinds_m,   only: wp
  use species_m, only: species_zval

  use TEMPLATE(geometry_m), only:       &
    geometry_t !=> TEMPLATE(geometry_t)

  use TEMPLATE(term_m), only:                     &
    term_set_energy => TEMPLATE(term_set_energy)

  use TEMPLATE(term_m), only:                                     &
    TEMPLATE(ionic_term_t)          => TEMPLATE(term_t),          &
    TEMPLATE(ionic_term_init)       => TEMPLATE(term_init),       &
    TEMPLATE(ionic_term_get)        => TEMPLATE(term_get),        &
    TEMPLATE(ionic_term_get_energy) => TEMPLATE(term_get_energy), &
    TEMPLATE(ionic_term_copy)       => TEMPLATE(term_copy),       &
    TEMPLATE(ionic_term_end)        => TEMPLATE(term_end)

  implicit none

  private
  public ::                               &
    TEMPLATE(ionic_term_t),               &
    TEMPLATE(ionic_term_init),            &
    TEMPLATE(ionic_term_update),          &
    TEMPLATE(ionic_term_get),             &
    TEMPLATE(ionic_term_get_energy),      &
    TEMPLATE(ionic_term_get_interaction), &
    TEMPLATE(ionic_term_copy),            &
    TEMPLATE(ionic_term_end)
  
contains
  
  ! ---------------------------------------------------------
  subroutine TEMPLATE(ionic_term_update)(this)
    type(TEMPLATE(ionic_term_t)), intent(inout) :: this
    !
    type(geometry_t), pointer :: geo
    real(kind=wp)             :: energy
    real(kind=wp)             :: rr, zi, zj
    integer                   :: i, j
    !
    nullify(geo)
    call TEMPLATE(ionic_term_get)(this, geo)
    ASSERT(associated(geo))
    energy=0.0_wp
    do i = 1, geo%natoms
      zi=species_zval(geo%atom(i)%spec)
      do j = i+1, geo%natoms
        rr=atom_distance(geo%atom(i), geo%atom(j))
        zj=species_zval(geo%atom(j)%spec)
        energy=energy+zi*zj/rr
      end do
    end do
    nullify(geo)
    call term_set_energy(this, energy)
    return
  end subroutine TEMPLATE(ionic_term_update)

  ! ---------------------------------------------------------
  subroutine TEMPLATE(ionic_term_get_interaction)(this, atom, energy, force)
    type(TEMPLATE(ionic_term_t)),          intent(in)  :: this
    type(atom_t),                          intent(in)  :: atom
    real(kind=wp),               optional, intent(out) :: energy
    real(kind=wp), dimension(:), optional, intent(out) :: force
    !
    type(geometry_t),         pointer :: geo
    real(kind=wp), dimension(MAX_DIM) :: r, f
    real(kind=wp)                     :: za, zi, rr, dd, e
    integer                           :: i
    !
    e=0.0_wp
    f=0.0_wp
    nullify(geo)
    call TEMPLATE(ionic_term_get)(this, geo)
    ASSERT(associated(geo))
    za=species_zval(atom%spec)
    do i = 1, geo%natoms
      zi=species_zval(geo%atom(i)%spec)
      r=atom%x-geo%atom(i)%x
      rr=sqrt(sum(r**2))
      dd=za*zi/rr
      f=f+(dd/rr**2)*r
      e=e+dd
    end do
    nullify(geo)
    if(present(energy))energy=e
    if(present(force))force=f
    return
  end subroutine TEMPLATE(ionic_term_get_interaction)

end module TEMPLATE(ionic_term_m)

!! Local Variables:
!! mode: f90
!! End:
