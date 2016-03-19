#include "global.h"

module ssys_external_oct_m

  use base_density_oct_m
  use base_potential_oct_m
  use global_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use storage_oct_m

  implicit none

  private

  public ::             &
    ssys_external_calc, &
    ssys_external_get

  interface ssys_external_get
    module procedure ssys_external_get_energy
    module procedure ssys_external_get_potential_1d
    module procedure ssys_external_get_potential_md
  end interface ssys_external_get

contains

  ! ---------------------------------------------------------
  subroutine ssys_external__calc__energy(this)
    type(base_potential_t), intent(inout) :: this

    type(base_density_t), pointer :: density
    type(storage_t),      pointer :: pdat, rdat
    real(kind=wp)                 :: energy

    PUSH_SUB(ssys_external__calc__energy)

    energy = 0.0_wp
    nullify(density, pdat, rdat)
    call base_potential_get(this, pdat)
    ASSERT(associated(pdat))
    call base_potential_get(this, density)
    ASSERT(associated(density))
    call base_density_get(density, rdat, total=.true.)
    ASSERT(associated(rdat))
    nullify(density)
    call storage_integrate(pdat, rdat, energy)
    nullify(pdat, rdat)
    call base_potential_set(this, energy=energy)

    POP_SUB(ssys_external__calc__energy)
  end subroutine ssys_external__calc__energy

  ! ---------------------------------------------------------
  subroutine ssys_external_calc(this)
    type(base_potential_t), intent(inout) :: this

    type(base_potential_iterator_t) :: iter
    type(base_potential_t), pointer :: subs
    integer                        :: ierr

    PUSH_SUB(ssys_external_calc)

    nullify(subs)
    call base_potential__reset__(this)
    call base_potential_init(iter, this)
    do
      nullify(subs)
      call base_potential_next(iter, subs, ierr)
      if(ierr/=BASE_POTENTIAL_OK)exit
      call ssys_external__calc__energy(subs)
      call base_potential__acc__(this, subs)
    end do
    call base_potential_end(iter)
    nullify(subs)
    call ssys_external__calc__energy(this)
    call base_potential__update__(this)

    POP_SUB(ssys_external_calc)
  end subroutine ssys_external_calc

  ! ---------------------------------------------------------
  subroutine ssys_external_get_energy(this, energy, except)
    type(base_potential_t),                   intent(in)  :: this
    real(kind=wp),                            intent(out) :: energy
    character(len=*), dimension(:), optional, intent(in)  :: except
    
    type(base_potential_t), pointer :: subs
    real(kind=wp)                   :: enrg
    integer                         :: indx
    
    PUSH_SUB(ssys_external_get_energy)
    
    nullify(subs)
    call base_potential_get(this, energy=energy)
    if(present(except))then
      do indx = 1, size(except)
        nullify(subs)
        call base_potential_gets(this, trim(adjustl(except(indx))), subs)
        ASSERT(associated(subs))
        call base_potential_get(subs, energy=enrg)
        energy = energy - enrg
      end do
      nullify(subs)
    end if
    
    POP_SUB(ssys_external_get_energy)
  end subroutine ssys_external_get_energy
  
  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_1d(this, that)
    type(base_potential_t),        intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    
    PUSH_SUB(ssys_external_get_potential_1d)
    
    call base_potential_get(this, that)
    
    POP_SUB(ssys_external_get_potential_1d)
  end subroutine ssys_external_get_potential_1d
  
  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_md(this, that)
    type(base_potential_t),          intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    
    PUSH_SUB(ssys_external_get_potential_md)
    
    call base_potential_get(this, that)
    
    POP_SUB(ssys_external_get_potential_md)
  end subroutine ssys_external_get_potential_md

end module ssys_external_oct_m

!! Local Variables:
!! mode: f90
!! End:
