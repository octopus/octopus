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

  use base_density_m, only: &
    base_density_t

  use base_density_m, only: &
    base_density_get

  use base_potential_m, only:            &
    ssys_external_t => base_potential_t

  use base_potential_m, only: &
    base_potential__update__, &
    base_potential__reset__,  &
    base_potential__acc__

  use base_potential_m, only: &
    base_potential_set,       &
    base_potential_get

  use base_potential_m, only:                      &
    ssys_external_new    => base_potential_new,    &
    ssys_external_del    => base_potential_del,    &
    ssys_external_init   => base_potential_init,   &
    ssys_external_start  => base_potential_start,  &
    ssys_external_update => base_potential_update, &
    ssys_external_stop   => base_potential_stop,   &
    ssys_external_next   => base_potential_next,   &
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

  interface ssys_external_get
    module procedure ssys_external_get_potential_by_name
    module procedure ssys_external_get_info
    module procedure ssys_external_get_potential_1d
    module procedure ssys_external_get_potential_md
  end interface ssys_external_get

contains

  ! ---------------------------------------------------------
  subroutine ssys_external__calc__energy(this)
    type(ssys_external_t), intent(inout) :: this

    type(base_density_t), pointer :: density
    type(storage_t),      pointer :: pdat, rdat
    real(kind=wp)                 :: energy

    PUSH_SUB(ssys_external__calc__energy)

    energy=0.0_wp
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
    type(ssys_external_t), intent(inout) :: this

    type(ssys_external_iterator_t) :: iter
    type(ssys_external_t), pointer :: subs
    integer                        :: ierr

    PUSH_SUB(ssys_external_calc)

    nullify(subs)
    call base_potential__reset__(this)
    call ssys_external_init(iter, this)
    do
      nullify(subs)
      call ssys_external_next(iter, subs, ierr)
      if(ierr/=SSYS_EXTERNAL_OK)exit
      call ssys_external__calc__energy(subs)
      call base_potential__acc__(this, subs)
    end do
    call ssys_external_end(iter)
    nullify(subs)
    call ssys_external__calc__energy(this)
    call base_potential__update__(this)

    POP_SUB(ssys_external_calc)
  end subroutine ssys_external_calc

  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_by_name(this, name, that)
    type(ssys_external_t),  intent(in) :: this
    character(len=*),       intent(in) :: name
    type(ssys_external_t), pointer     :: that
    
    PUSH_SUB(ssys_external_get_potential_by_name)
    
    !call base_potential_gets(this, name, that)
    
    POP_SUB(ssys_external_get_potential_by_name)
  end subroutine ssys_external_get_potential_by_name
  
  ! ---------------------------------------------------------
  subroutine ssys_external_get_energy(this, energy, except)
    type(ssys_external_t),                    intent(in)  :: this
    real(kind=wp),                            intent(out) :: energy
    character(len=*), dimension(:), optional, intent(in)  :: except
    
    type(ssys_external_t), pointer :: subs
    real(kind=wp)                  :: enrg
    integer                        :: indx
    
    PUSH_SUB(ssys_external_get_energy)
    
    nullify(subs)
    call base_potential_get(this, energy=energy)
    if(present(except))then
      do indx = 1, size(except)
        nullify(subs)
        !call base_potential_gets(this, trim(adjustl(except(indx))), subs)
        ASSERT(associated(subs))
        call base_potential_get(subs, energy=enrg)
        energy=energy-enrg
      end do
      nullify(subs)
    end if
    
    POP_SUB(ssys_external_get_energy)
  end subroutine ssys_external_get_energy
  
  ! ---------------------------------------------------------
  subroutine ssys_external_get_info(this, size, nspin, energy, except)
    type(ssys_external_t),                    intent(in)  :: this
    integer,                        optional, intent(out) :: size
    integer,                        optional, intent(out) :: nspin
    real(kind=wp),                  optional, intent(out) :: energy
    character(len=*), dimension(:), optional, intent(in)  :: except
    
    PUSH_SUB(ssys_external_get_info)
    
    if(present(energy))&
      call ssys_external_get_energy(this, energy, except)
    call base_potential_get(this, size=size, nspin=nspin)
    
    POP_SUB(ssys_external_get_info)
  end subroutine ssys_external_get_info
	
  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_1d(this, that)
    type(ssys_external_t),        intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    
    PUSH_SUB(ssys_external_get_potential_1d)
    
    call base_potential_get(this, that)
    
    POP_SUB(ssys_external_get_potential_1d)
  end subroutine ssys_external_get_potential_1d
  
  ! ---------------------------------------------------------
  subroutine ssys_external_get_potential_md(this, that)
    type(ssys_external_t),          intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    
    PUSH_SUB(ssys_external_get_potential_md)
    
    call base_potential_get(this, that)
    
    POP_SUB(ssys_external_get_potential_md)
  end subroutine ssys_external_get_potential_md

end module ssys_external_m

!! Local Variables:
!! mode: f90
!! End:
