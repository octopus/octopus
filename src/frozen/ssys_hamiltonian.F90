#include "global.h"

module ssys_hamiltonian_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: json_object_t
  use kinds_m, only: wp

  use simulation_m, only: &
    simulation_t

  use ssys_system_m, only: &
    ssys_system_t

  use ssys_external_m, only: &
    ssys_external_t

  use ssys_ionic_m, only: &
    ssys_ionic_t

  use ssys_tnadd_m, only: &
    ssys_tnadd_t

  use base_hamiltonian_m, only:               &
    ssys_hamiltonian_t => base_hamiltonian_t

  use base_hamiltonian_m, only: &
    base_hamiltonian__update__, &
    base_hamiltonian__reset__,  &
    base_hamiltonian__acc__

  use base_hamiltonian_m, only: &
    base_hamiltonian_get

  use base_hamiltonian_m, only:                         &
    ssys_hamiltonian_init   => base_hamiltonian_init,   &
    ssys_hamiltonian_start  => base_hamiltonian_start,  &
    ssys_hamiltonian_update => base_hamiltonian_update, &
    ssys_hamiltonian_stop   => base_hamiltonian_stop,   &
    ssys_hamiltonian_next   => base_hamiltonian_next,   &
    ssys_hamiltonian_copy   => base_hamiltonian_copy,   &
    ssys_hamiltonian_end    => base_hamiltonian_end

  use base_hamiltonian_m, only:                                 &
    ssys_hamiltonian_iterator_t => base_hamiltonian_iterator_t

  use base_hamiltonian_m, only:                                   &
    SSYS_HAMILTONIAN_OK          => BASE_HAMILTONIAN_OK,          &
    SSYS_HAMILTONIAN_KEY_ERROR   => BASE_HAMILTONIAN_KEY_ERROR,   &
    SSYS_HAMILTONIAN_EMPTY_ERROR => BASE_HAMILTONIAN_EMPTY_ERROR

  implicit none

  private
  public ::             &
    ssys_hamiltonian_t

  public ::                  &
    ssys_hamiltonian_init,   &
    ssys_hamiltonian_start,  &
    ssys_hamiltonian_update, &
    ssys_hamiltonian_stop,   &
    ssys_hamiltonian_next,   &
    ssys_hamiltonian_get,    &
    ssys_hamiltonian_copy,   &
    ssys_hamiltonian_end

  public ::                      &
    ssys_hamiltonian_iterator_t

  public ::                       &
    SSYS_HAMILTONIAN_OK,          &
    SSYS_HAMILTONIAN_KEY_ERROR,   &
    SSYS_HAMILTONIAN_EMPTY_ERROR

  interface ssys_hamiltonian_get
    module procedure ssys_hamiltonian_get_info
    module procedure ssys_hamiltonian_get_config
    module procedure ssys_hamiltonian_get_simulation
    module procedure ssys_hamiltonian_get_system
    module procedure ssys_hamiltonian_get_external
    module procedure ssys_hamiltonian_get_ionic
    module procedure ssys_hamiltonian_get_tnadd
  end interface ssys_hamiltonian_get

contains

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_acc(this)
    type(ssys_hamiltonian_t), intent(inout) :: this

    type(ssys_hamiltonian_iterator_t)  :: iter
    type(ssys_hamiltonian_t),  pointer :: subs
    integer                            :: ierr

    PUSH_SUB(ssys_hamiltonian_acc)

    nullify(subs)
    call base_hamiltonian__reset__(this)
    call ssys_hamiltonian_init(iter, this)
    do
      nullify(subs)
      call ssys_hamiltonian_next(iter, subs, ierr)
      if(ierr/=SSYS_HAMILTONIAN_OK)exit
      call base_hamiltonian__acc__(this, subs)
    end do
    call ssys_hamiltonian_end(iter)
    nullify(subs)
    call base_hamiltonian__update__(this)

    POP_SUB(ssys_hamiltonian_update)
  end subroutine ssys_hamiltonian_acc

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_info(this, size, nspin, energy)
    type(ssys_hamiltonian_t), intent(in)  :: this
    integer,        optional, intent(out) :: size
    integer,        optional, intent(out) :: nspin
    real(kind=wp),  optional, intent(out) :: energy

    PUSH_SUB(ssys_hamiltonian_get_info)

    call base_hamiltonian_get(this, size=size, nspin=nspin, energy=energy)

    POP_SUB(base_hamiltonian_get_info)
  end subroutine ssys_hamiltonian_get_info

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_config(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(json_object_t),     pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_config)

    call base_hamiltonian_get(this, that)

    POP_SUB(ssys_hamiltonian_get_config)
  end subroutine ssys_hamiltonian_get_config

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_system(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(ssys_system_t),     pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_system)

    call base_hamiltonian_get(this, that)

    POP_SUB(ssys_hamiltonian_get_system)
  end subroutine ssys_hamiltonian_get_system

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_simulation(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(simulation_t),      pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_simulation)

    call base_hamiltonian_get(this, that)

    POP_SUB(ssys_hamiltonian_get_simulation)
  end subroutine ssys_hamiltonian_get_simulation

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_external(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(ssys_external_t),   pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_external)

    call base_hamiltonian_get(this, "external", that)

    POP_SUB(ssys_hamiltonian_get_external)
  end subroutine ssys_hamiltonian_get_external

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_ionic(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(ssys_ionic_t),      pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_ionic)

    call base_hamiltonian_get(this, "ionic", that)

    POP_SUB(ssys_hamiltonian_get_ionic)
  end subroutine ssys_hamiltonian_get_ionic

  ! ---------------------------------------------------------
  subroutine ssys_hamiltonian_get_tnadd(this, that)
    type(ssys_hamiltonian_t), intent(in) :: this
    type(ssys_tnadd_t),      pointer     :: that

    PUSH_SUB(ssys_hamiltonian_get_tnadd)

    call base_hamiltonian_get(this, "tnadd", that)

    POP_SUB(ssys_hamiltonian_get_tnadd)
  end subroutine ssys_hamiltonian_get_tnadd

end module ssys_hamiltonian_m

!! Local Variables:
!! mode: f90
!! End:
