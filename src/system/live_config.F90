#include "global.h"

module live_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use live_handle_oct_m
  use messages_oct_m
  use profiling_oct_m
  use space_oct_m
  use states_oct_m
  use states_dim_oct_m
  use storage_oct_m

  implicit none

  private

  public ::           &
    live_config_parse

contains

  ! ---------------------------------------------------------
  subroutine live_config_parse_density(this, st)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st

    type(json_array_t), pointer :: list
    real(kind=wp)               :: chrg, qtot
    integer                     :: ispn, ik, ierr

    PUSH_SUB(live_config_parse_density)

    nullify(list)
    ASSERT(st%d%nspin>0)
    ASSERT(st%d%nspin<3)
    call json_get(this, "charge", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(st%d%nspin==json_len(list))
    qtot = 0.0_wp
    do ispn = 1, st%d%nspin
      chrg = 0.0_wp
      do ik = 1, st%d%nik
        if(ispn==states_dim_get_spin_index(st%d, ik))&
          chrg = chrg + st%d%kweights(ik) * sum(st%occ(:,ik))
      end do
      call json_set(list, ispn, chrg, ierr)
      ASSERT(ierr==JSON_OK)
      qtot = qtot + chrg
    end do
    ASSERT(.not.abs(1.0_wp-min(st%qtot,qtot)/max(st%qtot,qtot))>epsilon(qtot))
    nullify(list)

    POP_SUB(live_config_parse_density)
  end subroutine live_config_parse_density

  ! ---------------------------------------------------------
  subroutine live_config_parse_states(this, st)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_states)

    nullify(cnfg)
    call json_set(this, "charge", st%qtot)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_density(cnfg, st)
    nullify(cnfg)

    POP_SUB(live_config_parse_states)
  end subroutine live_config_parse_states

  ! ---------------------------------------------------------
  subroutine live_config_parse_system(this, st)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_system)

    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_states(cnfg, st)
    nullify(cnfg)

    POP_SUB(live_config_parse_system)
  end subroutine live_config_parse_system

  ! ---------------------------------------------------------
  subroutine live_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(live_config_parse_external)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_external)
  end subroutine live_config_parse_external

  ! ---------------------------------------------------------
  subroutine live_config_parse_ionic(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(live_config_parse_ionic)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(live_config_parse_ionic)
  end subroutine live_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine live_config_parse_hamiltonian(this)
    type(json_object_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(live_config_parse_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call live_config_parse_external(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call live_config_parse_ionic(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_hamiltonian)
  end subroutine live_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine live_config_parse_model(this, st)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse_model)

    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_system(cnfg, st)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_hamiltonian(cnfg)
    nullify(cnfg)

    POP_SUB(live_config_parse_model)
  end subroutine live_config_parse_model

  ! ---------------------------------------------------------
  subroutine live_config_parse(this, st, space)
    type(json_object_t), intent(out) :: this
    type(states_t),      intent(in)  :: st
    type(space_t),       intent(in)  :: space

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(live_config_parse)

    nullify(cnfg)
    call base_config_parse(this, st%d%nspin, space%dim)
    call json_set(this, "type", HNDL_TYPE_LIVE)
    call json_set(this, "name", "live")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call live_config_parse_model(cnfg, st)
    nullify(cnfg)

    POP_SUB(live_config_parse)
  end subroutine live_config_parse

end module live_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
