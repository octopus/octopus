#include "global.h"

module fio_handle_oct_m

  use base_density_oct_m
  use base_handle_oct_m
  use base_model_oct_m
  use fio_density_oct_m
  use fio_grid_oct_m
  use fio_model_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use space_oct_m

  implicit none

  private

  public ::         &
    HNDL_TYPE_FNIO

  public ::            &
    fio_handle_init,   &
    fio_handle_start,  &
    fio_handle_stop,   &
    fio_handle_copy,   &
    fio_handle_end

  integer, parameter :: HNDL_TYPE_FNIO = 1

contains

  ! ---------------------------------------------------------
  subroutine fio_handle__init__(this)
    type(base_handle_t), intent(inout) :: this

    type(base_density_t), pointer :: pdns

    PUSH_SUB(fio_handle__init__)

    nullify(pdns)
    call base_handle_get(this, pdns)
    ASSERT(associated(pdns))
    call fio_density__init__(pdns)
    nullify(pdns)

    POP_SUB(fio_handle__init__)
  end subroutine fio_handle__init__

  ! ---------------------------------------------------------
  subroutine fio_handle_init(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(fio_handle_init)

    call base_handle__init__(this, config)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FNIO)
    call fio_handle__init__(this)

    POP_SUB(fio_handle_init)
  end subroutine fio_handle_init

  ! ---------------------------------------------------------
  subroutine fio_handle__load__(this)
    type(base_handle_t), intent(inout) :: this

    type(base_model_t), pointer :: mdl

    PUSH_SUB(fio_handle__load__)

    nullify(mdl)
    call base_handle_get(this, mdl)
    ASSERT(associated(mdl))
    call fio_model__load__(mdl)
    nullify(mdl)

    POP_SUB(fio_handle__load__)
  end subroutine fio_handle__load__

  ! ---------------------------------------------------------
  subroutine fio_handle__start__(this, group)
    type(base_handle_t), intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: group

    PUSH_SUB(fio_handle__start__)

    call base_handle__start__(this, simstr)

    POP_SUB(fio_handle__start__)

  contains

    subroutine simstr(this)
      type(simulation_t), intent(inout) :: this

      PUSH_SUB(fio_handle__start__.simstr)

      call simulation_start(this, grdini)

      POP_SUB(fio_handle__start__.simstr)
    end subroutine simstr

    subroutine grdini(this, geo, space, config)
      type(grid_t),        intent(out) :: this
      type(geometry_t),    intent(in)  :: geo
      type(space_t),       intent(in)  :: space
      type(json_object_t), intent(in)  :: config

      PUSH_SUB(fio_handle__start__.grdini)

      call fio_grid_init(this, geo, space, group, config)

      POP_SUB(fio_handle__start__.grdini)
    end subroutine grdini

  end subroutine fio_handle__start__

  ! ---------------------------------------------------------
  subroutine fio_handle_start(this, group)
    type(base_handle_t), intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: group

    PUSH_SUB(fio_handle_start)

    call fio_handle__start__(this, group)
    call fio_handle__load__(this)

    POP_SUB(fio_handle_start)
  end subroutine fio_handle_start

  ! ---------------------------------------------------------
  subroutine fio_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(fio_handle_stop)

    call base_handle__stop__(this, simstp)

    POP_SUB(fio_handle_stop)

  contains

    subroutine simstp(this)
      type(simulation_t), intent(inout) :: this

      PUSH_SUB(fio_handle_stop.simstp)

      call simulation_stop(this, fio_grid_end)

      POP_SUB(fio_handle_stop.simstp)
    end subroutine simstp

  end subroutine fio_handle_stop

  ! ---------------------------------------------------------
  subroutine fio_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(fio_handle_copy)

    call fio_handle_end(this)
    call base_handle__copy__(this, that)

    POP_SUB(fio_handle_copy)
  end subroutine fio_handle_copy

  ! ---------------------------------------------------------
  subroutine fio_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    integer :: type
    logical :: start

    PUSH_SUB(fio_handle_end)

    call base_handle_get(this, type=type, started=start)
    ASSERT(type==HNDL_TYPE_FNIO)
    if(start) call fio_handle_stop(this)
    call base_handle_end(this)

    POP_SUB(fio_handle_end)
  end subroutine fio_handle_end

end module fio_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
