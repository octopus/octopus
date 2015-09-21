#include "global.h"

module fio_handle_m

  use base_handle_m
  use base_model_m
  use fio_model_m
  use global_m
  use json_m
  use messages_m
  use mpi_m
  use profiling_m

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
  subroutine fio_handle_init(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(fio_handle_init)

    call base_handle__init__(this, config)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FNIO)
    call base_handle__init__(this)

    POP_SUB(fio_handle_init)
  end subroutine fio_handle_init

  ! ---------------------------------------------------------
  subroutine fio_handle_start(this, mpi_grp)
    type(base_handle_t), intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: mpi_grp

    type(base_model_t), pointer :: mdl

    PUSH_SUB(fio_handle_start)

    nullify(mdl)
    call base_handle_get(this, mdl)
    ASSERT(associated(mdl))
    call fio_model_start(mdl, mpi_grp)
    nullify(mdl)
    call fio_handle__load__(this)

    POP_SUB(fio_handle_start)
  end subroutine fio_handle_start

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
  subroutine fio_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    type(base_model_t), pointer :: mdl

    PUSH_SUB(fio_handle_stop)

    nullify(mdl)
    call base_handle_get(this, mdl)
    ASSERT(associated(mdl))
    call fio_model_stop(mdl)
    nullify(mdl)

    POP_SUB(fio_handle_stop)
  end subroutine fio_handle_stop

  ! ---------------------------------------------------------
  subroutine fio_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(fio_handle_copy)

    call base_handle__copy__(this, that)

    POP_SUB(fio_handle_copy)
  end subroutine fio_handle_copy

  ! ---------------------------------------------------------
  subroutine fio_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    type(base_model_t), pointer :: mdl

    PUSH_SUB(fio_handle_end)

    nullify(mdl)
    call base_handle_get(this, mdl)
    ASSERT(associated(mdl))
    call fio_model_end(mdl)
    nullify(mdl)

    POP_SUB(fio_handle_end)
  end subroutine fio_handle_end

end module fio_handle_m

!! Local Variables:
!! mode: f90
!! End:
