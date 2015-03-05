#include "global.h"

module fio_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_get
  use mpi_m,  only: mpi_grp_t

  use base_handle_m, only: &
    base_handle__init__

  use base_handle_m, only:                  &
    fio_handle_copy => base_handle__copy__, &
    fio_handle_end  => base_handle__end__

  use base_handle_m, only:             &
    fio_handle_t   => base_handle_t,   &
    fio_handle_get => base_handle_get

  use fio_model_m, only: &
    fio_model_t,         &
    fio_model_start,     &
    fio_model_update,    &
    fio_model_stop,      &
    fio_model_end

  implicit none

  private

  public ::            &
    fio_handle_t,      &
    fio_handle_init,   &
    fio_handle_start,  &
    fio_handle_update, &
    fio_handle_stop,   &
    fio_handle_get,    &
    fio_handle_copy,   &
    fio_handle_end

  integer, public, parameter :: HNDL_TYPE_FNIO = 1

contains

  ! ---------------------------------------------------------
  subroutine fio_handle_init(this, config)
    type(fio_handle_t),  intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    integer :: type
    !
    PUSH_SUB(fio_handle_init)
    call base_handle__init__(this, config)
    call fio_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FNIO)
    call base_handle__init__(this)
    POP_SUB(fio_handle_init)
    return
  end subroutine fio_handle_init

  ! ---------------------------------------------------------
  subroutine fio_handle_start(this, mpi_grp)
    type(fio_handle_t), intent(inout) :: this
    type(mpi_grp_t),    intent(in)    :: mpi_grp
    !
    type(fio_model_t), pointer :: modl
    !
    PUSH_SUB(fio_handle_start)
    nullify(modl)
    call fio_handle_get(this, modl)
    ASSERT(associated(modl))
    call fio_model_start(modl, mpi_grp)
    nullify(modl)
    POP_SUB(fio_handle_start)
    return
  end subroutine fio_handle_start

  ! ---------------------------------------------------------
  subroutine fio_handle_update(this)
    type(fio_handle_t), intent(inout) :: this
    !
    type(fio_model_t), pointer :: modl
    !
    PUSH_SUB(fio_handle_update)
    nullify(modl)
    call fio_handle_get(this, modl)
    ASSERT(associated(modl))
    call fio_model_update(modl)
    nullify(modl)
    POP_SUB(fio_handle_update)
    return
  end subroutine fio_handle_update

  ! ---------------------------------------------------------
  subroutine fio_handle_stop(this)
    type(fio_handle_t), intent(inout) :: this
    !
    type(fio_model_t), pointer :: modl
    !
    PUSH_SUB(fio_handle_stop)
    nullify(modl)
    call fio_handle_get(this, modl)
    ASSERT(associated(modl))
    call fio_model_stop(modl)
    nullify(modl)
    POP_SUB(fio_handle_stop)
    return
  end subroutine fio_handle_stop

end module fio_handle_m

!! Local Variables:
!! mode: f90
!! End:
