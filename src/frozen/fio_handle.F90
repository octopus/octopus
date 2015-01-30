#include "global.h"

module fio_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_get
  use mpi_m,  only: mpi_grp_t

  use base_handle_m, only:               &
    handle_init   => base_handle_init,   &
    handle_start  => base_handle_start,  &
    handle_update => base_handle_update, & 
    handle_stop   => base_handle_stop,   &
    handle_copy   => base_handle_copy,   &
    handle_end    => base_handle_end

  use base_handle_m, only:               &
    fio_handle_t    => base_handle_t,    &
    fio_handle_get  => base_handle_get

  use fio_grid_m, only: &
    fio_grid_t

  use fio_external_m, only: &
    fio_external_t

  use fio_model_m, only: &
    fio_model_t,         &
    fio_model_init,      &
    fio_model_start,     &
    fio_model_update,    &
    fio_model_stop,      &
    fio_model_get,       &
    fio_model_copy,      &
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
    type(json_object_t), pointer :: cnfg
    type(fio_model_t),   pointer :: modl
    integer                      :: type, ierr
    !
    PUSH_SUB(fio_handle_init)
    nullify(cnfg, modl)
    call handle_init(this, config)
    call fio_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FNIO)
    call fio_handle_get(this, modl)
    ASSERT(associated(modl))
    call json_get(config, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_model_init(modl, cnfg)
    nullify(cnfg, modl)
    call handle_init(this)
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
    call handle_start(this)
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
    call handle_update(this)
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
    call handle_stop(this)
    call fio_handle_get(this, modl)
    ASSERT(associated(modl))
    call fio_model_stop(modl)
    nullify(modl)
    POP_SUB(fio_handle_stop)
    return
  end subroutine fio_handle_stop

  ! ---------------------------------------------------------
  subroutine fio_handle_copy(this, that)
    type(fio_handle_t), intent(inout) :: this
    type(fio_handle_t), intent(in)    :: that
    !
    type(fio_model_t), pointer :: omdl, imdl
    !
    PUSH_SUB(fio_handle_copy)
    nullify(omdl, imdl)
    call handle_copy(this, that)
    call fio_handle_get(that, imdl)
    ASSERT(associated(imdl))
    call fio_handle_get(this, omdl)
    ASSERT(associated(omdl))
    call fio_model_copy(omdl, imdl)
    nullify(omdl, imdl)
    POP_SUB(fio_handle_copy)
    return
  end subroutine fio_handle_copy

  ! ---------------------------------------------------------
  subroutine fio_handle_end(this)
    type(fio_handle_t), intent(inout) :: this
    !
    type(fio_model_t), pointer :: modl
    !
    PUSH_SUB(fio_handle_end)
    nullify(modl)
    call fio_handle_get(this, modl)
    if(associated(modl))then
      call fio_model_end(modl)
      nullify(modl)
    end if
    call handle_end(this)
    POP_SUB(fio_handle_end)
    return
  end subroutine fio_handle_end

end module fio_handle_m

!! Local Variables:
!! mode: f90
!! End:
