#include "global.h"

module frozen_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: JSON_OK, json_object_t, json_array_t, json_get
  use json_m, only: json_array_iterator_t, json_init, json_next, json_end
  use mpi_m,  only: mpi_world

  use grid_m, only: &
    grid_t

  use fio_external_m, only: &
    fio_external_t

  use fio_model_m, only: &
    fio_model_t,         &
    fio_model_get

  use fio_handle_m, only: &
    fio_handle_t,         &
    fio_handle_init,      &
    fio_handle_start,     &
    fio_handle_update,    &
    fio_handle_stop,      &
    fio_handle_get,       &
    fio_handle_copy,      &
    fio_handle_end

  use simulation_m, only:                &
    frozen_simulation_t => simulation_t

  use frozen_model_m, only: &
    frozen_model_t,         &
    frozen_model_init,      &
    frozen_model_start,     &
    frozen_model_update,    &
    frozen_model_copy,      &
    frozen_model_end

  use base_handle_m, only:               &
    handle_init   => base_handle_init,   &
    handle_update => base_handle_update, &
    handle_copy   => base_handle_copy,   &
    handle_end    => base_handle_end

  use base_handle_m, only:                    &
    frozen_handle_t     => base_handle_t,     &
    frozen_handle_start => base_handle_start, &
    frozen_handle_stop  => base_handle_stop,  &
    frozen_handle_next  => base_handle_next,  &
    frozen_handle_get   => base_handle_get

  use base_handle_m, only:                              &
    FROZEN_HANDLE_OK         => BASE_HANDLE_OK,         &
    frozen_handle_iterator_t => base_handle_iterator_t

  implicit none

  private

  public ::               &
    frozen_handle_t,      &
    frozen_handle_init,   &
    frozen_handle_start,  &
    frozen_handle_update, &
    frozen_handle_stop,   &
    frozen_handle_get,    &
    frozen_handle_copy,   &
    frozen_handle_end

  interface frozen_handle_init
    module procedure frozen_handle_init_handle
    module procedure frozen_handle_iterator_init
  end interface frozen_handle_init

  interface frozen_handle_copy
    module procedure frozen_handle_copy_handle
    module procedure frozen_handle_iterator_copy
  end interface frozen_handle_copy

  interface frozen_handle_end
    module procedure frozen_handle_end_handle
    module procedure frozen_handle_iterator_end
  end interface frozen_handle_end

  integer, public, parameter :: HNDL_TYPE_FRZN = 2

contains

  ! ---------------------------------------------------------
  subroutine frozen_handle_init_handle(this, config)
    type(frozen_handle_t), intent(out) :: this
    type(json_object_t),   intent(in)  :: config
    !
    type(json_array_iterator_t)   :: iter
    type(json_object_t),  pointer :: cnfg
    type(json_array_t),   pointer :: list
    type(frozen_model_t), pointer :: modl
    type(fio_handle_t),   pointer :: hndl
    integer                       :: type, ierr
    !
    PUSH_SUB(frozen_handle_init_handle)
    nullify(cnfg, list, modl, hndl)
    call handle_init(this, config)
    call frozen_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FRZN)
    call frozen_handle_get(this, modl)
    ASSERT(associated(modl))
    call json_get(config, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call frozen_model_init(modl, cnfg)
    nullify(modl)
    call json_get(config, "systems", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
      nullify(cnfg, hndl)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      SAFE_ALLOCATE(hndl)
      call fio_handle_init(hndl, cnfg)
      call handle_init(this, hndl, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list, hndl)
    call handle_init(this)
    POP_SUB(frozen_handle_init_handle)
    return
  end subroutine frozen_handle_init_handle

  ! ---------------------------------------------------------
  subroutine frozen_handle_update(this)
    type(frozen_handle_t), intent(inout) :: this
    !
    type(frozen_handle_iterator_t) :: iter
    type(frozen_model_t),  pointer :: mmdl
    type(fio_handle_t),    pointer :: shnd
    type(fio_model_t),     pointer :: smdl
    type(json_object_t),   pointer :: cnfg
    integer                        :: ierr
    !
    PUSH_SUB(frozen_handle_update)
    nullify(mmdl, shnd, smdl)
    call frozen_handle_get(this, mmdl)
    ASSERT(associated(mmdl))
    call frozen_handle_init(iter, this)
    do
      nullify(cnfg, shnd, smdl)
      call frozen_handle_next(iter, cnfg, shnd, ierr)
      if(ierr/=FROZEN_HANDLE_OK)exit
      call fio_handle_start(shnd, mpi_world)
      call fio_handle_update(shnd)
      call fio_handle_get(shnd, smdl)
      ASSERT(associated(smdl))
      call frozen_model_update(mmdl, smdl, cnfg)
      call fio_handle_stop(shnd)
    end do
    call frozen_handle_end(iter)
    call handle_update(this)
    nullify(mmdl, shnd, smdl)
    POP_SUB(frozen_handle_update)
    return
  end subroutine frozen_handle_update

  ! ---------------------------------------------------------
  subroutine frozen_handle_copy_handle(this, that)
    type(frozen_handle_t), intent(inout) :: this
    type(frozen_handle_t), intent(in)    :: that
    !
    type(frozen_model_t), pointer :: omdl, imdl
    !
    PUSH_SUB(frozen_handle_copy_handle)
    nullify(omdl, imdl)
    call handle_copy(this, that)
    call frozen_handle_get(that, imdl)
    ASSERT(associated(imdl))
    call frozen_handle_get(this, omdl)
    ASSERT(associated(omdl))
    call frozen_model_copy(omdl, imdl)
    nullify(omdl, imdl)
    POP_SUB(frozen_handle_copy_handle)
    return
  end subroutine frozen_handle_copy_handle

  ! ---------------------------------------------------------
  subroutine frozen_handle_end_handle(this)
    type(frozen_handle_t), intent(inout) :: this
    !
    type(frozen_handle_iterator_t) :: iter
    type(frozen_model_t),  pointer :: modl
    type(fio_handle_t),    pointer :: hndl
    integer                        :: ierr
    !
    PUSH_SUB(frozen_handle_end_handle)
    nullify(modl, hndl)
    call frozen_handle_init(iter, this)
    do
      nullify(hndl)
      call frozen_handle_next(iter, hndl, ierr)
      if(ierr/=FROZEN_HANDLE_OK)exit
      call fio_handle_end(hndl)
      SAFE_DEALLOCATE_P(hndl)
    end do
    call frozen_handle_end(iter)
    nullify(hndl)
    call frozen_handle_get(this, modl)
    ASSERT(associated(modl))
    call frozen_model_end(modl)
    nullify(modl)
    call handle_end(this)
    POP_SUB(frozen_handle_end_handle)
    return
  end subroutine frozen_handle_end_handle

  ! ---------------------------------------------------------
  subroutine frozen_handle_iterator_init(this, that)
    type(frozen_handle_iterator_t), intent(out) :: this
    type(frozen_handle_t),  target, intent(in)  :: that
    !
    PUSH_SUB(frozen_handle_iterator_init)
    call handle_init(this, that)
    POP_SUB(frozen_handle_iterator_init)
    return
  end subroutine frozen_handle_iterator_init

  ! ---------------------------------------------------------
  subroutine frozen_handle_iterator_copy(this, that)
    type(frozen_handle_iterator_t), intent(inout) :: this
    type(frozen_handle_iterator_t), intent(in)    :: that
    !
    PUSH_SUB(frozen_handle_iterator_copy)
    call handle_copy(this, that)
    POP_SUB(frozen_handle_iterator_copy)
    return
  end subroutine frozen_handle_iterator_copy

  ! ---------------------------------------------------------
  subroutine frozen_handle_iterator_end(this)
    type(frozen_handle_iterator_t), intent(inout) :: this
    !
    PUSH_SUB(frozen_handle_iterator_end)
    call handle_end(this)
    POP_SUB(frozen_handle_iterator_end)
    return
  end subroutine frozen_handle_iterator_end

end module frozen_handle_m

!! Local Variables:
!! mode: f90
!! End:
