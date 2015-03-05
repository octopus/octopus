#include "global.h"

module frozen_handle_m

  use global_m
  use messages_m
  use profiling_m

  use json_m, only: json_object_t
  use mpi_m,  only: mpi_world

  use fio_model_m, only: &
    fio_model_t

  use fio_handle_m, only: &
    fio_handle_t,         &
    fio_handle_init,      &
    fio_handle_start,     &
    fio_handle_update,    &
    fio_handle_stop,      &
    fio_handle_get

  use frozen_model_m, only: &
    frozen_model__update__

  use frozen_model_m, only: &
    frozen_model_t

  use base_handle_m, only: &
    base_handle__update__

  use base_handle_m, only: &
    base_handle_init,      &
    base_handle_next,      &
    base_handle_end

  use base_handle_m, only:                       &
    frozen_handle_start => base_handle__start__, &
    frozen_handle_stop  => base_handle__stop__

  use base_handle_m, only:                   &
    frozen_handle_t    => base_handle_t,     &
    frozen_handle_get  => base_handle_get,   &
    frozen_handle_copy => base_handle_copy,  &
    frozen_handle_end  => base_handle_end

  use base_handle_m, only:  &
    BASE_HANDLE_OK,         &
    base_handle_iterator_t

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

  integer, public, parameter :: HNDL_TYPE_FRZN = 2

contains

  ! ---------------------------------------------------------
  subroutine frozen_handle_init(this, config)
    type(frozen_handle_t), intent(out) :: this
    type(json_object_t),   intent(in)  :: config
    !
    integer :: type
    !
    PUSH_SUB(frozen_handle_init)
    call base_handle_init(this, config, fio_handle_init)
    call frozen_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FRZN)
    POP_SUB(frozen_handle_init)
    return
  end subroutine frozen_handle_init

  ! ---------------------------------------------------------
  subroutine frozen_handle__update__(this, that, config)
    type(frozen_handle_t), intent(inout) :: this
    type(fio_handle_t),    intent(in)    :: that
    type(json_object_t),   intent(in)    :: config
    !
    type(frozen_model_t), pointer :: mmdl
    type(fio_model_t),    pointer :: smdl
    !
    PUSH_SUB(frozen_handle__update__)
    nullify(mmdl, smdl)
    call frozen_handle_get(this, mmdl)
    ASSERT(associated(mmdl))
    call fio_handle_get(that, smdl)
    ASSERT(associated(smdl))
    call frozen_model__update__(mmdl, smdl, config)
    nullify(mmdl, smdl)
    POP_SUB(frozen_handle__update__)
    return
  end subroutine frozen_handle__update__

  ! ---------------------------------------------------------
  subroutine frozen_handle_update(this)
    type(frozen_handle_t), intent(inout) :: this
    !
    type(base_handle_iterator_t) :: iter
    type(fio_handle_t),  pointer :: hndl
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(frozen_handle_update)
    nullify(hndl, cnfg)
    call base_handle_init(iter, this)
    do
      nullify(hndl, cnfg)
      call base_handle_next(iter, cnfg, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call fio_handle_start(hndl, mpi_world)
      call fio_handle_update(hndl)
      call frozen_handle__update__(this, hndl, cnfg)
      call fio_handle_stop(hndl)
    end do
    call base_handle_end(iter)
    call base_handle__update__(this)
    nullify(hndl, cnfg)
    POP_SUB(frozen_handle_update)
    return
  end subroutine frozen_handle_update

end module frozen_handle_m

!! Local Variables:
!! mode: f90
!! End:
