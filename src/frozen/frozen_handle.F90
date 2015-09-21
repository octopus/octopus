#include "global.h"

module frozen_handle_m

  use base_handle_m
  use base_model_m
  use fio_handle_m
  use frozen_model_m
  use global_m
  use grid_m
  use json_m
  use messages_m
  use mpi_m
  use profiling_m

  use base_handle_m, only:            &
    frozen_handle_t => base_handle_t

  use base_handle_m, only:                     &
    frozen_handle_stop => base_handle__stop__

  use base_handle_m, only:                   &
    frozen_handle_get  => base_handle_get,   &
    frozen_handle_copy => base_handle_copy,  &
    frozen_handle_end  => base_handle_end

  implicit none

  private

  public ::          &
    frozen_handle_t

  public ::              &
    frozen_handle_init,  &
    frozen_handle_start, &
    frozen_handle_stop,  &
    frozen_handle_get,   &
    frozen_handle_copy,  &
    frozen_handle_end

  integer, public, parameter :: HNDL_TYPE_FRZN = 2

contains

  ! ---------------------------------------------------------
  subroutine frozen_handle_init(this, config)
    type(frozen_handle_t), intent(out) :: this
    type(json_object_t),   intent(in)  :: config

    integer :: type

    PUSH_SUB(frozen_handle_init)

    call base_handle_init(this, config, fio_handle_init)
    call frozen_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FRZN)

    POP_SUB(frozen_handle_init)
  end subroutine frozen_handle_init

  ! ---------------------------------------------------------
  subroutine frozen_handle_start(this, grid)
    type(frozen_handle_t), intent(inout) :: this
    type(grid_t),          intent(in)    :: grid

    type(base_handle_iterator_t)   :: iter
    type(frozen_handle_t), pointer :: hndl !> fio
    type(json_object_t),   pointer :: cnfg
    integer                        :: ierr

    PUSH_SUB(frozen_handle_start)

    nullify(hndl, cnfg)
    call base_handle__start__(this, grid)
    call base_handle__reset__(this)
    call base_handle_init(iter, this)
    do
      nullify(hndl, cnfg)
      call base_handle_next(iter, cnfg, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call fio_handle_start(hndl, mpi_world)
      call frozen_handle__acc__(this, hndl, cnfg)
      call fio_handle_stop(hndl)
    end do
    call frozen_handle_end(iter)
    call base_handle__update__(this)
    nullify(hndl, cnfg)

    POP_SUB(frozen_handle_start)
  end subroutine frozen_handle_start

  ! ---------------------------------------------------------
  subroutine frozen_handle__acc__(this, that, config)
    type(frozen_handle_t), intent(inout) :: this
    type(frozen_handle_t), intent(in)    :: that !> fio
    type(json_object_t),   intent(in)    :: config

    type(base_model_t), pointer :: mmdl
    type(base_model_t), pointer :: smdl

    PUSH_SUB(frozen_handle__acc__)

    nullify(mmdl, smdl)
    call frozen_handle_get(this, mmdl)
    ASSERT(associated(mmdl))
    call frozen_handle_get(that, smdl)
    ASSERT(associated(smdl))
    call frozen_model__acc__(mmdl, smdl, config)
    nullify(mmdl, smdl)

    POP_SUB(frozen_handle__acc__)
  end subroutine frozen_handle__acc__

end module frozen_handle_m

!! Local Variables:
!! mode: f90
!! End:
