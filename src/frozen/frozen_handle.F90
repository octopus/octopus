#include "global.h"

module frozen_handle_oct_m

  use base_handle_oct_m
  use base_model_oct_m
  use fio_handle_oct_m
  use frozen_model_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::         &
    HNDL_TYPE_FRZN

  public ::              &
    frozen_handle_init,  &
    frozen_handle_start, &
    frozen_handle_stop,  &
    frozen_handle_copy,  &
    frozen_handle_end

  integer, parameter :: HNDL_TYPE_FRZN = 2

contains

  ! ---------------------------------------------------------
  subroutine frozen_handle_init(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(frozen_handle_init)

    call base_handle_init(this, config, fio_handle_init)
    call base_handle_get(this, type)
    ASSERT(type==HNDL_TYPE_FRZN)

    POP_SUB(frozen_handle_init)
  end subroutine frozen_handle_init

  ! ---------------------------------------------------------
  subroutine frozen_handle__start__(this, mpi_group)
    type(base_handle_t), intent(inout) :: this !> frozen
    type(mpi_grp_t),     intent(in)    :: mpi_group

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl !> fio
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(frozen_handle__start__)

    nullify(hndl, cnfg)
    call base_handle__reset__(this)
    call base_handle_init(iter, this)
    do
      nullify(hndl, cnfg)
      call base_handle_next(iter, cnfg, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call fio_handle_start(hndl, mpi_group)
      call frozen_handle__acc__(this, hndl, cnfg)
      call fio_handle_stop(hndl)
    end do
    call base_handle_end(iter)
    call base_handle__update__(this)
    nullify(hndl, cnfg)

    POP_SUB(frozen_handle__start__)
  end subroutine frozen_handle__start__

  ! ---------------------------------------------------------
  subroutine frozen_handle_start(this, grid, mpi_group)
    type(base_handle_t), intent(inout) :: this !> frozen
    type(grid_t),        intent(in)    :: grid
    type(mpi_grp_t),     intent(in)    :: mpi_group

    PUSH_SUB(frozen_handle_start)

    call base_handle__start__(this, grid)
    call frozen_handle__start__(this, mpi_group)

    POP_SUB(frozen_handle_start)
  end subroutine frozen_handle_start

  ! ---------------------------------------------------------
  subroutine frozen_handle__acc__(this, that, config)
    type(base_handle_t), intent(inout) :: this !> frozen
    type(base_handle_t), intent(in)    :: that !> fio
    type(json_object_t), intent(in)    :: config

    type(base_model_t), pointer :: mmdl !> frozen
    type(base_model_t), pointer :: smdl !> fio

    PUSH_SUB(frozen_handle__acc__)

    nullify(mmdl, smdl)
    call base_handle_get(this, mmdl)
    ASSERT(associated(mmdl))
    call base_handle_get(that, smdl)
    ASSERT(associated(smdl))
    call frozen_model__acc__(mmdl, smdl, config)
    nullify(mmdl, smdl)

    POP_SUB(frozen_handle__acc__)
  end subroutine frozen_handle__acc__

  ! ---------------------------------------------------------
  subroutine frozen_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(frozen_handle_stop)

    call base_handle_stop(this)

    POP_SUB(frozen_handle_stop)
  end subroutine frozen_handle_stop

  ! ---------------------------------------------------------
  subroutine frozen_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(frozen_handle_copy)

    call base_handle_copy(this, that)

    POP_SUB(frozen_handle_copy)
  end subroutine frozen_handle_copy

  ! ---------------------------------------------------------
  subroutine frozen_handle__end__(this)
    type(base_handle_t), intent(inout) :: this

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl
    integer                      :: ierr

    PUSH_SUB(frozen_handle__end__)

    call base_handle_init(iter, this)
    do
      nullify(hndl)
      call base_handle_next(iter, hndl, ierr)
      if(ierr/=BASE_HANDLE_OK)exit
      call fio_handle_end(hndl)
    end do
    call base_handle_end(iter)
    nullify(hndl)

    POP_SUB(frozen_handle__end__)
  end subroutine frozen_handle__end__

  ! ---------------------------------------------------------
  subroutine frozen_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(frozen_handle_end)

    call frozen_handle__end__(this)
    call base_handle_end(this)

    POP_SUB(frozen_handle_end)
  end subroutine frozen_handle_end

end module frozen_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
