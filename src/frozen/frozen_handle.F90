#include "global.h"

module frozen_handle_oct_m

  use base_density_oct_m
  use base_geometry_oct_m
  use base_handle_oct_m
  use base_model_oct_m
  use fio_handle_oct_m
  use frozen_geometry_oct_m
  use frozen_model_oct_m
  use global_oct_m
  use grid_oct_m
  use json_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m
  use simulation_oct_m

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
  subroutine frozen_handle__init__(this, config)
    type(base_handle_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: config

    type(base_geometry_t), pointer :: pgeo

    PUSH_SUB(frozen_handle__init__)

    nullify(pgeo)
    call base_handle_get(this, pgeo)
    ASSERT(associated(pgeo))
    call frozen_geometry__init__(pgeo, config)
    nullify(pgeo)

    POP_SUB(frozen_handle__init__)
  end subroutine frozen_handle__init__

  ! ---------------------------------------------------------
  subroutine frozen_handle_init(this, config)
    type(base_handle_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    integer :: type

    PUSH_SUB(frozen_handle_init)

    call base_handle__init__(this, config)
    call base_handle_get(this, type=type)
    ASSERT(type==HNDL_TYPE_FRZN)
    call base_handle__init__(this, fio_handle_init)
    call frozen_handle__init__(this, config)

    POP_SUB(frozen_handle_init)
  end subroutine frozen_handle_init

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
  subroutine frozen_handle__load__(this, group)
    type(base_handle_t), intent(inout) :: this !> frozen
    type(mpi_grp_t),     intent(in)    :: group

    type(base_handle_iterator_t) :: iter
    type(base_handle_t), pointer :: hndl !> fio
    type(json_object_t), pointer :: cnfg

    PUSH_SUB(frozen_handle__load__)

    nullify(hndl, cnfg)
    call base_handle__reset__(this)
    call base_handle_init(iter, this)
    do
      nullify(hndl, cnfg)
      call base_handle_next(iter, hndl)
      if(.not.associated(hndl))exit
      call base_handle_get(hndl, cnfg)
      ASSERT(associated(cnfg))
      call fio_handle_start(hndl, group)
      call frozen_handle__acc__(this, hndl, cnfg)
      call fio_handle_stop(hndl)
    end do
    call base_handle_end(iter)
    call base_handle__update__(this)
    nullify(hndl, cnfg)

    POP_SUB(frozen_handle__load__)
  end subroutine frozen_handle__load__

  ! ---------------------------------------------------------
  subroutine frozen_handle_start(this, sim)
    type(base_handle_t), intent(inout) :: this !> frozen
    type(simulation_t),  intent(in)    :: sim
    
    type(base_density_t), pointer :: pdns
    type(mpi_grp_t),      pointer :: pgrp

    PUSH_SUB(frozen_handle_start)

    nullify(pdns, pgrp)
    call base_handle__start__(this, sim)
    call simulation_get(sim, pgrp)
    ASSERT(associated(pgrp))
    call frozen_handle__load__(this, pgrp)
    nullify(pgrp)
    call base_handle_get(this, pdns)
    ASSERT(associated(pdns))
    call base_density_set(pdns, static=.true.)
    nullify(pdns)

    POP_SUB(frozen_handle_start)
  end subroutine frozen_handle_start

  ! ---------------------------------------------------------
  subroutine frozen_handle_stop(this)
    type(base_handle_t), intent(inout) :: this

    PUSH_SUB(frozen_handle_stop)

    call base_handle__stop__(this)

    POP_SUB(frozen_handle_stop)
  end subroutine frozen_handle_stop

  ! ---------------------------------------------------------
  subroutine frozen_handle_copy(this, that)
    type(base_handle_t), intent(inout) :: this
    type(base_handle_t), intent(in)    :: that

    PUSH_SUB(frozen_handle_copy)

    call frozen_handle_end(this)
    call base_handle_copy(this, that)

    POP_SUB(frozen_handle_copy)
  end subroutine frozen_handle_copy

  ! ---------------------------------------------------------
  subroutine frozen_handle_end(this)
    type(base_handle_t), intent(inout) :: this

    integer :: type

    PUSH_SUB(frozen_handle_end)

    call base_handle_get(this, type=type)
    ASSERT(type==HNDL_TYPE_FRZN)
    call base_handle_end(this, fio_handle_end)

    POP_SUB(frozen_handle_end)
  end subroutine frozen_handle_end

end module frozen_handle_oct_m

!! Local Variables:
!! mode: f90
!! End:
