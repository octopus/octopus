#include "global.h"

module frozen_external_m

  use base_potential_m
  use basis_m
  use fio_external_m
  use global_m
  use json_m
  use kinds_m
  use mesh_m
  use messages_m
  use profiling_m
  use simulation_m
  use space_m

  use base_potential_m, only:              &
    frozen_external_t => base_potential_t

  use base_potential_m, only:                  &
    frozen_external_get => base_potential_get

  implicit none

  private

  public ::            &
    frozen_external_t

  public ::                 &
    frozen_external__acc__

  public ::              &
    frozen_external_get

contains

  ! ---------------------------------------------------------
  subroutine  frozen_external__acc__intrpl(this, intrpl, config)
    type(frozen_external_t),     intent(inout) :: this
    type(fio_external_intrpl_t), intent(in)    :: intrpl
    type(json_object_t),         intent(in)    :: config

    real(kind=wp), dimension(:),     pointer :: potn
    real(kind=wp), dimension(:), allocatable :: x
    type(simulation_t),              pointer :: sim
    type(space_t),                   pointer :: space
    type(mesh_t),                    pointer :: mesh
    type(basis_t)                            :: basis
    real(kind=wp)                            :: pot
    integer                                  :: indx, np

    PUSH_SUB(frozen_external__acc__intrpl)

    nullify(potn, sim, space, mesh)
    call frozen_external_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, space)
    ASSERT(associated(space))
    call basis_init(basis, space, config)
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call frozen_external_get(this, potn)
    ASSERT(associated(potn))
    call frozen_external_get(this, size=np)
    SAFE_ALLOCATE(x(space%dim))
    do indx = 1, np
      call basis_to_internal(basis, mesh%x(indx,1:space%dim), x)
      call fio_external_eval(intrpl, x, pot)
      potn(indx)=potn(indx)+pot
    end do
    SAFE_DEALLOCATE_A(x)
    nullify(potn, space, mesh)
    call basis_end(basis)

    POP_SUB(frozen_external__acc__intrpl)
  end subroutine frozen_external__acc__intrpl

  ! ---------------------------------------------------------
  subroutine frozen_external__acc__(this, that, config)
    type(frozen_external_t), intent(inout) :: this
    type(frozen_external_t), intent(in)    :: that !> fio
    type(json_object_t),     intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    type(fio_external_intrpl_t)  :: intrp
    integer                      :: type, ierr

    PUSH_SUB(frozen_external__update__)

    nullify(cnfg, list)
    call json_get(config, "type", type, ierr)
    if(ierr==JSON_OK)then
      call fio_external_init(intrp, that, type)
    else
      call fio_external_init(intrp, that)
    end if
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
       nullify(cnfg)
       call json_next(iter, cnfg, ierr)
       if(ierr/=JSON_OK)exit
       call frozen_external__acc__intrpl(this, intrp, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call fio_external_end(intrp)

    POP_SUB(frozen_external__acc__)
  end subroutine frozen_external__acc__

end module frozen_external_m

!! Local Variables:
!! mode: f90
!! End:
