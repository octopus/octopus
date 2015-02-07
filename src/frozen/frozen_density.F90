#include "global.h"

module frozen_density_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use json_m,  only: json_array_t, json_array_iterator_t, json_init, json_next, json_end
  use kinds_m, only: wp
  use mesh_m,  only: mesh_t
  use space_m, only: space_t

  use basis_m, only: basis_t, basis_init, basis_to_internal, basis_end

  use fio_density_m, only: &
    fio_density_t,         &
    fio_density_init,      &
    fio_density_eval,      &
    fio_density_get,       &
    fio_density_end

  use fio_density_m, only: &
    fio_density_intrpl_t

  use simulation_m, only:                    &
    simulation_t   => simulation_t,   &
    simulation_get => simulation_get

  use base_density_m, only:                     &
    frozen_density_t     => base_density_t,     &
    frozen_density_init  => base_density_init,  &
    frozen_density_start => base_density_start, &
    frozen_density_stop  => base_density_stop,  &
    frozen_density_get   => base_density_get,   &
    frozen_density_copy  => base_density_copy,  &
    frozen_density_end   => base_density_end

  implicit none

  private
  public ::                &
    frozen_density_t,      &
    frozen_density_init,   &
    frozen_density_start,  &
    frozen_density_update, &
    frozen_density_stop,   &
    frozen_density_get,    &
    frozen_density_copy,   &
    frozen_density_end

contains

  ! ---------------------------------------------------------
  subroutine  frozen_density_update_intrpl(this, intrpl, config)
    type(frozen_density_t),     intent(inout) :: this
    type(fio_density_intrpl_t), intent(in)    :: intrpl
    type(json_object_t),        intent(in)    :: config
    !
    real(kind=wp), dimension(:,:),   pointer :: dnst
    real(kind=wp), dimension(:), allocatable :: x, rho
    type(simulation_t),       pointer :: sim
    type(space_t),                   pointer :: space
    type(mesh_t),                    pointer :: mesh
    type(basis_t)                            :: basis
    integer                                  :: indx, np, nspin
    !
    PUSH_SUB(frozen_density_update_intrpl)
    nullify(dnst, sim, space, mesh)
    call frozen_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, space)
    ASSERT(associated(space))
    call basis_init(basis, space, config)
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call frozen_density_get(this, dnst)
    ASSERT(associated(dnst))
    call frozen_density_get(this, size=np, nspin=nspin)
    ASSERT(nspin>0)
    SAFE_ALLOCATE(x(space%dim))
    SAFE_ALLOCATE(rho(nspin))
    do indx = 1, np
      call basis_to_internal(basis, mesh%x(indx,1:space%dim), x)
      call fio_density_eval(intrpl, x, rho)
      dnst(indx,:)=dnst(indx,:)+rho
    end do
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(x)
    nullify(dnst, space, mesh)
    call basis_end(basis)
    POP_SUB(frozen_density_update_intrpl)
    return
  end subroutine frozen_density_update_intrpl

  ! ---------------------------------------------------------
  subroutine frozen_density_update(this, that, config)
    type(frozen_density_t), intent(inout) :: this
    type(fio_density_t),    intent(in)    :: that
    type(json_object_t),    intent(in)    :: config
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    type(fio_density_intrpl_t)   :: intrp
    integer                      :: type, ierr
    !
    PUSH_SUB(frozen_density_update)
    nullify(cnfg, list)
    call json_get(config, "type", type, ierr)
    if(ierr==JSON_OK)then
      call fio_density_init(intrp, that, type)
    else
      call fio_density_init(intrp, that)
    end if
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
       nullify(cnfg)
       call json_next(iter, cnfg, ierr)
       if(ierr/=JSON_OK)exit
       call frozen_density_update_intrpl(this, intrp, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call fio_density_end(intrp)
    POP_SUB(frozen_density_update)
    return
  end subroutine frozen_density_update

end module frozen_density_m
 
!! Local Variables:
!! mode: f90
!! End:
