#include "global.h"

module frozen_density_m

  use base_density_m
  use basis_m
  use fio_density_m
  use global_m
  use intrpl_m
  use json_m
  use kinds_m
  use mesh_m
  use messages_m
  use profiling_m
  use simulation_m
  use space_m

  use base_density_m, only:             &
    frozen_density_t => base_density_t

  use base_density_m, only:                 &
    frozen_density_get => base_density_get

  implicit none

  private
  public ::           &
    frozen_density_t

  public ::                &
    frozen_density__acc__

  public ::             &
    frozen_density_get

contains

  ! ---------------------------------------------------------
  subroutine  frozen_density__acc__intrpl(this, intrpl, config)
    type(frozen_density_t),     intent(inout) :: this
    type(fio_density_intrpl_t), intent(in)    :: intrpl
    type(json_object_t),        intent(in)    :: config
    !
    real(kind=wp), dimension(:,:),   pointer :: dnst
    real(kind=wp), dimension(:), allocatable :: x, rho
    type(simulation_t),              pointer :: sim
    type(space_t),                   pointer :: space
    type(mesh_t),                    pointer :: mesh
    type(basis_t)                            :: basis
    integer                                  :: indx, jndx, np, nspin

    PUSH_SUB(frozen_density__acc__intrpl)

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
    ASSERT(nspin<3)
    SAFE_ALLOCATE(x(space%dim))
    SAFE_ALLOCATE(rho(nspin))
    do indx = 1, np
      call basis_to_internal(basis, mesh%x(indx,1:space%dim), x)
      call fio_density_eval(intrpl, x, rho)
      do jndx = 1, nspin
        dnst(indx,jndx)=dnst(indx,jndx)+rho(jndx)
      end do
    end do
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(x)
    nullify(dnst, space, mesh)
    call basis_end(basis)

    POP_SUB(frozen_density__acc__intrpl)
  end subroutine frozen_density__acc__intrpl

  ! ---------------------------------------------------------
  subroutine frozen_density__acc__(this, that, config)
    type(frozen_density_t), intent(inout) :: this
    type(frozen_density_t), intent(in)    :: that !> fio
    type(json_object_t),    intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    type(fio_density_intrpl_t)   :: intrp
    integer                      :: type, ierr

    PUSH_SUB(frozen_density__acc__)

    nullify(cnfg, list)
    call json_get(config, "type", type, ierr)
    if(ierr/=JSON_OK)type=NEAREST
    call fio_density_init(intrp, that, type)
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)>0)
    call json_init(iter, list)
    do
       nullify(cnfg)
       call json_next(iter, cnfg, ierr)
       if(ierr/=JSON_OK)exit
       call frozen_density__acc__intrpl(this, intrp, cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call fio_density_end(intrp)

    POP_SUB(frozen_density__acc__)

  end subroutine frozen_density__acc__

end module frozen_density_m
 
!! Local Variables:
!! mode: f90
!! End:
