#include "global.h"

#undef SUBTEMPLATE_NAME
#undef SUBTEMPLATE_TYPE
#undef TEMPLATE_NAME
#undef TEMPLATE_TYPE
#define TEMPLATE_NAME frozen
#define SUBTEMPLATE_NAME fio
#include "tbase_density.F90"
#undef SUBTEMPLATE_NAME
#undef TEMPLATE_NAME

module frozen_density_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp
  use mesh_m,  only: mesh_t

  use frozen_simulation_m, only: &
    simulation_t

  use frozen_simulation_m, only:             &
    simulation_get => frozen_simulation_get

  use frozen_base_density_m, only:                     &
    base_density_update => frozen_base_density_update

  use frozen_base_density_m, only:                                             &
    frozen_density_t                 => frozen_base_density_t,                 &
    frozen_density_init              => frozen_base_density_init,              &
    frozen_density_start             => frozen_base_density_start,             &
    frozen_density_extend            => frozen_base_density_extend,            &
    frozen_density_get               => frozen_base_density_get,               &
    frozen_density_get_size          => frozen_base_density_get_size,          &
    frozen_density_get_nspin         => frozen_base_density_get_nspin,         &
    frozen_density_get_density       => frozen_base_density_get_density,       &
    frozen_density_get_total_density => frozen_base_density_get_total_density, &
    frozen_density_copy              => frozen_base_density_copy,              &
    frozen_density_end               => frozen_base_density_end

  use frozen_base_density_m, only:                                               &
    frozen_density_interpolation_t    => frozen_base_density_interpolation_t,    &
    frozen_density_interpolation_init => frozen_base_density_interpolation_init, &
    frozen_density_interpolation_eval => frozen_base_density_interpolation_eval, &
    frozen_density_interpolation_copy => frozen_base_density_interpolation_copy, &
    frozen_density_interpolation_end  => frozen_base_density_interpolation_end

  use fio_m, only:      &
    sub_t   => fio_t,   &
    sub_get => fio_get

  use fio_m, only:                                &
    interpolation_t    => fio_interpolation_t,    &
    interpolation_init => fio_interpolation_init, &
    interpolation_eval => fio_interpolation_eval, &
    interpolation_end  => fio_interpolation_end

  use fio_density_m, only:          &
    sub_density_t => fio_density_t

  implicit none

  private
  public ::                           &
    frozen_density_t,                 &
    frozen_density_init,              &
    frozen_density_start,             &
    frozen_density_extend,            &
    frozen_density_update,            &
    frozen_density_get,               &
    frozen_density_get_size,          &
    frozen_density_get_nspin,         &
    frozen_density_get_density,       &
    frozen_density_get_total_density, &
    frozen_density_copy,              &
    frozen_density_end

  public ::                            &
    frozen_density_interpolation_t,    &
    frozen_density_interpolation_init, &
    frozen_density_interpolation_eval, &
    frozen_density_interpolation_copy, &
    frozen_density_interpolation_end

  interface frozen_density_update
    module procedure base_density_update
    module procedure frozen_density_update_build
  end interface frozen_density_update

contains

  ! ---------------------------------------------------------
  subroutine frozen_density_update_build(this, that)
    type(frozen_density_t), intent(inout) :: this
    type(sub_t),            intent(in)    :: that
    !
    real(kind=wp), dimension(:,:),   pointer :: density
    real(kind=wp), dimension(:), allocatable :: rho
    type(simulation_t),              pointer :: sim
    type(mesh_t),                    pointer :: mesh
    type(sub_density_t),             pointer :: srho
    type(interpolation_t)                    :: intrp
    integer                                  :: i
    !
    nullify(density, sim, mesh, srho)
    call frozen_density_get_density(this, density)
    ASSERT(associated(density))
    call frozen_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    SAFE_ALLOCATE(rho(frozen_density_get_nspin(this)))
    call sub_get(that, srho)
    ASSERT(associated(srho))
    call interpolation_init(intrp, that, srho)
    do i = 1, frozen_density_get_size(this)
      call interpolation_eval(intrp, mesh%x(i,:), rho)
      density(i,:)=density(i,:)+rho
    end do
    call interpolation_end(intrp)
    nullify(density, sim, mesh, srho)
    SAFE_DEALLOCATE_A(rho)
    return
  end subroutine frozen_density_update_build

end module frozen_density_m

#define TEMPLATE_NAME frozen
#define SUBTEMPLATE_NAME fio
#include "tgeo.F90"
#include "tstates.F90"
#include "tsystem.F90"
#include "tterm.F90"
#include "tpotential.F90"
#undef SUBTEMPLATE_NAME
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
