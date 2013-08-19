#include "global.h"

#undef SUBTEMPLATE_NAME
#undef SUBTEMPLATE_TYPE
#undef TEMPLATE_NAME
#undef TEMPLATE_TYPE
#define TEMPLATE_NAME frozen
#define SUBTEMPLATE_NAME fio
#include "texternal_potential.F90"
#undef SUBTEMPLATE_NAME
#undef TEMPLATE_TYPE
#undef TEMPLATE_NAME

module frozen_external_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp
  use mesh_m,  only: mesh_t

  use frozen_simulation_m, only:             &
    simulation_t   !=> frozen_simulation_t,   &

  use frozen_simulation_m, only:             &
    simulation_get => frozen_simulation_get

  use frozen_external_potential_m, only:                           &
    external_potential_update => frozen_external_potential_update

  use frozen_external_potential_m, only:                                      &
    frozen_external_t             => frozen_external_potential_t,             &
    frozen_external_init          => frozen_external_potential_init,          &
    frozen_external_start         => frozen_external_potential_start,         &
    frozen_external_extend        => frozen_external_potential_extend,        &
    frozen_external_get           => frozen_external_potential_get,           &
    frozen_external_get_size      => frozen_external_potential_get_size,      &
    frozen_external_get_energy    => frozen_external_potential_get_energy,    &
    frozen_external_get_potential => frozen_external_potential_get_potential, &
    frozen_external_copy          => frozen_external_potential_copy,          &
    frozen_external_end           => frozen_external_potential_end

  use frozen_external_potential_m, only:                                                &
    frozen_external_interpolation_t    => frozen_external_potential_interpolation_t,    &
    frozen_external_interpolation_init => frozen_external_potential_interpolation_init, &
    frozen_external_interpolation_eval => frozen_external_potential_interpolation_eval, &
    frozen_external_interpolation_copy => frozen_external_potential_interpolation_copy, &
    frozen_external_interpolation_end  => frozen_external_potential_interpolation_end

  use fio_m, only:      &
    sub_t   => fio_t,   &
    sub_get => fio_get

  use fio_m, only:                                &
    interpolation_t    => fio_interpolation_t,    &
    interpolation_init => fio_interpolation_init, &
    interpolation_eval => fio_interpolation_eval, &
    interpolation_end  => fio_interpolation_end

  use fio_external_m, only:           &
    sub_external_t => fio_external_t

  implicit none

  private
  public ::                        &
    frozen_external_t,             &
    frozen_external_init,          &
    frozen_external_start,         &
    frozen_external_extend,        &
    frozen_external_update,        &
    frozen_external_get,           &
    frozen_external_get_size,      &
    frozen_external_get_energy,    &
    frozen_external_get_potential, &
    frozen_external_copy,          &
    frozen_external_end

  public ::                             &
    frozen_external_interpolation_t,    &
    frozen_external_interpolation_init, &
    frozen_external_interpolation_eval, &
    frozen_external_interpolation_copy, &
    frozen_external_interpolation_end

  interface frozen_external_update
    module procedure frozen_external_update_build
    module procedure frozen_external_update_finalize
  end interface frozen_external_update

contains

  ! ---------------------------------------------------------
  subroutine frozen_external_update_build(this, that)
    type(frozen_external_t), intent(inout) :: this
    type(sub_t),             intent(in)    :: that
    !
    real(kind=wp), dimension(:), pointer :: potential
    type(simulation_t),          pointer :: sim
    type(mesh_t),                pointer :: mesh
    type(sub_external_t),        pointer :: ep
    type(interpolation_t)                :: intrp
    real(kind=wp)                        :: pot
    integer                              :: i, ierr
    !
    nullify(potential, sim, mesh, ep)
    call frozen_external_get_potential(this, potential)
    ASSERT(associated(potential))
    call frozen_external_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    call sub_get(that, ep)
    ASSERT(associated(ep))
    print *, "***: frozen_external_update_build: interpolation_init"
    call interpolation_init(intrp, that, ep)
    do i = 1, frozen_external_get_size(this)
      call interpolation_eval(intrp, mesh%x(i,:), pot)
      potential(i)=potential(i)+pot
    end do
    print *, "***: frozen_external_update_build: interpolation_end"
    call interpolation_end(intrp)
    nullify(potential, sim, mesh, ep)
    return
  end subroutine frozen_external_update_build

  ! ---------------------------------------------------------
  subroutine frozen_external_update_finalize(this)
    type(frozen_external_t), intent(inout) :: this
    !
    call external_potential_update(this)
    return
  end subroutine frozen_external_update_finalize

end module frozen_external_m

!! Local Variables:
!! mode: f90
!! End:
