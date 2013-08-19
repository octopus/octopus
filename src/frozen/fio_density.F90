#include "global.h"
#include "template.h"

#undef SUBTEMPLATE_NAME
#undef SUBTEMPLATE_TYPE
#undef TEMPLATE_NAME
#undef TEMPLATE_TYPE
#define TEMPLATE_NAME fio
#include "tbase_density.F90"
#undef TEMPLATE_TYPE
#undef TEMPLATE_NAME

module fio_density_m

  use global_m
  use messages_m
  use profiling_m

  use datasets_m,  only: tmpdir
  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_array_t,  json_get
  use json_m,      only: json_array_iterator_t, json_init, json_end, json_next
  use kinds_m,     only: wp

  use fio_simulation_m, only:         &
    simulation_t !=> fio_simulation_t

  use fio_base_density_m, only:                   &
    base_density_start => fio_base_density_start 

  use fio_base_density_m, only:                                          &
    fio_density_t                 => fio_base_density_t,                 &
    fio_density_init              => fio_base_density_init,              &
    fio_density_update            => fio_base_density_update,            &
    fio_density_get               => fio_base_density_get,               &
    fio_density_get_size          => fio_base_density_get_size,          &
    fio_density_get_nspin         => fio_base_density_get_nspin,         &
    fio_density_get_density       => fio_base_density_get_density,       &
    fio_density_get_total_density => fio_base_density_get_total_density, &
    fio_density_copy              => fio_base_density_copy,              &
    fio_density_end               => fio_base_density_end

  use fio_base_density_m, only:                                            &
    fio_density_interpolation_t    => fio_base_density_interpolation_t,    &
    fio_density_interpolation_init => fio_base_density_interpolation_init, &
    fio_density_interpolation_eval => fio_base_density_interpolation_eval, &
    fio_density_interpolation_copy => fio_base_density_interpolation_copy, &
    fio_density_interpolation_end  => fio_base_density_interpolation_end

  implicit none

  private
  public ::                        &
    fio_density_t,                 &
    fio_density_init,              &
    fio_density_start,             &
    fio_density_update,            &
    fio_density_get,               &
    fio_density_get_size,          &
    fio_density_get_nspin,         &
    fio_density_get_density,       &
    fio_density_get_total_density, &
    fio_density_copy,              &
    fio_density_end

  public ::                         &
    fio_density_interpolation_t,    &
    fio_density_interpolation_init, &
    fio_density_interpolation_eval, &
    fio_density_interpolation_copy, &
    fio_density_interpolation_end

contains

  ! ---------------------------------------------------------
  subroutine fio_density_start(this, sim)
    type(fio_density_t), intent(inout) :: this
    type(simulation_t),  intent(in)    :: sim
    !
    real(kind=wp), dimension(:,:), pointer :: density
    type(json_object_t),           pointer :: cnfg
    type(json_array_t),            pointer :: list
    character(len=MAX_PATH_LEN)            :: dir, filename
    type(json_array_iterator_t)            :: iter
    integer                                :: np, isp, ierr
    !
    nullify(density, cnfg, list)
    call fio_density_get(this, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "dir", dir, ierr)
    if(ierr/=JSON_OK)dir='../'//trim(tmpdir)//GS_DIR
    call base_density_start(this, sim)
    call fio_density_get_density(this, density)
    ASSERT(associated(density))
    np=fio_density_get_size(this)
    isp=0
    nullify(list)
    call json_get(cnfg, "files", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
      call json_next(iter, filename, ierr)
      if(ierr/=JSON_OK)exit
      isp=isp+1
      call io_binary_read(trim(dir)//trim(filename), np, density(:,isp), ierr, offset=0)
      if(ierr/=0)then
        call fio_density_end(this)
        message(1)="Could not read the density file: '"//trim(dir)//trim(filename)//"'"
        write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
        call messages_fatal(2)
      end if
    end do
    call json_end(iter)
    ASSERT(isp==fio_density_get_nspin(this))
    nullify(density, cnfg, list)
    return
  end subroutine fio_density_start

end module fio_density_m

#define TEMPLATE_NAME fio
#include "tgeo.F90"
#include "tstates.F90"
#include "tsystem.F90"
#include "tterm.F90"
#include "tpotential.F90"
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
