#include "global.h"

#define TEMPLATE_NAME fio
#include "texternal_potential.F90"
#undef TEMPLATE_NAME

module fio_external_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_get
  use kinds_m,     only: wp

  use fio_simulation_m, only:         &
    simulation_t !=> fio_simulation_t

  use fio_external_potential_m, only:                         &
    external_potential_start => fio_external_potential_start

  use fio_external_potential_m, only:                                   &
    fio_external_t             => fio_external_potential_t,             &
    fio_external_init          => fio_external_potential_init,          &
    fio_external_update        => fio_external_potential_update,        &
    fio_external_get           => fio_external_potential_get,           &
    fio_external_get_size      => fio_external_potential_get_size,      &
    fio_external_get_energy    => fio_external_potential_get_energy,    &
    fio_external_get_potential => fio_external_potential_get_potential, &
    fio_external_copy          => fio_external_potential_copy,          &
    fio_external_end           => fio_external_potential_end

  use fio_external_potential_m, only:                                             &
    fio_external_interpolation_t    => fio_external_potential_interpolation_t,    &
    fio_external_interpolation_init => fio_external_potential_interpolation_init, &
    fio_external_interpolation_eval => fio_external_potential_interpolation_eval, &
    fio_external_interpolation_copy => fio_external_potential_interpolation_copy, &
    fio_external_interpolation_end  => fio_external_potential_interpolation_end

  implicit none

  private
  public ::                     &
    fio_external_t,             &
    fio_external_init,          &
    fio_external_start,         &
    fio_external_update,        &
    fio_external_get,           &
    fio_external_get_size,      &
    fio_external_get_energy,    &
    fio_external_get_potential, &
    fio_external_copy,          &
    fio_external_end

  public ::                          &
    fio_external_interpolation_t,    &
    fio_external_interpolation_init, &
    fio_external_interpolation_eval, &
    fio_external_interpolation_copy, &
    fio_external_interpolation_end

contains

  ! ---------------------------------------------------------
  subroutine fio_external_read(this)
    type(fio_external_t), intent(inout) :: this
    !
    real(kind=wp), dimension(:), pointer :: potential
    type(json_object_t),         pointer :: cnfg
    character(len=MAX_PATH_LEN)          :: file
    integer                              :: np, ierr
    !
    nullify(potential, cnfg)
    call fio_external_get(this, cnfg)
    call json_get(cnfg, "dir", file, ierr)
    if(ierr/=JSON_OK)file=STATIC_DIR
    file=trim(adjustl(file))//"v0.obf"
    call fio_external_get_potential(this, potential)
    ASSERT(associated(potential))
    np=fio_external_get_size(this)
    call io_binary_read(file, np, potential, ierr, offset=0)
    if(ierr/=0)then
      call fio_external_end(this)
      message(1)="Could not read the potential file: '"//trim(file)//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(potential)
    return
  end subroutine fio_external_read

  ! ---------------------------------------------------------
  subroutine fio_external_start(this, sim)
    type(fio_external_t), intent(inout) :: this
    type(simulation_t),   intent(in)    :: sim
    !
    real(kind=wp), dimension(:), pointer :: potential
    type(json_object_t), pointer :: cnfg
    logical                      :: read
    integer                      :: ierr
    !
    nullify(cnfg)
    call fio_external_get(this, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "allocate", read, ierr)
    if(ierr/=JSON_OK)read=.false.
    call external_potential_start(this, sim)
    if(read)&
      call fio_external_read(this)
    return
  end subroutine fio_external_start

end module fio_external_m

!! Local Variables:
!! mode: f90
!! End:
