#include "global.h"

module fio_external_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_get
  use kinds_m,     only: wp
  use path_m,      only: path_join

  use bepot_m, only:                                &
    fio_external_t          => bepot_t,             &
    fio_external_init       => bepot_init,          &
    fio_external_start      => bepot_start,         &
    fio_external_stop       => bepot_stop,          &
    fio_external_get        => bepot_get,           &
    fio_external_get_size   => bepot_get_size,      &
    fio_external_get_energy => bepot_get_energy,    &
    fio_external_copy       => bepot_copy,          &
    fio_external_end        => bepot_end

  use bepot_m, only:                         &
    fio_external_intrpl_t => bepot_intrpl_t

  implicit none

  private
  public ::                     &
    fio_external_t,             &
    fio_external_init,          &
    fio_external_start,         &
    fio_external_update,        &
    fio_external_stop,          &
    fio_external_get,           &
    fio_external_get_size,      &
    fio_external_get_energy,    &
    fio_external_copy,          &
    fio_external_end

  public ::                       &
    fio_external_intrpl_t

contains

  ! ---------------------------------------------------------
  subroutine fio_external_read(this, dir, file)
    type(fio_external_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    !
    real(kind=wp), dimension(:), pointer :: potential
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: np, ierr
    !
    PUSH_SUB(fio_external_read)
    nullify(potential)
    call fio_external_get(this, potential)
    ASSERT(associated(potential))
    np=fio_external_get_size(this)
    call path_join(dir, file, fpth)
    call io_binary_read(fpth, np, potential, ierr, offset=0)
    if(ierr/=0)then
      call fio_external_end(this)
      message(1)="Could not read the potential file: '"//trim(adjustl(fpth))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(potential)
    POP_SUB(fio_external_read)
    return
  end subroutine fio_external_read

  ! ---------------------------------------------------------
  subroutine fio_external_update(this)
    type(fio_external_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
    logical                      :: read
    integer                      :: ierr
    !
    PUSH_SUB(fio_external_update)
    nullify(cnfg)
    call fio_external_get(this, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(cnfg, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_external_read(this, trim(adjustl(dir)), trim(adjustl(file)))
    nullify(cnfg)
    POP_SUB(fio_external_update)
    return
  end subroutine fio_external_update

end module fio_external_m

!! Local Variables:
!! mode: f90
!! End:
