#include "global.h"

module fio_external_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_get
  use kinds_m,     only: wp
  use path_m,      only: path_join

  use base_external_m, only:                   &
    fio_external_t     => base_external_t,     &
    fio_external_init  => base_external_init,  &
    fio_external_start => base_external_start, &
    fio_external_stop  => base_external_stop,  &
    fio_external_eval  => base_external_eval,  &
    fio_external_get   => base_external_get,   &
    fio_external_copy  => base_external_copy,  &
    fio_external_end   => base_external_end

  use base_external_m, only:                         &
    fio_external_intrpl_t => base_external_intrpl_t

  implicit none

  private
  public ::              &
    fio_external_t,      &
    fio_external_init,   &
    fio_external_start,  &
    fio_external_update, &
    fio_external_stop,   &
    fio_external_eval,   &
    fio_external_get,    &
    fio_external_copy,   &
    fio_external_end

  public ::                &
    fio_external_intrpl_t

contains

  ! ---------------------------------------------------------
  subroutine fio_external_read(this, dir, file)
    type(fio_external_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    !
    real(kind=wp), dimension(:), pointer :: potn
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: np, ierr
    !
    PUSH_SUB(fio_external_read)
    nullify(potn)
    call fio_external_get(this, potn)
    ASSERT(associated(potn))
!    call fio_external_get(this, size=np)
    call path_join(dir, file, fpth)
    call io_binary_read(fpth, np, potn, ierr, offset=0)
    if(ierr/=0)then
      call fio_external_end(this)
      message(1)="Could not read the potential file: '"//trim(adjustl(fpth))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(potn)
    POP_SUB(fio_external_read)
    return
  end subroutine fio_external_read

  ! ---------------------------------------------------------
  subroutine fio_external_update(this)
    type(fio_external_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
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
