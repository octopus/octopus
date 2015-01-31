#include "global.h"

module fio_density_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_len, json_object_t, json_array_t,  json_get
  use json_m,      only: json_array_iterator_t, json_init, json_end, json_next
  use kinds_m,     only: wp
  use path_m,      only: path_join

  use fio_simulation_m, only: &
    fio_simulation_t

  use base_density_m, only:                  &
    fio_density_t     => base_density_t,     &
    fio_density_init  => base_density_init,  &
    fio_density_start => base_density_start, &
    fio_density_stop  => base_density_stop,  &
    fio_density_eval  => base_density_eval,  &
    fio_density_get   => base_density_get,   &
    fio_density_copy  => base_density_copy,  &
    fio_density_end   => base_density_end

  use base_density_m, only:                        &
    fio_density_intrpl_t => base_density_intrpl_t

  implicit none

  private
  public ::             &
    fio_density_t,      &
    fio_density_init,   &
    fio_density_start,  &
    fio_density_update, &
    fio_density_stop,   &
    fio_density_eval,   &
    fio_density_get,    &
    fio_density_copy,   &
    fio_density_end

  public ::               &
    fio_density_intrpl_t

contains

  ! ---------------------------------------------------------
  subroutine fio_density_read(this, dir, file, ispin)
    type(fio_density_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dir
    character(len=*),    intent(in)    :: file
    integer,             intent(in)    :: ispin
    !
    real(kind=wp), dimension(:,:), pointer :: dnst
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: np, ierr
    !
    PUSH_SUB(fio_density_read)
    nullify(dnst)
    call fio_density_get(this, dnst)
    ASSERT(associated(dnst))
    call fio_density_get(this, size=np)
    call path_join(dir, file, fpth)
    call io_binary_read(fpth, np, dnst(:,ispin), ierr, offset=0)
    if(ierr/=0)then
      call fio_density_end(this)
      message(1)="Could not read the density file: '"//trim(adjustl(fpth))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(dnst)
    POP_SUB(fio_density_read)
    return
  end subroutine fio_density_read

  ! ---------------------------------------------------------
  subroutine fio_density_update(this)
    type(fio_density_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: isp, nspin, ierr
    !
    PUSH_SUB(fio_density_update)
    nullify(cnfg, list)
    call fio_density_get(this, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(cnfg, "files", list, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_density_get(this, nspin=nspin)
    ASSERT(json_len(list)==nspin)
    isp=0
    call json_init(iter, list)
    do
      call json_next(iter, file, ierr)
      if(ierr/=JSON_OK)exit
      isp=isp+1
      call fio_density_read(this, trim(adjustl(dir)), trim(adjustl(file)), isp)
    end do
    call json_end(iter)
    ASSERT(isp==nspin)
    nullify(cnfg, list)
    POP_SUB(fio_density_update)
    return
  end subroutine fio_density_update

end module fio_density_m

!! Local Variables:
!! mode: f90
!! End:
