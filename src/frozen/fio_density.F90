#include "global.h"

#define EXTERNAL_LEVEL base
#define TEMPLATE_LEVEL fio
#define TEMPLATE_NAME density

module fio_density_m

  use base_density_m
  use global_m
  use io_binary_m
  use json_m
  use kinds_m
  use messages_m
  use path_m
  use profiling_m

#define INCLUDE_PREFIX
#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX

  implicit none

  private

  public ::              &
    fio_density__load__

  public ::           &
    fio_density_init, &
    fio_density_eval, &
    fio_density_get,  &
    fio_density_copy, &
    fio_density_end

#define INCLUDE_HEADER
#include "intrpl_inc.F90"
#undef INCLUDE_HEADER

  interface fio_density_init
    module procedure fio_density_init_type
    module procedure fio_density_init_copy
  end interface fio_density_init

  interface fio_density_eval
    module procedure fio_density_eval_1d
    module procedure fio_density_eval_md
  end interface fio_density_eval

  interface fio_density_copy
    module procedure fio_density_copy_type
  end interface fio_density_copy

  interface fio_density_end
    module procedure fio_density_end_type
  end interface fio_density_end

contains

  ! ---------------------------------------------------------
  subroutine fio_density_init_type(this, config)
    type(base_density_t), intent(out) :: this
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(fio_density_init_type)

    call base_density__init__(this, config)

    POP_SUB(fio_density_init_type)
  end subroutine fio_density_init_type

  ! ---------------------------------------------------------
  subroutine fio_density_init_copy(this, that)
    type(base_density_t), intent(out) :: this
    type(base_density_t), intent(in)  :: that

    PUSH_SUB(fio_density_init_copy)

    call base_density__init__(this, that)

    POP_SUB(fio_density_init_copy)
  end subroutine fio_density_init_copy

  ! ---------------------------------------------------------
  subroutine fio_density__read__(this, dir, file, ispin)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    integer,              intent(in)    :: ispin

    real(kind=wp), dimension(:,:), pointer :: dnst
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: np, ierr

    PUSH_SUB(fio_density__read__)

    nullify(dnst)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call base_density_get(this, size=np)
    call path_join(dir, file, fpth)
    call io_binary_read(fpth, np, dnst(:,ispin), ierr, offset=0)
    if(ierr/=0)then
      call fio_density_end(this)
      message(1) = "Could not read the density file: '"//trim(adjustl(fpth))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(dnst)

    POP_SUB(fio_density__read__)
  end subroutine fio_density__read__

  ! ---------------------------------------------------------
  subroutine fio_density__load__(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: isp, nspin, ierr
    logical                      :: load

    PUSH_SUB(fio_density__load__)

    nullify(cnfg, list)
    call base_density_get(this, use=load)
    if(load)then
      call base_density_get(this, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "dir", dir, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "files", list, ierr)
      ASSERT(ierr==JSON_OK)
      call base_density_get(this, nspin=nspin)
      ASSERT(json_len(list)==nspin)
      isp = 0
      call json_init(iter, list)
      do
        call json_next(iter, file, ierr)
        if(ierr/=JSON_OK)exit
        isp = isp + 1
        call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), isp)
      end do
      call json_end(iter)
      ASSERT(isp==nspin)
      nullify(cnfg, list)
      call base_density__update__(this)
    end if

    POP_SUB(fio_density__load__)
  end subroutine fio_density__load__

  ! ---------------------------------------------------------
  subroutine fio_density_copy_type(this, that)
    type(base_density_t), intent(inout) :: this
    type(base_density_t), intent(in)    :: that

    PUSH_SUB(fio_density_copy_type)

    call base_density__copy__(this, that)

    POP_SUB(fio_density_copy_type)
  end subroutine fio_density_copy_type

  ! ---------------------------------------------------------
  subroutine fio_density_end_type(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(fio_density_end_type)

    call base_density__end__(this)

    POP_SUB(fio_density_end_type)
  end subroutine fio_density_end_type

  ! ---------------------------------------------------------
  subroutine fio_density_eval_1d(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_density_eval_1d)

    call fio_density_intrpl__eval__(this, x, val, ierr)
    if(ierr/=FIO_DENSITY_INTRPL_OK) val = 0.0_wp

    POP_SUB(fio_density_eval_1d)
  end subroutine fio_density_eval_1d

  ! ---------------------------------------------------------
  subroutine fio_density_eval_md(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_density_eval_md)

    call fio_density_intrpl__eval__(this, x, val, ierr)
    if(ierr/=FIO_DENSITY_INTRPL_OK) val = 0.0_wp

    POP_SUB(fio_density_eval_md)
  end subroutine fio_density_eval_md

#define INCLUDE_BODY
#include "intrpl_inc.F90"
#undef INCLUDE_BODY

end module fio_density_m

#undef EXTERNAL_LEVEL
#undef TEMPLATE_LEVEL
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
