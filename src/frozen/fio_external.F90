#include "global.h"

module fio_external_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_get
  use kinds_m,     only: wp
  use path_m,      only: path_join

  use fio_system_m, only: &
    fio_system_t

  use base_external_m, only: &
    base_external__init__,   &
    base_external__update__, &
    base_external__copy__,   &
    base_external__end__

  use base_external_m, only:                      &
    fio_external_start => base_external__start__, &
    fio_external_stop  => base_external__stop__

  use base_external_m, only: &
    base_external_get

  use base_external_m, only:           &
    fio_external_t => base_external_t

#define TEMPLATE_NAME fio_external
#define INCLUDE_PREFIX
#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_NAME

  implicit none

  private
  public ::              &
    fio_external_t,      &
    fio_external_init,   &
    fio_external_start,  &
    fio_external_update, &
    fio_external_stop,   &
    fio_external_get,    &
    fio_external_copy,   &
    fio_external_end

#define TEMPLATE_NAME fio_external
#define INCLUDE_HEADER
#include "intrpl_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_NAME

  interface fio_external_init
    module procedure fio_external_init_potential
    module procedure fio_external_init_copy
  end interface fio_external_init

  interface fio_external_get
    module procedure fio_external_get_info
    module procedure fio_external_get_config
    module procedure fio_external_get_simulation
    module procedure fio_external_get_potential_1d
    module procedure fio_external_get_potential_2d
  end interface fio_external_get

  interface fio_external_copy
    module procedure fio_external_copy_potential
  end interface fio_external_copy

  interface fio_external_end
    module procedure fio_external_end_potential
  end interface fio_external_end

contains

  ! ---------------------------------------------------------
  subroutine fio_external_init_potential(this, sys, config)
    type(fio_external_t), intent(out) :: this
    type(fio_system_t),   intent(out) :: sys
    type(json_object_t),  intent(in)  :: config

    PUSH_SUB(fio_external_init_potential)

    call base_external__init__(this, sys, config)

    POP_SUB(fio_external_init_potential)
  end subroutine fio_external_init_potential

  ! ---------------------------------------------------------
  subroutine fio_external_init_copy(this, that)
    type(fio_external_t), intent(out) :: this
    type(fio_external_t), intent(in)  :: that

    PUSH_SUB(fio_external_init_copy)

    call base_external__init__(this, that)

    POP_SUB(fio_external_init_copy)
  end subroutine fio_external_init_copy

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
    call fio_external_get(this, size=np)
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
    call base_external__update__(this)
    POP_SUB(fio_external_update)
    return
  end subroutine fio_external_update

  ! ---------------------------------------------------------
  subroutine fio_external_get_info(this, size, nspin)
    type(fio_external_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin

    PUSH_SUB(fio_external_get_info)

    call base_external_get(this, size=size, nspin=nspin)

    POP_SUB(fio_external_get_info)
  end subroutine fio_external_get_info

  ! ---------------------------------------------------------
  subroutine fio_external_get_config(this, that)
    type(fio_external_t), target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(fio_external_get_config)

    call base_external_get(this, that)

    POP_SUB(fio_external_get_config)
  end subroutine fio_external_get_config

  ! ---------------------------------------------------------
  subroutine fio_external_get_simulation(this, that)
    type(fio_external_t), target, intent(in) :: this
    type(simulation_t),  pointer             :: that

    PUSH_SUB(fio_external_get_simulation)

    call base_external_get(this, that)

    POP_SUB(fio_external_get_simulation)
  end subroutine fio_external_get_simulation

  ! ---------------------------------------------------------
  subroutine fio_external_get_potential_1d(this, that, total)
    type(fio_external_t),                   intent(in) :: this
    real(kind=wp),           dimension(:), pointer     :: that
    real(kind=wp), optional, dimension(:), pointer     :: total

    PUSH_SUB(fio_external_get_potential_1d)

    call fio_external_get(this, that, total)

    POP_SUB(fio_external_get_potential_1d)
  end subroutine fio_external_get_potential_1d

  ! ---------------------------------------------------------
  subroutine fio_external_get_potential_2d(this, that)
    type(fio_external_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(external_get_fio_potential_2d)

    call base_external_get(this, that)

    POP_SUB(fio_external_get_potential_2d)
  end subroutine fio_external_get_potential_2d

  ! ---------------------------------------------------------
  subroutine fio_external_copy_potential(this, that)
    type(fio_external_t), intent(inout) :: this
    type(fio_external_t), intent(in)    :: that

    PUSH_SUB(fio_external_copy_potential)

    call base_external__copy__(this, that)

    POP_SUB(fio_external_copy_potential)
  end subroutine fio_external_copy_potential

  ! ---------------------------------------------------------
  subroutine fio_external_end_potential(this)
    type(fio_external_t), intent(inout) :: this

    PUSH_SUB(fio_external_end_potential)

    call base_external__end__(this)

    POP_SUB(fio_external_end_potential)
  end subroutine fio_external_end_potential

#define TEMPLATE_NAME fio_external
#define INCLUDE_BODY
#include "intrpl_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_NAME

end module fio_external_m

!! Local Variables:
!! mode: f90
!! End:
