#include "global.h"

module fio_external_m

  use atom_m
  use global_m
  use messages_m
  use profiling_m
  use io_binary_m
  use json_m
  use kinds_m
  use path_m
  use species_m
  use base_geom_m
  use fio_simulation_m
  use fio_system_m
  use base_potential_m, only:           &
    fio_external_t => base_potential_t
  use base_potential_m

#define TEMPLATE_NAME fio_external
#define INCLUDE_PREFIX
#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_NAME

  implicit none

  private
  public ::         &
    fio_external_t

  public ::               &
    fio_external__load__

  public ::            &
    fio_external_eval

  public ::            &
    fio_external_init, &
    fio_external_get,  &
    fio_external_copy, &
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
    module procedure fio_external_get_system
    module procedure fio_external_get_storage
    module procedure fio_external_get_potential_1d
    module procedure fio_external_get_potential_md
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

    call base_potential__init__(this, sys, config)

    POP_SUB(fio_external_init_potential)
  end subroutine fio_external_init_potential

  ! ---------------------------------------------------------
  subroutine fio_external_init_copy(this, that)
    type(fio_external_t), intent(out) :: this
    type(fio_external_t), intent(in)  :: that

    PUSH_SUB(fio_external_init_copy)

    call base_potential__init__(this, that)

    POP_SUB(fio_external_init_copy)
  end subroutine fio_external_init_copy

  ! ---------------------------------------------------------
  subroutine fio_external__read__(this, dir, file)
    type(fio_external_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    !
    real(kind=wp), dimension(:), pointer :: potn
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: np, ierr
    !
    PUSH_SUB(fio_external__read__)
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
    POP_SUB(fio_external__read__)
    return
  end subroutine fio_external__read__

  ! ---------------------------------------------------------
  subroutine fio_external__load__(this)
    type(fio_external_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ierr

    PUSH_SUB(fio_external__load__)

    nullify(cnfg)
    call fio_external_get(this, cnfg)
    ASSERT(associated(cnfg))
    call json_get(cnfg, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(cnfg, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_external__read__(this, trim(adjustl(dir)), trim(adjustl(file)))
    nullify(cnfg)
    call base_potential__update__(this)

    POP_SUB(fio_external__load__)
  end subroutine fio_external__load__

  ! ---------------------------------------------------------
  subroutine fio_external_get_info(this, size, nspin)
    type(fio_external_t), intent(in)  :: this
    integer,    optional, intent(out) :: size
    integer,    optional, intent(out) :: nspin

    PUSH_SUB(fio_external_get_info)

    call base_potential_get(this, size=size, nspin=nspin)

    POP_SUB(fio_external_get_info)
  end subroutine fio_external_get_info

  ! ---------------------------------------------------------
  subroutine fio_external_get_config(this, that)
    type(fio_external_t), intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(fio_external_get_config)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_config)
  end subroutine fio_external_get_config

  ! ---------------------------------------------------------
  subroutine fio_external_get_simulation(this, that)
    type(fio_external_t),    intent(in) :: this
    type(fio_simulation_t), pointer     :: that

    PUSH_SUB(fio_external_get_simulation)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_simulation)
  end subroutine fio_external_get_simulation

  ! ---------------------------------------------------------
  subroutine fio_external_get_system(this, that)
    type(fio_external_t), intent(in) :: this
    type(fio_system_t),  pointer     :: that

    PUSH_SUB(fio_external_get_system)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_system)
  end subroutine fio_external_get_system

  ! ---------------------------------------------------------
  subroutine fio_external_get_storage(this, that)
    type(fio_external_t), intent(in) :: this
    type(storage_t),     pointer     :: that

    PUSH_SUB(fio_external_get_storage)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_storage)
  end subroutine fio_external_get_storage

  ! ---------------------------------------------------------
  subroutine fio_external_get_potential_1d(this, that)
    type(fio_external_t),         intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that

    PUSH_SUB(fio_external_get_potential_1d)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_potential_1d)
  end subroutine fio_external_get_potential_1d

  ! ---------------------------------------------------------
  subroutine fio_external_get_potential_md(this, that)
    type(fio_external_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(fio_external_get_potential_md)

    call base_potential_get(this, that)

    POP_SUB(fio_external_get_potential_md)
  end subroutine fio_external_get_potential_md

  ! ---------------------------------------------------------
  subroutine fio_external_copy_potential(this, that)
    type(fio_external_t), intent(inout) :: this
    type(fio_external_t), intent(in)    :: that

    PUSH_SUB(fio_external_copy_potential)

    call base_potential__copy__(this, that)

    POP_SUB(fio_external_copy_potential)
  end subroutine fio_external_copy_potential

  ! ---------------------------------------------------------
  subroutine fio_external_end_potential(this)
    type(fio_external_t), intent(inout) :: this

    PUSH_SUB(fio_external_end_potential)

    call base_potential__end__(this)

    POP_SUB(fio_external_end_potential)
  end subroutine fio_external_end_potential

  ! ---------------------------------------------------------
  pure function fio_external_calc(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c

    real(kind=wp) :: v

    real(kind=wp) :: r

    r=sqrt(sum((x-y)**2))
    if(r<r_small) r=r_small
    v=-c/r

  end function fio_external_calc

  ! ---------------------------------------------------------
  subroutine fio_external_classical(this, x, v)
    type(fio_external_t),        intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v

    type(fio_system_t),     pointer :: sys
    type(base_geom_t),      pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(base_geom_iterator_t)      :: iter

    PUSH_SUB(fio_external_classical)

    v=0.0_wp
    nullify(sys, geom)
    call fio_external_get(this, sys)
    ASSERT(associated(sys))
    call fio_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geom_init(iter, geom)
    do
      nullify(atom)
      call base_geom_next(iter, atom)
      if(.not.associated(atom))exit
      v=v+fio_external_calc(x, atom%x, species_zval(atom%species))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call fio_geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+fio_external_calc(x, catom%x, catom%charge)
    !end do
    call base_geom_end(iter)
    !nullify(catom)
    nullify(sys, geom)

    POP_SUB(fio_external_classical)
  end subroutine fio_external_classical

  ! ---------------------------------------------------------
  subroutine fio_external_eval(this, x, v)
    type(fio_external_intrpl_t), intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: v

    type(fio_external_t), pointer :: epot
    integer                       :: ierr

    PUSH_SUB(fio_external_eval)

    nullify(epot)
    call fio_external_intrpl_eval(this, x, v, ierr)
    if(ierr/=FIO_EXTERNAL_INTRPL_OK)then
      call fio_external_get(this, epot)
      ASSERT(associated(epot))
      call fio_external_classical(epot, x, v)
      nullify(epot)
    end if

    POP_SUB(fio_external_intrpl_eval)
  end subroutine fio_external_eval

#define TEMPLATE_NAME fio_external
#define INCLUDE_BODY
#include "intrpl_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_NAME

end module fio_external_m

!! Local Variables:
!! mode: f90
!! End:
