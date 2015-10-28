#include "global.h"

#define EXTERNAL_LEVEL base
#define TEMPLATE_LEVEL fio
#define EXTERNAL_NAME potential
#define TEMPLATE_NAME external

module fio_external_m

  use atom_m
  use base_geometry_m
  use base_potential_m
  use base_system_m
  use global_m
  use io_binary_m
  use json_m
  use kinds_m
  use messages_m
  use path_m
  use profiling_m
  use species_m

#define INCLUDE_PREFIX
#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX

  implicit none

  private

  public ::               &
    fio_external__load__

  public ::            &
    fio_external_init, &
    fio_external_eval, &
    fio_external_get,  &
    fio_external_copy, &
    fio_external_end

#define INCLUDE_HEADER
#include "intrpl_inc.F90"
#undef INCLUDE_HEADER

  interface fio_external_init
    module procedure fio_external_init_type
    module procedure fio_external_init_copy
  end interface fio_external_init

  interface fio_external_copy
    module procedure fio_external_copy_type
  end interface fio_external_copy

  interface fio_external_end
    module procedure fio_external_end_type
  end interface fio_external_end

contains

  ! ---------------------------------------------------------
  subroutine fio_external_init_type(this, sys, config)
    type(base_potential_t), intent(out) :: this
    type(base_system_t),    intent(out) :: sys
    type(json_object_t),    intent(in)  :: config

    PUSH_SUB(fio_external_init_type)

    call base_potential__init__(this, sys, config)

    POP_SUB(fio_external_init_type)
  end subroutine fio_external_init_type

  ! ---------------------------------------------------------
  subroutine fio_external_init_copy(this, that)
    type(base_potential_t), intent(out) :: this
    type(base_potential_t), intent(in)  :: that

    PUSH_SUB(fio_external_init_copy)

    call base_potential__init__(this, that)

    POP_SUB(fio_external_init_copy)
  end subroutine fio_external_init_copy

  ! ---------------------------------------------------------
  subroutine fio_external__read__(this, dir, file)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: dir
    character(len=*),       intent(in)    :: file

    real(kind=wp), dimension(:), pointer :: potn
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: np, ierr

    PUSH_SUB(fio_external__read__)

    nullify(potn)
    call base_potential_get(this, potn)
    ASSERT(associated(potn))
    call base_potential_get(this, size=np)
    call path_join(dir, file, fpth)
    call io_binary_read(fpth, np, potn, ierr, offset=0)
    if(ierr/=0)then
      call fio_external_end(this)
      message(1) = "Could not read the potential file: '"//trim(adjustl(fpth))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    nullify(potn)

    POP_SUB(fio_external__read__)
  end subroutine fio_external__read__

  ! ---------------------------------------------------------
  subroutine fio_external__load__(this)
    type(base_potential_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ierr
    logical                      :: load

    PUSH_SUB(fio_external__load__)

    nullify(cnfg)
    call base_potential_get(this, use=load)
    if(load)then
      call base_potential_get(this, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "dir", dir, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "file", file, ierr)
      ASSERT(ierr==JSON_OK)
      call fio_external__read__(this, trim(adjustl(dir)), trim(adjustl(file)))
      nullify(cnfg)
      call base_potential__update__(this)
    end if

    POP_SUB(fio_external__load__)
  end subroutine fio_external__load__

  ! ---------------------------------------------------------
  subroutine fio_external_copy_type(this, that)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that

    PUSH_SUB(fio_external_copy_type)

    call base_potential__copy__(this, that)

    POP_SUB(fio_external_copy_type)
  end subroutine fio_external_copy_type

  ! ---------------------------------------------------------
  subroutine fio_external_end_type(this)
    type(base_potential_t), intent(inout) :: this

    PUSH_SUB(fio_external_end_type)

    call base_potential__end__(this)

    POP_SUB(fio_external_end_type)
  end subroutine fio_external_end_type

  ! ---------------------------------------------------------
  pure function fio_external_calc(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c

    real(kind=wp) :: v

    real(kind=wp) :: r

    r = sqrt(sum((x-y)**2))
    if(r<r_small) r = r_small
    v = -c / r

  end function fio_external_calc

  ! ---------------------------------------------------------
  subroutine fio_external_classical(this, x, v)
    type(base_potential_t),      intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v

    type(base_system_t),    pointer :: sys
    type(base_geometry_t),  pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(base_geometry_iterator_t)  :: iter

    PUSH_SUB(fio_external_classical)

    v = 0.0_wp
    nullify(sys, geom)
    call base_potential_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geometry_init(iter, geom)
    do
      nullify(atom)
      call base_geometry_next(iter, atom)
      if(.not.associated(atom))exit
      ASSERT(associated(atom%species))
      v = v + fio_external_calc(x, atom%x, species_zval(atom%species))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call fio_geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+fio_external_calc(x, catom%x, catom%charge)
    !end do
    call base_geometry_end(iter)
    !nullify(catom)
    nullify(sys, geom)

    POP_SUB(fio_external_classical)
  end subroutine fio_external_classical

  ! ---------------------------------------------------------
  subroutine fio_external_eval(this, x, v)
    type(fio_external_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v


    integer :: ierr

    PUSH_SUB(fio_external_eval)

    call fio_external_intrpl__eval__(this, x, v, ierr)
    if(ierr/=FIO_EXTERNAL_INTRPL_OK) call fio_external_classical(this%self, x, v)

    POP_SUB(fio_external_eval)
  end subroutine fio_external_eval

#define INCLUDE_BODY
#include "intrpl_inc.F90"
#undef INCLUDE_BODY

end module fio_external_m

#undef EXTERNAL_LEVEL
#undef TEMPLATE_LEVEL
#undef EXTERNAL_NAME
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
