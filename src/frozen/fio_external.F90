#include "global.h"

module fio_external_oct_m

  use atom_oct_m
  use base_geometry_oct_m
  use base_potential_oct_m
  use base_system_oct_m
  use global_oct_m
  use intrpl_oct_m
  use io_function_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use path_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use species_oct_m
  use storage_oct_m

  implicit none

  private

  public ::               &
    fio_external__load__

  public ::                &
    fio_external_intrpl_t

  public ::                   &
    fio_external_intrpl_init, &
    fio_external_intrpl_eval, &
    fio_external_intrpl_get,  &
    fio_external_intrpl_copy, &
    fio_external_intrpl_end

  type :: fio_external_intrpl_t
    private
    type(base_potential_t), pointer :: self =>null()
    type(intrpl_t)                  :: intrp
  end type fio_external_intrpl_t

  interface fio_external_intrpl_get
    module procedure fio_external_intrpl_get_info
    module procedure fio_external_intrpl_get_type
  end interface fio_external_intrpl_get

contains

  ! ---------------------------------------------------------
  subroutine fio_external__read__(this, dir, file)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: dir
    character(len=*),       intent(in)    :: file

    real(kind=wp), dimension(:), pointer :: potn
    type(simulation_t),          pointer :: sim
    type(mesh_t),                pointer :: mesh
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: ierr

    PUSH_SUB(fio_external__read__)

    nullify(potn, sim, mesh)
    call path_join(dir, file, fpth)
    call base_potential_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_potential_get(this, potn)
    ASSERT(associated(potn))
    call dio_function_input(fpth, mesh, potn, ierr)
    nullify(potn, mesh)
    if(ierr/=0)then
      call base_potential_end(this)
      message(1) = "Could not read the input file: '"//trim(adjustl(file))//".obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

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
  subroutine fio_external_intrpl_init(this, that, type)
    type(fio_external_intrpl_t),    intent(out) :: this
    type(base_potential_t), target, intent(in)  :: that
    integer,              optional, intent(in)  :: type

    type(storage_t), pointer :: data

    PUSH_SUB(fio_external_intrpl_init)

    nullify(data)
    this%self => that
    call base_potential_get(that, data)
    ASSERT(associated(data))
    call intrpl_init(this%intrp, data, type=type)
    nullify(data)

    POP_SUB(fio_external_intrpl_init)
  end subroutine fio_external_intrpl_init

  ! ---------------------------------------------------------
  pure function fio_external_calc(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c

    real(kind=wp) :: v

    real(kind=wp) :: r
    integer       :: n

    n = min(size(x), size(y))
    r = sqrt(sum((x(1:n)-y(1:n))**2))
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
  subroutine fio_external_intrpl_eval(this, x, val)
    type(fio_external_intrpl_t), intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_external_intrpl_eval)

    ierr = INTRPL_NI
    if(associated(this%self)) call intrpl_eval(this%intrp, x, val, ierr)
    if(ierr/=INTRPL_OK) call fio_external_classical(this%self, x, val)

    POP_SUB(fio_external_intrpl_eval)
  end subroutine fio_external_intrpl_eval

  ! ---------------------------------------------------------
  subroutine fio_external_intrpl_get_info(this, type, dim, default)
    type(fio_external_intrpl_t), intent(in)  :: this
    integer,           optional, intent(out) :: type
    integer,           optional, intent(out) :: dim
    real(kind=wp),     optional, intent(out) :: default

    PUSH_SUB(fio_external_intrpl_get_info)

    call intrpl_get(this%intrp, type, dim, default)

    POP_SUB(fio_external_intrpl_get_info)
  end subroutine fio_external_intrpl_get_info

  ! ---------------------------------------------------------
  subroutine fio_external_intrpl_get_type(this, that)
    type(fio_external_intrpl_t), intent(in) :: this
    type(base_potential_t),     pointer     :: that

    PUSH_SUB(fio_external_intrpl_get_type)

    nullify(that)
    if(associated(this%self)) that => this%self

    POP_SUB(fio_external_intrpl_get_type)
  end subroutine fio_external_intrpl_get_type

  ! ---------------------------------------------------------
  subroutine fio_external_intrpl_copy(this, that)
    type(fio_external_intrpl_t), intent(out) :: this
    type(fio_external_intrpl_t), intent(in)  :: that

    PUSH_SUB(fio_external_intrpl_copy)

    call fio_external_intrpl_end(this)
    if(associated(that%self))then
      this%self => that%self
      call intrpl_copy(this%intrp, that%intrp)
    end if

    POP_SUB(fio_external_intrpl_copy)
  end subroutine fio_external_intrpl_copy

  ! ---------------------------------------------------------
  subroutine fio_external_intrpl_end(this)
    type(fio_external_intrpl_t), intent(inout) :: this

    PUSH_SUB(fio_external_intrpl_end)

    if(associated(this%self)) call intrpl_end(this%intrp)
    nullify(this%self)

    POP_SUB(fio_external_intrpl_end)
  end subroutine fio_external_intrpl_end

end module fio_external_oct_m

!! Local Variables:
!! mode: f90
!! End:
