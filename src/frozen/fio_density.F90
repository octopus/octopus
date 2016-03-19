#include "global.h"

module fio_density_oct_m

  use base_density_oct_m
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
  use storage_oct_m

  implicit none

  private

  public ::              &
    fio_density__init__, &
    fio_density__load__

  public ::               &
    fio_density_intrpl_t

  public ::                  &
    fio_density_intrpl_init, &
    fio_density_intrpl_eval, &
    fio_density_intrpl_get,  &
    fio_density_intrpl_copy, &
    fio_density_intrpl_end

  type :: fio_density_intrpl_t
    private
    type(base_density_t), pointer :: self =>null()
    type(intrpl_t)                :: intrp
  end type fio_density_intrpl_t

  interface fio_density_intrpl_eval
    module procedure fio_density_intrpl_eval_1d
    module procedure fio_density_intrpl_eval_md
  end interface fio_density_intrpl_eval

  interface fio_density_intrpl_get
    module procedure fio_density_intrpl_get_info
    module procedure fio_density_intrpl_get_type
  end interface fio_density_intrpl_get

contains

  ! ---------------------------------------------------------
  subroutine fio_density__init__(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    real(kind=wp)                :: chrg
    integer                      :: ispn, nspin, ierr
    logical                      :: ipol

    PUSH_SUB(fio_density__init__)

    nullify(cnfg, list)
    call base_density_get(this, cnfg)
    ASSERT(associated(cnfg))
    call base_density_get(this, nspin=nspin)
    call json_get(cnfg, "invertpolarization", ipol, ierr)
    if(ierr/=JSON_OK) ipol = .false.
    if(ipol.and.(nspin>1))then
      call json_get(cnfg, "charge", list, ierr)
      ASSERT(ierr==JSON_OK)
      do ispn = 1, nspin
        call json_get(list, ispn, chrg, ierr)
        ASSERT(ierr==JSON_OK)
        call base_density_set(this, chrg, nspin-ispn+1)
      end do
      nullify(list)
    end if
    nullify(cnfg)

    POP_SUB(fio_density__init__)
  end subroutine fio_density__init__

  ! ---------------------------------------------------------
  subroutine fio_density__read__(this, dir, file, ispin)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    integer,              intent(in)    :: ispin

    real(kind=wp), dimension(:,:), pointer :: dnst
    type(simulation_t),            pointer :: sim
    type(mesh_t),                  pointer :: mesh
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: ierr

    PUSH_SUB(fio_density__read__)

    nullify(dnst, sim, mesh)
    call path_join(dir, file, fpth)
    call base_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh, fine=.true.)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dio_function_input(fpth, mesh, dnst(:,ispin), ierr)
    nullify(dnst, mesh)
    if(ierr/=0)then
      call base_density_end(this)
      message(1) = "Could not read the input file: '"//trim(adjustl(file))//".obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

    POP_SUB(fio_density__read__)
  end subroutine fio_density__read__

  ! ---------------------------------------------------------
  subroutine fio_density__load__(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ispn, nspin, ierr
    logical                      :: load, ipol

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
      call json_get(cnfg, "invertpolarization", ipol, ierr)
      if(ierr/=JSON_OK) ipol = .false.
      do ispn = 1, nspin
        call json_get(list, ispn, file, ierr)
        ASSERT(ierr==JSON_OK)
        if(ipol)then
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), nspin-ispn+1)
        else
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), ispn)
        end if
      end do
      nullify(cnfg, list)
      call base_density__update__(this)
    end if

    POP_SUB(fio_density__load__)
  end subroutine fio_density__load__

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_init(this, that, type)
    type(fio_density_intrpl_t),   intent(out) :: this
    type(base_density_t), target, intent(in)  :: that
    integer,            optional, intent(in)  :: type

    type(storage_t), pointer :: data

    PUSH_SUB(fio_density_intrpl_init)

    nullify(data)
    this%self => that
    call base_density_get(that, data)
    ASSERT(associated(data))
    call intrpl_init(this%intrp, data, type=type)
    nullify(data)

    POP_SUB(fio_density_intrpl_init)
  end subroutine fio_density_intrpl_init

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_eval_1d(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_density_intrpl_eval_1d)

    ierr = INTRPL_NI
    if(associated(this%self)) call intrpl_eval(this%intrp, x, val, ierr)
    if(ierr/=INTRPL_OK) val = 0.0_wp

    POP_SUB(fio_density_intrpl_eval_1d)
  end subroutine fio_density_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_eval_md(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val

    integer :: ierr

    PUSH_SUB(fio_density_intrpl_eval_md)

    ierr = INTRPL_NI
    if(associated(this%self)) call intrpl_eval(this%intrp, x, val, ierr)
    if(ierr/=INTRPL_OK) val = 0.0_wp

    POP_SUB(fio_density_intrpl_eval_md)
  end subroutine fio_density_intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_get_info(this, type, dim, default)
    type(fio_density_intrpl_t), intent(in)  :: this
    integer,          optional, intent(out) :: type
    integer,          optional, intent(out) :: dim
    real(kind=wp),    optional, intent(out) :: default

    PUSH_SUB(fio_density_intrpl_get_info)

    call intrpl_get(this%intrp, type, dim, default)

    POP_SUB(fio_density_intrpl_get_info)
  end subroutine fio_density_intrpl_get_info

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_get_type(this, that)
    type(fio_density_intrpl_t), intent(in) :: this
    type(base_density_t),      pointer     :: that

    PUSH_SUB(fio_density_intrpl_get_type)

    nullify(that)
    if(associated(this%self)) that => this%self

    POP_SUB(fio_density_intrpl_get_type)
  end subroutine fio_density_intrpl_get_type

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_copy(this, that)
    type(fio_density_intrpl_t), intent(out) :: this
    type(fio_density_intrpl_t), intent(in)  :: that

    PUSH_SUB(fio_density_intrpl_copy)

    call fio_density_intrpl_end(this)
    if(associated(that%self))then
      this%self => that%self
      call intrpl_copy(this%intrp, that%intrp)
    end if

    POP_SUB(fio_density_intrpl_copy)
  end subroutine fio_density_intrpl_copy

  ! ---------------------------------------------------------
  subroutine fio_density_intrpl_end(this)
    type(fio_density_intrpl_t), intent(inout) :: this

    PUSH_SUB(fio_density_intrpl_end)

    if(associated(this%self)) call intrpl_end(this%intrp)
    nullify(this%self)

    POP_SUB(fio_density_intrpl_end)
  end subroutine fio_density_intrpl_end

end module fio_density_oct_m

!! Local Variables:
!! mode: f90
!! End:
