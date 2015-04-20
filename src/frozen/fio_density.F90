#include "global.h"

module fio_density_m

  use global_m
  use messages_m
  use profiling_m

  use io_binary_m, only: io_binary_read
  use json_m,      only: JSON_OK, json_object_t, json_array_t, json_array_iterator_t
  use json_m,      only: json_len, json_init, json_next, json_get, json_end
  use kinds_m,     only: wp
  use path_m,      only: path_join

  use base_density_m, only:          &
    fio_density_t => base_density_t

  use base_density_m, only: &
    base_density__init__,   &
    base_density__update__, &
    base_density__copy__,   &
    base_density__end__

  use base_density_m, only: &
    base_density_get

#define TEMPLATE_NAME fio_density
#define INCLUDE_PREFIX
#include "intrpl_inc.F90"
#undef INCLUDE_PREFIX
#undef TEMPLATE_NAME

  implicit none

  private
  public ::        &
    fio_density_t

  public ::              &
    fio_density__load__

  public ::           &
    fio_density_eval

  public ::           &
    fio_density_init, &
    fio_density_get,  &
    fio_density_copy, &
    fio_density_end

#define TEMPLATE_NAME fio_density
#define INCLUDE_HEADER
#include "intrpl_inc.F90"
#undef INCLUDE_HEADER
#undef TEMPLATE_NAME

  interface fio_density_init
    module procedure fio_density_init_density
    module procedure fio_density_init_copy
  end interface fio_density_init

  interface fio_density_eval
    module procedure fio_density_eval_1d
    module procedure fio_density_eval_md
  end interface fio_density_eval

  interface fio_density_get
    module procedure fio_density_get_info
    module procedure fio_density_get_config
    module procedure fio_density_get_simulation
    module procedure fio_density_get_density_1d
    module procedure fio_density_get_density_2d
  end interface fio_density_get

  interface fio_density_copy
    module procedure fio_density_copy_density
  end interface fio_density_copy

  interface fio_density_end
    module procedure fio_density_end_density
  end interface fio_density_end

contains

  ! ---------------------------------------------------------
  subroutine fio_density_init_density(this, config)
    type(fio_density_t), intent(out) :: this
    type(json_object_t), intent(in)  :: config

    PUSH_SUB(fio_density_init_density)

    call base_density__init__(this, config)

    POP_SUB(fio_density_init_density)
  end subroutine fio_density_init_density

  ! ---------------------------------------------------------
  subroutine fio_density_init_copy(this, that)
    type(fio_density_t), intent(out) :: this
    type(fio_density_t), intent(in)  :: that

    PUSH_SUB(fio_density_init_copy)

    call base_density__init__(this, that)

    POP_SUB(fio_density_init_copy)
  end subroutine fio_density_init_copy

  ! ---------------------------------------------------------
  subroutine fio_density__read__(this, dir, file, ispin)
    type(fio_density_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dir
    character(len=*),    intent(in)    :: file
    integer,             intent(in)    :: ispin

    real(kind=wp), dimension(:,:), pointer :: dnst
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: np, ierr

    PUSH_SUB(fio_density__read__)

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

    POP_SUB(fio_density__read__)
  end subroutine fio_density__read__

  ! ---------------------------------------------------------
  subroutine fio_density__load__(this)
    type(fio_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: isp, nspin, ierr

    PUSH_SUB(fio_density__load__)

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
      call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), isp)
    end do
    call json_end(iter)
    ASSERT(isp==nspin)
    nullify(cnfg, list)
    call base_density__update__(this)

    POP_SUB(fio_density__load__)
  end subroutine fio_density__load__

  ! ---------------------------------------------------------
  subroutine fio_density_get_info(this, size, nspin, fine)
    type(fio_density_t), intent(in)  :: this
    integer,   optional, intent(out) :: size
    integer,   optional, intent(out) :: nspin
    logical,   optional, intent(out) :: fine

    PUSH_SUB(fio_density_get_info)

    call base_density_get(this, size=size, nspin=nspin, fine=fine)

    POP_SUB(fio_density_get_info)
  end subroutine fio_density_get_info

  ! ---------------------------------------------------------
  subroutine fio_density_get_config(this, that)
    type(fio_density_t),  intent(in) :: this
    type(json_object_t), pointer     :: that

    PUSH_SUB(fio_density_get_config)

    call base_density_get(this, that)

    POP_SUB(fio_density_get_config)
  end subroutine fio_density_get_config

  ! ---------------------------------------------------------
  subroutine fio_density_get_simulation(this, that)
    type(fio_density_t),  intent(in) :: this
    type(simulation_t),  pointer     :: that

    PUSH_SUB(fio_density_get_simulation)

    call base_density_get(this, that)

    POP_SUB(fio_density_get_simulation)
  end subroutine fio_density_get_simulation

  ! ---------------------------------------------------------
  subroutine fio_density_get_density_1d(this, that, total)
    type(fio_density_t),          intent(in) :: this
    real(kind=wp), dimension(:), pointer     :: that
    logical,            optional, intent(in) :: total

    PUSH_SUB(fio_density_get_density_1d)

    call base_density_get(this, that, total)

    POP_SUB(fio_density_get_density_1d)
  end subroutine fio_density_get_density_1d

  ! ---------------------------------------------------------
  subroutine fio_density_get_density_2d(this, that)
    type(fio_density_t),            intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that

    PUSH_SUB(density_get_fio_density_2d)

    call base_density_get(this, that)

    POP_SUB(fio_density_get_density_2d)
  end subroutine fio_density_get_density_2d

  ! ---------------------------------------------------------
  subroutine fio_density_copy_density(this, that)
    type(fio_density_t), intent(inout) :: this
    type(fio_density_t), intent(in)    :: that

    PUSH_SUB(fio_density_copy_density)

    call base_density__copy__(this, that)

    POP_SUB(fio_density_copy_density)
  end subroutine fio_density_copy_density

  ! ---------------------------------------------------------
  subroutine fio_density_end_density(this)
    type(fio_density_t), intent(inout) :: this

    PUSH_SUB(fio_density_end_density)

    call base_density__end__(this)

    POP_SUB(fio_density_end_density)
  end subroutine fio_density_end_density

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin_1_n(this, that)
    real(kind=wp),               intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(that))
    case(1)
      this=that(1)
    case(2)
      this=sum(that)
    case default
      this=-1.0_wp
    end select

  end subroutine fio_density_adjust_spin_1_n

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin_2_n(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(that))
    case(1)
      this=0.5_wp*that(1)
    case(2)
      this=that
    case default
      this=-1.0_wp
    end select

  end subroutine fio_density_adjust_spin_2_n

  ! ---------------------------------------------------------
  pure subroutine fio_density_adjust_spin(this, that)
    real(kind=wp), dimension(:), intent(out) :: this
    real(kind=wp), dimension(:), intent(in)  :: that

    select case(size(this))
    case(1)
      call fio_density_adjust_spin_1_n(this(1), that)
    case(2)
      call fio_density_adjust_spin_2_n(this, that)
    case default
      this=-1.0_wp
    end select

  end subroutine fio_density_adjust_spin

  ! ---------------------------------------------------------
  subroutine fio_density_eval_1d(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: val

    real(kind=wp), allocatable, dimension(:) :: tvl
    real(kind=wp),              dimension(1) :: tv1
    integer                                  :: ierr, nspin

    PUSH_SUB(fio_density_eval_1d)

    call fio_density_get(this%self, nspin=nspin)
    if(nspin==1)then
      call fio_density_intrpl_eval(this, x, val, ierr)
      if(ierr/=FIO_DENSITY_INTRPL_OK)val=0.0_wp
    else
      SAFE_ALLOCATE(tvl(1:nspin))
      call fio_density_intrpl_eval(this, x, tvl, ierr)
      if(ierr/=FIO_DENSITY_INTRPL_OK)tvl=0.0_wp
      call fio_density_adjust_spin(tv1, tvl)
      SAFE_DEALLOCATE_A(tvl)
      val=tv1(1)
    end if

    POP_SUB(fio_density_eval_1d)
  end subroutine fio_density_eval_1d

  ! ---------------------------------------------------------
  subroutine fio_density_eval_md(this, x, val)
    type(fio_density_intrpl_t),  intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(out) :: val

    real(kind=wp), allocatable, dimension(:) :: tvl
    integer                                  :: ierr, nspin

    PUSH_SUB(fio_density_eval_md)

    call fio_density_get(this%self, nspin=nspin)
    if(nspin==size(val))then
      call fio_density_intrpl_eval(this, x, val, ierr)
      if(ierr/=FIO_DENSITY_INTRPL_OK)val=0.0_wp
    else
      SAFE_ALLOCATE(tvl(1:nspin))
      call fio_density_intrpl_eval(this, x, tvl, ierr)
      if(ierr/=FIO_DENSITY_INTRPL_OK)tvl=0.0_wp
      call fio_density_adjust_spin(val, tvl)
      SAFE_DEALLOCATE_A(tvl)
    end if

    POP_SUB(fio_density_eval_md)
  end subroutine fio_density_eval_md

#define TEMPLATE_NAME fio_density
#define INCLUDE_BODY
#include "intrpl_inc.F90"
#undef INCLUDE_BODY
#undef TEMPLATE_NAME

end module fio_density_m

!! Local Variables:
!! mode: f90
!! End:
