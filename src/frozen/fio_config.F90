#include "global.h"

#define TEMPLATE_NAME fio
#define EXTERNAL
#include "thamiltonian.F90"
#include "tcalc.F90"
#include "twrap.F90"
#undef EXTERNAL
#undef TEMPLATE_NAME

module fio_config_m

  use global_m
  use messages_m
  use profiling_m

  use interpolation_m, only: NEAREST
  use io_m, only: io_open, io_close
  use json_parser_m, only: json_parser_t, json_parser_init, json_parser_end, &
    json_parser_parse, json_parser_error
  use parser_m, only: block_t, parse_block_string, parse_block_logical, parse_block_integer, parse_block_float
  use unit_m, only: units_to_atomic
  use unit_system_m, only: units_inp

  use basis_m, only: basis_t, basis_init, basis_end, basis_copy
  use json_m,  only: JSON_OK, json_object_t, json_isdef, json_init, json_set, json_get, json_write
  use kinds_m, only: wp

  implicit none

  private
  public ::           &
    fio_config_parse

contains

  ! ---------------------------------------------------------
  subroutine fio_config_basis_parse(this, block, ndim, line, cols)
    type(json_object_t), intent(out) :: this
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: line
    integer,             intent(in)  :: cols
    !
    real(kind=wp), dimension(ndim) :: array
    integer                        :: icol, col
    !
    icol=3
    array=0.0_wp
    call json_init(this)
    call json_set(this, "dim", ndim)
    do col=icol, min(ndim+icol-1, cols-1)
      call parse_block_float(block, line-1, col, array(col-icol+1))
    end do
    array=units_to_atomic(units_inp%length, array)
    call json_set(this, "r", array)
    array=0.0_wp
    do col=ndim+icol, min(2*ndim+icol-1, cols-1)
      call parse_block_float(block, line-1, col, array(col-ndim-icol+1))
    end do
    call json_set(this, "theta", array)
    return
  end subroutine fio_config_basis_parse

  ! ---------------------------------------------------------
  subroutine fio_config_calc_parse(this, file)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: file
    !
    type(json_parser_t) :: parser
    integer             :: iunit, ierr
    !
    call json_init(this)
    iunit=io_open(trim(file), action='read', status="old", is_tmp=.true.)
    if(iunit>0)then
      call json_parser_init(parser, iunit, ierr)
      if(ierr/=0)call json_parser_error(parser, fatal=.true.)
      call json_parser_parse(parser, this, ierr)
      if((ierr/=0).or.(.not.json_isdef(this)))call json_parser_error(parser, fatal=.true.)
      call json_parser_end(parser)
      call io_close(iunit)
    else
      message(1)="Could not open the config file: '"//trim(file)//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(2)
    end if
    return
  end subroutine fio_config_calc_parse

  ! ---------------------------------------------------------
  subroutine fio_config_calc_set_dir(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    character(len=MAX_PATH_LEN) :: dir
    integer                     :: ierr
    !
    call json_get(this, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    dir=trim(adjustl(dirname))//trim(adjustl(dir))
    call json_set(this, "dir", trim(dir))
    return
  end subroutine fio_config_calc_set_dir

  ! ---------------------------------------------------------
  subroutine fio_config_calc_update(this, type, usept, dirname)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: type
    logical,             intent(in)    :: usept
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: l1, l2, l3
    integer                      :: ierr
    !
    nullify(l1, l2, l3)
    call json_get(this, "simulation", l1, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(l1, "grid", l2, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(l2, "simul_box", l3, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_calc_set_dir(l3, dirname)
    nullify(l3)
    call json_get(l2, "mesh", l3, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_calc_set_dir(l3, dirname)
    nullify(l1, l2, l3)
    call json_get(this, "system", l1, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(l1, "states", l2, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(l2, "density", l3, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(l3, "interpolation", type)
    call fio_config_calc_set_dir(l3, dirname)
    nullify(l1, l2, l3)
    call json_get(this, "hamiltonian", l1, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(l1, "external", l2, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(l2, "interpolation", type)
    call json_set(l2, "allocate", usept)
    call fio_config_calc_set_dir(l2, dirname)
    nullify(l1, l2)
    return
  end subroutine fio_config_calc_update

  ! ---------------------------------------------------------
  subroutine fio_config_parse(this, block, ndim, line, cols)
    type(json_object_t), intent(out) :: this
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: ndim
    integer,             intent(in)  :: line
    integer,             intent(in)  :: cols
    !
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dirname
    logical                      :: usept
    integer                      :: icol, iunit, intrp, ierr
    !
    icol=0
    nullify(cnfg)
    call json_init(this)
    call json_set(this, "temporary", .true.)
    call json_set(this, "pass_simulation", .false.)
    call parse_block_string(block, line-1, icol, dirname)
    dirname=trim(adjustl(dirname))//"/"//STATIC_DIR
    call json_set(this, "dir", trim(dirname))
    intrp=NEAREST
    if(cols>(icol+1))call parse_block_integer(block, line-1, icol+1, intrp)
    call json_set(this, "interpolation", intrp)
    usept=.false.
    if(cols>(icol+2))call parse_block_logical(block, line-1, icol+2, usept)
    call json_set(this, "use_potential", usept)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_config_basis_parse(cnfg, block, ndim, line, cols)
    call json_set(this, "basis", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_config_calc_parse(cnfg, trim(adjustl(dirname))//"config.json")
    call fio_config_calc_update(cnfg, intrp, usept, trim(adjustl(dirname)))
    call json_set(this, "calculation", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_config_parse

end module fio_config_m

!! Local Variables:
!! mode: f90
!! End:
