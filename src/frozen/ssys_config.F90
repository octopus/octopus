#include "global.h"

module ssys_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  use base_handle_oct_m
  use fio_config_oct_m
  use fio_handle_oct_m
  use frozen_config_oct_m
  use frozen_handle_oct_m
  use functional_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ssys_handle_oct_m
  use storage_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::            &
    ssys_config_use,   &
    ssys_config_add,   &
    ssys_config_parse

contains

  ! ---------------------------------------------------------
  elemental function ssys_calc_nrot(this) result(that)
    integer, intent(in) :: this

    integer :: that

    that = 1
    if(this>2) that = (this * (this - 1)) / 2

  end function ssys_calc_nrot

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_coordinates(this, ndim, block, line, icol, cols)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: line
    integer,             intent(in)  :: icol
    integer,             intent(in)  :: cols

    real(kind=wp), allocatable, dimension(:) :: array
    integer                                  :: iclm, irot

    PUSH_SUB(ssys_config_parse_coordinates)

    SAFE_ALLOCATE(array(ndim))
    array = 0.0_wp
    call json_init(this)
    call json_set(this, "dim", ndim)
    do iclm = 1, min(ndim, cols-icol)
      call parse_block_float(block, line, iclm+icol-1, array(iclm))
    end do
    array(1:ndim) = units_to_atomic(units_inp%length, array(1:ndim))
    call json_set(this, "r", array)
    SAFE_DEALLOCATE_A(array)
    irot = ssys_calc_nrot(ndim)
    SAFE_ALLOCATE(array(irot))
    array = 0.0_wp
    do iclm = 1, min(irot, cols-icol-ndim)
      call parse_block_float(block, line, iclm+icol+ndim-1, array(iclm))
    end do
    call json_set(this, "theta", array)
    SAFE_DEALLOCATE_A(array)

    POP_SUB(ssys_config_parse_coordinates)
  end subroutine ssys_config_parse_coordinates

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_subsystems_coordinates_block(this, ndim)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: ndim

    type(json_object_t),        pointer :: cnfg, pstn
    type(json_array_t),         pointer :: list
    type(block_t)                       :: block
    character(len=BASE_HANDLE_NAME_LEN) :: name
    integer                             :: ilin, nlin, ncls, ierr

    !%Variable SubSystemCoordinates
    !%Type block
    !%Section System::Subsystems
    !%Description
    !% Lists the name of the subsystem, the coordinates and the rotation to apply uppon reading.
    !% A subsystem can figure multiple times.
    !%
    !% <tt>%SubSystemCoordinates
    !% <br>&nbsp;&nbsp;'name_1' | x | y | z | o_xy | o_xz | o_yz
    !% <br>&nbsp;&nbsp;'name_2' | x | y | z | o_xy
    !% <br>&nbsp;&nbsp;'name_2' | x | y | z 
    !% <br>&nbsp;&nbsp;'name_3' | x
    !% <br>%</tt>
    !%
    !%End

    PUSH_SUB(ssys_config_parse_subsystems_coordinates_block)

    nullify(cnfg, list, pstn)
    if(parse_block('SubSystemCoordinates',block)==0) then
      nlin = parse_block_n(block)
      if(nlin>0)then
        do ilin = 1, nlin
          ncls=parse_block_cols(block, ilin-1)
          ASSERT(ncls>0)
          call parse_block_string(block, ilin-1, 0, name)
          call json_get(this, trim(adjustl(name)), cnfg, ierr)
          if(ierr/=JSON_OK)then
            message(1) = "Subsystem '"//trim(adjustl(name))//"' was not specified in the SubSystems block."
            call messages_fatal(1)
          end if
          call json_get(cnfg, "positions", list, ierr)
          ASSERT(ierr==JSON_OK)
          SAFE_ALLOCATE(pstn)
          call ssys_config_parse_coordinates(pstn, ndim, block, ilin-1, 1, ncls)
          call json_append(list, pstn)
          nullify(cnfg, list, pstn)
        end do
      end if
      call parse_block_end(block)
    end if

    POP_SUB(ssys_config_parse_subsystems_coordinates_block)
  end subroutine ssys_config_parse_subsystems_coordinates_block

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_subsystems_block(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim

    type(json_object_t),        pointer :: cnfg
    type(block_t)                       :: block
    character(len=BASE_HANDLE_NAME_LEN) :: name
    integer                             :: ilin, nlin, ncls, type

    !%Variable SubSystems
    !%Type block
    !%Section System::Subsystems
    !%Description
    !% Lists the name, the subsystem type, the directory and the optional parameters to be used on the subsystem calculation.
    !%
    !% <tt>%SubSystems
    !% <br>&nbsp;&nbsp;'name_1' | frozen | 'directory_1' | nearest
    !% <br>&nbsp;&nbsp;'name_2' | frozen | 'directory_2' | nearest
    !% <br>&nbsp;&nbsp;'name_3' | frozen | 'directory_2' | nearest
    !% <br>%</tt>
    !%
    !% For frozen subsystems the optional parameter is the type of interpolation to use.
    !%
    !%Option frozen 2
    !% Frozen subsystems.
    !%Option nearest 1
    !% Uses the same value as the nearest grid point.
    !%Option qshep 2
    !% Uses Quadratic Shepard interpolation.
    !% 
    !%End

    PUSH_SUB(ssys_config_parse_subsystems_block)

    nullify(cnfg)
    call json_init(this)
    if(parse_block('SubSystems',block)==0) then
      nlin = parse_block_n(block)
      if(nlin>0)then
        do ilin = 1, nlin
          ncls=parse_block_cols(block, ilin-1)
          ASSERT(ncls>0)
          if(ncls<2)then
            message(1) = "In the SubSystems block at least subsystem name and type must be specified."
            call messages_fatal(1)
          end if
          call parse_block_string(block, ilin-1, 0, name)
          call parse_block_integer(block, ilin-1, 1, type)
          SAFE_ALLOCATE(cnfg)
          select case(type)
          case(HNDL_TYPE_FRZN)
            call fio_config_parse(cnfg, block, ilin-1, 2, ncls)
          case default
            message(1)="Unknown subsystems type."
            call messages_fatal(1)
          end select
          call json_set(cnfg, "name", trim(adjustl(name)))
          call json_set(this, trim(adjustl(name)), cnfg)
          nullify(cnfg)
        end do
      end if
      call parse_block_end(block)
      call ssys_config_parse_subsystems_coordinates_block(this, ndim)
    end if

    POP_SUB(ssys_config_parse_subsystems_block)
  end subroutine ssys_config_parse_subsystems_block

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_systems(this, nspin, ndim)
    type(json_array_t), intent(inout) :: this
    integer,            intent(in)    :: nspin
    integer,            intent(in)    :: ndim

    type(json_object_t)                 :: dict
    type(json_object_t),        pointer :: cnfg, ocfg, icfg
    type(json_array_t),         pointer :: frzn
    type(json_object_iterator_t)        :: iter
    character(len=BASE_HANDLE_NAME_LEN) :: name
    integer                             :: type, ierr

    PUSH_SUB(ssys_config_parse_systems)

    nullify(cnfg, ocfg, icfg, frzn)
    call ssys_config_parse_subsystems_block(dict, ndim)
    call json_init(iter, dict)
    do
      nullify(ocfg, icfg)
      call json_next(iter, name, icfg, ierr)
      if(ierr/=JSON_OK)exit
      SAFE_ALLOCATE(ocfg)
      call json_copy(ocfg, icfg)
      call json_get(ocfg, "type", type, ierr)
      ASSERT(ierr==JSON_OK)
      select case(type)
      case(HNDL_TYPE_FNIO)
        if(.not.associated(frzn))then
          SAFE_ALLOCATE(cnfg)
          call frozen_config_parse(cnfg, nspin, ndim)
          call json_append(this, cnfg)
          call json_get(cnfg, "systems", frzn, ierr)
          ASSERT(ierr==JSON_OK)
          nullify(cnfg)
        end if
        call json_append(frzn, ocfg)
      case default
        message(1) = "Unknown subsystems type."
        call messages_fatal(1)
      end select
    end do
    call json_end(iter)
    nullify(ocfg, icfg)
    call json_end(dict)

    POP_SUB(ssys_config_parse_systems)
  end subroutine ssys_config_parse_systems

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_external(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_external)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_external)
  end subroutine ssys_config_parse_external

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_ionic(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(ssys_config_parse_ionic)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_TERM)

    POP_SUB(ssys_config_parse_ionic)
  end subroutine ssys_config_parse_ionic

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_functional(this, id, factor, polarized)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: id
    real(kind=wp),       intent(in)  :: factor
    logical,             intent(in)  :: polarized

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_functional)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_FNCT)
    call json_set(this, "factor", factor)
    call json_set(this, "spin", polarized)
    SAFE_ALLOCATE(cnfg)
    call functional_init(cnfg, id=id)
    call json_set(this, "functional", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    if(id>FUNCT_XC_NONE)then
      call storage_init(cnfg, full=.false.)
    else
      call storage_init(cnfg, full=.false., allocate=.false.)
    end if
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_functional)
  end subroutine ssys_config_parse_functional

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_tnadd(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
 
    type(json_object_t), pointer :: cnfg
    real(kind=wp)                :: factor
    integer                      :: id
    logical                      :: plrz

    PUSH_SUB(ssys_config_parse_tnadd)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_HMLT)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false., allocate=.true.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)
    call parse_variable('TnaddFunctional', FUNCT_XC_NONE, id)
    !%Variable TnaddFactor
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian::Subsystems
    !%Description
    !% Chooses the Kinetic Functional amplification factor.
    !%End
    call parse_variable('TnaddFactor', 1.0_wp, factor)
    if(abs(factor)<1.0e-7_wp)then
      message(1) = "The 'TnaddFactor' value specified may be too small."
      call messages_warning(1)
    end if
    !%Variable TnaddPolarized
    !%Type logical
    !%Default yes
    !%Section Hamiltonian::Subsystems
    !%Description
    !% Calculates the Kinetic Functional with spin polarization or not.
    !%End
    call parse_variable('TnaddPolarized', (nspin>1), plrz)
    call json_set(this, "spin", plrz)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_functional(cnfg, id, factor, plrz)
    call json_set(this, "total", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_functional(cnfg, id, factor, plrz)
    call json_set(cnfg, "system", "live")
    call json_set(this, "live", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_tnadd)
  end subroutine ssys_config_parse_tnadd

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_hamiltonian(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(ssys_config_parse_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_external(cnfg)
    call json_set(this, "external", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_ionic(cnfg)
    call json_set(this, "ionic", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call ssys_config_parse_tnadd(cnfg, nspin)
    call json_set(this, "tnadd", cnfg)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_hamiltonian)
  end subroutine ssys_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine ssys_config_parse_model(this, nspin)
    type(json_object_t), intent(inout) :: this
    integer,             intent(in)    :: nspin

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse_model)

    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_hamiltonian(cnfg, nspin)
    nullify(cnfg)

    POP_SUB(ssys_config_parse_model)
  end subroutine ssys_config_parse_model

  ! ---------------------------------------------------------
  function ssys_config_use() result(that)

    logical :: that

    PUSH_SUB(ssys_config_use)

    that = parse_is_defined("SubSystems")

    POP_SUB(ssys_config_use)
  end function ssys_config_use

  ! ---------------------------------------------------------
  subroutine ssys_config_add(this, that)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that

    type(json_array_t), pointer :: list
    integer                     :: ierr

    PUSH_SUB(ssys_config_add)

    nullify(list)
    call json_get(this, "systems", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_append(list, that)
    nullify(list)

    POP_SUB(ssys_config_add)
  end subroutine ssys_config_add

  ! ---------------------------------------------------------
  subroutine ssys_config_parse(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin
    integer,             intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    integer                      :: ierr

    PUSH_SUB(ssys_config_parse)

    nullify(cnfg, list)
    call base_config_parse(this, nspin, ndim)
    call json_set(this, "type", HNDL_TYPE_SSYS)
    call json_set(this, "name", "main")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_model(cnfg, nspin)
    nullify(cnfg)
    call json_get(this, "systems", list, ierr)
    ASSERT(ierr==JSON_OK)
    call ssys_config_parse_systems(list, nspin, ndim)
    nullify(cnfg, list)

    POP_SUB(ssys_config_parse)
  end subroutine ssys_config_parse

end module ssys_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
