#include "global.h"

module base_config_oct_m

  use base_hamiltonian_oct_m
  use base_handle_oct_m
  use functional_oct_m
  use global_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use storage_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::               &
    BASE_CONFIG_NAME_LEN

  public ::     &
    PARSE_OK,   &
    PARSE_FAIL
  
  public ::                    &
    base_config_use,           &
    base_config_add,           &
    base_config_parse,         &
    base_config_parse_kinetic

  integer, parameter :: BASE_CONFIG_NAME_LEN = BASE_HANDLE_NAME_LEN
  
  integer, parameter :: PARSE_OK   =  0
  integer, parameter :: PARSE_FAIL = -1

  integer, parameter :: default_ndim  = 3
  integer, parameter :: default_nspin = 1

  interface base_config_parse
    module procedure base_config_parse_type
    module procedure base_config_parse_pass
  end interface base_config_parse

contains

  ! ---------------------------------------------------------
  function base_config_use() result(that)

    logical :: that

    PUSH_SUB(base_config_use)

    that = parse_is_defined("SubSystems")

    POP_SUB(base_config_use)
  end function base_config_use

  ! ---------------------------------------------------------
  subroutine base_config_add(this, name, that)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: name
    type(json_object_t), intent(in)    :: that

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(base_config_add)

    nullify(cnfg)
    call json_get(this, "subsystems", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_set(cnfg, name, that)
    nullify(cnfg)

    POP_SUB(base_config_add)
  end subroutine base_config_add

  ! ---------------------------------------------------------
  subroutine base_config_parse_space(this, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    integer :: idim

    PUSH_SUB(base_config_parse_space)

    idim = default_ndim
    if(present(ndim)) idim = ndim
    ASSERT(idim>0)
    call json_init(this)
    call json_set(this, "dimensions", idim)

    POP_SUB(base_config_parse_space)
  end subroutine base_config_parse_space

  ! ---------------------------------------------------------
  subroutine base_config_parse_geometry(this)
    type(json_object_t), intent(out) :: this

    PUSH_SUB(base_config_parse_geometry)

    call json_init(this)
    call json_set(this, "default", .true.)
    call json_set(this, "reduce", .false.)

    POP_SUB(base_config_parse_geometry)
  end subroutine base_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine base_config_parse_density(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    type(json_object_t),             pointer :: cnfg
    real(kind=wp), dimension(:), allocatable :: chrg
    integer                                  :: ispin

    PUSH_SUB(base_config_parse_density)

    nullify(cnfg)
    ispin = default_nspin
    if(present(nspin)) ispin = nspin
    ASSERT(ispin>0)
    SAFE_ALLOCATE(chrg(1:ispin))
    chrg = 0.0_wp
    call json_init(this)
    call json_set(this, "reduce", .false.)
    call json_set(this, "external", .false.)
    call json_set(this, "default", .true.)
    call json_set(this, "nspin", ispin)
    call json_set(this, "charge", chrg)
    SAFE_DEALLOCATE_A(chrg)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, ndim=ispin, fine=.true.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_density)
  end subroutine base_config_parse_density

  ! ---------------------------------------------------------
  subroutine base_config_parse_states(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_states)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_density(cnfg, nspin)
    call json_set(this, "density", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_states)
  end subroutine base_config_parse_states

  ! ---------------------------------------------------------
  subroutine base_config_parse_system(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_system)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_space(cnfg, ndim)
    call json_set(this, "space", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_geometry(cnfg)
    call json_set(this, "geometry", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_states(cnfg, nspin)
    call json_set(this, "states", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_system)
  end subroutine base_config_parse_system

  ! ---------------------------------------------------------
  subroutine base_config_parse_kinetic(this, nspin)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: nspin

    type(json_object_t), pointer :: cnfg
    real(kind=wp)                :: factor
    integer                      :: id
    logical                      :: plrz


    PUSH_SUB(base_config_parse_kinetic)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", TERM_TYPE_FNCT)
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
    call json_set(this, "factor", factor)
    !%Variable TnaddPolarized
    !%Type logical
    !%Default yes
    !%Section Hamiltonian::Subsystems
    !%Description
    !% Calculates the Kinetic Functional with spin polarization or not.
    !%End
    call parse_variable('TnaddPolarized', (nspin>1), plrz)
    call json_set(this, "spin", plrz)
    call parse_variable('TnaddFunctional', FUNCT_XC_NONE, id)
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

    POP_SUB(base_config_parse_kinetic)
  end subroutine base_config_parse_kinetic

  ! ---------------------------------------------------------
  subroutine base_config_parse_hamiltonian(this)
    type(json_object_t), intent(out) :: this

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_hamiltonian)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", TERM_TYPE_HMLT)
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false., allocate=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_hamiltonian)
  end subroutine base_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine base_config_parse_model(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_model)

    nullify(cnfg)
    call json_init(this)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_system(cnfg, nspin, ndim)
    call json_set(this, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_hamiltonian(cnfg)
    call json_set(this, "hamiltonian", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_model)
  end subroutine base_config_parse_model

  ! ---------------------------------------------------------
  elemental function base_config_calc_nrot(this) result(that)
    integer, intent(in) :: this

    integer :: that

    that = 1
    if(this>2) that = (this * (this - 1)) / 2

  end function base_config_calc_nrot

  ! ---------------------------------------------------------
  subroutine base_config_parse_position(this, ndim, block, line, icol, cols)
    type(json_object_t), intent(out) :: this
    integer,             intent(in)  :: ndim
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: line
    integer,             intent(in)  :: icol
    integer,             intent(in)  :: cols

    real(kind=wp), allocatable, dimension(:) :: array
    integer                                  :: iclm, irot

    PUSH_SUB(base_config_parse_position)

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
    irot = base_config_calc_nrot(ndim)
    SAFE_ALLOCATE(array(irot))
    array = 0.0_wp
    do iclm = 1, min(irot, cols-icol-ndim)
      call parse_block_float(block, line, iclm+icol+ndim-1, array(iclm))
    end do
    call json_set(this, "theta", array)
    SAFE_DEALLOCATE_A(array)

    POP_SUB(base_config_parse_position)
  end subroutine base_config_parse_position

  ! ---------------------------------------------------------
  subroutine base_config_parse_positions(this, name, ndim)
    type(json_array_t), intent(out) :: this
    character(len=*),   intent(in)  :: name
    integer,  optional, intent(in)  :: ndim

    character(len=BASE_CONFIG_NAME_LEN) :: inam
    type(json_object_t),        pointer :: cnfg
    type(block_t)                       :: block
    integer                             :: idim, ilin, nlin, ncls

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

    PUSH_SUB(base_config_parse_positions)

    nullify(cnfg)
    idim = default_ndim
    if(present(ndim)) idim = ndim
    ASSERT(idim>0)
    call json_init(this)
    if(parse_block('SubSystemCoordinates',block)==0) then
      nlin = parse_block_n(block)
      if(nlin>0)then
        do ilin = 1, nlin
          ncls=parse_block_cols(block, ilin-1)
          ASSERT(ncls>0)
          call parse_block_string(block, ilin-1, 0, inam)
          if(trim(adjustl(name))/=trim(adjustl(inam))) cycle
          SAFE_ALLOCATE(cnfg)
          call base_config_parse_position(cnfg, idim, block, ilin-1, 1, ncls)
          call json_append(this, cnfg)
          nullify(cnfg)
        end do
      end if
      call parse_block_end(block)
    end if

    POP_SUB(base_config_parse_positions)
  end subroutine base_config_parse_positions

  ! ---------------------------------------------------------
  subroutine base_config_parse_subsystems(this, ndim, subs_config_parse)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: ndim

    interface
      subroutine subs_config_parse(this, type, block, line, ierr)
        use json_oct_m
        use parser_oct_m
        type(json_object_t), intent(out) :: this
        integer,             intent(in)  :: type
        type(block_t),       intent(in)  :: block
        integer,             intent(in)  :: line
        integer,             intent(out) :: ierr
      end subroutine subs_config_parse
    end interface
    
    type(json_object_t),        pointer :: cnfg
    type(json_array_t),         pointer :: list
    type(block_t)                       :: block
    character(len=BASE_HANDLE_NAME_LEN) :: name
    integer                             :: ilin, nlin, ncls, type, ierr

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

    PUSH_SUB(base_config_parse_subsystems)

    nullify(cnfg, list)
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
          call subs_config_parse(cnfg, type, block, ilin, ierr)
          if(ierr==PARSE_OK)then
            ASSERT(json_isdef(cnfg))
            call json_set(cnfg, "name", trim(adjustl(name)))
            SAFE_ALLOCATE(list)
            call base_config_parse_positions(list, trim(adjustl(name)), ndim)
            if(json_len(list)>0)then
              call json_set(cnfg, "positions", list)
            else
              call json_end(list)
              SAFE_DEALLOCATE_P(list)
            end if
            nullify(list)
            call json_set(this, trim(adjustl(name)), cnfg)
          else
            call json_end(cnfg)
            SAFE_DEALLOCATE_P(cnfg)
          end if
          nullify(cnfg)
        end do
      end if
      call parse_block_end(block)
    end if

    POP_SUB(base_config_parse_subsystems)
  end subroutine base_config_parse_subsystems

  ! ---------------------------------------------------------
  subroutine base_config_parse_type(this, nspin, ndim)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(base_config_parse_type)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HNDL_TYPE_NONE)
    call json_set(this, "name", "base")
    SAFE_ALLOCATE(cnfg)
    call simulation_init(cnfg)
    call json_set(this, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_model(cnfg, nspin, ndim)
    call json_set(this, "model", cnfg)
    nullify(cnfg)

    POP_SUB(base_config_parse_type)
  end subroutine base_config_parse_type

  ! ---------------------------------------------------------
  subroutine base_config_parse_pass(this, nspin, ndim, subs_config_parse)
    type(json_object_t), intent(out) :: this
    integer,   optional, intent(in)  :: nspin
    integer,   optional, intent(in)  :: ndim

    interface
      subroutine subs_config_parse(this, type, block, line, ierr)
        use json_oct_m
        use parser_oct_m
        type(json_object_t), intent(out) :: this
        integer,             intent(in)  :: type
        type(block_t),       intent(in)  :: block
        integer,             intent(in)  :: line
        integer,             intent(out) :: ierr
      end subroutine subs_config_parse
    end interface
    
    type(json_object_t), pointer :: cnfg
    
    PUSH_SUB(base_config_parse_pass)

    nullify(cnfg)
    call base_config_parse(this, nspin, ndim)
    SAFE_ALLOCATE(cnfg)
    call base_config_parse_subsystems(cnfg, ndim, subs_config_parse)
    if(json_len(cnfg)>0)then
      call json_set(this, "subsystems", cnfg)
    else
      call json_end(cnfg)
      SAFE_DEALLOCATE_P(cnfg)
    end if
    nullify(cnfg)

    POP_SUB(base_config_parse_pass)
  end subroutine base_config_parse_pass

end module base_config_oct_m

!! Local Variables:
!! mode: f90
!! End:

