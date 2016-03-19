#include "global.h"

module fio_config_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m
  use curvilinear_oct_m
  use fio_handle_oct_m
  use global_oct_m
  use intrpl_oct_m
  use io_oct_m
  use json_oct_m
  use json_parser_oct_m
  use kinds_oct_m
  use loct_oct_m
  use messages_oct_m
  use parser_oct_m
  use path_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::              &
    input_frozen_dir,    &
    input_frozen_config

  public ::           &
    fio_config_parse

  interface fio_config_parse
    module procedure fio_config_parse_dir
    module procedure fio_config_parse_block
  end interface fio_config_parse

  character(len=*), parameter :: input_static_dir = STATIC_DIR

  character(len=*), parameter :: input_frozen_dir    = "frozen"
  character(len=*), parameter :: input_frozen_config = "config.json"

contains

  ! ---------------------------------------------------------
  subroutine fio_config_parse_get_dir(topdir, reldir, outdir)
    character(len=*), intent(in)  :: topdir
    character(len=*), intent(in)  :: reldir
    character(len=*), intent(out) :: outdir

    character(len=max_path_len) :: apth, bpth
    logical                     :: exst

    PUSH_SUB(fio_config_parse_get_dir)

    call path_realpath(reldir, outdir)
    if(.not.path_isabs(reldir))then
       call path_realpath(topdir, apth)
       call path_join(apth, reldir, bpth)
       call path_realpath(bpth, outdir)
       exst = loct_dir_exists(outdir)
       if(.not.exst)then
          call path_getcwd(apth)
          call path_join(apth, reldir, bpth)
          call path_realpath(bpth, outdir)
          exst = loct_dir_exists(outdir)
          if(.not.exst) call path_realpath(topdir, outdir)
       end if
    end if

    POP_SUB(fio_config_parse_get_dir)
  end subroutine fio_config_parse_get_dir

  ! ---------------------------------------------------------
  subroutine fio_config_parse_get_file(topdir, reldir, outdir, file, ierr)
    character(len=*), intent(in)  :: topdir
    character(len=*), intent(in)  :: reldir
    character(len=*), intent(out) :: outdir
    character(len=*), intent(in)  :: file
    integer,          intent(out) :: ierr

    character(len=max_path_len) :: fpth
    logical                     :: exst

    PUSH_SUB(fio_config_parse_get_file)

    ierr = 0
    call fio_config_parse_get_dir(topdir, reldir, outdir)
    call path_join(outdir, file, fpth)
    inquire(file=trim(adjustl(fpth)), exist=exst)
    if(.not.exst)then
      call path_realpath(topdir, outdir)
      call path_join(outdir, file, fpth)
      inquire(file=trim(adjustl(fpth)), exist=exst)
      if(.not.exst) ierr = -1
    end if

    POP_SUB(fio_config_parse_get_file)
  end subroutine fio_config_parse_get_file

  ! ---------------------------------------------------------
  subroutine fio_config_parse_storage(this, allocate)
    type(json_object_t), intent(inout) :: this
    logical,   optional, intent(in)    :: allocate

    PUSH_SUB(fio_config_parse_storage)

    if(present(allocate)) call json_set(this, "allocate", allocate)

    POP_SUB(fio_config_parse_storage)
  end subroutine fio_config_parse_storage

  ! ---------------------------------------------------------
  subroutine fio_config_parse_simul_box(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: ierr

    PUSH_SUB(fio_config_parse_simul_box)

    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_get_file(dirname, idir, odir, file, ierr)
    if(ierr/=0)then
      message(1) = "Could not find the input file: '"//trim(adjustl(file))//"',"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
      call messages_fatal(2)
    end if
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))

    POP_SUB(fio_config_parse_simul_box)
  end subroutine fio_config_parse_simul_box

  ! ---------------------------------------------------------
  subroutine fio_config_parse_index(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    character(len=MAX_PATH_LEN)  :: idir, odir, file
    integer :: ierr

    PUSH_SUB(fio_config_parse_index)

    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_get_file(dirname, idir, odir, file, ierr)
    if(ierr/=0)then
      message(1) = "Could not find the input file: '"//trim(adjustl(file))//"',"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
      call messages_fatal(2)
    end if
    call fio_config_parse_get_file(dirname, idir, odir, "lxyz.obf", ierr)
    if(ierr/=0)then
      message(1) = "Could not find the input file: 'lxyz.obf',"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
      call messages_fatal(2)
    end if
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))

    POP_SUB(fio_config_parse_index)
  end subroutine fio_config_parse_index

  ! ---------------------------------------------------------
  subroutine fio_config_parse_mesh(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: idir, odir, file
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_mesh)

    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_get_file(dirname, idir, odir, file, ierr)
    if(ierr/=0)then
      message(1) = "Could not find the input file: '"//trim(adjustl(file))//"',"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
      call messages_fatal(2)
    end if
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    nullify(cnfg)
    call json_get(this, "index", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_index(cnfg, dirname)
    nullify(cnfg)

    POP_SUB(fio_config_parse_mesh)
  end subroutine fio_config_parse_mesh

  ! ---------------------------------------------------------
  subroutine fio_config_parse_grid(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_grid)

    nullify(cnfg)
    call json_get(this, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simul_box(cnfg, dirname)
    nullify(cnfg)
    call json_get(this, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_mesh(cnfg, dirname)
    nullify(cnfg)

    POP_SUB(fio_config_parse_grid)
  end subroutine fio_config_parse_grid

  ! ---------------------------------------------------------
  subroutine fio_config_parse_simulation(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_simulation)

    nullify(cnfg)
    call json_get(this, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_grid(cnfg, dirname)
    nullify(cnfg)

    POP_SUB(fio_config_parse_simulation)
  end subroutine fio_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine fio_config_parse_density(this, dirname, usedensity, invertpol)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity
    logical,   optional, intent(in)    :: invertpol

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    character(len=MAX_PATH_LEN)  :: idir, odir, file
    integer                      :: indx, nspin, ierr

    PUSH_SUB(fio_config_parse_density)

    nullify(cnfg, list)
    call json_get(this, "nspin", nspin, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "files", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)==nspin)
    do indx = 1, nspin
      call json_get(list, indx, file, ierr)
      ASSERT(ierr==JSON_OK)
      file = trim(adjustl(file))//".obf"
      call fio_config_parse_get_file(dirname, idir, odir, file, ierr)
      if(ierr/=0)then
        message(1) = "Could not find the input file: '"//trim(adjustl(file))//"',"
        message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
        call messages_fatal(2)
      end if
      call json_set(list, indx, trim(adjustl(file)), ierr)
      ASSERT(ierr==JSON_OK)
    end do
    call json_set(this, "dir", trim(adjustl(odir)))
    nullify(list)
    call json_set(this, "invertpolarization", invertpol)
    call json_get(this, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_storage(cnfg, usedensity)
    nullify(cnfg)

    POP_SUB(fio_config_parse_density)
  end subroutine fio_config_parse_density

  ! ---------------------------------------------------------
  subroutine fio_config_parse_states(this, dirname, usedensity, invertpol)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity
    logical,   optional, intent(in)    :: invertpol

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_states)

    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_density(cnfg, dirname, usedensity, invertpol)
    nullify(cnfg)

    POP_SUB(fio_config_parse_states)
  end subroutine fio_config_parse_states

  ! ---------------------------------------------------------
  subroutine fio_config_parse_system(this, dirname, usedensity, invertpol)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity
    logical,   optional, intent(in)    :: invertpol

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_system)

    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_states(cnfg, dirname, usedensity, invertpol)
    nullify(cnfg)

    POP_SUB(fio_config_parse_system)
  end subroutine fio_config_parse_system

  ! ---------------------------------------------------------
  subroutine fio_config_parse_external(this, dirname, usepotential)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usepotential

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: idir, odir, file
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_external)

    nullify(cnfg)
    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    file = trim(adjustl(file))//".obf"
    call fio_config_parse_get_file(dirname, idir, odir, file, ierr)
    if(ierr/=0)then
      message(1) = "Could not find the input file: '"//trim(adjustl(file))//"',"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"',"
      call messages_fatal(2)
    end if
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    call json_get(this, "storage", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_storage(cnfg, usepotential)
    nullify(cnfg)

    POP_SUB(fio_config_parse_external)
  end subroutine fio_config_parse_external

  ! ---------------------------------------------------------
  subroutine fio_config_parse_hamiltonian(this, dirname, usepotential)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usepotential

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_hamiltonian)

    nullify(cnfg)
    call json_get(this, "external", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_external(cnfg, dirname, usepotential)
    nullify(cnfg)

    POP_SUB(fio_config_parse_hamiltonian)
  end subroutine fio_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine fio_config_parse_model(this, dirname, usepotential, usedensity, invertpol)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usepotential
    logical,   optional, intent(in)    :: usedensity
    logical,   optional, intent(in)    :: invertpol

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_model)

    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_system(cnfg, dirname, usedensity, invertpol)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_hamiltonian(cnfg, dirname, usepotential)
    nullify(cnfg)

    POP_SUB(fio_config_parse_model)
  end subroutine fio_config_parse_model

  ! ---------------------------------------------------------
  subroutine fio_config_parse_read(this, dirname, outdir)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    character(len=*),    intent(out) :: outdir

    type(json_parser_t)         :: parser
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: iunit, ierr

    PUSH_SUB(fio_config_parse_read)

    call fio_config_parse_get_file(dirname, ".", outdir, input_frozen_config, ierr)
    if(ierr/=0)then
      call fio_config_parse_get_file(dirname, input_frozen_dir, outdir, input_frozen_config, ierr)
      if(ierr/=0)then
        call path_join(input_static_dir, input_frozen_dir, dir)
        call fio_config_parse_get_file(dirname, dir, outdir, input_frozen_config, ierr)
        if(ierr/=0)then
          message(1) = "Could not find the file: '"//trim(adjustl(input_frozen_config))//"'"
          call messages_fatal(1)
        end if
      end if
    end if
    call path_join(outdir, input_frozen_config, file)
    iunit = io_open(file, action='read', status="old")
    ASSERT(iunit>0)
    call json_parser_init(parser, iunit, ierr)
    if(ierr/=0) call json_parser_error(parser, fatal=.true.)
    call json_parser_parse(parser, this, ierr)
    if((ierr/=0).or.(.not.json_isdef(this)))&
      call json_parser_error(parser, fatal=.true.)
    call json_parser_end(parser)
    call io_close(iunit)

    POP_SUB(fio_config_parse_read)
  end subroutine fio_config_parse_read

  ! ---------------------------------------------------------
  subroutine fio_config_parse_dir(this, dirname, interpolation, usepotential, usedensity, invertpol)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    integer,   optional, intent(in)  :: interpolation
    logical,   optional, intent(in)  :: usepotential
    logical,   optional, intent(in)  :: usedensity
    logical,   optional, intent(in)  :: invertpol

    character(len=MAX_PATH_LEN)  :: outdir
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_dir)

    nullify(cnfg)
    call fio_config_parse_read(this, dirname, outdir)
    call json_set(this, "type", HNDL_TYPE_FNIO)
    call json_set(this, "name", "fio")
    if(present(interpolation)) call json_set(this, "interpolation", interpolation)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simulation(cnfg, outdir)
    nullify(cnfg)
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_model(cnfg, outdir, usepotential, usedensity, invertpol)
    nullify(cnfg)

    POP_SUB(fio_config_parse_dir)
  end subroutine fio_config_parse_dir

  ! ---------------------------------------------------------
  subroutine fio_config_parse_block(this, block, line, icol, ncol)
    type(json_object_t), intent(out) :: this
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: line
    integer,             intent(in)  :: icol
    integer,             intent(in)  :: ncol

    character(len=MAX_PATH_LEN)  :: dirname
    integer                      :: intrp
    logical                      :: usepot, userho, invpol

    PUSH_SUB(fio_config_parse_block)

    call parse_block_string(block, line, icol, dirname)
    intrp = NEAREST
    usepot = .true.
    userho = .true.
    invpol = .false.
    if(ncol>(icol+1))call parse_block_integer(block, line, icol+1, intrp)
    if(ncol>(icol+2))call parse_block_logical(block, line, icol+2, usepot)
    if(ncol>(icol+3))call parse_block_logical(block, line, icol+3, userho)
    if(ncol>(icol+4))call parse_block_logical(block, line, icol+4, invpol)
    call fio_config_parse(this, dirname, intrp, usepot, userho, invpol)

    POP_SUB(fio_config_parse_block)
  end subroutine fio_config_parse_block

end module fio_config_oct_m

!! Local Variables:
!! mode: f90
!! End:
