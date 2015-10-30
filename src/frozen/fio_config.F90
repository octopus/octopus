#include "global.h"

module fio_config_m

  use base_config_m
  use base_hamiltonian_m
  use curvilinear_m
  use fio_handle_m
  use global_m
  use intrpl_m
  use io_m
  use json_m
  use json_parser_m
  use kinds_m
  use loct_m
  use messages_m
  use parser_m
  use path_m
  use profiling_m
  use species_m

  implicit none

  private

  public ::           &
    fio_config_parse

  interface fio_config_parse
    module procedure fio_config_parse_dir
    module procedure fio_config_parse_block
  end interface fio_config_parse

  character(len=*), public, parameter :: input_mesh_dir      = "./restart/"//GS_DIR
  character(len=*), public, parameter :: input_mesh_file     = "mesh"
  character(len=*), public, parameter :: input_static_dir    = "./"//STATIC_DIR
  character(len=*), public, parameter :: input_external_file = "v0.obf"

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
  subroutine fio_config_parse_get_file(topdir, reldir, outdir, file)
    character(len=*), intent(in)  :: topdir
    character(len=*), intent(in)  :: reldir
    character(len=*), intent(out) :: outdir
    character(len=*), intent(in)  :: file

    character(len=max_path_len) :: fpth
    logical                     :: exst

    PUSH_SUB(fio_config_parse_get_file)

    call fio_config_parse_get_dir(topdir, reldir, outdir)
    call path_join(outdir, file, fpth)
    inquire(file=trim(adjustl(fpth)), exist=exst)
    if(.not.exst)then
      call path_realpath(topdir, outdir)
      call path_join(outdir, file, fpth)
      inquire(file=trim(adjustl(fpth)), exist=exst)
      if(.not.exst)then
        message(1) = "Could not find the file: '"//trim(adjustl(file))//"',"
        message(2) = "at directory: '"//trim(adjustl(reldir))//"."
        call messages_fatal(2)
      end if
    end if

    POP_SUB(fio_config_parse_get_file)
  end subroutine fio_config_parse_get_file

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
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))

    POP_SUB(fio_config_parse_simul_box)
  end subroutine fio_config_parse_simul_box

  ! ---------------------------------------------------------
  subroutine fio_config_parse_mesh(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname

    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: ierr

    PUSH_SUB(fio_config_parse_mesh)

    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))

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
  subroutine fio_config_parse_density(this, dirname, usedensity)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity

    type(json_array_t), pointer :: list
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: indx, nspin, ierr

    PUSH_SUB(fio_config_parse_density)

    nullify(list)
    call json_get(this, "nspin", nspin, ierr)
    ASSERT(ierr==JSON_OK)
    if(present(usedensity)) call json_set(this, "allocate", usedensity)
    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "files", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)==nspin)
    do indx = 1, nspin
      call json_get(list, indx, file, ierr)
      ASSERT(ierr==JSON_OK)
      call fio_config_parse_get_file(dirname, idir, odir, file)
      call json_set(list, indx, trim(adjustl(file)), ierr)
      ASSERT(ierr==JSON_OK)
    end do
    call json_set(this, "dir", trim(adjustl(odir)))
    nullify(list)

    POP_SUB(fio_config_parse_density)
  end subroutine fio_config_parse_density

  ! ---------------------------------------------------------
  subroutine fio_config_parse_states(this, dirname, usedensity)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_states)

    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_density(cnfg, dirname, usedensity)
    nullify(cnfg)

    POP_SUB(fio_config_parse_states)
  end subroutine fio_config_parse_states

  ! ---------------------------------------------------------
  subroutine fio_config_parse_system(this, dirname, usedensity)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usedensity

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_system)

    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_states(cnfg, dirname, usedensity)
    nullify(cnfg)

    POP_SUB(fio_config_parse_system)
  end subroutine fio_config_parse_system

  ! ---------------------------------------------------------
  subroutine fio_config_parse_external(this, dirname, usepotential)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usepotential

    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: ierr

    PUSH_SUB(fio_config_parse_external)

    if(present(usepotential)) call json_set(this, "allocate", usepotential)
    call json_get(this, "dir", idir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(this, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))

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
  subroutine fio_config_parse_read(this, dirname)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname

    type(json_parser_t)         :: parser
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: iunit, ierr

    PUSH_SUB(fio_config_parse_read)

    call json_init(this)
    call fio_config_parse_get_file(dirname, input_static_dir, dir, "config.json")
    call path_join(dir, "config.json", file)
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
  subroutine fio_config_parse_model(this, dirname, usepotential, usedensity)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    logical,   optional, intent(in)    :: usepotential
    logical,   optional, intent(in)    :: usedensity

    type(json_object_t), pointer :: cnfg
    type(json_object_t)          :: modl
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_model)

    nullify(cnfg)
    call fio_config_parse_read(modl, dirname)
    call json_update(this, modl)
    call json_end(modl)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simulation(cnfg, dirname)
    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_system(cnfg, dirname, usedensity)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_hamiltonian(cnfg, dirname, usepotential)
    nullify(cnfg)

    POP_SUB(fio_config_parse_model)
  end subroutine fio_config_parse_model

  ! ---------------------------------------------------------
  subroutine fio_config_parse_dir(this, dirname, interpolation, usepotential, usedensity)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    integer,   optional, intent(in)  :: interpolation
    logical,   optional, intent(in)  :: usepotential
    logical,   optional, intent(in)  :: usedensity

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_config_parse_dir)

    nullify(cnfg)
    call base_config_parse(this)
    call json_set(this, "type", HNDL_TYPE_FNIO)
    call json_set(this, "name", "fio")
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_model(cnfg, dirname, usepotential, usedensity)
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
    logical                      :: usepot, userho

    PUSH_SUB(fio_config_parse_block)

    call parse_block_string(block, line, icol, dirname)
    intrp = NEAREST
    usepot = .true.
    userho = .true.
    if(ncol>(icol+1))call parse_block_integer(block, line, icol+1, intrp)
    if(ncol>(icol+2))call parse_block_logical(block, line, icol+2, usepot)
    if(ncol>(icol+3))call parse_block_logical(block, line, icol+3, userho)
    call fio_config_parse(this, dirname, intrp, usepot, userho)

    POP_SUB(fio_config_parse_block)
  end subroutine fio_config_parse_block

end module fio_config_m

!! Local Variables:
!! mode: f90
!! End:
