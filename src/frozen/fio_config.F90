#include "global.h"

module fio_config_m

  use global_m
  use messages_m
  use profiling_m
  use curvilinear_m
  use intrpl_m
  use io_m
  use json_m
  use json_m
  use json_m
  use json_parser_m
  use json_parser_m
  use kinds_m
  use loct_m
  use parser_m
  use path_m
  use species_m
  use base_hamiltonian_m
  use base_config_m
  use fio_handle_m

  implicit none

  private
  public ::           &
    fio_config_parse

  interface fio_config_parse
    module procedure fio_config_parse_dir
    module procedure fio_config_parse_block
  end interface fio_config_parse

  integer,          public, parameter :: default_ndim        = 3
  integer,          public, parameter :: default_nspin       = 1
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
    !
    character(len=max_path_len) :: apth, bpth
    logical                     :: exst
    !
    PUSH_SUB(fio_config_parse_get_dir)
    call path_realpath(reldir, outdir)
    if(.not.path_isabs(reldir))then
       call path_realpath(topdir, apth)
       call path_join(apth, reldir, bpth)
       call path_realpath(bpth, outdir)
       exst=loct_dir_exists(outdir)
       if(.not.exst)then
          call path_getcwd(apth)
          call path_join(apth, reldir, bpth)
          call path_realpath(bpth, outdir)
          exst=loct_dir_exists(outdir)
          if(.not.exst)call path_realpath(topdir, outdir)
       end if
    end if
    POP_SUB(fio_config_parse_get_dir)
    return
  end subroutine fio_config_parse_get_dir

  ! ---------------------------------------------------------
  subroutine fio_config_parse_get_file(topdir, reldir, outdir, file)
    character(len=*), intent(in)  :: topdir
    character(len=*), intent(in)  :: reldir
    character(len=*), intent(out) :: outdir
    character(len=*), intent(in)  :: file
    !
    character(len=max_path_len) :: fpth
    logical                     :: exst
    !
    PUSH_SUB(fio_config_parse_get_file)
    call fio_config_parse_get_dir(topdir, reldir, outdir)
    call path_join(outdir, file, fpth)
    inquire(file=trim(adjustl(fpth)), exist=exst)
    if(.not.exst)then
      call path_realpath(topdir, outdir)
      call path_join(outdir, file, fpth)
      inquire(file=trim(adjustl(fpth)), exist=exst)
      if(.not.exst)then
        message(1)="Could not find the file: '"//trim(adjustl(file))//"',"
        message(2)="at directory: '"//trim(adjustl(reldir))//"."
        call messages_fatal(2)
      end if
    end if
    POP_SUB(fio_config_parse_get_file)
    return
  end subroutine fio_config_parse_get_file

  ! ---------------------------------------------------------
  subroutine fio_config_parse_simul_box(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: idim, ierr
    !
    call json_get(this, "dimensions", idim, ierr)
    if(ierr/=JSON_OK)call json_set(this, "dimensions", default_ndim)
    call json_get(this, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_mesh_dir
    call json_get(this, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_mesh_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_simul_box

  ! ---------------------------------------------------------
  subroutine fio_config_parse_curvilinear(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: mthd, ierr
    !
    call json_get(this, "method", mthd, ierr)
    if(ierr/=JSON_OK)call json_set(this, "method", CURV_METHOD_UNIFORM)
    call json_get(this, "method", mthd, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(mthd==CURV_METHOD_UNIFORM)
    return
  end subroutine fio_config_parse_curvilinear

  ! ---------------------------------------------------------
  subroutine fio_config_parse_mesh(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_array_t), pointer :: list
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: ierr
    !
    nullify(list)
    call json_get(this, "spacing", list, ierr)
    ASSERT(ierr==JSON_OK)
    nullify(list)
    call json_get(this, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_mesh_dir
    call json_get(this, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_mesh_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_mesh

  ! ---------------------------------------------------------
  subroutine fio_config_parse_grid(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(this, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simul_box(cnfg, dirname)
    nullify(cnfg)
    call json_get(this, "curvilinear", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_curvilinear(cnfg)
    nullify(cnfg)
    call json_get(this, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_mesh(cnfg, dirname)
    nullify(cnfg)
    return
  end subroutine fio_config_parse_grid

  ! ---------------------------------------------------------
  subroutine fio_config_parse_simulation(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(this, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_grid(cnfg, dirname)
    nullify(cnfg)
    return
  end subroutine fio_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine fio_config_parse_space(this)
    type(json_object_t), intent(inout) :: this
    !
    integer :: idim, ierr
    !
    call json_get(this, "dimensions", idim, ierr)
    if(ierr/=JSON_OK)call json_set(this, "dimensions", default_ndim)
    return
  end subroutine fio_config_parse_space

  ! ---------------------------------------------------------
  subroutine fio_config_parse_species(this)
    type(json_object_t), intent(inout) :: this
    !
    call json_set(this, "type", SPECIES_FROZEN)
    return
  end subroutine fio_config_parse_species

  ! ---------------------------------------------------------
  subroutine fio_config_parse_geometry(this)
    type(json_object_t), intent(inout) :: this
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    integer                      :: n, ierr
    !
    nullify(cnfg, list)
    call json_get(this, "nspecies", n, ierr)
    ASSERT((ierr==JSON_OK).and.(n>0))
    call json_get(this, "species", list, ierr)
    ASSERT((ierr==JSON_OK).and.(json_len(list)==n))
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call fio_config_parse_species(cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call json_get(this, "natoms", n, ierr)
    ASSERT((ierr==JSON_OK).and.(n>0))
    call json_get(this, "atom", list, ierr)
    ASSERT((ierr==JSON_OK).and.(json_len(list)==n))
    nullify(list)
    call json_get(this, "ncatoms", n, ierr)
    if((ierr==JSON_OK).and.(n>0))then
      call json_get(this, "catom", list, ierr)
      ASSERT((ierr==JSON_OK).and.(json_len(list)==n))
      nullify(list)
    end if
    return
  end subroutine fio_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine fio_config_parse_density(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_array_t), pointer :: list
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: indx, nspin, ierr
    !
    nullify(list)
    call json_get(this, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=default_nspin
    ASSERT(nspin>0)
    call json_set(this, "nspin", nspin)
    call json_get(this, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_static_dir
    call json_get(this, "files", list, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(list)
      call json_init(list)
      call json_set(this, "files", list)
      if(nspin>1)then
        do indx=1, nspin
          write(unit=file, fmt="(a,i1,a)") "density-sp", indx, ".obf"
          call json_append(list, trim(adjustl(file)))
        end do
      else
        call json_append(list, "density.obf")
      end if
    end if
    ASSERT(json_len(list)==nspin)
    do indx=1, nspin
      call json_get(list, indx, file, ierr)
      ASSERT(ierr==JSON_OK)
      call fio_config_parse_get_file(dirname, idir, odir, file)
      call json_set(list, indx, trim(adjustl(file)), ierr)
      ASSERT(ierr==JSON_OK)
    end do
    call json_set(this, "dir", trim(adjustl(odir)))
    nullify(list)
    return
  end subroutine fio_config_parse_density

  ! ---------------------------------------------------------
  subroutine fio_config_parse_states(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    real(kind=wp)                :: charge
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(this, "charge", charge, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(charge>0.0_wp)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_density(cnfg, dirname)
    nullify(cnfg)
    return
  end subroutine fio_config_parse_states

  ! ---------------------------------------------------------
  subroutine fio_config_parse_system(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    nullify(cnfg)
    call json_get(this, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_space(cnfg)
    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_geometry(cnfg)
    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_states(cnfg, dirname)
    nullify(cnfg)
    return
  end subroutine fio_config_parse_system

  ! ---------------------------------------------------------
  subroutine fio_config_parse_external(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: type, ierr
    !
    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)
    call json_get(this, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_static_dir
    call json_get(this, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_external_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_external

  ! ---------------------------------------------------------
  subroutine fio_config_parse_hamiltonian(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr
    !
    nullify(cnfg)
    call json_get(this, "type", type, ierr)
    if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_HMLT)
    call json_get(this, "external", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "external", cnfg)
    end if
    call fio_config_parse_external(cnfg, dirname)
    nullify(cnfg)
    return
  end subroutine fio_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine fio_config_parse_model(this, dirname)
    type(json_object_t), intent(inout) :: this
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(fio_config_parse_model)
    nullify(cnfg)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simulation(cnfg, dirname)
    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_system(cnfg, dirname)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_hamiltonian(cnfg, dirname)
    nullify(cnfg)
    POP_SUB(fio_config_parse_model)
    return
  end subroutine fio_config_parse_model

  ! ---------------------------------------------------------
  subroutine fio_config_parse_read(this, dirname)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    !
    type(json_parser_t)         :: parser
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: iunit, ierr
    !
    PUSH_SUB(fio_config_parse_read)
    call json_init(this)
    call fio_config_parse_get_file(dirname, input_static_dir, dir, "config.json")
    call path_join(dir, "config.json", file)
    iunit=io_open(file, action='read', status="old")
    ASSERT(iunit>0)
    call json_parser_init(parser, iunit, ierr)
    if(ierr/=0)call json_parser_error(parser, fatal=.true.)
    call json_parser_parse(parser, this, ierr)
    if((ierr/=0).or.(.not.json_isdef(this)))&
      call json_parser_error(parser, fatal=.true.)
    call json_parser_end(parser)
    call io_close(iunit)
    POP_SUB(fio_config_parse_read)
    return
  end subroutine fio_config_parse_read

  ! ---------------------------------------------------------
  subroutine fio_config_parse_dir(this, dirname)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    !
    type(json_object_t)          :: that
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(fio_config_parse_dir)
    nullify(cnfg)
    call base_config_parse(this)
    call json_set(this, "name", "fio")
    call json_set(this, "type", HNDL_TYPE_FNIO)
    call fio_config_parse_read(that, dirname)
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_update(cnfg, that)
    call json_end(that)
    call fio_config_parse_model(cnfg, dirname)
    nullify(cnfg)
    POP_SUB(fio_config_parse_dir)
    return
  end subroutine fio_config_parse_dir

  ! ---------------------------------------------------------
  subroutine fio_config_parse_block(this, block, line, icol, ncol)
    type(json_object_t), intent(out) :: this
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: line
    integer,             intent(in)  :: icol
    integer,             intent(in)  :: ncol
    !
    character(len=MAX_PATH_LEN)  :: dirname
    integer                      :: intrp
    !
    PUSH_SUB(fio_config_parse_block)
    call parse_block_string(block, line, icol, dirname)
    call fio_config_parse(this, dirname)
    if(ncol>(icol+1))then
      call parse_block_integer(block, line, icol+1, intrp)
      call json_set(this, "interpolation", intrp)
    end if
    POP_SUB(fio_config_parse_block)
    return
  end subroutine fio_config_parse_block

end module fio_config_m

!! Local Variables:
!! mode: f90
!! End:
