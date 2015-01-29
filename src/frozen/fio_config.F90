#include "global.h"

module fio_config_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m, only: CURV_METHOD_UNIFORM
  use intrpl_m,      only: NEAREST
  use io_m,          only: io_open, io_close
  use json_m,        only: JSON_OK, json_object_t, json_array_t
  use json_m,        only: json_isdef, json_len, json_init, json_set, json_get, json_append, json_copy, json_end
  use json_parser_m, only: json_parser_t, json_parser_init, json_parser_end
  use json_parser_m, only: json_parser_parse, json_parser_error
  use kinds_m,       only: wp
  use loct_m,        only: loct_dir_exists
  use parser_m,      only: block_t, parse_block_string, parse_block_integer

  use bhmlt_m, only: &
    POTN_TYPE,       &
    HMLT_TYPE

  use bcnfg_m, only: &
    FNIO_TYPE,       &
    bcnfg_parse

  implicit none

  private
  public ::           &
    fio_config_parse

  integer,          parameter :: default_ndim        = 3
  integer,          parameter :: default_nspin       = 1
  character(len=*), parameter :: input_mesh_dir      = "./restart/"//GS_DIR
  character(len=*), parameter :: input_mesh_file     = "mesh"
  character(len=*), parameter :: input_static_dir    = "./"//STATIC_DIR
  character(len=*), parameter :: input_density_file  = "mesh"
  character(len=*), parameter :: input_external_file = "v0.obf"

contains

  ! ---------------------------------------------------------
  subroutine fio_config_parse_get_dir(topdir, reldir, outdir)
    use path_m
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
    use path_m
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
  subroutine fio_config_parse_simul_box(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: idim, ierr
    !
    call json_get(that, "dimensions", idim, ierr)
    if(ierr/=JSON_OK)idim=default_ndim
    call json_set(this, "dimensions", idim)
    call json_get(that, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_mesh_dir
    call json_get(that, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_mesh_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_simul_box

  ! ---------------------------------------------------------
  subroutine fio_config_parse_curvilinear(this, that)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    !
    integer :: mthd, ierr
    !
    call json_get(that, "method", mthd, ierr)
    if(ierr/=JSON_OK)mthd=CURV_METHOD_UNIFORM
    call json_set(this, "method", mthd)
    return
  end subroutine fio_config_parse_curvilinear

  ! ---------------------------------------------------------
  subroutine fio_config_parse_mesh(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_array_t), pointer :: olst, ilst
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: ierr
    !
    nullify(olst, ilst)
    call json_get(this, "spacing", olst, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "spacing", ilst, ierr)
    ASSERT(ierr==JSON_OK)
    call json_copy(olst, ilst)
    nullify(olst, ilst)
    call json_get(that, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_mesh_dir
    call json_get(that, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_mesh_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_mesh

  ! ---------------------------------------------------------
  subroutine fio_config_parse_grid(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: ierr
    !
    nullify(ocfg, icfg)
    call json_get(this, "simul_box", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "simul_box", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simul_box(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    call json_get(this, "curvilinear", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "curvilinear", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_curvilinear(ocfg, icfg)
    nullify(ocfg, icfg)
    call json_get(this, "mesh", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "mesh", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_mesh(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    return
  end subroutine fio_config_parse_grid

  ! ---------------------------------------------------------
  subroutine fio_config_parse_simulation(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: ierr
    !
    nullify(ocfg, icfg)
    call json_get(this, "grid", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "grid", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_grid(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    return
  end subroutine fio_config_parse_simulation

  ! ---------------------------------------------------------
  subroutine fio_config_parse_space(this, that)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    !
    integer :: idim, ierr
    !
    call json_get(that, "dimensions", idim, ierr)
    if(ierr/=JSON_OK)idim=default_ndim
    call json_set(this, "dimensions", idim)
    return
  end subroutine fio_config_parse_space

  ! ---------------------------------------------------------
  subroutine fio_config_parse_geometry(this, that)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    !
    type(json_array_t), pointer :: olst, ilst
    integer                     :: n, ierr
    !
    nullify(olst, ilst)
    call json_get(that, "nspecies", n, ierr)
    if((ierr==JSON_OK).and.(n>0))then
      call json_set(this, "nspecies", n)
      call json_get(this, "species", olst, ierr)
      if(ierr/=JSON_OK)then
        SAFE_ALLOCATE(olst)
        call json_init(olst)
        call json_set(this, "species", olst)
      end if
      call json_get(that, "species", ilst, ierr)
      ASSERT(ierr==JSON_OK)
      call json_copy(olst, ilst)
      nullify(olst, ilst)
    end if
    call json_get(that, "natoms", n, ierr)
    if((ierr==JSON_OK).and.(n>0))then
      call json_set(this, "natoms", n)
      call json_get(this, "atom", olst, ierr)
      if(ierr/=JSON_OK)then
        SAFE_ALLOCATE(olst)
        call json_init(olst)
        call json_set(this, "atom", olst)
      end if
      call json_get(that, "atom", ilst, ierr)
      ASSERT(ierr==JSON_OK)
      call json_copy(olst, ilst)
      nullify(olst, ilst)
    end if
    call json_get(that, "ncatoms", n, ierr)
    if((ierr==JSON_OK).and.(n>0))then
      call json_set(this, "ncatoms", n)
      call json_get(this, "catom", olst, ierr)
      if(ierr/=JSON_OK)then
        SAFE_ALLOCATE(olst)
        call json_init(olst)
        call json_set(this, "catom", olst)
      end if
      call json_get(that, "catom", ilst, ierr)
      ASSERT(ierr==JSON_OK)
      call json_copy(olst, ilst)
      nullify(olst, ilst)
    end if
    return
  end subroutine fio_config_parse_geometry

  ! ---------------------------------------------------------
  subroutine fio_config_parse_density(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_array_t), pointer :: ilst, olst
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: indx, nspin, ierr
    !
    nullify(ilst, olst)
    call json_get(that, "nspin", nspin, ierr)
    if(ierr/=JSON_OK)nspin=default_nspin
    ASSERT(nspin>0)
    call json_set(this, "nspin", nspin)
    call json_get(that, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_static_dir
    call json_get(that, "files", ilst, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(ilst)
      call json_init(ilst)
      do indx=1, nspin
        if(nspin>1)then
          write(unit=file, fmt="(a,i1,a)") "density-sp", indx, ".obf"
        else
          file="density.obf"
        end if
        call json_append(ilst, trim(adjustl(file)))
      end do
    end if
    call json_get(this, "files", olst, ierr)
    ASSERT(ierr/=JSON_OK)
    SAFE_ALLOCATE(olst)
    call json_init(olst)
    call json_set(this, "files", olst)
    do indx=1, nspin
      call json_get(ilst, indx, file, ierr)
      ASSERT(ierr==JSON_OK)
      call fio_config_parse_get_file(dirname, idir, odir, file)
      call json_append(olst, trim(adjustl(file)))
    end do
    call json_set(this, "dir", trim(adjustl(odir)))
    nullify(ilst, olst)
    return
  end subroutine fio_config_parse_density

  ! ---------------------------------------------------------
  subroutine fio_config_parse_states(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: ierr
    !
    nullify(ocfg, icfg)
    call json_get(this, "density", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "density", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_density(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    return
  end subroutine fio_config_parse_states

  ! ---------------------------------------------------------
  subroutine fio_config_parse_system(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: ierr
    !
    nullify(ocfg, icfg)
    call json_get(this, "space", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "space", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_space(ocfg, icfg)
    nullify(ocfg, icfg)
    call json_get(this, "geometry", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "geometry", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_geometry(ocfg, icfg)
    nullify(ocfg, icfg)
    call json_get(this, "states", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "states", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_states(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    return
  end subroutine fio_config_parse_system

  ! ---------------------------------------------------------
  subroutine fio_config_parse_external(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    character(len=MAX_PATH_LEN) :: idir, odir, file
    integer                     :: type, ierr
    !
    call json_get(that, "type", type, ierr)
    if(ierr/=JSON_OK)type=POTN_TYPE
    call json_set(this, "type", type)
    call json_get(that, "dir", idir, ierr)
    if(ierr/=JSON_OK)idir=input_static_dir
    call json_get(that, "file", file, ierr)
    if(ierr/=JSON_OK)file=input_external_file
    call fio_config_parse_get_file(dirname, idir, odir, file)
    call json_set(this, "dir", trim(adjustl(odir)))
    call json_set(this, "file", trim(adjustl(file)))
    return
  end subroutine fio_config_parse_external

  ! ---------------------------------------------------------
  subroutine fio_config_parse_hamiltonian(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: type, ierr
    !
    nullify(ocfg, icfg)
    call json_get(that, "type", type, ierr)
    if(ierr/=JSON_OK)type=HMLT_TYPE
    call json_set(this, "type", type)
    call json_get(this, "external", ocfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(ocfg)
      call json_init(ocfg)
      call json_set(this, "external", ocfg)
    end if
    call json_get(that, "external", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_external(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    return
  end subroutine fio_config_parse_hamiltonian

  ! ---------------------------------------------------------
  subroutine fio_config_parse_model(this, that, dirname)
    type(json_object_t), intent(inout) :: this
    type(json_object_t), intent(in)    :: that
    character(len=*),    intent(in)    :: dirname
    !
    type(json_object_t), pointer :: ocfg, icfg
    integer                      :: ierr
    !
    PUSH_SUB(fio_config_parse_model)
    nullify(ocfg, icfg)
    call json_get(this, "simulation", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "simulation", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_simulation(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    call json_get(this, "system", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "system", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_system(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    call json_get(this, "hamiltonian", ocfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(that, "hamiltonian", icfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_hamiltonian(ocfg, icfg, dirname)
    nullify(ocfg, icfg)
    POP_SUB(fio_config_parse_model)
    return
  end subroutine fio_config_parse_model

  ! ---------------------------------------------------------
  subroutine fio_config_parse_read(this, dirname)
    type(json_object_t), intent(out) :: this
    character(len=*),    intent(in)  :: dirname
    !
    type(json_parser_t)         :: parser
    character(len=MAX_PATH_LEN) :: dir
    integer                     :: iunit, ierr
    !
    PUSH_SUB(fio_config_parse_read)
    call json_init(this)
    call fio_config_parse_get_file(dirname, input_static_dir, dir, "config.json")
    iunit=io_open(trim(adjustl(dir))//"/"//"config.json", action='read', status="old", is_tmp=.true.)
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
  subroutine fio_config_parse(this, block, line, icol, ncol)
    type(json_object_t), intent(out) :: this
    type(block_t),       intent(in)  :: block
    integer,             intent(in)  :: line
    integer,             intent(in)  :: icol
    integer,             intent(in)  :: ncol
    !
    type(json_object_t)          :: that
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dirname
    integer                      :: intrp, ierr
    !
    PUSH_SUB(fio_config_parse)
    nullify(cnfg)
    call bcnfg_parse(this)
    call json_set(this, "name", "fio")
    call json_set(this, "type", FNIO_TYPE)
    call parse_block_string(block, line, icol, dirname)
    call fio_config_parse_read(that, dirname)
    call json_get(this, "model", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_config_parse_model(cnfg, that, dirname)
    nullify(cnfg)
    if(ncol>(icol+1))then
      call parse_block_integer(block, line, icol+1, intrp)
      call json_set(this, "interpolation", intrp)
    end if
    call json_end(that)
    POP_SUB(fio_config_parse)
    return
  end subroutine fio_config_parse

end module fio_config_m

!! Local Variables:
!! mode: f90
!! End:
