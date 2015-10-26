#include "global.h"

module fio_mesh_m

  use curvilinear_m
  use checksum_interface_m
  use global_m
  use hypercube_m
  use index_m
  use json_m
  use kinds_m
  use mesh_m
  use mesh_cube_map_m
  use messages_m
  use mpi_m
  use path_m
  use profiling_m
  use simul_box_m

  implicit none

  private

  public ::         &
    fio_mesh_init,  &
    fio_mesh_copy,  &
    fio_mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_index_init(this, np, ndim, mpi_grp, config)
    type(index_t),         intent(inout) :: this
    integer,               intent(in)    :: np
    integer,               intent(in)    :: ndim
    type(mpi_grp_t),       intent(in)    :: mpi_grp
    type(json_object_t),   intent(in)    :: config

    character(len=MAX_PATH_LEN)       :: dir, fpth
    integer(kind=kind(this%checksum)) :: chksm
    integer                           :: ierr
    integer                           :: i11, i21, i12, i22, i13, i23

    PUSH_SUB(fio_index_init)

    SAFE_ALLOCATE(this%lxyz(np,MAX_DIM))
    i11 = this%nr(1,1)
    i21 = this%nr(2,1)
    i12 = this%nr(1,2)
    i22 = this%nr(2,2)
    i13 = this%nr(1,3)
    i23 = this%nr(2,3)
    SAFE_ALLOCATE(this%lxyz_inv(i11:i21,i12:i22,i13:i23))
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call index_load_lxyz(this, np, dir, mpi_grp, ierr)
    if(ierr/=0) then
      call path_join(dir, "lxyz.obf", fpth)
      message(1) = "Error: failed to read file: '"//trim(adjustl(fpth))//"'."
      call messages_fatal(1)
    end if
    call checksum_calculate(1, np*ndim, this%lxyz(1,1), chksm)
    ASSERT(chksm==this%checksum)

    POP_SUB(fio_index_init)
  end subroutine fio_index_init

  ! ---------------------------------------------------------
  subroutine fio_index_end(this)
    type(index_t), intent(inout) :: this

    PUSH_SUB(fio_index_end)

    if(this%is_hypercube) call hypercube_end(this%hypercube)
    this%is_hypercube = .false.
    this%dim = 0
    this%nr = 0
    this%ll = 0
    SAFE_DEALLOCATE_A(this%lxyz)
    SAFE_DEALLOCATE_A(this%lxyz_inv)
    this%enlarge = 0
    this%checksum = int(0, kind=kind(this%checksum))

    POP_SUB(fio_index_end)
  end subroutine fio_index_end

  ! ---------------------------------------------------------
  subroutine fio_mesh_init(this, sb, cv, mpi_grp, config)
    type(mesh_t),                intent(out) :: this
    type(simul_box_t),   target, intent(in)  :: sb
    type(curvilinear_t), target, intent(in)  :: cv
    type(mpi_grp_t),             intent(in)  :: mpi_grp
    type(json_object_t),         intent(in)  :: config

    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: i, ia, ib, ierr

    PUSH_SUB(fio_mesh_init)

    ASSERT(.not.sb%mr_flag)
    ASSERT(cv%method==CURV_METHOD_UNIFORM)
    this%sb => sb
    this%cv => cv
    this%use_curvilinear = .false.
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call mesh_load(this, dir, file, mpi_grp, ierr)
    if(ierr==0)then
      ASSERT(.not.this%idx%is_hypercube)
      call fio_index_init(this%idx, this%np_part_global, this%sb%dim, mpi_grp, config)
      this%np = this%np_global
      this%np_part = this%np_part_global
      this%spacing = 0.0_wp
      call json_get(config, "spacing", this%spacing(1:this%sb%dim), ierr=ierr)
      ASSERT(ierr==JSON_OK)
      SAFE_ALLOCATE(this%vol_pp(1:1))
      this%vol_pp(1) = product(this%spacing(1:this%sb%dim))
      this%volume_element = this%vol_pp(1)
      SAFE_ALLOCATE(this%x(1:this%np_part_global,1:MAX_DIM))
      do i=1, this%np_part_global
        this%x(i,:) = mesh_x_global(this, i, .true.)
      end do
      call mesh_cube_map_init(this%cube_map, this%idx, this%np_global)
      call mesh_read_fingerprint(this, dir, "grid", mpi_grp, ia, ib, ierr)
      ASSERT((ia==0).and.(ib==0).and.(ierr==0))
    else
      message(1)="Could not open the mesh info file: '"//trim(adjustl(file))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if

    POP_SUB(fio_mesh_init)
  end subroutine fio_mesh_init

  ! ---------------------------------------------------------
  subroutine fio_mesh_copy(this, that)
    type(mesh_t), intent(inout) :: this
    type(mesh_t), intent(in)    :: that

    PUSH_SUB(fio_mesh_copy)

    ASSERT(.false.)
    this=that

    POP_SUB(fio_mesh_copy)
  end subroutine fio_mesh_copy

  ! ---------------------------------------------------------
  subroutine fio_mesh_end(this)
    type(mesh_t), intent(inout) :: this

    PUSH_SUB(fio_mesh_end)

    call mesh_end(this)

    POP_SUB(fio_mesh_end)
  end subroutine fio_mesh_end

end module fio_mesh_m

!! Local Variables:
!! mode: f90
!! End:
