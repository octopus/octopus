#include "global.h"

module fio_mesh_oct_m

  use curvilinear_oct_m
  use checksum_interface_oct_m
  use fio_index_oct_m
  use global_oct_m
  use hypercube_oct_m
  use index_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use mesh_cube_map_oct_m
  use messages_oct_m
  use mpi_oct_m
  use path_oct_m
  use profiling_oct_m
  use simul_box_oct_m

  implicit none

  private

  public ::         &
    fio_mesh_init,  &
    fio_mesh_copy,  &
    fio_mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_mesh_init(this, sb, cv, group, config)
    type(mesh_t),                intent(out) :: this
    type(simul_box_t),   target, intent(in)  :: sb
    type(curvilinear_t), target, intent(in)  :: cv
    type(mpi_grp_t),             intent(in)  :: group
    type(json_object_t),         intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: i, ia, ib, ierr

    PUSH_SUB(fio_mesh_init)

    nullify(cnfg)
    ASSERT(.not.sb%mr_flag)
    ASSERT(cv%method==CURV_METHOD_UNIFORM)
    this%sb => sb
    this%cv => cv
    this%use_curvilinear = .false.
    call json_get(config, "globalsize", this%np_global, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "globalfullsize", this%np_part_global, ierr)
    ASSERT(ierr==JSON_OK)
    this%np = this%np_global
    this%np_part = this%np_part_global
    this%parallel_in_domains = .false.
    call json_get(config, "spacing", this%spacing(1:this%sb%dim), ierr=ierr)
    if(size(this%spacing)>this%sb%dim) this%spacing(this%sb%dim+1:) = 0.0_wp
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(this%vol_pp(1:1))
    this%vol_pp(1) = product(this%spacing(1:this%sb%dim))
    this%volume_element = this%vol_pp(1)
    call json_get(config, "index", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(associated(cnfg))
    call fio_index_init(this%idx, group, cnfg)
    ASSERT(.not.this%idx%is_hypercube)
    ASSERT(this%idx%dim==this%sb%dim)
    nullify(cnfg)
    SAFE_ALLOCATE(this%x(1:this%np_part_global,1:MAX_DIM))
    do i=1, this%np_part_global
      this%x(i,:) = mesh_x_global(this, i, .true.)
    end do
    call mesh_cube_map_init(this%cube_map, this%idx, this%np_global)
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call mesh_read_fingerprint(this, dir, file, group, ia, ib, ierr)
    ASSERT((ia==0).and.(ib==0).and.(ierr==0))

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

end module fio_mesh_oct_m

!! Local Variables:
!! mode: f90
!! End:
