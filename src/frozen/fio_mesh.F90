#include "global.h"

module fio_mesh_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m,        only: CURV_METHOD_UNIFORM, curvilinear_t
  use checksum_interface_m, only: checksum_calculate
  use hypercube_m,          only: hypercube_end
  use index_m,              only: index_t, index_load_lxyz
  use json_m,               only: JSON_OK, json_object_t, json_get
  use kinds_m,              only: wp
  use mesh_cube_map_m,      only: mesh_cube_map_init
  !use mesh_init_m,          only: mesh_get_vol_pp
  use mpi_m,                only: mpi_grp_t

  use mesh_m, only:        &
    mesh_load,             &
    mesh_read_fingerprint, &
    mesh_x_global

  use imesh_m, only:            &
    fio_mesh_t    => mesh_t,    &
    fio_mesh_copy => mesh_copy, &
    fio_mesh_end  => mesh_end

  use fio_simul_box_m, only: &
    fio_simul_box_t

  implicit none

  private
  public ::         &
    fio_mesh_t,     &
    fio_mesh_init,  &
    fio_mesh_start, &
    fio_mesh_stop,  &
    fio_mesh_copy,  &
    fio_mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_index_init(this, np, ndim, mpi_grp, config)
    type(index_t),         intent(out) :: this
    integer,               intent(in)  :: np
    integer,               intent(in)  :: ndim
    type(mpi_grp_t),       intent(in)  :: mpi_grp
    type(json_object_t),   intent(in)  :: config
    !
    character(len=MAX_PATH_LEN) :: dir
    integer                     :: indx, ierr
    integer                     :: i11, i21, i12, i22, i13, i23
    !
    PUSH_SUB(fio_index_init)
    SAFE_ALLOCATE(this%lxyz(np,MAX_DIM))
    i11=this%nr(1,1)
    i21=this%nr(2,1)
    i12=this%nr(1,2)
    i22=this%nr(2,2)
    i13=this%nr(1,3)
    i23=this%nr(2,3)
    SAFE_ALLOCATE(this%lxyz_inv(i11:i21,i12:i22,i13:i23))
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call index_load_lxyz(this, np, dir, mpi_grp, ierr)
    if(ierr/=0) then
      message(1) = "Error: failed to read file: '"//trim(adjustl(dir))//"/lxyz.obf'."
      call messages_fatal(1)
    end if
    call checksum_calculate(1, np*ndim, this%lxyz(1,1), this%checksum)
    POP_SUB(fio_index_init)
    return
  end subroutine fio_index_init

  ! ---------------------------------------------------------
  subroutine fio_index_end(this)
    type(index_t), intent(inout) :: this
    !
    character(len=MAX_PATH_LEN) :: dir
    integer                     :: indx, ierr
    integer                     :: i11, i21, i12, i22, i13, i23
    !
    PUSH_SUB(fio_index_end)
    if(this%is_hypercube)call hypercube_end(this%hypercube)
    this%is_hypercube=.false.
    this%dim=0
    this%nr=0
    this%ll=0
    SAFE_DEALLOCATE_P(this%lxyz)
    nullify(this%lxyz)
    SAFE_DEALLOCATE_P(this%lxyz_inv)
    nullify(this%lxyz_inv)
    this%enlarge=0
    this%checksum=0
    POP_SUB(fio_index_end)
    return
  end subroutine fio_index_end

  ! ---------------------------------------------------------
  subroutine fio_mesh_init(this, sb, cv, config)
    type(fio_mesh_t),              intent(out) :: this
    type(fio_simul_box_t), target, intent(in)  :: sb
    type(curvilinear_t),   target, intent(in)  :: cv
    type(json_object_t),           intent(in)  :: config
    !
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: i, ia, ib, ierr
    !
    PUSH_SUB(fio_mesh_init)
    ASSERT(.not.sb%mr_flag)
    this%sb=>sb
    this%cv=>cv
    call json_get(config, "spacing", this%spacing(1:sb%dim), ierr=ierr)
    ASSERT(ierr==JSON_OK)
    this%spacing=(/this%spacing(1:sb%dim),(0.0_wp,i=sb%dim+1,MAX_DIM)/)
    POP_SUB(fio_mesh_init)
    return
  end subroutine fio_mesh_init

  ! ---------------------------------------------------------
  subroutine fio_mesh_start(this, mpi_grp, config)
    type(fio_mesh_t),    intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: mpi_grp
    type(json_object_t), intent(in)    :: config
    !
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: i, ia, ib, ierr
    !
    PUSH_SUB(fio_mesh_start)
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call mesh_load(this, dir, file, mpi_grp, ierr)
    if(ierr==0)then
      this%np=this%np_global
      this%np_part=this%np_part_global
      call fio_index_init(this%idx, this%np_part_global, this%sb%dim, mpi_grp, config)
      call mesh_read_fingerprint(this, dir, "grid", mpi_grp, ia, ib, ierr)
      ASSERT((ia==0).and.(ib==0).and.(ierr==0))
      this%use_curvilinear=(this%cv%method/=CURV_METHOD_UNIFORM).or.this%sb%mr_flag
      nullify(this%resolution)
      SAFE_ALLOCATE(this%x(this%np_part_global,MAX_DIM))
      do i=1, this%np_part_global
        this%x(i,:) = mesh_x_global(this, i, .true.)
      end do
      call mesh_cube_map_init(this%cube_map, this%idx, this%np_global)
      !call mesh_get_vol_pp(this, sb)
    else
      message(1)="Could not open the mesh info file: '"//trim(adjustl(file))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(2)
    end if
    POP_SUB(fio_mesh_start)
    return
  end subroutine fio_mesh_start

  ! ---------------------------------------------------------
  subroutine fio_mesh_stop(this)
    type(fio_mesh_t), intent(inout) :: this
    !
    PUSH_SUB(fio_mesh_stop)
    call fio_index_end(this%idx)
    SAFE_DEALLOCATE_P(this%x)
    SAFE_DEALLOCATE_P(this%resolution)
    SAFE_DEALLOCATE_P(this%vol_pp)
    nullify(this%x, this%resolution, this%vol_pp)
    POP_SUB(fio_mesh_stop)
    return
  end subroutine fio_mesh_stop

end module fio_mesh_m

!! Local Variables:
!! mode: f90
!! End:
