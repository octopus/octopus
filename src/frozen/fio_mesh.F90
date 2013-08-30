#include "global.h"

module fio_mesh_m

  use global_m
  use messages_m
  use profiling_m

  use datasets_m,      only: tmpdir
  use curvilinear_m,   only: CURV_METHOD_UNIFORM, curvilinear_t
  use checksum_m,      only: checksum_calculate
  use io_binary_m,     only: io_binary_read
  use io_m,            only: io_open, io_close
  use json_m,          only: JSON_OK, json_object_t, json_get
  use kinds_m,         only: wp
  use mesh_cube_map_m, only: mesh_cube_map_init
  !use mesh_init_m,     only: mesh_get_vol_pp

  use mesh_m, only:           &
    mesh_init_from_file,      &
    mesh_lxyz_init_from_file, &
    mesh_read_fingerprint,    &
    mesh_x_global

  use mesh_m, only:           &
    fio_mesh_t   =>mesh_t,    &
    fio_mesh_end =>mesh_end

  use fio_simul_box_m, only: &
    fio_simul_box_t

  implicit none

  private
  public ::        &
    fio_mesh_t,    &
    fio_mesh_init, &
    fio_mesh_copy, &
    fio_mesh_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_mesh_init(this, sb, cv, config)
    type(fio_mesh_t),              intent(out) :: this
    type(fio_simul_box_t), target, intent(in)  :: sb
    type(curvilinear_t),   target, intent(in)  :: cv
    type(json_object_t),           intent(in)  :: config
    !
    character(len=MAX_PATH_LEN) :: dir
    integer                     :: i, ia, ib, iunit, ierr
    integer                     :: i11, i21, i12, i22, i13, i23
    !
    PUSH_SUB(fio_mesh_init)
    ASSERT(.not.sb%mr_flag)
    call json_get(config, "dir", dir, ierr)
    if(ierr/=JSON_OK)dir="./"//trim(tmpdir)//GS_DIR
    iunit=io_open(trim(dir)//"mesh", action="read", status="old")
    if(iunit>0)then
      this%sb=>sb
      this%cv=>cv
      call json_get(config, "spacing", this%spacing(1:sb%dim), ierr=ierr)
      ASSERT(ierr==JSON_OK)
      this%spacing=(/this%spacing(1:sb%dim),(0.0_wp,i=sb%dim+1,MAX_DIM)/)
      call mesh_init_from_file(this, iunit)
      call io_close(iunit)
      this%np=this%np_global
      this%np_part=this%np_part_global
      this%idx%sb=>sb
      SAFE_ALLOCATE(this%idx%lxyz(this%np_part_global,MAX_DIM))
      i11=this%idx%nr(1,1)
      i21=this%idx%nr(2,1)
      i12=this%idx%nr(1,2)
      i22=this%idx%nr(2,2)
      i13=this%idx%nr(1,3)
      i23=this%idx%nr(2,3)
      SAFE_ALLOCATE(this%idx%lxyz_inv(i11:i21,i12:i22,i13:i23))
      call io_binary_read(trim(dir)//'lxyz.obf', this%np_part_global*this%sb%dim, this%idx%lxyz, ierr)
      if(ierr>0) then
        message(1) = "Error: failed to read file: '"//trim(dir)//"lxyz.obf'"
        call messages_fatal(1)
      end if
      this%idx%lxyz(:,this%sb%dim+1:MAX_DIM)=0
      forall(i = 1:this%np_part_global)
        this%idx%lxyz_inv(this%idx%lxyz(i,1),this%idx%lxyz(i,2),this%idx%lxyz(i,3))=i
      end forall
      call checksum_calculate(1, this%np_part_global*this%sb%dim, this%idx%lxyz(1,1), this%idx%checksum)
      call mesh_read_fingerprint(this, trim(dir)//"grid", ia, ib)
      ASSERT((ia==0).and.(ib==0))
      this%use_curvilinear=(cv%method/=CURV_METHOD_UNIFORM).or.sb%mr_flag
      nullify(this%resolution)
      SAFE_ALLOCATE(this%x(this%np_part_global,MAX_DIM))
      do i=1, this%np_part_global
        this%x(i,:) = mesh_x_global(this, i, .true.)
      end do
      call mesh_cube_map_init(this%cube_map, this%idx, this%np_global)
      !call mesh_get_vol_pp(this, sb)
    else
      message(1)="Could not open the mesh info file: '"//trim(dir)//"mesh'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(2)
    end if
    POP_SUB(fio_mesh_init)
    return
  end subroutine fio_mesh_init

  ! ---------------------------------------------------------
  subroutine fio_mesh_copy(this, that)
    type(fio_mesh_t), intent(out) :: this
    type(fio_mesh_t), intent(in)  :: that
    !
    PUSH_SUB(fio_mesh_copy)
    ASSERT(.false.)
    POP_SUB(fio_mesh_copy)
    return
  end subroutine fio_mesh_copy

end module fio_mesh_m

!! Local Variables:
!! mode: f90
!! End:
