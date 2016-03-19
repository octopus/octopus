#include "global.h"

module fio_index_oct_m

  use checksum_interface_oct_m
  use global_oct_m
  use hypercube_oct_m
  use index_oct_m
  use json_oct_m
  use messages_oct_m
  use mpi_oct_m
  use profiling_oct_m

  implicit none

  private

  public ::         &
    fio_index_init, &
    fio_index_copy, &
    fio_index_end

contains
  
  ! ---------------------------------------------------------
  subroutine fio_index_init(this, group, config)
    type(index_t),         intent(out) :: this
    type(mpi_grp_t),       intent(in)  :: group
    type(json_object_t),   intent(in)  :: config

    character(len=MAX_PATH_LEN)       :: dir, file
    integer(kind=kind(this%checksum)) :: chksm
    integer                           :: size, ierr
    integer                           :: i11, i21, i12, i22, i13, i23

    PUSH_SUB(fio_index_init)

    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call index_load(this, dir, file, group, ierr)
    if(ierr/=0)then
      message(1) = "Could not read from the input file."
      call messages_fatal(1)
    end if
    call json_get(config, "size", size, ierr)
    ASSERT(ierr==JSON_OK)
    SAFE_ALLOCATE(this%lxyz(size,MAX_DIM))
    i11 = this%nr(1,1)
    i21 = this%nr(2,1)
    i12 = this%nr(1,2)
    i22 = this%nr(2,2)
    i13 = this%nr(1,3)
    i23 = this%nr(2,3)
    SAFE_ALLOCATE(this%lxyz_inv(i11:i21,i12:i22,i13:i23))
    call index_load_lxyz(this, size, dir, group, ierr)
    if(ierr/=0) then
      message(1) = "Could not read from the input file."
      call messages_fatal(1)
    end if
    call checksum_calculate(1, size*this%dim, this%lxyz(1,1), chksm)
    ASSERT(chksm==this%checksum)

    POP_SUB(fio_index_init)
  end subroutine fio_index_init

  ! ---------------------------------------------------------
  subroutine fio_index_copy(this, that)
    type(index_t), intent(inout) :: this
    type(index_t), intent(in)    :: that

    PUSH_SUB(fio_index_copy)

    ASSERT(.false.)
    this = that

    POP_SUB(fio_index_copy)
  end subroutine fio_index_copy

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

end module fio_index_oct_m

!! Local Variables:
!! mode: f90
!! End:
