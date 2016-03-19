#include "global.h"

module fio_simul_box_oct_m

  use geometry_oct_m
  use global_oct_m
  use json_oct_m
  use kpoints_oct_m
  use mpi_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use symmetries_oct_m

  implicit none

  private

  public ::             &
    fio_simul_box_init, &
    fio_simul_box_copy, &
    fio_simul_box_end
  
contains

  ! ---------------------------------------------------------
  subroutine fio_simul_box_init(this, geo, space, config)
    type(simul_box_t),   intent(out) :: this
    type(geometry_t),    intent(in)  :: geo
    type(space_t),       intent(in)  :: space
    type(json_object_t), intent(in)  :: config

    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: ierr

    PUSH_SUB(fio_simul_box_init)

    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call simul_box_load(this, dir, file, mpi_world, ierr)
    if(ierr==0)then
      ASSERT(this%dim==space%dim)
      ASSERT(this%box_shape/=HYPERCUBE)
      call simul_box_lookup_init(this, geo)
      ASSERT(.not.this%mr_flag)
      nullify(this%hr_area%radius, this%hr_area%interp%posi, this%hr_area%interp%ww)
      ASSERT(this%periodic_dim==0)
      call symmetries_init(this%symm, geo, this%dim, this%periodic_dim, this%rlattice)
      call kpoints_init(this%kpoints, this%symm, this%dim, this%rlattice, this%klattice, .true.)
    else
      message(1) = "Could not read from the input file."
      call messages_fatal(1)
    end if

    POP_SUB(fio_simul_box_init)
  end subroutine fio_simul_box_init

  ! ---------------------------------------------------------
  subroutine fio_simul_box_copy(this, that)
    type(simul_box_t), intent(inout) :: this
    type(simul_box_t), intent(in)    :: that

    PUSH_SUB(fio_simul_box_copy)

    call simul_box_copy(this, that)

    POP_SUB(fio_simul_box_copy)
  end subroutine fio_simul_box_copy

  ! ---------------------------------------------------------
  subroutine fio_simul_box_end(this)
    type(simul_box_t), intent(inout) :: this

    PUSH_SUB(fio_simul_box_end)

    call simul_box_end(this)

    POP_SUB(fio_simul_box_end)
  end subroutine fio_simul_box_end

end module fio_simul_box_oct_m

!! Local Variables:
!! mode: f90
!! End:
