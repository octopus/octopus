#include "global.h"

module fio_simul_box_m

  use global_m
  use messages_m
  use profiling_m

  use datasets_m,   only: tmpdir
  use io_m,         only: io_open, io_close
  use json_m,       only: JSON_OK, json_object_t, json_get
  use kinds_m,      only: wp
  use kpoints_m,    only: kpoints_init
  use symmetries_m, only: symmetries_init

  use simul_box_m, only:      &
    HYPERCUBE,                &
    simul_box_init_from_file, &
    simul_box_lookup_init

  use simul_box_m, only:                  &
    fio_simul_box_t    => simul_box_t,    &
    fio_simul_box_copy => simul_box_copy, &
    fio_simul_box_end  => simul_box_end

  use fio_geometry_m, only: &
    fio_geometry_t

  implicit none

  private
  public ::             &
    fio_simul_box_t,    &
    fio_simul_box_init, &
    fio_simul_box_copy, &
    fio_simul_box_end
  
contains

  ! ---------------------------------------------------------
  subroutine fio_simul_box_init(this, geo, config)
    type(fio_simul_box_t), intent(out) :: this
    type(fio_geometry_t),  intent(in)  :: geo
    type(json_object_t),   intent(in)  :: config
    !
    character(len=MAX_PATH_LEN)  :: dir
    integer                      :: ierr, iunit, order
    !
    PUSH_SUB(fio_simul_box_init)
    call json_get(config, "dir", dir, ierr)
    if(ierr/=JSON_OK)dir="./"//trim(tmpdir)//GS_DIR
    iunit=io_open(trim(dir)//"mesh", action="read", status="old")
    if(iunit>0)then
      call simul_box_init_from_file(this, iunit)
      call io_close(iunit)
      ASSERT(this%box_shape/=HYPERCUBE)
      this%complex_boundaries=.false.
      call simul_box_lookup_init(this, geo)
      ASSERT(.not.this%mr_flag)
      nullify(this%hr_area%radius, this%hr_area%interp%posi, this%hr_area%interp%ww)
      ASSERT(this%periodic_dim==0)
      call symmetries_init(this%symm, geo, this%dim, this%periodic_dim, this%rlattice, this%lsize)
      call kpoints_init(this%kpoints, this%symm, this%dim, this%rlattice, this%klattice, .true.)
    else
      message(1)="Could not open the simulation box info file: '"//trim(dir)//"mesh'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(2)
    end if
    POP_SUB(fio_simul_box_init)
    return
  end subroutine fio_simul_box_init

end module fio_simul_box_m

!! Local Variables:
!! mode: f90
!! End:
