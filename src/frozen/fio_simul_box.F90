#include "global.h"

module fio_simul_box_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,       only: JSON_OK, json_object_t, json_get
  use kpoints_m,    only: kpoints_init
  use mpi_m,        only: mpi_world
  use symmetries_m, only: symmetries_init

  use geometry_m, only: &
    geometry_t

  use simul_box_m, only:   &
    HYPERCUBE,             &
    simul_box_load,        &
    simul_box_lookup_init

  use bgeom_m, only:           &
    fio_geom_t   => bgeom_t,   &
    fio_geom_get => bgeom_get

  use simul_box_m, only:                  &
    fio_simul_box_t    => simul_box_t,    &
    fio_simul_box_copy => simul_box_copy, &
    fio_simul_box_end  => simul_box_end

  implicit none

  private
  public ::             &
    fio_simul_box_t,    &
    fio_simul_box_init, &
    fio_simul_box_copy, &
    fio_simul_box_end
  
contains

  ! ---------------------------------------------------------
  subroutine fio_simul_box_init(this, geom, config)
    type(fio_simul_box_t), intent(out) :: this
    type(fio_geom_t),      intent(in)  :: geom
    type(json_object_t),   intent(in)  :: config
    !
    type(geometry_t),   pointer :: geo
    character(len=MAX_PATH_LEN) :: dir, file
    integer                     :: ierr, order
    !
    PUSH_SUB(fio_simul_box_init)
    nullify(geo)
    call json_get(config, "dir", dir, ierr)
    ASSERT(ierr==JSON_OK)
    call json_get(config, "file", file, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_geom_get(geom, geo)
    ASSERT(associated(geo))
    call simul_box_load(this, dir, file, mpi_world, ierr)
    if(ierr==0)then
      ASSERT(this%box_shape/=HYPERCUBE)
      call simul_box_lookup_init(this, geo)
      ASSERT(.not.this%mr_flag)
      nullify(this%hr_area%radius, this%hr_area%interp%posi, this%hr_area%interp%ww)
      ASSERT(this%periodic_dim==0)
      call symmetries_init(this%symm, geo, this%dim, this%periodic_dim, this%rlattice)
      call kpoints_init(this%kpoints, this%symm, this%dim, this%rlattice, this%klattice, .true.)
    else
      message(1)="Error reading the simulation box info file: '"//trim(adjustl(file))//"'"
      message(2)="from the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i10)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if
    nullify(geo)
    POP_SUB(fio_simul_box_init)
    return
  end subroutine fio_simul_box_init

end module fio_simul_box_m

!! Local Variables:
!! mode: f90
!! End:
