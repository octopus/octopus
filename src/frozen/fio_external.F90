#include "global.h"

module fio_external_oct_m

  use base_potential_oct_m
  use global_oct_m
  use io_function_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use path_oct_m
  use profiling_oct_m
  use simulation_oct_m

  implicit none

  private

  public ::               &
    fio_external__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_external__read__(this, dir, file)
    type(base_potential_t), intent(inout) :: this
    character(len=*),       intent(in)    :: dir
    character(len=*),       intent(in)    :: file

    real(kind=wp), dimension(:), pointer :: potn
    type(simulation_t),          pointer :: sim
    type(mesh_t),                pointer :: mesh
    character(len=MAX_PATH_LEN)          :: fpth
    integer                              :: ierr

    PUSH_SUB(fio_external__read__)

    nullify(potn, sim, mesh)
    call path_join(dir, file, fpth)
    call base_potential_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_potential_get(this, potn)
    ASSERT(associated(potn))
    call dio_function_input(fpth, mesh, potn, ierr)
    nullify(potn, mesh)
    if(ierr/=0)then
      call base_potential_end(this)
      message(1) = "Could not read the input file: '"//trim(adjustl(file))//".obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

    POP_SUB(fio_external__read__)
  end subroutine fio_external__read__

  ! ---------------------------------------------------------
  subroutine fio_external__load__(this)
    type(base_potential_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ierr
    logical                      :: load

    PUSH_SUB(fio_external__load__)

    nullify(cnfg)
    call base_potential_get(this, use=load)
    if(load)then
      call base_potential_get(this, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "dir", dir, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "file", file, ierr)
      ASSERT(ierr==JSON_OK)
      call fio_external__read__(this, trim(adjustl(dir)), trim(adjustl(file)))
      nullify(cnfg)
      call base_potential_update(this)
    end if

    POP_SUB(fio_external__load__)
  end subroutine fio_external__load__

end module fio_external_oct_m

!! Local Variables:
!! mode: f90
!! End:
