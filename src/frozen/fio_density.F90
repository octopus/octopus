#include "global.h"

module fio_density_oct_m

  use base_density_oct_m
  use dnst_oct_m
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

  public ::              &
    fio_density__init__, &
    fio_density__load__

contains

  ! ---------------------------------------------------------
  subroutine fio_density__init__(this)
    type(base_density_t), intent(inout) :: this

    PUSH_SUB(fio_density__init__)

    call base_density__init__(this, init)

    POP_SUB(fio_density__init__)

  contains
    
    subroutine init(this, config)
      type(dnst_t),        intent(out) :: this
      type(json_object_t), intent(in)  :: config
      
      real(kind=wp) :: chrb, chrt
      integer       :: ispn, nspin, ierr
      logical       :: load, ipol

      PUSH_SUB(fio_density__init__.init)

      call dnst_init(this, config)
      call json_get(config, "load", load, ierr)
      if(ierr/=JSON_OK) load = .true.
      if(load)then
        call dnst_get(this, nspin=nspin)
        call json_get(config, "invertpolarization", ipol, ierr)
        if(ierr/=JSON_OK) ipol = .false.
        if(ipol.and.(nspin>1))then
          do ispn = 1, nspin/2
            call dnst_get(this, chrb, spin=ispn)
            call dnst_get(this, chrt, spin=nspin-ispn+1)
            call dnst_set(this, chrt, spin=ispn)
            call dnst_set(this, chrb, spin=nspin-ispn+1)
          end do
        end if
      else
        call dnst_reset(this)
      end if

      POP_SUB(fio_density__init__.init)
    end subroutine init
    
  end subroutine fio_density__init__

  ! ---------------------------------------------------------
  subroutine fio_density__read__(this, dir, file, ispin)
    type(base_density_t), intent(inout) :: this
    character(len=*),     intent(in)    :: dir
    character(len=*),     intent(in)    :: file
    integer,              intent(in)    :: ispin

    real(kind=wp), dimension(:,:), pointer :: dnst
    type(simulation_t),            pointer :: sim
    type(mesh_t),                  pointer :: mesh
    character(len=MAX_PATH_LEN)            :: fpth
    integer                                :: ierr

    PUSH_SUB(fio_density__read__)

    nullify(dnst, sim, mesh)
    call path_join(dir, file, fpth)
    call base_density_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, mesh, fine=.true.)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_density_get(this, dnst)
    ASSERT(associated(dnst))
    call dio_function_input(fpth, mesh, dnst(:,ispin), ierr)
    nullify(dnst, mesh)
    if(ierr/=0)then
      call base_density_end(this)
      message(1) = "Could not read the input file: '"//trim(adjustl(file))//".obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

    POP_SUB(fio_density__read__)
  end subroutine fio_density__read__

  ! ---------------------------------------------------------
  subroutine fio_density__load__(this)
    type(base_density_t), intent(inout) :: this

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    character(len=MAX_PATH_LEN)  :: dir, file
    integer                      :: ispn, nspin, ierr
    logical                      :: load, ipol

    PUSH_SUB(fio_density__load__)

    nullify(cnfg, list)
    call base_density_get(this, use=load)
    if(load)then
      call base_density_get(this, cnfg)
      ASSERT(associated(cnfg))
      call json_get(cnfg, "dir", dir, ierr)
      ASSERT(ierr==JSON_OK)
      call json_get(cnfg, "files", list, ierr)
      ASSERT(ierr==JSON_OK)
      call base_density_get(this, nspin=nspin)
      ASSERT(json_len(list)==nspin)
      call json_get(cnfg, "invertpolarization", ipol, ierr)
      if(ierr/=JSON_OK) ipol = .false.
      do ispn = 1, nspin
        call json_get(list, ispn, file, ierr)
        ASSERT(ierr==JSON_OK)
        if(ipol)then
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), nspin-ispn+1)
        else
          call fio_density__read__(this, trim(adjustl(dir)), trim(adjustl(file)), ispn)
        end if
      end do
      nullify(cnfg, list)
      call base_density__update__(this)
    end if
    !call base_density_set(this, static=.true.)

    POP_SUB(fio_density__load__)
  end subroutine fio_density__load__

end module fio_density_oct_m

!! Local Variables:
!! mode: f90
!! End:
