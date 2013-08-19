#include "global.h"

module output_fio_m

  use global_m
  use messages_m
  use profiling_m

  use datasets_m,    only: tmpdir
  use epot_m,        only: epot_t
  use geometry_m,    only: geometry_t, geometry_create_data_object
  use grid_m,        only: grid_t
  use hamiltonian_m, only: hamiltonian_t
  use io_m,          only: io_open, io_close
  use json_m,        only: json_init, json_object_t, json_array_t, json_set, json_append, json_write, json_end
  use kinds_m,       only: wp
  use mesh_m,        only: mesh_t
  use simul_box_m,   only: simul_box_t
  use simulation_m,  only: simulation_t
  use space_m,       only: space_create_data_object
  use states_m,      only: states_t

  implicit none

  private
  public ::     &
    output_fio
  
contains

  ! ---------------------------------------------------------
  subroutine fio_simul_box_create_config(this, config)
    type(simul_box_t),   intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    call json_init(config)
    call json_set(config, "dir", "../"//trim(tmpdir)//GS_DIR)
    return
  end subroutine fio_simul_box_create_config

  ! ---------------------------------------------------------
  subroutine fio_mesh_create_config(this, config)
    type(mesh_t),        intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    call json_init(config)
    call json_set(config, "dir", "../"//trim(tmpdir)//GS_DIR)
    call json_set(config, "spacing", this%spacing)
    return
  end subroutine fio_mesh_create_config

  ! ---------------------------------------------------------
  subroutine fio_grid_create_config(this, config)
    type(grid_t),        intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call fio_simul_box_create_config(this%sb, cnfg)
    call json_set(config, "simul_box", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_mesh_create_config(this%mesh, cnfg)
    call json_set(config, "mesh", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_grid_create_config

  ! ---------------------------------------------------------
  subroutine fio_simulation_create_config(this, config)
    type(grid_t),        intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call fio_grid_create_config(this, cnfg)
    call json_set(config, "grid", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_simulation_create_config

  ! ---------------------------------------------------------
  subroutine fio_density_create_config(this, config)
    type(states_t),      intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    character(len=MAX_PATH_LEN) :: filename
    type(json_array_t), pointer :: list
    integer                     :: isp
    !
    nullify(list)
    call json_init(config)
    call json_set(config, "dir", "./")
    call json_set(config, "SpinComponents", this%d%nspin)
    SAFE_ALLOCATE(list)
    call json_init(list)
    do isp = 1, this%d%nspin
      if(this%d%nspin==1) then
        write(unit=filename, fmt='(a)') 'density'
      else
        write(unit=filename, fmt='(a,i1)') 'density-sp', isp
      endif
      call json_append(list, trim(adjustl(filename))//'.obf')
    end do
    call json_set(config, "files", list)
    nullify(list)
    return
  end subroutine fio_density_create_config

  ! ---------------------------------------------------------
  subroutine fio_states_create_config(this, config)
    type(states_t),      intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call fio_density_create_config(this, cnfg)
    call json_set(config, "density", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_states_create_config

  ! ---------------------------------------------------------
  subroutine fio_system_create_config(geo, st, config)
    type(geometry_t),    intent(in)  :: geo
    type(states_t),      intent(in)  :: st
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call space_create_data_object(geo%space, cnfg)
    call json_set(config, "space", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call geometry_create_data_object(geo, cnfg)
    call json_set(config, "geometry", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_states_create_config(st, cnfg)
    call json_set(config, "states", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_system_create_config

  ! ---------------------------------------------------------
  subroutine fio_external_create_config(this, config)
    type(epot_t),        intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    call json_init(config)
    call json_set(config, "dir", "./")
    return
  end subroutine fio_external_create_config

  ! ---------------------------------------------------------
  subroutine fio_hamiltonian_create_config(this, config)
    type(hamiltonian_t), intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call fio_external_create_config(this%ep, cnfg)
    call json_set(config, "external", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_hamiltonian_create_config

  ! ---------------------------------------------------------
  subroutine fio_calc_create_config(gr, geo, st, hm, config)
    type(grid_t),        intent(in)  :: gr
    type(geometry_t),    intent(in)  :: geo
    type(states_t),      intent(in)  :: st
    type(hamiltonian_t), intent(in)  :: hm
    type(json_object_t), intent(out) :: config
    !
    type(json_object_t), pointer :: cnfg
    !
    nullify(cnfg)
    call json_init(config)
    SAFE_ALLOCATE(cnfg)
    call fio_simulation_create_config(gr, cnfg)
    call json_set(config, "simulation", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_system_create_config(geo, st, cnfg)
    call json_set(config, "system", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call fio_hamiltonian_create_config(hm, cnfg)
    call json_set(config, "hamiltonian", cnfg)
    nullify(cnfg)
    return
  end subroutine fio_calc_create_config

  ! ---------------------------------------------------------
  subroutine output_fio(gr, geo, st, hm, file)
    type(grid_t),        intent(in)  :: gr
    type(geometry_t),    intent(in)  :: geo
    type(states_t),      intent(in)  :: st
    type(hamiltonian_t), intent(in)  :: hm
    character(len=*),    intent(in)  :: file
    !
    type(json_object_t) :: config
    integer             :: iunit
    !
    iunit = io_open(file, action='write')
    if(iunit>0)then
      call fio_calc_create_config(gr, geo, st, hm, config)
      call json_write(config, iunit)
      call json_end(config)
      call io_close(iunit)
    else
      message(1)="Could not write the output file: '"//trim(file)//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(2)
    end if
    return
  end subroutine output_fio

end module output_fio_m

!! Local Variables:
!! mode: f90
!! End:
