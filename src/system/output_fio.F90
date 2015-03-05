#include "global.h"

module output_fio_m

  use global_m
  use messages_m
  use profiling_m

  use curvilinear_m, only: curvilinear_t
  use epot_m,        only: epot_t
  use geometry_m,    only: geometry_t, geometry_create_data_object
  use grid_m,        only: grid_t
  use hamiltonian_m, only: hamiltonian_t
  use io_m,          only: io_open, io_close
  use json_m,        only: JSON_OK, json_object_t, json_array_t, json_array_iterator_t
  use json_m,        only: json_init, json_next, json_append, json_set, json_get, json_write, json_end
  use kinds_m,       only: wp
  use mesh_m,        only: mesh_t
  use path_m,        only: path_join
  use simul_box_m,   only: simul_box_t
  use space_m,       only: space_t, space_create_data_object
  use species_m,     only: SPECIES_FROZEN
  use states_m,      only: states_t

  ! use base_hamiltonian_m, only: &
  !   HMLT_TYPE_POTN,             &
  !   HMLT_TYPE_HMLT

  ! use base_config_m, only: &
  !   base_config_parse

  ! use fio_config_m, only: &
  !   input_mesh_dir,       &
  !   input_mesh_file,      &
  !   input_external_file

  implicit none

  private
  public ::     &
    output_fio
  
contains

  ! ---------------------------------------------------------
  subroutine output_parse_config_simul_box(this, sb)
    type(json_object_t), intent(inout) :: this
    type(simul_box_t),   intent(in)    :: sb
    !
    integer :: idim, ierr
    !
    !This routine is not compatible anymore with the current restart machinery and should be rewritten.
    PUSH_SUB(output_parse_config_simul_box)
    call json_get(this, "dimensions", idim, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(idim==sb%dim)
    !call json_set(this, "dir", input_mesh_dir)
    !call json_set(this, "file", input_mesh_file)
    POP_SUB(output_parse_config_simul_box)
    return
  end subroutine output_parse_config_simul_box

  ! ---------------------------------------------------------
  subroutine output_parse_config_curvilinear(this, cv)
    type(json_object_t), intent(inout) :: this
    type(curvilinear_t), intent(in)    :: cv
    !
    PUSH_SUB(output_parse_config_curvilinear)
    call json_set(this, "method", cv%method)
    POP_SUB(output_parse_config_curvilinear)
    return
  end subroutine output_parse_config_curvilinear

  ! ---------------------------------------------------------
  subroutine output_parse_config_mesh(this, mesh)
    type(json_object_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    !
    type(json_array_t), pointer :: list
    integer                     :: indx, ierr
    !
    !This routine is not compatible anymore with the current restart machinery and should be rewritten.
    PUSH_SUB(output_parse_config_mesh)
    nullify(list)
    call json_get(this, "spacing", list, ierr)
    ASSERT(ierr==JSON_OK)
    do indx=1, mesh%sb%dim
      call json_append(list, mesh%spacing(indx))
    end do
    !call json_set(this, "dir", input_mesh_dir)
    !call json_set(this, "file", input_mesh_file)
    POP_SUB(output_parse_config_mesh)
    return
  end subroutine output_parse_config_mesh

  ! ---------------------------------------------------------
  subroutine output_parse_config_grid(this, grid)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_grid)
    nullify(cnfg)
    call json_get(this, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_simul_box(cnfg, grid%sb)
    nullify(cnfg)
    call json_get(this, "curvilinear", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_curvilinear(cnfg, grid%cv)
    nullify(cnfg)
    call json_get(this, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_mesh(cnfg, grid%mesh)
    nullify(cnfg)
    POP_SUB(output_parse_config_grid)
    return
  end subroutine output_parse_config_grid

  ! ---------------------------------------------------------
  subroutine output_parse_config_simulation(this, grid)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_simulation)
    nullify(cnfg)
    call json_get(this, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_grid(cnfg, grid)
    nullify(cnfg)
    POP_SUB(output_parse_config_simulation)
    return
  end subroutine output_parse_config_simulation

  ! ---------------------------------------------------------
  subroutine output_parse_config_space(this, space)
    type(json_object_t), intent(inout) :: this
    type(space_t),       intent(in)    :: space
    !
    PUSH_SUB(output_parse_config_space)
    call space_create_data_object(space, this)
    POP_SUB(output_parse_config_space)
    return
  end subroutine output_parse_config_space

  ! ---------------------------------------------------------
  subroutine output_parse_config_species(this)
    type(json_object_t), intent(inout) :: this
    !
    PUSH_SUB(output_parse_config_species)
    call json_set(this, "type", SPECIES_FROZEN)
    POP_SUB(output_parse_config_species)
    return
  end subroutine output_parse_config_species

  ! ---------------------------------------------------------
  subroutine output_parse_config_geometry(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    !
    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_geometry)
    nullify(cnfg, list)
    call geometry_create_data_object(geo, this)
    call json_get(this, "species", list, ierr)
    ASSERT(ierr==JSON_OK)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call output_parse_config_species(cnfg)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    POP_SUB(output_parse_config_geometry)
    return
  end subroutine output_parse_config_geometry

  ! ---------------------------------------------------------
  subroutine output_parse_config_density(this, st, dir)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    character(len=*),    intent(in)    :: dir
    !
    character(len=MAX_PATH_LEN) :: file
    type(json_array_t), pointer :: list
    integer                     :: ispn, nspin, ierr
    !
    PUSH_SUB(output_parse_config_density)
    nullify(list)
    call json_get(this, "nspin", nspin, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(nspin==st%d%nspin)
    call json_set(this, "dir", trim(adjustl(dir)))
    SAFE_ALLOCATE(list)
    call json_init(list)
    call json_set(this, "files", list)
    if(st%d%nspin>1) then
      do ispn = 1, st%d%nspin
        write(unit=file, fmt="(a,i1,a)") "density-sp", ispn, ".obf"
        call json_append(list, trim(adjustl(file)))
      end do
    else
      call json_append(list, "density.obf")
    end if
    nullify(list)
    POP_SUB(output_parse_config_density)
    return
  end subroutine output_parse_config_density

  ! ---------------------------------------------------------
  subroutine output_parse_config_states(this, st, dir)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    character(len=*),    intent(in)    :: dir
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_states)
    nullify(cnfg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_density(cnfg, st, dir)
    nullify(cnfg)
    POP_SUB(output_parse_config_states)
    return
  end subroutine output_parse_config_states

  ! ---------------------------------------------------------
  subroutine output_parse_config_system(this, geo, st, dir)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    character(len=*),    intent(in)    :: dir
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_system)
    nullify(cnfg)
    call json_get(this, "space", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_end(cnfg)
    call output_parse_config_space(cnfg, geo%space)
    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call json_end(cnfg)
    call output_parse_config_geometry(cnfg, geo)
    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_states(cnfg, st, dir)
    nullify(cnfg)
    POP_SUB(output_parse_config_system)
    return
  end subroutine output_parse_config_system

  ! ---------------------------------------------------------
  subroutine output_parse_config_external(this, epot, dir)
    type(json_object_t), intent(inout) :: this
    type(epot_t),        intent(in)    :: epot
    character(len=*),    intent(in)    :: dir
    !
    integer :: type, ierr
    !
    PUSH_SUB(output_parse_config_external)
    call json_get(this, "type", type, ierr)
    !if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_POTN)
    call json_set(this, "dir", trim(adjustl(dir)))
    !call json_set(this, "file", input_external_file)
    POP_SUB(output_parse_config_external)
    return
  end subroutine output_parse_config_external

  ! ---------------------------------------------------------
  subroutine output_parse_config_hamiltonian(this, hm, dir)
    type(json_object_t), intent(inout) :: this
    type(hamiltonian_t), intent(in)    :: hm
    character(len=*),    intent(in)    :: dir
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: type, ierr
    !
    PUSH_SUB(output_parse_config_hamiltonian)
    nullify(cnfg)
    call json_get(this, "type", type, ierr)
    !if(ierr/=JSON_OK)call json_set(this, "type", HMLT_TYPE_HMLT)
    call json_get(this, "external", cnfg, ierr)
    if(ierr/=JSON_OK)then
      SAFE_ALLOCATE(cnfg)
      call json_init(cnfg)
      call json_set(this, "external", cnfg)
    end if
    call output_parse_config_external(cnfg, hm%ep, dir)
    nullify(cnfg)
    POP_SUB(output_parse_config_hamiltonian)
    return
  end subroutine output_parse_config_hamiltonian

  ! ---------------------------------------------------------
  subroutine output_parse_config_model(this, gr, geo, st, hm, dir)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: hm
    character(len=*),    intent(in)    :: dir
    !
    type(json_object_t), pointer :: cnfg
    integer                      :: ierr
    !
    PUSH_SUB(output_parse_config_model)
    nullify(cnfg)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_simulation(cnfg, gr)
    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_system(cnfg, geo, st, dir)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_hamiltonian(cnfg, hm, dir)
    nullify(cnfg)
    POP_SUB(output_parse_config_model)
    return
  end subroutine output_parse_config_model

  ! ---------------------------------------------------------
  subroutine output_fio(gr, geo, st, hm, dir, file)
    type(grid_t),        intent(in) :: gr
    type(geometry_t),    intent(in) :: geo
    type(states_t),      intent(in) :: st
    type(hamiltonian_t), intent(in) :: hm
    character(len=*),    intent(in) :: dir
    character(len=*),    intent(in) :: file
    !
    type(json_object_t)          :: config
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: fnam
    integer                      :: iunit, ierr
    !
    PUSH_SUB(output_fio)
    nullify(cnfg)
    call path_join(dir, file, fnam)
    iunit=io_open(fnam, action='write')
    if(iunit>0)then
      !call base_config_parse(config, geo%space%dim, st%d%nspin)
      call json_get(config, "model", cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call output_parse_config_model(cnfg, gr, geo, st, hm, dir)
      call json_write(cnfg, iunit)
      nullify(cnfg)
      call json_end(config)
      call io_close(iunit)
    else
      message(1)="Could not write the output file: '"//trim(adjustl(fnam))//"'"
      write(unit=message(2), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(2)
    end if
    POP_SUB(output_fio)
    return
  end subroutine output_fio

end module output_fio_m

!! Local Variables:
!! mode: f90
!! End:
