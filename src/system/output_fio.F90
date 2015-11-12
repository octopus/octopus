#include "global.h"

module output_fio_m

  use base_config_m
  use base_hamiltonian_m  
  use curvilinear_m
  use epot_m
  use fio_config_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_m
  use io_function_m
  use json_m
  use json_m
  use kinds_m
  use mesh_m
  use messages_m
  use mpi_m
  use path_m
  use profiling_m
  use simul_box_m
  use space_m
  use species_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m
  
  implicit none

  private

  public ::     &
    output_fio

  integer, parameter :: output_format = OPTION__OUTPUTFORMAT__BINARY

  character(len=*), parameter :: grid_dump_file = "grid"

contains

  ! ---------------------------------------------------------
  subroutine output_parse_config_simul_box(this, sb, dir, file, group)
    type(json_object_t), intent(out) :: this
    type(simul_box_t),   intent(in)  :: sb
    character(len=*),    intent(in)  :: dir
    character(len=*),    intent(in)  :: file
    type(mpi_grp_t),     intent(in)  :: group

    integer :: ierr

    PUSH_SUB(output_parse_config_simul_box)

    ASSERT(sb%dim>0)
    call json_init(this)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", trim(adjustl(file)))
    call simul_box_dump(sb, dir, file, group, ierr)
    if(ierr/=0)then
      message(1) = "Could not write to the output file."
      call messages_fatal(1)
    end if

    POP_SUB(output_parse_config_simul_box)
  end subroutine output_parse_config_simul_box

  ! ---------------------------------------------------------
  subroutine output_parse_config_curvilinear(this, cv)
    type(json_object_t), intent(out) :: this
    type(curvilinear_t), intent(in)  :: cv

    PUSH_SUB(output_parse_config_curvilinear)

    ASSERT(cv%method==CURV_METHOD_UNIFORM)
    call json_init(this)
    call json_set(this, "method", cv%method)

    POP_SUB(output_parse_config_curvilinear)
  end subroutine output_parse_config_curvilinear

  ! ---------------------------------------------------------
  subroutine output_parse_config_index(this, index, np, dir, file, group)
    type(json_object_t), intent(out) :: this
    type(index_t),       intent(in)  :: index
    integer,             intent(in)  :: np
    character(len=*),    intent(in)  :: dir
    character(len=*),    intent(in)  :: file
    type(mpi_grp_t),     intent(in)  :: group

    integer :: ierr

    PUSH_SUB(output_parse_config_index)

    call json_init(this)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", trim(adjustl(file)))
    call index_dump(index, dir, file, group, ierr)
    if(ierr/=0)then
      message(1) = "Could not write to the output file."
      call messages_fatal(1)
    end if
    call json_set(this, "size", np)
    call index_dump_lxyz(index, np, dir, group, ierr)
    if(ierr/=0)then
      message(1) = "Could not write to the output file."
      call messages_fatal(1)
    end if

    POP_SUB(output_parse_config_index)
  end subroutine output_parse_config_index

  ! ---------------------------------------------------------
  subroutine output_parse_config_mesh(this, mesh, dir, file, group)
    type(json_object_t), intent(out) :: this
    type(mesh_t),        intent(in)  :: mesh
    character(len=*),    intent(in)  :: dir
    character(len=*),    intent(in)  :: file
    type(mpi_grp_t),     intent(in)  :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_mesh)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", trim(adjustl(file)))
    call json_set(this, "localsize", mesh%np)
    call json_set(this, "localfullsize", mesh%np_part)
    call json_set(this, "globalsize", mesh%np_global)
    call json_set(this, "globalfullsize", mesh%np_part_global)
    call json_set(this, "spacing", mesh%spacing)
    call mesh_write_fingerprint(mesh, dir, file, group, ierr)
    if(ierr/=0)then
      message(1) = "Could not write to the output file."
      call messages_fatal(1)
    end if
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_index(cnfg, mesh%idx, mesh%np_part_global, dir, file, group)
    call json_set(this, "index", cnfg)
    nullify(cnfg)

    POP_SUB(output_parse_config_mesh)
  end subroutine output_parse_config_mesh

  ! ---------------------------------------------------------
  subroutine output_parse_config_grid(this, grid, dir, file, group)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    character(len=*),    intent(in)    :: dir
    character(len=*),    intent(in)    :: file
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(output_parse_config_grid)

    nullify(cnfg)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", trim(adjustl(file)))
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_mesh(cnfg, grid%mesh, dir, file, group)
    call json_set(this, "mesh", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_simul_box(cnfg, grid%sb, dir, file, group)
    call json_set(this, "simul_box", cnfg)
    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_curvilinear(cnfg, grid%cv)
    call json_set(this, "curvilinear", cnfg)
    nullify(cnfg)

    POP_SUB(output_parse_config_grid)
  end subroutine output_parse_config_grid

  ! ---------------------------------------------------------
  subroutine output_parse_config_simulation(this, grid, dir, group)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_simulation)

    nullify(cnfg)
    call json_get(this, "grid", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_grid(cnfg, grid, dir, grid_dump_file, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_simulation)
  end subroutine output_parse_config_simulation

  ! ---------------------------------------------------------
  subroutine output_parse_config_species(this)
    type(json_object_t), intent(inout) :: this

    PUSH_SUB(output_parse_config_species)

    call json_set(this, "type", SPECIES_FROZEN)

    POP_SUB(output_parse_config_species)
  end subroutine output_parse_config_species

  ! ---------------------------------------------------------
  subroutine output_parse_config_molecule(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(json_array_iterator_t)  :: iter
    integer                      :: ierr

    PUSH_SUB(output_parse_config_molecule)

    nullify(cnfg, list)
    call json_end(this)
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

    POP_SUB(output_parse_config_molecule)
  end subroutine output_parse_config_molecule

  ! ---------------------------------------------------------
  subroutine output_parse_config_geometry(this, geo)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_geometry)

    nullify(cnfg)
    call json_get(this, "molecule", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_molecule(cnfg, geo)
    nullify(cnfg)
    
    POP_SUB(output_parse_config_geometry)
  end subroutine output_parse_config_geometry

  ! ---------------------------------------------------------
  subroutine output_parse_config_density_charge(this, st)
    type(json_array_t), intent(inout) :: this
    type(states_t),     intent(in)    :: st

    real(kind=wp) :: chrg, qtot
    integer       :: ispn, ik, ierr

    PUSH_SUB(output_parse_config_density_charge)

    qtot = 0.0_wp
    do ispn = 1, st%d%nspin
      chrg = 0.0_wp
      do ik = 1, st%d%nik
        if(ispn==states_dim_get_spin_index(st%d, ik))&
          chrg = chrg + st%d%kweights(ik) * sum(st%occ(:,ik))
      end do
      call json_set(this, ispn, chrg, ierr)
      ASSERT(ierr==JSON_OK)
      qtot = qtot + chrg
    end do
    ASSERT(.not.abs(1.0_wp-min(st%qtot,qtot)/max(st%qtot,qtot))>epsilon(qtot))

    POP_SUB(output_parse_config_density_charge)
  end subroutine output_parse_config_density_charge

  ! ---------------------------------------------------------
  subroutine output_parse_config_density_files(this, st, gr, geo, dir, group)
    type(json_array_t), intent(out) :: this
    type(states_t),     intent(in)  :: st
    type(grid_t),       intent(in)  :: gr
    type(geometry_t),   intent(in)  :: geo
    character(len=*),   intent(in)  :: dir
    type(mpi_grp_t),    intent(in)  :: group

    character(len=MAX_PATH_LEN) :: file
    type(unit_t)                :: unts
    integer                     :: ispn, ierr

    PUSH_SUB(output_parse_config_density_files)

    call json_init(this)
    unts = units_out%length**(-gr%mesh%sb%dim)
    do ispn = 1, st%d%nspin
      if(st%d%nspin>1) then
        write(unit=file, fmt="(a,i1)") "density-sp", ispn
      else
        file="density"
      end if
      call json_append(this, trim(adjustl(file)))
      call dio_function_output(output_format, dir, file, gr%fine%mesh, &
        st%rho(:,ispn), unts, ierr, geo = geo, grp = group)
      if(ierr/=0)then
        message(1) = "Could not write to the output file: '"//trim(adjustl(file))//".obf'"
        message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
        write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
        call messages_fatal(3)
      end if
    end do

    POP_SUB(output_parse_config_density_files)
  end subroutine output_parse_config_density_files

  ! ---------------------------------------------------------
  subroutine output_parse_config_density(this, st, gr, geo, dir, group)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_array_t), pointer :: list
    integer                     :: ierr

    PUSH_SUB(output_parse_config_density)

    nullify(list)
    ASSERT(st%d%nspin>0)
    ASSERT(st%d%nspin<3)
    call json_get(this, "charge", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(st%d%nspin==json_len(list))
    call output_parse_config_density_charge(list, st)
    nullify(list)
    call json_set(this, "dir", trim(adjustl(dir)))
    SAFE_ALLOCATE(list)
    call output_parse_config_density_files(list, st, gr, geo, dir, group)
    call json_set(this, "files", list)
    nullify(list)

    POP_SUB(output_parse_config_density)
  end subroutine output_parse_config_density

  ! ---------------------------------------------------------
  subroutine output_parse_config_states(this, st, gr, geo, dir, group)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_states)

    nullify(cnfg)
    ASSERT(st%qtot>0.0_wp)
    call json_set(this, "charge", st%qtot)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_density(cnfg, st, gr, geo, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_states)
  end subroutine output_parse_config_states

  ! ---------------------------------------------------------
  subroutine output_parse_config_system(this, geo, st, gr, dir, group)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    type(grid_t),        intent(in)    :: gr
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_system)

    nullify(cnfg)
    call json_get(this, "geometry", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_geometry(cnfg, geo)
    nullify(cnfg)
    call json_get(this, "states", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_states(cnfg, st, gr, geo, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_system)
  end subroutine output_parse_config_system

  ! ---------------------------------------------------------
  subroutine output_parse_config_external(this, epot, gr, geo, dir, group)
    type(json_object_t), intent(out) :: this
    type(epot_t),        intent(in)  :: epot
    type(grid_t),        intent(in)  :: gr
    type(geometry_t),    intent(in)  :: geo
    character(len=*),    intent(in)  :: dir
    type(mpi_grp_t),     intent(in)  :: group

    integer :: ierr

    PUSH_SUB(output_parse_config_external)

    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", "v0")
    call dio_function_output(output_format, dir, "v0", &
      gr%mesh, epot%vpsl, units_out%energy, ierr, geo = geo, grp = group)
    if(ierr/=0)then
      message(1) = "Could not write to the output file: 'v0.obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if

    POP_SUB(output_parse_config_external)
  end subroutine output_parse_config_external

  ! ---------------------------------------------------------
  subroutine output_parse_config_hamiltonian(this, hm, gr, geo, dir, group)
    type(json_object_t), intent(inout) :: this
    type(hamiltonian_t), intent(in)    :: hm
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(output_parse_config_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_external(cnfg, hm%ep, gr, geo, dir, group)
    call json_set(this, "external", cnfg)
    nullify(cnfg)

    POP_SUB(output_parse_config_hamiltonian)
  end subroutine output_parse_config_hamiltonian

  ! ---------------------------------------------------------
  subroutine output_parse_config_model(this, gr, geo, st, hm, dir, group)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: gr
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: hm
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_model)

    nullify(cnfg)
    call json_get(this, "simulation", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_simulation(cnfg, gr, dir, group)
    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_system(cnfg, geo, st, gr, dir, group)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_hamiltonian(cnfg, hm, gr, geo, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_model)
  end subroutine output_parse_config_model

  ! ---------------------------------------------------------
  subroutine output_fio(gr, geo, st, hm, dir, group)
    type(grid_t),        intent(in) :: gr
    type(geometry_t),    intent(in) :: geo
    type(states_t),      intent(in) :: st
    type(hamiltonian_t), intent(in) :: hm
    character(len=*),    intent(in) :: dir
    type(mpi_grp_t),     intent(in) :: group

    type(json_object_t)          :: config
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: dnam, fnam
    integer                      :: iunit, ierr

    PUSH_SUB(output_fio)

    nullify(cnfg)
    call path_join(dir, input_frozen_dir, dnam)
    call io_mkdir(dnam)
    call path_join(dnam, input_frozen_config, fnam)
    iunit = io_open(fnam, action='write')
    if(iunit>0)then
      call base_config_parse(config, geo%space%dim, st%d%nspin)
      call json_get(config, "model", cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call output_parse_config_model(cnfg, gr, geo, st, hm, dnam, group)
      call json_write(cnfg, iunit)
      nullify(cnfg)
      call json_end(config)
      call io_close(iunit)
    else
      message(1) = "Could not write to the output file: '"//trim(adjustl(input_frozen_config))//"'"
      message(2) = "in the directory: '"//trim(adjustl(dnam))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(3)
    end if

    POP_SUB(output_fio)
  end subroutine output_fio

end module output_fio_m

!! Local Variables:
!! mode: f90
!! End:
