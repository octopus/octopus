#include "global.h"

module output_fio_oct_m

  use base_config_oct_m
  use base_hamiltonian_oct_m  
  use batch_oct_m
  use curvilinear_oct_m
  use density_oct_m
  use epot_oct_m
  use fio_config_oct_m
  use fio_handle_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use index_oct_m
  use io_oct_m
  use io_function_oct_m
  use json_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use path_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use storage_oct_m
  use unit_oct_m
  use unit_system_oct_m
  
  implicit none

  private

  public ::     &
    output_fio

  integer, parameter :: output_format = OPTION__OUTPUTFORMAT__BINARY

  character(len=*), parameter :: grid_dump_file = "grid"

  integer, parameter :: TYPE_NULL = 0
  integer, parameter :: TYPE_ASSC = 1
  integer, parameter :: TYPE_ALLC = 2

  type :: density_batch_t
    private
    integer                                         :: type = TYPE_NULL
    type(batch_t),                          pointer :: psib =>null()
    real(kind=wp),    dimension(:,:,:), allocatable :: dpsi
    complex(kind=wp), dimension(:,:,:), allocatable :: zpsi
  end type density_batch_t

  type :: density_acc_t
    private
    integer                 :: nst  = 0
    type(states_t), pointer :: st   =>null()
    type(grid_t),   pointer :: grid =>null()
    type(density_calc_t)    :: calc
  end type density_acc_t

contains

  ! ---------------------------------------------------------
  subroutine density_batch_init(this, st, grid, ik, ib, nst)
    type(density_batch_t), intent(out) :: this
    type(states_t),        intent(in)  :: st
    type(grid_t),          intent(in)  :: grid
    integer,               intent(in)  :: ik
    integer,               intent(in)  :: ib
    integer,               intent(in)  :: nst

    integer :: ist, nblk

    PUSH_SUB(density_batch_init)

    this%type = TYPE_NULL
    nullify(this%psib)
    if(states_block_min(st,ib)<=nst)then
      if(states_block_max(st,ib)<=nst)then
        this%type = TYPE_ASSC
        this%psib => st%group%psib(ib,ik)
      else
        this%type = TYPE_ALLC
        SAFE_ALLOCATE(this%psib)
        nblk = nst - states_block_min(st,ib) + 1
        call batch_init(this%psib, st%d%dim, nblk)
        if(states_are_real(st))then
          SAFE_ALLOCATE(this%dpsi(1:grid%mesh%np,1:st%d%dim,states_block_min(st,ib):nst))
          do ist = states_block_min(st,ib), nst
            call states_get_state(st, grid%mesh, ist, ik, this%dpsi(:,:,ist))
            call batch_add_state(this%psib, ist, this%dpsi(:,:,ist))
          end do
        else
          SAFE_ALLOCATE(this%zpsi(1:grid%mesh%np,1:st%d%dim,states_block_min(st,ib):nst))
          do ist = states_block_min(st,ib), nst
            call states_get_state(st, grid%mesh, ist, ik, this%zpsi(:,:,ist))
            call batch_add_state(this%psib, ist, this%zpsi(:,:,ist))
          end do
        end if
      end if
    end if

    POP_SUB(density_batch_init)
  end subroutine density_batch_init

  ! ---------------------------------------------------------
  subroutine density_batch_get(this, that)
    type(density_batch_t), intent(in) :: this
    type(batch_t),        pointer     :: that

    PUSH_SUB(density_batch_get)

    nullify(that)
    if(associated(this%psib)) that => this%psib

    POP_SUB(density_batch_get)
  end subroutine density_batch_get

  ! ---------------------------------------------------------
  subroutine density_batch_end(this)
    type(density_batch_t), intent(inout) :: this

    PUSH_SUB(density_batch_end)

    if(this%type==TYPE_ALLC)then
      call batch_end(this%psib)
      SAFE_DEALLOCATE_P(this%psib)
      if(allocated(this%dpsi))then
        SAFE_DEALLOCATE_A(this%dpsi)
      else
        SAFE_DEALLOCATE_A(this%zpsi)
      end if
    end if
    nullify(this%psib)
    this%type = TYPE_NULL

    POP_SUB(density_batch_end)
  end subroutine density_batch_end

  ! ---------------------------------------------------------
  subroutine density_acc_init(this, st, grid, nst, density)
    type(density_acc_t),           intent(out) :: this
    type(states_t),        target, intent(in)  :: st
    type(grid_t),          target, intent(in)  :: grid
    integer,                       intent(in)  :: nst
    real(kind=wp), dimension(:,:), intent(out) :: density

    PUSH_SUB(density_acc_init)

    this%nst = nst
    this%st => st
    this%grid => grid
    call density_calc_init(this%calc, this%st, this%grid, density)

    POP_SUB(density_acc_init)
  end subroutine density_acc_init

  ! ---------------------------------------------------------
  subroutine density_acc(this, ik, ib)
    type(density_acc_t), intent(inout) :: this
    integer,             intent(in)    :: ik
    integer,             intent(in)    :: ib

    type(density_batch_t)  :: btch
    type(batch_t), pointer :: psib

    PUSH_SUB(density_acc)

    nullify(psib)
    call density_batch_init(btch, this%st, this%grid, ik, ib, this%nst)
    call density_batch_get(btch, psib)
    if(associated(psib)) call density_calc_accumulate(this%calc, ik, psib)
    nullify(psib)
    call density_batch_end(btch)

    POP_SUB(density_acc)
  end subroutine density_acc

  ! ---------------------------------------------------------
  subroutine density_acc_end(this)
    type(density_acc_t), intent(inout) :: this

    PUSH_SUB(density_acc_end)

    this%nst = 0
    nullify(this%st, this%grid)
    call density_calc_end(this%calc)

    POP_SUB(density_acc_end)
  end subroutine density_acc_end

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
    call json_set(this, "spacing", mesh%spacing(1:mesh%sb%dim))
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
  subroutine output_parse_config_density_charge(this, st, nst, charge)
    type(json_array_t), intent(inout) :: this
    type(states_t),     intent(in)    :: st
    integer,            intent(in)    :: nst
    real(kind=wp),      intent(in)    :: charge

    real(kind=wp) :: chrg, qtot
    integer       :: ispn, ik, ierr

    PUSH_SUB(output_parse_config_density_charge)

    qtot = 0.0_wp
    do ispn = 1, st%d%nspin
      chrg = 0.0_wp
      do ik = 1, st%d%nik
        if(ispn==states_dim_get_spin_index(st%d, ik))&
          chrg = chrg + st%d%kweights(ik) * sum(st%occ(1:nst,ik))
      end do
      call json_set(this, ispn, chrg, ierr)
      ASSERT(ierr==JSON_OK)
      qtot = qtot + chrg
    end do
    ASSERT(.not.abs(1.0_wp-min(charge,qtot)/max(charge,qtot))>epsilon(qtot))

    POP_SUB(output_parse_config_density_charge)
  end subroutine output_parse_config_density_charge

  ! ---------------------------------------------------------
  subroutine output_parse_config_density_calc(st, grid, nst, density)
    type(states_t),                intent(inout) :: st
    type(grid_t),                  intent(in)    :: grid
    integer,                       intent(in)    :: nst
    real(kind=wp), dimension(:,:), intent(out)   :: density

    type(density_acc_t) :: calc
    integer             :: ik, ib

    PUSH_SUB(output_parse_config_density_calc)

    call density_acc_init(calc, st, grid, nst, density)
    if((st%st_start<=nst).and.(nst<=st%st_end))then
      do ik = st%d%kpt%start, st%d%kpt%end
        do ib = st%group%block_start, st%group%block_end
          call density_acc(calc, ik, ib)
        end do
      end do
    end if
    call density_acc_end(calc)

    POP_SUB(output_parse_config_density_calc)
  end subroutine output_parse_config_density_calc

  ! ---------------------------------------------------------
  subroutine output_parse_config_density_files(this, st, grid, geo, nst, dir, group)
    type(json_array_t), intent(out)   :: this
    type(states_t),     intent(inout) :: st
    type(grid_t),       intent(in)    :: grid
    type(geometry_t),   intent(in)    :: geo
    integer,            intent(in)    :: nst
    character(len=*),   intent(in)    :: dir
    type(mpi_grp_t),    intent(in)    :: group

    real(kind=wp), dimension(:,:), allocatable :: density
    character(len=MAX_PATH_LEN)                :: file
    type(unit_t)                               :: unts
    integer                                    :: ispn, ierr

    PUSH_SUB(output_parse_config_density_files)

    call json_init(this)
    unts = units_out%length**(-grid%mesh%sb%dim)
    SAFE_ALLOCATE(density(1:grid%fine%mesh%np_part,1:st%d%nspin))
    call output_parse_config_density_calc(st, grid, nst, density)
    do ispn = 1, st%d%nspin
      if(st%d%nspin>1) then
        write(unit=file, fmt="(a,i1)") "density-sp", ispn
      else
        file="density"
      end if
      call json_append(this, trim(adjustl(file)))
      call dio_function_output(output_format, dir, file, grid%fine%mesh, &
        density(:,ispn), unts, ierr, geo = geo, grp = group)
      if(ierr/=0)then
        message(1) = "Could not write to the output file: '"//trim(adjustl(file))//".obf'"
        message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
        write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
        call messages_fatal(3)
      end if
    end do
    SAFE_DEALLOCATE_A(density)

    POP_SUB(output_parse_config_density_files)
  end subroutine output_parse_config_density_files

  ! ---------------------------------------------------------
  subroutine output_parse_config_density(this, st, grid, geo, nst, charge, dir, group)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: grid
    type(geometry_t),    intent(in)    :: geo
    integer,             intent(in)    :: nst
    real(kind=wp),       intent(in)    :: charge
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
    call output_parse_config_density_charge(list, st, nst, charge)
    nullify(list)
    call json_set(this, "dir", trim(adjustl(dir)))
    SAFE_ALLOCATE(list)
    call output_parse_config_density_files(list, st, grid, geo, nst, dir, group)
    call json_set(this, "files", list)
    nullify(list)

    POP_SUB(output_parse_config_density)
  end subroutine output_parse_config_density

  ! ---------------------------------------------------------
  subroutine output_parse_config_states(this, st, grid, geo, dir, group)
    type(json_object_t), intent(inout) :: this
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: grid
    type(geometry_t),    intent(in)    :: geo
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    real(kind=wp)                :: chrg
    integer                      :: fnst, ierr
    integer                      :: indx, jndx

    PUSH_SUB(output_parse_config_states)

    nullify(cnfg)
    ASSERT(st%qtot>0.0_wp)
    !%Variable FrozenStates
    !%Type integer
    !%Default Number of States
    !%Section Output::Subsystems
    !%Description
    !% Subdirectory of static/frozen where <tt>Octopus</tt> will write the frozen system.
    !%End
    call parse_variable('FrozenStates', st%nst, fnst)
    if(fnst<0) fnst = fnst + st%nst
    if((fnst>st%nst).or.(fnst<1))then
      write(message(1),'(a,i6,a)') "Can not freeze ", fnst, " states."
      write(message(2),'(a,i6,a)') "Can not freeze less than 1 or more than ", st%nst, "states."
      call messages_fatal(2)
    end if
    chrg = 0.0_wp
    do indx = 1, fnst
      do jndx = 1, st%d%nik
        chrg = chrg + st%occ(indx,jndx) * st%d%kweights(jndx)
      end do
    end do
    call json_set(this, "charge", chrg)
    call json_get(this, "density", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_density(cnfg, st, grid, geo, fnst, chrg, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_states)
  end subroutine output_parse_config_states

  ! ---------------------------------------------------------
  subroutine output_parse_config_system(this, geo, st, grid, dir, group)
    type(json_object_t), intent(inout) :: this
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: grid
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
    call output_parse_config_states(cnfg, st, grid, geo, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_system)
  end subroutine output_parse_config_system

  ! ---------------------------------------------------------
  subroutine output_parse_config_external(this, epot, grid, geo, dir, group)
    type(json_object_t), intent(out) :: this
    type(epot_t),        intent(in)  :: epot
    type(grid_t),        intent(in)  :: grid
    type(geometry_t),    intent(in)  :: geo
    character(len=*),    intent(in)  :: dir
    type(mpi_grp_t),     intent(in)  :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_external)

    nullify(cnfg)
    call json_init(this)
    call json_set(this, "type", HMLT_TYPE_POTN)
    call json_set(this, "dir", trim(adjustl(dir)))
    call json_set(this, "file", "v0")
    call dio_function_output(output_format, dir, "v0", &
      grid%mesh, epot%vpsl, units_out%energy, ierr, geo = geo, grp = group)
    if(ierr/=0)then
      message(1) = "Could not write to the output file: 'v0.obf'"
      message(2) = "in the directory: '"//trim(adjustl(dir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", ierr
      call messages_fatal(3)
    end if
    SAFE_ALLOCATE(cnfg)
    call storage_init(cnfg, full=.false.)
    call json_set(this, "storage", cnfg)
    nullify(cnfg)

    POP_SUB(output_parse_config_external)
  end subroutine output_parse_config_external

  ! ---------------------------------------------------------
  subroutine output_parse_config_hamiltonian(this, hm, grid, geo, dir, group)
    type(json_object_t), intent(inout) :: this
    type(hamiltonian_t), intent(in)    :: hm
    type(grid_t),        intent(in)    :: grid
    type(geometry_t),    intent(in)    :: geo
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg

    PUSH_SUB(output_parse_config_hamiltonian)

    nullify(cnfg)
    SAFE_ALLOCATE(cnfg)
    call output_parse_config_external(cnfg, hm%ep, grid, geo, dir, group)
    call json_set(this, "external", cnfg)
    nullify(cnfg)

    POP_SUB(output_parse_config_hamiltonian)
  end subroutine output_parse_config_hamiltonian

  ! ---------------------------------------------------------
  subroutine output_parse_config_model(this, grid, geo, st, hm, dir, group)
    type(json_object_t), intent(inout) :: this
    type(grid_t),        intent(in)    :: grid
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(in)    :: hm
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(output_parse_config_model)

    nullify(cnfg)
    call json_get(this, "system", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_system(cnfg, geo, st, grid, dir, group)
    nullify(cnfg)
    call json_get(this, "hamiltonian", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call output_parse_config_hamiltonian(cnfg, hm, grid, geo, dir, group)
    nullify(cnfg)

    POP_SUB(output_parse_config_model)
  end subroutine output_parse_config_model

  ! ---------------------------------------------------------
  subroutine output_fio(grid, geo, st, hm, dir, group)
    type(grid_t),        intent(in)    :: grid
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(in)    :: hm
    character(len=*),    intent(in)    :: dir
    type(mpi_grp_t),     intent(in)    :: group

    type(json_object_t)          :: config
    type(json_object_t), pointer :: cnfg
    character(len=MAX_PATH_LEN)  :: tdir, bdir, fdir, odir, fnam
    integer                      :: iunit, ierr

    PUSH_SUB(output_fio)

    nullify(cnfg)
    call path_realpath(dir, bdir)
    !%Variable FrozenDir
    !%Type string
    !%Default "static/frozen"
    !%Section Output::Subsystems
    !%Description
    !% Subdirectory of static/ where <tt>Octopus</tt> will write the frozen system.
    !%End
    call parse_variable('FrozenDir', input_frozen_dir, tdir)
    if(path_isabs(tdir))then
      fdir = tdir
    else
      call path_join(bdir, tdir, fdir)
    end if
    call io_mkdir(fdir)
    call path_realpath(fdir, odir)
    call path_join(odir, input_frozen_config, fnam)
    iunit = io_open(fnam, action='write')
    if(iunit>0)then
      call base_config_parse(config, st%d%nspin, geo%space%dim)
      call json_set(config, "type", HNDL_TYPE_FNIO)
      call json_set(config, "name", "fio")
      call json_get(config, "simulation", cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call output_parse_config_simulation(cnfg, grid, odir, group)
      nullify(cnfg)
      call json_get(config, "model", cnfg, ierr)
      ASSERT(ierr==JSON_OK)
      call output_parse_config_model(cnfg, grid, geo, st, hm, odir, group)
      nullify(cnfg)
      call json_write(config, iunit)
      call json_end(config)
      call io_close(iunit)
    else
      message(1) = "Could not write to the output file: '"//trim(adjustl(input_frozen_config))//"'"
      message(2) = "in the directory: '"//trim(adjustl(odir))//"'"
      write(unit=message(3), fmt="(a,i3)") "I/O Error: ", iunit
      call messages_fatal(3)
    end if

    POP_SUB(output_fio)
  end subroutine output_fio

end module output_fio_oct_m

!! Local Variables:
!! mode: f90
!! End:
