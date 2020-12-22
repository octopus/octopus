!! Copyright (C) 2014 M. Oliveira, J. Jornet-Somoza
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

program oct_local_multipoles
  use atom_oct_m
  use basins_oct_m
  use box_oct_m
  use box_union_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use kick_oct_m
  use loct_oct_m
  use local_write_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use simul_box_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none
  
  type local_domain_t
    type(box_union_t)    :: domain         !< boxes that form each domain.
    character(len=500)   :: clist          !< list of centers for each domain.
    character(len=15)    :: lab            !< declared name for each domain.
    integer              :: dshape         !< shape of box for each domain.
    logical, allocatable :: inside(:)      !< relation of mesh points on each domain.
    logical, allocatable :: ions_inside(:) !< relation of ions inside each domain.
    FLOAT,   allocatable :: dcm(:)         !< store the center of mass of each domain on the real space.
    type(local_write_t)  :: writ           !< write option for local domains analysis.
  end type local_domain_t

  type(electrons_t), pointer :: sys
  integer, parameter    :: BADER = 512
  FLOAT                 :: BaderThreshold

  ! Initialize stuff
  call global_init(is_serial = .false.)
  call calc_mode_par_init()

  call parser_init()
  
  call messages_init()

  call io_init()
  call profiling_init(global_namespace)
 
  call print_header()
  call messages_print_stress(stdout, "Local Domains mode")
  call messages_print_stress(stdout)
    
  call messages_experimental("oct-local_multipoles utility")

  call unit_system_init(global_namespace)
  call restart_module_init(global_namespace)

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  sys => electrons_t(global_namespace)
  call sys%init_parallelization(mpi_world)

  call local_domains()

  SAFE_DEALLOCATE_P(sys)
  call profiling_end(global_namespace)
  call io_end()
  call print_date("Calculation ended on ")
  call messages_end()
  call parser_end()
  call global_end()

contains
  ! -------------
  !> Reads a global grid and density and make local domains
  !! This is a high-level interface that reads the input file and
  !! calls the proper function.
  subroutine local_domains()
    integer                        :: nd !< number of local domains.
    type(local_domain_t), allocatable :: loc_domains(:)
    integer                        :: err, iter, l_start, l_end, l_step, id
    integer                        :: ia, n_spec_def, read_data, iunit, ispec
    integer                        :: length, folder_index
    FLOAT                          :: default_dt, dt
    character(len=MAX_PATH_LEN)    :: filename, folder, folder_default, radiifile
    character(len=MAX_PATH_LEN)    :: in_folder, restart_folder, basename, frmt, ldrestart_folder
    character(len=15)              :: lab
    logical                        :: iterate, ldupdate, ldoverwrite, ldrestart
    type(kick_t)                   :: kick
    type(restart_t)                :: restart, restart_ld

    PUSH_SUB(local_domains)

    message(1) = 'Info: Creating local domains'
    message(2) = ''
    call messages_info(2)

    write(folder_default,'(a)') 'restart/gs/'

    !%Variable LDFolder
    !%Type string
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The folder name where the density used as input file is.
    !%End
    call parse_variable(global_namespace, 'LDFolder', folder_default, folder)

    ! Check if the folder is finished by an /
    if (index(folder, '/', back = .true.) /= len_trim(folder)) then
      folder = trim(folder)//'/'
    end if

    default_dt = M_ZERO
    call parse_variable(global_namespace, 'TDTimeStep', default_dt, dt, unit = units_inp%time)
    if (dt <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
      write(message(2),'(a)') 'Input: TDTimeStep reset to 0. Check input file.'
      call messages_info(2)
    end if

    !%Variable LDFilename
    !%Type string
    !%Default 'density'
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Input filename. The original filename for the density which is going to be 
    !% fragmented into domains.
    !%End
    call parse_variable(global_namespace, 'LDFilename', 'density', basename)
    if ( basename == " " ) basename = ""
    ! Delete the extension if present
    length = len_trim(basename)
    if ( length > 4) then
      if ( basename(length-3:length) == '.obf' ) then
        basename = trim(basename(1:length-4))
      end if
    end if

    !%Variable LDBaderThreshold
    !%Type float
    !%Default 0.01
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable sets the threshold for the basins calculations. Recommended values: 
    !% 0.01 -> intramolecular volumes; 0.2 -> intermolecular volumes.
    !%End
    call parse_variable(global_namespace, 'LDBaderThreshold', CNST(0.01), BaderThreshold)

    !%Variable LDUpdate
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Controls if the calculation of the local domains is desired at each iteration.
    !%End
    call parse_variable(global_namespace, 'LDUpdate', .false., ldupdate)

    !%Variable LDOverWrite
    !%Type logical
    !%Default true
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Controls whether to over-write existing files.
    !%End
    call parse_variable(global_namespace, 'LDOverWrite', .true., ldoverwrite)                       

    !%Variable LDRadiiFile
    !%Type string
    !%Default 'default'
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Full path for the radii file. If set, def_rsize will be reset to the new values. 
    !% This file should have the same format as share/PP/default.
    !%End
    call parse_variable(global_namespace, 'LDRadiiFile', 'default', radiifile)

    if(trim(radiifile) /= "default") then
      n_spec_def = max(0, loct_number_of_lines(radiifile))
      if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment
     
      ! get parameters from file
      do ia = 1, sys%geo%nspecies
        read_data = 0
        iunit = io_open(radiifile, global_namespace, action='read', status='old', die=.false.)
        if(iunit > 0) then
          if(ia == 1) then
            write(message(1),'(a,a)') 'Redefining def_rsize from file:', trim(radiifile)
            call messages_info(1)
          end if
          read(iunit,*)
          default_file: do ispec = 1, n_spec_def
            read(iunit,*) lab
            if (trim(lab) == trim(species_label(sys%geo%species(ia)))) then
              call read_from_default_file(iunit, read_data, sys%geo%species(ia))
              exit default_file
            end if
          end do default_file
          call io_close(iunit)
        else
          write(message(1),'(a,a)') 'Octopus could not open then file:', trim(radiifile)
          call messages_warning(1)
        end if
      end do
    end if

    !TODO: use standart FromSratch and %RestartOptions 
    !%Variable LDRestart
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Restart information will be read from <tt>LDRestartFolder</tt>.
    !%End
    call parse_variable(global_namespace, 'LDRestart', .false., ldrestart)

    if (ldrestart) then
      write(folder_default,'(a)')'ld.general'

      !%Variable LDRestartFolder
      !%Type string
      !%Default "ld.general"
      !%Section Utilities::oct-local_multipoles
      !%Description
      !% The folder name where the density used as input file is.
      !%End
      call parse_variable(global_namespace, 'LDRestartFolder', folder_default, ldrestart_folder)

      ! Check if the folder is finished by an /
      if (index(ldrestart_folder, '/', .true.) /= len_trim(ldrestart_folder)) then
        write(ldrestart_folder,'(a,a1)') trim(ldrestart_folder), '/'
      end if
    end if


    !%Variable LDIterateFolder
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable decides if a folder is going to be iterated.
    !%End
    call parse_variable(global_namespace, 'LDIterateFolder', .false., iterate)

    !%Variable LDStart
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The starting number of the filename or folder.
    !%End
    call parse_variable(global_namespace, 'LDStart', 0, l_start)

    !%Variable LDEnd
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The last number of the filename or folder.
    !%End
    call parse_variable(global_namespace, 'LDEnd', 0, l_end)

    !%Variable LDStep
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The padding between the filenames or folder.
    !%End
    call parse_variable(global_namespace, 'LDStep', 1, l_step)

    message(1) = 'Info: Computing local multipoles'
    message(2) = ''
    call messages_info(2)

    call local_init(sys%space, sys%gr%mesh, sys%geo, nd, loc_domains)

    ! Starting loop over selected densities.
    if (any(loc_domains(:)%dshape == BADER)) then
      call messages_experimental('Bader volumes in oct-local_multipoles')
    end if

    call kick_init(kick, global_namespace, sys%gr%sb, sys%st%d%ispin)
    do id = 1, nd
      call local_write_init(loc_domains(id)%writ, global_namespace, loc_domains(id)%lab, 0, dt)
    end do

    !TODO: initialize hamiltonian if needed: check for LDOuput = energy or potential, using local_write_check_hm(local%writ)

    if (ldrestart) then
      !TODO: check for domains & mesh compatibility 
      call restart_init(restart_ld, global_namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, sys%mc, err, &
                        dir=trim(ldrestart_folder), mesh = sys%gr%mesh)
      call local_restart(sys%gr%mesh, sys%geo, nd, loc_domains, restart_ld)
      call restart_end(restart_ld)
    end if

    ! Initialize the restart directory from <tt>LDFolder</tt> value.
    ! This directory has to have the files 'grid', 'mesh' and 'lxyz.obf'
    ! and the files that are going to be used, must be inside this folder
    if (iterate) then
      ! Delete the last / and find the previous /, if any
      in_folder = folder(1:len_trim(folder)-1)
      folder_index = index(in_folder, '/', .true.)
      restart_folder = in_folder(1:folder_index)
    else 
      restart_folder = folder
    end if
    call restart_init(restart, global_namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, sys%mc, err, &
                      dir=trim(restart_folder), mesh = sys%gr%mesh)

!!$    call loct_progress_bar(-1, l_end-l_start)
    do iter = l_start, l_end, l_step
      if (iterate) then
        ! Delete the last / and add the corresponding folder number
        write(in_folder,'(a,i0.7,a)') folder(folder_index+1:len_trim(folder)-1),iter, "/"
        write(filename, '(a,a,a)') trim(in_folder), trim(basename)
      else
        in_folder = "."
        if ( l_start /= l_end ) then
          ! Here, we are only considering 10 character long filenames.
          ! Subtract the initial part given at 'ConvertFilename' from the format and pad
          ! with zeros.
          write(frmt,'(a,i0,a)')"(a,i0.",10-len_trim(basename),")"
          write(filename, fmt=trim(frmt)) trim(basename), iter
        else 
          ! Assuming filename is given complete in the 'LDFilename'
          write(filename, '(a,a,a,a)') trim(in_folder),"/", trim(basename)
          filename = basename
        end if
      end if
      ! FIXME: there is a special function for reading the density. Why not use that?
      ! TODO: Find the function that reads the density. Which one?
      ! FIXME: why only real functions? Please generalize.
      ! TODO: up to know the local_multipoles utlity acts over density functions, which are real.
      if (err == 0) then
        call drestart_read_mesh_function(restart, trim(filename), sys%gr%mesh, sys%st%rho(:,1), err)
      end if
      if (err /= 0 ) then
        write(message(1),*) 'While reading density: "', trim(filename), '", error code:', err
        call messages_fatal(1)
      end if

      ! Look for the mesh points inside local domains
      if ((iter == l_start .and. .not. ldrestart) .or. ldupdate) then
        call local_inside_domain(sys%gr%mesh, sys%geo, sys%mc, nd, loc_domains, global_namespace, sys%st%rho(:,1))
      end if

      do id = 1, nd
        call local_write_iter(loc_domains(id)%writ, global_namespace, loc_domains(id)%lab, loc_domains(id)%ions_inside, &
          loc_domains(id)%inside, loc_domains(id)%dcm, sys%gr, sys%st, sys%hm, sys%ks, sys%geo, kick, iter, l_start, ldoverwrite)
      end do
      call loct_progress_bar(iter-l_start, l_end-l_start) 
    end do

    call restart_end(restart)
    call kick_end(kick)
    call local_end(nd, loc_domains)

    message(1) = 'Info: Exiting local domains'
    message(2) = ''
    call messages_info(2)

    POP_SUB(local_domains)
  end subroutine local_domains

  ! ---------------------------------------------------------
  !> Initialize local_domain_t variable, allocating variable 
  !! and reading parameters from input file. 
  ! ---------------------------------------------------------
  subroutine local_init(space, mesh, geo, nd, loc_domains)
    type(space_t),    intent(in)  :: space
    type(mesh_t),     intent(in)  :: mesh
    type(geometry_t), intent(in)  :: geo
    integer,          intent(out) :: nd
    type(local_domain_t), allocatable, intent(out) :: loc_domains(:)

    integer           :: id
    type(block_t)     :: blk

    PUSH_SUB(local_init)

    !TODO: use same strategy used in %Species block to define variables
    !%Variable LocalDomains
    !%Type block
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The LocalDomains are by definition part of the global grid. The domains are defined by 
    !% selecting a type shape. The domain box will be constructed using the given parameters. 
    !% A local domain could be construct by addition of several box centered on the ions.
    !% The grid points inside this box will belong to the local domain. 
    !% 
    !% The format of this block is the following:<br>
    !% <tt> 'Label' | Shape | %< | Shape dependencies >% </tt>
    !% <br>The first field is the label of the domain. 
    !% Label = string with the name of the new local domain.
    !% The second is the shape type of the box used to define the domain.
    !% Shape = SPHERE, CYLINDER, PARALLELEPIPED, MINIMUM, BADER.
    !% Some types may need some parameters given in the remaining fields of the row.
    !% (the valid options are detailed below). 
    !%
    !% <tt>%LocalDomains
    !% <br>case(SPHERE):         | rsize | %<dim origin coordinates>
    !% <br>case(CYLINDER):       | rsize | xsize | %<origin coordinates>
    !% <br>case(PARALLELEPIPED): | %<lsize> | %<origin coordinates>
    !% <br>case(MINIMUM):        | rsize | 'center_list' 
    !% <br>case(BADER):          | 'center_list' 
    !% <br>%</tt>
    !% <br>rsize: Radius in input length units
    !% <br>xsize: the length of the cylinder in the x-direction 
    !% <br>origin coordinates: in input length units separated by |, where the box is centered.
    !% <br>lsize: half of the length of the parallelepiped in each direction.
    !% <br>center_list: string containing the list of atoms in xyz file for each domain in the form "2,16-23"
    !%End

    ! First, find out if there is a LocalDomains block.
    nd = 0
    if (parse_block(global_namespace, 'LocalDomains', blk) == 0) then
      nd = parse_block_n(blk)
    end if

    SAFE_ALLOCATE(loc_domains(1:nd))

    block: do id = 1, nd
      SAFE_ALLOCATE(loc_domains(id)%inside(1:mesh%np))
      SAFE_ALLOCATE(loc_domains(id)%ions_inside(1:geo%natoms))
      SAFE_ALLOCATE(loc_domains(id)%dcm(1:geo%space%dim))
      call parse_block_string(blk, id-1, 0, loc_domains(id)%lab)
      call local_read_from_block(space, geo, blk, id-1, loc_domains(id)%domain, loc_domains(id)%dshape, loc_domains(id)%clist, &
        global_namespace)
    end do block
    call parse_block_end(blk)
    message(1) = ''
    call messages_info(1)

    POP_SUB(local_init)
  end subroutine local_init

  ! ---------------------------------------------------------
  !> Ending local_domain_t variable, allocating variable 
  !! and reading parameters from input file. 
  ! ---------------------------------------------------------
  subroutine local_end(nd, loc_domains)
    integer,              intent(in)     :: nd
    type(local_domain_t), allocatable, intent(inout) :: loc_domains(:)

    integer :: id

    PUSH_SUB(local_end)

    do id = 1, nd
      call box_union_end(loc_domains(id)%domain)
      call local_write_end(loc_domains(id)%writ)
      SAFE_DEALLOCATE_A(loc_domains(id)%inside)
      SAFE_DEALLOCATE_A(loc_domains(id)%ions_inside)
      SAFE_DEALLOCATE_A(loc_domains(id)%dcm)
    end do
    SAFE_DEALLOCATE_A(loc_domains)

    POP_SUB(local_end)
  end subroutine local_end

  ! ---------------------------------------------------------
  subroutine local_read_from_block(space, geo, blk, row, dom, shape, clist, namespace)
    type(space_t),     intent(in)        :: space
    type(geometry_t),  intent(in)        :: geo
    type(block_t),     intent(in)        :: blk
    integer,           intent(in)        :: row
    type(box_union_t), intent(inout)     :: dom
    integer,           intent(out)       :: shape
    character(len=*),  intent(out)       :: clist
    type(namespace_t), intent(in)        :: namespace
    
    integer                     :: ic, nb
    FLOAT                       :: rsize, xsize
    FLOAT                       :: center(MAX_DIM), lsize(MAX_DIM)
   
    PUSH_SUB(local_read_from_block)

    ! Initializing variables in dom
    nb = 1
    rsize = -M_ONE
    xsize = M_ZERO
    lsize(:) = M_ZERO
    center(:) = M_ZERO

    call parse_block_integer(blk, row, 1, shape)

    select case (shape)
    case (MINIMUM)
      if(geo%reduced_coordinates) then
        message(1) = "The 'minimum' box shape cannot be used if atomic positions"
        message(2) = "are given as reduced coordinates."
        call messages_fatal(2)
      end if
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if (rsize < M_ZERO) call messages_input_error(namespace, 'radius', row=row, column=2)
      call parse_block_string(blk, row, 3, clist)
      nb = 0
      do ic = 1, geo%natoms
        if(loct_isinstringlist(ic, clist)) nb = nb + 1
      end do

    case (SPHERE)
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if(rsize < M_ZERO) call messages_input_error(namespace, 'radius', row=row, column=2)
      do ic = 1, geo%space%dim 
        call parse_block_float(blk, row, 2 + ic, center(ic), unit = units_inp%length)
      end do
      lsize(1:space%dim) = rsize

    case (CYLINDER)
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if(rsize < M_ZERO) call messages_input_error(namespace, 'radius', row=row, column=2)
      call parse_block_float(blk, row, 3, xsize, unit = units_inp%length)
      do ic = 1, geo%space%dim 
        call parse_block_float(blk, row, 3 + ic, center(ic), unit = units_inp%length)
      end do
      lsize(1)     = xsize
      lsize(2:space%dim) = rsize

    case (PARALLELEPIPED)
      do ic = 1, geo%space%dim
        call parse_block_float(blk, row, 1 + ic, lsize(ic), unit = units_inp%length)
      end do
      do ic = 1, geo%space%dim
        call parse_block_float(blk, row, 1 + space%dim + ic, center(ic), unit = units_inp%length)
      end do

    case (BADER)
      ! FIXME: when input error exists --> segmentation fault appears
      call parse_block_string(blk, row, 2, clist)
      nb = 0
      do ic = 1, geo%natoms
        if(loct_isinstringlist(ic, clist)) nb = nb + 1
      end do
    end select

    call local_domains_init(geo, dom, space%dim, shape, center, rsize, lsize, nb, clist, namespace)

    POP_SUB(local_read_from_block)

  end subroutine local_read_from_block

  !!---------------------------------------------------------------------------^
  subroutine local_domains_init(geo, dom, dim, shape, center, rsize, lsize, nb, clist, namespace)
    type(geometry_t),  intent(in)    :: geo
    type(box_union_t), intent(inout) :: dom
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: shape
    FLOAT,             intent(in)    :: center(dim)
    FLOAT,             intent(in)    :: rsize
    FLOAT,             intent(in)    :: lsize(MAX_DIM)
    integer,           intent(in)    :: nb
    character(len=*),  intent(in)    :: clist
    type(namespace_t), intent(in)    :: namespace

    integer                  :: ia, ibox, ic, bshape
    FLOAT                    :: bcenter(MAX_DIM), bsize(MAX_DIM)
    type(box_t), allocatable :: boxes(:)

    PUSH_SUB(local_domains_init)

    SAFE_ALLOCATE(boxes(1:nb))
    bsize(:) = M_ZERO
    ibox = 1

    select case (shape)
    case (SPHERE, CYLINDER)
      call box_create(boxes(ibox), shape, dim, lsize, center, namespace)
    case (PARALLELEPIPED)
      bshape         = 3
      call box_create(boxes(ibox), bshape, dim, lsize, center, namespace)
    case (MINIMUM) 
      bshape         = SPHERE
      do ia = 1, geo%natoms
        if(loct_isinstringlist(ia, clist))then
          bcenter(1:geo%space%dim) = geo%atom(ia)%x(1:geo%space%dim)
          bsize(:) = rsize
          if (bsize(1) < M_EPSILON) bsize(1) = species_def_rsize(geo%atom(ia)%species)
          call box_create(boxes(ibox), bshape, dim, bsize, bcenter, namespace)
          ibox = ibox + 1
        end if
      end do
    case (BADER) 
      bshape         = SPHERE
      do ia = 1, geo%natoms
        if(loct_isinstringlist(ia, clist))then
          bcenter(1:geo%space%dim) = geo%atom(ia)%x(1:geo%space%dim)
          bsize(:) = species_def_h(geo%atom(ia)%species)
          call box_create(boxes(ibox), bshape, dim, bsize, bcenter, namespace)
          ibox = ibox + 1
        end if
      end do
    end select

    call box_union_init(dom, nb, boxes)

    ! TODO: Check for a conflict between box_union and clist for Bader Volumes
    ic = 0
    if (shape == MINIMUM) then
      do ia = 1, geo%natoms
        if (box_union_inside(dom, geo%atom(ia)%x) .and. .not.loct_isinstringlist(ia, clist) ) then
          ic = ic + 1
          if(ic <= 20) write(message(ic),'(a,a,I0,a,a)')'Atom: ',trim(species_label(geo%atom(ia)%species)),ia, &
                                 ' is inside the union box BUT not in list: ',trim(clist)
        end if
      end do
      if (ic > 0) then 
        if (ic > 20) ic = 20
        call messages_warning(ic)
        message(1) = ' THIS COULD GIVE INCORRECT RESULTS '
        call messages_warning(1)
        if (ic == 20) then
          message(1) = ' AT LEAST 19 ATOMS ARE NOT PRESENT IN THE LIST'
          call messages_warning(1)
        end if
      end if
    end if

    do ibox = 1, nb 
      call box_end(boxes(ibox))
    end do
    SAFE_DEALLOCATE_A(boxes)

    POP_SUB(local_domains_init)
  end subroutine local_domains_init

  ! ---------------------------------------------------------
  subroutine local_restart(mesh, geo, nd, loc_domains, restart)
    type(mesh_t),         intent(in)    :: mesh
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: nd
    type(local_domain_t), intent(inout) :: loc_domains(:)
    type(restart_t),      intent(in)    :: restart

    integer                     :: id, ip, iunit, ierr
    character(len=MAX_PATH_LEN) :: filename, tmp
    character(len=31) :: line(1)
    FLOAT, allocatable          :: inside(:)

    PUSH_SUB(local_restart)

    message(1) = 'Info: Reading mesh points inside each local domain'
    call messages_info(1)

    SAFE_ALLOCATE(inside(1:mesh%np))
    !Read local domain information from ldomains.info 
    filename = "ldomains.info"    
    iunit = restart_open(restart, filename, status='old')
    call restart_read(restart, iunit, line, 1, ierr)    
    read(line(1),'(a25,1x,i5)') tmp, ierr
    call restart_close(restart, iunit)

    filename = "ldomains"
    call drestart_read_mesh_function(restart, trim(filename), mesh, inside, ierr) 

    do id = 1, nd
      loc_domains(id)%inside = .false.
      do ip = 1 , mesh%np
        if (bitand(int(inside(ip)), 2**id) /= 0) loc_domains(id)%inside(ip) = .true.
      end do

      !Check for atom list inside each domain
      call local_ions_inside(loc_domains(id)%inside, geo, mesh, loc_domains(id)%ions_inside, loc_domains(id)%clist)

      !Compute center of mass of each domain
      call local_center_of_mass(loc_domains(id)%ions_inside, geo, loc_domains(id)%dcm)
    end do

    SAFE_DEALLOCATE_A(inside)

    POP_SUB(local_restart)
  end subroutine local_restart

  ! ---------------------------------------------------------
  subroutine local_inside_domain(mesh, geo, mc, nd, loc_domains, namespace, ff)
    type(mesh_t),           intent(in)    :: mesh
    type(geometry_t),       intent(in)    :: geo
    type(multicomm_t),      intent(in)    :: mc
    integer,                intent(in)    :: nd
    type(local_domain_t),   intent(inout) :: loc_domains(:)
    type(namespace_t),      intent(in)    :: namespace
    FLOAT,                  intent(in)    :: ff(:)
    
    integer(8)          :: how
    integer             :: id, ip, iunit, ierr
    type(basins_t)      :: basins
    FLOAT, allocatable  :: ff2(:,:)
    character(len=MAX_PATH_LEN)   :: filename, base_folder, folder, frmt
    character(len=140), allocatable  :: lines(:)
    type(restart_t)     :: restart
    FLOAT, allocatable :: dble_domain_map(:), domain_mesh(:)

    PUSH_SUB(local_inside_domain)

    !TODO: find a parallel algorithm to perform the charge-density fragmentation.
    message(1) = 'Info: Assigning mesh points inside each local domain'
    call messages_info(1)

    if (debug%info) then
      call parse_variable(global_namespace, 'LDOutputFormat', 0, how)
      if (.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
        call messages_input_error(global_namespace, 'LDOutputFormat')
      end if
    end if

    SAFE_ALLOCATE(ff2(1:mesh%np,1))
    ff2(1:mesh%np,1) = ff(1:mesh%np)

    ! First we do the Bader domains
    if (any(loc_domains(:)%dshape == BADER)) then
      if (mesh%parallel_in_domains) then
        write(message(1),'(a)') 'Bader volumes can only be computed in serial'
        call messages_fatal(1)
      end if

      call add_dens_to_ion_x(ff2, namespace, mesh, geo)
      call basins_init(basins, mesh)
      call parse_variable(namespace, 'LDBaderThreshold', CNST(0.01), BaderThreshold)
      call basins_analyze(basins, mesh, ff2(:,1), ff2, BaderThreshold)
      call bader_union_inside(basins, mesh, geo, nd, loc_domains)

      if (debug%info) then
        filename = 'basinsmap'
        call dio_function_output(how, trim('local.general'), trim(filename), global_namespace, mesh, &
          DBLE(basins%map(1:mesh%np)), unit_one, ierr, geo = geo)
        call dio_function_output(how, trim('local.general'), 'dens_ff2', global_namespace, mesh, ff2(:,1), unit_one, ierr, &
          geo = geo)
        call io_close(iunit)
      end if
      call basins_end(basins)
    end if

    ! Then all the remaining domains
    do id = 1, nd
      if (loc_domains(id)%dshape == BADER) cycle
      call box_union_inside_vec(loc_domains(id)%domain, mesh%np, mesh%x, loc_domains(id)%inside)
    end do

    if (debug%info) then
      ! Write extra information about domains
      SAFE_ALLOCATE(dble_domain_map(1:mesh%np))
      SAFE_ALLOCATE(domain_mesh(1:mesh%np))

      domain_mesh = M_ZERO
      do id = 1, nd
        dble_domain_map = M_ZERO
        do ip = 1, mesh%np
          if (loc_domains(id)%inside(ip)) then
            dble_domain_map(ip) = TOFLOAT(id)
            domain_mesh(ip) = domain_mesh(ip) + dble_domain_map(ip)
          end if
        end do

        write(filename,'(a,a,a)') 'domain.', trim(loc_domains(id)%lab)
        call dio_function_output(how, trim('local.general'), trim(filename), global_namespace, mesh, &
          dble_domain_map(1:mesh%np), unit_one, ierr, geo = geo)
      end do

      write(filename,'(a,a,a)') 'domain.mesh'
      call dio_function_output(how, trim('local.general'), trim(filename), global_namespace, mesh, domain_mesh(1:mesh%np), &
        unit_one, ierr, geo = geo)

      SAFE_DEALLOCATE_A(dble_domain_map)
      SAFE_DEALLOCATE_A(domain_mesh)
    end if

    do id = 1, nd
      !Check for atom list inside each domain
      call local_ions_inside(loc_domains(id)%inside, geo, mesh, loc_domains(id)%ions_inside, loc_domains(id)%clist)

      !Compute center of mass of each domain
      call local_center_of_mass(loc_domains(id)%ions_inside, geo, loc_domains(id)%dcm)
    end do

    !Write restart file for local domains
    base_folder = "./restart/"
    folder = "ld/"
    filename = "ldomains"
    write(message(1),'(a,a)')'Info: Writing restart info to ', trim(filename)
    call messages_info(1)
    call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_DUMP, mc, ierr, mesh=mesh, &
      dir=trim(base_folder)//trim(folder))
    ff2 = M_ZERO
    SAFE_ALLOCATE(lines(1:nd+2))
    write(lines(1),'(a,1x,i5)') 'Number of local domains =', nd
    write(lines(2),'(a3,1x,a15,1x,a5,1x,71x,a9,1x,a14)') '#id', 'label', 'shape', 'Atom list', 'center of mass'
    do id = 1, nd
      write(frmt,'(a,i0,a)') '(i3,1x,a15,1x,i5,1x,a80,1x', geo%space%dim, '(f10.6,1x))'
      write(lines(id+2), fmt=trim(frmt))id, trim(loc_domains(id)%lab), loc_domains(id)%dshape, trim(loc_domains(id)%clist), &
        loc_domains(id)%dcm(1:geo%space%dim)
      do ip = 1, mesh%np
        if (loc_domains(id)%inside(ip)) ff2(ip, 1) = ff2(ip, 1) + 2**DBLE(id)
      end do
    end do
    call drestart_write_mesh_function(restart, filename, mesh, ff2(1:mesh%np, 1), ierr)

    filename = "ldomains.info"
    iunit = restart_open(restart, filename)
    call restart_write(restart, iunit, lines, nd+2, ierr)
    call restart_end(restart)
    call io_close(iunit)
    SAFE_DEALLOCATE_A(lines)
    
    SAFE_DEALLOCATE_A(ff2)
    
    POP_SUB(local_inside_domain)
  end subroutine local_inside_domain

  ! ---------------------------------------------------------
  subroutine bader_union_inside(basins, mesh, geo, nd, loc_domains)
    type(basins_t),       intent(inout) :: basins
    type(mesh_t),         intent(in)    :: mesh
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: nd 
    type(local_domain_t), intent(inout) :: loc_domains(:)

    integer               :: ia, ib, id, ip, ix, rankmin, n_boxes
    integer, allocatable  :: domain_map(:), ion_map(:)
    FLOAT                 :: dmin, dd
    FLOAT, allocatable    :: xi(:)
    logical               :: lduseatomicradii

    PUSH_SUB(bader_union_inside)

    SAFE_ALLOCATE(xi(1:geo%space%dim))

    !%Variable LDUseAtomicRadii
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% If set, atomic radii will be used to assign lone pairs to ion. 
    !%End
    call parse_variable(global_namespace, 'LDUseAtomicRadii', .false., lduseatomicradii)

    SAFE_ALLOCATE(ion_map(geo%natoms))

    do id = 1, nd
      if (loc_domains(id)%dshape /= BADER) cycle

      n_boxes = box_union_get_nboxes(loc_domains(id)%domain)
      SAFE_ALLOCATE(domain_map(1:n_boxes))

      ! Assign basins%map to ions
      if (lduseatomicradii .and. (id == 1)) then
        ion_map = 0
        do ia = 1, geo%natoms
          ion_map(ia) = basins%map(mesh_nearest_point(mesh, geo%atom(ia)%x, dmin, rankmin))
        end do
        ! Assign lonely pair to ion in a atomic radii distance
        do ip = 1, mesh%np
          if (all(ion_map(:) /= basins%map(ip)) ) then
            do ia = 1, geo%natoms
              dd = sum((mesh%x(ip, 1:geo%space%dim) - geo%atom(ia)%x(1:geo%space%dim))**2)
              dd = sqrt(dd)
              if (dd <= species_vdw_radius(geo%atom(ia)%species)) basins%map(ip) = ion_map(ia)
            end do
          end if
        end do
      end if
      domain_map = 0
      do ib = 1, n_boxes
        xi = box_union_get_center(loc_domains(id)%domain, ib)
        ix = mesh_nearest_point(mesh, xi, dmin, rankmin)
        domain_map(ib) = basins%map(ix)
      end do
      do ip = 1, mesh%np
        loc_domains(id)%inside(ip) = any(domain_map(1:n_boxes) == basins%map(ip))
      end do

      SAFE_DEALLOCATE_A(domain_map)
    end do

    SAFE_DEALLOCATE_A(ion_map)
    SAFE_DEALLOCATE_A(xi)

    POP_SUB(bader_union_inside)
  end subroutine bader_union_inside

  ! ---------------------------------------------------------
  subroutine add_dens_to_ion_x(ff, namespace, mesh, geo)
    FLOAT,              intent(inout) :: ff(:,:)
    type(namespace_t),  intent(in)    :: namespace
    type(mesh_t),       intent(in)    :: mesh
    type(geometry_t),   intent(in)    :: geo

    integer :: ia, is
    FLOAT, allocatable :: ffs(:)

    PUSH_SUB(add_dens_to_ion_x)

    SAFE_ALLOCATE(ffs(1:mesh%np))
    do ia = 1, geo%natoms
      call species_get_long_range_density(geo%atom(ia)%species, namespace, geo%atom(ia)%x, mesh, ffs)
      do is = 1, ubound(ff, dim=2)
        ff(1:mesh%np, is) = ff(1:mesh%np, is) - ffs(1:mesh%np)
      end do
    end do
    SAFE_DEALLOCATE_A(ffs)

    POP_SUB(add_dens_to_ion_x)
  end subroutine add_dens_to_ion_x

  ! ---------------------------------------------------------
  !Check for the ions inside each local domain.
  subroutine local_ions_inside(inside, geo, mesh, ions_inside, ions_list)
    logical,            intent(in)  :: inside(:)
    type(geometry_t),   intent(in)  :: geo
    type(mesh_t),       intent(in)  :: mesh
    logical,            intent(out) :: ions_inside(:)
    character(len=500), intent(out) :: ions_list

    integer              :: ia, ix, ic
    integer, allocatable :: inside_tmp(:)
    integer              :: rankmin
    FLOAT                :: dmin
    character(len=500)   :: chtmp

    PUSH_SUB(local_ions_inside)

    SAFE_ALLOCATE(inside_tmp(geo%natoms))
    inside_tmp = 0
   
    do ia = 1, geo%natoms
      ix = mesh_nearest_point(mesh, geo%atom(ia)%x, dmin, rankmin )
      if (rankmin /= mesh%mpi_grp%rank) cycle
      if (inside(ix)) inside_tmp(ia) = 1
    end do

    if(mesh%parallel_in_domains) then
      call comm_allreduce(mesh%mpi_grp%comm, inside_tmp, geo%natoms)
    end if                               
    ions_inside = inside_tmp == 1

    SAFE_DEALLOCATE_A(inside_tmp)

    !print list of atoms
    ic = 0
    ix = 0
    ions_list = ""
    chtmp = ""
    do ia = 1, geo%natoms-1
      if (ions_inside(ia)) then
        if (ic == 0 .or. .not.ions_inside(ia+1)) then
          write(ions_list, '(a,i0)') trim(chtmp), ia
        end if
        if (ions_inside(ia + 1)) then 
          ic = ic + 1
          chtmp = trim(ions_list)//"-"
        else
          ic = 0
          chtmp = trim(ions_list)//","
        end if

      else
        cycle 
      end if
    end do

    if (ions_inside(geo%natoms)) then
      write(ions_list, '(a,i0)') trim(chtmp), ia
    end if

    POP_SUB(local_ions_inside)
  end subroutine local_ions_inside

  ! ---------------------------------------------------------
  subroutine local_center_of_mass(ions_inside, geo, center)
    logical,           intent(in)  :: ions_inside(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(out) :: center(:)

    integer :: ia
    FLOAT :: sumw

    PUSH_SUB(local_center_of_mass)

    sumw = M_ZERO
    center = M_ZERO
    do ia = 1, geo%natoms
      if (ions_inside(ia)) then
        center(1:geo%space%dim) = center(1:geo%space%dim) + geo%atom(ia)%x(1:geo%space%dim)*species_mass(geo%atom(ia)%species)
        sumw = sumw + species_mass(geo%atom(ia)%species)
      end if
    end do
    center(1:geo%space%dim) = center(1:geo%space%dim)/sumw

    POP_SUB(local_center_of_mass)
  end subroutine local_center_of_mass

end program oct_local_multipoles

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
