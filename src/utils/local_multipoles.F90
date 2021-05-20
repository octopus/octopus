!! Copyright (C) 2014 M. Oliveira, J. Jornet-Somoza
!! Copyright (C) 2020-2021 M. Oliveira
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
  use box_cylinder_oct_m
  use box_oct_m
  use box_parallelepiped_oct_m
  use box_sphere_oct_m
  use box_union_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
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
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none
  
  type local_domain_t
    class(box_t), pointer :: box            !< box defining the domain
    character(len=500)    :: clist          !< list of centers
    character(len=15)     :: lab            !< declared name
    integer               :: dshape         !< shape of domain
    logical,  allocatable :: mesh_mask(:)   !< mesh points inside the domain
    logical,  allocatable :: ions_mask(:)   !< ions inside the domain
    type(local_write_t)   :: writ           !< write handler for local domains analysis
  end type local_domain_t

  type(electrons_t), pointer :: sys
  integer, parameter         :: BADER = 9  ! should be = OPTION__OUTPUT__BADER

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
    if (dt < M_ZERO) then
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
    if (basename == " ") basename = ""
    ! Delete the extension if present
    length = len_trim(basename)
    if (length > 4) then
      if (basename(length-3:length) == '.obf') then
        basename = trim(basename(1:length-4))
      end if
    end if

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
      do ia = 1, sys%ions%nspecies
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
            if (trim(lab) == trim(species_label(sys%ions%species(ia)))) then
              call read_from_default_file(iunit, read_data, sys%ions%species(ia))
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

    call local_init(sys%space, sys%gr%mesh, sys%ions, nd, loc_domains)

    ! Starting loop over selected densities.
    if (any(loc_domains(:)%dshape == BADER)) then
      call messages_experimental('Bader volumes in oct-local_multipoles')
    end if

    call kick_init(kick, global_namespace, sys%space, sys%kpoints, sys%st%d%ispin)
    do id = 1, nd
      call local_write_init(loc_domains(id)%writ, global_namespace, loc_domains(id)%lab, 0, dt)
    end do

    !TODO: initialize hamiltonian if needed: check for LDOuput = energy or potential, using local_write_check_hm(local%writ)

    if (ldrestart) then
      !TODO: check for domains & mesh compatibility 
      call restart_init(restart_ld, global_namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, sys%mc, err, &
                        dir=trim(ldrestart_folder), mesh = sys%gr%mesh)
      call local_restart_read(sys%space, sys%gr%mesh, sys%ions, nd, loc_domains, restart_ld)
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
        call drestart_read_mesh_function(restart, sys%space, trim(filename), sys%gr%mesh, sys%st%rho(:,1), err)
      end if
      if (err /= 0 ) then
        write(message(1),*) 'While reading density: "', trim(filename), '", error code:', err
        call messages_fatal(1)
      end if

      ! Look for the mesh points inside local domains
      if ((iter == l_start .and. .not. ldrestart) .or. ldupdate) then
        call local_inside_domain(sys%space, sys%gr%mesh, sys%ions, nd, loc_domains, global_namespace, sys%st%rho(:,1))
        call local_restart_write(global_namespace, sys%space, sys%gr%mesh, sys%mc, nd, loc_domains)
      end if

      do id = 1, nd
        call local_write_iter(loc_domains(id)%writ, global_namespace, sys%space, loc_domains(id)%lab, loc_domains(id)%ions_mask, &
          loc_domains(id)%mesh_mask, sys%gr%mesh, sys%st, sys%hm, sys%ks, sys%ions, kick, iter, l_start, ldoverwrite)
      end do
      call loct_progress_bar(iter-l_start, l_end-l_start) 
    end do

    call restart_end(restart)
    call kick_end(kick)

    do id = 1, nd
      call local_end(loc_domains(id))
    end do
    SAFE_DEALLOCATE_A(loc_domains)

    message(1) = 'Info: Exiting local domains'
    message(2) = ''
    call messages_info(2)

    POP_SUB(local_domains)
  end subroutine local_domains

  ! ---------------------------------------------------------
  !> Initialize local_domain_t variable, allocating variable 
  !! and reading parameters from input file. 
  ! ---------------------------------------------------------
  subroutine local_init(space, mesh, ions, nd, loc_domains)
    type(space_t),    intent(in)  :: space
    type(mesh_t),     intent(in)  :: mesh
    type(ions_t),     intent(in)  :: ions
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
      SAFE_ALLOCATE(loc_domains(id)%mesh_mask(1:mesh%np))
      SAFE_ALLOCATE(loc_domains(id)%ions_mask(1:ions%natoms))
      call parse_block_string(blk, id-1, 0, loc_domains(id)%lab)
      call local_read_from_block(loc_domains(id), space, ions, blk, id-1, global_namespace)
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
  subroutine local_end(domain)
    type(local_domain_t), intent(inout) :: domain

    class(box_t), pointer :: box

    PUSH_SUB(local_end)

    if (domain%dshape /= BADER) then
      box => domain%box
      SAFE_DEALLOCATE_P(box)
    end if
    call local_write_end(domain%writ)
    SAFE_DEALLOCATE_A(domain%mesh_mask)
    SAFE_DEALLOCATE_A(domain%ions_mask)

    POP_SUB(local_end)
  end subroutine local_end

  ! ---------------------------------------------------------
  subroutine local_read_from_block(domain, space, ions, blk, row, namespace)
    type(local_domain_t), intent(inout) :: domain 
    type(space_t),        intent(in)    :: space
    type(ions_t),         intent(in)    :: ions
    type(block_t),        intent(in)    :: blk
    integer,              intent(in)    :: row
    type(namespace_t),    intent(in)    :: namespace

    integer :: ic, ia
    FLOAT   :: radius, center(space%dim), length(space%dim)
    class(box_union_t), pointer :: minimum_box

    PUSH_SUB(local_read_from_block)

    call parse_block_integer(blk, row, 1, domain%dshape)

    select case (domain%dshape)
    case (MINIMUM)
      call parse_block_float(blk, row, 2, radius, unit = units_inp%length)
      if (radius < M_ZERO) call messages_input_error(namespace, 'radius', row=row, column=2)
      call parse_block_string(blk, row, 3, domain%clist)

      minimum_box => box_union_t(space%dim)
      do ia = 1, ions%natoms
        if (loct_isinstringlist(ia, domain%clist)) then
          if (radius < M_EPSILON) then
            radius = species_def_rsize(ions%atom(ia)%species)
          end if
          call minimum_box%add_box(box_sphere_t(space%dim, ions%pos(:, ia), radius))
        end if
      end do
      domain%box => minimum_box

      ic = 0
      do ia = 1, ions%natoms
        if (domain%box%contains_point(ions%pos(:, ia)) .and. .not. loct_isinstringlist(ia, domain%clist)) then
          ic = ic + 1
          if(ic <= 20) write(message(ic),'(a,a,I0,a,a)')'Atom: ',trim(species_label(ions%atom(ia)%species)), ia, &
            ' is inside the union box BUT not in list: ', trim(domain%clist)
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

    case (SPHERE)
      call parse_block_float(blk, row, 2, radius, unit = units_inp%length)
      if (radius < M_ZERO) call messages_input_error(namespace, 'radius', row=row, column=2)
      do ic = 1, space%dim
        call parse_block_float(blk, row, 2 + ic, center(ic), unit = units_inp%length)
      end do

      domain%box => box_sphere_t(space%dim, center, radius)

    case (CYLINDER)
      call parse_block_float(blk, row, 2, radius, unit = units_inp%length)
      if (radius < M_ZERO) then
        call messages_input_error(namespace, 'radius', row=row, column=2)
      end if
      call parse_block_float(blk, row, 3, length(1), unit = units_inp%length)
      do ic = 1, space%dim 
        call parse_block_float(blk, row, 3 + ic, center(ic), unit = units_inp%length)
      end do

      domain%box => box_cylinder_t(space%dim, center, radius, 1, length(1)*M_TWO, namespace)

    case (PARALLELEPIPED)
      do ic = 1, space%dim
        call parse_block_float(blk, row, 1 + ic, length(ic), unit = units_inp%length)
      end do
      do ic = 1, space%dim
        call parse_block_float(blk, row, 1 + space%dim + ic, center(ic), unit = units_inp%length)
      end do

      domain%box => box_parallelepiped_t(space%dim, center, length*M_TWO)

    case (BADER)
      call parse_block_string(blk, row, 2, domain%clist)
      POP_SUB(local_read_from_block)
      return
    end select

    POP_SUB(local_read_from_block)
  end subroutine local_read_from_block

  ! ---------------------------------------------------------
  !> Write restart files for local domains
  subroutine local_restart_write(namespace, space, mesh, mc, nd, loc_domains)
    type(namespace_t),    intent(in) :: namespace
    type(space_t),        intent(in) :: space
    type(mesh_t),         intent(in) :: mesh
    type(multicomm_t),    intent(in) :: mc
    integer,              intent(in) :: nd
    type(local_domain_t), intent(in) :: loc_domains(:)

    integer :: id, ip, iunit, ierr
    character(len=MAX_PATH_LEN) :: dirname
    character(len=140), allocatable :: lines(:)
    FLOAT, allocatable  :: ff(:,:)
    type(restart_t) :: restart

    PUSH_SUB(local_restart_write)

    dirname = "./restart/ld"
    write(message(1),'(a,a)')'Info: Writing restart info to ', trim(dirname)
    call messages_info(1)

    call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_DUMP, mc, ierr, mesh=mesh, dir=trim(dirname))

    SAFE_ALLOCATE(lines(1:nd+2))
    write(lines(1),'(a,1x,i5)') 'Number of local domains =', nd
    write(lines(2),'(a3,1x,a15,1x,a5,1x)') '#id', 'label', 'shape'
    do id = 1, nd
      write(lines(id+2), '(i3,1x,a15,1x,i5)') id, trim(loc_domains(id)%lab), loc_domains(id)%dshape
    end do
    iunit = restart_open(restart, "ldomains.info")
    call restart_write(restart, iunit, lines, nd+2, ierr)
    call io_close(iunit)
    SAFE_DEALLOCATE_A(lines)

    SAFE_ALLOCATE(ff(1:mesh%np, 1))
    ff = M_ZERO
    do id = 1, nd
      do ip = 1, mesh%np
        if (loc_domains(id)%mesh_mask(ip)) ff(ip, 1) = ff(ip, 1) + 2**TOFLOAT(id)
      end do
    end do
    call drestart_write_mesh_function(restart, space, "ldomains", mesh, ff(1:mesh%np, 1), ierr)
    SAFE_DEALLOCATE_A(ff)

    call restart_end(restart)

    POP_SUB(local_restart_write)
  end subroutine local_restart_write

  ! ---------------------------------------------------------
  subroutine local_restart_read(space, mesh, ions, nd, loc_domains, restart)
    type(space_t),        intent(in)    :: space
    type(mesh_t),         intent(in)    :: mesh
    type(ions_t),         intent(in)    :: ions
    integer,              intent(in)    :: nd
    type(local_domain_t), intent(inout) :: loc_domains(:)
    type(restart_t),      intent(in)    :: restart

    integer            :: id, ip, iunit, ierr
    character(len=31)  :: line(1), tmp
    FLOAT, allocatable :: mask(:)

    PUSH_SUB(local_restart_read)

    message(1) = 'Info: Reading mesh points inside each local domain'
    call messages_info(1)

    SAFE_ALLOCATE(mask(1:mesh%np))

    !Read local domain information from ldomains.info 
    iunit = restart_open(restart, "ldomains.info", status='old')
    call restart_read(restart, iunit, line, 1, ierr)    
    read(line(1),'(a25,1x,i5)') tmp, ierr
    call restart_close(restart, iunit)

    call drestart_read_mesh_function(restart, space, "ldomains", mesh, mask, ierr) 

    do id = 1, nd
      loc_domains(id)%mesh_mask = .false.
      do ip = 1 , mesh%np
        if (bitand(int(mask(ip)), 2**id) /= 0) loc_domains(id)%mesh_mask(ip) = .true.
      end do

      !Check for atom list inside each domain
      call local_ions_mask(loc_domains(id)%mesh_mask, ions, mesh, loc_domains(id)%ions_mask)
    end do

    SAFE_DEALLOCATE_A(mask)

    POP_SUB(local_restart_read)
  end subroutine local_restart_read

  ! ---------------------------------------------------------
  subroutine local_inside_domain(space, mesh, ions, nd, loc_domains, namespace, ff)
    type(space_t),          intent(in)    :: space
    type(mesh_t),           intent(in)    :: mesh
    type(ions_t),           intent(in)    :: ions
    integer,                intent(in)    :: nd
    type(local_domain_t),   intent(inout) :: loc_domains(:)
    type(namespace_t),      intent(in)    :: namespace
    FLOAT,                  intent(in)    :: ff(:)

    integer(8)          :: how
    integer             :: id, ip, ierr
    type(basins_t)      :: basins
    character(len=MAX_PATH_LEN) :: filename
    FLOAT, allocatable :: dble_domain_map(:), domain_mesh(:)

    PUSH_SUB(local_inside_domain)

    !TODO: find a parallel algorithm to perform the charge-density fragmentation.
    message(1) = 'Info: Assigning mesh points inside each local domain'
    call messages_info(1)

    ! If we have Bader domains, then we need to get the basins
    if (any(loc_domains(:)%dshape == BADER)) then
      if (mesh%parallel_in_domains) then
        write(message(1),'(a)') 'Bader volumes can only be computed in serial'
        call messages_fatal(1)
      end if
      call create_basins(namespace, space, mesh, ions, ff, basins)
    end if

    do id = 1, nd
      ! Create the mask that tells which mesh points are inside the domain
      select case (loc_domains(id)%dshape)
      case (BADER)
        call bader_domain_create_mask(loc_domains(id), basins, mesh, ions)
      case default
        call box_domain_create_mask(loc_domains(id), mesh)
      end select

      !Check for atom list inside each domain
      call local_ions_mask(loc_domains(id)%mesh_mask, ions, mesh, loc_domains(id)%ions_mask)
    end do

    if (any(loc_domains(:)%dshape == BADER)) then
      call basins_end(basins)
    end if

    if (debug%info) then
      call parse_variable(namespace, 'LDOutputFormat', 0, how)
      if (.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
        call messages_input_error(namespace, 'LDOutputFormat')
      end if

      ! Write extra information about domains
      SAFE_ALLOCATE(dble_domain_map(1:mesh%np))
      SAFE_ALLOCATE(domain_mesh(1:mesh%np))

      domain_mesh = M_ZERO
      do id = 1, nd
        dble_domain_map = M_ZERO
        do ip = 1, mesh%np
          if (loc_domains(id)%mesh_mask(ip)) then
            dble_domain_map(ip) = TOFLOAT(id)
            domain_mesh(ip) = domain_mesh(ip) + dble_domain_map(ip)
          end if
        end do

        write(filename,'(a,a)') 'domain.', trim(loc_domains(id)%lab)
        call dio_function_output(how, 'local.general', trim(filename), namespace, space, mesh, dble_domain_map, unit_one, ierr, &
          ions = ions)
      end do

      call dio_function_output(how, 'local.general', 'domain.mesh', namespace, space, mesh, domain_mesh, unit_one, ierr, &
        ions = ions)

      SAFE_DEALLOCATE_A(dble_domain_map)
      SAFE_DEALLOCATE_A(domain_mesh)
    end if
    
    POP_SUB(local_inside_domain)
  end subroutine local_inside_domain

  ! ---------------------------------------------------------
  subroutine create_basins(namespace, space, mesh, ions, ff, basins)
    type(namespace_t),  intent(in)    :: namespace
    type(space_t),      intent(in)    :: space
    type(mesh_t),       intent(in)    :: mesh
    type(ions_t),       intent(in)    :: ions
    FLOAT,              intent(in)    :: ff(:)
    type(basins_t),     intent(inout) :: basins

    logical               :: use_atomic_radii
    integer               :: ip, ierr, ia, is, rankmin
    integer(8)            :: how
    FLOAT                 :: bader_threshold, dmin, dd
    integer, allocatable  :: ion_map(:)
    FLOAT,   allocatable  :: ff2(:,:), ffs(:)

    PUSH_SUB(create_basins)

    !%Variable LDBaderThreshold
    !%Type float
    !%Default 0.01
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable sets the threshold for the basins calculations. Recommended values:
    !% 0.01 -> intramolecular volumes; 0.2 -> intermolecular volumes.
    !%End
    call parse_variable(namespace, 'LDBaderThreshold', CNST(0.01), bader_threshold)

    ! Copy function used to construct the basins, as we need to modify it
    SAFE_ALLOCATE(ff2(1:mesh%np,1))
    ff2(1:mesh%np,1) = ff(1:mesh%np)

    ! Add long range density from atoms
    SAFE_ALLOCATE(ffs(1:mesh%np))
    do ia = 1, ions%natoms
      call species_get_long_range_density(ions%atom(ia)%species, ions%namespace, ions%space, ions%latt, &
        ions%pos(:, ia), mesh, ffs)
      do is = 1, ubound(ff2, dim=2)
        ff2(1:mesh%np, is) = ff2(1:mesh%np, is) - ffs(1:mesh%np)
      end do
    end do
    SAFE_DEALLOCATE_A(ffs)

    call basins_init(basins, mesh)
    call basins_analyze(basins, mesh, ff2(:,1), ff2, bader_threshold)

    if (debug%info) then
      call parse_variable(namespace, 'LDOutputFormat', 0, how)
      if (.not. varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
        call messages_input_error(namespace, 'LDOutputFormat')
      end if

      call dio_function_output(how, 'local.general', 'basinsmap', namespace, space, mesh, TOFLOAT(basins%map(1:mesh%np)), &
        unit_one, ierr, ions = ions)
      call dio_function_output(how, 'local.general', 'dens_ff2', namespace, space, mesh, ff2(:,1), unit_one, ierr, ions = ions)
    end if
    SAFE_DEALLOCATE_A(ff2)

    !%Variable LDUseAtomicRadii
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% If set, atomic radii will be used to assign lone pairs to ion. 
    !%End
    call parse_variable(global_namespace, 'LDUseAtomicRadii', .false., use_atomic_radii)

    if (use_atomic_radii) then
      SAFE_ALLOCATE(ion_map(ions%natoms))
      ion_map = 0
      do ia = 1, ions%natoms
        ion_map(ia) = basins%map(mesh_nearest_point(mesh, ions%pos(:, ia), dmin, rankmin))
      end do

      do ip = 1, mesh%np
        ! Check if lonely pair has no ions assigned to it
        if (all(ion_map(:) /= basins%map(ip)) ) then
          ! Assign lonely pair to ion in a atomic radii distance
          do ia = 1, ions%natoms
            dd = norm2(mesh%x(ip, :) - ions%pos(:, ia))
            if (dd <= species_vdw_radius(ions%atom(ia)%species)) basins%map(ip) = ion_map(ia)
          end do
        end if
      end do

      SAFE_DEALLOCATE_A(ion_map)
    end if

    POP_SUB(create_basins)
  end subroutine create_basins

  ! ---------------------------------------------------------
  subroutine box_domain_create_mask(domain, mesh)
    class(local_domain_t), intent(inout) :: domain
    type(mesh_t),          intent(in)    :: mesh

    PUSH_SUB(box_domain_create_mask)

    domain%mesh_mask = domain%box%contains_points(mesh%np, mesh%x)

    POP_SUB(box_domain_create_mask)
  end subroutine box_domain_create_mask
  
  ! ---------------------------------------------------------
  subroutine bader_domain_create_mask(domain, basins, mesh, ions)
    class(local_domain_t), intent(inout) :: domain
    type(basins_t),        intent(in)    :: basins
    type(mesh_t),          intent(in)    :: mesh
    type(ions_t),          intent(in)    :: ions

    integer               :: ia, ib, ip, ix, n_basins, rankmin
    FLOAT                 :: dmin, xi(MAX_DIM)
    integer, allocatable  :: domain_map(:)

    PUSH_SUB(bader_domain_create_mask)

    ! Count then number of basins to consider. That is the number of ions specified in the input by the user
    n_basins = 0
    do ia = 1, ions%natoms
      if (loct_isinstringlist(ia, domain%clist)) n_basins = n_basins + 1
    end do

    SAFE_ALLOCATE(domain_map(1:n_basins))

    ! Assign basins to ions
    domain_map = 0
    ib = 0
    do ia = 1, ions%natoms
      if (.not. loct_isinstringlist(ia, domain%clist)) cycle
      ib = ib + 1
      xi(1:ions%space%dim) = ions%pos(:, ia)
      ix = mesh_nearest_point(mesh, xi, dmin, rankmin)
      domain_map(ib) = basins%map(ix)
    end do

    ! Now construct the mask: a point belongs to this Bader domain if it is part
    ! of at least one basin that is part of this domain
    do ip = 1, mesh%np
      domain%mesh_mask(ip) = any(domain_map(1:n_basins) == basins%map(ip))
    end do

    SAFE_DEALLOCATE_A(domain_map)

    POP_SUB(bader_domain_create_mask)
  end subroutine bader_domain_create_mask

  ! ---------------------------------------------------------
  !Check for the ions inside each local domain.
  subroutine local_ions_mask(mesh_mask, ions, mesh, ions_mask)
    logical,            intent(in)  :: mesh_mask(:)
    type(ions_t),       intent(in)  :: ions
    type(mesh_t),       intent(in)  :: mesh
    logical,            intent(out) :: ions_mask(:)

    integer              :: ia, ix
    integer, allocatable :: mask_tmp(:)
    integer              :: rankmin
    FLOAT                :: dmin, xi(MAX_DIM)

    PUSH_SUB(local_ions_mask)

    SAFE_ALLOCATE(mask_tmp(ions%natoms))
    mask_tmp = 0

    do ia = 1, ions%natoms
      xi(1:ions%space%dim) = ions%pos(:, ia)
      ix = mesh_nearest_point(mesh, xi, dmin, rankmin )
      if (rankmin /= mesh%mpi_grp%rank) cycle
      if (mesh_mask(ix)) mask_tmp(ia) = 1
    end do

    if(mesh%parallel_in_domains) then
      call mesh%allreduce(mask_tmp, ions%natoms)
    end if                               
    ions_mask = mask_tmp == 1

    SAFE_DEALLOCATE_A(mask_tmp)

    POP_SUB(local_ions_mask)
  end subroutine local_ions_mask

end program oct_local_multipoles

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
