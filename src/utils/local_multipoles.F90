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
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use simul_box_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  type local_domain_t
    integer                         :: nd              !< number of local domains.
    type(box_union_t), allocatable  :: domain(:)       !< boxes that form each domain.
    character(len=500), allocatable :: clist(:)        !< list of centers for each domain.
    character(len=15), allocatable  :: lab(:)          !< declared name for each domain.
    integer, allocatable            :: dshape(:)       !< shape of box for each domain.
    logical, allocatable            :: inside(:,:)     !< relation of mesh points on each domain.
    logical, allocatable            :: ions_inside(:,:)!< relation of ions inside each domain.
    FLOAT, allocatable              :: dcm(:,:)        !< store the center of mass of each domain on the real space.
    type(local_write_t)             :: writ            !< write option for local domains analysis.
  end type local_domain_t

  type(system_t)        :: sys
  type(simul_box_t)     :: sb
  integer, parameter    :: BADER = 512
  FLOAT                 :: BaderThreshold
  type(namespace_t)     :: default_namespace

  ! Initialize stuff
  call global_init(is_serial = .false.)
  call calc_mode_par_init()

  call parser_init()
  default_namespace = namespace_t("")

  call messages_init(default_namespace)

  call io_init(default_namespace)
  call profiling_init(default_namespace)

  call print_header()
  call messages_print_stress(stdout, "Local Domains mode")
  call messages_print_stress(stdout)

  call messages_experimental("oct-local_multipoles utility")

  call unit_system_init(default_namespace)
  call restart_module_init(default_namespace)

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  call system_init(sys, default_namespace)
  call simul_box_init(sb, default_namespace, sys%geo, sys%space)

  call local_domains()

  call simul_box_end(sb)
  call system_end(sys)
  call profiling_end(default_namespace)
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
    type(local_domain_t)           :: local
    integer                        :: err, iter, l_start, l_end, l_step
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

    write(folder_default,'(a)')'restart/gs/'

    !%Variable LDFolder
    !%Type string
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The folder name where the density used as input file is.
    !%End
    call parse_variable(default_namespace, 'LDFolder', folder_default, folder)

    ! Check if the folder is finished by an /
    if(index(folder, '/', back = .true.) /= len_trim(folder)) then
      folder = trim(folder)//'/'
    end if

    default_dt = M_ZERO
    call parse_variable(default_namespace, 'TDTimeStep', default_dt, dt, unit = units_inp%time)
    if (dt <= M_ZERO) then
      write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
      write(message(2),'(a)') 'Input: TDTimeStep reset to 0. Check input file'
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
    call parse_variable(default_namespace, 'LDFilename', 'density', basename)
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
    !% 0.01 -> intramolecular volumes; 0.2 -> intermolecular volumes
    !%End
    call parse_variable(default_namespace, 'LDBaderThreshold', CNST(0.01), BaderThreshold)

    !%Variable LDUpdate
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Controls if the calculation of the local domains is desired at each iteration.
    !%End
    call parse_variable(default_namespace, 'LDUpdate', .false., ldupdate)

    !%Variable LDOverWrite
    !%Type logical
    !%Default true
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Controls whether to over-write existing files.
    !%End
    call parse_variable(default_namespace, 'LDOverWrite', .true., ldoverwrite)

    !%Variable LDRadiiFile
    !%Type string
    !%Default 'default'
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Full path for the radii file. If set, def_rsize will be reset to the new values.
    !% This file should have the same format as share/PP/default.
    !%End
    call parse_variable(default_namespace, 'LDRadiiFile', 'default', radiifile)

    if(trim(radiifile) /= "default") then
      n_spec_def = max(0, loct_number_of_lines(radiifile))
      if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

      ! get parameters from file
      do ia = 1, sys%geo%nspecies
        read_data = 0
        iunit = io_open(radiifile, default_namespace, action='read', status='old', die=.false.)
        if(iunit > 0) then
          if(ia == 1) then
            write(message(1),'(a,a)')'Redefining def_rsize from file:', trim(radiifile)
            call messages_info(1)
          end if
          read(iunit,*)
          default_file: do ispec = 1, n_spec_def
            read(iunit,*) lab
            if(trim(lab) == trim(species_label(sys%geo%species(ia)))) then
              call read_from_default_file(iunit, read_data, sys%geo%species(ia))
              exit default_file
            end if
          end do default_file
          call io_close(iunit)
        else
          write(message(1),'(a,a)')' Octopus could not open then file:',trim(radiifile)
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
    call parse_variable(default_namespace, 'LDRestart', .false., ldrestart)

    if (ldrestart) then
      write(folder_default,'(a)')'ld.general'

      !%Variable LDRestartFolder
      !%Type string
      !%Default "ld.general"
      !%Section Utilities::oct-local_multipoles
      !%Description
      !% The folder name where the density used as input file is.
      !%End
      call parse_variable(default_namespace, 'LDRestartFolder', folder_default, ldrestart_folder)

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
    call parse_variable(default_namespace, 'LDIterateFolder', .false., iterate)

    !%Variable LDStart
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The starting number of the filename or folder.
    !%End
    call parse_variable(default_namespace, 'LDStart', 0, l_start)

    !%Variable LDEnd
    !%Type integer
    !%Default 0
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The last number of the filename or folder.
    !%End
    call parse_variable(default_namespace, 'LDEnd', 0, l_end)

    !%Variable LDStep
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The padding between the filenames or folder.
    !%End
    call parse_variable(default_namespace, 'LDStep', 1, l_step)

    message(1) = 'Info: Computing local multipoles'
    message(2) = ''
    call messages_info(2)

    call local_init(local)

    ! Starting loop over selected densities.
    if ( any( local%dshape(:) == BADER )) then
      call messages_experimental('Bader volumes in oct-local_multipoles')
    end if

    call kick_init(kick, default_namespace, sys%gr%mesh%sb, sys%st%d%ispin)
    call local_write_init(local%writ, default_namespace, local%nd, local%lab, 0, dt)

    !TODO: initialize hamiltonian if needed: check for LDOuput = energy or potential, using local_write_check_hm(local%writ)

    if (ldrestart) then
      !TODO: check for domains & mesh compatibility
      call restart_init(restart_ld, default_namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, sys%mc, err, &
        dir=trim(ldrestart_folder), mesh = sys%gr%mesh)
      call local_restart(local, restart_ld)
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
    call restart_init(restart, default_namespace, RESTART_UNDEFINED, RESTART_TYPE_LOAD, sys%mc, err, &
      dir=trim(restart_folder), mesh = sys%gr%mesh)

!!$    call loct_progress_bar(-1, l_end-l_start)
    do iter = l_start, l_end, l_step
      if (iterate) then
        ! Delete the last / and add the corresponding folder number
        write(in_folder,'(a,i0.7,a)') folder(folder_index+1:len_trim(folder)-1),iter,"/"
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
      if(err == 0) &
        call drestart_read_mesh_function(restart, trim(filename), sys%gr%mesh, sys%st%rho(:,1), err)
      if (err /= 0 ) then
        write(message(1),*) 'While reading density: "', trim(filename), '", error code:', err
        call messages_fatal(1)
      end if

      ! Look for the mesh points inside local domains
      if(iter == l_start .and. .not.ldrestart) then
        call local_inside_domain(local, default_namespace, sys%st%rho(:,1))
      else
        if (ldupdate) call local_inside_domain(local, default_namespace, sys%st%rho(:,1))
      end if

      call local_write_iter(local%writ, default_namespace, local%nd, local%lab, local%ions_inside, local%inside, local%dcm, &
        sys%gr, sys%st, sys%hm, sys%ks, sys%geo, kick, iter, l_start, ldoverwrite)
      call loct_progress_bar(iter-l_start, l_end-l_start)
    end do

    call restart_end(restart)
    call kick_end(kick)
    call local_end(local)

    message(1) = 'Info: Exiting local domains'
    message(2) = ''
    call messages_info(2)

    POP_SUB(local_domains)
  end subroutine local_domains

  ! ---------------------------------------------------------
  !> Initialize local_domain_t variable, allocating variable
  !! and reading parameters from input file.
  ! ---------------------------------------------------------
  subroutine local_init(local)
    type(local_domain_t), intent(inout) :: local

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
    local%nd = 0
    if(parse_block(default_namespace, 'LocalDomains', blk) == 0) then
      local%nd = parse_block_n(blk)
    end if

    SAFE_ALLOCATE(local%domain(1:local%nd))
    SAFE_ALLOCATE(local%clist(1:local%nd))
    SAFE_ALLOCATE(local%dshape(1:local%nd))
    SAFE_ALLOCATE(local%lab(1:local%nd))
    SAFE_ALLOCATE(local%inside(1:sys%gr%mesh%np, 1:local%nd))
    SAFE_ALLOCATE(local%ions_inside(1:sys%geo%natoms, 1:local%nd))
    SAFE_ALLOCATE(local%dcm(1:sys%space%dim, 1:local%nd))

    block: do id = 1, local%nd
      call parse_block_string(blk, id-1, 0, local%lab(id))
      call local_read_from_block(blk, id-1, local%domain(id), local%dshape(id), local%clist(id), default_namespace)
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
  subroutine local_end(local)
    type(local_domain_t), intent(inout) :: local
    integer :: id

    PUSH_SUB(local_end)
    do id = 1, local%nd
      call box_union_end(local%domain(id))
    end do
    call local_write_end(local%writ, local%nd)
    SAFE_DEALLOCATE_A(local%lab)
    SAFE_DEALLOCATE_A(local%domain)
    SAFE_DEALLOCATE_A(local%clist)
    SAFE_DEALLOCATE_A(local%dshape)
    SAFE_DEALLOCATE_A(local%inside)
    SAFE_DEALLOCATE_A(local%dcm)

    POP_SUB(local_end)
  end subroutine local_end

  ! ---------------------------------------------------------
  subroutine local_read_from_block(blk, row, dom, shape, clist, namespace)
    type(block_t),     intent(in)        :: blk
    integer,           intent(in)        :: row
    type(box_union_t), intent(inout)     :: dom
    integer,           intent(out)       :: shape
    character(len=*),  intent(out)       :: clist
    type(namespace_t), intent(in)        :: namespace

    integer                     :: dim, ic, idir, nb
    FLOAT                       :: lgst, val, rsize, xsize
    FLOAT                       :: center(MAX_DIM), lsize(MAX_DIM)

    PUSH_SUB(local_read_from_block)

    ! Initializing variables in dom
    shape = 1
    nb = 1
    rsize = -M_ONE
    xsize = M_ZERO
    lsize(:) = M_ZERO
    center(:) = M_ZERO
    dim = sys%space%dim

    call parse_block_integer(blk, row, 1, shape)

    select case(shape)
    case(MINIMUM)
      if(sys%geo%reduced_coordinates) then
        message(1) = "The 'minimum' box shape cannot be used if atomic positions"
        message(2) = "are given as reduced coordinates."
        call messages_fatal(2)
      end if
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if(rsize < M_ZERO) call messages_input_error('radius')
      call parse_block_string(blk, row, 3, clist)
      nb = 0
      do ic = 1, sys%geo%natoms
        if(loct_isinstringlist(ic, clist)) nb = nb + 1
      end do
    case(SPHERE)
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if(rsize < M_ZERO) call messages_input_error('radius')
      do ic = 1, dim
        call parse_block_float(blk, row, 2 + ic, center(ic), unit = units_inp%length)
      end do
    case(CYLINDER)
      call parse_block_float(blk, row, 2, rsize, unit = units_inp%length)
      if(rsize < M_ZERO) call messages_input_error('radius')
      call parse_block_float(blk, row, 3, xsize, unit = units_inp%length)
      do ic = 1, dim
        call parse_block_float(blk, row, 3 + ic, center(ic), unit = units_inp%length)
      end do
    case(PARALLELEPIPED)
      do ic = 1, dim
        call parse_block_float(blk, row, 2 + ic, lsize(ic), unit = units_inp%length)
      end do
      do ic = 1, dim
        call parse_block_float(blk, row, 2 + dim + ic, center(ic), unit = units_inp%length)
      end do
    case(BADER)
      ! FIXME: when input error exists --> segmentation fault appears
      call parse_block_string(blk, row, 2, clist)
      nb = 0
      do ic = 1, sys%geo%natoms
        if(loct_isinstringlist(ic, clist)) nb = nb + 1
      end do
    end select
    ! fill in lsize structure
    select case(shape)
    case(SPHERE)
      lsize(1:dim) = rsize
    case(CYLINDER)
      lsize(1)     = xsize
      lsize(2:dim) = rsize
    case(MINIMUM, BADER)
      do idir = 1, dim
        lgst = M_ZERO; val = M_ZERO
        do ic = 1, sys%geo%natoms
          if(loct_isinstringlist(ic, clist)) val = abs(sys%geo%atom(ic)%x(idir))
          if(lgst <= val) lgst = val
        end do
        lsize(idir) =  lgst + rsize
      end do
    end select
    call local_domains_init(dom, dim, shape, center, rsize, lsize, nb, clist, namespace)

    POP_SUB(local_read_from_block)

  end subroutine local_read_from_block

  !!---------------------------------------------------------------------------^
  subroutine local_domains_init(dom, dim, shape, center, rsize, lsize, nb, clist, namespace)
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
    FLOAT                    :: bcenter(dim), bsize(dim)
    type(box_t), allocatable :: boxes(:)

    PUSH_SUB(local_domains_init)

    SAFE_ALLOCATE(boxes(1:nb))
    bsize(:) = M_ZERO
    ibox = 1

    select case(shape)
    case(SPHERE, CYLINDER)
      call box_create(boxes(ibox), shape, dim, lsize, center, namespace)
    case(PARALLELEPIPED)
      bshape         = 3
      call box_create(boxes(ibox), bshape, dim, lsize, center, namespace)
    case(MINIMUM)
      bshape         = SPHERE
      do ia = 1, sys%geo%natoms
        if(loct_isinstringlist(ia, clist))then
          bcenter(1:dim) = sys%geo%atom(ia)%x(1:dim)
          bsize(:) = rsize
          if (bsize(1) < M_EPSILON) bsize(1) = species_def_rsize(sys%geo%atom(ia)%species)
          call box_create(boxes(ibox), bshape, dim, bsize, bcenter, namespace)
          ibox = ibox + 1
        end if
      end do
    case(BADER)
      bshape         = SPHERE
      do ia = 1, sys%geo%natoms
        if(loct_isinstringlist(ia, clist))then
          bcenter(1:dim) = sys%geo%atom(ia)%x(1:dim)
          bsize(:) = species_def_h(sys%geo%atom(ia)%species)
          call box_create(boxes(ibox), bshape, dim, bsize, bcenter, namespace)
          ibox = ibox + 1
        end if
      end do
    end select

    call box_union_init(dom, nb, boxes)

    ! TODO: Check for a conflict between box_union and clist for Bader Volumes
    ic = 0
    if (shape == MINIMUM ) then
      do ia = 1, sys%geo%natoms
        if (box_union_inside(dom, sys%geo%atom(ia)%x).and. .not.loct_isinstringlist(ia, clist) ) then
          ic = ic + 1
          if(ic <= 20) write(message(ic),'(a,a,I0,a,a)')'Atom: ',trim(species_label(sys%geo%atom(ia)%species)),ia, &
            ' is inside the union box BUT not in list: ',trim(clist)
        end if
      end do
      if (ic > 0) then
        if( ic > 20 ) ic = 20
        call messages_warning(ic)
        message(1) = ' THIS COULD GIVE INCORRECT RESULTS '
        call messages_warning(1)
        if ( ic == 20 ) then
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
  subroutine local_restart(lcl, restart)
    type(local_domain_t), intent(inout) :: lcl
    type(restart_t),      intent(in)    :: restart

    integer                     :: id, ip, iunit, ierr
    character(len=MAX_PATH_LEN) :: filename, tmp
    character(len=31) :: line(1)
    FLOAT, allocatable          :: inside(:)

    PUSH_SUB(local_restart)

    message(1) = 'Info: Reading mesh points inside each local domain'
    call messages_info(1)
    SAFE_ALLOCATE(inside(1:sys%gr%mesh%np))
    !Read local domain information from ldomains.info
    filename = "ldomains.info"
    iunit = restart_open(restart, filename, status='old')
    call restart_read(restart, iunit, line, 1, ierr)
    read(line(1),'(a25,1x,i5)') tmp, ierr
    call restart_close(restart, iunit)

    filename = "ldomains"
    call drestart_read_mesh_function(restart, trim(filename), sys%gr%mesh, inside, ierr)

    lcl%inside = .false.
    do ip = 1 , sys%gr%mesh%np
      do id = 1, lcl%nd
        if (bitand(int(inside(ip)), 2**id) /= 0) lcl%inside(ip,id) = .true.
      end do
    end do

    !Check for atom list inside each domain
    call local_ions_inside(lcl%nd, lcl%inside, sys%geo, sys%gr%mesh, lcl%ions_inside, lcl%clist)

    !Compute center of mass of each domain
    call local_center_of_mass(lcl%nd, lcl%ions_inside, sys%geo, lcl%dcm)

    SAFE_DEALLOCATE_A(inside)

    POP_SUB(local_restart)
  end subroutine local_restart
  ! ---------------------------------------------------------
  subroutine local_inside_domain(lcl, namespace, ff)
    type(local_domain_t),   intent(inout) :: lcl
    type(namespace_t),      intent(in)    :: namespace
    FLOAT,                  intent(in)    :: ff(:)

    integer(8)          :: how
    integer             :: id, ip, iunit, ierr
    type(basins_t)      :: basins
    FLOAT, allocatable  :: ff2(:,:)
    logical             :: extra_write
    character(len=MAX_PATH_LEN)   :: filename, base_folder, folder, frmt
    character(len=140), allocatable  :: lines(:)
    type(restart_t)     :: restart

    PUSH_SUB(local_inside_domain)

    !TODO: find a parallel algorithm to perform the charge-density fragmentation.
    message(1) = 'Info: Assigning mesh points inside each local domain'
    call messages_info(1)

    SAFE_ALLOCATE(ff2(1:sys%gr%mesh%np,1))
    ff2(1:sys%gr%mesh%np,1) = ff(1:sys%gr%mesh%np)

    if (any(lcl%dshape(:) == BADER)) then

      if (sys%gr%mesh%parallel_in_domains) then
        write(message(1),'(a)') 'Bader volumes can only be computed in serial'
        call messages_fatal(1)
      end if

      call add_dens_to_ion_x(ff2, namespace, sys%geo)
      call basins_init(basins, sys%gr%mesh)
      call parse_variable(namespace, 'LDBaderThreshold', CNST(0.01), BaderThreshold)
      call basins_analyze(basins, sys%gr%mesh, ff2(:,1), ff2, BaderThreshold)
      call parse_variable(namespace, 'LDExtraWrite', .false., extra_write)
      call bader_union_inside(basins, lcl%nd, lcl%domain, lcl%lab, lcl%dshape, lcl%inside)

      if (extra_write) then
        call messages_obsolete_variable(namespace, 'LDOutputHow', 'LDOutputFormat')
        call parse_variable(namespace, 'LDOutputFormat', 0, how)
        if(.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
          call messages_input_error('LDOutputFormat')
        end if
        filename = 'basinsmap'
        call dio_function_output(how, trim('local.general'), trim(filename), default_namespace, &
          sys%gr%mesh, DBLE(basins%map(1:sys%gr%mesh%np)), unit_one, ierr, geo = sys%geo)
        call dio_function_output(how, trim('local.general'), 'dens_ff2', default_namespace, &
          sys%gr%mesh, ff2(:,1), unit_one, ierr, geo = sys%geo)
        call io_close(iunit)
      end if
      call basins_end(basins)
    else
      do id = 1, lcl%nd
        call box_union_inside_vec(lcl%domain(id), sys%gr%mesh%np, sys%gr%mesh%x, lcl%inside(:,id))
      end do
    end if

    !Check for atom list inside each domain
    call local_ions_inside(lcl%nd, lcl%inside, sys%geo, sys%gr%mesh, lcl%ions_inside, lcl%clist)

    !Compute center of mass of each domain
    call local_center_of_mass(lcl%nd, lcl%ions_inside, sys%geo, lcl%dcm)

    !Write restart file for local domains
    base_folder = "./restart/"
    folder = "ld/"
    filename = "ldomains"
    write(message(1),'(a,a)')'Info: Writing restart info to ', trim(filename)
    call messages_info(1)
    call restart_init(restart, namespace, RESTART_UNDEFINED, RESTART_TYPE_DUMP, sys%mc, ierr, &
      mesh=sys%gr%mesh, dir=trim(base_folder)//trim(folder))
    ff2 = M_ZERO
    SAFE_ALLOCATE(lines(1:lcl%nd+2))
    write(lines(1),'(a,1x,i5)')'Number of local domains =', lcl%nd
    write(lines(2),'(a3,1x,a15,1x,a5,1x,71x,a9,1x,a14)')'#id','label','shape','Atom list','center of mass'
    do id = 1, lcl%nd
      write(frmt,'(a,i0,a)')'(i3,1x,a15,1x,i5,1x,a80,1x',sys%space%dim,'(f10.6,1x))'
      write(lines(id+2),fmt=trim(frmt))id, trim(lcl%lab(id)), lcl%dshape(id), trim(lcl%clist(id)),lcl%dcm(1:sys%space%dim, id)
      do ip = 1, sys%gr%mesh%np
        if (lcl%inside(ip, id)) ff2(ip,1) = ff2(ip,1) + 2**DBLE(id)
      end do
    end do
    call drestart_write_mesh_function(restart, filename, sys%gr%mesh, ff2(1:sys%gr%mesh%np, 1), ierr)

    filename = "ldomains.info"
    iunit = restart_open(restart, filename)
    call restart_write(restart, iunit, lines, lcl%nd+2, ierr)
    call restart_end(restart)
    call io_close(iunit)
    SAFE_DEALLOCATE_A(lines)

    SAFE_DEALLOCATE_A(ff2)

    POP_SUB(local_inside_domain)

  end subroutine local_inside_domain

  ! ---------------------------------------------------------
  subroutine bader_union_inside(basins, nd, dom, lab, dsh, inside)
    type(basins_t),    intent(inout) :: basins
    integer,           intent(in)    :: nd
    type(box_union_t), intent(in)    :: dom(:)
    character(len=15), intent(in)    :: lab(:)
    integer,           intent(in)    :: dsh(:)
    logical,           intent(out)   :: inside(:,:)

    integer(8)            :: how
    integer               :: ia, ib, id, ierr, ip, ix, rankmin
    integer               :: max_check
    integer, allocatable  :: dunit(:), domain_map(:,:), ion_map(:)
    FLOAT                 :: dmin, dd
    FLOAT, allocatable    :: xi(:), dble_domain_map(:,:), domain_mesh(:)
    logical               :: extra_write, lduseatomicradii
    character(len=64)     :: filename

    PUSH_SUB(bader_union_inside)

    SAFE_ALLOCATE(xi(1:sys%space%dim))
    inside = .false.


    !%Variable LDUseAtomicRadii
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% If set, atomic radii will be used to assign lone pairs to ion.
    !%End
    call parse_variable(default_namespace, 'LDUseAtomicRadii', .false., lduseatomicradii)

    SAFE_ALLOCATE(ion_map(1:sys%geo%natoms))

    max_check = 1
    do id = 1, nd
      if (box_union_get_nboxes(dom(id)) > max_check) max_check = box_union_get_nboxes(dom(id))
    end do
    SAFE_ALLOCATE(domain_map(1:nd, 1:max_check))

    do id = 1, nd
      if( dsh(id) /= BADER ) then
        call box_union_inside_vec(dom(id), sys%gr%mesh%np, sys%gr%mesh%x, inside(:,id))
      else
        ! Assign basins%map to ions
        if (lduseatomicradii .and. (id == 1)) then
          ion_map = 0
          do ia = 1, sys%geo%natoms
            ion_map(ia) = basins%map(mesh_nearest_point(sys%gr%mesh, sys%geo%atom(ia)%x, dmin, rankmin))
          end do
          ! Assign lonely pair to ion in a atomic radii distance
          do ip = 1, sys%gr%mesh%np
            if( all(ion_map(:) /= basins%map(ip)) ) then
              do ia = 1, sys%geo%natoms
                dd = sum((sys%gr%mesh%x(ip, 1:sys%gr%mesh%sb%dim) &
                  - sys%geo%atom(ia)%x(1:sys%gr%mesh%sb%dim))**2)
                dd = sqrt(dd)
                if ( dd <= species_vdw_radius(sys%geo%atom(ia)%species) ) basins%map(ip) = ion_map(ia)
              end do
            end if
          end do
        end if
        domain_map = 0
        do ib = 1, box_union_get_nboxes(dom(id))
          xi = box_union_get_center(dom(id), ib)
          ix = mesh_nearest_point(sys%gr%mesh, xi, dmin, rankmin)
          domain_map(id,ib) = basins%map(ix)
        end do
        do ip = 1, sys%gr%mesh%np
          if(any( domain_map(id, 1:box_union_get_nboxes(dom(id))) == basins%map(ip) ) ) inside(ip,id) = .true.
        end do
      end if
    end do

    SAFE_DEALLOCATE_A(domain_map)
    SAFE_DEALLOCATE_A(ion_map)

    !%Variable LDExtraWrite
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Writes additional information to files, when computing local multipoles. For
    !% example, it writes coordinates of each local domain.
    !%End
    call parse_variable(default_namespace, 'LDExtraWrite', .false., extra_write)

    SAFE_ALLOCATE(dble_domain_map(1:nd, 1:sys%gr%mesh%np))
    SAFE_ALLOCATE(domain_mesh(1:sys%gr%mesh%np))

    if (extra_write) then
      call parse_variable(default_namespace, 'LDOutputFormat', 0, how)
      if(.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
        call messages_input_error('LDOutputFormat')
      end if
      dble_domain_map(1:nd, 1:sys%gr%mesh%np) = M_ZERO
      domain_mesh(1:sys%gr%mesh%np) = M_ZERO
      do ip = 1, sys%gr%mesh%np
        do id = 1, nd
          if (inside(ip, id)) then
            dble_domain_map(id, ip) = DBLE(id)
            domain_mesh(ip) = domain_mesh(ip) + dble_domain_map(id, ip)
          end if
        end do
      end do

      write(filename,'(a,a,a)')'domain.mesh'
      call dio_function_output(how, trim('local.general'), trim(filename), default_namespace, &
        sys%gr%mesh, domain_mesh(1:sys%gr%mesh%np) , unit_one, ierr, geo = sys%geo)

      do id = 1, nd
        write(filename,'(a,a,a)')'domain.',trim(lab(id))
        call dio_function_output(how, trim('local.general'), trim(filename), default_namespace, &
          sys%gr%mesh, dble_domain_map(id, 1:sys%gr%mesh%np) , unit_one, ierr, geo = sys%geo)
      end do
    end if

    SAFE_DEALLOCATE_A(dble_domain_map)
    SAFE_DEALLOCATE_A(domain_mesh)

    SAFE_DEALLOCATE_A(xi)
    SAFE_DEALLOCATE_A(dunit)

    POP_SUB(bader_union_inside)

  end subroutine bader_union_inside

  ! ---------------------------------------------------------
  subroutine add_dens_to_ion_x(ff, namespace, geo)
    FLOAT,              intent(inout) :: ff(:,:)
    type(namespace_t),  intent(in)    :: namespace
    type(geometry_t),   intent(inout) :: geo

    integer :: ia, is
    FLOAT, allocatable :: ffs(:)

    PUSH_SUB(add_dens_to_ion_x)

    SAFE_ALLOCATE(ffs(1:sys%gr%mesh%np))
    do ia = 1, geo%natoms
      call species_get_density(geo%atom(ia)%species, namespace, geo%atom(ia)%x, sys%gr%mesh, ffs)
      do is = 1, sys%st%d%nspin
        ff(1:sys%gr%mesh%np,is) = ff(1:sys%gr%mesh%np, is) - ffs(1:sys%gr%mesh%np)
      end do
    end do

    SAFE_DEALLOCATE_A(ffs)
    POP_SUB(add_dens_to_ion_x)
  end subroutine add_dens_to_ion_x
  ! ---------------------------------------------------------
  !Check for the ions inside each local domain.
  subroutine local_ions_inside(nd, inside, geo, mesh, ions_inside, ions_list)
    integer,           intent(in)  :: nd
    logical,           intent(in)  :: inside(:,:)
    type(geometry_t),  intent(in)  :: geo
    type(mesh_t),      intent(in)  :: mesh
    logical,           intent(out) :: ions_inside(:,:)
    character(len=500),intent(out) :: ions_list(:)

    integer              :: ia, id, ix, ic
    integer, allocatable :: inside_tmp(:,:)
    integer              :: rankmin
    FLOAT                :: dmin
    character(len=500)   :: chtmp

    PUSH_SUB(local_ions_inside)

    SAFE_ALLOCATE(inside_tmp(geo%natoms,nd))
    inside_tmp = 0
    ions_list = ""

    ions_inside = .false.
    do ia = 1, geo%natoms
      ix = mesh_nearest_point(mesh, geo%atom(ia)%x, dmin, rankmin )
      if (rankmin /= mesh%mpi_grp%rank) cycle
      do id = 1, nd
        if (inside(ix, id)) inside_tmp(ia, id) = 1
      end do
    end do

    if(mesh%parallel_in_domains) then
      call comm_allreduce(mesh%mpi_grp%comm, inside_tmp, dim = (/geo%natoms,nd/))
    end if

    do ia = 1, geo%natoms
      do id = 1, nd
        if (inside_tmp(ia, id) == 1) ions_inside(ia, id) = .true.
      end do
    end do

    SAFE_DEALLOCATE_A(inside_tmp)

    !print list of atoms
    do id = 1, nd
      ic = 0
      ix = 0
      chtmp = ions_list(id)
      do ia = 1, geo%natoms-1
        if (ions_inside(ia, id)) then
          if ( ic == 0 .or. .not.ions_inside(ia+1, id)) then
            write(ions_list(id), '(a,i0)')trim(chtmp),ia
          end if
          if (ions_inside(ia+1,id)) then
            ic = ic + 1
            chtmp = trim(ions_list(id))//"-"
          else
            ic = 0
            chtmp = trim(ions_list(id))//","
          end if

        else
          cycle
        end if
      end do

      if (ions_inside(geo%natoms,id)) then
        write(ions_list(id), '(a,i0)')trim(chtmp),ia
      end if

!      write(message(1),'(a,1x,i0,1x,a,1x,a)')'Atoms inside domain',id,':',trim(ions_list(id))
!      call messages_warning(1)
    end do

    POP_SUB(local_ions_inside)
  end subroutine local_ions_inside
  ! ---------------------------------------------------------
  subroutine local_center_of_mass(nd, ions_inside, geo, center)
    integer,           intent(in)  :: nd
    logical,           intent(in)  :: ions_inside(:,:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(out) :: center(:,:)

    integer            :: ia, id
    FLOAT, allocatable :: sumw(:)

    PUSH_SUB(local_center_of_mass)

    SAFE_ALLOCATE(sumw(1:nd))
    sumw(1:nd) = M_ZERO

    center(:,:) = M_ZERO
    do ia = 1, geo%natoms
      do  id = 1, nd
        if ( ions_inside(ia, id) ) then
          center(1:geo%space%dim,id) = center(1:geo%space%dim,id) &
            + geo%atom(ia)%x(1:geo%space%dim)*species_mass(geo%atom(ia)%species)
          sumw(id) = sumw(id) + species_mass(geo%atom(ia)%species)
        end if
      end do
    end do
    do id = 1, nd
      center(1:geo%space%dim,id) = center(1:geo%space%dim,id) / sumw(id)
    end do

    SAFE_DEALLOCATE_A(sumw)

    POP_SUB(local_center_of_mass)
  end subroutine local_center_of_mass

end program oct_local_multipoles

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
