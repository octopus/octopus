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
!! $Id$

#include "global.h"

program oct_local_multipoles
  use atom_m
  use basins_m
  use box_m
  use box_union_m
  use datasets_m
  use command_line_m
  use geometry_m
  use global_m
  use hamiltonian_m
  use io_m
  use io_binary_m
  use io_function_m
  use kick_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use local_write_m
  use parser_m
  use profiling_m
  use restart_m
  use space_m
  use species_m
  use simul_box_m
  use system_m    
  use unit_m
  use unit_system_m
  use utils_m

  implicit none
  
  type local_domain_t
    integer                         :: nd
    type(box_union_t), allocatable  :: domain(:)
    character(len=15), allocatable  :: lab(:) 
    integer, allocatable            :: dshape(:)
    logical, allocatable            :: inside(:,:)
    FLOAT, allocatable              :: dcm(:,:)         !< store the center of mass of each domain on the real space.
    type(local_write_t)             :: writ       
  end type local_domain_t

  type(system_t)        :: sys
  type(simul_box_t)     :: sb
  integer, parameter    :: BADER = 512
  FLOAT                 :: BaderThreshold


  ! Initialize stuff
  call global_init(is_serial = .true.)

  call messages_init()

  call messages_experimental("oct-local_multipoles utility")

  call datasets_init(1)
  call io_init()
  call profiling_init()

  call unit_system_init()
  call system_init(sys)
  call simul_box_init(sb, sys%geo, sys%space)

  call local_domains()

  call simul_box_end(sb)
  call geometry_end(sys%geo)
  call space_end(sys%space)
  call profiling_output()
  call profiling_end()
  call io_end()
  call datasets_end()
  call messages_end()
  call global_end()

contains
  ! -------------
  !> Reads a global grid and density and make local domains
  !! This is a high-level interface that reads the input file and
  !! calls the proper function.
  subroutine local_domains()
    type(local_domain_t)           :: local
    integer                        :: err, iter, l_start, l_end, l_step, last_slash
    integer                        :: length
    FLOAT                          :: default_dt, dt
    character(64)                  :: filename, folder, folder_default, aux, base_folder
    logical                        :: iterate, ldupdate
    type(kick_t)                   :: kick

    PUSH_SUB(local_domains)
    
    message(1) = 'Info: Creating local domains'
    message(2) = ''
    call messages_info(2)

    write(folder_default,'(a)')'restart/gs/'

    !%Variable GlobalDensityFolder
    !%Type string
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The folder name where the input files are.
    !%End
    call parse_string(datasets_check('GlobalDensityFolder'), folder_default, folder)

    ! Check if the folder is finished by an /
    if (index(folder, '/', .true.) /= len_trim(folder)) then
      write(folder,'(a,a1)') trim(folder), '/'
    end if

    ! Guess the base folder (which should change while iterating)
    base_folder = ""
    last_slash = index(folder(1:len_trim(folder)-1), '/', .true.)
    if ( last_slash > 0 ) then
      base_folder = folder(1:last_slash)
      folder = folder(last_slash+1:len_trim(folder))
    end if
    
    aux = folder(1:3)
    if (aux == 'td.') then
      aux = trim(folder(4:len_trim(folder)-1))
      read(aux,'(I10.0)')iter
      default_dt = M_ZERO
      call parse_float(datasets_check('TDTimeStep'), default_dt, dt, unit = units_inp%time)
      if (dt <= M_ZERO) then
        write(message(1),'(a)') 'Input: TDTimeStep must be positive.'
        write(message(2),'(a)') 'Input: TDTimeStep reset to 0. Check input file'
        call messages_info(2)
      end if
    else
      iter = 0
      dt = M_ZERO
    end if

    !%Variable GlobalDensityFilename
    !%Type string
    !%Default 'density'
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Input filename. The original filename for the density which is going to be 
    !% fragmented into domains.
    !%End
    call parse_string(datasets_check('GlobalDensityFilename'), 'density', filename)
    if ( filename == " " ) filename = ""
    ! Delete the extension if present
    length = len_trim(filename)
    if ( length > 4) then
      if ( filename(length-3:length) == '.obf' ) then
        filename = trim(filename(1:length-4))
      end if
    end if

    !%Variable LocalBaderThreshold
    !%Type float
    !%Default 0.01
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable sets the threshold for the basins calculations. Recommended values: 
    !% Recommended values: 0.01 -> intramolecular volumes; 0.2 -> intermolecular volumes
    !%End
    call parse_float(datasets_check('LocalBaderThreshold'), CNST(0.01), BaderThreshold)

    !%Variable LDUpdate
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Controls if the calculation of the local domains is desired at each iteration.
    !%End
    call parse_logical(datasets_check('LDUpdate'), .false., ldupdate)

    ! Documentation in convert.F90
    call parse_logical(datasets_check('ConvertIterateFolder'), .false., iterate)
    call parse_integer(datasets_check('ConvertStart'), iter, l_start)
    call parse_integer(datasets_check('ConvertEnd'), iter, l_end)
    call parse_integer(datasets_check('ConvertStep'), 1, l_step)

    message(1) = 'Info: Computing local multipoles'
    message(2) = ''
    call messages_info(2)
    
    call local_init(local)

    ! Starting loop over selected densities.
    if ( any( local%dshape(:) == BADER )) then
      call messages_experimental('Bader volumes in oct-local_multipoles')
    end if

    call kick_init(kick,  sys%st%d%ispin, sys%gr%mesh%sb%dim, sys%gr%mesh%sb%periodic_dim )
    call local_write_init(local%writ, local%nd, local%lab, 0, dt)
    call loct_progress_bar(-1, l_end-l_start) 
    do iter = l_start, l_end, l_step
      if (iterate) then
        write(folder,'(a,i0.7,a)') folder(1:3),iter,"/"
      end if
      call drestart_read_function(trim(base_folder) // trim(folder), trim(filename), sys%gr%mesh, sys%st%rho(:,1), err)
      if (err /= 0 ) then
        write(message(1),*) 'While reading density: "', trim(base_folder) // trim(folder), trim(filename), '", error code:', err
        call messages_fatal(1)
      end if
      ! Look for the mesh points inside local domains
      if(iter == l_start) then
        call local_inside_domain(local, sys%st%rho(:,1), .true.)
      else
        call local_inside_domain(local, sys%st%rho(:,1), ldupdate)
      end if

      call local_write_iter(local%writ, local%nd, local%domain, local%lab, local%inside, local%dcm, & 
                              sys%gr, sys%st, sys%geo, kick, iter, dt)
      call loct_progress_bar(iter-l_start, l_end-l_start) 
    end do

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

    !%Variable LocalDomain
    !%Type block
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% The LocalDomains are by definition part of the global grid. The domains are defined by 
    !% selecting a type shape. The domain box will be constructed using the given parameters. 
    !% A local domain could be construct by addition of several box centered on the ions.
    !% The grid points inside this box will belong to the local domain. 
    !% 
    !% The format of this block is the following: 
    !% 'Label' | Shape | %< | Shape dependencies >%
    !%  The first field is the label of the domain. 
    !% Label = string with the name of the new local domain.
    !% The second is the shape type of the box used to define the domain.
    !% Shape = SPHERE, CYLINDER, PARALLELEPIPED, MINIMUM, BADER
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
    !%
    !% rsize < Radius in input length units
    !% xsize < the length of the cylinder in the x-direction 
    !% origin coordinates < in input length units separated by | . where is the box centered.
    !% lsize <  half of the length of the parallelepiped in each direction.
    !% center_list < string containing the list of atoms in xyz file for each domain in the form "2,16-23"
    !% 
    !%End

    ! First, find out if there is a LocalDomains block.
    local%nd = 0
    if(parse_block(datasets_check('LocalDomains'), blk) == 0) then
      local%nd = parse_block_n(blk)
    end if

    SAFE_ALLOCATE(local%domain(1:local%nd))
    SAFE_ALLOCATE(local%dshape(1:local%nd))
    SAFE_ALLOCATE(local%lab(1:local%nd))
    SAFE_ALLOCATE(local%inside(1:sys%gr%mesh%np, 1:local%nd))
    SAFE_ALLOCATE(local%dcm(1:sys%space%dim, 1:local%nd))

    block: do id = 1, local%nd
      call parse_block_string(blk, id-1, 0, local%lab(id))
      call local_read_from_block(blk, id-1, local%domain(id), local%dshape(id))
    end do block
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
    call local_write_end(local%writ)
    SAFE_DEALLOCATE_A(local%lab)
    SAFE_DEALLOCATE_A(local%domain)
    SAFE_DEALLOCATE_A(local%dshape)
    SAFE_DEALLOCATE_A(local%inside)
    SAFE_DEALLOCATE_A(local%dcm)

    POP_SUB(local_end)
  end subroutine local_end

  ! ---------------------------------------------------------
  subroutine local_read_from_block(blk, row, dom, shape)
    type(block_t),     intent(in)        :: blk
    integer,           intent(in)        :: row
    type(box_union_t), intent(inout)     :: dom
    integer,           intent(out)       :: shape
    
    integer           :: dim, ic, idir, nb
    character(len=80) :: clist
    FLOAT             :: lgst, val, rsize, xsize
    FLOAT             :: center(MAX_DIM), lsize(MAX_DIM)
   
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
        call parse_block_float(blk, row, 2, rsize)
        rsize = units_from_atomic(units_inp%length**(-1), rsize)
        if(rsize < M_ZERO) call input_error('radius')
        call parse_block_string(blk, row, 3, clist)
        nb = 0
        do ic = 1, sys%geo%natoms
          if(loct_isinstringlist(ic, clist)) nb = nb + 1
        end do
      case(SPHERE)
        call parse_block_float(blk, row, 2, rsize)
        rsize = units_from_atomic(units_inp%length**(-1), rsize)
        if(rsize < M_ZERO) call input_error('radius')
        do ic = 1, dim 
          call parse_block_float(blk, row, 2 + ic, center(ic))
          center(ic) = units_from_atomic(units_inp%length**(-1), center(ic))
        end do
      case(CYLINDER)
        call parse_block_float(blk, row, 2, rsize)
        rsize = units_from_atomic(units_inp%length**(-1), rsize)
        if(rsize < M_ZERO) call input_error('radius')
        call parse_block_float(blk, row, 3, xsize)
        do ic = 1, dim 
          call parse_block_float(blk, row, 3 + ic, center(ic))
          center(ic) = units_from_atomic(units_inp%length**(-1), center(ic))
        end do
      case(PARALLELEPIPED)
        do ic = 1, dim 
          call parse_block_float(blk, row, 2 + ic, lsize(ic))
          lsize(ic) = units_from_atomic(units_inp%length**(-1), lsize(ic))
        end do
        do ic = 1, dim 
          call parse_block_float(blk, row, 2 + dim + ic, center(ic))
          center(ic) = units_from_atomic(units_inp%length**(-1), center(ic))
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
    call local_domains_init(dom, dim, shape, center, rsize, lsize, nb, clist)

    POP_SUB(local_read_from_block)

  end subroutine local_read_from_block

  !!---------------------------------------------------------------------------^
  subroutine local_domains_init(dom, dim, shape, center, rsize, lsize, nb, clist)
    type(box_union_t), intent(inout) :: dom
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: shape
    FLOAT,             intent(in)    :: center(dim)
    FLOAT,             intent(in)    :: rsize
    FLOAT,             intent(in)    :: lsize(MAX_DIM)
    integer,           intent(in)    :: nb
    character(len=80), intent(in)    :: clist

    integer                  :: ia, ibox, ic, bshape
    FLOAT                    :: bcenter(dim), bsize(dim)
    type(box_t), allocatable :: boxes(:)

    PUSH_SUB(local_domains_init)

    SAFE_ALLOCATE(boxes(1:nb))
    bsize(:) = M_ZERO
    ibox = 1

    select case(shape)
      case(SPHERE, CYLINDER) 
        call box_create(boxes(ibox), shape, dim, lsize, center)
      case(PARALLELEPIPED) 
        bshape         = 3
        call box_create(boxes(ibox), bshape, dim, lsize, center)
      case(MINIMUM) 
        bshape         = SPHERE
        do ia = 1, sys%geo%natoms
          if(loct_isinstringlist(ia, clist))then
            bcenter(1:dim) = sys%geo%atom(ia)%x(1:dim)
            bsize(:) = rsize
            if (bsize(1) < M_EPSILON) bsize(1) = species_def_rsize(sys%geo%atom(ia)%spec)
            call box_create(boxes(ibox), bshape, dim, bsize, bcenter)
            ibox = ibox + 1
          end if
        end do
      case(BADER) 
        bshape         = SPHERE
        do ia = 1, sys%geo%natoms
          if(loct_isinstringlist(ia, clist))then
            bcenter(1:dim) = sys%geo%atom(ia)%x(1:dim)
            bsize(:) = species_def_rsize(sys%geo%atom(ia)%spec)
            call box_create(boxes(ibox), bshape, dim, bsize, bcenter)
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
          if ( shape == BADER ) then
            write(message(1),'(a,a,I0,a,a)')'Atom: ',trim(species_label(sys%geo%atom(ia)%spec)),ia, &
                                             ' is inside the union box BUT not in list: ',trim(clist)
            message(2) = ' BUG: This could not be possible for BADER SHAPED BOXES'
            call messages_fatal(2)
          end if
          if( ic <= 20 ) write(message(ic),'(a,a,I0,a,a)')'Atom: ',trim(species_label(sys%geo%atom(ia)%spec)),ia, & 
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
  subroutine local_inside_domain(lcl, ff, update)
    type(local_domain_t),   intent(inout) :: lcl
    FLOAT,                  intent(in)    :: ff(:)
    logical,                intent(in)    :: update
    
    integer             :: id, ip, ix, iunit
    type(basins_t)      :: basins
    FLOAT, allocatable  :: ff2(:,:)
    logical             :: extra_write
    character(len=64)   :: filename
    
    PUSH_SUB(local_inside_domain)

    if (update) then
      if (any(lcl%dshape(:) == BADER)) then
        SAFE_ALLOCATE(ff2(1:sys%gr%mesh%np_part,1)); ff2(1:sys%gr%mesh%np_part,1) = ff(1:sys%gr%mesh%np_part)
        call add_dens_to_ion_x(ff2,sys%geo)
        call basins_init(basins, sys%gr%mesh)
        call parse_float(datasets_check('LocalBaderThreshold'), CNST(0.01), BaderThreshold)
        call basins_analyze(basins, sys%gr%mesh, ff2(:,1), ff2, BaderThreshold)
        call parse_logical(datasets_check('LocalMultipolesExtraWrite'), .false., extra_write)
        call bader_union_inside(basins, lcl%nd, lcl%domain, lcl%lab, lcl%dshape, lcl%inside) 
        if (extra_write) then
          filename = 'basinsmap'
          iunit = io_open(file=trim(filename), action='write')
          do ip = 1, sys%gr%mesh%np
            write(iunit, '(3(f21.12,1x),i9)') (units_from_atomic(units_out%length,sys%gr%mesh%x(ip,ix)),ix=1,3), &
             basins%map(ip)
          end do
          call io_close(iunit)
        end if
        call local_center_of_mass(lcl%nd, lcl%domain, sys%geo, lcl%dcm)
        SAFE_DEALLOCATE_A(ff2)
      else
        do id = 1, lcl%nd
          call box_union_inside_vec(lcl%domain(id), sys%gr%mesh%np, sys%gr%mesh%x, lcl%inside(:,id))
        end do
        call local_center_of_mass(lcl%nd, lcl%domain, sys%geo, lcl%dcm)
      end if
    end if 
    
    POP_SUB(local_inside_domain)

  end subroutine local_inside_domain

  ! ---------------------------------------------------------
  subroutine bader_union_inside(basins, nd, dom, lab, dsh, inside)
    type(basins_t),    intent(in)  :: basins
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    character(len=15), intent(in)  :: lab(:)
    integer,           intent(in)  :: dsh(:)
    logical,           intent(out) :: inside(:,:)

    integer           :: ib, id, ip, ix, nb, rankmin
    integer           :: max_check
    integer, allocatable :: dunit(:), domain_map(:,:)
    FLOAT             :: dmin
    FLOAT, allocatable :: xi(:)
    logical           :: extra_write
    character(len=64) :: filename

    PUSH_SUB(bader_union_inside)

    SAFE_ALLOCATE(xi(1:sys%space%dim))
    inside = .false.

    do id = 1, nd
      if( dsh(id) /= BADER ) then
        call box_union_inside_vec(dom(id), sys%gr%mesh%np, sys%gr%mesh%x, inside(:,id))
      else
        nb = box_union_get_nboxes(dom(id))
        max_check = nb
        SAFE_ALLOCATE(domain_map(nd,max_check))
        do ib = 1, nb
          xi = box_union_get_center(dom(id), ib)
          ix = mesh_nearest_point(sys%gr%mesh, xi, dmin, rankmin)
          domain_map(id,ib) = basins%map(ix)
        end do
        do ip = 1, sys%gr%mesh%np
          if(any( domain_map(id,:) == basins%map(ip) ) ) inside(ip,id) = .true.
        end do
        SAFE_DEALLOCATE_A(domain_map)
      end if
    end do

    !%Variable LocalMultipolesExtraWrite
    !%Type logical
    !%Default false
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Writes additional information to files, when computing local multipoles. For 
    !% example, it writes coordinates of each local domain.
    !%End
    call parse_logical(datasets_check('LocalMultipolesExtraWrite'), .false., extra_write)

    if (extra_write) then
      SAFE_ALLOCATE(dunit(1:nd))
      do ip = 1, sys%gr%mesh%np
        do id = 1, nd
          if(ip == 1)then
            write(filename,'(a,a,a)')'local.general/domain.',trim(lab(id)),'.xyz'
            dunit(id) = io_open(file=trim(filename), action='write')
          end if
          if (inside(ip,id)) then
            write(dunit(id),*)(units_from_atomic(units_out%length, sys%gr%mesh%x(ip,ib)),ib=1,3)
          end if
        end do
      end do
      
      do id = 1, nd
        call io_close(dunit(id))
      end do
    end if

    SAFE_DEALLOCATE_A(xi)
    SAFE_DEALLOCATE_A(dunit)

    POP_SUB(bader_union_inside)

  end subroutine bader_union_inside

  ! ---------------------------------------------------------
  subroutine add_dens_to_ion_x(ff, geo)
    FLOAT,              intent(inout)   :: ff(:,:)
    type(geometry_t),   intent(inout)   :: geo

    integer :: ia, ix, rankmin
    FLOAT   :: dmin

    PUSH_SUB(add_dens_to_ion_x)

    do ia = 1, geo%natoms
      ix = mesh_nearest_point(sys%gr%mesh, geo%atom(ia)%x, dmin, rankmin)
      ff(ix,1) = ff(ix,1) + species_z(geo%atom(ia)%spec)
    end do

    POP_SUB(add_dens_to_ion_x)
  end subroutine add_dens_to_ion_x

  ! ---------------------------------------------------------
  subroutine local_center_of_mass(nd, dom, geo, center)
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(out) :: center(:,:)

    integer            :: ia, id
    FLOAT, allocatable :: sumw(:)

    PUSH_SUB(local_center_of_mass)

    center(:,:) = M_ZERO
    SAFE_ALLOCATE(sumw(1:nd)); sumw(:) = M_ZERO
    do ia = 1, geo%natoms
      do  id = 1, nd
        if (box_union_inside(dom(id),geo%atom(ia)%x)) then
          center(1:geo%space%dim,id) = center(1:geo%space%dim,id) &
                   + geo%atom(ia)%x(1:geo%space%dim)*species_weight(geo%atom(ia)%spec)     
          sumw(id) = sumw(id) + species_weight(geo%atom(ia)%spec)
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
