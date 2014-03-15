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
  use batch_m
  use box_m
  use box_union_m
  use datasets_m
  use command_line_m
  use geometry_m
  use global_m
  use io_m
  use io_binary_m
  use io_function_m
  use lalg_adv_m
  use loct_m
  use messages_m
  use mesh_function_m
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
  
  integer           :: ierr
  character*256     :: config_str
  type(system_t)    :: sys
  type(simul_box_t) :: sb

  ! Initialize stuff
  call global_init(is_serial = .true.)

  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(config_str)
  call getopt_end()

  call messages_init()

  call datasets_init(1)
  call io_init()

  call unit_system_init()
  call system_init(sys)
  call simul_box_init(sb, sys%geo, sys%space)

  call local_domains()

  call simul_box_end(sb)
  call geometry_end(sys%geo)
  call space_end(sys%space)
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
    type(box_union_t), allocatable :: domain(:)
    integer                        :: err, id, nd, lmax, iter
    FLOAT                          :: default_dt, dt, spacing
    FLOAT, allocatable             :: read_ff(:)
    character(64)                  :: filename, folder, folder_default, aux
    character(len=15), allocatable :: lab(:)
    logical                        :: wrt_multipoles, wrt_local_densities

    PUSH_SUB(local_domains)

    message(1) = 'Info: Creating local domains'
    message(2) = ''
    call messages_info(2)

    write(folder_default,'(a)')'restart'

    !%Variable GlobalDensityFolder
    !%Type string
    !%Section Utilities::oct-convert
    !%Description
    !% The folder name where the input files are.
    !%End
    call parse_string(datasets_check('GlobalDensityFolder'), folder_default, folder)
    aux = trim(folder(1:3))
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
    
    SAFE_ALLOCATE(read_ff(1:sys%gr%mesh%np)); read_ff(:) = M_ZERO
    call drestart_read_function(folder, filename, sys%gr%mesh, read_ff, err)
    if (err /= 0 ) then
     write(message(1),*) 'While reading density: "',trim(folder),trim(filename),'", error code:',err
     call messages_fatal(1)
    end if
    
    !%Variable ConvertMultipole
    !%Type logical
    !%Default false
    !%Section Utilities::oct-convert
    !%Description
    !% This variable decides if a folder is going to be iterated or the
    !% filename is going to be iterated.
    !%End
    call parse_logical(datasets_check('LocalMultipoles'), .false., wrt_multipoles)

    !%Variable LocalMultipoleLmax 
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable means the maximum electric multipole of the density output.
    !%End
    call parse_integer(datasets_check('LocalMultipoleLmax '), 1, lmax)

    !%Variable domain
    !%Type union_box_t
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% This variable stores the information about the local boxes 
    !% for each new local domain. Information is read in a block format.
    !% <tt>
    !% %LocalDomains
    !% 'Label' | Shape | rsize  %< | Shape dependencies >%
    !% Label = string with the name of the new local domain.
    !% Shape = SPHERE, CYLINDER, PARALLELEPIPED, MINIMUM
    !% Shape dependencies:
    !% case(SPHERE):         | rsize | %<dim origin coordinates>
    !% case(CYLINDER):       | rsize | xsize | %<origin coordinates>
    !% case(PARALLELEPIPED): | %<lsize> | %<origin coordinates>
    !% case(MINIMUM):        | 'center_list' 
    !% rsize < Radius in input length units
    !% xsize < the length of the cylinder in the x-direction 
    !% origin coordinates < in input length units separated by | . where is the box centered.
    !% lsize <  half of the length of the parallelepiped in each direction.
    !% center_list < string containing the list of atoms in xyz file for each domain in the form "2,16-23"
    !% </tt>
    !%End
    call local_domains_read(domain, nd, lab)

    if ( wrt_multipoles ) then
      message(1) = 'Info: Computing local multipoles'
      message(2) = ''
      call messages_info(2)
!      call calc_local_multipoles(sys, nd, domain, lab, lmax, read_ff, sys%gr%mesh%np, iter, dt) 
      call calc_local_multipoles(nd, domain, lab, lmax, read_ff, sys%gr%mesh%np, iter, dt) 
    end if
    SAFE_DEALLOCATE_A(lab)
    do id = 1, nd
      call box_union_end(domain(id))
    end do

    message(1) = 'Info: Exiting local domains'
    message(2) = ''
    call messages_info(2)
    POP_SUB(local_domains)
  end subroutine local_domains

  ! ---------------------------------------------------------
  !> Reads the information (from the input file) about a local_t variable, initializing
  !! part of it (it has to be completed later with "local_init").
  ! ---------------------------------------------------------
  subroutine local_domains_read(domain, ndomain, lab)
    type(box_union_t), allocatable, intent(out) :: domain(:)
    integer,                        intent(out) :: ndomain
    character(len=15), allocatable, intent(out) :: lab(:)

    integer           :: ib, id
    type(block_t)     :: blk

    PUSH_SUB(local_domains_read)

    !%Variable LocalDomain
    !%Type block
    !%Section Execution::Local_multipoles utility
    !%Description
    !% Specifies the shape and properties of each local box.
    !%End

    ! First, find out if there is a Species block.
    ndomain = 0
    if(parse_block(datasets_check('LocalDomains'), blk) == 0) then
      ndomain = parse_block_n(blk)
    end if
    SAFE_ALLOCATE(domain(1:ndomain))
    SAFE_ALLOCATE(lab(1:ndomain))

    block: do id = 1, ndomain
      call parse_block_string(blk, id-1, 0, lab(id))
      write(message(1),'(a,a)')' Reading Local Domain: ',trim(lab(id))
      call messages_info(1)
      call read_from_domain_block(blk, id-1, domain(id))
    end do block
    message(1) = ''
    call messages_info(1)

    POP_SUB(local_domains_read)
  end subroutine local_domains_read

  ! ---------------------------------------------------------
  subroutine read_from_domain_block(blk, row, dom)
    type(block_t),     intent(in)        :: blk
    integer,           intent(in)        :: row
    type(box_union_t), intent(inout)     :: dom
    
    integer           :: dim, ic, idir, nb, shape
    character(len=80) :: clist, default
    character(len=15) :: lab
    FLOAT             :: lgst, val, rsize, xsize
    FLOAT             :: center(MAX_DIM), lsize(MAX_DIM)
   
    PUSH_SUB(read_from_domain_block)

! > Initializing variables in dom
    shape = 1
    nb = 1
    rsize = -M_ONE
    xsize = M_ZERO
    lsize(:) = M_ZERO
    center(:) = M_ZERO
    dim = sys%space%dim
!
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
        call parse_block_string(blk, row, 3, clist)
        nb = 0
        do ic = 1, sys%geo%natoms
          if(loct_isinstringlist(ic, clist)) nb = nb + 1
        end do
        message(1) = 'Bader shape is not yet implemented'
        call messages_fatal(1)
    end select
      ! fill in lsize structure
      select case(shape)
      case(SPHERE)
        lsize(1:dim) = rsize
      case(CYLINDER)       
        lsize(1)     = xsize
        lsize(2:dim) = rsize
      case(MINIMUM)
        do idir = 1, dim
          lgst = M_ZERO; val = M_ZERO
          do ic = 1, sys%geo%natoms
            if(loct_isinstringlist(ic, clist)) val = abs(sys%geo%atom(ic)%x(idir))
            if(lgst <= val) lgst = val
          end do
          lsize(idir) =  lgst + rsize
        end do
      end select

    call local_domains_init(dom, dim, shape, center, rsize, xsize, lsize, nb, clist)

    POP_SUB(read_from_domain_block)

  end subroutine read_from_domain_block

  ! ---------------------------------------------------------
  subroutine local_domains_init(dom, dim, shape, center, rsize, xsize, lsize, nb, clist)
    type(box_union_t), intent(inout) :: dom
    integer,           intent(in)    :: dim
    integer,           intent(in)    :: shape
    FLOAT,             intent(in)    :: center(dim)
    FLOAT,             intent(in)    :: rsize
    FLOAT,             intent(in)    :: xsize
    FLOAT,             intent(in)    :: lsize(MAX_DIM)
    integer,           intent(in)    :: nb
    character(len=80), intent(out)    :: clist

    integer                  :: ia, ibox, ic, id, nboxes, bshape
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
    end select
    call box_union_init(dom, nb, boxes)

   !< Check for a conflict between box_unon and clist
    ic = 0
    if (shape == MINIMUM ) then
      do ia = 1, sys%geo%natoms
        if (box_union_inside(dom, sys%geo%atom(ia)%x).and. .not.loct_isinstringlist(ia, clist) ) then
          ic = ic + 1
          if( ic <= 20 ) write(message(ic),'(a,a,I0,a,a)')'Atom: ',trim(species_label(sys%geo%atom(ia)%spec)),ia, & 
                                 ' is inside the union box BUT not in list: ',trim(clist)
        end if
      end do
      if (ic > 0) then 
         if( ic > 20 ) ic = 20
         call messages_info(ic)
         message(1) = '!'
         message(2) = '! WARNING: THIS FACT COULD GIVE INCORRECT RESULTS '
         message(3) = '! '
         call messages_info(3)
         if ( ic == 20 ) then
           message(1) = 'WARNING: AT LIST 19 ATOMS ARE NOT PRESENT IN THE LIST'
           call messages_info(1)
         end if
         message(1) = ''
           call messages_info(1)
      end if
    end if

    do ibox = 1, nb 
      call box_end(boxes(ibox))
    end do
    SAFE_DEALLOCATE_A(boxes)

    POP_SUB(local_domains_init)
  end subroutine local_domains_init

  ! ---------------------------------------------------------
  !> Computes the local multipoles and write them on a files.  
  ! ---------------------------------------------------------
  subroutine calc_local_multipoles(nd, dom, lab, lmax, ff, np, iter, dt)
    integer,                intent(in) :: nd
    type(box_union_t),      intent(in) :: dom(:)
    character(len=15),      intent(in) :: lab(:)
    integer,                intent(in) :: lmax
    FLOAT,                  intent(in) :: ff(1:np)
    integer,                intent(in) :: np
    integer,                intent(in) :: iter
    FLOAT,                  intent(in) :: dt   

    FLOAT, allocatable :: multipoles(:,:), ion_dipole(:,:), dcenter(:,:)
    integer            :: id, is, nspin

    PUSH_SUB(calc_local_multipoles)

    SAFE_ALLOCATE(multipoles(1:(lmax + 1)**2, nd)); multipoles(:,:) = M_ZERO

!> TODO: For instance spin are not included and just work with real densities.

    nspin = sys%st%d%nspin
    call local_center_of_mass(nd, dom, sys%geo, dcenter)
    call local_geometry_dipole(nd, dom, sys%geo, ion_dipole)
    call dmf_local_multipoles(sys%gr%mesh, nd, dom, ff, lmax, multipoles)
    do id = 1, nd
      multipoles(2:sys%space%dim+1, id) = -ion_dipole(1:sys%space%dim, id)/nspin - multipoles(2:sys%space%dim+1, id)
      call wrt_local_multipoles(multipoles(:,id), dcenter(:,id), lmax, sys%st%d%nspin, lab(id), iter, dt, &
                                 local_geometry_charge(dom(id), sys%geo))
    end do
    SAFE_DEALLOCATE_A(dcenter)
    SAFE_DEALLOCATE_A(ion_dipole)
    SAFE_DEALLOCATE_A(multipoles)

    POP_SUB(calc_local_multipoles)
  end subroutine calc_local_multipoles

  ! ---------------------------------------------------------
  subroutine local_center_of_mass(nd, dom, geo, center)
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,allocatable, intent(out) :: center(:,:)

    integer            :: ia, ibox, id
    FLOAT, allocatable :: sumw(:)

    PUSH_SUB(local_center_of_mass)

    SAFE_ALLOCATE(center(1:sys%space%dim, nd)); center(:,:) = M_ZERO
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

  ! ---------------------------------------------------------
  subroutine local_geometry_dipole(nd, dom, geo, dipole)
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,allocatable, intent(out) :: dipole(:,:)

    integer :: ia, ibox, id

    PUSH_SUB(local_geometry_dipole)

    SAFE_ALLOCATE(dipole(1:sys%space%dim, nd))
    dipole(:,:) = M_ZERO
    do ia = 1, geo%natoms
      do  id = 1, nd
        if (box_union_inside(dom(id), geo%atom(ia)%x)) then
          dipole(1:geo%space%dim, id) = dipole(1:geo%space%dim, id) + &
          species_zval(geo%atom(ia)%spec)*(geo%atom(ia)%x(1:geo%space%dim))
        end if
      end do
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(local_geometry_dipole)
  end subroutine local_geometry_dipole

  ! ---------------------------------------------------------
  real*8 function local_geometry_charge(dom, geo) result(charge)
    type(box_union_t), intent(in)  :: dom
    type(geometry_t),  intent(in)  :: geo

    integer :: ia, ibox, id

    PUSH_SUB(local_geometry_charge)

    charge = M_ZERO
    do ia = 1, geo%natoms
      if (box_union_inside(dom, geo%atom(ia)%x)) then
        charge = charge + species_zval(geo%atom(ia)%spec)
      end if
    end do

    POP_SUB(local_geometry_charge)
  end function local_geometry_charge

  ! ---------------------------------------------------------
  subroutine wrt_local_multipoles(multipoles, center, lmax, nspin, label, iter, dt, loc_charge)
    FLOAT,         intent(in) :: multipoles(:)
    FLOAT,         intent(in) :: center(:)
    integer,       intent(in) :: lmax
    integer,       intent(in) :: nspin
    character(15), intent(in) :: label
    integer,       intent(in) :: iter
    FLOAT,         intent(in) :: dt   
    FLOAT,         intent(in) :: loc_charge
   
    integer             :: add_lm, ll, mm, out_multipoles
    character(len=64)   :: filename
    character(len=21)   :: aux
    logical             :: file_exists
    
    PUSH_SUB(wrt_local_multipoles)

    call io_mkdir('local.multipoles')
    write(filename,'(a,a,a)')'local.multipoles/',trim(label),'.multipoles'
    inquire(file=filename, exist=file_exists)
    if (file_exists) then 
      out_multipoles = io_open(file=trim(filename), action='write', position='append', status='old')
    else
      out_multipoles = io_open(file=trim(filename), action='write', status='new')
      call write_multipoles_header(out_multipoles,nspin,lmax,center)
    end if
    write(aux,'(I10)') iter
    write(out_multipoles,'(a10)', advance='no')aux
    write(aux,'(es21.12)') iter*(units_from_atomic(units_inp%time, dt))
    write(out_multipoles,'(a)', advance='no')aux
    add_lm = 1
    do ll = 0, lmax
      do mm = -ll, ll
        write(aux,'(es21.12)')units_from_atomic(units_out%length**ll, multipoles(add_lm))
        write(out_multipoles,'(a)', advance='no')aux
        add_lm = add_lm + 1
      end do
    end do
    write(out_multipoles,'(es21.12)')loc_charge
!    write(out_multipoles,'(a)')""
    
    if(iand(sys%outp%how, C_OUTPUT_HOW_BILD) /= 0 .AND. sys%space%dim == 3)then
      call out_bld_multipoles(multipoles(2:sys%space%dim+1), center, lmax, label, iter)
    else
      message(1) = "Error. Output Format not available."
      message(2) = ''
      call messages_info(2)
    end if  
    
    call io_close(out_multipoles)

    POP_SUB(wrt_local_multipoles)
  end subroutine wrt_local_multipoles

  ! ---------------------------------------------------------
  subroutine out_bld_multipoles(multipoles, center, lmax, label, iter)
    FLOAT,         intent(in) :: multipoles(:)
    FLOAT,         intent(in) :: center(:)
    integer,       intent(in) :: lmax
    character(15), intent(in) :: label
    integer,       intent(in) :: iter
   
    integer             :: add_lm, ll, mm, out_bld
    character(len=80)   :: filename, folder
    FLOAT               :: dipolearrow(3,2)

    PUSH_SUB(out_bld_multipoles)
    
    write(folder,'(a,a)')'local.multipoles/',trim(label)
    call io_mkdir(folder)
    write(filename,'(a,a,a,a,i7.7,a)')trim(folder),'/',trim(label),'.',iter,'.bld'
!    out_bld = io_open(file=trim(filename), action='write', status='new')
    out_bld = io_open(file=trim(filename), action='write')

    write(out_bld,'(a,a,a,i7)')'.comment ** Arrow for the dipole moment centered at the center of mass for ', &
                        trim(label), ' domain and iteration number: ',iter
    write(out_bld,'(a)')''
    write(out_bld,'(a)')'.color red'
    write(out_bld,'(a,3(f12.6,2x),a)')'.sphere ',(units_from_atomic(units_out%length,center(ll)), ll= 1, 3),' 0.2' 
    do ll = 1, 3
      dipolearrow(ll,1) = units_from_atomic(units_out%length, center(ll) - multipoles(ll))
      dipolearrow(ll,2) = units_from_atomic(units_out%length, center(ll) + multipoles(ll))
    end do
    write(out_bld,'(a,6(f12.6,2x),a)')'.arrow ',(dipolearrow(ll,1), ll= 1, 3), &
                                     (dipolearrow(ll,2), ll= 1, 3), ' 0.1 0.5 0.90'
    POP_SUB(out_bld_multipoles)
  end subroutine out_bld_multipoles

  ! ---------------------------------------------------------
  !> Write header of multipoles files
  ! ---------------------------------------------------------
  subroutine write_multipoles_header(iunit, nspin, lmax, center)
    integer,           intent(in)    :: iunit
    integer,           intent(in)    :: nspin
    integer,           intent(in)    :: lmax
    FLOAT,             intent(in)    :: center(:)

    integer        :: is, ll, mm
    character(64)  :: space,frmt
    character(21)  :: aux

    PUSH_SUB(write_multipoles_header)

    call messages_print_stress(iunit)
      write(iunit, '(a)', advance='no')'# Center of mass: ('
    do ll = 1, sys%space%dim-1
      write(iunit, '(f21.12,a2)', advance='no') units_from_atomic(units_out%length,center(ll)),' ,'
    enddo
      write(iunit, '(f21.12,a2)', advance='yes') units_from_atomic(units_out%length,center(sys%space%dim)),' )'
    write(aux,'(a)')"# Iter";write(iunit,'(a20)',advance='no')aux
    write(aux,'(a)')"t";write(iunit,'(a12)',advance='no')aux!;write(iunit,'(a3)',advance='no')" "
    do is=1,nspin
      write(aux,'(a19,i1,a1)')"Electronic charges(",is,")";write(iunit,'(a)',advance='no')aux; write(iunit,'(a6)',advance='no')" "
      write(aux,'(a4,i1,a1)')"<x>(",is,")";write(iunit,'(a)',advance='no')aux
      write(aux,'(a4,i1,a1)')"<y>(",is,")";write(iunit,'(a)',advance='no')aux
      write(aux,'(a4,i1,a1)')"<z>(",is,")";write(iunit,'(a)',advance='no')aux
      write(aux,'(a12,i1,a1)')"Ion charges(",is,")";write(iunit,'(a)',advance='no')aux; write(iunit,'(a6)',advance='no')" "
    end do
    write(iunit,'(a)',advance='yes')""
    write(aux,'(a)')"#[Iter n.]";write(iunit,'(a16)',advance='no')aux
    write(aux,'(a)')'[' // trim(units_abbrev(units_out%time)) // ']';write(iunit,'(a)',advance='no')aux
      do is = 1, nspin
        do ll = 0, lmax
          do mm = -ll, ll
            select case(ll)
            case(0)
              write(aux,'(a)') 'Electrons';write(iunit,'(a)', advance='no') aux; write(iunit,'(a2)',advance='no')" "
            case(1)
              write(aux,'(a)') '[' // trim(units_abbrev(units_out%length)) // ']'
              write(iunit,'(a)', advance='no') aux
            case default
              write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
              write(aux,'(a)')'[' // trim(aux) // ']'
              write(iunit,'(a)', advance='no') aux
            end select
          end do
        end do
      end do
    write(aux,'(a)') 'Electrons';write(iunit,'(a)', advance='no') aux; write(iunit,'(a2)',advance='no')" "
    write(iunit,'(a)')""
    call messages_print_stress(iunit)

  POP_SUB(write_multipoles_header)
   
  end subroutine write_multipoles_header
  ! ---------------------------------------------------------

end program oct_local_multipoles

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
