!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module io_function_m
  use cube_function_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use index_m
  use loct_parser_m
  use mesh_m
  use messages_m
  use profiling_m
  use simul_box_m
  use units_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif
  use par_vec_m
  use mpi_m
  use mpi_debug_m
  use varinfo_m

  implicit none

  private
  public ::                       &
    io_function_read_how,         &
    io_function_fill_how,         &
    write_xsf_geometry,           &
    write_xsf_geometry_file,      &
    dinput_function,              &
    zinput_function,              &
    doutput_function,             &
    zoutput_function,             &
    dio_function_out_text,        &
    zio_function_out_text,        &
    io_output_tensor,             &
    io_output_dipole

  integer, parameter, public ::   &
    output_axis_x     =     1,    &
    output_axis_y     =     2,    &
    output_axis_z     =     4,    &
    output_plane_x    =     8,    &
    output_plane_y    =    16,    &
    output_plane_z    =    32,    &
    output_dx         =    64,    &
    output_netcdf     =   128,    &
    output_plain      =   256,    &
    output_mesh_index =   512,    &
    output_xcrysden   =  1024,    &
    output_matlab     =  2048,    &
    output_meshgrid   =  4096,    &
    boundary_points   =  8192,    &
    output_binary     = 16384,    &
    output_etsf       = 32768

  ! doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private ::  &
    doutput_kind      =    1,     &
    zoutput_kind      =   -1

  ! index to axis mapping
  character, parameter, public :: &
    index2axis(3) = (/ 'x', 'y', 'z' /)

  ! index to label mapping
  character(len=3), parameter ::  &
    index2label(3) = (/ 're ', 'im ', 'abs' /)

  type(profile_t), save :: read_prof, write_prof

contains

  ! ---------------------------------------------------------
  subroutine io_function_read_how(sb, how)
    type(simul_box_t), intent(in)  :: sb
    integer,           intent(out) :: how

    call push_sub('out.output_init')

    how = 0

    !%Variable OutputHow
    !%Type flag
    !%Default 0
    !%Section Output
    !%Description
    !% Describes the format of the output files (see <tt>Output</tt>).
    !% Example: "axis_x + plane_x + dx"
    !%Option axis_x 1
    !% The values of the function on the <math>x</math> axis are printed. The string ".y=0,z=0" is appended
    !% to previous file names.
    !%Option axis_y 2
    !% The values of the function on the <math>y</math> axis are printed. The string ".x=0,z=0" is appended
    !% to previous file names.
    !%Option axis_z 4
    !% The values of the function on the <math>z</math> axis are printed. The string ".x=0,y=0" is appended
    !% to previous file names.
    !%Option plane_x 8
    !% A plane slice at <math>x=0</math> is printed. The string ".x=0'' is appended
    !% to previous file names.
    !%Option plane_y 16
    !% A plane slice at <math>y=0</math> is printed. The string ".y=0" is appended
    !% to previous file names.
    !%Option plane_z 32
    !% A plane slice at <math>y=0</math> is printed. The string ".z=0" is appended to
    !% previous file names.
    !%Option dx 64
    !% For printing three-dimensional information, the open-source program
    !% visualization tool OpenDX (http://www.opendx.org/) can be used. The string
    !% ".dx" is appended to previous file names.
    !%Option netcdf 128
    !% Outputs in NetCDF (http://www.unidata.ucar.edu/packages/netcdf/) format. This file
    !% can then be read, for example, by OpenDX. The string ".ncdf" is appended to previous file names.
    !% Requires the NetCDF library.
    !%Option plain 256
    !% Restart files are output in plain binary.
    !%Option mesh_index 512
    !% Generates output files of a given quantity (density, wfs, ...) which include
    !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
    !% useful for small meshes and debugging purposes.
    !% The output can also be used to display the mesh directly. A gnuplot script for mesh vizualization
    !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>
    !%Option xcrysden 1024
    !% A format for printing structures and three-dimensional information, which can be visualized by
    !% the open-source program XCrySDen (http://www.xcrysden.org/). The string
    !% ".xsf" is appended to previous file names. Note that lattice vectors and coordinates are as
    !% specified by UnitsOutput.
    !%Option matlab 2048
    !% In combination with plane_x, plane_y and plane_z this option produces output files 
    !% which are suitable for 2D Matlab functions like mesh(), surf() or waterfall(). To load 
    !% these files into Matlab you can use, e.g.
    !%
    !%   >> density = load('static/density-1.x=0.matlab.abs');
    !%   >> mesh(density);
    !%
    !%Option meshgrid 4096
    !% Outputs in Matlab mode the internal mesh in a format similar to e.g.
    !%
    !%   >> [x,y] = meshgrid(-2:.2:2,-1:.15:1)
    !%
    !% The x meshgrid is contained in a file *.meshgrid.x and the y-grid can be found in
    !% *.meshgrid.y
    !%Option boundary_points 8192
    !% This option includes the output of the mesh enlargement. Default is without.
    !%Option binary 16384
    !% Plain binary, new format.
    !%Option etsf 32768
    !% ETSF file format (http://www.etsf.eu/resources/software/standardization_project).
    !% Requires the ETSF_IO library.
    !%End
    call loct_parse_int(datasets_check('OutputHow'), 0, how)
    if(.not.varinfo_valid_option('OutputHow', how, is_flag=.true.)) then
      call input_error('OutputHow')
    end if

    if(how .eq. 0) then
      write(message(1), '(a)') 'Must specify output method with variable OutputHow.'
      call write_fatal(1)
     endif

    ! some modes are not available in some circumstances, so we reset how
    if(sb%dim == 1) how = iand(how, not(output_axis_y + output_plane_z))
    if(sb%dim <= 2) how = iand(how, not(output_axis_z + &
      output_plane_x + output_plane_y + output_dx))
#if !defined(HAVE_NETCDF)
    how = iand(how, not(output_netcdf))
#endif
#if !defined(HAVE_ETSF_IO)
    how = iand(how, not(output_etsf))
#endif

    call pop_sub()
  end subroutine io_function_read_how

  ! -------------------------------------------------------------------
  ! Use this function to quickly plot functions for debugging purposes:
  ! call doutput_function(io_function_fill_how("AxisX_and_PlaneX_and_DX"), &
  !                       ".", "func", m, sb, func, M_ONE, ierr)
  ! -------------------------------------------------------------------
  integer function io_function_fill_how(where) result(how)
    character(len=*), intent(in) :: where

    how = 0
    if(index(where, "AxisX").ne.0)     how = ior(how, output_axis_x)
    if(index(where, "AxisY").ne.0)     how = ior(how, output_axis_y)
    if(index(where, "AxisZ").ne.0)     how = ior(how, output_axis_z)
    if(index(where, "PlaneX").ne.0)    how = ior(how, output_plane_x)
    if(index(where, "PlaneY").ne.0)    how = ior(how, output_plane_y)
    if(index(where, "PlaneZ").ne.0)    how = ior(how, output_plane_z)
    if(index(where, "DX").ne.0)        how = ior(how, output_dx)
    if(index(where, "XCrySDen").ne.0)  how = ior(how, output_xcrysden)
    if(index(where, "Plain").ne.0)     how = ior(how, output_plain)
    if(index(where, "Binary").ne.0)    how = ior(how, output_binary)
    if(index(where, "MeshIndex").ne.0) how = ior(how, output_mesh_index)
#if defined(HAVE_NETCDF)
    if(index(where, "NETCDF").ne.0)    how = ior(how, output_netcdf)
#endif

  end function io_function_fill_how

  ! ---------------------------------------------------------
  subroutine write_xsf_geometry_file(dir, fname, geo, sb, offset)
    character(len=*),   intent(in) :: dir, fname
    type(geometry_t),   intent(in) :: geo
    type(simul_box_t),  intent(in) :: sb
    FLOAT,    optional, intent(in) :: offset(:)

    integer iunit
    character(len=6) position
    FLOAT, allocatable:: offset_(:)

    call push_sub('io_function.write_xsf_geometry_file')

    if( .not. mpi_grp_is_root(mpi_world)) return

    call io_mkdir(dir)
    position = 'asis'
    iunit = io_open(trim(dir)//'/'//trim(fname)//'.xsf', action='write', position=position)

    SAFE_ALLOCATE(offset_(1:sb%dim))
    if(.not. present(offset)) then
      offset_(1:sb%dim) = 0
    else
      ASSERT(ubound(offset, 1) >= sb%dim)
      offset_(1:sb%dim) = offset(1:sb%dim)
    endif

    call write_xsf_geometry(iunit, geo, sb, offset_(:))

    call io_close(iunit)

    call pop_sub()
  end subroutine write_xsf_geometry_file

  ! ---------------------------------------------------------
! for format specification see:
! http://www.xcrysden.org/doc/XSF.html#__toc__11
  subroutine write_xsf_geometry(iunit, geo, sb, offset)
    integer,            intent(in) :: iunit
    type(geometry_t),   intent(in) :: geo
    type(simul_box_t),  intent(in) :: sb
    FLOAT,    optional, intent(in) :: offset(:)

    integer idir, iatom

    call push_sub('io_function.write_xsf_geometry')

    select case(sb%periodic_dim)
      case(3)
        write(iunit, '(a)') 'CRYSTAL'
      case(2)
        write(iunit, '(a)') 'SLAB'
      case(1)
        write(iunit, '(a)') 'WIRE'  ! could also be POLYMER
      case(0)
        write(iunit, '(a)') 'MOLECULE'
    end select

    write(iunit, '(a)') 'PRIMVEC'

    do idir = 1, sb%dim
      write(iunit, '(3f12.6)') 2 * sb%lsize(idir) / units_out%length%factor * sb%rlattice(1:3, idir)
    enddo

    write(iunit, '(a)') 'PRIMCOORD'

    ! BoxOffset should be considered here
    write(iunit, '(i10, a)') geo%natoms, ' 1'
    do iatom = 1, geo%natoms
      write(iunit, '(a10, 3f12.6)') trim(geo%atom(iatom)%label), &
        ((geo%atom(iatom)%x(1:sb%dim) - offset(1:sb%dim)) / units_out%length%factor)
    enddo

    call pop_sub()
  end subroutine write_xsf_geometry


#if defined(HAVE_NETCDF)
  ! ---------------------------------------------------------
  subroutine ncdf_error(func, status, filename, ierr)
    character(len=*), intent(in)    :: func
    integer,          intent(in)    :: status
    character(len=*), intent(in)    :: filename
    integer,          intent(inout) :: ierr
    if(status .eq. NF90_NOERR) return
    write(message(1),'(3a)') "NETCDF error in function '" , trim(func) , "'"
    write(message(2),'(3a)') "(reading/writing ", trim(filename) , ")"
    write(message(3), '(6x,a,a)')'Error code = ', trim(nf90_strerror(status))
    call write_warning(3)
    ierr = 5
  end subroutine ncdf_error
#endif

  subroutine transpose3(in, out)
    FLOAT, intent(in)  :: in(:, :, :)
    FLOAT, intent(out) :: out(:, :, :)
    integer :: ix, iy, iz

    do ix = lbound(in, 1), ubound(in, 1)
      do iy = lbound(in, 2), ubound(in, 2)
        do iz = lbound(in, 3), ubound(in, 3)
          out(iz, iy, ix) = in(ix, iy, iz)
        end do
      end do
    end do
  end subroutine transpose3


  ! ---------------------------------------------------------
  subroutine io_output_tensor(iunit, tensor, ndim, factor)
    integer, intent(in) :: iunit
    FLOAT,   intent(in) :: tensor(:,:)
    integer, intent(in) :: ndim
    FLOAT,   intent(in) :: factor
    
    FLOAT :: trace
    integer :: j

    trace = M_z0
    do j = 1, ndim
      write(iunit, '(3f20.6)') tensor(j, 1:ndim)/factor
      trace = trace + tensor(j, j)
    end do
    trace = trace/TOFLOAT(ndim)

    write(iunit, '(a, f20.6)')  'Isotropic average', trace/factor
      
  end subroutine io_output_tensor


  ! ---------------------------------------------------------
  subroutine io_output_dipole(iunit, dipole, ndim) 
    integer, intent(in) :: iunit
    FLOAT,   intent(in) :: dipole(:)
    integer, intent(in) :: ndim
    
    integer :: idir
    FLOAT, parameter :: ATOMIC_TO_DEBYE = CNST(2.5417462)

    write(iunit, '(3a)') 'Dipole [', trim(units_out%length%abbrev), ']:                    [Debye]'
    do idir = 1, ndim
      write(iunit, '(6x,a,i1,a,es14.5,3x,2es14.5)') '<x', idir, '> = ', &
        dipole(idir) / units_out%length%factor, dipole(idir) * ATOMIC_TO_DEBYE
    end do

  end subroutine io_output_dipole


  ! ---------------------------------------------------------
#include "undef.F90"
#include "real.F90"
#include "io_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "io_function_inc.F90"
#include "undef.F90"

end module io_function_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
