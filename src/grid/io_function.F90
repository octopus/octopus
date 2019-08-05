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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module io_function_oct_m
  use comm_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use distributed_oct_m
  use geometry_oct_m
  use global_oct_m
  use index_oct_m
  use io_oct_m
  use io_binary_oct_m
  use io_csv_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_debug_oct_m
  use namespace_oct_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use vtk_oct_m

  implicit none

  private
  public ::                       &
    io_function_read_how,         &
    io_function_fill_how,         &
    write_bild_forces_file,       &
    write_canonicalized_xyz_file, &
    write_xsf_geometry,           &
    write_xsf_geometry_file,      &
    dio_function_input,           &
    zio_function_input,           &
    dio_function_output,          &
    zio_function_output,          &
    io_function_output_vector,    &
    dio_function_output_global,   &
    zio_function_output_global,   &
    io_function_output_vector_BZ, &
    io_function_output_global_BZ

#if defined(HAVE_NETCDF)
 public ::                        &
    dout_cf_netcdf,               &
    zout_cf_netcdf
#endif

  !> doutput_kind => real variables; zoutput_kind => complex variables.
  integer, parameter, private ::  &
    DOUTPUT_KIND      =    1,     &
    ZOUTPUT_KIND      =   -1

  !> index to label mapping
  character(len=3), parameter ::  &
    index2label(3) = (/ 're ', 'im ', 'abs' /)

  type(profile_t), save :: read_prof, write_prof

  interface io_function_output_vector
    module procedure dio_function_output_vector, zio_function_output_vector
  end interface io_function_output_vector

  interface io_function_output_vector_BZ
    module procedure dio_function_output_vector_BZ, zio_function_output_vector_BZ
  end interface io_function_output_vector_BZ

  interface io_function_output_global_BZ
    module procedure dio_function_output_global_BZ, zio_function_output_global_BZ
  end interface io_function_output_global_BZ


contains

  ! ---------------------------------------------------------
  subroutine io_function_read_how(sb, namespace, how, ignore_error)
    type(simul_box_t), intent(in)  :: sb
    type(namespace_t), intent(in)  :: namespace
    integer(8),        intent(out) :: how
    logical, optional, intent(in)  :: ignore_error !> Ignore error check. Used when called from some external utility.

    PUSH_SUB(io_function_read_how)

    how = 0_8
    
    call messages_obsolete_variable(namespace, 'OutputHow', 'OutputFormat')
    
    !%Variable OutputFormat
    !%Type flag
    !%Default 0
    !%Section Output
    !%Description
    !% Describes the format of the output files (see <tt>Output</tt>).
    !% Example: <tt>axis_x + plane_x + dx</tt>
    !%Option axis_x bit(0) 
    !% The values of the function on the <i>x</i> axis are printed. The string <tt>.y=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_y bit(1)
    !% The values of the function on the <i>y</i> axis are printed. The string <tt>.x=0,z=0</tt> is appended
    !% to previous file names.
    !%Option axis_z bit(2)
    !% The values of the function on the <i>z</i> axis are printed. The string <tt>.x=0,y=0</tt> is appended
    !% to previous file names.
    !%Option plane_x bit(3)
    !% A plane slice at <i>x</i> = 0 is printed. The string <tt>.x=0</tt> is appended
    !% to previous file names.
    !%Option plane_y bit(4)
    !% A plane slice at <i>y</i> = 0 is printed. The string <tt>.y=0</tt> is appended
    !% to previous file names.
    !%Option plane_z bit(5)
    !% A plane slice at <i>z</i> = 0 is printed. The string <tt>.z=0</tt> is appended to
    !% previous file names.
    !%Option dx bit(6)
    !% For printing three-dimensional information, the open-source program
    !% visualization tool <a href=http://www.opendx.org>OpenDX</a> can be used. The string
    !% <tt>.dx</tt> is appended to previous file names. Available only in 3D.
    !%Option netcdf bit(7)
    !% Outputs in <a href=http://www.unidata.ucar.edu/packages/netcdf>NetCDF</a> format. This file
    !% can then be read, for example, by OpenDX. The string <tt>.ncdf</tt> is appended to previous file names.
    !% Requires the NetCDF library. Only writes the real part of complex functions.
    !%Option mesh_index bit(8)
    !% Generates output files of a given quantity (density, wavefunctions, ...) which include
    !% the internal numbering of mesh points. Since this mode produces large datafiles this is only 
    !% useful for small meshes and debugging purposes.
    !% The output can also be used to display the mesh directly. A Gnuplot script for mesh visualization
    !% can be found under <tt>PREFIX/share/octopus/util/display_mesh_index.gp</tt>.
    !%Option xcrysden bit(9)
    !% A format for printing structures and three-dimensional information, which can be visualized by
    !% the free open-source program <a href=http://www.xcrysden.org>XCrySDen</a> and others. The string
    !% <tt>.xsf</tt> is appended to previous file names. Note that lattice vectors and coordinates are as
    !% specified by <tt>UnitsOutput</tt>. Available in 2D and 3D.
    !%Option matlab bit(10)
    !% In combination with <tt>plane_x</tt>, <tt>plane_y</tt> and
    !% <tt>plane_z</tt>, this option produces output files which are
    !% suitable for 2D Matlab functions like <tt>mesh()</tt>,
    !% <tt>surf()</tt>, or <tt>waterfall()</tt>. To load these files
    !% into Matlab you can use, <i>e.g.</i>
    !%<tt>
    !%   >> density = load('static/density-1.x=0.matlab.abs');
    !%   >> mesh(density);
    !%</tt>
    !%Option meshgrid bit(11)
    !% Outputs in Matlab mode the internal mesh in a format similar to
    !%<tt>
    !%   >> [x,y] = meshgrid(-2:.2:2,-1:.15:1)
    !%</tt>
    !% The <i>x</i> meshgrid is contained in a file <tt>*.meshgrid.x</tt> and the <i>y</i>-grid can be found in
    !% <tt>*.meshgrid.y</tt>.
    !%Option boundary_points bit(12)
    !% This option includes the output of the mesh enlargement. Default is without.
    !% Supported only by <tt>binary</tt>, <tt>axis</tt>, <tt>plane</tt>, <tt>mesh_index</tt>,
    !% and <tt>matlab</tt> formats.
    !% Not all types of <tt>Output</tt> will have this information available. Not supported when parallel in domains.
    !%Option binary bit(13)
    !% Plain binary, new format.
    !%Option etsf bit(14)
    !% <a href=http://www.etsf.eu/resources/software/standardization_project>ETSF file format</a>.
    !% Requires the ETSF_IO library. Applies only to <tt>Output = density</tt>, <tt>geometry</tt>,
    !% <tt>wfs</tt>, and/or <tt>wfs_fourier</tt>.
    !%Option xyz bit(15)
    !% Geometry will be output in XYZ format. Does not affect other outputs.
    !%Option cube bit(16)
    !% Generates output in the <a href=http://paulbourke.net/dataformats/cube>cube file format</a>.
    !% Available only in 3D. Only writes the real part of complex functions.
    !%Option bild bit(19)
    !% Generates output in <a href=http://plato.cgl.ucsf.edu/chimera/docs/UsersGuide/bild.html>BILD format</a>.
    !%Option vtk bit(20)
    !% Generates output in <a href=http://www.vtk.org/VTK/img/file-formats.pdf>VTK legacy format</a>.
    !%Option integrate_xy bit(21)
    !% Integrates the function in the x-y plane and the result on the <i>z</i> axis is printed.
    !%Option integrate_xz bit(22)
    !% Integrates the function in the x-z plane and the result on the <i>y</i> axis is printed
    !%Option integrate_yz bit(23)
    !% Integrates the function in the y-z plane and the result on the <i>x</i> axis is printed
    !%End
    call parse_variable(namespace, 'OutputFormat', 0, how)
    if(.not.varinfo_valid_option('OutputFormat', how, is_flag=.true.)) then
      call messages_input_error('OutputFormat')
    end if

    if(how  ==  0 .and. .not. optional_default(ignore_error, .false.)) then
      write(message(1), '(a)') 'Must specify output method with variable OutputFormat.'
      call messages_fatal(1, only_root_writes = .true.)
     end if

    ! some modes are not available in some circumstances
    if(sb%dim == 1) then
      if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Y) /= 0) then
        message(1) = "OutputFormat = axis_y not available with Dimensions = 1."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Z) /= 0) then
        message(1) = "OutputFormat = plane_z not available with Dimensions = 1."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__XCRYSDEN) /= 0) then
        message(1) = "OutputFormat = xcrysden not available with Dimensions = 1."
        call messages_fatal(1)
      end if
    end if

    if(sb%dim <= 2) then
      if(bitand(how, OPTION__OUTPUTFORMAT__AXIS_Z) /= 0) then
        message(1) = "OutputFormat = axis_z not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_X) /= 0) then
        message(1) = "OutputFormat = plane_x not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__PLANE_Y) /= 0) then
        message(1) = "OutputFormat = plane_y not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__INTEGRATE_XY) /= 0) then
        message(1) = "OutputFormat = integrate_xy not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__INTEGRATE_XZ) /= 0) then
        message(1) = "OutputFormat = integrate_xz not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__INTEGRATE_YZ) /= 0) then
        message(1) = "OutputFormat = integrate_yz not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__DX) /= 0) then
        message(1) = "OutputFormat = dx not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
      if(bitand(how, OPTION__OUTPUTFORMAT__CUBE) /= 0) then
        message(1) = "OutputFormat = cube not available with Dimensions <= 2."
        call messages_fatal(1)
      end if
    end if

#if !defined(HAVE_NETCDF)
    if (bitand(how, OPTION__OUTPUTFORMAT__NETCDF) /= 0) then
      message(1) = 'Octopus was compiled without NetCDF support.'
      message(2) = 'It is not possible to write output in NetCDF format.'
      call messages_fatal(2)
    end if
#endif
#if !defined(HAVE_ETSF_IO)
    if (bitand(how, OPTION__OUTPUTFORMAT__ETSF) /= 0) then
      message(1) = 'Octopus was compiled without ETSF_IO support.'
      message(2) = 'It is not possible to write output in ETSF format.'
      call messages_fatal(2)
    end if
#endif

    POP_SUB(io_function_read_how)
  end subroutine io_function_read_how

  ! -------------------------------------------------------------------
  !> Use this function to quickly plot functions for debugging purposes:
  !! call dio_function_output(io_function_fill_how("AxisX_and_PlaneX_and_DX"), &
  !                       ".", "func", mesh, sb, func, M_ONE, ierr)
  ! -------------------------------------------------------------------
  integer(8) function io_function_fill_how(where) result(how)
    character(len=*), intent(in) :: where

    PUSH_SUB(io_function_fill_how)

    how = 0
    if(index(where, "AxisX")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_X)
    if(index(where, "AxisY")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_Y)
    if(index(where, "AxisZ")     /= 0) how = ior(how, OPTION__OUTPUTFORMAT__AXIS_Z)
    IF(index(where, "PlaneX")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_X)
    if(index(where, "PlaneY")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Y)
    if(index(where, "PlaneZ")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Z)
    if(index(where, "IntegrateXY") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_XY)
    if(index(where, "IntegrateXZ") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_XZ)
    if(index(where, "IntegrateYZ") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__INTEGRATE_YZ)
    if(index(where, "PlaneZ")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__PLANE_Z)
    if(index(where, "DX")        /= 0) how = ior(how, OPTION__OUTPUTFORMAT__DX)
    if(index(where, "XCrySDen")  /= 0) how = ior(how, OPTION__OUTPUTFORMAT__XCRYSDEN)
    if(index(where, "Binary")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__BINARY)
    if(index(where, "MeshIndex") /= 0) how = ior(how, OPTION__OUTPUTFORMAT__MESH_INDEX)
    if(index(where, "XYZ")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__XYZ)
#if defined(HAVE_NETCDF)
    if(index(where, "NETCDF")    /= 0) how = ior(how, OPTION__OUTPUTFORMAT__NETCDF)
#endif
    if(index(where, "Cube")      /= 0) how = ior(how, OPTION__OUTPUTFORMAT__CUBE)
    if(index(where, "VTK")       /= 0) how = ior(how, OPTION__OUTPUTFORMAT__VTK)

    POP_SUB(io_function_fill_how)
  end function io_function_fill_how

  ! ---------------------------------------------------------
  subroutine write_bild_forces_file(dir, fname, geo, mesh)
    character(len=*),   intent(in) :: dir, fname
    type(geometry_t),   intent(in) :: geo
    type(mesh_t),       intent(in) :: mesh

    integer :: iunit, iatom, idir
    FLOAT, allocatable :: forces(:,:), center(:,:)
    character(len=20) frmt

    if( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_bild_forces_file)

    call io_mkdir_old(dir)
    iunit = io_open_old(trim(dir)//'/'//trim(fname)//'.bild', action='write', position='asis')

    write(frmt,'(a,i0,a)')'(a,2(', mesh%sb%dim,'f16.6,1x))'

    SAFE_ALLOCATE(forces(1:geo%natoms, 1:mesh%sb%dim))
    SAFE_ALLOCATE(center(1:geo%natoms, 1:mesh%sb%dim))
    forall(iatom = 1:geo%natoms, idir = 1:mesh%sb%dim)
      forces(iatom, idir) = units_from_atomic(units_out%force, geo%atom(iatom)%f(idir))
      center(iatom, idir) = units_from_atomic(units_out%length, geo%atom(iatom)%x(idir))
    end forall
    write(iunit, '(a)')'.comment : force vectors in ['//trim(units_abbrev(units_out%force))//']'
    write(iunit, *)
    write(iunit, '(a)')'.color red'
    write(iunit, *)
    do iatom = 1, geo%natoms
      write(iunit, '(a,1x,i4,1x,a2,1x,a6,1x,f10.6,a)')'.comment :', iatom, trim(geo%atom(iatom)%label), & 
                         'force:', sqrt(sum(forces(iatom,:)**2)),'['//trim(units_abbrev(units_out%force))//']'
      write(iunit,fmt=trim(frmt))'.arrow',(center(iatom, idir), idir = 1, mesh%sb%dim), &
                                 (center(iatom, idir) + forces(iatom, idir), idir = 1, mesh%sb%dim)
      write(iunit,*)
    end do

    SAFE_DEALLOCATE_A(forces)
    SAFE_DEALLOCATE_A(center)

    call io_close(iunit)

    POP_SUB(write_bild_forces_file)
  end subroutine write_bild_forces_file

  ! ---------------------------------------------------------
  !> Write canonicalized xyz file with atom labels and positions in Angstroms.
  !> Includes information about simulation box and periodicity when applicable.
  !> This differs from a normal xyz file by including information about box
  !> shape and always using Angstroms.
  subroutine write_canonicalized_xyz_file(dir, fname, geo, mesh)
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: fname
    type(geometry_t), intent(in) :: geo
    type(mesh_t),     intent(in) :: mesh

    integer :: iunit
    integer :: idir
    FLOAT :: position
    integer :: iatom

    PUSH_SUB(write_canonicalized_xyz_file)

    call io_mkdir_old(dir)
    iunit = io_open_old(trim(dir)//'/'//trim(fname)//'.xyz', action='write', position='asis')

    write(iunit, '(i6)') geo%natoms
    call simul_box_write_short_info(mesh%sb, iunit)

    ! xyz-style labels and positions:
    do iatom=1, geo%natoms
      write(iunit, '(10a)', advance='no') geo%atom(iatom)%label
      do idir=1, 3
        if(idir <= mesh%sb%dim) then
          position = geo%atom(iatom)%x(idir)
        else
          position = M_ZERO
        end if
        write(iunit, '(xf11.6)', advance='no') units_from_atomic(unit_angstrom, position)
      end do
      write(iunit, '()')
    end do

    call io_close(iunit)

    POP_SUB(write_canonicalized_xyz_file)
  end subroutine write_canonicalized_xyz_file

  ! ---------------------------------------------------------
  subroutine write_xsf_geometry_file(dir, fname, geo, mesh, write_forces)
    character(len=*),   intent(in) :: dir, fname
    type(geometry_t),   intent(in) :: geo
    type(mesh_t),       intent(in) :: mesh
    logical,  optional, intent(in) :: write_forces

    integer :: iunit, iatom, idir
    FLOAT, allocatable :: forces(:,:)
    logical :: write_forces_

    if( .not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(write_xsf_geometry_file)

    call io_mkdir_old(dir)
    iunit = io_open_old(trim(dir)//'/'//trim(fname)//'.xsf', action='write', position='asis')

    if(.not. present(write_forces)) then
      write_forces_ = .false.
    else
      write_forces_ = write_forces
    end if

    if(write_forces_) then
      SAFE_ALLOCATE(forces(1:geo%natoms, 1:mesh%sb%dim))
      forall(iatom = 1:geo%natoms, idir = 1:mesh%sb%dim)
        forces(iatom, idir) = units_from_atomic(units_out%force, geo%atom(iatom)%f(idir))
      end forall
      call write_xsf_geometry(iunit, geo, mesh, forces = forces)
      SAFE_DEALLOCATE_A(forces)
    else
      call write_xsf_geometry(iunit, geo, mesh)
    end if

    call io_close(iunit)

    POP_SUB(write_xsf_geometry_file)
  end subroutine write_xsf_geometry_file

  ! ---------------------------------------------------------
  !> for format specification see:
  !! http://www.xcrysden.org/doc/XSF.html#__toc__11
  subroutine write_xsf_geometry(iunit, geo, mesh, forces, index)
    integer,           intent(in) :: iunit
    type(geometry_t),  intent(in) :: geo
    type(mesh_t),      intent(in) :: mesh
    FLOAT,   optional, intent(in) :: forces(:, :)
    integer, optional, intent(in) :: index !< for use in writing animated files

    integer :: idir, idir2, iatom, index_
    character(len=7) :: index_str
    FLOAT :: offset(3)

    PUSH_SUB(write_xsf_geometry)

    if(present(index)) then
      write(index_str, '(a,i6)') ' ', index
      index_ = index
    else
      write(index_str, '(a)') ''
      index_ = 1
    end if

    offset = M_ZERO
    ! The corner of the cell is always (0,0,0) to XCrySDen
    ! so the offset is applied to the atomic coordinates.
    ! Offset in periodic directions:
    offset(1:3) = -matmul(mesh%sb%rlattice_primitive(1:3,1:3), mesh%sb%lsize(1:3))
    ! Offset in aperiodic directions:
    do idir = mesh%sb%periodic_dim + 1, 3
      offset(idir) = -(mesh%idx%ll(idir) - 1)/2 * mesh%spacing(idir)
    end do

    if(simul_box_is_periodic(mesh%sb)) then
      if(index_ == 1) then
        select case(mesh%sb%periodic_dim)
          case(3)
            write(iunit, '(a)') 'CRYSTAL'
          case(2)
            write(iunit, '(a)') 'SLAB'
          case(1)
            write(iunit, '(a)') 'POLYMER'
        end select
      end if

      write(iunit, '(a)') 'PRIMVEC'//trim(index_str)

      do idir = 1, mesh%sb%dim
        write(iunit, '(3f12.6)') (units_from_atomic(units_out%length, &
          mesh%sb%rlattice(idir2, idir)), idir2 = 1, mesh%sb%dim)
      end do

      write(iunit, '(a)') 'PRIMCOORD'//trim(index_str)
      write(iunit, '(i10, a)') geo%natoms, ' 1'
    else
      write(iunit, '(a)') 'ATOMS'//trim(index_str)
    end if

    ! BoxOffset should be considered here
    do iatom = 1, geo%natoms
      write(iunit, '(a10, 3f12.6)', advance='no') trim(geo%atom(iatom)%label), &
        (units_from_atomic(units_out%length, geo%atom(iatom)%x(idir) - offset(idir)), idir = 1, mesh%sb%dim)
      if(present(forces)) then
        write(iunit, '(5x, 3f12.6)', advance='no') (forces(iatom, idir), idir = 1, mesh%sb%dim)
      end if
      write(iunit, '()')
    end do

    POP_SUB(write_xsf_geometry)
  end subroutine write_xsf_geometry


#if defined(HAVE_NETCDF)
  ! ---------------------------------------------------------
  subroutine ncdf_error(func, status, filename, ierr)
    character(len=*), intent(in)    :: func
    integer,          intent(in)    :: status
    character(len=*), intent(in)    :: filename
    integer,          intent(inout) :: ierr

    PUSH_SUB(ncdf_error)

    if(status  ==  NF90_NOERR) then
    POP_SUB(ncdf_error)
      return
    end if

    write(message(1),'(3a)') "NETCDF error in function '" , trim(func) , "'"
    write(message(2),'(3a)') "(reading/writing ", trim(filename) , ")"
    write(message(3), '(6x,a,a)')'Error code = ', trim(nf90_strerror(status))
    call messages_warning(3)
    ierr = 5

    POP_SUB(ncdf_error)
  end subroutine ncdf_error
#endif

  ! ---------------------------------------------------------
  subroutine transpose3(in, out)
    FLOAT, intent(in)  :: in(:, :, :)
    FLOAT, intent(out) :: out(:, :, :)
    integer :: ix, iy, iz

    PUSH_SUB(transpose3)

    do ix = lbound(in, 1), ubound(in, 1)
      do iy = lbound(in, 2), ubound(in, 2)
        do iz = lbound(in, 3), ubound(in, 3)
          out(iz, iy, ix) = in(ix, iy, iz)
        end do
      end do
    end do

    POP_SUB(transpose3)
  end subroutine transpose3

#include "undef.F90"
#include "real.F90"
#include "io_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "io_function_inc.F90"
#include "undef.F90"

end module io_function_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
