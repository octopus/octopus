!! Copyright (C) 2009-2011 X. Andrade, M. Oliveira
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
!! $Id: output_etsf_inc.F90 5880 2009-09-03 23:44:44Z dstrubbe $

! ---------------------------------------------------------

subroutine output_etsf(st, gr, geo, dir, outp)
  type(states_t),         intent(in) :: st
  type(grid_t),           intent(in) :: gr
  type(geometry_t),       intent(in) :: geo
  character(len=*),       intent(in) :: dir
  type(output_t),         intent(in) :: outp

  type(cube_t) :: cube
  type(cube_function_t) :: cf
#ifdef HAVE_ETSF_IO
  type(fourier_shell_t) :: shell
  integer :: ncid
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  type(etsf_dims) :: geometry_dims, density_dims, wfs_dims, pw_dims
  type(etsf_groups_flags) :: geometry_flags, density_flags, wfs_flags, pw_flags
#endif

  PUSH_SUB(output_etsf)

#ifndef HAVE_ETSF_IO
  ASSERT(.false.)
#endif

  !Create a cube
  call cube_init(cube, gr%mesh%idx%ll, gr%sb, fft=.true.)
  call cube_function_null(cf)
  
  ! To create an etsf file one has to do the following:
  !
  ! * Calculate the dimensions and the flags with the _dims functions
  !   for all sets of values.
  ! * Init the data file with the flags and the dims.
  ! * Call the _write functions for all sets.
  ! * Close the file.
  !
  ! Note: to keep things clean, new data MUST be added following this
  ! scheme and using functions.


#ifdef HAVE_ETSF_IO

  ! geometry
  if (iand(outp%what, C_OUTPUT_GEOMETRY).ne.0) then

    call output_etsf_geometry_dims(geo, gr%sb, geometry_dims, geometry_flags)

    call output_etsf_file_init(dir//"/geometry-etsf.nc", "Crystallographic_data file", geometry_dims, geometry_flags, ncid)

    call output_etsf_geometry_write(geo, gr%sb, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
  end if

  ! density
  if (iand(outp%what, C_OUTPUT_DENSITY).ne.0) then
    call dcube_function_alloc_RS(cube, cf)

    call output_etsf_geometry_dims(geo, gr%sb, density_dims, density_flags)
    call output_etsf_density_dims(st, gr%mesh, cube, density_dims, density_flags)

    call output_etsf_file_init(dir//"/density-etsf.nc", "Density file", density_dims, density_flags, ncid)

    call output_etsf_density_write(st, gr%mesh, cube, cf, ncid)
    call output_etsf_geometry_write(geo, gr%sb, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call dcube_function_free_rs(cube, cf)
  end if

  ! wave-functions
  if (iand(outp%what, C_OUTPUT_WFS).ne.0) then
    call dcube_function_alloc_RS(cube, cf)

    call output_etsf_geometry_dims(geo, gr%sb, wfs_dims, wfs_flags)
    call output_etsf_kpoints_dims(gr%sb, wfs_dims, wfs_flags)
    call output_etsf_electrons_dims(st, wfs_dims, wfs_flags)
    call output_etsf_wfs_rsp_dims(st, gr%mesh, cube, wfs_dims, wfs_flags)

    call output_etsf_file_init(dir//"/wfs-etsf.nc", "Wavefunctions file", wfs_dims, wfs_flags, ncid)

    call output_etsf_electrons_write(st, ncid)
    call output_etsf_geometry_write(geo, gr%sb, ncid)
    call output_etsf_kpoints_write(gr%sb, ncid)
    call output_etsf_wfs_rsp_write(st, gr%mesh, cube, cf, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
    
    call dcube_function_free_rs(cube, cf)
  end if

  ! wave-functions in fourier space
  if (iand(outp%what, C_OUTPUT_WFS_FOURIER).ne.0) then
    call zcube_function_alloc_RS(cube, cf)
    call fourier_shell_init(shell, cube, gr%mesh)

    call output_etsf_geometry_dims(geo, gr%sb, pw_dims, pw_flags)
    call output_etsf_kpoints_dims(gr%sb, pw_dims, pw_flags)
    call output_etsf_electrons_dims(st, pw_dims, pw_flags)
    call output_etsf_basisdata_dims(st, gr%mesh, cube, shell, pw_dims, pw_flags)
    call output_etsf_wfs_pw_dims(st, gr%mesh, cube, shell, pw_dims, pw_flags)

    call output_etsf_file_init(dir//"/wfs-pw-etsf.nc", "Wavefunctions file", pw_dims, pw_flags, ncid)

    call output_etsf_electrons_write(st, ncid)
    call output_etsf_geometry_write(geo, gr%sb, ncid)
    call output_etsf_kpoints_write(gr%sb, ncid)
    call output_etsf_basisdata_write(st, gr%mesh, shell, ncid)
    call output_etsf_wfs_pw_write(st, gr%mesh, cube, cf, shell, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call fourier_shell_end(shell)
    call zcube_function_free_rs(cube, cf)
  end if
#endif

  call cube_end(cube)

  POP_SUB(output_etsf)
end subroutine output_etsf

#ifdef HAVE_ETSF_IO
! --------------------------------------------------------

subroutine output_etsf_file_init(filename, filetype, dims, flags, ncid)
  character(len=*),        intent(in)    :: filename
  character(len=*),        intent(in)    :: filetype
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  integer,                 intent(out)   :: ncid

  logical :: lstat
  type(etsf_io_low_error) :: error_data

  PUSH_SUB(output_etsf_file_init)

  call etsf_io_data_init(filename, flags, dims, filetype, "Created by "//PACKAGE_STRING, &
    lstat, error_data, overwrite = .true., k_dependent = .false.)
  if (.not. lstat) call output_etsf_error(error_data)

  call etsf_io_low_open_modify(ncid, filename, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  POP_SUB(output_etsf_file_init)
end subroutine output_etsf_file_init

! --------------------------------------------------------

subroutine output_etsf_error(error_data)
  type(etsf_io_low_error), intent(in) :: error_data

  PUSH_SUB(output_etsf_error)
  
  call etsf_io_low_error_handle(error_data)
  message(1) = "ETSF_IO returned a fatal error. See message above."
  call messages_fatal(1)

  POP_SUB(output_etsf_error)
end subroutine output_etsf_error

! --------------------------------------------------------

subroutine output_etsf_geometry_dims(geo, sb, dims, flags)
  type(geometry_t),        intent(in)    :: geo
  type(simul_box_t),       intent(in)    :: sb
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  
  PUSH_SUB(output_etsf_geometry_dims)

  flags%geometry = etsf_geometry_all - etsf_geometry_valence_charges - etsf_geometry_pseudo_types

  dims%number_of_symmetry_operations = symmetries_number(sb%symm)
  dims%number_of_atom_species = geo%nspecies
  dims%number_of_atoms = geo%natoms

  POP_SUB(output_etsf_geometry_dims)
end subroutine output_etsf_geometry_dims

! --------------------------------------------------------

subroutine output_etsf_geometry_write(geo, sb, ncid)
  type(geometry_t),       intent(in)    :: geo
  type(simul_box_t),      intent(in)    :: sb
  integer,                intent(in)    :: ncid

  type(etsf_geometry) :: geometry
  integer :: idir, isymm, ispecies, i, j
  FLOAT :: offset(1:3)
  type(etsf_io_low_error)  :: error_data
  logical :: lstat

  PUSH_SUB(output_etsf_geometry_write)

  ! Primitive vectors
  SAFE_ALLOCATE(geometry%primitive_vectors(1:3, 1:3))
  do idir = 1, sb%dim
    geometry%primitive_vectors(1:3, idir) = sb%rlattice(1:3, idir)
  end do

  ! The symmetries
  SAFE_ALLOCATE(geometry%space_group)
  geometry%space_group = symmetries_space_group_number(sb%symm)
  SAFE_ALLOCATE(geometry%reduced_symmetry_matrices(1:3, 1:3, 1:symmetries_number(sb%symm)))
  SAFE_ALLOCATE(geometry%reduced_symmetry_translations(1:3, 1:symmetries_number(sb%symm)))

  do isymm = 1, symmetries_number(sb%symm)
    geometry%reduced_symmetry_matrices(1:3, 1:3, isymm) = symm_op_rotation_matrix(sb%symm%ops(isymm))
    geometry%reduced_symmetry_translations(1:3, isymm) = symm_op_translation_vector(sb%symm%ops(isymm))
  end do

  ! The species
  SAFE_ALLOCATE(geometry%atomic_numbers(geo%nspecies))
  SAFE_ALLOCATE(geometry%chemical_symbols(geo%nspecies))
  SAFE_ALLOCATE(geometry%atom_species_names(geo%nspecies))

  do ispecies = 1, geo%nspecies
    geometry%atomic_numbers(ispecies) = species_z(geo%species(ispecies))
    geometry%chemical_symbols(ispecies) = trim(species_label(geo%species(ispecies)))
    ! according to the specification atomic_numbers is enough, but
    ! v_sim wants atom_species_name, so we use the label as name
    geometry%atom_species_names(ispecies) = trim(species_label(geo%species(ispecies)))
  end do

  ! The atoms
  SAFE_ALLOCATE(geometry%atom_species(geo%natoms))

  do i = 1, geo%natoms
    do j = 1, geo%nspecies
      if (species_z(geo%atom(i)%spec) == species_z(geo%species(j))) then
        geometry%atom_species(i) = j
        exit
      end if
    end do
  end do

  ! The coordinates
  SAFE_ALLOCATE(geometry%reduced_atom_positions(3, geo%natoms))

  offset = M_ZERO
  offset(1:geo%space%dim) = -matmul(sb%rlattice_primitive(1:geo%space%dim, 1:geo%space%dim), sb%lsize(1:geo%space%dim))

  do i = 1, geo%natoms
    ! this is only valid if the primitive vectors are along the x, y, and z directions.
    do idir = 1, 3
      geometry%reduced_atom_positions(idir, i) = (geo%atom(i)%x(idir) - offset(idir))/geometry%primitive_vectors(idir, idir)
    end do
  end do

  call etsf_io_geometry_put(ncid, geometry, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)
  
  ! Free the geometry container
  SAFE_DEALLOCATE_P(geometry%primitive_vectors)
  SAFE_DEALLOCATE_P(geometry%reduced_symmetry_matrices)
  SAFE_DEALLOCATE_P(geometry%reduced_symmetry_translations)
  SAFE_DEALLOCATE_P(geometry%space_group)
  SAFE_DEALLOCATE_P(geometry%reduced_atom_positions)
  SAFE_DEALLOCATE_P(geometry%atom_species)
  SAFE_DEALLOCATE_P(geometry%atomic_numbers)
  SAFE_DEALLOCATE_P(geometry%chemical_symbols)
  SAFE_DEALLOCATE_P(geometry%atom_species_names)

  POP_SUB(output_etsf_geometry_write)
end subroutine output_etsf_geometry_write

! --------------------------------------------------------

subroutine output_etsf_kpoints_dims(sb, dims, flags)
  type(simul_box_t),       intent(in)    :: sb
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  
  PUSH_SUB(output_etsf_kpoints_dims)

  flags%kpoints = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights

  dims%number_of_kpoints = sb%kpoints%reduced%npoints

  POP_SUB(output_etsf_kpoints_dims)
end subroutine output_etsf_kpoints_dims

! --------------------------------------------------------

subroutine output_etsf_kpoints_write(sb, ncid)
  type(simul_box_t),      intent(in)    :: sb
  integer,                intent(in)    :: ncid
  
  type(etsf_kpoints), target :: kpoints
  integer  :: nkpoints, ikpoint
  type(etsf_io_low_error)  :: error_data
  logical :: lstat

  PUSH_SUB(output_etsf_kpoints_write)

  nkpoints = sb%kpoints%reduced%npoints
  
  !Create the kpoints container
  SAFE_ALLOCATE(kpoints%reduced_coordinates_of_kpoints(1:3, 1:nkpoints))
  SAFE_ALLOCATE(kpoints%kpoint_weights(1:nkpoints))
  
  do ikpoint = 1, nkpoints
    kpoints%reduced_coordinates_of_kpoints(1:3, ikpoint) = M_ZERO
    kpoints%reduced_coordinates_of_kpoints(1:sb%dim, ikpoint) = sb%kpoints%reduced%red_point(1:sb%dim, ikpoint)
    kpoints%kpoint_weights(ikpoint) = kpoints_get_weight(sb%kpoints, ikpoint)
  end do
  
  call etsf_io_kpoints_put(ncid, kpoints, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  SAFE_DEALLOCATE_P(kpoints%reduced_coordinates_of_kpoints)
  SAFE_DEALLOCATE_P(kpoints%kpoint_weights)

  POP_SUB(output_etsf_kpoints_write)
end subroutine output_etsf_kpoints_write

! --------------------------------------------------------

subroutine output_etsf_electrons_dims(st, dims, flags)
  type(states_t),          intent(in)    :: st
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_electrons_dims)

  flags%electrons = etsf_electrons_eigenvalues + etsf_electrons_occupations + &
    etsf_electrons_number_of_electrons
  
  !Set the dimensions
  dims%number_of_spins = 1
  if (st%d%ispin == SPIN_POLARIZED) dims%number_of_spins = 2

  dims%max_number_of_states = st%nst
  dims%number_of_spinor_components = st%d%dim

  POP_SUB(output_etsf_electrons_dims)
end subroutine output_etsf_electrons_dims

! --------------------------------------------------------

subroutine output_etsf_electrons_write(st, ncid)
  type(states_t),      intent(in)    :: st
  integer,             intent(in)    :: ncid

  type(etsf_electrons) :: electrons
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  integer :: nspin, ist, ik, ispin, nkpoints

  PUSH_SUB(output_etsf_electrons_write)

  SAFE_ALLOCATE(electrons%number_of_electrons)
  electrons%number_of_electrons = st%qtot

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin

  SAFE_ALLOCATE(electrons%eigenvalues%data3D(1:st%nst, 1:nkpoints, 1:nspin))
  SAFE_ALLOCATE(electrons%occupations%data3D(1:st%nst, 1:nkpoints, 1:nspin))

  do ist = 1, st%nst
    do ik = 1, nkpoints
      do ispin = 1, nspin 
        electrons%eigenvalues%data3D(ist, ik, ispin) = st%eigenval(ist, nspin*(ik-1) + ispin)
        electrons%occupations%data3D(ist, ik, ispin) = st%occ(ist, nspin*(ik-1) + ispin)
      end do
    end do
  end do

  call etsf_io_electrons_put(ncid, electrons, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  SAFE_DEALLOCATE_P(electrons%number_of_electrons)
  SAFE_DEALLOCATE_P(electrons%eigenvalues%data3D)
  SAFE_DEALLOCATE_P(electrons%occupations%data3D)

  POP_SUB(output_etsf_electrons_write)
end subroutine output_etsf_electrons_write

! --------------------------------------------------------

subroutine output_etsf_density_dims(st, mesh, cube, dims, flags)
  type(states_t),          intent(in)    :: st
  type(mesh_t),            intent(in)    :: mesh
  type(cube_t),            intent(in)    :: cube
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_density_dims)

  flags%main = etsf_main_density

  dims%number_of_components = st%d%nspin
  dims%number_of_grid_points_vector1 = cube%n(1)
  dims%number_of_grid_points_vector2 = cube%n(2)
  dims%number_of_grid_points_vector3 = cube%n(3)
  dims%real_or_complex_density = 1
  
  POP_SUB(output_etsf_density_dims)
end subroutine output_etsf_density_dims

! --------------------------------------------------------

subroutine output_etsf_density_write(st, mesh, cube, cf, ncid)
  type(states_t),        intent(in)    :: st
  type(mesh_t),          intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  integer,               intent(in)    :: ncid

  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  integer :: ispin
  FLOAT, allocatable :: d(:), md(:,:)

  PUSH_SUB(output_etsf_density_write)

  SAFE_ALLOCATE(main%density%data4D(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1:st%d%nspin))

  if (st%d%ispin /= SPINORS) then
    do ispin = 1, st%d%nspin
      call dmesh_to_cube(mesh, st%rho(:, ispin), cube, cf, local = .true.)
      main%density%data4D(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), ispin) = cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
    end do
  else
    SAFE_ALLOCATE(md(1:mesh%np, 1:3))
    SAFE_ALLOCATE(d(1:mesh%np_part))

    d = st%rho(:, 1) + st%rho(:, 2)
    call magnetic_density(mesh, st, st%rho, md)

    call dmesh_to_cube(mesh, d, cube, cf, local = .true.)
    main%density%data4D(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1) = cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
    do ispin = 1, 3
      call dmesh_to_cube(mesh, md(:, ispin), cube, cf, local = .true.)
      main%density%data4D(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), ispin + 1) = cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
    end do
    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(md)
  end if

  call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  SAFE_DEALLOCATE_P(main%density%data4D)

  POP_SUB(output_etsf_density_write)
end subroutine output_etsf_density_write


! --------------------------------------------------------

subroutine output_etsf_wfs_rsp_dims(st, mesh, cube, dims, flags)
  type(states_t),          intent(in)    :: st
  type(mesh_t),            intent(in)    :: mesh
  type(cube_t),            intent(in)    :: cube
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_wfs_rsp_dims)

  if (states_are_real(st)) then
    dims%real_or_complex_wavefunctions = 1
  elseif(states_are_complex(st)) then
    dims%real_or_complex_wavefunctions = 2
  end if

  dims%number_of_grid_points_vector1 = cube%n(1)
  dims%number_of_grid_points_vector2 = cube%n(2)
  dims%number_of_grid_points_vector3 = cube%n(3)

  flags%main = etsf_main_wfs_rsp

  POP_SUB(output_etsf_wfs_rsp_dims)
end subroutine output_etsf_wfs_rsp_dims

! --------------------------------------------------------

subroutine output_etsf_wfs_rsp_write(st, mesh, cube, cf, ncid)
  type(states_t),        intent(in)    :: st
  type(mesh_t),          intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  integer,               intent(in)    :: ncid

  integer :: ist, ispin, ik, idim, nspin, zdim, nkpoints
  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  REAL_DOUBLE, allocatable, target :: local_wfs(:,:,:,:,:,:,:)
  FLOAT, allocatable :: dpsi(:)
  CMPLX, allocatable :: zpsi(:)

  PUSH_SUB(output_etsf_wfs_rsp_write)

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin
  if (states_are_real(st)) then
    zdim = 1
  elseif(states_are_complex(st)) then
    zdim = 2
  end if

  if (states_are_real(st)) then
    SAFE_ALLOCATE(dpsi(1:mesh%np))
  else
    SAFE_ALLOCATE(zpsi(1:mesh%np))
  end if

  SAFE_ALLOCATE(local_wfs(1:zdim, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1:st%d%dim, 1:st%nst, 1:st%d%nik))
  do ispin = 1, nspin
    do ik = 1, st%d%nik, nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          if (states_are_real(st)) then
            call states_get_state(st, mesh, idim, ist, ik + ispin - 1, dpsi)

            call dmesh_to_cube(mesh, dpsi, cube, cf, local = .true.)
            local_wfs(1, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, ist, ik+(ispin-1)*nkpoints) = &
              cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))

          else if(states_are_complex(st)) then
            call states_get_state(st, mesh, idim, ist, ik + ispin - 1, zpsi)

            call dmesh_to_cube(mesh, real(zpsi, REAL_PRECISION), cube, cf, local = .true.)
            local_wfs(1, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, ist, ik+(ispin-1)*nkpoints) = &
              cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))

            call dmesh_to_cube(mesh, aimag(zpsi), cube, cf, local = .true.)
            local_wfs(2, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, ist, ik+(ispin-1)*nkpoints) = &
              cf%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))

          end if
        end do
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(dpsi)
  SAFE_DEALLOCATE_A(zpsi)

  main%real_space_wavefunctions%data7D => local_wfs

  call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)
    
  SAFE_DEALLOCATE_A(local_wfs)

  POP_SUB(output_etsf_wfs_rsp_write)
end subroutine output_etsf_wfs_rsp_write

! --------------------------------------------------

subroutine output_etsf_basisdata_dims(st, mesh, cube, shell, dims, flags)
  type(states_t),          intent(in)    :: st
  type(mesh_t),            intent(in)    :: mesh
  type(cube_t),            intent(in)    :: cube
  type(fourier_shell_t),   intent(in)    :: shell
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_basisdata_dims)

  flags%basisdata = etsf_basisdata_basis_set + etsf_basisdata_kin_cutoff + &
    etsf_basisdata_red_coord_pw

  POP_SUB(output_etsf_basisdata_dims)

end subroutine output_etsf_basisdata_dims

! --------------------------------------------------------

subroutine output_etsf_basisdata_write(st, mesh, shell, ncid)
  type(states_t),        intent(in)    :: st
  type(mesh_t),          intent(in)    :: mesh
  type(fourier_shell_t), intent(in)    :: shell
  integer,               intent(in)    :: ncid

  type(etsf_basisdata) :: basisdata
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  integer :: ig, ng, ix, iy, iz, ixx(1:3)
 
  PUSH_SUB(output_etsf_basisdata_write)

  if((maxval(mesh%spacing(1:3)) - minval(mesh%spacing(1:3))) > CNST(1e-10)) then
    message(1) = 'Cannot generate a ETSF plane wave wave-functions file,'
    message(2) = 'spacing is not the same for each direction.'
    call messages_fatal(2)
  end if

  SAFE_ALLOCATE(basisdata%basis_set)

  basisdata%basis_set = "plane_waves"

  SAFE_ALLOCATE(basisdata%kinetic_energy_cutoff)

  basisdata%kinetic_energy_cutoff = shell%ekin_cutoff

  ng = shell%ngvectors

  SAFE_ALLOCATE(basisdata%reduced_coordinates_of_plane_waves%data2D(1:3, 1:ng))

  basisdata%reduced_coordinates_of_plane_waves%data2D(1:3, 1:ng) = shell%red_gvec(1:3, 1:ng)

  call etsf_io_basisdata_put(ncid, basisdata, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  SAFE_DEALLOCATE_P(basisdata%basis_set)
  SAFE_DEALLOCATE_P(basisdata%kinetic_energy_cutoff)
  SAFE_DEALLOCATE_P(basisdata%reduced_coordinates_of_plane_waves%data2D)

  POP_SUB(output_etsf_basisdata_write)

end subroutine output_etsf_basisdata_write

! --------------------------------------------------------

subroutine output_etsf_wfs_pw_dims(st, mesh, cube, shell, dims, flags)
  type(states_t),          intent(in)    :: st
  type(mesh_t),            intent(in)    :: mesh
  type(cube_t),            intent(in)    :: cube
  type(fourier_shell_t),   intent(in)    :: shell
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags

  PUSH_SUB(output_etsf_wfs_pw_dims)

  flags%main = etsf_main_wfs_coeff

  dims%max_number_of_coefficients = shell%ngvectors
  dims%real_or_complex_coefficients = 2

  POP_SUB(output_etsf_wfs_pw_dims)

end subroutine output_etsf_wfs_pw_dims

! --------------------------------------------------------

subroutine output_etsf_wfs_pw_write(st, mesh, cube, cf, shell, ncid)
  type(states_t),        intent(in)    :: st
  type(mesh_t),          intent(in)    :: mesh
  type(cube_t),          intent(inout) :: cube
  type(cube_function_t), intent(inout) :: cf
  type(fourier_shell_t), intent(in)    :: shell
  integer,               intent(in)    :: ncid

  type(etsf_main) :: main
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  REAL_DOUBLE, allocatable, target :: local_wfs(:, :, :, :, :, :)
  CMPLX, allocatable :: zpsi(:)
  integer :: nkpoints, nspin, zdim
  integer :: idim, ist, iq, ikpoint, ispin
  integer :: ig, ng, ix, iy, iz, ixx(1:3)

  PUSH_SUB(output_etsf_wfs_pw_write)

  nspin = 1
  if (st%d%ispin == SPIN_POLARIZED) nspin = 2

  nkpoints = st%d%nik/nspin
  zdim = 2

  ng = shell%ngvectors

  !Write the wavefunctions to the file
  SAFE_ALLOCATE(local_wfs(1:zdim, 1:ng, 1:st%d%dim, 1:st%nst, 1:nkpoints, 1:nspin))

  SAFE_ALLOCATE(zpsi(1:mesh%np))

  do iq = 1, st%d%nik
    ispin = states_dim_get_spin_index(st%d, iq)
    ikpoint = states_dim_get_kpoint_index(st%d, iq)
    do ist = 1, st%nst
      do idim = 1, st%d%dim

        ! for the moment we treat all functions as complex
        call states_get_state(st, mesh, idim, ist, iq, zpsi)
        call zmesh_to_cube(mesh, zpsi, cube, cf, local = .true.)
        call zcube_function_rs2fs(cube, cf)

        do ig = 1, ng
          ix = shell%coords(1, ig)
          iy = shell%coords(2, ig)
          iz = shell%coords(3, ig)
          
          local_wfs(1, ig, idim, ist, ikpoint, ispin) = real(cf%fs(ix, iy, iz), 8)
          local_wfs(2, ig, idim, ist, ikpoint, ispin) = aimag(cf%fs(ix, iy, iz))
        end do

      end do
    end do
  end do

  SAFE_DEALLOCATE_A(zpsi)

  main%coefficients_of_wavefunctions%data6D => local_wfs

  call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
  if (.not. lstat) call output_etsf_error(error_data)

  call etsf_io_tools_set_time_reversal_symmetry(ncid, .false., lstat, error_data)

  SAFE_DEALLOCATE_A(local_wfs)

  POP_SUB(output_etsf_wfs_pw_write)

end subroutine output_etsf_wfs_pw_write

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
