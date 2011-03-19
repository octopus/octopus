!! Copyright (C) 2009
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

subroutine h_sys_output_etsf(st, gr, geo, dir, outp)
  type(states_t),         intent(in) :: st
  type(grid_t),           intent(in) :: gr
  type(geometry_t),       intent(in) :: geo
  character(len=*),       intent(in) :: dir
  type(h_sys_output_t),   intent(in) :: outp

  type(cube_function_t) :: cube
#ifdef HAVE_ETSF_IO
  logical :: lstat
  integer :: i, j, idir, is, ik, idim, ix, iy, iz, nkpoints, nspin, zdim, isymm, ikpoint, ispecies
  FLOAT   :: offset(MAX_DIM)
  FLOAT, allocatable :: d(:), md(:,:)
  REAL_DOUBLE, allocatable, target :: local_rho(:,:,:,:), local_wfs(:,:,:,:,:,:,:), local_ev(:, :, :)
  REAL_DOUBLE, allocatable, target :: local_occ(:,:,:), local_red_coord_kpt(:,:), local_kpoint_weights(:)
  type(etsf_io_low_error)  :: error_data
  type(etsf_dims) :: dims
  type(etsf_groups) :: groups
  type(etsf_groups_flags), target :: flags
  type(etsf_main), target :: main
  type(etsf_geometry), target :: geometry
  type(etsf_electrons), target :: electrons
  type(etsf_kpoints), target :: kpoints
#endif

  PUSH_SUB(h_sys_output_etsf)

#ifndef HAVE_ETSF_IO
  ASSERT(.false.)
#endif

  !Create a cube
  call cube_function_init(cube, gr%mesh%idx%ll)
  call dcube_function_alloc_RS(cube)

#ifdef HAVE_ETSF_IO

  ! First initialize the geometry, included in all files
  flags%geometry = etsf_geometry_all - etsf_geometry_valence_charges - etsf_geometry_pseudo_types
  groups%geometry => geometry

  ! Primitive vectors
  SAFE_ALLOCATE(geometry%primitive_vectors(1:3, 1:3))
  do idir = 1, gr%sb%dim
    geometry%primitive_vectors(1:3, idir) = gr%sb%rlattice(1:3, idir)
  end do

  ! The symmetries
  SAFE_ALLOCATE(geometry%space_group)
  geometry%space_group = symmetries_space_group_number(gr%sb%symm)
  dims%number_of_symmetry_operations = symmetries_number(gr%sb%symm)
  SAFE_ALLOCATE(geometry%reduced_symmetry_matrices(1:3, 1:3, 1:symmetries_number(gr%sb%symm)))
  SAFE_ALLOCATE(geometry%reduced_symmetry_translations(1:3, 1:symmetries_number(gr%sb%symm)))

  do isymm = 1, symmetries_number(gr%sb%symm)
    geometry%reduced_symmetry_matrices(1:3, 1:3, isymm) = symm_op_rotation_matrix(gr%sb%symm%ops(isymm))
    geometry%reduced_symmetry_translations(1:3, isymm) = symm_op_translation_vector(gr%sb%symm%ops(isymm))
  end do

  ! The species
  dims%number_of_atom_species = geo%nspecies

  SAFE_ALLOCATE(geometry%atomic_numbers(geo%nspecies))
  SAFE_ALLOCATE(geometry%chemical_symbols(geo%nspecies))
  SAFE_ALLOCATE(geometry%atom_species_names(geo%nspecies))

  do ispecies = 1, geo%nspecies
    geometry%atomic_numbers(ispecies) = species_z(geo%species(ispecies))
    geometry%chemical_symbols(ispecies) = trim(species_label(geo%species(ispecies)))
    ! according to the specification atomic_numbers is enough, but
    ! v_sim wants atom_species_name, so we use the label
    geometry%atom_species_names(ispecies) = trim(species_label(geo%species(ispecies)))
  end do

  ! The atoms
  dims%number_of_atoms = geo%natoms
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
  SAFE_ALLOCATE(geometry%reduced_atom_positions(3,geo%natoms))

  offset = -matmul(gr%sb%rlattice_primitive, gr%sb%lsize)
  do i = gr%sb%periodic_dim+1, 3
    offset(i)=-(cube%n(i) - 1)/2 * gr%mesh%spacing(i)
  end do

  do i = 1, geo%natoms
    ! this is only valid if the primitive vectors are along the x, y, and z directions.
    do idir = 1, 3
      geometry%reduced_atom_positions(idir, i) = (geo%atom(i)%x(idir) - offset(idir))/geometry%primitive_vectors(idir, idir)
    end do
  end do

  if (iand(outp%what, output_geometry).ne.0) then
    call etsf_io_data_init(dir//"/geometry-etsf.nc", flags, dims, &
      "Crystallographic_data file", &
      "Created by "//PACKAGE_STRING, &
      lstat, error_data, overwrite=.true.)
    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if

    !Write the geometry to the file
    call etsf_io_data_write(dir//"/geometry-etsf.nc", &
      groups, lstat, error_data)

    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if
  end if


  if (iand(outp%what, output_density).ne.0) then
    !Set the dimensions
    dims%number_of_components = st%d%nspin
    dims%number_of_grid_points_vector1 = cube%n(1)
    dims%number_of_grid_points_vector2 = cube%n(2)
    dims%number_of_grid_points_vector3 = cube%n(3)
    dims%real_or_complex_density = 1

    !Open the file
    flags%main = etsf_main_density
    call etsf_io_data_init(dir//"/density-etsf.nc", flags, dims, &
      & "Density file", &
      & "Created by "//PACKAGE_STRING, &
      & lstat, error_data, overwrite=.true.)
    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if

    !Write the density to the file
    SAFE_ALLOCATE(local_rho(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1:st%d%nspin))
    if (st%d%ispin /= SPINORS) then
      do i = 1, st%d%nspin
        call dmesh_to_cube(gr%mesh, st%rho(:,i), cube, local = .true.)
        local_rho(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), i) = cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
      end do
    else
      SAFE_ALLOCATE(md(1:gr%mesh%np, 1:3))
      SAFE_ALLOCATE(d(1:gr%mesh%np_part))

      d = st%rho(:, 1) + st%rho(:, 2)
      call magnetic_density(gr%mesh, st, st%rho, md)

      call dmesh_to_cube(gr%mesh, d, cube, local = .true.)
      local_rho(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1) = cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
      do i = 1, 3
        call dmesh_to_cube(gr%mesh, md(:,i), cube, local = .true.)
        local_rho(1:cube%n(1), 1:cube%n(2), 1:cube%n(3), i+1) = cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
      end do
      SAFE_DEALLOCATE_A(d)
      SAFE_DEALLOCATE_A(md)
    end if
    main%density%data4D => local_rho
    groups%main => main
    call etsf_io_data_write(dir//"/density-etsf.nc", &
      & groups, lstat, error_data)
    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if

    !Free the main container
    nullify(groups%main)
    nullify(main%density%data4D)
    SAFE_DEALLOCATE_A(local_rho)

    !Reset the dimensions
    dims%number_of_components = 1
    dims%number_of_grid_points_vector1 = 1
    dims%number_of_grid_points_vector2 = 1
    dims%number_of_grid_points_vector3 = 1
    dims%real_or_complex_density = 1
  end if


  if (iand(outp%what, output_wfs).ne.0) then
    !Set the dimensions
    nspin = 1
    if (st%d%ispin == SPIN_POLARIZED) nspin = 2
    nkpoints = st%d%nik/nspin
    if (states_are_real(st)) then
      zdim = 1
    elseif(states_are_complex(st)) then
      zdim = 2
    end if
    dims%max_number_of_states = st%nst
    dims%number_of_kpoints = nkpoints
    dims%number_of_spins = nspin
    dims%number_of_spinor_components = st%d%dim
    dims%number_of_grid_points_vector1 = cube%n(1)
    dims%number_of_grid_points_vector2 = cube%n(2)
    dims%number_of_grid_points_vector3 = cube%n(3)
    dims%real_or_complex_wavefunctions = zdim

    !Create the electrons container
    flags%electrons = etsf_electrons_eigenvalues + etsf_electrons_occupations + etsf_electrons_number_of_states
    SAFE_ALLOCATE(local_ev(1:st%nst, 1:nkpoints, 1:nspin))
    SAFE_ALLOCATE(local_occ(1:st%nst, 1:nkpoints, 1:nspin))
    do i = 1, st%nst
      do ik = 1, nkpoints
        do is = 1, nspin 
          local_ev(i, ik, is) = st%eigenval(i, nspin*(ik-1) + is)
          local_occ(i, ik, is) = st%occ(i, nspin*(ik-1) + is)
        end do
      end do
    end do
    electrons%eigenvalues%data3D => local_ev
    electrons%occupations%data3D => local_occ
    groups%electrons => electrons

    !Create the kpoints container
    flags%kpoints = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
    SAFE_ALLOCATE(local_red_coord_kpt(1:3, 1:nkpoints))
    SAFE_ALLOCATE(local_kpoint_weights(1:nkpoints))

    do ikpoint = 1, nkpoints
      local_red_coord_kpt(1:3, ikpoint) = M_ZERO
      local_red_coord_kpt(1:gr%sb%dim, ikpoint) = gr%sb%kpoints%reduced%red_point(1:gr%sb%dim, ikpoint)
      local_kpoint_weights(ikpoint) = kpoints_get_weight(gr%sb%kpoints, ikpoint)
    end do

    kpoints%reduced_coordinates_of_kpoints => local_red_coord_kpt
    kpoints%kpoint_weights => local_kpoint_weights
    groups%kpoints => kpoints

    !Open the file
    flags%main = etsf_main_wfs_rsp
    call etsf_io_data_init(dir//"/wfs-etsf.nc", flags, dims, &
      & "Wavefunctions file", &
      & "Created by "//PACKAGE_STRING, &
      & lstat, error_data, overwrite=.true.)
    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if

    !Write the wavefunctions to the file
    SAFE_ALLOCATE(local_wfs(1:zdim, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), 1:st%d%dim, 1:st%nst, 1:st%d%nik))
    do is = 1, nspin
      do ik = 1, st%d%nik, nspin
        do i = 1, st%nst
          do idim = 1, st%d%dim
            if (states_are_real(st)) then
              call dmesh_to_cube(gr%mesh, st%dpsi(1:gr%mesh%np_part, idim, i, ik+is-1), cube, local = .true.)
              local_wfs(1, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, i, ik+(is-1)*nkpoints) = &
                cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))
              
            else if(states_are_complex(st)) then
              call dmesh_to_cube(gr%mesh, &
                real(st%zpsi(1:gr%mesh%np_part, idim, i, ik+is-1), REAL_PRECISION), cube, local = .true.)
              local_wfs(1, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, i, ik+(is-1)*nkpoints) = &
                cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))

              call dmesh_to_cube(gr%mesh, aimag(st%zpsi(1:gr%mesh%np_part, idim, i, ik+is-1)), cube, local = .true.)
              local_wfs(2, 1:cube%n(1), 1:cube%n(2), 1:cube%n(3), idim, i, ik+(is-1)*nkpoints) = &
                cube%dRS(1:cube%n(1), 1:cube%n(2), 1:cube%n(3))

            end if
          end do
        end do
      end do
    end do

    main%real_space_wavefunctions%data7D => local_wfs
    groups%main => main
    call etsf_io_data_write(dir//"/wfs-etsf.nc", &
      & groups, lstat, error_data)
    if (.not. lstat) then
      call etsf_io_low_error_handle(error_data)
      message(1) = "ETSF_IO returned a fatal error. See message above."
      call write_fatal(1)
    end if

    !Free the main container
    nullify(groups%main)
    nullify(main%real_space_wavefunctions%data7D)
    SAFE_DEALLOCATE_A(local_wfs)

    !Free the k-points container
    SAFE_DEALLOCATE_A(local_red_coord_kpt)
    SAFE_DEALLOCATE_A(local_kpoint_weights)
    nullify(groups%kpoints)
    nullify(kpoints%reduced_coordinates_of_kpoints)
    nullify(kpoints%kpoint_weights)
    flags%kpoints = etsf_kpoints_none

    !Free the electrons container
    SAFE_DEALLOCATE_A(local_occ)
    SAFE_DEALLOCATE_A(local_ev)
    nullify(groups%electrons)
    nullify(electrons%eigenvalues%data3D)
    nullify(electrons%occupations%data3D)
    flags%electrons = etsf_electrons_none

    !Reset the dimensions
    dims%max_number_of_states = 1
    dims%number_of_kpoints = 1
    dims%number_of_spins = 1
    dims%number_of_spinor_components = 1
    dims%number_of_grid_points_vector1 = 1
    dims%number_of_grid_points_vector2 = 1
    dims%number_of_grid_points_vector3 = 1
    dims%real_or_complex_wavefunctions = 1
  end if

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
#endif

  call cube_function_end(cube)

  POP_SUB(h_sys_output_etsf)
end subroutine h_sys_output_etsf

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
