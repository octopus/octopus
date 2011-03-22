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
  integer :: i, j, idir, is, ik, idim, ix, iy, iz, nspin, zdim, isymm, nkpoints, ikpoint, ispecies, ncid
  FLOAT, allocatable :: d(:), md(:,:)
  REAL_DOUBLE, allocatable, target :: local_rho(:,:,:,:), local_wfs(:,:,:,:,:,:,:), local_ev(:, :, :)
  REAL_DOUBLE, allocatable, target :: local_occ(:,:,:)
  type(etsf_io_low_error)  :: error_data
  logical :: lstat
  type(etsf_dims) :: dims
  type(etsf_groups_flags), target :: flags
  type(etsf_main), target :: main
  type(etsf_electrons), target :: electrons
#endif

  PUSH_SUB(h_sys_output_etsf)

#ifndef HAVE_ETSF_IO
  ASSERT(.false.)
#endif

  !Create a cube
  call cube_function_init(cube, gr%mesh%idx%ll)
  call dcube_function_alloc_RS(cube)

#ifdef HAVE_ETSF_IO
  if (iand(outp%what, output_geometry).ne.0) then

    call output_etsf_geometry_dims(geo, gr%sb, dims, flags)

    call etsf_io_data_init(dir//"/geometry-etsf.nc", flags, dims, &
      "Crystallographic_data file", &
      "Created by "//PACKAGE_STRING, &
      lstat, error_data, overwrite=.true.)
    if (.not. lstat) call output_etsf_error(error_data)

    call etsf_io_low_open_modify(ncid, dir//"/geometry-etsf.nc", lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
    
    call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call output_etsf_geometry_write(geo, gr%sb, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
  end if


  if (iand(outp%what, output_density).ne.0) then

    !Set the dimensions
    call output_etsf_geometry_dims(geo, gr%sb, dims, flags)

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
    if (.not. lstat) call output_etsf_error(error_data)

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

    ! write the geometry
    call etsf_io_low_open_modify(ncid, dir//"/density-etsf.nc", lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
    
    call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call output_etsf_geometry_write(geo, gr%sb, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    !Free the main container
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
    call output_etsf_geometry_dims(geo, gr%sb, dims, flags)
    call output_etsf_kpoints_dims(gr%sb, dims, flags)

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
    dims%number_of_spins = nspin
    dims%number_of_spinor_components = st%d%dim
    dims%number_of_grid_points_vector1 = cube%n(1)
    dims%number_of_grid_points_vector2 = cube%n(2)
    dims%number_of_grid_points_vector3 = cube%n(3)
    dims%real_or_complex_wavefunctions = zdim

    !Create the electrons container
    flags%electrons = etsf_electrons_eigenvalues + etsf_electrons_occupations + &
      etsf_electrons_number_of_electrons

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
    SAFE_ALLOCATE(electrons%number_of_electrons)
    electrons%number_of_electrons = st%qtot
    electrons%eigenvalues%data3D => local_ev
    electrons%occupations%data3D => local_occ

    !Open the file
    flags%main = etsf_main_wfs_rsp
    call etsf_io_data_init(dir//"/wfs-etsf.nc", flags, dims, &
      "Wavefunctions file", &
      "Created by "//PACKAGE_STRING, &
      lstat, error_data, overwrite=.true.)
    if (.not. lstat) call output_etsf_error(error_data)

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

    call etsf_io_low_open_modify(ncid, dir//"/wfs-etsf.nc", lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
    
    call etsf_io_low_set_write_mode(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call etsf_io_main_put(ncid, main, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)
    
    call etsf_io_electrons_put(ncid, electrons, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    call output_etsf_geometry_write(geo, gr%sb, ncid)
    call output_etsf_kpoints_write(gr%sb, ncid)

    call etsf_io_low_close(ncid, lstat, error_data = error_data)
    if (.not. lstat) call output_etsf_error(error_data)

    !Free the main container
    nullify(main%real_space_wavefunctions%data7D)
    SAFE_DEALLOCATE_A(local_wfs)

    !Free the electrons container
    SAFE_DEALLOCATE_A(local_occ)
    SAFE_DEALLOCATE_A(local_ev)
    SAFE_DEALLOCATE_P(electrons%number_of_electrons)

    nullify(electrons%eigenvalues%data3D)
    nullify(electrons%occupations%data3D)
    flags%electrons = etsf_electrons_none

    !Reset the dimensions
    dims%max_number_of_states = 1
    dims%number_of_spins = 1
    dims%number_of_spinor_components = 1
    dims%number_of_grid_points_vector1 = 1
    dims%number_of_grid_points_vector2 = 1
    dims%number_of_grid_points_vector3 = 1
    dims%real_or_complex_wavefunctions = 1
  end if
#endif

  call cube_function_end(cube)

  POP_SUB(h_sys_output_etsf)
end subroutine h_sys_output_etsf

#ifdef HAVE_ETSF_IO
! --------------------------------------------------------

subroutine output_etsf_error(error_data)
  type(etsf_io_low_error), intent(in) :: error_data
  
  call etsf_io_low_error_handle(error_data)
  message(1) = "ETSF_IO returned a fatal error. See message above."
  call write_fatal(1)

end subroutine output_etsf_error

! --------------------------------------------------------

subroutine output_etsf_geometry_dims(geo, sb, dims, flags)
  type(geometry_t),        intent(in)    :: geo
  type(simul_box_t),       intent(in)    :: sb
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  
  flags%geometry = etsf_geometry_all - etsf_geometry_valence_charges - etsf_geometry_pseudo_types

  dims%number_of_symmetry_operations = symmetries_number(sb%symm)
  dims%number_of_atom_species = geo%nspecies
  dims%number_of_atoms = geo%natoms

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

  offset = -matmul(sb%rlattice_primitive, sb%lsize)

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

end subroutine output_etsf_geometry_write

! --------------------------------------------------------

subroutine output_etsf_kpoints_dims(sb, dims, flags)
  type(simul_box_t),       intent(in)    :: sb
  type(etsf_dims),         intent(inout) :: dims
  type(etsf_groups_flags), intent(inout) :: flags
  
  flags%kpoints = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights

  dims%number_of_kpoints = sb%kpoints%reduced%npoints

end subroutine output_etsf_kpoints_dims

! --------------------------------------------------------

subroutine output_etsf_kpoints_write(sb, ncid)
  type(simul_box_t),      intent(in)    :: sb
  integer,                intent(in)    :: ncid
  
  type(etsf_kpoints), target :: kpoints
  integer  :: nkpoints, ikpoint
  type(etsf_io_low_error)  :: error_data
  logical :: lstat

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

end subroutine output_etsf_kpoints_write

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
