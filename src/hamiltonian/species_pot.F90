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

module species_pot_oct_m
  use atom_oct_m
  use curvilinear_oct_m
  use global_oct_m
  use io_function_oct_m
  use index_oct_m
  use lattice_vectors_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use ps_oct_m
  use root_solver_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use splines_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use volume_oct_m

  implicit none

  private
  public ::                         &
    species_get_long_range_density, &
    species_get_nlcc,               &
    species_get_nlcc_grad,          &
    species_get_local,              &
    species_atom_density,           &
    species_atom_density_np,        &
    species_atom_density_derivative,&
    species_atom_density_derivative_np, & 
    species_atom_density_grad

  type(mesh_t), pointer :: mesh_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT, pointer :: pos_p(:)

contains


  ! ---------------------------------------------------------
  subroutine species_atom_density(mesh, space, namespace, sb, atom, spin_channels, rho)
    type(mesh_t),         intent(in)    :: mesh
    type(space_t),        intent(in)    :: space
    type(namespace_t),    intent(in)    :: namespace
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: rho(:, :) !< (mesh%np, spin_channels)

    integer :: isp, ip, in_points, icell
    FLOAT :: rr, x, pos(space%dim), nrm, rmax
    FLOAT :: xx(space%dim), yy(space%dim), rerho, imrho
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps
    type(volume_t) :: volume

#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(lattice_iterator_t) :: latt_iter

    PUSH_SUB(species_atom_density)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    species => atom%species
    rho = M_ZERO

    ! build density ...
    select case (species_type(species))
    case (SPECIES_FROM_FILE, SPECIES_USDEF, SPECIES_SOFT_COULOMB, SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN) ! ... from userdef
      do isp = 1, spin_channels
        rho(1:mesh%np, isp) = M_ONE
        x = (species_zval(species)/TOFLOAT(spin_channels)) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPECIES_CHARGE_DENSITY, SPECIES_JELLIUM_CHARGE_DENSITY)
      ! We put, for the electron density, the same as the positive density that 
      ! creates the external potential.
      ! This code is repeated in get_density, and should therefore be cleaned!!!!!

      if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
        call volume_init(volume)
        call volume_read_from_block(volume, namespace, trim(species_rho_string(species)))
      end if

      latt_iter = lattice_iterator_t(sb%latt, maxval(norm2(sb%latt%rlattice, dim=1)))
      rho = M_ZERO
      do icell = 1, latt_iter%n_cells
        yy = latt_iter%get(icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = atom%x, coords = xx)
          xx = xx + yy
          rr = norm2(xx)
          
          rerho = M_ZERO
          if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
            if(volume_in_volume(space, volume, xx)) rerho = M_ONE
          else
            call parse_expression(rerho, imrho, space%dim, xx, rr, M_ZERO, trim(species_rho_string(species)))
          end if
          rho(ip, 1) = rho(ip, 1) + rerho
       end do
      end do

      if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
         call volume_end(volume)
      end if

      if(spin_channels > 1) then
        rho(:, 1) = M_HALF*rho(:, 1)
        rho(:, 2) = rho(:, 1)
      end if

      ! rescale to match the valence charge
      do isp = 1, spin_channels
        x = species_zval(species) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPECIES_JELLIUM) ! ... from jellium
      in_points = 0
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = atom%x)
        if(rr <= species_jradius(species)) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_oct_m module.
 
        if (mesh%use_curvilinear) then
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(species)) then
              rho(ip, 1:spin_channels) = species_zval(species) /   &
                (mesh%vol_pp(ip)*TOFLOAT(in_points*spin_channels))
            end if
          end do
        else
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(species)) then
              rho(ip, 1:spin_channels) = species_zval(species) /   &
                (mesh%vol_pp(1)*TOFLOAT(in_points*spin_channels))
            end if
          end do
        end if
      end if

    case (SPECIES_JELLIUM_SLAB) ! ... from jellium slab
      in_points = 0
      do ip = 1, mesh%np
        rr = abs( mesh%x( ip, 3 ) )
        if( rr <= species_jthick(species)/M_TWO ) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, mesh%mpi_grp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_oct_m module.

        if (mesh%use_curvilinear) then
          do ip = 1, mesh%np
            rr = abs( mesh%x( ip, 3 ) )
            if( rr <= species_jthick(species)/M_TWO ) then
              rho(ip, 1:spin_channels) = species_zval(species) /   &
                (mesh%vol_pp(ip)*TOFLOAT(in_points*spin_channels))
            end if
          end do
        else
          do ip = 1, mesh%np
            rr = abs( mesh%x( ip, 3 ) )
            if( rr <= species_jthick(species)/M_TWO ) then
              rho(ip, 1:spin_channels) = species_zval(species) /   &
                (mesh%vol_pp(1)*TOFLOAT(in_points*spin_channels))
            end if
          end do
        end if
      end if

    case (SPECIES_PSEUDO, SPECIES_PSPIO)
      ! ...from pseudopotentials
      
      ps => species_ps(species)

      if(ps_has_density(ps)) then

        ASSERT(allocated(ps%density))

        rmax = CNST(0.0)
        do isp = 1, spin_channels
          rmax = max(rmax, spline_cutoff_radius(ps%density(isp), ps%projectors_sphere_threshold))
        end do
        
        latt_iter = lattice_iterator_t(sb%latt, rmax)
        do icell = 1, latt_iter%n_cells
          pos = atom%x(1:space%dim) + latt_iter%get(icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, R_SMALL)
            
            do isp = 1, spin_channels
              if(rr >= spline_range_max(ps%density(isp))) cycle
              rho(ip, isp) = rho(ip, isp) + spline_eval(ps%density(isp), rr)
            end do
            
          end do
        end do

      else 

        !we use the square root of the short-range local potential, just to put something that looks like a density

        latt_iter = lattice_iterator_t(sb%latt, spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold))
        do icell = 1, latt_iter%n_cells
          pos = atom%x(1:space%dim) + latt_iter%get(icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, R_SMALL)

            if(rr >= spline_range_max(ps%vl)) cycle
            
            do isp = 1, spin_channels
              rho(ip, isp) = rho(ip, isp) + sqrt(abs(spline_eval(ps%vl, rr)))
            end do
              
          end do
        end do

        ! normalize
        nrm = CNST(0.0)
        do isp = 1, spin_channels
          nrm = nrm + dmf_integrate(mesh, rho(:, isp))
        end do

        rho(1:mesh%np, 1:spin_channels) = rho(1:mesh%np, 1:spin_channels)*species_zval(species)/nrm

      end if

    end select

    POP_SUB(species_atom_density)
  end subroutine species_atom_density

  ! ---------------------------------------------------------
  ! A non periodized version of the routine species_atom_density
  ! This is used for the Hirshfeld routines
  ! TODO: implement it for other approaches than pseudo potentials.
 subroutine species_atom_density_np(mesh, atom, namespace, pos,  spin_channels, rho)
    type(mesh_t),         intent(in)    :: mesh
    type(atom_t), target, intent(in)    :: atom
    type(namespace_t),    intent(in)    :: namespace
    FLOAT,                intent(in)    :: pos(:)
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: rho(:, :) !< (mesh%np, spin_channels)
    integer :: isp, ip
    FLOAT :: rr, nrm, rmax
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps


    PUSH_SUB(species_atom_density_np)
 
    rho = M_ZERO
    species => atom%species
    ps => species_ps(species)
    select case (species_type(species))
    case (SPECIES_PSEUDO, SPECIES_PSPIO)
      ! ...from pseudopotentials

      if(ps_has_density(ps)) then

        ASSERT(allocated(ps%density))

        rmax = CNST(0.0)

        do isp = 1, spin_channels
          rmax = max(rmax, spline_cutoff_radius(ps%density(isp), ps%projectors_sphere_threshold))
        end do
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos)

          rr = max(rr, R_SMALL) 
           
          do isp = 1, spin_channels
            if(rr >= spline_range_max(ps%density(isp))) cycle
            rho(ip, isp) = rho(ip, isp) + spline_eval(ps%density(isp), rr)
          end do

        end do

      else

        !we use the square root of the short-range local potential, just to put something that looks like a density

        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos)
          rr = max(rr, R_SMALL)

          if(rr >= spline_range_max(ps%vl)) cycle

          do isp = 1, spin_channels
            rho(ip, isp) = rho(ip, isp) + sqrt(abs(spline_eval(ps%vl, rr)))
          end do

        end do

        ! normalize
        nrm = CNST(0.0)
        do isp = 1, spin_channels
          nrm = nrm + dmf_integrate(mesh, rho(:, isp))
        end do

        rho(1:mesh%np, 1:spin_channels) = rho(1:mesh%np, 1:spin_channels)*species_zval(species)/nrm

      end if
    case default
      call messages_not_implemented('species_atom_density_np for non-pseudopotential species', namespace=namespace)

    end select

    POP_SUB(species_atom_density_np)
  end subroutine species_atom_density_np


  ! ---------------------------------------------------------

  subroutine species_atom_density_derivative(mesh, space, sb, atom, namespace, spin_channels, drho)
    type(mesh_t),         intent(in)    :: mesh
    type(space_t),        intent(in)    :: space
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    type(namespace_t),    intent(in)    :: namespace
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: drho(:, :) !< (mesh%np, spin_channels)

    integer :: icell
    FLOAT :: pos(space%dim), range
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps
    type(lattice_iterator_t) :: latt_iter

    PUSH_SUB(species_atom_density_derivative)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    species => atom%species
    drho = M_ZERO

    ! build density ...
    select case (species_type(species))

    case (SPECIES_PSEUDO, SPECIES_PSPIO)
      ! ...from pseudopotentials
      
      ps => species_ps(species)

      if(ps_has_density(ps)) then

        range = spline_cutoff_radius(ps%density_der(1), ps%projectors_sphere_threshold)
        if (spin_channels == 2) range = max(range, spline_cutoff_radius(ps%density_der(2), ps%projectors_sphere_threshold))
        latt_iter = lattice_iterator_t(sb%latt, range)

        do icell = 1, latt_iter%n_cells
          pos = atom%x(1:space%dim) + latt_iter%get(icell)
          call species_atom_density_derivative_np(mesh, atom, namespace, pos, spin_channels,  drho)
        end do
      end if
       
    case default
      call messages_not_implemented('species_atom_density_derivative for non-pseudopotential species', namespace=namespace)

    end select

    POP_SUB(species_atom_density_derivative)
  end subroutine species_atom_density_derivative

  subroutine species_atom_density_derivative_np(mesh, atom, namespace, pos, spin_channels,  drho)
    type(mesh_t),         intent(in)    :: mesh
    type(atom_t),         intent(in)    :: atom
    type(namespace_t),    intent(in)    :: namespace
    FLOAT,                intent(in)    :: pos(:)
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: drho(:, :) !< (mesh%np, spin_channels)

    integer :: isp, ip
    FLOAT :: rr
    type(ps_t), pointer :: ps

    PUSH_SUB(species_atom_density_derivative_np)


    ps => species_ps(atom%species)

    if(ps_has_density(ps)) then

      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = pos)
        rr = max(rr, R_SMALL)

        do isp = 1, spin_channels
          if(rr >= spline_range_max(ps%density_der(isp))) cycle
          drho(ip, isp) = drho(ip, isp) + spline_eval(ps%density_der(isp), rr)
        end do

      end do
    else
      call messages_write('The pseudopotential for')
      call messages_write(species_label(atom%species))
      call messages_write(' does not contain the density.')
      call messages_fatal(namespace=namespace)
    end if

    POP_SUB(species_atom_density_derivative_np)
  end subroutine species_atom_density_derivative_np
  

  ! ---------------------------------------------------------
  ! Gradient of the atomic density, if available
  subroutine species_atom_density_grad(mesh, space, sb, atom, namespace, spin_channels, drho)
    type(mesh_t),         intent(in)    :: mesh
    type(space_t),        intent(in)    :: space
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    type(namespace_t),    intent(in)    :: namespace
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: drho(:, :, :) !< (mesh%np, spin_channels, dim)

    integer :: isp, ip, icell, idir
    FLOAT :: rr, pos(space%dim), range, spline
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps
    type(lattice_iterator_t) :: latt_iter

    PUSH_SUB(species_atom_density_grad)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    species => atom%species
    drho = M_ZERO

    ! build density ...
    select case (species_type(species))

    case (SPECIES_PSEUDO, SPECIES_PSPIO)
      ! ...from pseudopotentials
      
      ps => species_ps(species)

      if(ps_has_density(ps)) then

        range = spline_cutoff_radius(ps%density_der(1), ps%projectors_sphere_threshold)
        if (spin_channels == 2) range = max(range, spline_cutoff_radius(ps%density_der(2), ps%projectors_sphere_threshold))
        latt_iter = lattice_iterator_t(sb%latt, range)

        do icell = 1, latt_iter%n_cells
          pos = atom%x(1:space%dim) + latt_iter%get(icell)
        
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, R_SMALL)

            do isp = 1, spin_channels
              if(rr >= spline_range_max(ps%density_der(isp))) cycle
              spline = spline_eval(ps%density_der(isp), rr)

              do idir = 1, space%dim
                drho(ip, isp, idir) = drho(ip, isp, idir) - spline*(mesh%x(ip, idir) - pos(idir))/rr
              end do
           end do
          end do
        end do
  
      else 
        call messages_write('The pseudopotential for')
        call messages_write(species_label(species))
        call messages_write(' does not contain the density.')
        call messages_fatal(namespace=namespace)
      end if
      
    case default
      call messages_not_implemented('species_atom_density_grad for non-pseudopotential species', namespace=namespace)

    end select

    POP_SUB(species_atom_density_grad)
  end subroutine species_atom_density_grad

  ! ---------------------------------------------------------

  subroutine species_get_long_range_density(species, space, namespace, pos, mesh, rho)
    type(species_t),    target, intent(in)  :: species
    type(space_t),              intent(in)  :: space
    type(namespace_t),          intent(in)  :: namespace
    FLOAT,              target, intent(in)  :: pos(1:space%dim)
    type(mesh_t),       target, intent(in)  :: mesh
    FLOAT,                      intent(out) :: rho(:)

    type(root_solver_t) :: rs
    logical :: conv
    FLOAT   :: x(space%dim+1), startval(space%dim + 1)
    FLOAT   :: delta, alpha, beta, xx(space%dim), yy(space%dim), rr, imrho1, rerho
    FLOAT   :: dist2, dist2_min
    integer :: icell, ipos, ip
    type(lattice_iterator_t) :: latt_iter
    type(ps_t), pointer :: ps
    type(volume_t) :: volume
    logical :: have_point
#ifdef HAVE_MPI
    FLOAT   :: local_min(2), global_min(2)
#endif
    type(submesh_t)       :: sphere
    type(profile_t), save :: prof
    FLOAT,    allocatable :: rho_sphere(:)
    FLOAT, parameter      :: threshold = CNST(1e-6)
    FLOAT                 :: norm_factor
    
    PUSH_SUB(species_get_long_range_density)

    call profiling_in(prof, "SPECIES_LONG_RANGE_DENSITY")

    select case(species_type(species))

    case(SPECIES_PSEUDO, SPECIES_PSPIO)
      ps => species_ps(species)

      call submesh_init(sphere, space, mesh%sb, mesh, pos, spline_cutoff_radius(ps%nlr, threshold))
      SAFE_ALLOCATE(rho_sphere(1:sphere%np))

      do ip = 1, sphere%np
        rho_sphere(ip) = sphere%x(ip, 0)
      end do
      if(sphere%np > 0) call spline_eval_vec(ps%nlr, sphere%np, rho_sphere)

      rho(1:mesh%np) = M_ZERO

      ! A small amount of charge is missing with the cutoff, we
      ! renormalize so that the long range potential is exact
      norm_factor = abs(species_zval(species)/dsm_integrate(mesh, sphere, rho_sphere))
      
      do ip = 1, sphere%np
        rho(sphere%map(ip)) = rho(sphere%map(ip)) + norm_factor*rho_sphere(ip)
      end do

      SAFE_DEALLOCATE_A(rho_sphere)
      call submesh_end(sphere)
      nullify(ps)

    case(SPECIES_FULL_DELTA)

      dist2_min = huge(delta)
      ipos = 0

      do ip = 1, mesh%np

        rho(ip) = M_ZERO

        dist2 = sum((mesh%x(ip, :) - pos)**2)
        if (dist2 < dist2_min) then
          ipos = ip
          dist2_min = dist2
        end if

      end do

      have_point = .true.
#ifdef HAVE_MPI
      ! in parallel we have to find the minimum of the whole grid
      if(mesh%parallel_in_domains) then

        local_min = (/ dist2_min, TOFLOAT(mesh%mpi_grp%rank)/)
        call MPI_Allreduce(local_min, global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)

        if(mesh%mpi_grp%rank /= nint(global_min(2))) have_point = .false.

        dist2_min = global_min(1)
      end if
#endif

      if(have_point) then
        if(mesh%use_curvilinear) then
          rho(ipos) = -species_z(species)/mesh%vol_pp(ipos)
        else
          rho(ipos) = -species_z(species)/mesh%vol_pp(1)
        end if
      end if

      write(message(1), '(3a,f5.2,3a)') &
        "Info: species_full_delta species ", trim(species_label(species)), &
        " atom displaced ", units_from_atomic(units_out%length, sqrt(dist2_min)), &
        " [ ", trim(units_abbrev(units_out%length)), " ]"
      call messages_info(1)

    case(SPECIES_FULL_GAUSSIAN)

      ! periodic copies are not considered in this routine
      if (space%is_periodic()) then
        call messages_experimental("species_full_gaussian for periodic systems")
      end if

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      SAFE_ALLOCATE(rho_p(1:mesh%np))
      SAFE_ALLOCATE(grho_p(1:mesh%np, 1:space%dim+1))

      mesh_p => mesh
      pos_p => pos

      ! Initial guess.
      delta   = mesh%spacing(1)
      alpha   = sqrt(M_TWO)*species_sigma(species)*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:space%dim) = CNST(1.0)

      ! the dim+1 variable is the normalization of the delta function
      startval(space%dim+1) = beta

      ! get a better estimate for beta
      call getrho(space%dim, startval)
      beta = M_ONE / dmf_integrate(mesh, rho_p)
      startval(space%dim + 1) = beta

      ! solve equation
      call root_solver_init(rs, namespace, space%dim+1, solver_type=ROOT_NEWTON, maxiter=500, abs_tolerance=CNST(1.0e-10))
      call droot_solver_run(rs, func, x, conv, startval=startval)

      if(.not.conv) then
        write(message(1),'(a)') 'Internal error in species_get_density.'
        call messages_fatal(1, namespace=namespace)
      end if

      ! we want a charge of -Z
      rho = -species_z(species)*rho_p

      nullify(mesh_p)
      nullify(pos_p)
      SAFE_DEALLOCATE_A(grho_p)
      SAFE_DEALLOCATE_A(rho_p)


    case(SPECIES_CHARGE_DENSITY, SPECIES_JELLIUM_CHARGE_DENSITY)

      if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
        call volume_init(volume)
        call volume_read_from_block(volume, namespace, trim(species_rho_string(species)))
      end if
       
      latt_iter = lattice_iterator_t(mesh%sb%latt, maxval(norm2(mesh%sb%latt%rlattice, dim=1)))

      rho = M_ZERO
      do icell = 1, latt_iter%n_cells
        yy = latt_iter%get(icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos, coords = xx)
          xx = xx + yy
          rr = norm2(xx)

          rerho = M_ZERO
          if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
            if(volume_in_volume(space, volume, xx)) rerho = M_ONE
          else
            call parse_expression(rerho, imrho1, space%dim, xx, rr, M_ZERO, trim(species_rho_string(species)))
          end if
          rho(ip) = rho(ip) - rerho
        end do
      end do

      if(species_type(species) == SPECIES_JELLIUM_CHARGE_DENSITY) then
         call volume_end(volume)
      end if

      rr = species_zval(species) / abs(dmf_integrate(mesh, rho(:)))
      rho(1:mesh%np) = rr * rho(1:mesh%np)

    end select

    call profiling_out(prof)
    POP_SUB(species_get_long_range_density)
  end subroutine species_get_long_range_density


  ! ---------------------------------------------------------
  subroutine func(xin, ff, jacobian)
    FLOAT, intent(in)  :: xin(:)
    FLOAT, intent(out) :: ff(:), jacobian(:,:)

    FLOAT, allocatable :: xrho(:)
    integer :: idir, jdir, dim

    PUSH_SUB(func)

    dim = mesh_p%sb%dim

    call getrho(dim, xin(1:dim+1))
    SAFE_ALLOCATE(xrho(1:mesh_p%np))

    ! First, we calculate the function ff.
    do idir = 1, dim
      xrho(1:mesh_p%np) = rho_p(1:mesh_p%np) * mesh_p%x(1:mesh_p%np, idir)
      ff(idir) = dmf_integrate(mesh_p, xrho) - pos_p(idir)
    end do
    ff(dim+1) = dmf_integrate(mesh_p, rho_p) - M_ONE

    ! Now the jacobian.
    do idir = 1, dim
      do jdir = 1, dim+1
        xrho(1:mesh_p%np) = grho_p(1:mesh_p%np, jdir) * mesh_p%x(1:mesh_p%np, idir)
        jacobian(idir, jdir) = dmf_integrate(mesh_p, xrho)
      end do
    end do
    do jdir = 1, dim+1
      xrho(1:mesh_p%np) = grho_p(1:mesh_p%np, jdir)
      jacobian(dim+1, jdir) = dmf_integrate(mesh_p, xrho)
    end do

    SAFE_DEALLOCATE_A(xrho)
    POP_SUB(func)
  end subroutine func

  ! ---------------------------------------------------------
  subroutine species_get_nlcc(species, space, pos, mesh, rho_core, accumulate)
    type(species_t), target, intent(in)    :: species
    type(space_t),           intent(in)    :: space
    FLOAT,                   intent(in)    :: pos(1:space%dim)
    type(mesh_t),            intent(in)    :: mesh
    FLOAT,                   intent(inout) :: rho_core(:)
    logical, optional,       intent(in)    :: accumulate

    FLOAT :: center(space%dim), rr
    integer :: icell, ip
    type(lattice_iterator_t) :: latt_iter
    type(ps_t), pointer :: ps

    PUSH_SUB(species_get_nlcc)

    ! only for 3D pseudopotentials, please
    if(species_is_ps(species)) then
      ps => species_ps(species)
      if(.not. optional_default(accumulate, .false.)) rho_core = M_ZERO

      latt_iter = lattice_iterator_t(mesh%sb%latt, spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold))
      do icell = 1, latt_iter%n_cells
        center = pos + latt_iter%get(icell)
        do ip = 1, mesh%np
          rr = norm2(mesh%x(ip, 1:space%dim) - center)
          if(rr < spline_range_max(ps%core)) then
            rho_core(ip) = rho_core(ip) + spline_eval(ps%core, rr)
          end if
        end do
      end do
    else
      if(.not. optional_default(accumulate, .false.)) rho_core = M_ZERO
    end if

    POP_SUB(species_get_nlcc)
  end subroutine species_get_nlcc

  ! ---------------------------------------------------------
  subroutine species_get_nlcc_grad(species, space, pos, mesh, rho_core_grad)
    type(species_t), target, intent(in)  :: species
    type(space_t),           intent(in)  :: space
    FLOAT,                   intent(in)  :: pos(1:space%dim)
    type(mesh_t),            intent(in)  :: mesh
    FLOAT,                   intent(out) :: rho_core_grad(:,:)

    FLOAT :: center(space%dim), rr, spline
    integer :: icell, ip, idir
    type(lattice_iterator_t) :: latt_iter
    type(ps_t), pointer :: ps

    PUSH_SUB(species_get_nlcc_grad)

    ! only for 3D pseudopotentials, please
    if(species_is_ps(species)) then
      ps => species_ps(species)
      rho_core_grad = M_ZERO
      if(ps_has_nlcc(ps)) then

        latt_iter = lattice_iterator_t(mesh%sb%latt, spline_cutoff_radius(ps%core_der, ps%projectors_sphere_threshold))
        do icell = 1, latt_iter%n_cells
          center = pos + latt_iter%get(icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = center)
            rr = max(rr, R_SMALL)
            if(rr >= spline_range_max(ps%core_der)) cycle
            spline = spline_eval(ps%core_der, rr)

            do idir = 1, space%dim
              rho_core_grad(ip, idir) = rho_core_grad(ip, idir) - spline*(mesh%x(ip, idir)-center(idir))/rr
            end do
          end do
        end do
      end if
    else
      rho_core_grad = M_ZERO
    end if

    POP_SUB(species_get_nlcc_grad)
  end subroutine species_get_nlcc_grad

  ! ---------------------------------------------------------
  subroutine getrho(dim, xin)
    integer, intent(in) :: dim
    FLOAT,   intent(in) :: xin(1:dim+1)

    integer :: ip, idir, idx(dim)
    FLOAT   :: r, chi(dim)

    PUSH_SUB(getrho)

    rho_p = M_ZERO
    do ip = 1, mesh_p%np
      call mesh_local_index_to_coords(mesh_p, ip, idx)
      chi(1:dim) = idx(1:dim) * mesh_p%spacing(1:dim)

      r = norm2(chi - xin(1:dim))

      if( (r/alpha_p)**2 < CNST(10.0)) then
        grho_p(ip, dim+1) = exp(-(r/alpha_p)**2)
        rho_p(ip)         = xin(dim+1) * grho_p(ip, dim+1)
      else
        grho_p(ip, dim+1) = M_ZERO
        rho_p(ip)         = M_ZERO
      end if

      do idir = 1, dim
        grho_p(ip, idir) = (M_TWO/alpha_p**2) * (chi(idir) - xin(idir)) * rho_p(ip)
      end do
    end do

    POP_SUB(getrho)
  end subroutine getrho 


  ! ---------------------------------------------------------
  !> used when the density is not available, or otherwise the Poisson eqn would be used instead
  subroutine species_get_local(species, space, mesh, namespace, x_atom, vl)
    type(species_t), target, intent(in)  :: species
    type(space_t),           intent(in)  :: space
    type(mesh_t),            intent(in)  :: mesh
    type(namespace_t),       intent(in)  :: namespace
    FLOAT,                   intent(in)  :: x_atom(:)
    FLOAT,                   intent(out) :: vl(:)

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(space%dim), x_atom_per(space%dim), r, r2, threshold
    integer :: ip, err, icell
    type(ps_t), pointer :: ps
    CMPLX :: zpot
    type(lattice_iterator_t) :: latt_iter

    type(profile_t), save :: prof

    PUSH_SUB(species_get_local)

    call profiling_in(prof, "SPECIES_GET_LOCAL")

      select case(species_type(species))

      case(SPECIES_SOFT_COULOMB)

        call parse_variable(namespace, 'SpeciesProjectorSphereThreshold', CNST(0.001), threshold)

        !Assuming that we want to take the contribution from all replica that contributes up to 0.001
        ! to the center of the cell, we arrive to a range of 1000 a.u.. 
        latt_iter = lattice_iterator_t(mesh%sb%latt, species_zval(species) / threshold)
        vl = M_ZERO
        do icell = 1, latt_iter%n_cells
          x_atom_per = x_atom(1:space%dim) + latt_iter%get(icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, r, origin = x_atom_per)
            r2 = r*r
            vl(ip) = vl(ip) -species_zval(species)/sqrt(r2+species_sc_alpha(species))
          end do
        end do

      case(SPECIES_USDEF)
        !TODO: we should control the value of 5 by a variable.
        latt_iter = lattice_iterator_t(mesh%sb%latt, CNST(5.0) * maxval(norm2(mesh%sb%latt%rlattice, dim=1)))
        vl = M_ZERO
        do icell = 1, latt_iter%n_cells
          x_atom_per = x_atom(1:space%dim) + latt_iter%get(icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, r, origin = x_atom_per, coords = xx)

            ! Note that as the spec%user_def is in input units, we have to convert
            ! the units back and forth
            xx = units_from_atomic(units_inp%length, xx)
            r = units_from_atomic(units_inp%length, r)
            zpot = species_userdef_pot(species, space%dim, xx, r)
            vl(ip) = vl(ip) + units_to_atomic(units_inp%energy, TOFLOAT(zpot))
          end do
        end do

      case(SPECIES_FROM_FILE)

        ASSERT(.not. space%is_periodic())

        call dio_function_input(trim(species_filename(species)), namespace, space, mesh, vl, err)
        if(err /= 0) then
          write(message(1), '(a)')    'Error loading file '//trim(species_filename(species))//'.'
          write(message(2), '(a,i4)') 'Error code returned = ', err
          call messages_fatal(2, namespace=namespace)
        end if

      case(SPECIES_JELLIUM)

        ASSERT(.not. space%is_periodic())

        a1 = species_z(species)/(M_TWO*species_jradius(species)**3)
        a2 = species_z(species)/species_jradius(species)
        Rb2= species_jradius(species)**2
        
        do ip = 1, mesh%np
          
          xx = mesh%x(ip, :) - x_atom(1:space%dim)
          r = norm2(xx)
          
          if(r <= species_jradius(species)) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = -species_z(species)/r
          end if
          
        end do
      
      case(SPECIES_JELLIUM_SLAB)
        a1 = M_TWO *M_PI * species_z(species)/ (M_FOUR *mesh%sb%lsize(1) *mesh%sb%lsize(2) )

        do ip = 1, mesh%np

          r = abs( mesh%x(ip, 3 ) )

          if(r <= species_jthick(species)/M_TWO ) then
            vl(ip) = a1 *( r*r/species_jthick(species) + species_jthick(species)/M_FOUR )
          else
            vl(ip) = a1 *r
          end if

        end do

      case(SPECIES_PSEUDO, SPECIES_PSPIO)
       
        ASSERT(.not. space%is_periodic())

        ps => species_ps(species)

        do ip = 1, mesh%np
          r2 = sum((mesh%x(ip, 1:space%dim) - x_atom(1:space%dim))**2)
          if(r2 < spline_range_max(ps%vlr_sq)) then
            vl(ip) = spline_eval(ps%vlr_sq, r2)
          else
            vl(ip) = P_PROTON_CHARGE*species_zval(species)/sqrt(r2)
          end if
        end do

        nullify(ps)
        
      case(SPECIES_FULL_DELTA, SPECIES_FULL_GAUSSIAN, SPECIES_CHARGE_DENSITY, SPECIES_JELLIUM_CHARGE_DENSITY)
        vl(1:mesh%np) = M_ZERO
        
      end select

      call profiling_out(prof)
    POP_SUB(species_get_local)
  end subroutine species_get_local

end module species_pot_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
