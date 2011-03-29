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
!! $Id: grid.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module species_pot_m
  use curvilinear_m
  use datasets_m
  use double_grid_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use io_function_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use poisson_m
  use profiling_m
  use ps_m
  use root_solver_m
  use simul_box_m
  use solids_m
  use species_m
  use splines_m
  use submesh_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                      &
    guess_density,               &
    species_get_density,         &
    species_get_nlcc,            &
    species_get_orbital,         &
    species_get_orbital_submesh, &
    species_get_local

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 77

  type(mesh_t), pointer :: mesh_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT :: pos_p(MAX_DIM)

contains


  ! ---------------------------------------------------------
  subroutine atom_density(mesh, sb, atom, spin_channels, rho)
    type(mesh_t),      intent(in)    :: mesh
    type(simul_box_t), intent(in)    :: sb
    type(atom_t),      intent(in)    :: atom
    integer,           intent(in)    :: spin_channels
    FLOAT,             intent(inout) :: rho(:, :) ! (mesh%np, spin_channels)

    integer :: isp, ip, in_points, nn, icell
    FLOAT :: rr, x, pos(1:MAX_DIM)
    FLOAT :: psi1, psi2, xx(MAX_DIM), yy(MAX_DIM), rerho, imrho
    type(species_t), pointer :: spec
    type(ps_t), pointer :: ps

#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(periodic_copy_t) :: pp

    PUSH_SUB(atom_density)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    spec => atom%spec
    rho = M_ZERO

    ! build density ...
    select case (species_type(spec))
    case (SPEC_FROM_FILE, SPEC_USDEF, SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_PS_CPI, SPEC_PS_FHI) ! ... from userdef
      do isp = 1, spin_channels
        rho(1:mesh%np, isp) = M_ONE
        x = (species_zval(spec)/real(spin_channels, REAL_PRECISION)) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPEC_CHARGE_DENSITY)
      ! We put, for the electron density, the same as the positive density that 
      ! creates the external potential.

      call periodic_copy_init(pp, sb, spread(M_ZERO, dim=1, ncopies = sb%dim), &
        range = M_TWO * maxval(sb%lsize(1:sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy = periodic_copy_position(pp, sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = atom%x, coords = xx)
          xx(1:sb%dim) = xx(1:sb%dim) + yy(1:sb%dim)
          rr = sqrt(dot_product(xx(1:sb%dim), xx(1:sb%dim)))
          call parse_expression(rerho, imrho, sb%dim, xx, rr, M_ZERO, trim(species_rho_string(spec)))
          rho(ip, 1) = rho(ip, 1) + rerho
        end do
      end do
      call periodic_copy_end(pp)
      if(spin_channels > 1) then
        rho(:, 1) = M_HALF*rho(:, 1)
        rho(:, 2) = rho(:, 1)
      end if

    case (SPEC_POINT, SPEC_JELLI) ! ... from jellium
      in_points = 0
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin = atom%x)
        if(rr <= species_jradius(spec)) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, mesh%vp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_m module.
 
        if (mesh%use_curvilinear) then
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(spec)) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(ip)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        else
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = atom%x)
            if(rr <= species_jradius(spec)) then
              rho(ip, 1:spin_channels) = species_zval(spec) /   &
                (mesh%vol_pp(1)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        end if
      end if

    case (SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_UPF) ! ...from pseudopotential
      ! the outer loop sums densities over atoms in neighbour cells

      ps => species_ps(spec)

      call periodic_copy_init(pp, sb, atom%x, &
        range = spline_cutoff_radius(ps%Ur(1, 1), ps%projectors_sphere_threshold))

      do icell = 1, periodic_copy_num(pp)
        pos = periodic_copy_position(pp, sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos)
          rr = max(rr, r_small)
          do nn = 1, ps%conf%p
            select case(spin_channels)
            case(1)
              psi1 = spline_eval(ps%Ur(nn, 1), rr)
              rho(ip, 1) = rho(ip, 1) + ps%conf%occ(nn, 1)*psi1*psi1 /(M_FOUR*M_PI)
            case(2)
              psi1 = spline_eval(ps%Ur(nn, 1), rr)
              psi2 = spline_eval(ps%Ur(nn, 2), rr)
              rho(ip, 1) = rho(ip, 1) + ps%conf%occ(nn, 1)*psi1*psi1 /(M_FOUR*M_PI)
              rho(ip, 2) = rho(ip, 2) + ps%conf%occ(nn, 2)*psi2*psi2 /(M_FOUR*M_PI)
            end select
          end do
        end do
      end do
  
      call periodic_copy_end(pp)
      nullify(ps)

    end select

    POP_SUB(atom_density)
  end subroutine atom_density
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  !> builds a density which is the sum of the atomic densities
  subroutine guess_density(mesh, sb, geo, qtot, nspin, spin_channels, rho)
    type(mesh_t),      intent(in)  :: mesh
    type(simul_box_t), intent(in)  :: sb
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(in)  :: qtot  ! the total charge of the system
    integer,           intent(in)  :: nspin, spin_channels
    FLOAT,             intent(out) :: rho(:, :)

    integer :: ia, is, idir, gmd_opt
    integer, save :: iseed = 321
    type(block_t) :: blk
    FLOAT :: rr, rnd, phi, theta, mag(1:3), lmag, n1, n2
    FLOAT, allocatable :: atom_rho(:,:)
    logical :: parallelized_in_atoms


    PUSH_SUB(guess_density)

    parallelized_in_atoms = .false.

    if (spin_channels == 1) then
      gmd_opt = INITRHO_PARAMAGNETIC
    else
      !%Variable GuessMagnetDensity
      !%Type integer
      !%Default ferromagnetic
      !%Section SCF
      !%Description
      !% The guess density for the SCF cycle is just the sum of all the atomic densities.
      !% When performing spin-polarized or non-collinear-spin calculations this option sets 
      !% the guess magnetization density.
      !%
      !% For anti-ferromagnetic configurations, the <tt>user_defined</tt> option should be used.
      !%
      !% Note that if the <tt>paramagnetic</tt> option is used, the final ground state will also be
      !% paramagnetic, but the same is not true for the other options.
      !%Option paramagnetic 1
      !% Magnetization density is zero.
      !%Option ferromagnetic 2
      !% Magnetization density is the sum of the atomic magnetization densities.
      !%Option random 3
      !% Each atomic magnetization density is randomly rotated.
      !%Option user_defined 77
      !% The atomic magnetization densities are rotated so that the magnetization 
      !% vector has the same direction as a vector provided by the user. In this case,
      !% the <tt>AtomsMagnetDirection</tt> block has to be set.
      !%End
      call parse_integer(datasets_check('GuessMagnetDensity'), INITRHO_FERROMAGNETIC, gmd_opt)
      if(.not.varinfo_valid_option('GuessMagnetDensity', gmd_opt)) call input_error('GuessMagnetDensity')
      call messages_print_var_option(stdout, 'GuessMagnetDensity', gmd_opt)
    end if

    rho = M_ZERO
    select case (gmd_opt)
    case (INITRHO_PARAMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:mesh%np, 1:1))

      parallelized_in_atoms = .true.

      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call atom_density(mesh, sb, geo%atom(ia), 1, atom_rho)
        rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + atom_rho(1:mesh%np, 1)
      end do

      if (spin_channels == 2) then
        rho(1:mesh%np, 1) = M_HALF*rho(1:mesh%np, 1)
        rho(1:mesh%np, 2) = rho(1:mesh%np, 1)
      end if

    case (INITRHO_FERROMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:mesh%np, 1:2))

      parallelized_in_atoms = .true.

      atom_rho = M_ZERO
      rho = M_ZERO
      do ia = geo%atoms_dist%start, geo%atoms_dist%end
        call atom_density(mesh, sb, geo%atom(ia), 2, atom_rho(1:mesh%np, 1:2))
        rho(1:mesh%np, 1:2) = rho(1:mesh%np, 1:2) + atom_rho(1:mesh%np, 1:2)
      end do

    case (INITRHO_RANDOM) ! Randomly oriented spins
      SAFE_ALLOCATE(atom_rho(1:mesh%np, 1:2))
      do ia = 1, geo%natoms
        call atom_density(mesh, sb, geo%atom(ia), 2, atom_rho)

        if (nspin == 2) then
          call quickrnd(iseed, rnd)
          rnd = rnd - M_HALF
          if (rnd > M_ZERO) then
            rho(1:mesh%np, 1:2) = rho(1:mesh%np, 1:2) + atom_rho(1:mesh%np, 1:2)
          else
            rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + atom_rho(1:mesh%np, 2)
            rho(1:mesh%np, 2) = rho(1:mesh%np, 2) + atom_rho(1:mesh%np, 1)
          end if
        elseif (nspin == 4) then
          call quickrnd(iseed, phi)
          call quickrnd(iseed, theta)
          phi = phi*M_TWO*M_PI
          theta = theta*M_PI
          rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:mesh%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:mesh%np, 2)
          rho(1:mesh%np, 2) = rho(1:mesh%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:mesh%np, 1) &
            + cos(theta/M_TWO)**2*atom_rho(1:mesh%np, 2)
          rho(1:mesh%np, 3) = rho(1:mesh%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:mesh%np, 1) - atom_rho(1:mesh%np, 2))
          rho(1:mesh%np, 4) = rho(1:mesh%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:mesh%np, 1) - atom_rho(1:mesh%np, 2))
        end if
      end do

    case (INITRHO_USERDEF) ! User-defined
      
      !%Variable AtomsMagnetDirection
      !%Type block
      !%Section Hamiltonian
      !%Description
      !% This option is only used when <tt>GuessMagnetDensity</tt> is
      !% set to <tt>user_defined</tt>. It provides a direction for the
      !% magnetization vector of each atom when building the guess
      !% density. In order to do that, the user should specify the
      !% coordinates of a vector that has the desired direction and
      !% norm.  Note that it is necessary to maintain the ordering in
      !% which the species were defined in the coordinates
      !% specifications.
      !%
      !% For spin-polarized calculations, the vectors should have only
      !% one component; for non-collinear-spin calculations, they
      !% should have three components.
      !%End
      if(parse_block(datasets_check('AtomsMagnetDirection'), blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined."
        call messages_fatal(1)
      end if

      if (parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows."
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(atom_rho(1:mesh%np, 1:2))
      do ia = 1, geo%natoms
        !Read from AtomsMagnetDirection block 
        if (nspin == 2) then
          call parse_block_float(blk, ia-1, 0, mag(1))
          lmag = abs(mag(1))
        elseif (nspin == 4) then
          do idir = 1, 3
            call parse_block_float(blk, ia-1, idir-1, mag(idir))
            if (abs(mag(idir)) < CNST(1.0e-20)) mag(idir) = M_ZERO
          end do
          lmag = sqrt(dot_product(mag(1:3), mag(1:3)))
        end if

        !Get atomic density
        call atom_density(mesh, sb, geo%atom(ia), 2, atom_rho)

        !Scale magnetization density
        n1 = dmf_integrate(mesh, atom_rho(:, 1))
        n2 = dmf_integrate(mesh, atom_rho(:, 2))
        if (lmag > n1 + n2) then
          mag = mag*(n1 + n2)/lmag
          lmag = n1 + n2
        elseif (lmag == M_ZERO) then
          if (n1 - n2 == M_ZERO) then
            rho(1:mesh%np, 1:2) = rho(1:mesh%np, 1:2) + atom_rho(1:mesh%np, 1:2)
          else
            atom_rho(:, 1) = (atom_rho(:, 1) + atom_rho(:, 2))/M_TWO
            rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + atom_rho(1:mesh%np, 1)
            rho(1:mesh%np, 2) = rho(1:mesh%np, 2) + atom_rho(1:mesh%np, 1)
          end if
          cycle
        end if
        if (n1 - n2 /= lmag .and. n2 /= M_ZERO) then
          if (n1 - n2 < lmag) then
            atom_rho(:, 1) = atom_rho(:, 1) + (lmag - n1 + n2)/M_TWO/n2*atom_rho(:, 2)
            atom_rho(:, 2) = (n1 + n2 - lmag)/M_TWO/n2*atom_rho(:, 2)
          elseif (n1 - n2 > lmag) then
            atom_rho(:, 2) = atom_rho(:, 2) + (n1 - n2 - lmag)/M_TWO/n1*atom_rho(:, 1)
            atom_rho(:, 1) = (lmag + n1 + n2)/M_TWO/n1*atom_rho(:, 1)
          end if
        end if

        !Rotate magnetization density
        if (nspin == 2) then
          if (mag(1) > M_ZERO) then
            rho(1:mesh%np, 1:2) = rho(1:mesh%np, 1:2) + atom_rho(1:mesh%np, 1:2)
          else
            rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + atom_rho(1:mesh%np, 2)
            rho(1:mesh%np, 2) = rho(1:mesh%np, 2) + atom_rho(1:mesh%np, 1)
          end if

        elseif (nspin == 4) then
          theta = acos(mag(3)/lmag)
          if (mag(1) == M_ZERO) then
            if (mag(2) == M_ZERO) then
              phi = M_ZERO
            elseif (mag(2) < M_ZERO) then
              phi = M_PI*M_TWOTHIRD
            elseif (mag(2) > M_ZERO) then
              phi = M_PI*M_HALF
            end if
          else
            if (mag(2) < M_ZERO) then
              phi = M_TWO*M_PI - acos(mag(1)/sin(theta)/lmag)
            elseif (mag(2) >= M_ZERO) then
              phi = acos(mag(1)/sin(theta)/lmag)
            end if
          end if

          rho(1:mesh%np, 1) = rho(1:mesh%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:mesh%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:mesh%np, 2)
          rho(1:mesh%np, 2) = rho(1:mesh%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:mesh%np, 1) &
               + cos(theta/M_TWO)**2*atom_rho(1:mesh%np, 2)
          rho(1:mesh%np, 3) = rho(1:mesh%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:mesh%np, 1) - atom_rho(1:mesh%np, 2))
          rho(1:mesh%np, 4) = rho(1:mesh%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:mesh%np, 1) - atom_rho(1:mesh%np, 2))
        end if
      end do

      call parse_block_end(blk)

    end select


#ifdef HAVE_MPI
    if(geo%atoms_dist%parallel .and. parallelized_in_atoms) then
      do is = 1, spin_channels
        atom_rho(1:mesh%np, 1) = rho(1:mesh%np, is)
        call MPI_Allreduce(atom_rho(1, 1), rho(1, is), mesh%np, MPI_FLOAT, MPI_SUM, geo%atoms_dist%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    ! we now renormalize the density (necessary if we have a charged system)
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(mesh, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', rr
    call messages_info(1)

    rr = qtot / rr
    rho = rr * rho
    rr = M_ZERO
    do is = 1, spin_channels
      rr = rr + dmf_integrate(mesh, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', rr
    call messages_info(1)

    SAFE_DEALLOCATE_A(atom_rho)
    POP_SUB(guess_density)
  end subroutine guess_density


  ! ---------------------------------------------------------
  subroutine species_get_density(spec, pos, mesh, geo, rho)
    type(species_t),            intent(in)  :: spec
    FLOAT,                      intent(in)  :: pos(:)
    type(mesh_t),       target, intent(in)  :: mesh
    type(geometry_t),           intent(in)  :: geo
    FLOAT,                      intent(out) :: rho(:)

    type(root_solver_t) :: rs
    logical :: conv
    integer :: dim
    FLOAT   :: x(1:MAX_DIM+1), chi0(MAX_DIM), startval(MAX_DIM + 1)
    FLOAT   :: delta, alpha, beta, xx(MAX_DIM), yy(MAX_DIM), rr, imrho, rerho
    FLOAT   :: dist2, dist2_min
    integer :: icell, ipos, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps
    logical :: have_point
#ifdef HAVE_MPI
    real(8) :: local_min(2), global_min(2)
#endif

    PUSH_SUB(species_get_density)

    select case(species_type(spec))

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
      ps => species_ps(spec)
      rho = M_ZERO
      call periodic_copy_init(pp, mesh%sb, pos, range = spline_cutoff_radius(ps%nlr, ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        call dmf_put_radial_spline(mesh, ps%nlr, periodic_copy_position(pp, mesh%sb, icell), rho, add = .true.)
      end do
      call periodic_copy_end(pp)
      nullify(ps)

    case(SPEC_FULL_DELTA)

      dist2_min = huge(delta)
      ipos = 0

      do ip = 1, mesh%np

        rho(ip) = M_ZERO

        dist2 = sum((mesh%x(ip, 1:mesh%sb%dim) - pos(1:mesh%sb%dim))**2)
        if (dist2 < dist2_min) then
          ipos = ip
          dist2_min = dist2
        end if

      end do

      write(message(1), '(3a,f5.2,3a)') &
        "Info: spec_full_delta species ", trim(species_label(spec)), &
        " atom displaced ", units_from_atomic(units_out%length, sqrt(dist2_min)), &
        " [ ", trim(units_abbrev(units_out%length)), " ]"
      call messages_info(1)

      have_point = .true.
#ifdef HAVE_MPI
      ! in parallel we have to find the minimum of the whole grid
      if(mesh%parallel_in_domains) then

        local_min = (/ dist2_min, dble(mesh%mpi_grp%rank)/)
        call MPI_Allreduce(local_min, global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)

        if(mesh%mpi_grp%rank /= nint(global_min(2))) have_point = .false.

      end if
#endif
      if(have_point) then
        if(mesh%use_curvilinear) then
          rho(ipos) = -species_z(spec)/mesh%vol_pp(ipos)
        else
          rho(ipos) = -species_z(spec)/mesh%vol_pp(1)
        end if
      end if

    case(SPEC_FULL_GAUSSIAN)

      ! periodic copies are not considered in this routine
      if(simul_box_is_periodic(mesh%sb)) then
        call messages_experimental("spec_full_gaussian for periodic systems")
      end if

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      dim = mesh%sb%dim

      SAFE_ALLOCATE(rho_p(1:mesh%np))
      SAFE_ALLOCATE(grho_p(1:mesh%np, 1:dim+1))

      mesh_p => mesh
      pos_p = pos

      ! Initial guess.
      call curvilinear_x2chi(mesh%sb, mesh%cv, pos, chi0)
      delta   = mesh%spacing(1)
      alpha   = sqrt(M_TWO)*species_sigma(spec)*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:dim) = CNST(1.0)

      ! the dim+1 variable is the normalization of the delta function
      startval(dim+1) = beta

      ! get a better estimate for beta
      call getrho(startval)
      beta = M_ONE / dmf_integrate(mesh, rho_p)
      startval(dim+1) = beta

      ! solve equation
      call root_solver_init(rs, dim+1, &
        solver_type=ROOT_NEWTON, maxiter=500, abs_tolerance=CNST(1.0e-10))
      call droot_solver_run(rs, func, x, conv, startval=startval)

      if(.not.conv) then
        write(message(1),'(a)') 'Internal error in get_species_density.'
        call messages_fatal(1)
      end if

      ! we want a charge of -Z
      rho = -species_z(spec)*rho_p

      nullify(mesh_p)
      SAFE_DEALLOCATE_A(grho_p)
      SAFE_DEALLOCATE_A(rho_p)


    case(SPEC_CHARGE_DENSITY)

      call periodic_copy_init(pp, mesh%sb, spread(M_ZERO, dim=1, ncopies = mesh%sb%dim), &
        range = M_TWO * maxval(mesh%sb%lsize(1:mesh%sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy = periodic_copy_position(pp, mesh%sb, icell)
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, origin = pos, coords = xx)
          xx(1:mesh%sb%dim) = xx(1:mesh%sb%dim) + yy(1:mesh%sb%dim)
          rr = sqrt(dot_product(xx(1:mesh%sb%dim), xx(1:mesh%sb%dim)))
          call parse_expression(rerho, imrho, mesh%sb%dim, xx, rr, M_ZERO, trim(species_rho_string(spec)))
          rho(ip) = rho(ip) - rerho
        end do
      end do

      call periodic_copy_end(pp)

    end select

    POP_SUB(species_get_density)
  end subroutine species_get_density


  ! ---------------------------------------------------------
  subroutine func(xin, ff, jacobian)
    FLOAT, intent(in)  :: xin(:)
    FLOAT, intent(out) :: ff(:), jacobian(:,:)

    FLOAT, allocatable :: xrho(:)
    integer :: idir, jdir, dim

    PUSH_SUB(func)

    dim = mesh_p%sb%dim

    call getrho(xin)
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
  subroutine species_get_nlcc(spec, pos, mesh, geo, rho_core)
    type(species_t),  intent(in)  :: spec
    FLOAT,            intent(in)  :: pos(MAX_DIM)
    type(mesh_t),     intent(in)  :: mesh
    type(geometry_t), intent(in)  :: geo
    FLOAT,            intent(out) :: rho_core(:)

    integer :: icell
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps

    PUSH_SUB(species_get_density)

    ! only for 3D pseudopotentials, please
    if(species_is_ps(spec)) then
      ps => species_ps(spec)
      rho_core = M_ZERO
      call periodic_copy_init(pp, mesh%sb, pos, range = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        call dmf_put_radial_spline(mesh, ps%core, periodic_copy_position(pp, mesh%sb, icell), rho_core, add = .true.)
      end do
      call periodic_copy_end(pp)
    else
      rho_core = M_ZERO
    end if

    POP_SUB(species_get_density)
  end subroutine species_get_nlcc

  ! ---------------------------------------------------------
  subroutine getrho(xin)
    FLOAT, intent(in) :: xin(:)

    integer :: ip, jp, idir, dim
    FLOAT   :: r, chi(MAX_DIM)

    PUSH_SUB(getrho)

    dim = mesh_p%sb%dim
    rho_p = M_ZERO
    do ip = 1, mesh_p%np

      jp = ip
      if(mesh_p%parallel_in_domains) &
        jp = mesh_p%vp%local(mesh_p%vp%xlocal(mesh_p%vp%partno)+ip-1)

      chi(1:dim) = mesh_p%idx%lxyz(jp, 1:dim) * mesh_p%spacing(1:dim) + mesh_p%sb%box_offset(1:dim) 

      r = sqrt( sum( (chi(1:dim) - xin(1:dim))**2 ) )

      if( (r/alpha_p)**2 < CNST(10.0)) then
        grho_p(ip, dim+1) = exp(-(r/alpha_p)**2)
        rho_p(ip)         = xin(dim+1) * grho_p(ip, dim+1)
      else
        grho_p(ip, dim+1) = M_ZERO
        rho_p(ip)         = M_ZERO
      end if

      do idir = 1, dim
        grho_p(ip, idir) = (M_TWO/alpha_p**2) * (chi(idir)-xin(idir)) * rho_p(ip)
      end do
    end do

    POP_SUB(getrho)
  end subroutine getrho 


  ! ---------------------------------------------------------
  subroutine species_get_local(spec, mesh, x_atom, vl, time)
    type(species_t), intent(in) :: spec
    type(mesh_t),    intent(in) :: mesh
    FLOAT,           intent(in) :: x_atom(:)
    FLOAT,           intent(out):: vl(:)
    FLOAT, optional, intent(in) :: time

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, time_
    integer :: ip, err, idim
    type(ps_t), pointer :: ps

    type(profile_t), save :: prof

    PUSH_SUB(species_get_local)
    call profiling_in(prof, "SPECIES_GET_LOCAL")

    time_ = M_ZERO

    if (present(time)) time_ = time

      select case(species_type(spec))
      case(SPEC_USDEF)

        do ip = 1, mesh%np
          
          xx = M_ZERO
          xx(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim) - x_atom(1:mesh%sb%dim)
          r = sqrt(sum(xx(1:mesh%sb%dim)**2))
          
          ! Note that as the spec%user_def is in input units, we have to convert
          ! the units back and forth
          forall(idim = 1:mesh%sb%dim) xx(idim) = units_from_atomic(units_inp%length, xx(idim))
          r = units_from_atomic(units_inp%length, r)
          vl(ip) = units_to_atomic(units_inp%energy, species_userdef_pot(spec, mesh%sb%dim, xx, r, time_))

        end do


      case(SPEC_FROM_FILE)

        call dio_function_input(trim(species_filename(spec)), mesh, vl, err)
        if(err .ne. 0) then
          write(message(1), '(a)')    'File '//trim(species_filename(spec))//'not found.'
          write(message(2), '(a,i4)') 'Error code returned = ', err
          call messages_fatal(1)
        end if

      case(SPEC_POINT, SPEC_JELLI)
        a1 = species_z(spec)/(M_TWO*species_jradius(spec)**3)
        a2 = species_z(spec)/species_jradius(spec)
        Rb2= species_jradius(spec)**2
        
        do ip = 1, mesh%np
          
          xx(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim) - x_atom(1:mesh%sb%dim)
          r = sqrt(sum(xx(1:mesh%sb%dim)**2))
          
          if(r <= species_jradius(spec)) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = -species_z(spec)/r
          end if
          
        end do
        
      case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
        ps => species_ps(spec)
        do ip = 1, mesh%np
          vl(ip) = sum((mesh%x(ip, 1:mesh%sb%dim) - x_atom(1:mesh%sb%dim))**2)
        end do
        call spline_eval_vec(ps%vlr_sq, mesh%np, vl)
        nullify(ps)
        
      case(SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_CHARGE_DENSITY)
        vl(1:mesh%np) = M_ZERO
        
      end select

      call profiling_out(prof)
    POP_SUB(species_get_local)
  end subroutine species_get_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Places, in the function phi (defined in each point of the mesh), the
  !! iorb-th atomic orbital. The orbitals are obtained from the species data
  !! type, and are numbered from one to species_niwfs(spec). It may happen
  !! that there are different orbitals for each spin-polarization direction,
  !! and therefore the orbital is also characterized by the label "is".
  !!
  !! In order to put the orbital in the mesh, it is necessary to know where
  !! the species is, and this is given by the vector "pos".
  !!
  !! \todo Most of this work should be done inside the species
  !! module, and we should get rid of species_iwf_i, species_ifw_l, etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine species_get_orbital(spec, mesh, iorb, ispin, pos, phi)
    type(species_t),   intent(in)     :: spec
    type(mesh_t),      intent(in)     :: mesh
    integer,           intent(in)     :: iorb
    integer,           intent(in)     :: ispin   !< The spin index.
    FLOAT,             intent(in)     :: pos(:)  !< The position of the atom.
    FLOAT,             intent(out)    :: phi(:)  !< The function defined in the mesh where the orbitals is returned.

    integer :: i, l, m, ip
    FLOAT :: r2, x(1:MAX_DIM)
    FLOAT, allocatable :: xf(:, :), ylm(:), sphi(:)
    type(ps_t), pointer :: ps
    type(submesh_t) :: sphere

    PUSH_SUB(species_get_orbital)

    call species_iwf_ilm(spec, iorb, ispin, i, l, m)

    if(species_is_ps(spec)) then

      ps => species_ps(spec)
      SAFE_ALLOCATE(xf(1:mesh%np, 1:mesh%sb%dim))
      SAFE_ALLOCATE(ylm(1:mesh%np))

      do ip = 1, mesh%np
        x(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim) - pos(1:mesh%sb%dim)
        phi(ip) = sum(x(1:mesh%sb%dim)**2)
        xf(ip, 1:mesh%sb%dim) = x(1:mesh%sb%dim)
      end do

      call spline_eval_vec(ps%ur_sq(i, ispin), mesh%np, phi)
      call loct_ylm(mesh%np, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))

      do ip = 1, mesh%np
        phi(ip) = phi(ip)*ylm(ip)
      end do

      SAFE_DEALLOCATE_A(xf)
      SAFE_DEALLOCATE_A(ylm)

      nullify(ps)
    else

      do ip = 1, mesh%np
        x(1:mesh%sb%dim) = mesh%x(ip, 1:mesh%sb%dim) - pos(1:mesh%sb%dim)
        r2 = sum(x(1:mesh%sb%dim)**2)
        select case(mesh%sb%dim)
        case(1)
          phi(ip) = exp(-species_omega(spec)*r2/M_TWO) * hermite(i - 1, x(1)*sqrt(species_omega(spec)))
        case(2)
          phi(ip) = exp(-species_omega(spec)*r2/M_TWO) * &
            hermite(i - 1, x(1)*sqrt(species_omega(spec))) * hermite(l - 1, x(2)*sqrt(species_omega(spec)))
        case(3)
          phi(ip) = exp(-species_omega(spec)*r2/M_TWO) * &
               hermite(i - 1, x(1) * sqrt(species_omega(spec))) * &
               hermite(l - 1, x(2) * sqrt(species_omega(spec))) * &
               hermite(m - 1, x(3) * sqrt(species_omega(spec)))
        end select
      end do
    end if

    POP_SUB(species_get_orbital)
  end subroutine species_get_orbital


  ! ---------------------------------------------------------
  subroutine species_get_orbital_submesh(spec, submesh, iorb, ispin, pos, phi, derivative)
    type(species_t),   intent(in)  :: spec       !< The species.
    type(submesh_t),   intent(in)  :: submesh    !< The submesh descriptor where the orbital will be calculated.
    integer,           intent(in)  :: iorb       !< The index of the orbital to return.
    integer,           intent(in)  :: ispin      !< The spin index.
    FLOAT,             intent(in)  :: pos(:)     !< The position of the atom.
    FLOAT,             intent(out) :: phi(:)     !< The function defined in the mesh where the orbitals is returned.
    logical, optional, intent(in)  :: derivative !< If present and .true. returns the derivative of the orbital.

    integer :: i, l, m, ip
    FLOAT :: r2, x(1:MAX_DIM), sqrtw, ww
    FLOAT, allocatable :: xf(:, :), ylm(:)
    type(ps_t), pointer :: ps
    type(spline_t) :: dur
    logical :: derivative_
    
    if(submesh%ns == 0) return

    PUSH_SUB(species_get_orbital_submesh)

    derivative_ = optional_default(derivative, .false.)

    ASSERT(ubound(phi, dim = 1) >= submesh%ns)

    call species_iwf_ilm(spec, iorb, ispin, i, l, m)

    if(species_is_ps(spec)) then
      ps => species_ps(spec)
      
      if(derivative_) then
        call spline_init(dur)
        call spline_der(ps%ur(i, ispin), dur)
      end if

      SAFE_ALLOCATE(xf(1:submesh%ns, 1:submesh%mesh%sb%dim))
      SAFE_ALLOCATE(ylm(1:submesh%ns))
      do ip = 1, submesh%ns
        x(1:submesh%mesh%sb%dim) = submesh%mesh%x(submesh%jxyz(ip), 1:submesh%mesh%sb%dim) - pos(1:submesh%mesh%sb%dim)
        phi(ip) = sum(x(1:submesh%mesh%sb%dim)**2)
        if(derivative_) phi(ip) = sqrt(phi(ip))
        xf(ip, 1:submesh%mesh%sb%dim) = x(1:submesh%mesh%sb%dim)
      end do

      if(.not. derivative_) then
        call spline_eval_vec(ps%ur_sq(i, ispin), submesh%ns, phi)
      else
        call spline_eval_vec(dur, submesh%ns, phi)
      end if

      call loct_ylm(submesh%ns, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))

      do ip = 1, submesh%ns
        phi(ip) = phi(ip)*ylm(ip)
      end do

      SAFE_DEALLOCATE_A(xf)
      SAFE_DEALLOCATE_A(ylm)

      if(derivative_) call spline_end(dur)
      nullify(ps)
    else
      
      ASSERT(.not. derivative_)

      ww = species_omega(spec)
      sqrtw = sqrt(ww)

      select case(submesh%mesh%sb%dim)
      case(1)
        do ip = 1, submesh%ns
          x(1:submesh%mesh%sb%dim) = submesh%mesh%x(submesh%jxyz(ip), 1:submesh%mesh%sb%dim) - pos(1:submesh%mesh%sb%dim)
          r2 = sum(x(1:submesh%mesh%sb%dim)**2)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)
        end do
      case(2)
        do ip = 1, submesh%ns
          x(1:submesh%mesh%sb%dim) = submesh%mesh%x(submesh%jxyz(ip), 1:submesh%mesh%sb%dim) - pos(1:submesh%mesh%sb%dim)
          r2 = sum(x(1:submesh%mesh%sb%dim)**2)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)*hermite(l - 1, x(2)*sqrtw)
        end do
      case(3)
        do ip = 1, submesh%ns
          x(1:submesh%mesh%sb%dim) = submesh%mesh%x(submesh%jxyz(ip), 1:submesh%mesh%sb%dim) - pos(1:submesh%mesh%sb%dim)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)*hermite(l - 1, x(2)*sqrtw)*hermite(m - 1, x(3)*sqrtw)
        end do
      end select
      
    end if

    POP_SUB(species_get_orbital_submesh)
  end subroutine species_get_orbital_submesh

end module species_pot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
