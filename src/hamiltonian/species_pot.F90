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
  use parser_m
  use loct_math_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use root_solver_m
  use simul_box_m
  use species_m
  use solids_m
  use splines_m
  use submesh_m
  use poisson_m
  use profiling_m
  use ps_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use io_function_m

  implicit none

  private
  public ::                      &
    guess_density,               &
    species_get_density,         &
    species_get_orbital,         &
    species_get_orbital_submesh, &
    species_get_local

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 77

  type(mesh_t),       pointer :: m_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT :: pos_p(MAX_DIM)

contains


  ! ---------------------------------------------------------
  subroutine atom_density(m, sb, atom, spin_channels, rho)
    type(mesh_t),      intent(in)    :: m
    type(simul_box_t), intent(in)    :: sb
    type(atom_t),      intent(in)    :: atom
    integer,           intent(in)    :: spin_channels
    FLOAT,             intent(inout) :: rho(:, :) ! (m%np, spin_channels)

    integer :: i, in_points, n, icell
    FLOAT :: r, x, pos(1:MAX_DIM)
    FLOAT :: psi1, psi2, xx(MAX_DIM), yy(MAX_DIM), rerho, imrho
    type(species_t), pointer :: s
    type(ps_t), pointer :: ps

#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(periodic_copy_t) :: pp

    call push_sub('species_pot.atom_density')

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    s => atom%spec
    rho = M_ZERO

    ! build density ...
    select case (species_type(s))
    case (SPEC_FROM_FILE, SPEC_USDEF, SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_PS_CPI, SPEC_PS_FHI) ! ... from userdef
      do i = 1, spin_channels
        rho(1:m%np, i) = M_ONE
        x = (species_zval(s)/real(spin_channels, REAL_PRECISION)) / dmf_integrate(m, rho(:, i))
        rho(1:m%np, i) = x * rho(1:m%np, i)
      end do

    case (SPEC_CHARGE_DENSITY)
      ! We put, for the electron density, the same as the positive density that 
      ! creates the external potential.

      call periodic_copy_init(pp, sb, spread(M_ZERO, dim=1, ncopies = MAX_DIM), &
        range = M_TWO * maxval(sb%lsize(1:sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy = periodic_copy_position(pp, sb, icell)
        do i = 1, m%np
          call mesh_r(m, i, r, x = xx, a = atom%x)
          xx(1:sb%dim) = xx(1:sb%dim) + yy(1:sb%dim)
          r = sqrt(dot_product(xx(1:sb%dim), xx(1:sb%dim)))
          call parse_expression(rerho, imrho, sb%dim, xx, r, M_ZERO, trim(species_rho_string(s)))
          rho(i, 1) = rho(i, 1) + rerho
        end do
      end do
      call periodic_copy_end(pp)
      if(spin_channels > 1) then
        rho(:, 1) = M_HALF*rho(:, 1)
        rho(:, 2) = rho(:, 1)
      end if

    case (SPEC_POINT, SPEC_JELLI) ! ... from jellium
      in_points = 0
      do i = 1, m%np
        call mesh_r(m, i, r, a=atom%x)
        if(r <= species_jradius(s)) then
          in_points = in_points + 1
        end if
      end do

#if defined(HAVE_MPI)
      if(m%parallel_in_domains) then
        call MPI_Allreduce(in_points, in_points_red, 1, MPI_INTEGER, MPI_SUM, m%vp%comm, mpi_err)
        in_points = in_points_red
      end if
#endif

      if(in_points > 0) then
        ! This probably should be done inside the mesh_function_m module.
 
        if (m%use_curvilinear) then
          do i = 1, m%np
            call mesh_r(m, i, r, a=atom%x)
            if(r <= species_jradius(s)) then
              rho(i, 1:spin_channels) = species_zval(s) /   &
                (m%vol_pp(i)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        else
          do i = 1, m%np
            call mesh_r(m, i, r, a=atom%x)
            if(r <= species_jradius(s)) then
              rho(i, 1:spin_channels) = species_zval(s) /   &
                (m%vol_pp(1)*real(in_points*spin_channels, REAL_PRECISION))
            end if
          end do
        end if
      end if

    case (SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_UPF) ! ...from pseudopotential
      ! the outer loop sums densities over atoms in neighbour cells

      ps => species_ps(s)

      call periodic_copy_init(pp, sb, atom%x, &
        range = spline_cutoff_radius(ps%Ur(1, 1), ps%projectors_sphere_threshold))

      do icell = 1, periodic_copy_num(pp)
        pos = periodic_copy_position(pp, sb, icell)
        do i = 1, m%np
          call mesh_r(m, i, r, pos)
          r = max(r, r_small)
          do n = 1, ps%conf%p
            select case(spin_channels)
            case(1)
              psi1 = spline_eval(ps%Ur(n, 1), r)
              rho(i, 1) = rho(i, 1) + ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
            case(2)
              psi1 = spline_eval(ps%Ur(n, 1), r)
              psi2 = spline_eval(ps%Ur(n, 2), r)
              rho(i, 1) = rho(i, 1) + ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
              rho(i, 2) = rho(i, 2) + ps%conf%occ(n, 2)*psi2*psi2 /(M_FOUR*M_PI)
            end select
          end do
        end do
      end do
  
      call periodic_copy_end(pp)
      nullify(ps)

    end select

    call pop_sub()
  end subroutine atom_density
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  ! builds a density which is the sum of the atomic densities
  subroutine guess_density(m, sb, geo, qtot, nspin, spin_channels, rho)
    type(mesh_t),      intent(in)  :: m
    type(simul_box_t), intent(in)  :: sb
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(in)  :: qtot  ! the total charge of the system
    integer,           intent(in)  :: nspin, spin_channels
    FLOAT,             intent(out) :: rho(:, :)

    integer :: ia, is, gmd_opt, i
    integer, save :: iseed = 321
    type(block_t) :: blk
    FLOAT :: r, rnd, phi, theta, mag(MAX_DIM), lmag, n1, n2
    FLOAT, allocatable :: atom_rho(:,:)

    call push_sub('species_pot.guess_density')

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
      !% Note that if the <tt>paramagnetic</tt> option is used the final ground-state will also be
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
      SAFE_ALLOCATE(atom_rho(1:m%np, 1:1))
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 1, atom_rho(1:m%np, 1:1))
        rho(1:m%np, 1:1) = rho(1:m%np, 1:1) + atom_rho(1:m%np, 1:1)
      end do
      if (spin_channels == 2) then
        rho(1:m%np, 1) = M_HALF*rho(1:m%np, 1)
        rho(1:m%np, 2) = rho(1:m%np, 1)
      end if

    case (INITRHO_FERROMAGNETIC)
      SAFE_ALLOCATE(atom_rho(1:m%np, 1:2))
      atom_rho = M_ZERO
      rho = M_ZERO
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 2, atom_rho(1:m%np, 1:2))
        rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
      end do

    case (INITRHO_RANDOM) ! Randomly oriented spins
      SAFE_ALLOCATE(atom_rho(1:m%np, 1:2))
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 2, atom_rho)

        if (nspin == 2) then
          call quickrnd(iseed, rnd)
          rnd = rnd - M_HALF
          if (rnd > M_ZERO) then
            rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
          else
            rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 2)
            rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
          end if
        elseif (nspin == 4) then
          call quickrnd(iseed, phi)
          call quickrnd(iseed, theta)
          phi = phi*M_TWO*M_PI
          theta = theta*M_PI
          rho(1:m%np, 1) = rho(1:m%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 2) = rho(1:m%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
            + cos(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 3) = rho(1:m%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
          rho(1:m%np, 4) = rho(1:m%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
        end if
      end do

    case (INITRHO_USERDEF) ! User-defined
      
      !%Variable AtomsMagnetDirection
      !%Type block
      !%Section Hamiltonian
      !%Description
      !% This option is only used when <tt>GuessMagnetDensity</tt> is set to 
      !% <tt>user_defined</tt>. It provides a direction for each atoms magnetization 
      !% vector when building the guess density. In order to do that the user should
      !% specify the coordinates of a vector that has the desired direction and norm.
      !% Note that it is necessary to maintain the ordering in which the species
      !% were defined in the coordinates specifications.
      !%
      !% For spin-polarized calculations the vectors should have only one component and
      !% for non-collinear-spin calculations they should have three components.
      !%End
      if(parse_block(datasets_check('AtomsMagnetDirection'), blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined "
        call write_fatal(1)
      end if

      if (parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows"
        call write_fatal(1)
      end if

      SAFE_ALLOCATE(atom_rho(1:m%np, 1:2))
      do ia = 1, geo%natoms
        !Read from AtomsMagnetDirection block 
        if (nspin == 2) then
          call parse_block_float(blk, ia-1, 0, mag(1))
          lmag = abs(mag(1))
        elseif (nspin == 4) then
          do i = 1, 3
            call parse_block_float(blk, ia-1, i-1, mag(i))
            if (abs(mag(i)) < CNST(1.0e-20)) mag(i) = M_ZERO
          end do
          lmag = sqrt(dot_product(mag, mag))
        end if

        !Get atomic density
        call atom_density(m, sb, geo%atom(ia), 2, atom_rho)

        !Scale magnetization density
        n1 = dmf_integrate(m, atom_rho(:, 1))
        n2 = dmf_integrate(m, atom_rho(:, 2))
        if (lmag > n1 + n2) then
          mag = mag*(n1 + n2)/lmag
          lmag = n1 + n2
        elseif (lmag == M_ZERO) then
          if (n1 - n2 == M_ZERO) then
            rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
          else
            atom_rho(:, 1) = (atom_rho(:, 1) + atom_rho(:, 2))/M_TWO
            rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 1)
            rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
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
            rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
          else
            rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 2)
            rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
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

          rho(1:m%np, 1) = rho(1:m%np, 1) + cos(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
            + sin(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 2) = rho(1:m%np, 2) + sin(theta/M_TWO)**2*atom_rho(1:m%np, 1) &
               + cos(theta/M_TWO)**2*atom_rho(1:m%np, 2)
          rho(1:m%np, 3) = rho(1:m%np, 3) + cos(theta/M_TWO)*sin(theta/M_TWO)*cos(phi)* &
            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
          rho(1:m%np, 4) = rho(1:m%np, 4) - cos(theta/M_TWO)*sin(theta/M_TWO)*sin(phi)* &
            (atom_rho(1:m%np, 1) - atom_rho(1:m%np, 2))
        end if
      end do

      call parse_block_end(blk)

    end select

    ! we now renormalize the density (necessary if we have a charged system)
    r = M_ZERO
    do is = 1, spin_channels
      r = r + dmf_integrate(m, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Unnormalized total charge = ', r
    call write_info(1)

    r = qtot/r
    rho = r*rho
    r = M_ZERO
    do is = 1, spin_channels
      r = r + dmf_integrate(m, rho(:, is))
    end do

    write(message(1),'(a,f13.6)')'Info: Renormalized total charge = ', r
    call write_info(1)

    SAFE_DEALLOCATE_A(atom_rho)
    call pop_sub()
  end subroutine guess_density
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine species_get_density(s, pos, gr, geo, rho)
    type(species_t),             intent(in)  :: s
    FLOAT,                      intent(in)  :: pos(MAX_DIM)
    type(grid_t),       target, intent(in)  :: gr
    type(geometry_t),           intent(in)  :: geo
    FLOAT,                      intent(out) :: rho(:)

    type(root_solver_t) :: rs
    logical :: conv
    integer :: dim, i
    FLOAT   :: x(1:MAX_DIM+1), chi0(MAX_DIM), startval(MAX_DIM + 1)
    FLOAT   :: delta, alpha, beta, xx(MAX_DIM), yy(MAX_DIM), r, imrho, rerho
    FLOAT   :: dist2, dist2_min
    integer :: icell, ipos, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps
    logical :: have_point
#ifdef HAVE_MPI
    real(8) :: local_min(2), global_min(2)
#endif

    call push_sub('species_pot.species_get_density')

    select case(species_type(s))

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
      ps => species_ps(s)
      rho = M_ZERO
      call periodic_copy_init(pp, gr%sb, pos, range = spline_cutoff_radius(ps%nlr, ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        call dmf_put_radial_spline(gr%mesh, ps%nlr, periodic_copy_position(pp, gr%sb, icell), rho, add = .true.)
      end do
      call periodic_copy_end(pp)
      nullify(ps)

    case(SPEC_FULL_DELTA)

      dist2_min = huge(delta)
      ipos = 0

      do ip = 1, gr%mesh%np

        rho(ip) = M_ZERO

        dist2 = sum((gr%mesh%x(ip, 1:MAX_DIM) - pos(1:MAX_DIM))**2)
        if (dist2 < dist2_min) then
          ipos = ip
          dist2_min = dist2
        end if

      end do

      write(message(1), '(3a,f5.2,3a)') &
        "Info: spec_full_delta species ", trim(species_label(s)), &
        " atom displaced ", units_from_atomic(units_out%length, sqrt(dist2_min)), &
        " [ ", trim(units_abbrev(units_out%length)), " ]"
      call write_info(1)

      have_point = .true.
#ifdef HAVE_MPI
      ! in parallel we have to find the minimum of the whole grid
      if(gr%mesh%parallel_in_domains) then

        local_min = (/ dist2_min, dble(gr%mesh%mpi_grp%rank)/)
        call MPI_Allreduce(local_min, global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, gr%mesh%mpi_grp%comm, mpi_err)

        if(gr%mesh%mpi_grp%rank /= nint(global_min(2))) have_point = .false.

      end if
#endif
      if(have_point) then
        if(gr%mesh%use_curvilinear) then
          rho(ipos) = -species_z(s)/gr%mesh%vol_pp(ipos)
        else
          rho(ipos) = -species_z(s)/gr%mesh%vol_pp(1)
        end if
      end if

    case(SPEC_FULL_GAUSSIAN)

      ! periodic copies are not considered in this routine
      if(simul_box_is_periodic(gr%mesh%sb)) then
        call messages_devel_version("spec_full_gaussian for periodic systems")
      end if

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      dim = gr%mesh%sb%dim

      SAFE_ALLOCATE(rho_p(1:gr%mesh%np))
      SAFE_ALLOCATE(grho_p(1:gr%mesh%np, 1:dim+1))

      m_p   => gr%mesh
      pos_p = pos

      ! Initial guess.
      call curvilinear_x2chi(gr%mesh%sb, gr%cv, pos, chi0)
      delta   = gr%mesh%h(1)
      alpha   = sqrt(M_TWO)*species_sigma(s)*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:dim) = CNST(1.0)

      ! the dim+1 variable is the normalization of the delta function
      startval(dim+1) = beta

      ! get a better estimate for beta
      call getrho(startval)
      beta = M_ONE / dmf_integrate(gr%mesh, rho_p)
      startval(dim+1) = beta

      ! solve equation
      call root_solver_init(rs, dim+1, &
        solver_type=ROOT_NEWTON, maxiter=500, abs_tolerance=CNST(1.0e-10))
      call droot_solver_run(rs, func, x, conv, startval=startval)

      if(.not.conv) then
        write(message(1),'(a)') 'Internal error in get_species_density.'
        call write_fatal(1)
      end if

      ! we want a charge of -Z
      rho = -species_z(s)*rho_p

      nullify(m_p)
      SAFE_DEALLOCATE_A(grho_p)
      SAFE_DEALLOCATE_A(rho_p)


    case(SPEC_CHARGE_DENSITY)

      call periodic_copy_init(pp, gr%sb, spread(M_ZERO, dim=1, ncopies = MAX_DIM), &
        range = M_TWO * maxval(gr%sb%lsize(1:gr%sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy = periodic_copy_position(pp, gr%sb, icell)
        do i = 1, gr%mesh%np
          call mesh_r(gr%mesh, i, r, x = xx, a = pos)
          xx(1:gr%sb%dim) = xx(1:gr%sb%dim) + yy(1:gr%sb%dim)
          r = sqrt(dot_product(xx(1:gr%sb%dim), xx(1:gr%sb%dim)))
          call parse_expression(rerho, imrho, gr%sb%dim, xx, r, M_ZERO, trim(species_rho_string(s)))
          rho(i) = rho(i) - rerho
        end do
      end do

      call periodic_copy_end(pp)

    end select

    call pop_sub()
  end subroutine species_get_density

  subroutine func(xin, f, jacobian)
    FLOAT, intent(in)  :: xin(:)
    FLOAT, intent(out) :: f(:), jacobian(:,:)

    FLOAT, allocatable :: xrho(:)
    integer :: i, j, dim

    call push_sub('species_pot.func')

    dim = m_p%sb%dim

    call getrho(xin)
    SAFE_ALLOCATE(xrho(1:m_p%np))

    ! First, we calculate the function f.
    do i = 1, dim
      xrho(1:m_p%np) = rho_p(1:m_p%np) * m_p%x(1:m_p%np, i)
      f(i) = dmf_integrate(m_p, xrho) - pos_p(i)
    end do
    f(dim+1) = dmf_integrate(m_p, rho_p) - M_ONE

    ! Now the jacobian.
    do i = 1, dim
      do j = 1, dim+1
        xrho(1:m_p%np) = grho_p(1:m_p%np, j) * m_p%x(1:m_p%np, i)
        jacobian(i, j) = dmf_integrate(m_p, xrho)
      end do
    end do
    do j = 1, dim+1
      xrho(1:m_p%np) = grho_p(1:m_p%np, j)
      jacobian(dim+1, j) = dmf_integrate(m_p, xrho)
    end do

    SAFE_DEALLOCATE_A(xrho)
    call pop_sub()
  end subroutine func


  ! ---------------------------------------------------------
  subroutine getrho(xin)
    FLOAT, intent(in) :: xin(:)

    integer :: i, j, dim
    FLOAT   :: r, chi(MAX_DIM)

    call push_sub('species_pot.getrho')

    dim = m_p%sb%dim
    rho_p = M_ZERO
    do i = 1, m_p%np

      j = i
      if(m_p%parallel_in_domains) &
        j  = m_p%vp%local(m_p%vp%xlocal(m_p%vp%partno)+i-1)

      chi(1:dim) = m_p%idx%Lxyz(j, 1:dim) * m_p%h(1:dim) + m_p%sb%box_offset(1:dim) 

      r = sqrt( sum( (chi(1:dim) - xin(1:dim))**2 ) )

      if( (r/alpha_p)**2 < CNST(10.0)) then
        grho_p(i, dim+1) = exp(-(r/alpha_p)**2)
        rho_p(i)         = xin(dim+1) * grho_p(i, dim+1)
      else
        grho_p(i, dim+1) = M_ZERO
        rho_p(i)         = M_ZERO
      end if

      do j = 1, dim
        grho_p(i, j) = (M_TWO/alpha_p**2) * (chi(j)-xin(j)) * rho_p(i)
      end do
    end do

    call pop_sub()
  end subroutine getrho 

  ! ---------------------------------------------------------
  subroutine species_get_local(s, mesh, x_atom, vl, time)
    type(species_t),  intent(in) :: s
    type(mesh_t),    intent(in) :: mesh
    FLOAT,           intent(in) :: x_atom(MAX_DIM)
    FLOAT,           intent(out):: vl(:)
    FLOAT, optional, intent(in) :: time

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, time_
    integer :: ip, err, idim
    type(ps_t), pointer :: ps

    type(profile_t), save :: prof

    call push_sub('species_pot.species_get_local')
    call profiling_in(prof, "SPECIES_GET_LOCAL")

    time_ = M_ZERO

    if (present(time)) time_ = time

      select case(species_type(s))
      case(SPEC_USDEF)

        do ip = 1, mesh%np
          
          xx = M_ZERO
          xx(1:mesh%sb%dim) = mesh%x(ip,1:mesh%sb%dim)-x_atom(1:mesh%sb%dim)
          r = sqrt(sum(xx(:)**2))
          
          ! Note that as the s%user_def is in input units, we have to convert
          ! the units back and forth
          forall(idim = 1:mesh%sb%dim) xx(idim) = units_from_atomic(units_inp%length, xx(idim))
          r = units_from_atomic(units_inp%length, r)
          vl(ip) = units_to_atomic(units_inp%energy, species_userdef_pot(s, mesh%sb%dim, xx, r, time_))

        end do


      case(SPEC_FROM_FILE)

        call dinput_function(trim(species_filename(s)), mesh, vl, err)
        if(err .ne. 0) then
          write(message(1), '(a)')    'File '//trim(species_filename(s))//'not found.'
          write(message(2), '(a,i4)') 'Error code returned = ', err
          call write_fatal(1)
        end if

      case(SPEC_POINT, SPEC_JELLI)
        a1 = species_z(s)/(M_TWO*species_jradius(s)**3)
        a2 = species_z(s)/species_jradius(s)
        Rb2= species_jradius(s)**2
        
        do ip = 1, mesh%np
          
          xx(:) = mesh%x(ip,:)-x_atom(:)
          r = sqrt(sum(xx(:)**2))
          
          if(r <= species_jradius(s)) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = -species_z(s)/r
          end if
          
        end do
        
      case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
        ps => species_ps(s)
        do ip = 1, mesh%np
          vl(ip) = sum((mesh%x(ip, 1:MAX_DIM) - x_atom(1:MAX_DIM))**2)
        end do
        call spline_eval_vec(ps%vlr_sq, mesh%np, vl)
        nullify(ps)
        
      case(SPEC_FULL_DELTA, SPEC_FULL_GAUSSIAN, SPEC_CHARGE_DENSITY)
        vl(1:mesh%np) = M_ZERO
        
      end select

      call profiling_out(prof)
      call pop_sub()
  end subroutine species_get_local

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Places, in the function phi (defined in each point of the mesh), the
  ! j-th atomic orbital. The orbitals are obtained from the species data
  ! type, and are numbered from one to species_niwfs(spec). It may happen
  ! that there are different orbitals for each spin-polarization direction,
  ! and therefore the orbital is also characterized by the label "is".
  !
  ! In order to put the orbital in the mesh, it is necessary to know where
  ! the species is, and this is given by the vector "pos".
  !
  ! WARNING/TODO: Most of this work should be done inside the species
  ! module, and we should get rid of species_iwf_i, species_ifw_l, etc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine species_get_orbital(spec, mesh, j, dim, is, pos, phi)
    type(species_t),   intent(in)     :: spec
    type(mesh_t),      intent(in)     :: mesh
    integer,           intent(in)     :: j
    integer,           intent(in)     :: dim
    integer,           intent(in)     :: is
    FLOAT,             intent(in)     :: pos(:)
    FLOAT,             intent(out)    :: phi(:)

    integer :: i, l, m, ip
    FLOAT :: r2, x(1:MAX_DIM)
    FLOAT, allocatable :: xf(:, :), ylm(:)
    type(ps_t), pointer :: ps

    call push_sub('species_pot.species_get_orbital')

    call species_iwf_ilm(spec, j, is, i, l, m)

    if(species_is_ps(spec)) then

      ps => species_ps(spec)
      SAFE_ALLOCATE(xf(1:mesh%np, 1:MAX_DIM))
      SAFE_ALLOCATE(ylm(1:mesh%np))

      do ip = 1, mesh%np
        x(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM) - pos(1:MAX_DIM)
        phi(ip) = sum(x(1:MAX_DIM)**2)
        xf(ip, 1:MAX_DIM) = x(1:MAX_DIM)
      end do

      call spline_eval_vec(ps%ur_sq(i, is), mesh%np, phi)
      call loct_ylm(mesh%np, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))

      do ip = 1, mesh%np
        phi(ip) = phi(ip)*ylm(ip)
      end do

      SAFE_DEALLOCATE_A(xf)
      SAFE_DEALLOCATE_A(ylm)

      nullify(ps)
    else

      do ip = 1, mesh%np
        x(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM) - pos(1:MAX_DIM)
        r2 = sum(x(1:MAX_DIM)**2)
        select case(dim)
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

    call pop_sub()
  end subroutine species_get_orbital

  subroutine species_get_orbital_submesh(spec, submesh, j, dim, is, pos, phi)
    type(species_t),   intent(in)  :: spec
    type(submesh_t),   intent(in)  :: submesh
    integer,           intent(in)  :: j
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: is
    FLOAT,             intent(in)  :: pos(:)
    FLOAT,             intent(out) :: phi(:)

    integer :: i, l, m, ip
    FLOAT :: r2, x(1:MAX_DIM), sqrtw, ww
    FLOAT, allocatable :: xf(:, :), ylm(:)
    type(ps_t), pointer :: ps

    if(submesh%ns == 0) return

    call push_sub('species_pot.species_get_orbital_submesh')

    call species_iwf_ilm(spec, j, is, i, l, m)

    if(species_is_ps(spec)) then

      ps => species_ps(spec)
      SAFE_ALLOCATE(xf(1:submesh%ns, 1:MAX_DIM))
      SAFE_ALLOCATE(ylm(1:submesh%ns))
      do ip = 1, submesh%ns
        x(1:MAX_DIM) = submesh%mesh%x(submesh%jxyz(ip), 1:MAX_DIM) - pos(1:MAX_DIM)
        phi(ip) = sum(x(1:MAX_DIM)**2)
        xf(ip, 1:MAX_DIM) = x(1:MAX_DIM)
      end do
      call spline_eval_vec(ps%ur_sq(i, is), submesh%ns, phi)
      call loct_ylm(submesh%ns, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))

      do ip = 1, submesh%ns
        phi(ip) = phi(ip)*ylm(ip)
      end do

      SAFE_DEALLOCATE_A(xf)
      SAFE_DEALLOCATE_A(ylm)

      nullify(ps)
    else
      ww = species_omega(spec)
      sqrtw = sqrt(ww)

      select case(dim)
      case(1)
        do ip = 1, submesh%ns
          x(1:MAX_DIM) = submesh%mesh%x(submesh%jxyz(ip), 1:MAX_DIM) - pos(1:MAX_DIM)
          r2 = sum(x(1:MAX_DIM)**2)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)
        end do
      case(2)
        do ip = 1, submesh%ns
          x(1:MAX_DIM) = submesh%mesh%x(submesh%jxyz(ip), 1:MAX_DIM) - pos(1:MAX_DIM)
          r2 = sum(x(1:MAX_DIM)**2)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)*hermite(l - 1, x(2)*sqrtw)
        end do
      case(3)
        do ip = 1, submesh%ns
          x(1:MAX_DIM) = submesh%mesh%x(submesh%jxyz(ip), 1:MAX_DIM) - pos(1:MAX_DIM)
          phi(ip) = exp(-ww*r2/M_TWO)*hermite(i - 1, x(1)*sqrtw)*hermite(l - 1, x(2)*sqrtw)*hermite(m - 1, x(3)*sqrtw)
        end do
      end select
      
    end if

    call pop_sub()
  end subroutine species_get_orbital_submesh

end module species_pot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
