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
  use curvlinear_m
  use datasets_m
  use double_grid_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use loct_parser_m
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
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                    &
    guess_density,             &
    species_pot_init,          &
    species_pot_end,           &
    species_get_density,       &
    species_get_orbital,      &
    species_get_local

  integer, parameter :: INITRHO_PARAMAGNETIC  = 1, &
                        INITRHO_FERROMAGNETIC = 2, &
                        INITRHO_RANDOM        = 3, &
                        INITRHO_USERDEF       = 123

  type(mesh_t),       pointer :: m_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT :: pos_p(MAX_DIM)

contains

 ! ---------------------------------------------------------
  subroutine species_pot_init(this, gr, filter)
    type(species_t),      intent(inout) :: this
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: filter

    character(len=256) :: dirname

    call push_sub('species_pot.species_pot_init')
    
    if(species_is_ps(this)) then
      call ps_separate(this%ps)

      call ps_getradius(this%ps)

      if(filter .ne. PS_FILTER_NONE) then 
        call ps_filter(this%ps, filter, mesh_gcutoff(gr%mesh))
        call ps_getradius(this%ps) ! radius may have changed
      end if

      call ps_derivatives(this%ps)
      
      if(in_debug_mode) then
        write(dirname, '(a)') 'debug/geometry'
        call io_mkdir(dirname)
        call species_debug(trim(dirname), this)
      end if

    end if

    call pop_sub()

  end subroutine species_pot_init

  ! ---------------------------------------------------------
  subroutine species_pot_end(this)
    type(species_t),      intent(inout) :: this
    
    ! for the moment there is nothing to destroy
  end subroutine species_pot_end

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
#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(periodic_copy_t) :: pp

    call push_sub('species_pot.atom_density')

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    s => atom%spec
    rho = M_ZERO

    ! build density ...
    select case (s%type)
    case (SPEC_USDEF, SPEC_ALL_E, SPEC_PS_CPI, SPEC_PS_FHI) ! ... from userdef
      do i = 1, spin_channels
        rho(1:m%np, i) = M_ONE
        x = (real(s%z_val, REAL_PRECISION)/real(spin_channels, REAL_PRECISION)) / dmf_integrate(m, rho(:, i))
        rho(1:m%np, i) = x * rho(1:m%np, i)
      end do

    case (SPEC_CHARGE_DENSITY)
      ! We put, for the electron density, the same as the positive density that create the external potential.

      call periodic_copy_init(pp, sb, spread(M_ZERO, dim=1, ncopies = MAX_DIM), &
        range = M_TWO * maxval(sb%lsize(1:sb%dim)))

      rho = M_ZERO
      do icell = 1, periodic_copy_num(pp)
        yy = periodic_copy_position(pp, sb, icell)
        do i = 1, m%np
          call mesh_r(m, i, r, x = xx, a = atom%x)
          xx(1:sb%dim) = xx(1:sb%dim) + yy(1:sb%dim)
          r = sqrt(dot_product(xx(1:sb%dim), xx(1:sb%dim)))
          call loct_parse_expression(rerho, imrho, sb%dim, xx, r, M_ZERO, trim(s%rho))
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
        if(r <= s%jradius) then
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
        do i = 1, m%np
          call mesh_r(m, i, r, a=atom%x)
          if(r <= s%jradius) then
            rho(i, 1:spin_channels) = real(s%z_val, REAL_PRECISION) /   &
              (m%vol_pp(i)*real(in_points*spin_channels, REAL_PRECISION))
          end if
        end do
      end if

    case (SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_UPF) ! ...from pseudopotential
      ! the outer loop sums densities over atoms in neighbour cells

      call periodic_copy_init(pp, sb, atom%x, &
        range = spline_cutoff_radius(s%ps%Ur(1, 1), s%ps%projectors_sphere_threshold))

      do icell = 1, periodic_copy_num(pp)
        pos = periodic_copy_position(pp, sb, icell)
        do i = 1, m%np
          call mesh_r(m, i, r, pos)
          r = max(r, r_small)
          do n = 1, s%ps%conf%p
            select case(spin_channels)
            case(1)
              psi1 = spline_eval(s%ps%Ur(n, 1), r)
              rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
            case(2)
              psi1 = spline_eval(s%ps%Ur(n, 1), r)
              psi2 = spline_eval(s%ps%Ur(n, 2), r)
              rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
              rho(i, 2) = rho(i, 2) + s%ps%conf%occ(n, 2)*psi2*psi2 /(M_FOUR*M_PI)
            end select
          end do
        end do
      end do

      call periodic_copy_end(pp)

    end select

    call pop_sub()
  end subroutine atom_density

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
      !% When performing spin-polarized or non-collinear spin calculations this option sets 
      !% the guess magnetization density.
      !%
      !% For anti-ferromagnetic configurations the <tt>user_defined</tt> option should be used.
      !%
      !% Note that if the <tt>paramagnetic</tt> option is used the final ground-state will also be
      !% paramagnetic, but the same is not true for the other options.
      !%Option paramagnetic 1
      !% Magnetization density is zero.
      !%Option ferromagnetic 2
      !% Magnetization density is the sum of the atomic magnetization densities.
      !%Option random 3
      !% Each atomic magnetization density is randomly rotated.
      !%Option user_defined 123
      !% The atomic magnetization densities are rotated so that the magnetization 
      !% vector has the same direction as a vector provided by the user. In this case,
      !% the <tt>AtomsMagnetDirection</tt> block has to be set.
      !%End
      call loct_parse_int(datasets_check('GuessMagnetDensity'), INITRHO_FERROMAGNETIC, gmd_opt)
      if(.not.varinfo_valid_option('GuessMagnetDensity', gmd_opt)) call input_error('GuessMagnetDensity')
      call messages_print_var_option(stdout, 'GuessMagnetDensity', gmd_opt)
    end if

    rho = M_ZERO
    select case (gmd_opt)
    case (INITRHO_PARAMAGNETIC)
      ALLOCATE(atom_rho(m%np, 1), m%np*1)
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 1, atom_rho(1:m%np, 1:1))
        rho(1:m%np, 1:1) = rho(1:m%np, 1:1) + atom_rho(1:m%np, 1:1)
      end do
      if (spin_channels == 2) then
        rho(1:m%np, 1) = M_HALF*rho(1:m%np, 1)
        rho(1:m%np, 2) = rho(1:m%np, 1)
      end if

    case (INITRHO_FERROMAGNETIC)
      ALLOCATE(atom_rho(m%np, 2), m%np*2)
      atom_rho = M_ZERO
      rho = M_ZERO
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 2, atom_rho(1:m%np, 1:2))
        rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
      end do

    case (INITRHO_RANDOM) ! Random oriented spins
      ALLOCATE(atom_rho(m%np, 2), m%np*2)
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
      !% for non-collinear spin calculations they should have three components.
      !%End
      if(loct_parse_block(datasets_check('AtomsMagnetDirection'), blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined "
        call write_fatal(1)
      end if

      if (loct_parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows"
        call write_fatal(1)
      end if

      ALLOCATE(atom_rho(m%np, 2), m%np*2)
      do ia = 1, geo%natoms
        !Read from AtomsMagnetDirection block 
        if (nspin == 2) then
          call loct_parse_block_float(blk, ia-1, 0, mag(1))
          lmag = abs(mag(1))
        elseif (nspin == 4) then
          do i = 1, 3
            call loct_parse_block_float(blk, ia-1, i-1, mag(i))
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

      call loct_parse_block_end(blk)

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

    deallocate(atom_rho)
    call pop_sub()
  end subroutine guess_density


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
    integer :: icell
    type(periodic_copy_t) :: pp

    call push_sub('species_pot.species_get_density')

    select case(s%type)

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
      rho = M_ZERO
      call periodic_copy_init(pp, gr%sb, pos, range = spline_cutoff_radius(s%ps%nlr, s%ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        call dmf_put_radial_spline(gr%mesh, s%ps%nlr, periodic_copy_position(pp, gr%sb, icell), rho, add = .true.)
      end do
      call periodic_copy_end(pp)
      
    case(SPEC_ALL_E)

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      dim = gr%mesh%sb%dim

      ALLOCATE(rho_p(gr%mesh%np), gr%mesh%np)
      ALLOCATE(grho_p(gr%mesh%np, dim+1), 4*gr%mesh%np)

      m_p   => gr%mesh
      pos_p = pos

      ! Initial guess.
      call curvlinear_x2chi(gr%mesh%sb, geo, gr%cv, pos, chi0)
      delta   = gr%mesh%h(1)
      alpha   = sqrt(M_TWO)*s%sigma*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:dim) = chi0(1:dim)
      
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
      rho = - s%z * rho_p

      nullify(m_p)
      deallocate(grho_p, rho_p)

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
          call loct_parse_expression(rerho, imrho, gr%sb%dim, xx, r, M_ZERO, trim(s%rho))
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

    dim = m_p%sb%dim

    call getrho(xin)
    ALLOCATE(xrho(m_p%np), m_p%np)

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

    deallocate(xrho)
  end subroutine func


  ! ---------------------------------------------------------
  subroutine getrho(xin)
    FLOAT, intent(in) :: xin(:)

    integer :: i, j, dim
    FLOAT   :: r, chi(MAX_DIM)

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

  end subroutine getrho 

  ! ---------------------------------------------------------
  subroutine species_get_local(s, mesh, x_atom, vl, time)
    type(species_t),  intent(in) :: s
    type(mesh_t),    intent(in) :: mesh
    FLOAT,           intent(in) :: x_atom(MAX_DIM)
    FLOAT,           intent(out):: vl(:)
    FLOAT, optional, intent(in) :: time

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, pot_re, pot_im, time_
    integer :: ip

    type(profile_t), save :: prof

    call profiling_in(prof, "SPECIES_GET_LOCAL")

    time_ = M_ZERO

    if (present(time)) time_ = time

      select case(s%type)
      case(SPEC_USDEF)

        do ip = 1, mesh%np
          
          xx(:) = mesh%x(ip,:)-x_atom(:)
          r = sqrt(sum(xx(:)**2))
          
          ! Note that as the s%user_def is in input units, we have to convert
          ! the units back and forth
          xx(:) = xx(:)/units_inp%length%factor ! convert from a.u. to input units
          r = r/units_inp%length%factor
          
          call loct_parse_expression(                            &
               pot_re, pot_im, mesh%sb%dim, xx, r, time_, s%user_def)
          vl(ip) = pot_re * units_inp%energy%factor  ! convert from input units to a.u.

        end do

      case(SPEC_POINT, SPEC_JELLI)
        a1 = s%Z/(M_TWO*s%jradius**3)
        a2 = s%Z/s%jradius
        Rb2= s%jradius**2
        
        do ip = 1, mesh%np
          
          xx(:) = mesh%x(ip,:)-x_atom(:)
          r = sqrt(sum(xx(:)**2))
          
          if(r <= s%jradius) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = - s%Z/r
          end if
          
        end do
        
      case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
        do ip = 1, mesh%np
          vl(ip) = sum((mesh%x(ip, 1:MAX_DIM) - x_atom(1:MAX_DIM))**2)
        end do
        call spline_eval_vec(s%ps%vlr_sq, mesh%np, vl)
        
      case(SPEC_ALL_E, SPEC_CHARGE_DENSITY)
        vl(1:mesh%np) = M_ZERO
        
      end select

      call profiling_out(prof)
      
  end subroutine species_get_local

  subroutine species_get_orbital(s, mesh, j, dim, is, pos, phi)
    type(species_t),   intent(in)  :: s
    type(mesh_t),      intent(in)  :: mesh
    integer,           intent(in)  :: j
    integer,           intent(in)  :: dim
    integer,           intent(in)  :: is
    FLOAT,             intent(in)  :: pos(:)
    FLOAT,             intent(out) :: phi(:)

    integer :: i, l, m, ip
    FLOAT :: r2, x(1:MAX_DIM)
    FLOAT, allocatable :: xf(:, :), ylm(:)

    i = s%iwf_i(j, is)
    l = s%iwf_l(j, is)
    m = s%iwf_m(j, is)

    if(species_is_ps(s)) then

      ALLOCATE(xf(1:mesh%np, 1:MAX_DIM), mesh%np*MAX_DIM)
      ALLOCATE(ylm(1:mesh%np), mesh%np)

      do ip = 1, mesh%np
        x(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM) - pos(1:MAX_DIM)
        phi(ip) = sqrt(sum(x(1:MAX_DIM)**2))
        xf(ip, 1:MAX_DIM) = x(1:MAX_DIM)
      end do

      call spline_eval_vec(s%ps%ur(i, is), mesh%np, phi)
      call loct_ylm(mesh%np, xf(1, 1), xf(1, 2), xf(1, 3), l, m, ylm(1))

      do ip = 1, mesh%np
        phi(ip) = phi(ip)*ylm(ip)
      end do

    else

      do ip = 1, mesh%np
        x(1:MAX_DIM) = mesh%x(ip, 1:MAX_DIM) - pos(1:MAX_DIM)
        r2 = sum(x(1:MAX_DIM)**2)
        select case(dim)
        case(1)
          phi(ip) = exp(-s%omega*r2/M_TWO)*hermite(i - 1, x(1)*sqrt(s%omega))
        case(2)
          phi(ip) = exp(-s%omega*r2/M_TWO)*hermite(i - 1, x(1)*sqrt(s%omega))*hermite(l - 1, x(2)*sqrt(s%omega))
        case(3)
          phi(ip) = exp(-s%omega*r2/M_TWO)*&
               hermite(i - 1, x(1)*sqrt(s%omega))*hermite(l - 1, x(2)*sqrt(s%omega))*hermite(m - 1, x(3)*sqrt(s%omega))
        end select
      end do
    end if

  end subroutine species_get_orbital
  
end module species_pot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
