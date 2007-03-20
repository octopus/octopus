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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: grid.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module specie_pot_m
  use curvlinear_m
  use datasets_m
  use double_grid_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use lib_oct_gsl_spline_m
  use lib_oct_parser_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use root_solver_m
  use simul_box_m
  use specie_m
  use poisson_m
  use units_m
  use varinfo_m

  implicit none

  private
  public ::                   &
    guess_density,            &
    specie_pot_init,          &
    specie_pot_end,           &
    specie_get_density,       & 
    specie_get_gdensity,      &
    specie_get_local,         &
    specie_get_glocal,        &
    specie_get_g2local


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
  subroutine specie_pot_init(this, gr)
    type(specie_t),      intent(inout) :: this
    type(grid_t),        intent(in)    :: gr

    type(loct_spline_t) :: vlc
    integer :: l, k
    FLOAT :: qmax, hmax
    
    call push_sub('specie_pot.specie_pot_init')
    
    if(dg_add_localization_density(gr%dgrid) .and. specie_is_ps(this) ) then

      call loct_spline_init(vlc)
      call dg_get_potential_correction(gr%dgrid, vlc)
      call loct_spline_times(this%z_val, vlc)

      call loct_spline_init(this%ps%vll)
      call loct_spline_sum(vlc, this%ps%vl, this%ps%vll)

      call loct_spline_init(this%ps%dvll)
      call loct_spline_der(this%ps%vll, this%ps%dvll)

      call loct_spline_end(vlc)

      if ( dg_filter_potential(gr%dgrid) ) then 
        
        hmax = dg_get_hmax(gr%dgrid, gr%m)
        qmax = M_PI/hmax
        
        write(message(1),'(a, f13.10, a, f13.10)')  'Info: hmax', hmax, ' qmax ', qmax
        call write_info(1)
        
        call loct_spline_filter(this%ps%vll, fs = (/ qmax, CNST(18.0) /), rs = (/ CNST(3.0), CNST(10.0) /))
        
        do l = 0, this%ps%l_max
          do k = 1, this%ps%kbc
            call loct_spline_filter(this%ps%kb(l, k), l, fs = (/ qmax, CNST(18.0) /), rs = (/ CNST(3.0), CNST(10.0) /))
          end do
        end do
        
      end if
      
    end if

    call pop_sub()

  end subroutine specie_pot_init

  ! ---------------------------------------------------------

  subroutine specie_pot_end(this, gr)
    type(specie_t),      intent(inout) :: this
    type(grid_t),        intent(in)    :: gr

    if(dg_add_localization_density(gr%dgrid) .and. specie_is_ps(this) ) then
      
      call loct_spline_end(this%ps%vll)
      call loct_spline_end(this%ps%dvll)
      
    end if

  end subroutine specie_pot_end

  ! ---------------------------------------------------------
  subroutine atom_density(m, sb, atom, spin_channels, rho)
    type(mesh_t),      intent(in)    :: m
    type(simul_box_t), intent(in)    :: sb
    type(atom_t),      intent(in)    :: atom
    integer,           intent(in)    :: spin_channels
    FLOAT,             intent(inout) :: rho(:, :) ! (m%np, spin_channels)

    integer :: i, in_points, k, n
    FLOAT :: r, x
    FLOAT :: psi1, psi2
    type(specie_t), pointer :: s
#if defined(HAVE_MPI)
    integer :: in_points_red
#endif

    call push_sub('specie_pot.atom_density')

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

      do k = 1, 3**sb%periodic_dim
        do i = 1, m%np
          call mesh_r(m, i, r, a=atom%x + sb%shift(k,:))
          r = max(r, r_small)
          do n = 1, s%ps%conf%p
            select case(spin_channels)
            case(1)
              psi1 = loct_splint(s%ps%Ur(n, 1), r)
              rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
            case(2)
              psi1 = loct_splint(s%ps%Ur(n, 1), r)
              psi2 = loct_splint(s%ps%Ur(n, 2), r)
              rho(i, 1) = rho(i, 1) + s%ps%conf%occ(n, 1)*psi1*psi1 /(M_FOUR*M_PI)
              rho(i, 2) = rho(i, 2) + s%ps%conf%occ(n, 2)*psi2*psi2 /(M_FOUR*M_PI)
            end select
          end do
        end do
      end do
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
    C_POINTER :: blk
    FLOAT :: r, rnd, phi, theta, mag(MAX_DIM)
    FLOAT, allocatable :: atom_rho(:,:)

    call push_sub('specie_pot.guess_density')

    if (spin_channels == 1) then
      gmd_opt = 1
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
      call loct_parse_int(check_inp('GuessMagnetDensity'), INITRHO_FERROMAGNETIC, gmd_opt)
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
      !% specify the coordinates of a vector that has the desired direction. The norm 
      !% of the vector is ignored. Note that it is necessaty to maintain the 
      !% ordering in which the species were defined in the coordinates specifications.
      !%
      !% For spin-polarized calculations the vectors should have only one component and
      !% for non-collinear spin calculations they should have three components.
      !%End
      if(loct_parse_block(check_inp('AtomsMagnetDirection'), blk) < 0) then
        message(1) = "AtomsMagnetDirection block is not defined "
        call write_fatal(1)
      end if

      if (loct_parse_block_n(blk) /= geo%natoms) then
        message(1) = "AtomsMagnetDirection block has the wrong number of rows"
        call write_fatal(1)
      end if

      ALLOCATE(atom_rho(m%np, 2), m%np*2)
      do ia = 1, geo%natoms
        call atom_density(m, sb, geo%atom(ia), 2, atom_rho)

        if (nspin == 2) then
          call loct_parse_block_float(blk, ia-1, 0, mag(1))
          if (mag(1) > M_ZERO) then
            rho(1:m%np, 1:2) = rho(1:m%np, 1:2) + atom_rho(1:m%np, 1:2)
          else
            rho(1:m%np, 1) = rho(1:m%np, 1) + atom_rho(1:m%np, 2)
            rho(1:m%np, 2) = rho(1:m%np, 2) + atom_rho(1:m%np, 1)
          end if

        elseif (nspin == 4) then
          do i = 1, 3
            call loct_parse_block_float(blk, ia-1, i-1, mag(i))
            if (abs(mag(i)) < CNST(1.0e-20)) mag(i) = M_ZERO
          end do

          theta = acos(mag(3)/sqrt(dot_product(mag, mag)))
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
              phi = M_TWO*M_PI - acos(mag(1)/sin(theta)/sqrt(dot_product(mag, mag)))
            elseif (mag(2) >= M_ZERO) then
              phi = acos(mag(1)/sin(theta)/sqrt(dot_product(mag, mag)))
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


  subroutine specie_get_density(s, pos, gr, geo, rho)
    type(specie_t),             intent(in)  :: s
    FLOAT,                      intent(in)  :: pos(MAX_DIM)
    type(grid_t),       target, intent(in)  :: gr
    type(geometry_t),   target, intent(in)  :: geo
    FLOAT,                      intent(out) :: rho(:)

    type(root_solver_t) :: rs
    logical :: conv
    integer :: dim
    FLOAT   :: x(1:MAX_DIM+1), chi0(MAX_DIM), startval(MAX_DIM + 1)
    FLOAT   :: delta, alpha, beta
    FLOAT   :: r, correction
    integer :: ip
    type(loct_spline_t) :: rho_corr

    call push_sub('specie_grid.specie_get_density')

    select case(s%type)

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)

      call loct_spline_init(rho_corr)
      call dg_get_density_correction(gr%dgrid, rho_corr)

      do ip = 1, gr%m%np
        r = sqrt(sum((pos(1:MAX_DIM)-gr%m%x(ip, 1:MAX_DIM))**2))
        rho(ip) = -s%z_val*loct_splint(rho_corr, r)
      end do

      call loct_spline_end(rho_corr)

      correction = -s%z_val / dmf_integrate(gr%m, rho)
      rho(1:gr%m%np) = correction * rho(1:gr%m%np)
      
      write(message(1),'(a, f13.10)')  'Info: Localization charge correction ', correction
      call write_info(1)

    case(SPEC_ALL_E)

      ! --------------------------------------------------------------
      ! Constructs density for an all-electron atom with the procedure
      ! sketched in Modine et al. [Phys. Rev. B 55, 10289 (1997)],
      ! section II.B
      ! --------------------------------------------------------------
      dim = gr%m%sb%dim

      ALLOCATE(rho_p(gr%m%np), gr%m%np)
      ALLOCATE(grho_p(gr%m%np, dim+1), 4*gr%m%np)

      m_p   => gr%m
      pos_p = pos

      ! Initial guess.
      call curvlinear_x2chi(gr%m%sb, geo, gr%cv, pos, chi0)
      delta   = gr%m%h(1)
      alpha   = sqrt(M_TWO)*s%sigma*delta
      alpha_p = alpha  ! global copy of alpha
      beta    = M_ONE

      ! the first dim variables are the position of the delta function
      startval(1:dim) = chi0(1:dim)
      
      ! the dim+1 variable is the normalization of the delta function
      startval(dim+1) = beta

      ! get a better estimate for beta
      call getrho(startval)
      beta = M_ONE / dmf_integrate(gr%m, rho_p)
      startval(dim+1) = beta

      ! solve equation
      call root_solver_init(rs, solver_type=ROOT_NEWTON, maxiter=500, abs_tolerance=CNST(1.0e-10))
      call droot_solver_run(rs, func, x, conv, startval=startval)

      if(.not.conv) then
        write(message(1),'(a)') 'Internal error in get_specie_density.'
        call write_fatal(1)
      end if

      ! we want a charge of -Z
      rho = - s%z * rho_p

      nullify(m_p)
      deallocate(grho_p, rho_p)
    end select

    call pop_sub()
  end subroutine specie_get_density

  subroutine specie_get_gdensity(s, pos, gr, rho)
    type(specie_t),             intent(in)  :: s
    FLOAT,                      intent(in)  :: pos(MAX_DIM)
    type(grid_t),       target, intent(in)  :: gr
    FLOAT,                      intent(out) :: rho(:, :)

    FLOAT   :: x(1:MAX_DIM), r
    integer :: ip
    type(loct_spline_t) :: rho_corr, drho_corr

    call push_sub('specie_grid.specie_get_density')

    select case(s%type)

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)

      call loct_spline_init(rho_corr)
      call loct_spline_init(drho_corr)
      call dg_get_density_correction(gr%dgrid, rho_corr)
      call loct_spline_der(rho_corr, drho_corr)

      do ip = 1, gr%m%np
        x(1:MAX_DIM) = pos(1:MAX_DIM) - gr%m%x(ip, 1:MAX_DIM)
        r = sqrt(sum(x(1:MAX_DIM)**2))
        if ( r > CNST(1e-10) ) then 
          rho(ip, 1:MAX_DIM) = -s%z_val*loct_splint(drho_corr, r)*x(1:MAX_DIM)/r
        else
          rho(ip, 1:MAX_DIM) = M_ZERO
        end if
      end do

    end select

    call pop_sub()
  end subroutine specie_get_gdensity


  ! ---------------------------------------------------------
  subroutine func(xin, f, jacobian)
    FLOAT :: xin(:), f(:), jacobian(:,:)

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
    FLOAT   :: r, chi(MAX_DIM), x(MAX_DIM)

    dim = m_p%sb%dim
    rho_p = M_ZERO; x = M_ZERO
    do i = 1, m_p%np

      j = i
      if(m_p%parallel_in_domains) &
        j  = m_p%vp%local(m_p%vp%xlocal(m_p%vp%partno)+i-1)

      chi(1:dim) = m_p%Lxyz(j, 1:dim) * m_p%h(1:dim) + m_p%sb%box_offset(1:dim) 

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
  subroutine specie_get_local(s, gr, x_atom, vl, time)
    type(specie_t),  intent(in) :: s
    type(grid_t),    intent(in) :: gr
    FLOAT,           intent(in) :: x_atom(MAX_DIM)
    FLOAT,           intent(out):: vl(:)
    FLOAT, optional, intent(in) :: time

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, pot_re, pot_im, time_
    integer :: ip


    time_ = M_ZERO

    if (present(time)) time_ = time

      select case(s%type)
      case(SPEC_USDEF)

        do ip = 1, gr%m%np
          
          xx(:) = gr%m%x(ip,:)-x_atom(:)
          r = sqrt(sum(xx(:)**2))
          
          ! Note that as the s%user_def is in input units, we have to convert
          ! the units back and forth
          xx(:) = xx(:)/units_inp%length%factor ! convert from a.u. to input units
          r = r/units_inp%length%factor
          
          call loct_parse_expression(                            &
               pot_re, pot_im, xx(1), xx(2), xx(3), r, time_, s%user_def)
          vl(ip) = pot_re * units_inp%energy%factor  ! convert from input units to a.u.

        end do

      case(SPEC_POINT, SPEC_JELLI)
        a1 = s%Z/(M_TWO*s%jradius**3)
        a2 = s%Z/s%jradius
        Rb2= s%jradius**2
        
        do ip = 1, gr%m%np
          
          xx(:) = gr%m%x(ip,:)-x_atom(:)
          r = sqrt(sum(xx(:)**2))
          
          if(r <= s%jradius) then
            vl(ip) = (a1*(r*r - Rb2) - a2)
          else
            vl(ip) = - s%Z/r
          end if
          
        end do
        
      case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
        call double_grid_apply_local(gr%dgrid, s, gr%m, x_atom, vl(:))
        
      case(SPEC_ALL_E)
        vl(1:gr%m%np) = M_ZERO
        
      end select
      
  end subroutine specie_get_local

  ! ---------------------------------------------------------
  ! returns the gradient of the external potential
  ! ---------------------------------------------------------
  subroutine specie_get_glocal(s, gr, x_atom, gv, time)
    type(specie_t),  intent(in) :: s
    type(grid_t),    intent(inout) :: gr
    FLOAT,           intent(in) :: x_atom(MAX_DIM)
    FLOAT,          intent(out) :: gv(:, :)
    FLOAT, optional, intent(in) :: time

    FLOAT, parameter :: Delta = CNST(1e-4)
    FLOAT :: x(MAX_DIM), r, l1, l2, pot_re, pot_im, time_
    FLOAT, allocatable :: grho(:, :), gpot(:)
    integer :: i, ip

    gv    = M_ZERO

    time_ = M_ZERO
    if (present(time)) time_ = time

    select case(s%type)
    case(SPEC_USDEF)

      do ip = 1, gr%m%np
        x(:) = gr%m%x(ip, :) - x_atom(:)
        r = sqrt(sum(x(:)**2))
        ! Note that as the s%user_def is in input units, we have to convert
        ! the units back and forth
        x(:) = x(:)/units_inp%length%factor      ! convert from a.u. to input units
        r = r / units_inp%length%factor 
        do i = 1, 3
          x(i) = x(i) - Delta/units_inp%length%factor
          call loct_parse_expression(pot_re, pot_im,           &
               x(1), x(2), x(3), r, time_, s%user_def)
          l1 = pot_re * units_inp%energy%factor     ! convert from input units to a.u.
          
          x(i) = x(i) + M_TWO*Delta/units_inp%length%factor
          call loct_parse_expression(pot_re, pot_im,           &
               x(1), x(2), x(3), r, time_, s%user_def)
          l2 = pot_re * units_inp%energy%factor     ! convert from input units to a.u.
          
          gv(ip, i) = (l2 - l1)/(M_TWO*Delta)
        end do
      end do

    case(SPEC_POINT, SPEC_JELLI)
      l1 = s%Z/(M_TWO*s%jradius**3)

      do ip = 1, gr%m%np
        x(:) = gr%m%x(ip, :) - x_atom(:)
        r = sqrt(sum(x(:)**2))
        if(r <= s%jradius) then
          gv(ip, 1:3) = l1*x(1:3)
        else
          gv(ip, 1:3) = s%Z*x(1:3)/r**3
        end if
      end do

    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)
#if 0
      do ip = 1, gr%m%np

        x(:) = gr%m%x(ip, :) - x_atom(:)
        r = sqrt(sum(x(:)**2))
        if(r>CNST(1e-5)) then
          if (dg_add_localization_density(gr%dgrid)) then 
            dvl_r = loct_splint(s%ps%dvll, r)
          else 
            dvl_r = loct_splint(s%ps%dvl, r)
          end if
          gv(ip, 1:3) = -dvl_r*x(1:3)/r
        else 
          gv(ip, 1:3) = M_ZERO
        end if

      end do

#else 
      call double_grid_apply_glocal(gr%dgrid, s, gr%m, x_atom, gv(:, :))
#endif 

      if (dg_add_localization_density(gr%dgrid)) then
        
        ALLOCATE(grho(NP, MAX_DIM), NP*MAX_DIM)
        ALLOCATE(gpot(NP_PART), NP_PART)
        
        call specie_get_gdensity(s, x_atom, gr, grho)

        do i = 1, NDIM
          call dpoisson_solve(gr, gpot(:), grho(:, i))
          gv(1:gr%m%np, i) = gv(1:gr%m%np, i) + gpot(1:gr%m%np)
        end do

        deallocate(grho)
        deallocate(gpot)

      end if

    case(SPEC_ALL_E)
      gv(1:gr%m%np, 1:3) = M_ZERO
        
    end select

  end subroutine specie_get_glocal

  subroutine specie_get_g2local(s, gr, x_atom, g2v, time)
    type(specie_t),  intent(in) :: s
    type(grid_t),    intent(in) :: gr
    FLOAT,           intent(in) :: x_atom(MAX_DIM)
    FLOAT,          intent(out) :: g2v(:, :, :)
    FLOAT, optional, intent(in) :: time

    FLOAT, parameter :: Delta = CNST(1e-4)
    FLOAT :: x(MAX_DIM), r, dvl_rr, d2vl_r
    integer :: ii, jj, ip


    select case(s%type)
    case(SPEC_PS_PSF, SPEC_PS_HGH, SPEC_PS_CPI, SPEC_PS_FHI, SPEC_PS_UPF)

      do ip = 1, gr%m%np
        x(:) = gr%m%x(ip, :) - x_atom(:)
        r = sqrt(sum(x(:)**2))

        if(r>CNST(0.00001)) then
          dvl_rr = loct_splint(s%ps%dvl,r)/r
          d2vl_r = loct_splint(s%ps%d2vl, r)
          
          do ii= 1, 3
            do jj = 1, 3
              g2v(ip, ii, jj) = ddelta(ii, jj)*dvl_rr  + x(ii)*x(jj)/(r*r)*(d2vl_r - dvl_rr)
            end do
          end do
          
        end if
      end do

    case default
      write(message(1),'(a)')    'Second derivative of the potential not implemented.'
      call write_fatal(1)
      
    end select
    
  end subroutine specie_get_g2local

end module specie_pot_m
