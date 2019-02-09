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
  use double_grid_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use io_function_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use periodic_copy_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use species_oct_m
  use splines_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                         &
    species_get_density,            &
    species_get_nlcc,               &
    species_get_local,              &
    species_atom_density,           &
    species_atom_density_grad

  type(mesh_t), pointer :: mesh_p
  FLOAT, allocatable :: rho_p(:)
  FLOAT, allocatable :: grho_p(:, :)
  FLOAT :: alpha_p
  FLOAT :: pos_p(MAX_DIM)

contains


  ! ---------------------------------------------------------
  subroutine species_atom_density(mesh, sb, atom, spin_channels, rho)
    type(mesh_t),         intent(in)    :: mesh
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: rho(:, :) !< (mesh%np, spin_channels)

    integer :: isp, ip, in_points, icell
    FLOAT :: rr, x, pos(1:MAX_DIM), nrm, rmax
    FLOAT :: xx(MAX_DIM), yy(MAX_DIM), rerho, imrho
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps

#if defined(HAVE_MPI)
    integer :: in_points_red
#endif
    type(periodic_copy_t) :: pp

    PUSH_SUB(species_atom_density)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    species => atom%species
    rho = M_ZERO

    ! build density ...
    select case (species_type(species))
    case (SPECIES_FULL_DELTA) ! ... from userdef
      do isp = 1, spin_channels
        rho(1:mesh%np, isp) = M_ONE
        x = (species_zval(species)/real(spin_channels, REAL_PRECISION)) / dmf_integrate(mesh, rho(:, isp))
        rho(1:mesh%np, isp) = x * rho(1:mesh%np, isp)
      end do

    case (SPECIES_PSEUDO)
      ! ...from pseudopotentials
      
      pos(1:MAX_DIM) = M_ZERO
      ps => species_ps(species)

      if(ps_has_density(ps)) then

        ASSERT(associated(ps%density))

        rmax = CNST(0.0)
        do isp = 1, spin_channels
          rmax = max(rmax, spline_cutoff_radius(ps%density(isp), ps%projectors_sphere_threshold))
        end do
        
        call periodic_copy_init(pp, sb, atom%x, rmax)

        do icell = 1, periodic_copy_num(pp)
          pos(1:sb%dim) = periodic_copy_position(pp, sb, icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, r_small)
            
            do isp = 1, spin_channels
              if(rr >= spline_range_max(ps%density(isp))) cycle
              rho(ip, isp) = rho(ip, isp) + spline_eval(ps%density(isp), rr)
            end do
            
          end do
        end do
  
        call periodic_copy_end(pp)

      else 

        !we use the square root of the short-range local potential, just to put something that looks like a density

        call periodic_copy_init(pp, sb, atom%x, range = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold))
      
        do icell = 1, periodic_copy_num(pp)
          pos(1:sb%dim) = periodic_copy_position(pp, sb, icell)
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, r_small)

            if(rr >= spline_range_max(ps%vl)) cycle
            
            do isp = 1, spin_channels
              rho(ip, isp) = rho(ip, isp) + sqrt(abs(spline_eval(ps%vl, rr)))
            end do
              
          end do
        end do
  
        call periodic_copy_end(pp)

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
  ! Gradient of the atomic density, if available
  subroutine species_atom_density_grad(mesh, sb, atom, spin_channels, drho)
    type(mesh_t),         intent(in)    :: mesh
    type(simul_box_t),    intent(in)    :: sb
    type(atom_t), target, intent(in)    :: atom
    integer,              intent(in)    :: spin_channels
    FLOAT,                intent(inout) :: drho(:, :, :) !< (mesh%np, spin_channels, dim)

    integer :: isp, ip, icell, idir
    FLOAT :: rr, pos(1:MAX_DIM), range
    type(species_t), pointer :: species
    type(ps_t), pointer :: ps
    type(periodic_copy_t) :: pp

    PUSH_SUB(species_atom_density_grad)

    ASSERT(spin_channels == 1 .or. spin_channels == 2)

    species => atom%species
    drho = M_ZERO

    ! build density ...
    select case (species_type(species))

    case (SPECIES_PSEUDO)
      ! ...from pseudopotentials
      
      pos(1:MAX_DIM) = M_ZERO
      ps => species_ps(species)

      if(ps_has_density(ps)) then

        range = spline_cutoff_radius(ps%density_der(1), ps%projectors_sphere_threshold)
        if (spin_channels == 2) range = max(range, spline_cutoff_radius(ps%density_der(2), ps%projectors_sphere_threshold))
        call periodic_copy_init(pp, sb, atom%x, range = range)

        do icell = 1, periodic_copy_num(pp)
          pos(1:sb%dim) = periodic_copy_position(pp, sb, icell)
        
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, origin = pos)
            rr = max(rr, r_small)

            do idir = 1, mesh%sb%dim
              do isp = 1, spin_channels
                if(rr >= spline_range_max(ps%density_der(isp))) cycle
                drho(ip, isp, idir) = drho(ip, isp, idir) - spline_eval(ps%density_der(isp), rr)*(mesh%x(ip, idir)-pos(idir))/rr
              end do
           end do
          end do
        end do
  
        call periodic_copy_end(pp)

      else 
        call messages_write('The pseudopotential for')
        call messages_write(species_label(species))
        call messages_write(' does not contain the density.')
        call messages_fatal()
      end if
      
    case default
      call messages_not_implemented('species_atom_density_grad for non-pseudopotential species')

    end select

    POP_SUB(species_atom_density_grad)
  end subroutine species_atom_density_grad

  ! ---------------------------------------------------------

  subroutine species_get_density(species, pos, mesh, rho)
    type(species_t),    target, intent(in)  :: species
    FLOAT,                      intent(in)  :: pos(:)
    type(mesh_t),       target, intent(in)  :: mesh
    FLOAT,                      intent(out) :: rho(:)

    logical :: conv
    integer :: dim
    FLOAT   :: x(1:MAX_DIM+1), startval(MAX_DIM + 1)
    FLOAT   :: delta, alpha, beta, xx(MAX_DIM), yy(MAX_DIM), rr, imrho1, rerho
    FLOAT   :: dist2, dist2_min
    integer :: icell, ipos, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps
    logical :: have_point
#ifdef HAVE_MPI
    real(8) :: local_min(2), global_min(2)
#endif
    type(submesh_t)       :: sphere
    type(profile_t), save :: prof
    FLOAT,    allocatable :: rho_sphere(:)
    FLOAT, parameter      :: threshold = CNST(1e-6)
    FLOAT                 :: norm_factor
    
    PUSH_SUB(species_get_density)

    call profiling_in(prof, "SPECIES_DENSITY")

    select case(species_type(species))

    case(SPECIES_PSEUDO)
      ps => species_ps(species)

      call submesh_init(sphere, mesh%sb, mesh, pos, spline_cutoff_radius(ps%nlr, threshold))
      SAFE_ALLOCATE(rho_sphere(1:sphere%np))
      
      forall(ip = 1:sphere%np) rho_sphere(ip) = sphere%x(ip, 0)
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

        dist2 = sum((mesh%x(ip, 1:mesh%sb%dim) - pos(1:mesh%sb%dim))**2)
        if (dist2 < dist2_min) then
          ipos = ip
          dist2_min = dist2
        end if

      end do

      have_point = .true.
#ifdef HAVE_MPI
      ! in parallel we have to find the minimum of the whole grid
      if(mesh%parallel_in_domains) then

        local_min = (/ dist2_min, dble(mesh%mpi_grp%rank)/)
        call MPI_Allreduce(local_min, global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, mesh%mpi_grp%comm, mpi_err)

        if(mesh%mpi_grp%rank /= nint(global_min(2))) have_point = .false.

        dist2_min = global_min(1)
      end if
#endif

      if(have_point) rho(ipos) = -species_z(species)/mesh%volume_element

      write(message(1), '(3a,f5.2,3a)') &
        "Info: species_full_delta species ", trim(species_label(species)), &
        " atom displaced ", units_from_atomic(units_out%length, sqrt(dist2_min)), &
        " [ ", trim(units_abbrev(units_out%length)), " ]"
      call messages_info(1)

    end select

    call profiling_out(prof)
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
  subroutine species_get_nlcc(species, pos, mesh, rho_core)
    type(species_t), target, intent(in)  :: species
    FLOAT,                   intent(in)  :: pos(MAX_DIM)
    type(mesh_t),            intent(in)  :: mesh
    FLOAT,                   intent(out) :: rho_core(:)

    FLOAT :: center(MAX_DIM), rr
    integer :: icell, ip
    type(periodic_copy_t) :: pp
    type(ps_t), pointer :: ps

    PUSH_SUB(species_get_nlcc)

    ! only for 3D pseudopotentials, please
    if(species_is_ps(species)) then
      ps => species_ps(species)
      rho_core = M_ZERO
      call periodic_copy_init(pp, mesh%sb, pos, range = spline_cutoff_radius(ps%core, ps%projectors_sphere_threshold))
      do icell = 1, periodic_copy_num(pp)
        center(1:mesh%sb%dim) = periodic_copy_position(pp, mesh%sb, icell)
        do ip = 1, mesh%np
          rr = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
          if(rr < spline_range_max(ps%core)) then
            rho_core(ip) = rho_core(ip) + spline_eval(ps%core, rr)
          end if
        end do
      end do
      call periodic_copy_end(pp)
    else
      rho_core = M_ZERO
    end if

    POP_SUB(species_get_nlcc)
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
        jp = mesh_p%vp%local(mesh_p%vp%xlocal+ip-1)

      chi(1:dim) = mesh_p%idx%lxyz(jp, 1:dim) * mesh_p%spacing(1:dim)

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
  !> used when the density is not available, or otherwise the Poisson eqn would be used instead
  subroutine species_get_local(species, mesh, x_atom, vl)
    type(species_t), target, intent(in)  :: species
    type(mesh_t),            intent(in)  :: mesh
    FLOAT,                   intent(in)  :: x_atom(:)
    FLOAT,                   intent(out) :: vl(:)

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(MAX_DIM), r, r2
    integer :: ip, err, idim
    type(ps_t), pointer :: ps
    CMPLX :: zpot

    type(profile_t), save :: prof

    PUSH_SUB(species_get_local)

    call profiling_in(prof, "SPECIES_GET_LOCAL")

      select case(species_type(species))
      case(SPECIES_PSEUDO)
       
        ps => species_ps(species)

        do ip = 1, mesh%np
          r2 = sum((mesh%x(ip, 1:mesh%sb%dim) - x_atom(1:mesh%sb%dim))**2)
          if(r2 < spline_range_max(ps%vlr_sq)) then
            vl(ip) = spline_eval(ps%vlr_sq, r2)
          else
            vl(ip) = P_PROTON_CHARGE*species_zval(species)/sqrt(r2)
          end if
        end do

        nullify(ps)
        
      case(SPECIES_FULL_DELTA)
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
