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
!! $Id$

#include "global.h"

module external_pot_m
  use datasets_m
  use derivatives_m
  use double_grid_m
  use functions_m
  use gauge_field_m
  use global_m
  use grid_m
  use io_m
  use ion_interaction_m
  use lalg_basic_m
  use loct_parser_m
  use splines_m
  use magnetic_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use simul_box_m
  use units_m
  use logrid_m
  use poisson_cutoffs_m
  use ps_m
  use specie_m
  use specie_pot_m
  use spline_filter_m
  use solids_m
  use geometry_m
  use states_m
  use submesh_m
  use lasers_m
  use profiling_m
  use mpi_m
  use mpi_debug_m
  use varinfo_m
  use poisson_m
  use projector_m

  implicit none

  private
  public ::                    &
    epot_t,                    &
    epot_init,                 &
    epot_end,                  &
    epot_generate,             &
    epot_forces,               &
    epot_local_potential


  type epot_t
    ! Classic charges:
    integer :: classic_pot        ! How to include the classic charges
    FLOAT, pointer :: vclassic(:) ! We use it to store the potential of the classic charges

    ! Ions
    FLOAT,       pointer :: vpsl(:)       ! the local part of the pseudopotentials
    type(projector_t), pointer :: p(:)    ! non-local projectors
    logical :: non_local
    integer :: natoms

    ! External e-m fields
    integer :: no_lasers                   ! number of laser pulses used
    type(laser_t), pointer :: lasers(:)    ! lasers stuff
    FLOAT, pointer :: E_field(:)           ! static electric field
    FLOAT, pointer :: v_static(:)          ! static scalar potential
    FLOAT, pointer :: B_field(:)           ! static magnetic field
    FLOAT, pointer :: A_static(:,:)        ! static vector potential
    type(gauge_field_t) :: gfield
    integer :: reltype            ! type of relativistic correction to use

    ! The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    ! *effective* electrons in a quantum dot. It affects the spin Zeeman term.
    FLOAT :: gyromagnetic_ratio

    ! SO prefactor (1.0 = normal SO, 0.0 = no SO)
    FLOAT :: so_strength
    
    FLOAT :: eii
    FLOAT, pointer :: fii(:, :)

  end type epot_t

  integer, public, parameter :: &
    NOREL      = 0,             &
    SPIN_ORBIT = 1,             &
    APP_ZORA   = 2,             &
    ZORA       = 3

contains

  ! ---------------------------------------------------------
  subroutine epot_init(ep, gr, geo, ispin)
    type(epot_t),     intent(out)   :: ep
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(inout) :: geo
    integer,          intent(in)    :: ispin

    integer :: i
    type(block_t) :: blk
    FLOAT, allocatable :: x(:)
    integer :: filter

    call push_sub('epot.epot_init')

    !%Variable FilterPotentials
    !%Type integer
    !%Default 1
    !%Section Hamiltonian
    !%Description
    !% Octopus filters the pseudopotentials so that they no
    !% longer contain Fourier components larger than the mesh itself. This is
    !% very useful to decrease the egg-box effect, and so should be used in
    !% all instances where atoms move.
    !%Option filter_none 1
    !% Do not filter
    !%Option filter_TS 2
    !% The filter of M. Tafipolsky and R. Schmid, J. Chem. Phys. 124, 174102 (2006)
    !%Option filter_BSB 3
    !% The filter of E. L. Briggs, D. J. Sullivan, and J. Bernholc, Phys. Rev. B 54, 14362 (1996)
    !%End
    call loct_parse_int(check_inp('FilterPotentials'), PS_FILTER_NONE, filter)
    if(.not.varinfo_valid_option('FilterPotentials', filter)) call input_error('FilterPotentials')
    call messages_print_var_option(stdout, "FilterPotentials", filter)

    if(filter == PS_FILTER_TS) call spline_filter_mask_init()
    do i = 1, geo%nspecies
      call specie_pot_init(geo%specie(i), gr, filter)
    end do

    ! Local part of the pseudopotentials
    ALLOCATE(ep%vpsl(NP), NP)
    !$omp parallel workshare
    ep%vpsl(1:NP) = M_ZERO
    !$omp end parallel workshare

    ep%classic_pot = 0
    if(geo%ncatoms > 0) then

      !%Variable ClassicPotential
      !%Type integer
      !%Default 0
      !%Section Hamiltonian
      !%Description
      !% If <tt>true</tt>, add to the external potential the potential generated by 
      !% the point charges read from the PDB input (see <tt>PBDCoordinates</tt>).
      !%End
      call loct_parse_int(check_inp('ClassicPotential'), 0, ep%classic_pot)
      if(ep%classic_pot > 0) then
        message(1) = 'Info: generating classic external potential'
        call write_info(1)

        ALLOCATE(ep%Vclassic(NP), NP)
        call epot_generate_classic(ep, gr%m, geo)
      end if
    end if

    ! lasers
    call laser_init(ep%no_lasers, ep%lasers, gr%m)

    ! No more "UserDefinedTDPotential" from this version on.
    call obsolete_variable('UserDefinedTDPotential', 'TDExternalFields')

    !%Variable StaticElectricField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% A static constant electrical field may be added to the usual Hamiltonian,
    !% by setting the block StaticElectricField.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the electrical field vector.
    !%End
    nullify(ep%E_field, ep%v_static)
    if(loct_parse_block(check_inp('StaticElectricField'), blk)==0) then
      ALLOCATE(ep%E_field(NDIM), NDIM)
      do i = 1, NDIM
        call loct_parse_block_float(blk, 0, i-1, ep%E_field(i))
      end do
      call loct_parse_block_end(blk)

      ep%E_field(:) = ep%E_field(:) * units_inp%energy%factor/units_inp%length%factor

      ! Compute the scalar potential
      ALLOCATE(ep%v_static(NP), NP)
      do i = 1, NP
        ep%v_static(i) = sum(gr%m%x(i,:)*ep%E_field(:))
      end do
    end if

    !%Variable StaticMagneticField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% A static constant magnetic field may be added to the usual Hamiltonian,
    !% by setting the block StaticMagneticField. 
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the magnetic field vector. Note that
    !% if you are running the code in 1D mode this will not work, and if you
    !% are running the code in 2D mode the magnetic field will have to be in
    !% the z-direction, so that the first two columns should be zero.
    !%
    !% The magnetic field should always be entered in atomic units, regardless
    !% of the "Units" variable. Note that we use the "Gaussian" system
    !% meaning 1 au[B] = 1.7152553 * 10^7 gauss, which corresponds to
    !% 1.7152553 * 10^3 Tesla.
    !%End
    nullify(ep%B_field, ep%A_static)
    if(loct_parse_block(check_inp('StaticMagneticField'), blk)==0) then

      ALLOCATE(ep%B_field(3), 3)
      do i = 1, 3
        call loct_parse_block_float(blk, 0, i-1, ep%B_field(i))
      end do
      select case(calc_dim)
      case(1)
        call input_error('StaticMagneticField')
      case(2)
        if(ep%B_field(1)**2+ep%B_field(2)**2 > M_ZERO) call input_error('StaticMagneticField')
      end select
      call loct_parse_block_end(blk)

      ! Compute the vector potential
      ALLOCATE(ep%A_static(NP, NDIM), NP*NDIM)
      ALLOCATE(x(NDIM), NDIM)
      do i = 1, NP
        x(1:NDIM) = gr%m%x(i, 1:NDIM)
        select case (NDIM)
        case (2)
          ep%A_static(i, :) = (/x(2), -x(1)/)*ep%B_field(3)
        case (3)
          ep%A_static(i, :) = (/x(2)*ep%B_field(3) - x(3)*ep%B_field(2), &
            x(3)*ep%B_field(1) - x(1)*ep%B_field(3), x(1)*ep%B_field(2) - x(2)*ep%B_field(1)/)
        end select
      end do
      deallocate(x)
      ep%A_static = -M_HALF/P_c*ep%A_static

    end if
    
    !%Variable GyromagneticRatio
    !%Type float
    !%Default 2.0023193043768
    !%Section Hamiltonian
    !%Description
    !% The gyromagnetic ratio of the electron. This is of course a physical 
    !% constant, and the default value is the exact one that you should not 
    !% touch, unless : 
    !% 
    !% (i)  You want to disconnect the anomalous Zeeman term in the Hamiltonian 
    !% (then set it to zero, this number only affects this term);
    !% 
    !% (ii) You are using an effective Hamiltonian, as it is the case when
    !% you calculate a 2D electron gas, in which case you have an effective
    !% gyromagnetic factor that depends on the material.
    !%End
    call loct_parse_float(check_inp('GyromagneticRatio'), P_g, ep%gyromagnetic_ratio)

    !%Variable RelativisticCorrection
    !%Type integer
    !%Default non_relativistic
    !%Section Hamiltonian
    !%Description
    !% The default value means that <i>no</i> relativistic correction is used. To
    !% include spin-orbit coupling turn <tt>RelativisticCorrection</tt> to <tt>spin_orbit</tt> 
    !% (this will only work if <tt>SpinComponents</tt> has been set to <tt>non_collinear</tt>, which ensures
    !% the use of spinors).
    !%Option non_relativistic 0
    !% No relativistic corrections.
    !%Option spin_orbit 1
    !% Spin-Orbit.
    !%Option app_zora 2
    !% Approximated ZORA (Not implemented)
    !%Option zora 3
    !% ZORA (Not implemented)
    !%End
    call loct_parse_int(check_inp('RelativisticCorrection'), NOREL, ep%reltype)
    if(.not.varinfo_valid_option('RelativisticCorrection', ep%reltype)) call input_error('RelativisticCorrection')
    if (ispin /= SPINORS .and. ep%reltype == SPIN_ORBIT) then
      message(1) = "The Spin-orbit term can only be applied when using Spinors."
      call write_fatal(1)
    end if
    ! This is temporary...
    if(ep%reltype > SPIN_ORBIT) then
      message(1) = 'Error: ZORA corrections not implemented.'
      call write_fatal(1)
    end if
    call messages_print_var_option(stdout, "RelativisticCorrection", ep%reltype)

    !%Variable SOStrength
    !%Type float
    !%Default M_ONE
    !%Section Hamiltonian
    !%Description
    !% Tuning of the spin-orbit coupling strength: setting this value to zero turns off spin-orbit terms in
    !% the hamiltonian, and setting it to one corresponds to full spin-orbit.
    !%End
    call loct_parse_float(check_inp('SOStrength'), M_ONE, ep%so_strength)
    
    ALLOCATE(ep%p(geo%natoms), geo%natoms)

    do i = 1, geo%natoms
      call projector_null(ep%p(i))
    end do
    ep%natoms = geo%natoms
    ep%non_local = .false.

    ALLOCATE(ep%fii(1:MAX_DIM, 1:geo%natoms), MAX_DIM*geo%natoms)

    call gauge_field_init(ep%gfield, gr%sb)

    call pop_sub()
  end subroutine epot_init


  ! ---------------------------------------------------------
  subroutine epot_end(ep, gr, geo)
    type(epot_t),      intent(inout) :: ep
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(inout) :: geo

    integer :: i, iproj

    call push_sub('epot.epot_end')

    deallocate(ep%fii)

    do i = 1, geo%nspecies
      call specie_pot_end(geo%specie(i), gr)
    end do

    if(associated(ep%vpsl)) then
      deallocate(ep%vpsl)
      nullify(ep%vpsl)
    end if

    if(ep%classic_pot > 0) then
      ep%classic_pot = 0
      ! sanity check
      ASSERT(associated(ep%Vclassic)) 
      deallocate(ep%Vclassic)         ! and clean up
      nullify(ep%Vclassic)
    end if

    ! the external laser
    call laser_end(ep%no_lasers, ep%lasers)

    ! the macroscopic fields
    if(associated(ep%E_field))  deallocate(ep%E_field)
    if(associated(ep%v_static)) deallocate(ep%v_static)
    if(associated(ep%B_field))  deallocate(ep%B_field)
    if(associated(ep%A_static)) deallocate(ep%A_static)

    do iproj = 1, geo%natoms
      if(.not. specie_is_ps(geo%atom(iproj)%spec)) cycle
      call projector_end(ep%p(iproj))
    end do
    
    ASSERT(associated(ep%p))
    deallocate(ep%p)

    call pop_sub()

  end subroutine epot_end

  ! ---------------------------------------------------------

  subroutine epot_generate(ep, gr, geo, st, time)
    type(epot_t),      intent(inout) :: ep
    type(grid_t), target,  intent(inout) :: gr
    type(geometry_t),  intent(inout) :: geo
    type(states_t),    intent(inout) :: st
    FLOAT,   optional, intent(in)    :: time

    FLOAT   :: time_
    integer :: ia
    type(atom_t),     pointer :: atm

    type(mesh_t),      pointer :: m
    type(simul_box_t), pointer :: sb
    type(profile_t), save :: epot_generate_prof

#ifdef HAVE_MPI
    FLOAT,    allocatable :: vpsltmp(:)
    type(profile_t), save :: epot_reduce
#endif

    call profiling_in(epot_generate_prof, "EPOT_GENERATE")
    call push_sub('epot.epot_generate')

    sb  => gr%sb
    m   => gr%m

    time_ = M_ZERO
    if (present(time)) time_ = time

    ! Local.
    ep%vpsl = M_ZERO
    do ia = geo%atoms_start, geo%atoms_end
      call epot_local_potential(ep, gr, geo, geo%atom(ia), ep%vpsl, time_, st%rho_core)
    end do

#ifdef HAVE_MPI
    call profiling_in(epot_reduce, "EPOT_REDUCE")
    if(geo%parallel_in_atoms) then
      ALLOCATE(vpsltmp(1:NP), NP)
      call MPI_Allreduce(ep%vpsl, vpsltmp, NP, MPI_FLOAT, MPI_SUM, geo%mpi_grp%comm, mpi_err)
      call lalg_copy(NP, vpsltmp, ep%vpsl)
      deallocate(vpsltmp)
    end if
    call profiling_out(epot_reduce)
#endif

    ! we assume that we need to recalculate the ion ion energy
    call ion_interaction_calculate(geo, sb, ep%eii, ep%fii)

    ! the pseudo potential part.
    do ia = 1, geo%natoms
      atm => geo%atom(ia)
      if(.not. specie_is_ps(atm%spec)) cycle
      ep%non_local = .true.
      call projector_end(ep%p(ia))
      call projector_init(ep%p(ia), gr%m, sb, atm, st%d%dim, ep%reltype)
      if(simul_box_is_periodic(sb)) call projector_init_phases(ep%p(ia), gr%m, st%d%nik, st%d%kpoints)
      call projector_build(ep%p(ia), gr, atm, ep%so_strength)
    end do

    if (ep%classic_pot > 0) then
      ep%vpsl(1:m%np) = ep%vpsl(1:m%np) + ep%vclassic(1:m%np)
    end if

    call pop_sub()
    call profiling_out(epot_generate_prof)

  end subroutine epot_generate

  subroutine epot_local_potential(ep, gr, geo, a, vpsl, time, rho_core)
    type(epot_t),             intent(in)    :: ep
    type(grid_t),             intent(inout) :: gr
    type(geometry_t),         intent(in)    :: geo
    type(atom_t),             intent(inout) :: a
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,                    intent(in)    :: time
    FLOAT,          optional, pointer       :: rho_core(:)

    integer :: i
    FLOAT :: x(MAX_DIM), radius
    FLOAT, allocatable  :: rho(:), vl(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    
    call push_sub('epot.epot_local_potential')
    call profiling_in(prof, "EPOT_LOCAL")

    ALLOCATE(vl(1:NP_PART), NP_PART)
#ifdef USE_OMP
    !$omp parallel workshare
    vl(1:NP) = M_ZERO
    !$omp end parallel workshare
#endif

    !Local potential, we can get it by solving the poisson equation
    !(for all electron species or pseudopotentials in periodic
    !systems) or by applying it directly to the grid

    if(a%spec%has_density .or. (specie_is_ps(a%spec) .and. simul_box_is_periodic(gr%sb))) then

      ALLOCATE(rho(1:NP), NP)

      !this has to be optimized so the poisson solution is made once
      !for all species, perhaps even include it in the hartree term
      call specie_get_density(a%spec, a%x, gr, geo, rho)

      vl(1:NP) = M_ZERO   ! vl has to be initialized before entering routine
      ! and our best guess for the potential is zero
      call dpoisson_solve(gr, vl, rho)

      deallocate(rho)

    else

      !Local potential
      call specie_get_local(a%spec, gr, a%x(1:NDIM), vl, time)

    end if

    !$omp parallel workshare
    vpsl(1:NP) = vpsl(1:NP) + vl(1:NP)
    !$omp end parallel workshare

    !the localized part
    if(specie_is_ps(a%spec)) then

      radius = double_grid_get_rmax(gr%dgrid, a%spec, gr%m) + gr%m%h(1)

      call submesh_init_sphere(sphere, gr%sb, gr%m, a%x, radius)
      call double_grid_apply_local(gr%dgrid, a%spec, gr%m, sphere, a%x, vl(1:sphere%ns))

      vpsl(sphere%jxyz(1:sphere%ns)) = vpsl(sphere%jxyz(1:sphere%ns)) + vl(1:sphere%ns)
      call submesh_end(sphere)

    end if

    deallocate(vl)

    !Non-local core corrections
    if(present(rho_core) .and. a%spec%nlcc .and. specie_is_ps(a%spec)) then
      do i = 1, NP
        x(1:NDIM) = gr%m%x(i, 1:NDIM) - a%x(1:NDIM)
        rho_core(i) = rho_core(i) + specie_get_nlcc(a%spec, x)
      end do
    end if

    call profiling_out(prof)
    call pop_sub()
  end subroutine epot_local_potential


  ! ---------------------------------------------------------
  subroutine epot_generate_classic(ep, m, geo)
    type(epot_t),     intent(inout) :: ep
    type(mesh_t),     intent(in)    :: m
    type(geometry_t), intent(in)    :: geo

    integer i, ia
    FLOAT :: r, rc

    call push_sub('epot.epot_generate_classic')

    ep%Vclassic = M_ZERO
    do ia = 1, geo%ncatoms
      do i = 1, m%np
        call mesh_r(m, i, r, a=geo%catom(ia)%x)
        select case(ep%classic_pot)
        case(1) ! point charge
          if(r < r_small) r = r_small
          ep%Vclassic(i) = ep%Vclassic(i) - geo%catom(ia)%charge/r
        case(2) ! gaussion smeared charge
          select case(geo%catom(ia)%label(1:1)) ! covalent radii
          case('H')
            rc = CNST(0.4)*P_Ang
          case('C')
            rc = CNST(0.8)*P_Ang
          case default
            rc = CNST(0.7)*P_Ang
          end select
          if(abs(r - rc) < r_small) r = rc + sign(r_small, r-rc)
          ep%Vclassic(i) = ep%Vclassic(i) - geo%catom(ia)%charge*(r**4 - rc**4)/(r**5 - rc**5)
        end select
      end do
    end do

    call pop_sub()
  end subroutine epot_generate_classic


  ! ---------------------------------------------------------
  subroutine epot_forces(gr, geo, ep, st, t)
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(inout) :: geo
    type(epot_t),     intent(in)    :: ep
    type(states_t),   intent(inout) :: st
    FLOAT,     optional, intent(in) :: t

    integer :: i, j
    FLOAT :: x(MAX_DIM), time
    
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    call push_sub('epot.epot_forces')

    time = M_ZERO
    if(present(t)) time = t

    ! the ion-ion term is already calculated
    do i = 1, geo%natoms
      geo%atom(i)%f(1:MAX_DIM) = ep%fii(1:MAX_DIM, i)
    end do
    
    if (wfs_are_real(st) ) then 
      call dcalc_forces_from_potential(gr, geo, ep, st, time)
    else
      call zcalc_forces_from_potential(gr, geo, ep, st, time)
    end if
    
    !TODO: forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, ep%no_lasers
        select case(ep%lasers(j)%field)
        case(E_FIELD_ELECTRIC)
          call laser_field(gr%sb, ep%lasers(j), x, t)
          do i = 1, geo%natoms
            geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + geo%atom(i)%spec%z_val*x(1:NDIM)
          end do

        case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL, E_FIELD_SCALAR_POTENTIAL)
          write(message(1),'(a)') 'The forces are currently not properly calculated if time-dependent'
          write(message(2),'(a)') 'magnetic fields are present.'
          call write_fatal(2)
        end select
      end do
    end if

    if(associated(ep%E_field)) then
      do i = 1, geo%natoms
        geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + geo%atom(i)%spec%Z_val * ep%E_field(1:NDIM)
      end do
    end if
    
    call pop_sub()
    call profiling_out(forces_prof)

  end subroutine epot_forces

#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot_m



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
