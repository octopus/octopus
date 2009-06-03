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
  use gauge_field_m
  use global_m
  use grid_m
  use io_m
  use ion_interaction_m
  use lalg_basic_m
  use lalg_adv_m
  use linear_response_m
  use loct_parser_m
  use splines_m
  use magnetic_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use multigrid_m
  use simul_box_m
  use units_m
  use logrid_m
  use poisson_cutoff_m
  use ps_m
  use smear_m
  use species_m
  use species_pot_m
  use spline_filter_m
  use solids_m
  use geometry_m
  use states_m
  use states_dim_m
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
  public ::                        &
    epot_t,                        &
    epot_init,                     &
    epot_end,                      &
    epot_generate,                 &
    epot_forces,                   &
    epot_local_potential,          &
    epot_precalc_local_potential,  &
    epot_dipole_periodic,          &
    dcalc_forces_from_potential,   &
    zcalc_forces_from_potential

  type epot_t
    ! Classical charges:
    integer        :: classical_pot ! how to include the classical charges
    FLOAT, pointer :: Vclassical(:) ! We use it to store the potential of the classical charges

    ! Ions
    FLOAT,             pointer :: vpsl(:)       ! the local part of the pseudopotentials
                                                ! plus the potential from static electric fields
    type(projector_t), pointer :: proj(:)       ! non-local projectors
    type(projector_t), pointer :: proj_fine(:)  ! non-local projectors in the fine grid
    logical                    :: non_local
    integer                    :: natoms

    ! External e-m fields
    integer                :: no_lasers            ! number of laser pulses used
    type(laser_t), pointer :: lasers(:)            ! lasers stuff
    FLOAT,         pointer :: E_field(:)           ! static electric field
    FLOAT, pointer         :: v_static(:)          ! static scalar potential
    FLOAT, pointer         :: B_field(:)           ! static magnetic field
    FLOAT, pointer         :: A_static(:,:)        ! static vector potential
    type(gauge_field_t)    :: gfield               ! the time-dependent gauge field
    integer                :: reltype              ! type of relativistic correction to use

    ! The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    ! *effective* electrons in a quantum dot. It affects the spin Zeeman term.
    FLOAT :: gyromagnetic_ratio

    ! SO prefactor (1.0 = normal SO, 0.0 = no SO)
    FLOAT :: so_strength
    
    ! the ion-ion energy and force
    FLOAT          :: eii
    FLOAT, pointer :: fii(:, :)
    
    real(4), pointer :: local_potential(:,:)
    logical        :: local_potential_precalculated
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
    call loct_parse_int(datasets_check('FilterPotentials'), PS_FILTER_NONE, filter)
    if(.not.varinfo_valid_option('FilterPotentials', filter)) call input_error('FilterPotentials')
    call messages_print_var_option(stdout, "FilterPotentials", filter)

    if(filter == PS_FILTER_TS) call spline_filter_mask_init()
    do i = 1, geo%nspecies
      call species_pot_init(geo%species(i), mesh_gcutoff(gr%mesh), filter)
    end do

    ! Local part of the pseudopotentials
    SAFE_ALLOCATE(ep%vpsl(1:gr%mesh%np))

    ep%vpsl(1:gr%mesh%np) = M_ZERO

    ep%classical_pot = 0
    if(geo%ncatoms > 0) then

      !%Variable ClassicalPotential
      !%Type integer
      !%Default 0
      !%Section Hamiltonian
      !%Description
      !% If <tt>true</tt>, add to the external potential the potential generated by 
      !% the point charges read from the PDB input (see <tt>PBDCoordinates</tt>).
      !%End
      call loct_parse_int(datasets_check('ClassicalPotential'), 0, ep%classical_pot)
      if(ep%classical_pot > 0) then
        message(1) = 'Info: generating classical external potential'
        call write_info(1)

        SAFE_ALLOCATE(ep%Vclassical(1:gr%mesh%np))
        call epot_generate_classical(ep, gr%mesh, geo)
      end if
    end if

    ! lasers
    call laser_init(ep%no_lasers, ep%lasers, gr%mesh)

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
    if(loct_parse_block(datasets_check('StaticElectricField'), blk)==0) then
      SAFE_ALLOCATE(ep%E_field(1:gr%mesh%sb%dim))
      do i = 1, gr%mesh%sb%dim
        call loct_parse_block_float(blk, 0, i-1, ep%E_field(i))
      end do
      call loct_parse_block_end(blk)

      ep%E_field(:) = ep%E_field(:) * units_inp%energy%factor/units_inp%length%factor

      ! Compute the scalar potential
      SAFE_ALLOCATE(ep%v_static(1:gr%mesh%np))
      do i = 1, gr%mesh%np
        ep%v_static(i) = sum(gr%mesh%x(i,:)*ep%E_field(:))
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
    if(loct_parse_block(datasets_check('StaticMagneticField'), blk)==0) then

      SAFE_ALLOCATE(ep%B_field(1:3))
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
      SAFE_ALLOCATE(ep%A_static(1:gr%mesh%np, 1:gr%mesh%sb%dim))
      SAFE_ALLOCATE(x(1:gr%mesh%sb%dim))
      do i = 1, gr%mesh%np
        x(1:gr%mesh%sb%dim) = gr%mesh%x(i, 1:gr%mesh%sb%dim)
        select case (gr%mesh%sb%dim)
        case (2)
          ep%A_static(i, :) = (/x(2), -x(1)/)*ep%B_field(3)
        case (3)
          ep%A_static(i, :) = (/x(2)*ep%B_field(3) - x(3)*ep%B_field(2), &
            x(3)*ep%B_field(1) - x(1)*ep%B_field(3), x(1)*ep%B_field(2) - x(2)*ep%B_field(1)/)
        end select
      end do
      SAFE_DEALLOCATE_A(x)
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
    call loct_parse_float(datasets_check('GyromagneticRatio'), P_g, ep%gyromagnetic_ratio)

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
    call loct_parse_int(datasets_check('RelativisticCorrection'), NOREL, ep%reltype)
    if(.not.varinfo_valid_option('RelativisticCorrection', ep%reltype)) call input_error('RelativisticCorrection')
    if (ispin /= SPINORS .and. ep%reltype == SPIN_ORBIT) then
      message(1) = "The spin-orbit term can only be applied when using spinors."
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
    if (ep%reltype == SPIN_ORBIT) then
      call loct_parse_float(datasets_check('SOStrength'), M_ONE, ep%so_strength)
    else
      ep%so_strength = M_ONE
    end if

    SAFE_ALLOCATE(ep%proj(1:geo%natoms))
    do i = 1, geo%natoms
      call projector_null(ep%proj(i))
    end do

    if(gr%have_fine_mesh) then
      SAFE_ALLOCATE(ep%proj_fine(1:geo%natoms))
      do i = 1, geo%natoms
        call projector_null(ep%proj_fine(i))
      end do
    else
      ep%proj_fine => ep%proj
    end if

    ep%natoms = geo%natoms
    ep%non_local = .false.

    SAFE_ALLOCATE(ep%fii(1:MAX_DIM, 1:geo%natoms))

    call gauge_field_init(ep%gfield, gr%sb)

    nullify(ep%local_potential)
    ep%local_potential_precalculated = .false.

    call pop_sub()
  end subroutine epot_init


  ! ---------------------------------------------------------
  subroutine epot_end(ep, gr, geo)
    type(epot_t),      intent(inout) :: ep
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(inout) :: geo

    integer :: i, iproj

    call push_sub('epot.epot_end')

    SAFE_DEALLOCATE_P(ep%local_potential)

    SAFE_DEALLOCATE_P(ep%fii)

    if(associated(ep%vpsl)) then
      SAFE_DEALLOCATE_P(ep%vpsl)
      nullify(ep%vpsl)
    end if

    if(ep%classical_pot > 0) then
      ep%classical_pot = 0
      ! sanity check
      ASSERT(associated(ep%Vclassical)) 
      SAFE_DEALLOCATE_P(ep%Vclassical)         ! and clean up
      nullify(ep%Vclassical)
    end if

    ! the external laser
    call laser_end(ep%no_lasers, ep%lasers)

    ! the macroscopic fields
    SAFE_DEALLOCATE_P(ep%E_field)
    SAFE_DEALLOCATE_P(ep%v_static)
    SAFE_DEALLOCATE_P(ep%B_field)
    SAFE_DEALLOCATE_P(ep%A_static)

    do iproj = 1, geo%natoms
      if(.not. species_is_ps(geo%atom(iproj)%spec)) cycle
      call projector_end(ep%proj(iproj))
    end do

    ASSERT(associated(ep%proj))
    SAFE_DEALLOCATE_P(ep%proj)

    if(gr%have_fine_mesh) then
      do iproj = 1, geo%natoms
        if(.not. species_is_ps(geo%atom(iproj)%spec)) cycle
        call projector_end(ep%proj_fine(iproj))
      end do
      SAFE_DEALLOCATE_P(ep%proj_fine)
    end if

    call pop_sub()

  end subroutine epot_end

  ! ---------------------------------------------------------

  subroutine epot_generate(ep, gr, geo, st, time)
    type(epot_t),          intent(inout) :: ep
    type(grid_t), target,  intent(inout) :: gr
    type(geometry_t),      intent(inout) :: geo
    type(states_t),        intent(inout) :: st
    FLOAT,       optional, intent(in)    :: time

    FLOAT   :: time_
    integer :: ia
    type(atom_t),      pointer :: atm
    type(mesh_t),      pointer :: mesh
    type(simul_box_t), pointer :: sb
    type(profile_t), save :: epot_generate_prof

#ifdef HAVE_MPI
    FLOAT,    allocatable :: vpsltmp(:)
    type(profile_t), save :: epot_reduce
#endif

    call profiling_in(epot_generate_prof, "EPOT_GENERATE")
    call push_sub('epot.epot_generate')

    sb   => gr%sb
    mesh => gr%mesh

    time_ = M_ZERO
    if (present(time)) time_ = time

    ! Local.
    ep%vpsl = M_ZERO
    if(st%nlcc) st%rho_core = M_ZERO
    do ia = geo%atoms%start, geo%atoms%end
      if(st%nlcc) then
        call epot_local_potential(ep, gr, gr%mesh, geo, ia, ep%vpsl, time_, st%rho_core)
      else
        call epot_local_potential(ep, gr, gr%mesh, geo, ia, ep%vpsl, time_)
      end if
    end do

#ifdef HAVE_MPI
    call profiling_in(epot_reduce, "EPOT_REDUCE")
    if(geo%atoms%parallel) then
      SAFE_ALLOCATE(vpsltmp(1:gr%mesh%np))
      call MPI_Allreduce(ep%vpsl, vpsltmp, gr%mesh%np, MPI_FLOAT, MPI_SUM, geo%atoms%mpi_grp%comm, mpi_err)
      call lalg_copy(gr%mesh%np, vpsltmp, ep%vpsl)
      if(associated(st%rho_core)) then
        call MPI_Allreduce(st%rho_core, vpsltmp, gr%mesh%np, MPI_FLOAT, MPI_SUM, geo%atoms%mpi_grp%comm, mpi_err)
        call lalg_copy(gr%mesh%np, vpsltmp, st%rho_core)
      end if
      SAFE_DEALLOCATE_A(vpsltmp)
    end if
    call profiling_out(epot_reduce)
#endif

    ! we assume that we need to recalculate the ion-ion energy
    call ion_interaction_calculate(geo, sb, ep%eii, ep%fii)

    ! the pseudopotential part.
    do ia = 1, geo%natoms
      atm => geo%atom(ia)
      if(.not. species_is_ps(atm%spec)) cycle
      ep%non_local = .true.
      call projector_end(ep%proj(ia))
      call projector_init(ep%proj(ia), gr%mesh, sb, atm, st%d%dim, ep%reltype)
      if(simul_box_is_periodic(sb) .or. associated(ep%a_static)) then
        call projector_init_phases(ep%proj(ia), sb, st%d%nik, st%d%kpoints, vec_pot_var = ep%a_static)
      end if
      call projector_build(ep%proj(ia), gr, atm, ep%so_strength)

      ! the projectors in the fine grid
      if(gr%have_fine_mesh) then
        call projector_end(ep%proj_fine(ia))
        call projector_init(ep%proj_fine(ia), gr%fine%mesh, sb, atm, st%d%dim, ep%reltype)
        if(simul_box_is_periodic(sb) .or. associated(ep%a_static)) then
          call projector_init_phases(ep%proj_fine(ia), sb, st%d%nik, st%d%kpoints, vec_pot_var = ep%a_static)
        end if
        call projector_build(ep%proj_fine(ia), gr, atm, ep%so_strength)
      end if

    end do

    ! add static electric fields
    if (ep%classical_pot > 0)     ep%vpsl(1:mesh%np) = ep%vpsl(1:mesh%np) + ep%Vclassical(1:mesh%np)
    if (associated(ep%e_field)) ep%vpsl(1:mesh%np) = ep%vpsl(1:mesh%np) + ep%v_static(1:mesh%np)

    call pop_sub()
    call profiling_out(epot_generate_prof)
  end subroutine epot_generate

  subroutine epot_local_potential(ep, gr, mesh, geo, iatom, vpsl, time, rho_core)
    type(epot_t),             intent(in)    :: ep
    type(grid_t),             intent(inout) :: gr
    type(mesh_t),             intent(inout) :: mesh
    type(geometry_t),         intent(in)    :: geo
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,                    intent(in)    :: time
    FLOAT,          optional, pointer       :: rho_core(:)

    integer :: i, ip
    FLOAT :: x(MAX_DIM), radius
    FLOAT, allocatable  :: rho(:), vl(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    
    call push_sub('epot.epot_local_potential')
    call profiling_in(prof, "EPOT_LOCAL")

    if(ep%local_potential_precalculated) then

      forall(ip = 1:mesh%np) vpsl(ip) = vpsl(ip) + ep%local_potential(ip, iatom)

    else
      
      SAFE_ALLOCATE(vl(1:mesh%np_part))
      
      !Local potential, we can get it by solving the Poisson equation
      !(for all-electron species or pseudopotentials in periodic
      !systems) or by applying it directly to the grid
      
      if( species_has_density(geo%atom(iatom)%spec) .or. &
          (species_is_ps(geo%atom(iatom)%spec) .and. simul_box_is_periodic(gr%sb)) ) then

        SAFE_ALLOCATE(rho(1:mesh%np))

        !this has to be optimized so the Poisson solution is made once
        !for all species, perhaps even include it in the Hartree term
        call species_get_density(geo%atom(iatom)%spec, geo%atom(iatom)%x, gr, geo, rho)

        vl(1:mesh%np) = M_ZERO   ! vl has to be initialized before entering routine
        ! and our best guess for the potential is zero
        call dpoisson_solve(gr, vl, rho)

        SAFE_DEALLOCATE_A(rho)
      else

        !Local potential
        call species_get_local(geo%atom(iatom)%spec, mesh, geo%atom(iatom)%x(1:gr%mesh%sb%dim), vl, time)
      end if
      
      vpsl(1:mesh%np) = vpsl(1:mesh%np) + vl(1:mesh%np)

      !the localized part
      if(species_is_ps(geo%atom(iatom)%spec)) then
        
        radius = double_grid_get_rmax(gr%dgrid, geo%atom(iatom)%spec, mesh) + mesh%h(1)
        
        call submesh_init_sphere(sphere, gr%sb, mesh, geo%atom(iatom)%x, radius)
        call double_grid_apply_local(gr%dgrid, geo%atom(iatom)%spec, mesh, sphere, geo%atom(iatom)%x, vl(1:sphere%ns))
        
        vpsl(sphere%jxyz(1:sphere%ns)) = vpsl(sphere%jxyz(1:sphere%ns)) + vl(1:sphere%ns)
        call submesh_end(sphere)
        
      end if

      SAFE_DEALLOCATE_A(vl)
    end if

    !Non-local core corrections
    if( present(rho_core) .and. &
        species_has_nlcc(geo%atom(iatom)%spec) .and. &
        species_is_ps(geo%atom(iatom)%spec)) then
      do i = 1, mesh%np
        x(1:gr%mesh%sb%dim) = mesh%x(i, 1:gr%mesh%sb%dim) - geo%atom(iatom)%x(1:gr%mesh%sb%dim)
        rho_core(i) = rho_core(i) + species_get_nlcc(geo%atom(iatom)%spec, x)
      end do
    end if

    call profiling_out(prof)
    call pop_sub()
  end subroutine epot_local_potential


  ! ---------------------------------------------------------
  subroutine epot_generate_classical(ep, m, geo)
    type(epot_t),     intent(inout) :: ep
    type(mesh_t),     intent(in)    :: m
    type(geometry_t), intent(in)    :: geo

    integer i, ia
    FLOAT :: r, rc

    call push_sub('epot.epot_generate_classical')

    ep%Vclassical = M_ZERO
    do ia = 1, geo%ncatoms
      do i = 1, m%np
        call mesh_r(m, i, r, a=geo%catom(ia)%x)
        select case(ep%classical_pot)
        case(1) ! point charge
          if(r < r_small) r = r_small
          ep%Vclassical(i) = ep%Vclassical(i) - geo%catom(ia)%charge/r
        case(2) ! Gaussian smeared charge
          select case(geo%catom(ia)%label(1:1)) ! covalent radii
          case('H')
            rc = CNST(0.4)*P_Ang
          case('C')
            rc = CNST(0.8)*P_Ang
          case default
            rc = CNST(0.7)*P_Ang
          end select
          if(abs(r - rc) < r_small) r = rc + sign(r_small, r-rc)
          ep%Vclassical(i) = ep%Vclassical(i) - geo%catom(ia)%charge*(r**4 - rc**4)/(r**5 - rc**5)
        end select
      end do
    end do

    call pop_sub()
  end subroutine epot_generate_classical


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
    
    if (states_are_real(st) ) then 
      call dcalc_forces_from_potential(gr, geo, ep, st, time)
    else
      call zcalc_forces_from_potential(gr, geo, ep, st, time)
    end if
    
    !TODO: forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, ep%no_lasers
        select case(laser_kind(ep%lasers(j)))
        case(E_FIELD_ELECTRIC)
          call laser_field(gr%sb, ep%lasers(j), x, t)
          do i = 1, geo%natoms
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) + &
              P_PROTON_CHARGE * species_zval(geo%atom(i)%spec) * x(1:gr%mesh%sb%dim)
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
        geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) + &
          P_PROTON_CHARGE * species_zval(geo%atom(i)%spec) * ep%E_field(1:gr%mesh%sb%dim)
      end do
    end if
    
    call pop_sub()
    call profiling_out(forces_prof)

  end subroutine epot_forces

  ! ---------------------------------------------------------
  ! Uses the single-point Berry`s phase method to calculate dipole moment in a periodic system
  ! This is only accurate in the limit of a large supercell.
  ! It is implemented only for an orthogonal unit cell.
  ! mu = - eL/2*pi Im ln <Psi|exp(-i(2*pi/L)x)|Psi>
  ! E Yaschenko, L Fu, L Resca, R Resta, Phys. Rev. B 58, 1222-1229 (1998)
  ! Single-point Berry`s phase method for dipole should not be used when there is more than one k-point.
  ! in this case, finite differences should be used to construct derivatives with respect to k
  FLOAT function epot_dipole_periodic(st, gr, dir) result(dipole)
    type(states_t), intent(in) :: st
    type(grid_t),   intent(in) :: gr
    integer,        intent(in) :: dir

    integer ik, ist, ist2, idim, ip
    CMPLX, allocatable :: matrix(:, :), tmp(:)
    CMPLX :: det, phase

    call push_sub('epot.epot_dipole_periodic')

    if(.not. smear_is_semiconducting(st%smear)) then
      message(1) = "Warning: single-point Berry's phase dipole calculation not correct without integer occupations."
      call write_warning(1)
    endif

    SAFE_ALLOCATE(matrix(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(tmp(1:gr%mesh%np))

    ! TODO: add in sum over k-points in orthogonal directions here

    phase = exp(-M_zI*(M_PI/gr%sb%lsize(dir)))

    do ik = st%d%kpt%start, st%d%kpt%end ! determinants for different spins multiply since matrix is block-diagonal
      do ist = 1, st%nst
        do ist2 = 1, st%nst

          if (ist .le. ist2) then
            matrix(ist, ist2) = M_Z0
            do idim = 1, st%d%dim ! spinor components

              if(states_are_complex(st)) then
                forall(ip = 1:gr%mesh%np)
                  tmp(ip) = conjg(st%zpsi(ip, idim, ist, ik))*&
                       exp(-M_zI*(M_PI/gr%sb%lsize(dir))*gr%mesh%x(ip, dir))*st%zpsi(ip, idim, ist2, ik)
                end forall
              else
                forall(ip = 1:gr%mesh%np)
                  tmp(ip) = st%dpsi(ip, idim, ist, ik)*&
                       exp(-M_zI*(M_PI/gr%sb%lsize(dir))*gr%mesh%x(ip, dir))*st%dpsi(ip, idim, ist2, ik)
                end forall
              end if

              matrix(ist, ist2) = matrix(ist, ist2) + zmf_integrate(gr%mesh, tmp)

              ! factor of two removed from exp since actual lattice vector is 2 * lsize
            end do
          else
            ! enforce Hermiticity of matrix
            matrix(ist, ist2) = conjg(matrix(ist2, ist))
          end if
        enddo !ist2
      enddo !ist
    enddo !ik
      
    det = lalg_determinant(st%nst, matrix(1:st%nst, 1:st%nst), invert = .false.)
    dipole = -(2 * gr%sb%lsize(dir) / (2 * M_Pi)) * aimag(log(det)) * st%smear%el_per_state

    SAFE_DEALLOCATE_A(matrix)
    call pop_sub()
  end function epot_dipole_periodic

  subroutine epot_precalc_local_potential(ep, gr, mesh, geo, time)
    type(epot_t),             intent(inout) :: ep
    type(grid_t),             intent(inout) :: gr
    type(mesh_t),             intent(inout) :: mesh
    type(geometry_t),         intent(inout) :: geo
    FLOAT,                    intent(in)    :: time

    integer :: iatom
    FLOAT, allocatable :: tmp(:)

    if(.not. associated(ep%local_potential)) then
      SAFE_ALLOCATE(ep%local_potential(1:mesh%np, 1:geo%natoms))
    end if

    ep%local_potential_precalculated = .false.

    SAFE_ALLOCATE(tmp(1:mesh%np))

    do iatom = 1, geo%natoms
      tmp(1:mesh%np) = M_ZERO
      call epot_local_potential(ep, gr, mesh, geo, iatom, tmp, time)
      ep%local_potential(1:mesh%np, iatom) = tmp(1:mesh%np)
    end do

    ep%local_potential_precalculated = .true.

  end subroutine epot_precalc_local_potential

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
