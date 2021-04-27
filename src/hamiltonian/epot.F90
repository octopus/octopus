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

module epot_oct_m
  use atom_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use gauge_field_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use ion_interaction_oct_m
  use kick_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use species_pot_oct_m
  use splines_oct_m
  use spline_filter_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use tdfunction_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  private
  public ::                        &
    epot_t,                        &
    epot_init,                     &
    epot_end,                      &
    epot_generate,                 &
    epot_local_potential,          &
    epot_precalc_local_potential,  &
    epot_have_kick,                &
    epot_have_external_potentials, &
    local_potential_has_density

  integer, public, parameter :: &
    CLASSICAL_NONE     = 0, & !< no classical charges
    CLASSICAL_POINT    = 1, & !< classical charges treated as point charges
    CLASSICAL_GAUSSIAN = 2    !< classical charges treated as Gaussian distributions
     
  integer, public, parameter :: &
    NOREL      = 0,             &
    SPIN_ORBIT = 1

  type epot_t
    ! Components are public by default

    ! Classical charges:
    integer            :: classical_pot !< how to include the classical charges
    FLOAT, allocatable :: Vclassical(:) !< We use it to store the potential of the classical charges

    ! Ions
    FLOAT,             allocatable :: vpsl(:)       !< the local part of the pseudopotentials
                                                !< plus the potential from static electric fields
    type(projector_t), allocatable :: proj(:)       !< non-local projectors
    logical                        :: non_local
    integer                        :: natoms

    ! External e-m fields
    FLOAT, allocatable     :: E_field(:)           !< static electric field
    FLOAT, allocatable     :: v_ext(:)             !< static scalar potential - 1:gr%mesh%np_part
    FLOAT, allocatable     :: B_field(:)           !< static magnetic field
    FLOAT, allocatable     :: A_static(:,:)        !< static vector potential
    type(gauge_field_t)    :: gfield               !< the time-dependent gauge field
    integer                :: reltype              !< type of relativistic correction to use

    !> The possible kick
    type(kick_t) :: kick

    !> The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    !! *effective* electrons in a quantum dot. It affects the spin Zeeman term.)
    FLOAT :: gyromagnetic_ratio

    !> SO prefactor (1.0 = normal SO, 0.0 = no SO)
    FLOAT, private :: so_strength
    
    !> the ion-ion energy and force
    FLOAT              :: eii
    FLOAT, allocatable :: fii(:, :)
    FLOAT, allocatable :: vdw_forces(:, :)
    
    FLOAT, allocatable, private :: local_potential(:,:)
    logical,            private :: local_potential_precalculated

    logical,                  private :: have_density
    type(poisson_t), pointer, private :: poisson_solver

    logical              :: nlcc = .false.   !< does any species have non-local core corrections?
  end type epot_t

contains

  ! ---------------------------------------------------------
  subroutine epot_init(ep, namespace, gr, geo, psolver, ispin, xc_family, mc, kpoints)
    type(epot_t),                       intent(out)   :: ep
    type(namespace_t),                  intent(in)    :: namespace
    type(grid_t),                       intent(in)    :: gr
    type(geometry_t),                   intent(inout) :: geo
    type(poisson_t),  target,           intent(in)    :: psolver
    integer,                            intent(in)    :: ispin
    integer,                            intent(in)    :: xc_family
    type(multicomm_t),                  intent(in)    :: mc
    type(kpoints_t),                    intent(in)    :: kpoints


    integer :: ispec, ia
    integer :: filter

    PUSH_SUB(epot_init)

    !%Variable FilterPotentials
    !%Type integer
    !%Default filter_ts
    !%Section Hamiltonian
    !%Description
    !% <tt>Octopus</tt> can filter the pseudopotentials so that they no
    !% longer contain Fourier components larger than the mesh itself. This is
    !% very useful to decrease the egg-box effect, and so should be used in
    !% all instances where atoms move (<i>e.g.</i> geometry optimization,
    !% molecular dynamics, and vibrational modes).
    !%Option filter_none 0
    !% Do not filter.
    !%Option filter_TS 2
    !% The filter of M. Tafipolsky and R. Schmid, <i>J. Chem. Phys.</i> <b>124</b>, 174102 (2006).
    !%Option filter_BSB 3
    !% The filter of E. L. Briggs, D. J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b>, 14362 (1996).
    !%End
    call parse_variable(namespace, 'FilterPotentials', PS_FILTER_TS, filter)
    if(.not.varinfo_valid_option('FilterPotentials', filter)) call messages_input_error(namespace, 'FilterPotentials')
    call messages_print_var_option(stdout, "FilterPotentials", filter)

    if(family_is_mgga(xc_family) .and. filter /= PS_FILTER_NONE) &
      call messages_not_implemented("FilterPotentials different from filter_none with MGGA", namespace=namespace)

    if(filter == PS_FILTER_TS) call spline_filter_mask_init()
    do ispec = 1, geo%nspecies
      call species_pot_init(geo%species(ispec), namespace, mesh_gcutoff(gr%mesh), filter)
    end do

    SAFE_ALLOCATE(ep%vpsl(1:gr%mesh%np))

    ep%vpsl(1:gr%mesh%np) = M_ZERO

    ep%classical_pot = 0
    if(geo%ncatoms > 0) then

      if (geo%space%is_periodic()) then
        call messages_not_implemented("classical atoms in periodic systems", namespace=namespace)
      end if

      !%Variable ClassicalPotential
      !%Type integer
      !%Default no
      !%Section Hamiltonian
      !%Description
      !% Whether and how to add to the external potential the potential generated by 
      !% the classical charges read from block <tt>PDBClassical</tt>, for QM/MM calculations.
      !% Not available in periodic systems.
      !%Option no 0
      !%  No classical charges.
      !%Option point_charges 1
      !%  Classical charges are treated as point charges.
      !%Option gaussian_smeared 2
      !%  Classical charges are treated as Gaussian distributions. 
      !%  Smearing widths are hard-coded by species (experimental).
      !%End
      call parse_variable(namespace, 'ClassicalPotential', 0, ep%classical_pot)
      if(ep%classical_pot  ==  CLASSICAL_GAUSSIAN) then
        call messages_experimental("Gaussian smeared classical charges")
        ! This method probably works but definitely needs to be made user-friendly:
        ! i.e. telling the user what widths are used and letting them be set somehow.
      end if

      if(ep%classical_pot > 0) then
        message(1) = 'Info: generating classical external potential.'
        call messages_info(1)

        SAFE_ALLOCATE(ep%Vclassical(1:gr%mesh%np))
        call epot_generate_classical(ep, gr%mesh, geo)
      end if
    end if

    call kick_init(ep%kick, namespace, gr%sb, kpoints, ispin)

    ! No more "UserDefinedTDPotential" from this version on.
    call messages_obsolete_variable(namespace, 'UserDefinedTDPotential', 'TDExternalFields')

    !%Variable GyromagneticRatio
    !%Type float
    !%Default 2.0023193043768
    !%Section Hamiltonian
    !%Description
    !% The gyromagnetic ratio of the electron. This is of course a physical 
    !% constant, and the default value is the exact one that you should not 
    !% touch, unless: 
    !% (i)  You want to disconnect the anomalous Zeeman term in the Hamiltonian 
    !% (then set it to zero; this number only affects that term);
    !% (ii) You are using an effective Hamiltonian, as is the case when
    !% you calculate a 2D electron gas, in which case you have an effective
    !% gyromagnetic factor that depends on the material.
    !%End
    call parse_variable(namespace, 'GyromagneticRatio', P_g, ep%gyromagnetic_ratio)

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
    !% Spin-orbit.
    !%End
    call parse_variable(namespace, 'RelativisticCorrection', NOREL, ep%reltype)
    if(.not.varinfo_valid_option('RelativisticCorrection', ep%reltype)) then
      call messages_input_error(namespace, 'RelativisticCorrection')
    end if
    if (ispin /= SPINORS .and. ep%reltype == SPIN_ORBIT) then
      message(1) = "The spin-orbit term can only be applied when using spinors."
      call messages_fatal(1, namespace=namespace)
    end if

    call messages_print_var_option(stdout, "RelativisticCorrection", ep%reltype)

    !%Variable SOStrength
    !%Type float
    !%Default 1.0
    !%Section Hamiltonian
    !%Description
    !% Tuning of the spin-orbit coupling strength: setting this value to zero turns off spin-orbit terms in
    !% the Hamiltonian, and setting it to one corresponds to full spin-orbit.
    !%End
    if (ep%reltype == SPIN_ORBIT) then
      call parse_variable(namespace, 'SOStrength', M_ONE, ep%so_strength)
    else
      ep%so_strength = M_ONE
    end if

    SAFE_ALLOCATE(ep%proj(1:geo%natoms))

    ep%natoms = geo%natoms
    ep%non_local = .false.

    ep%eii = M_ZERO
    SAFE_ALLOCATE(ep%fii(1:gr%sb%dim, 1:geo%natoms))
    ep%fii = M_ZERO

    SAFE_ALLOCATE(ep%vdw_forces(1:geo%space%dim, 1:geo%natoms))
    ep%vdw_forces = M_ZERO

    ep%local_potential_precalculated = .false.
    

    ep%have_density = .false.
    do ia = 1, geo%natoms
      if (local_potential_has_density(geo%space, geo%atom(ia))) then
        ep%have_density = .true.
        exit
      end if
    end do

    if(ep%have_density) then
      ep%poisson_solver => psolver
    else
      nullify(ep%poisson_solver)
    end if

    ! find out if we need non-local core corrections
    ep%nlcc = .false.
    do ia = 1, geo%nspecies
      ep%nlcc = (ep%nlcc.or.species_has_nlcc(geo%species(ia)))
    end do

    POP_SUB(epot_init)
  end subroutine epot_init

  ! ---------------------------------------------------------
  subroutine epot_end(ep)
    type(epot_t), intent(inout) :: ep

    integer :: iproj

    PUSH_SUB(epot_end)

    if(ep%have_density) then
      nullify(ep%poisson_solver)
    end if

    SAFE_DEALLOCATE_A(ep%local_potential)
    SAFE_DEALLOCATE_A(ep%fii)
    SAFE_DEALLOCATE_A(ep%vdw_forces)
    SAFE_DEALLOCATE_A(ep%vpsl)

    if(ep%classical_pot > 0) then
      ep%classical_pot = 0
      ! sanity check
      ASSERT(allocated(ep%Vclassical)) 
      SAFE_DEALLOCATE_A(ep%Vclassical)         ! and clean up
    end if


    call kick_end(ep%kick)

    ! the macroscopic fields
    SAFE_DEALLOCATE_A(ep%E_field)
    SAFE_DEALLOCATE_A(ep%v_ext)
    SAFE_DEALLOCATE_A(ep%B_field)
    SAFE_DEALLOCATE_A(ep%A_static)

    do iproj = 1, ep%natoms
      if (projector_is_null(ep%proj(iproj))) cycle
      call projector_end(ep%proj(iproj))
    end do

    ASSERT(allocated(ep%proj))
    SAFE_DEALLOCATE_A(ep%proj)

    POP_SUB(epot_end)

  end subroutine epot_end

  ! ---------------------------------------------------------
  subroutine epot_generate(ep, namespace, gr, geo, st_d)
    type(epot_t),             intent(inout) :: ep
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),     target, intent(in)    :: gr
    type(geometry_t), target, intent(inout) :: geo
    type(states_elec_dim_t),  intent(inout) :: st_d

    integer :: ia
    type(mesh_t),      pointer :: mesh
    type(simul_box_t), pointer :: sb
    type(profile_t), save :: epot_generate_prof
    type(ps_t), pointer :: ps
    
    call profiling_in(epot_generate_prof, "EPOT_GENERATE")
    PUSH_SUB(epot_generate)

    sb   => gr%sb
    mesh => gr%mesh

    ! Local part
    ep%vpsl = M_ZERO

    ! we assume that we need to recalculate the ion-ion energy
    call ion_interaction_calculate(geo%ion_interaction, geo%space, geo%latt, geo%atom, geo%natoms, geo%catom, geo%ncatoms, &
      sb%lsize, ep%eii, ep%fii)

    ! the pseudopotential part.
    do ia = 1, geo%natoms
      if (.not. species_is_ps(geo%atom(ia)%species)) cycle
      call projector_end(ep%proj(ia))
      call projector_init(ep%proj(ia), geo%atom(ia), namespace, st_d%dim, ep%reltype)
    end do

    do ia = geo%atoms_dist%start, geo%atoms_dist%end
      if(ep%proj(ia)%type == PROJ_NONE) cycle
      ps => species_ps(geo%atom(ia)%species)
      call submesh_init(ep%proj(ia)%sphere, geo%space, mesh%sb, mesh, geo%atom(ia)%x, ps%rc_max + mesh%spacing(1))
    end do

    if(geo%atoms_dist%parallel) then
      do ia = 1, geo%natoms
        if(ep%proj(ia)%type == PROJ_NONE) cycle
        ps => species_ps(geo%atom(ia)%species)
        call submesh_broadcast(ep%proj(ia)%sphere, mesh, geo%atom(ia)%x, ps%rc_max + mesh%spacing(1), &
          geo%atoms_dist%node(ia), geo%atoms_dist%mpi_grp)
      end do
    end if
    
    do ia = 1, geo%natoms
      call projector_build(ep%proj(ia), geo%atom(ia), ep%so_strength)
      if(.not. projector_is(ep%proj(ia), PROJ_NONE)) ep%non_local = .true.
    end do

    POP_SUB(epot_generate)
    call profiling_out(epot_generate_prof)
  end subroutine epot_generate

  ! ---------------------------------------------------------

  logical pure function local_potential_has_density(space, atom) result(has_density)
    type(space_t),        intent(in)    :: space
    type(atom_t),         intent(in)    :: atom
    
    has_density = species_has_density(atom%species) .or. (species_is_ps(atom%species) .and. space%is_periodic())

  end function local_potential_has_density
  
  ! ---------------------------------------------------------
  subroutine epot_local_potential(ep, namespace, space, mesh, atom, iatom, vpsl)
    type(epot_t),             intent(in)    :: ep
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(mesh_t),             intent(in)    :: mesh
    type(atom_t),             intent(in)    :: atom
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)

    integer :: ip
    FLOAT :: radius, r
    FLOAT, allocatable :: vl(:), rho(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    type(ps_t), pointer :: ps

    PUSH_SUB(epot_local_potential)
    call profiling_in(prof, "EPOT_LOCAL")

    if(ep%local_potential_precalculated) then

      do ip = 1, mesh%np
        vpsl(ip) = vpsl(ip) + ep%local_potential(ip, iatom)
      end do
    else

      !Local potential, we can get it by solving the Poisson equation
      !(for all-electron species or pseudopotentials in periodic
      !systems) or by applying it directly to the grid

      if(local_potential_has_density(space, atom)) then
        SAFE_ALLOCATE(rho(1:mesh%np))

        call species_get_long_range_density(atom%species, space, namespace, atom%x, mesh, rho)

        SAFE_ALLOCATE(vl(1:mesh%np))
          
        if(poisson_solver_is_iterative(ep%poisson_solver)) then
          ! vl has to be initialized before entering routine
          ! and our best guess for the potential is zero
          vl(1:mesh%np) = M_ZERO
        end if

        call dpoisson_solve(ep%poisson_solver, vl, rho, all_nodes = .false.)

        SAFE_DEALLOCATE_A(rho)

      else

        SAFE_ALLOCATE(vl(1:mesh%np))
        call species_get_local(atom%species, space, mesh, namespace, atom%x(1:mesh%sb%dim), vl)

      end if

      if(allocated(vl)) then
        call lalg_axpy(mesh%np, M_ONE, vl, vpsl)
        SAFE_DEALLOCATE_A(vl)
      end if

      !the localized part
      if(species_is_ps(atom%species)) then

        ps => species_ps(atom%species)

        radius = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold) + mesh%spacing(1)

        call submesh_init(sphere, space, mesh%sb, mesh, atom%x, radius)
        SAFE_ALLOCATE(vl(1:sphere%np))

        do ip = 1, sphere%np
          r = sphere%x(ip, 0)
          vl(ip) = spline_eval(ps%vl, r)
        end do

        call submesh_add_to_mesh(sphere, vl, vpsl)

        SAFE_DEALLOCATE_A(vl)
        call submesh_end(sphere)
        nullify(ps)

      end if

    end if

    call profiling_out(prof)
    POP_SUB(epot_local_potential)
  end subroutine epot_local_potential


  ! ---------------------------------------------------------
  subroutine epot_generate_classical(ep, mesh, geo)
    type(epot_t),     intent(inout) :: ep
    type(mesh_t),     intent(in)    :: mesh
    type(geometry_t), intent(in)    :: geo

    integer ip, ia
    FLOAT :: rr, rc

    PUSH_SUB(epot_generate_classical)

    ep%Vclassical = M_ZERO
    do ia = 1, geo%ncatoms
      do ip = 1, mesh%np
        call mesh_r(mesh, ip, rr, origin=geo%catom(ia)%x)
        select case(ep%classical_pot)
        case(CLASSICAL_POINT)
          if(rr < R_SMALL) rr = R_SMALL
          ep%Vclassical(ip) = ep%Vclassical(ip) - geo%catom(ia)%charge/rr
        case(CLASSICAL_GAUSSIAN)
          select case(geo%catom(ia)%label(1:1)) ! covalent radii
          case('H')
            rc = CNST(0.4) * P_Ang
          case('C')
            rc = CNST(0.8) * P_Ang
          case default
            rc = CNST(0.7) * P_Ang
          end select
          if(abs(rr - rc) < R_SMALL) rr = rc + sign(R_SMALL, rr - rc)
          ep%Vclassical(ip) = ep%Vclassical(ip) - geo%catom(ia)%charge * (rr**4 - rc**4) / (rr**5 - rc**5)
        case default
          message(1) = 'Unknown type of classical potential in epot_generate_classical'
          call messages_fatal(1)
        end select
      end do
    end do

    POP_SUB(epot_generate_classical)
  end subroutine epot_generate_classical

  ! ---------------------------------------------------------
  subroutine epot_precalc_local_potential(ep, namespace, gr, geo)
    type(epot_t),      intent(inout) :: ep
    type(namespace_t), intent(in)    :: namespace
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(in)    :: geo

    integer :: iatom
    
    PUSH_SUB(epot_precalc_local_potential)
    
    if(.not. allocated(ep%local_potential)) then
      SAFE_ALLOCATE(ep%local_potential(1:gr%mesh%np, 1:geo%natoms))
    end if

    ep%local_potential_precalculated = .false.

    do iatom = 1, geo%natoms
      ep%local_potential(1:gr%mesh%np, iatom) = M_ZERO 
      call epot_local_potential(ep, namespace, geo%space, gr%mesh, geo%atom(iatom), iatom, &
                ep%local_potential(1:gr%mesh%np, iatom))!, time)
    end do
    ep%local_potential_precalculated = .true.

    POP_SUB(epot_precalc_local_potential)
  end subroutine epot_precalc_local_potential

  ! ---------------------------------------------------------

  logical function epot_have_kick(ep)
    type(epot_t), intent(in)  :: ep

    PUSH_SUB(epot_have_kick)

    epot_have_kick = .false.

    if( ep%kick%delta_strength /= M_ZERO ) epot_have_kick = .true.

    POP_SUB(epot_have_kick)

  end function epot_have_kick

  ! ---------------------------------------------------------

  logical function epot_have_external_potentials(ep)
    type(epot_t), intent(in)  :: ep

    PUSH_SUB(epot_have_external_potentials)

    epot_have_external_potentials =  .false.

    if( allocated(ep%E_field) ) epot_have_external_potentials = .true.

    POP_SUB(epot_have_external_potentials)

  end function epot_have_external_potentials

end module epot_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
