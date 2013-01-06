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

module epot_m
  use atom_m
  use comm_m
  use datasets_m
  use derivatives_m
  use double_grid_m
  use gauge_field_m
  use geometry_m
  use global_m
  use grid_m
  use index_m
  use io_m
  use kick_m
  use lalg_adv_m
  use lalg_basic_m
  use lasers_m
  use linear_response_m
  use loct_math_m
  use logrid_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multigrid_m
  use parser_m
  use poisson_m
  use poisson_cutoff_m
  use poisson_sete_m
  use profiling_m
  use projector_m
  use ps_m
  use simul_box_m
  use solids_m
  use species_m
  use species_pot_m
  use splines_m
  use spline_filter_m
  use states_m
  use states_dim_m
  use submesh_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                        &
    epot_t,                        &
    epot_init,                     &
    epot_end,                      &
    epot_generate,                 &
    epot_local_potential,          &
    epot_precalc_local_potential

  integer, public, parameter :: &
    CLASSICAL_NONE     = 0, & !< no classical charges
    CLASSICAL_POINT    = 1, & !< classical charges treated as point charges
    CLASSICAL_GAUSSIAN = 2    !< classical charges treated as Gaussian distributions
     
  integer, public, parameter :: &
    NOREL      = 0,             &
    SPIN_ORBIT = 1

  type epot_t
    ! Classical charges:
    integer        :: classical_pot !< how to include the classical charges
    FLOAT, pointer :: Vclassical(:) !< We use it to store the potential of the classical charges

    ! Ions
    FLOAT,             pointer :: vpsl(:)       !< the local part of the pseudopotentials
                                                !< plus the potential from static electric fields
    FLOAT,             pointer :: Imvpsl(:)     !< cmplxscl: imaginary part of vpsil          
    FLOAT,             pointer :: vpsl_lead(:, :) !< (np, NLEADS) the local part of the leads
    type(projector_t), pointer :: proj(:)       !< non-local projectors
    logical                    :: non_local
    integer                    :: natoms

    ! External e-m fields
    integer                :: no_lasers            !< number of laser pulses used
    type(laser_t), pointer :: lasers(:)            !< lasers stuff
    FLOAT,         pointer :: E_field(:)           !< static electric field
    FLOAT, pointer         :: v_static(:)          !< static scalar potential
    FLOAT, pointer         :: B_field(:)           !< static magnetic field
    FLOAT, pointer         :: A_static(:,:)        !< static vector potential
    type(gauge_field_t)    :: gfield               !< the time-dependent gauge field
    integer                :: reltype              !< type of relativistic correction to use

    !> The possible kick
    type(kick_t) :: kick

    !> The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    !! *effective* electrons in a quantum dot. It affects the spin Zeeman term.)
    FLOAT :: gyromagnetic_ratio

    !> SO prefactor (1.0 = normal SO, 0.0 = no SO)
    FLOAT :: so_strength
    
    !> the ion-ion energy and force
    FLOAT          :: eii
    FLOAT, pointer :: fii(:, :)
    
    real(4), pointer :: local_potential(:,:)
    real(4), pointer :: Imlocal_potential(:,:) !cmplxscl
    logical          :: local_potential_precalculated

    logical          :: ignore_external_ions
    logical                  :: have_density
    type(poisson_t), pointer :: poisson_solver
  end type epot_t

contains

  ! ---------------------------------------------------------
  subroutine epot_init( ep, gr, geo, ispin, nik, cmplxscl)
    type(epot_t),       intent(out)   :: ep
    type(grid_t),       intent(in)    :: gr
    type(geometry_t),   intent(inout) :: geo
    integer,            intent(in)    :: ispin
    integer,            intent(in)    :: nik
    logical,            intent(in)    :: cmplxscl


    integer :: ispec, ip, idir, ia, gauge_2d
    type(block_t) :: blk
    FLOAT, allocatable :: grx(:)
    integer :: filter

    PUSH_SUB(epot_init)

    !%Variable FilterPotentials
    !%Type integer
    !%Default filter_none
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
    call parse_integer(datasets_check('FilterPotentials'), PS_FILTER_NONE, filter)
    if(.not.varinfo_valid_option('FilterPotentials', filter)) call input_error('FilterPotentials')
    call messages_print_var_option(stdout, "FilterPotentials", filter)

    if(filter == PS_FILTER_TS) call spline_filter_mask_init()
    do ispec = 1, geo%nspecies
      call species_pot_init(geo%species(ispec), mesh_gcutoff(gr%mesh), filter)
    end do

    ! Local part of the pseudopotentials
    SAFE_ALLOCATE(ep%vpsl(1:gr%mesh%np))
    ep%vpsl(1:gr%mesh%np) = M_ZERO

    nullify(ep%Imvpsl)

    if(cmplxscl) then
      SAFE_ALLOCATE(ep%Imvpsl(1:gr%mesh%np))
      ep%Imvpsl(1:gr%mesh%np) = M_ZERO
    end if
    

    nullify(ep%vpsl_lead)
    if (gr%ob_grid%open_boundaries) then
      SAFE_ALLOCATE(ep%vpsl_lead(1:maxval(gr%ob_grid%lead(1:NLEADS)%mesh%np), 1:NLEADS))
      ep%vpsl_lead(1:maxval(gr%ob_grid%lead(1:NLEADS)%mesh%np), 1:NLEADS) = M_ZERO
    end if
    
    ep%classical_pot = 0
    nullify(ep%Vclassical)
    if(geo%ncatoms > 0) then

      !%Variable ClassicalPotential
      !%Type integer
      !%Default 0
      !%Section Hamiltonian
      !%Description
      !% Whether and how to add to the external potential the potential generated by 
      !% the classical charges read from the PDB input (see <tt>PBDCoordinates</tt>).
      !%Option no 0
      !%  No classical charges.
      !%Option point_charges 1
      !%  Classical charges are treated as point charges.
      !%Option gaussian_smeared 2
      !%  Classical charges are treated as Gaussian distributions. 
      !%  Smearing widths are hard-coded by species (experimental).
      !%End
      call parse_integer(datasets_check('ClassicalPotential'), 0, ep%classical_pot)
      if(ep%classical_pot .eq. CLASSICAL_GAUSSIAN) then
        call messages_experimental("Gaussian smeared classical charges")
        ! This method probably works but definitely needs to be made user-friendly:
        ! i.e. telling the user what widths are used and letting them be set somehow.
      endif

      if(ep%classical_pot > 0) then
        message(1) = 'Info: generating classical external potential.'
        call messages_info(1)

        SAFE_ALLOCATE(ep%Vclassical(1:gr%mesh%np))
        call epot_generate_classical(ep, gr%mesh, geo)
      end if
    end if

    ! lasers
    call laser_init(ep%no_lasers, ep%lasers, gr%mesh)

    select case(ispin)
    case(UNPOLARIZED)
      call kick_init(ep%kick, 1, gr%mesh%sb%dim)
    case(SPIN_POLARIZED)
      call kick_init(ep%kick, 2, gr%mesh%sb%dim)
    case(SPINORS)
      call kick_init(ep%kick, 4, gr%mesh%sb%dim)
    end select

    ! No more "UserDefinedTDPotential" from this version on.
    call messages_obsolete_variable('UserDefinedTDPotential', 'TDExternalFields')

    !%Variable StaticElectricField
    !%Type block
    !%Default 0
    !%Section Hamiltonian
    !%Description
    !% A static constant electric field may be added to the usual Hamiltonian,
    !% by setting the block <tt>StaticElectricField</tt>.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the electric field vector.
    !% It can be applied in a periodic direction of a large supercell via
    !% the single-point Berry phase.
    !%End
    nullify(ep%E_field, ep%v_static)
    if(parse_block(datasets_check('StaticElectricField'), blk)==0) then
      SAFE_ALLOCATE(ep%E_field(1:gr%sb%dim))
      do idir = 1, gr%sb%dim
        call parse_block_float(blk, 0, idir - 1, ep%E_field(idir), units_inp%energy / units_inp%length)

        if(idir <= gr%sb%periodic_dim .and. abs(ep%E_field(idir)) > M_EPSILON) then
          message(1) = "Applying StaticElectricField in a periodic direction is only accurate for large supercells."
          if(nik == 1) then
            call messages_warning(1)
          else
            message(2) = "Single-point Berry phase is not appropriate when k-point sampling is needed."
            call messages_warning(2)
          endif
        endif
      end do
      call parse_block_end(blk)

      if(gr%sb%periodic_dim < gr%sb%dim) then
        ! Compute the scalar potential
        !
        ! Note that the -1 sign is missing. This is because we
        ! consider the electrons with +1 charge. The electric field
        ! however retains the sign because we also consider protons to
        ! have +1 charge when calculating the force.
        SAFE_ALLOCATE(ep%v_static(1:gr%mesh%np))
        forall(ip = 1:gr%mesh%np)
          ep%v_static(ip) = sum(gr%mesh%x(ip, gr%sb%periodic_dim + 1:gr%sb%dim) * ep%E_field(gr%sb%periodic_dim + 1:gr%sb%dim))
        end forall
      endif
    endif


    !%Variable StaticMagneticField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% A static constant magnetic field may be added to the usual Hamiltonian,
    !% by setting the block <tt>StaticMagneticField</tt>.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the magnetic field vector. Note that
    !% if you are running the code in 1D mode, this will not work, and if you
    !% are running the code in 2D mode the magnetic field will have to be in
    !% the <i>z</i>-direction, so that the first two columns should be zero.
    !%
    !% The magnetic field should always be entered in atomic units, regardless
    !% of the <tt>Units</tt> variable. Note that we use the "Gaussian" system
    !% meaning 1 au[B] = 1.7152553 * 10^7 gauss, which corresponds to
    !% 1.7152553 * 10^3 Tesla.
    !%End
    nullify(ep%B_field, ep%A_static)
    if(parse_block(datasets_check('StaticMagneticField'), blk) == 0) then

      !%Variable StaticMagneticField2DGauge
      !%Type integer
      !%Default 0
      !%Section Hamiltonian
      !%Description
      !% The gauge of the static vector potential A when a magnetic field B = (0,0,B_z) is applied onto a 2D-system.
      !%Option linear_xy 0
      !% Linear gauge with A = ((1/2)/P_c)*(-y,x)*B_z. This is the default.
      !%Option linear_y 1
      !% Linear gauge with A = (1/P_c)*(-y,0)*B_z
      !%End
      call parse_integer(datasets_check('StaticMagneticField2DGauge'), 0, gauge_2d)
      if(.not.varinfo_valid_option('StaticMagneticField2DGauge', gauge_2d)) &
        call input_error('StaticMagneticField2DGauge')


      SAFE_ALLOCATE(ep%B_field(1:3))
      do idir = 1, 3
        call parse_block_float(blk, 0, idir - 1, ep%B_field(idir))
      end do
      select case(gr%sb%dim)
      case(1)
        call input_error('StaticMagneticField')
      case(2)
        if(ep%B_field(1)**2 + ep%B_field(2)**2 > M_ZERO) call input_error('StaticMagneticField')
      end select
      call parse_block_end(blk)

      if(gr%sb%dim > 3) call messages_not_implemented('Magnetic field for dim > 3')

      ! Compute the vector potential
      SAFE_ALLOCATE(ep%A_static(1:gr%mesh%np, 1:gr%sb%dim))
      SAFE_ALLOCATE(grx(1:gr%sb%dim))

      select case(gr%sb%dim)
      case(2)
        select case(gauge_2d)
        case(0) ! linear_xy
          do ip = 1, gr%mesh%np
            grx(1:gr%sb%dim) = gr%mesh%x(ip, 1:gr%sb%dim)
            ep%A_static(ip, :) = M_HALF/P_C*(/grx(2), -grx(1)/) * ep%B_field(3)
          end do
        case(1) ! linear y
          do ip = 1, gr%mesh%np
            grx(1:gr%sb%dim) = gr%mesh%x(ip, 1:gr%sb%dim)
            ep%A_static(ip, :) = M_ONE/P_C*(/grx(2), M_ZERO/) * ep%B_field(3)
          end do
        end select
      case(3)
        do ip = 1, gr%mesh%np
          grx(1:gr%sb%dim) = gr%mesh%x(ip, 1:gr%sb%dim)
          ep%A_static(ip, :) = M_HALF/P_C*(/grx(2) * ep%B_field(3) - grx(3) * ep%B_field(2), &
                               grx(3) * ep%B_field(1) - grx(1) * ep%B_field(3), &
                               grx(1) * ep%B_field(2) - grx(2) * ep%B_field(1)/)
        end do
      end select

      SAFE_DEALLOCATE_A(grx)

    end if
    
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
    call parse_float(datasets_check('GyromagneticRatio'), P_g, ep%gyromagnetic_ratio)

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
    call parse_integer(datasets_check('RelativisticCorrection'), NOREL, ep%reltype)
    if(.not.varinfo_valid_option('RelativisticCorrection', ep%reltype)) call input_error('RelativisticCorrection')
    if (ispin /= SPINORS .and. ep%reltype == SPIN_ORBIT) then
      message(1) = "The spin-orbit term can only be applied when using spinors."
      call messages_fatal(1)
    end if

    call messages_print_var_option(stdout, "RelativisticCorrection", ep%reltype)

    !%Variable SOStrength
    !%Type float
    !%Default 1
    !%Section Hamiltonian
    !%Description
    !% Tuning of the spin-orbit coupling strength: setting this value to zero turns off spin-orbit terms in
    !% the Hamiltonian, and setting it to one corresponds to full spin-orbit.
    !%End
    if (ep%reltype == SPIN_ORBIT) then
      call parse_float(datasets_check('SOStrength'), M_ONE, ep%so_strength)
    else
      ep%so_strength = M_ONE
    end if

    !%Variable IgnoreExternalIons
    !%Type logical
    !%Default no
    !%Section Hamiltonian
    !%Description
    !% If this variable is set to "yes", then the ions that are outside the simulation box do not feel any
    !% external force (and therefore progress at constant velocity), and do not originate any force on other
    !% ions, or any potential on the electronic system.
    !%
    !% This feature is only available for finite systems; if the system is periodic in any dimension, 
    !% this variable cannot be set to "yes".
    !%End
    call parse_logical(datasets_check('IgnoreExternalIons'), .false., ep%ignore_external_ions)
    if(ep%ignore_external_ions) then
      if(gr%sb%periodic_dim > 0) call input_error('IgnoreExternalIons')
    end if

    SAFE_ALLOCATE(ep%proj(1:geo%natoms))
    do ia = 1, geo%natoms
      call projector_null(ep%proj(ia))
    end do

    ep%natoms = geo%natoms
    ep%non_local = .false.

    ep%eii = M_ZERO
    SAFE_ALLOCATE(ep%fii(1:gr%sb%dim, 1:geo%natoms))
    ep%fii = M_ZERO

    call gauge_field_nullify(ep%gfield)

    nullify(ep%local_potential)
    nullify(ep%Imlocal_potential)
    ep%local_potential_precalculated = .false.
    

    ep%have_density = .false.
    do ia = 1, geo%natoms
      if(local_potential_has_density(gr%mesh%sb, geo%atom(ia))) then
        ep%have_density = .true.
        exit
      end if
    end do

    if(ep%have_density) then
      ep%poisson_solver => psolver

      if (poisson_get_solver(ep%poisson_solver) == POISSON_SETE) then
        SAFE_ALLOCATE(rho_nuc(1:gr%mesh%np))
        rho_nuc(1:gr%mesh%np) = M_ZERO
        SAFE_ALLOCATE(v_es(1:gr%mesh%np))
        v_es(1:gr%mesh%np) = M_ZERO
      end if
      
    else
      nullify(ep%poisson_solver)
    end if

    POP_SUB(epot_init)
  end subroutine epot_init


  ! ---------------------------------------------------------
  subroutine epot_end(ep, geo)
    type(epot_t),      intent(inout) :: ep
    type(geometry_t),  intent(inout) :: geo

    integer :: iproj

    PUSH_SUB(epot_end)

    if(ep%have_density) then

      if (poisson_get_solver(ep%poisson_solver) == POISSON_SETE) then 
        SAFE_DEALLOCATE_A(rho_nuc)
        SAFE_DEALLOCATE_A(v_es)
        SAFE_DEALLOCATE_A(v_es3)
      end if

      nullify(ep%poisson_solver)
    end if

    SAFE_DEALLOCATE_P(ep%local_potential)
    SAFE_DEALLOCATE_P(ep%Imlocal_potential)!cmplxscl 
    SAFE_DEALLOCATE_P(ep%fii)

    SAFE_DEALLOCATE_P(ep%vpsl)
    SAFE_DEALLOCATE_P(ep%vpsl_lead)
    SAFE_DEALLOCATE_P(ep%Imvpsl)!cmplxscl

    if(ep%classical_pot > 0) then
      ep%classical_pot = 0
      ! sanity check
      ASSERT(associated(ep%Vclassical)) 
      SAFE_DEALLOCATE_P(ep%Vclassical)         ! and clean up
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

    POP_SUB(epot_end)

  end subroutine epot_end

  ! ---------------------------------------------------------
  subroutine epot_generate(ep, gr, geo, st, cmplxscl)
    type(epot_t),             intent(inout) :: ep
    type(grid_t),     target, intent(in)    :: gr
    type(geometry_t), target, intent(in)    :: geo
    type(states_t),           intent(inout) :: st
    logical,                  intent(in)    :: cmplxscl  

    integer :: ia, ip
    type(atom_t),      pointer :: atm
    type(mesh_t),      pointer :: mesh
    type(simul_box_t), pointer :: sb
    type(profile_t), save :: epot_generate_prof
    FLOAT,    allocatable :: density(:)
    FLOAT,    allocatable :: tmp(:)
    type(profile_t), save :: epot_reduce

    call profiling_in(epot_generate_prof, "EPOT_GENERATE")
    PUSH_SUB(epot_generate)

    sb   => gr%sb
    mesh => gr%mesh

    SAFE_ALLOCATE(density(1:mesh%np))
    density = M_ZERO

    ! Local part
    ep%vpsl = M_ZERO
    if(associated(ep%Imvpsl)) ep%Imvpsl = M_ZERO
    if(geo%nlcc) st%rho_core = M_ZERO

    do ia = geo%atoms_dist%start, geo%atoms_dist%end
      if(.not.simul_box_in_box(sb, geo, geo%atom(ia)%x) .and. ep%ignore_external_ions) cycle
      if(geo%nlcc) then
        if(cmplxscl) then
        call epot_local_potential(ep, gr%der, gr%dgrid, geo, ia, ep%vpsl, ep%Imvpsl, &
          rho_core = st%rho_core, density = density)
        else
        call epot_local_potential(ep, gr%der, gr%dgrid, geo, ia, ep%vpsl, &
          rho_core = st%rho_core, density = density)
        end if
      else
        if(cmplxscl) then
          call epot_local_potential(ep, gr%der, gr%dgrid, geo, ia, ep%vpsl, ep%Imvpsl, density = density)
        else
          call epot_local_potential(ep, gr%der, gr%dgrid, geo, ia, ep%vpsl, density = density)
        end if
      end if
    end do

    ! reduce over atoms if required
    if(geo%atoms_dist%parallel) then
      call profiling_in(epot_reduce, "EPOT_REDUCE")

      call comm_allreduce(geo%atoms_dist%mpi_grp%comm, ep%vpsl, dim = gr%mesh%np)
      if(associated(st%rho_core)) call comm_allreduce(geo%atoms_dist%mpi_grp%comm, st%rho_core, dim = gr%mesh%np)
      if(ep%have_density) call comm_allreduce(geo%atoms_dist%mpi_grp%comm, density, dim = gr%mesh%np)

      call profiling_out(epot_reduce)
    end if

    if(ep%have_density) then
      ! now we solve the poisson equation with the density of all nodes
      SAFE_ALLOCATE(tmp(1:gr%mesh%np_part))

      if(poisson_solver_is_iterative(ep%poisson_solver)) tmp(1:mesh%np) = M_ZERO

      call dpoisson_solve(ep%poisson_solver, tmp, density)

      forall(ip = 1:mesh%np) ep%vpsl(ip) = ep%vpsl(ip) + tmp(ip)

      SAFE_DEALLOCATE_A(tmp)
    end if
    SAFE_DEALLOCATE_A(density)

    ! we assume that we need to recalculate the ion-ion energy
    call ion_interaction_calculate(geo, sb, gr, ep, ep%eii, ep%fii)

    ! the pseudopotential part.
    do ia = 1, geo%natoms
      atm => geo%atom(ia)
      if(.not. species_is_ps(atm%spec)) cycle
      if(.not.simul_box_in_box(sb, geo, geo%atom(ia)%x) .and. ep%ignore_external_ions) cycle
      call projector_end(ep%proj(ia))
      call projector_init(ep%proj(ia), gr%mesh, atm, st%d%dim, ep%reltype)
      call projector_build(ep%proj(ia), gr, atm, ep%so_strength)
      if(.not. projector_is(ep%proj(ia), M_NONE)) ep%non_local = .true.
    end do

    ! add static electric fields
    if (ep%classical_pot > 0)   ep%vpsl(1:mesh%np) = ep%vpsl(1:mesh%np) + ep%Vclassical(1:mesh%np)
    if (associated(ep%e_field) .and. sb%periodic_dim < sb%dim) ep%vpsl(1:mesh%np) = ep%vpsl(1:mesh%np) + ep%v_static(1:mesh%np)

    POP_SUB(epot_generate)
    call profiling_out(epot_generate_prof)
  end subroutine epot_generate

  ! ---------------------------------------------------------

  logical pure function local_potential_has_density(sb, atom) result(has_density)
    type(simul_box_t),        intent(in)    :: sb
    type(atom_t),             intent(in)    :: atom
    
    
    has_density = &
      species_has_density(atom%spec) .or. (species_is_ps(atom%spec) .and. simul_box_is_periodic(sb)) &
      .or. (species_is_ps(atom%spec) .and. simul_box_complex_boundaries(sb))

  end function local_potential_has_density
  
  ! ---------------------------------------------------------
  subroutine epot_local_potential(ep, der, dgrid, geo, iatom, vpsl, Imvpsl, rho_core, density)
    type(epot_t),             intent(inout) :: ep
    type(derivatives_t),      intent(in)    :: der
    type(double_grid_t),      intent(in)    :: dgrid
    type(geometry_t),         intent(in)    :: geo
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,          optional, intent(inout) :: Imvpsl(:)
    FLOAT,          optional, pointer       :: rho_core(:)
    FLOAT,          optional, intent(inout) :: density(:) !< If present, the ionic density will be added here.

    integer :: ip
    FLOAT :: radius
    FLOAT, allocatable :: vl(:), Imvl(:), rho(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    integer :: counter, conversion(3), nx_half, ny_half, nz_half !ROA
    integer :: nx, ny, nz !ROA

    PUSH_SUB(epot_local_potential)
    call profiling_in(prof, "EPOT_LOCAL")

    if(ep%local_potential_precalculated) then

      forall(ip = 1:der%mesh%np) vpsl(ip) = vpsl(ip) + ep%local_potential(ip, iatom)
      if(present(Imvpsl)) then !cmplxscl
        forall(ip = 1:der%mesh%np) Imvpsl(ip) = Imvpsl(ip) + ep%Imlocal_potential(ip, iatom)
      end if
    else

      !Local potential, we can get it by solving the Poisson equation
      !(for all-electron species or pseudopotentials in periodic
      !systems) or by applying it directly to the grid

      if(local_potential_has_density(der%mesh%sb, geo%atom(iatom))) then
        SAFE_ALLOCATE(rho(1:der%mesh%np))

        call species_get_density(geo%atom(iatom)%spec, geo%atom(iatom)%x, der%mesh, rho)

        if(present(density)) then
          forall(ip = 1:der%mesh%np) density(ip) = density(ip) + rho(ip)
        else

          SAFE_ALLOCATE(vl(1:der%mesh%np))
          
          if(poisson_solver_is_iterative(ep%poisson_solver)) then
            ! vl has to be initialized before entering routine
            ! and our best guess for the potential is zero
            vl(1:der%mesh%np) = M_ZERO
          end if

          call dpoisson_solve(ep%poisson_solver, vl, rho, all_nodes = .false.)
          
        end if

        count_atoms=iatom

        if (poisson_get_solver(ep%poisson_solver) == POISSON_SETE) then  !SEC
          write(68,*) "Calling rhonuc iatom", iatom
          rho_nuc(1:der%mesh%np) = rho_nuc(1:der%mesh%np) + rho(1:der%mesh%np)
          if (iatom==geo%natoms.and.calc_gate_energy==0) then
             write(68,*) "Entering the zone"
             nx = der%mesh%idx%nr(2,1) - der%mesh%idx%nr(1,1) + 1 - 2*der%mesh%idx%enlarge(1)
             ny = der%mesh%idx%nr(2,2) - der%mesh%idx%nr(1,2) + 1 - 2*der%mesh%idx%enlarge(2)
             nz = der%mesh%idx%nr(2,3) - der%mesh%idx%nr(1,3) + 1 - 2*der%mesh%idx%enlarge(3)
             nx_half = (der%mesh%idx%nr(2,1) - der%mesh%idx%nr(1,1) - 2*der%mesh%idx%enlarge(1))/2 + 1
             ny_half = (der%mesh%idx%nr(2,2) - der%mesh%idx%nr(1,2) - 2*der%mesh%idx%enlarge(2))/2 + 1
             nz_half = (der%mesh%idx%nr(2,3) - der%mesh%idx%nr(1,3) - 2*der%mesh%idx%enlarge(3))/2 + 1
             do counter = 1, der%mesh%np
               call  index_to_coords(der%mesh%idx,der%mesh%sb%dim,counter,conversion)
               conversion(1) = conversion(1) + nx_half
               conversion(2) = conversion(2) + ny_half
               conversion(3) = conversion(3) + nz_half
               v_es(counter)=v_es3(conversion(1),conversion(2),conversion(3))
            enddo
            es_energy = dmf_dotp(der%mesh, rho_nuc, v_es)
            write(68,*) "Calculating ion energy due to bias:", es_energy*CNST(2.0*13.60569193)
          end if
         end if
         SAFE_DEALLOCATE_A(rho)

      else

        SAFE_ALLOCATE(vl(1:der%mesh%np))
        if(present(Imvpsl)) then !cmplxscl
          SAFE_ALLOCATE(Imvl(1:der%mesh%np))
          call species_get_local(geo%atom(iatom)%spec, der%mesh, geo%atom(iatom)%x(1:der%mesh%sb%dim), vl, Imvl)
        else
          call species_get_local(geo%atom(iatom)%spec, der%mesh, geo%atom(iatom)%x(1:der%mesh%sb%dim), vl)                      
        end if
      end if

      if(allocated(vl)) then
        forall(ip = 1:der%mesh%np) vpsl(ip) = vpsl(ip) + vl(ip)
        SAFE_DEALLOCATE_A(vl)
        if(allocated(Imvl) ) then !cmplxscl
          forall(ip = 1:der%mesh%np) Imvpsl(ip) = Imvpsl(ip) + Imvl(ip)
          SAFE_DEALLOCATE_A(Imvl)
        end if
      end if

      !the localized part
      if(species_is_ps(geo%atom(iatom)%spec)) then

        radius = double_grid_get_rmax(dgrid, geo%atom(iatom)%spec, der%mesh) + der%mesh%spacing(1)

        call submesh_init_sphere(sphere, der%mesh%sb, der%mesh, geo%atom(iatom)%x, radius)
        SAFE_ALLOCATE(vl(1:sphere%np))

        call double_grid_apply_local(dgrid, geo%atom(iatom)%spec, der%mesh, sphere, geo%atom(iatom)%x, vl)
        vpsl(sphere%map(1:sphere%np)) = vpsl(sphere%map(1:sphere%np)) + vl(1:sphere%np)

        SAFE_DEALLOCATE_A(vl)
        call submesh_end(sphere)

      end if

    end if

    !Non-local core corrections
    if(present(rho_core) .and. &
      species_has_nlcc(geo%atom(iatom)%spec) .and. &
      species_is_ps(geo%atom(iatom)%spec)) then
      SAFE_ALLOCATE(rho(1:der%mesh%np))
      call species_get_nlcc(geo%atom(iatom)%spec, geo%atom(iatom)%x, der%mesh, rho)
      forall(ip = 1:der%mesh%np) rho_core(ip) = rho_core(ip) + rho(ip)
      SAFE_DEALLOCATE_A(rho)
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
          if(rr < r_small) rr = r_small
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
          if(abs(rr - rc) < r_small) rr = rc + sign(r_small, rr - rc)
          ep%Vclassical(ip) = ep%Vclassical(ip) - geo%catom(ia)%charge * (rr**4 - rc**4) / (rr**5 - rc**5)
        case default
          call input_error('ClassicalPotential')
        end select
      end do
    end do

    POP_SUB(epot_generate_classical)
  end subroutine epot_generate_classical

  ! ---------------------------------------------------------
  subroutine epot_precalc_local_potential(ep, gr, geo)
    type(epot_t),     intent(inout) :: ep
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo

    integer :: iatom
    FLOAT, allocatable :: tmp(:), Imtmp(:)
    logical  :: cmplxscl
    
    PUSH_SUB(epot_precalc_local_potential)
    cmplxscl =.false.
    if(associated(ep%Imvpsl)) cmplxscl = .true.
    
    if(.not. associated(ep%local_potential)) then
      SAFE_ALLOCATE(ep%local_potential(1:gr%mesh%np, 1:geo%natoms))
    end if

    if(.not. associated(ep%Imlocal_potential) .and. cmplxscl) then
      SAFE_ALLOCATE(ep%Imlocal_potential(1:gr%mesh%np, 1:geo%natoms))
      SAFE_ALLOCATE(Imtmp(1:gr%mesh%np))
    end if


    ep%local_potential_precalculated = .false.

    SAFE_ALLOCATE(tmp(1:gr%mesh%np))
    
    if(cmplxscl) then
      do iatom = 1, geo%natoms
        tmp(1:gr%mesh%np) = M_ZERO
        Imtmp(1:gr%mesh%np) = M_ZERO
        call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, tmp, Imtmp)!, time)
        ep%local_potential(1:gr%mesh%np, iatom) = tmp(1:gr%mesh%np) 
        ep%Imlocal_potential(1:gr%mesh%np, iatom) = Imtmp(1:gr%mesh%np)
      end do
    else
      do iatom = 1, geo%natoms
        tmp(1:gr%mesh%np) = M_ZERO
        call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, tmp)!, time)
        ep%local_potential(1:gr%mesh%np, iatom) = tmp(1:gr%mesh%np)
      end do
    end if
    ep%local_potential_precalculated = .true.

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(Imtmp)
    POP_SUB(epot_precalc_local_potential)
  end subroutine epot_precalc_local_potential

  ! ---------------------------------------------------------
  !> For details about this routine, see
  !! http://www.tddft.org/programs/octopus/wiki/index.php/Developers:Ion-Ion_interaction
  subroutine ion_interaction_calculate(geo, sb, gr, ep, energy, force)
    type(geometry_t), target, intent(in)    :: geo
    type(simul_box_t),        intent(in)    :: sb
    type(epot_t),             intent(inout) :: ep
    type(grid_t),             intent(in)    :: gr
    FLOAT,                    intent(out)   :: energy
    FLOAT,    dimension(:,:), intent(out)   :: force
    !
    FLOAT, parameter :: alpha = CNST(1.1313708)
    !
    FLOAT, dimension(sb%dim) :: r, f
    FLOAT :: rr, dd, zi, zj, epsilon, sigma
    integer :: iatom, jatom, natom, iindex, jindex
    type(species_t), pointer :: spci, spcj
    type(profile_t), save :: ion_ion_prof
    logical,  allocatable :: in_box(:)

    call profiling_in(ion_ion_prof, "ION_ION_INTERACTION")
    PUSH_SUB(ion_interaction_calculate)

    energy = M_ZERO
    force(1:sb%dim,1:geo%natoms) = M_ZERO

    if(simul_box_is_periodic(sb)) then
      ASSERT(geo%ncatoms==0)
      ! This depends on the area, but we should check if it is fully consistent.        
      spci => geo%atom(1)%spec
      if( species_type(spci).eq.SPEC_JELLI_SLAB ) then
        energy = energy +M_PI *species_zval(spci)**2 /( M_FOUR *sb%lsize(1) *sb%lsize(2) ) &
          & *( sb%lsize(3) - species_jthick(spci) /M_THREE ) 
      else
        call ion_interaction_periodic(geo, sb, energy, force)
      end if

    else if(simul_box_complex_boundaries(sb)) then
      ASSERT(geo%ncatoms==0)
      ! only interaction inside the cell
      write(68,'(a)') "Calling SETE interaction"
      write(6,'(a)') "Calling SETE interaction"
      !ep%eii= M_ZERO
      write(68,*) "energy, ep%eii before ion interaction sete" , energy, ep%eii
      call ion_interaction_sete(gr, sb, geo, ep) 
      write(*,*) "energy, ep%eii", energy, ep%eii
    else
      natom=geo%natoms+geo%ncatoms
      if(ep%ignore_external_ions) then
        SAFE_ALLOCATE(in_box(natom))
        do iatom = 1, geo%natoms
          in_box(iatom) = simul_box_in_box(sb, geo, geo%atom(iatom)%x)
        end do
        do iatom = 1, geo%ncatoms
          in_box(geo%natoms+iatom) = simul_box_in_box(sb, geo, geo%catom(iatom)%x)
        end do
      end if
      ! only interaction inside the cell
      do iatom = 1, geo%natoms
        if(ep%ignore_external_ions) then
          if(.not. in_box(iatom)) cycle
        end if
        spci=>geo%atom(iatom)%spec
        zi=species_zval(spci)
        select case(species_type(spci))
        case(SPEC_JELLI)
          energy=energy+(M_THREE/M_FIVE)*zi**2/species_jradius(spci)
          ! The part depending on the simulation sphere is neglected
        case(SPEC_JELLI_SLAB)
          energy=energy-M_PI*zi**2/(M_FOUR*sb%lsize(1)*sb%lsize(2)) &
            *species_jthick(spci)/M_THREE
          ! The part depending on the simulation box transverse dimension is neglected
        end select
        do jatom = iatom+1, geo%natoms
          if(ep%ignore_external_ions) then
            if(.not. in_box(jatom)) cycle
          end if
          spcj=>geo%atom(jatom)%spec
          r=geo%atom(iatom)%x(1:sb%dim)-geo%atom(jatom)%x(1:sb%dim)
          rr=sqrt(sum(r**2))
          iindex=species_index(spci)
          jindex=species_index(spcj)
          select case(geo%ionic_interaction_type(iindex,jindex))
          case(INTERACTION_COULOMB)
            zj=species_zval(spcj)
            !the force
            dd=zi*zj/rr
            f=(dd/rr**2)*r
            force(1:sb%dim,iatom)=force(1:sb%dim,iatom)+f
            force(1:sb%dim,jatom)=force(1:sb%dim,jatom)-f
            !energy
            energy=energy+dd
          case(INTERACTION_LJ)
            epsilon=geo%ionic_interaction_parameter(LJ_EPSILON, iindex, jindex)
            sigma  =geo%ionic_interaction_parameter(LJ_SIGMA,   iindex, jindex)
            dd = (sigma/rr)**6
            !the force
            f=(CNST(24.0)*epsilon*(dd/rr**2)*(CNST(2.0)*dd-M_ONE))*r
            force(1:sb%dim,iatom)=force(1:sb%dim,iatom)+f
            force(1:sb%dim,jatom)=force(1:sb%dim,jatom)-f
            !energy
            energy=energy+CNST(4.0)*epsilon*dd*(dd-M_ONE)
          end select
        end do !jatom
      end do !iatom
      do iatom = 1, geo%natoms
        if(ep%ignore_external_ions) then
          if(.not. in_box(iatom)) cycle
        end if
        do jatom = 1, geo%ncatoms
          if(ep%ignore_external_ions) then
            if(.not. in_box(geo%natoms+jatom)) cycle
          end if
          r=geo%atom(iatom)%x(1:sb%dim)-geo%catom(jatom)%x(1:sb%dim)
          rr=sqrt(sum(r**2))
          !INTERACTION_COULOMB
          zi=species_zval(geo%atom(iatom)%spec)
          zj=geo%catom(jatom)%charge
          !the force
          dd=zi*zj/rr
          force(1:sb%dim,iatom)=force(1:sb%dim,iatom)+(dd/rr*2)*r
          !energy
          energy=energy+dd
        end do !jatom
      end do !iatom
    end if
    SAFE_DEALLOCATE_A(in_box)
    call profiling_out(ion_ion_prof)
    POP_SUB(ion_interaction_calculate)
    return
  end subroutine ion_interaction_calculate

  ! ---------------------------------------------------------
  subroutine ion_interaction_periodic(geo, sb, energy, force)
    type(geometry_t),  target, intent(in)    :: geo
    type(simul_box_t),         intent(in)    :: sb
    FLOAT,                     intent(out)   :: energy
    FLOAT,                     intent(out)   :: force(:, :) !< sb%dim, geo%natoms

    type(species_t), pointer :: spec
    FLOAT :: rr, xi(1:MAX_DIM), zi, zj, ereal, efourier, eself, erfc, rcut
    integer :: iatom, jatom, icopy
    type(periodic_copy_t) :: pc
    integer :: ix, iy, iz, isph, ss, idim
    FLOAT   :: gg(1:MAX_DIM), gg2, gx
    FLOAT   :: factor, charge
    CMPLX   :: sumatoms, tmp(1:MAX_DIM)
    FLOAT, parameter :: alpha = CNST(1.1313708)
    CMPLX, allocatable :: phase(:)

    PUSH_SUB(ion_interaction_periodic)

    ! Check http://www.duke.edu/~kz10/file/ewald.pdf for the equations
    ! implemented here.

    if(any(geo%ionic_interaction_type /= INTERACTION_COULOMB)) then
      message(1) = "Cannot calculate non-Coulombic interaction for periodic systems."
      call messages_fatal(1)
    end if

    ereal = M_ZERO

    force(1:sb%dim, 1:geo%natoms) = M_ZERO

    ! if the system is periodic we have to add the energy of the
    ! interaction with the copies
    
    rcut = CNST(6.0)/alpha

    ! the short-range part is calculated directly
    do iatom = 1, geo%natoms
      spec => geo%atom(iatom)%spec
      if (.not. species_is_ps(spec)) cycle
      zi = species_zval(geo%atom(iatom)%spec)

      call periodic_copy_init(pc, sb, geo%atom(iatom)%x, rcut)
      
      do icopy = 1, periodic_copy_num(pc)
        xi(1:sb%dim) = periodic_copy_position(pc, sb, icopy)
        
        do jatom = 1, geo%natoms
          zj = species_zval(geo%atom(jatom)%spec)
          rr = sqrt( sum( (xi(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))**2 ) )
          
          if(rr < CNST(1e-5)) cycle
          
          erfc = M_ONE - loct_erf(alpha*rr)

          ! energy
          ereal = ereal + M_HALF*zj*zi*erfc/rr
          
          ! force
          force(1:sb%dim, jatom) = force(1:sb%dim, jatom) + &
            M_HALF*zj*zi*(xi(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))*&
            (erfc/rr + M_TWO*alpha/sqrt(M_PI)*exp(-(alpha*rr)**2))/rr**2
        end do
        
      end do
      
      call periodic_copy_end(pc)
    end do

    ! self-interaction
    eself = M_ZERO
    charge = M_ZERO
    do iatom = 1, geo%natoms
      zi = species_zval(geo%atom(iatom)%spec)
      charge = charge + zi
      eself = eself - alpha/sqrt(M_PI)*zi**2
    end do

    ! And the long-range part, using an Ewald sum
    SAFE_ALLOCATE(phase(1:geo%natoms))

    ! First the G = 0 term (charge was calculated previously)
    efourier = -M_PI*charge**2/(M_TWO*alpha**2*sb%rcell_volume)

    ! get a converged value for the cutoff in g
    rcut = sum(sb%klattice(1:sb%dim, 1))**2
    do idim = 2, sb%dim
      rcut = min(rcut, sum(sb%klattice(1:sb%dim, idim))**2)
    end do

    isph = ceiling(CNST(6.0)*alpha/rcut)

    do ix = -isph, isph
      do iy = -isph, isph
        do iz = -isph, isph
          
          ss = ix**2 + iy**2 + iz**2
          
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:sb%dim) = ix*sb%klattice(1:sb%dim, 1) + iy*sb%klattice(1:sb%dim, 2) + iz*sb%klattice(1:sb%dim, 3)
          gg2 = sum(gg(1:sb%dim)**2)
          
          ! g=0 must be removed from the sum
          if(gg2 == M_ZERO) cycle

          gx = -CNST(0.25)*gg2/alpha**2

          if(gx < CNST(-36.0)) cycle

          factor = M_TWO*M_PI/sb%rcell_volume*exp(gx)/gg2

          if(factor < epsilon(factor)) cycle

          sumatoms = M_Z0
          do iatom = 1, geo%natoms
            zi = species_zval(geo%atom(iatom)%spec)
            xi(1:sb%dim) = geo%atom(iatom)%x(1:sb%dim)
            gx = sum(gg(1:sb%dim)*xi(1:sb%dim))
            phase(iatom) = zi*TOCMPLX(cos(gx), sin(gx))
            sumatoms = sumatoms + phase(iatom)
          end do
          
          efourier = efourier + factor*sumatoms*conjg(sumatoms)
          
          do iatom = 1, geo%natoms
            tmp(1:sb%dim) = M_ZI*gg(1:sb%dim)*phase(iatom)
            force(1:sb%dim, iatom) = force(1:sb%dim, iatom) &
              - factor*(conjg(tmp(1:sb%dim))*sumatoms + tmp(1:sb%dim)*conjg(sumatoms))
          end do
          
        end do
      end do
    end do
    
    energy = ereal + efourier + eself
    
    SAFE_DEALLOCATE_A(phase)

    POP_SUB(ion_interaction_periodic)
  end subroutine ion_interaction_periodic


  ! ------------------------------------------------------------
  subroutine ion_interaction_sete(gr, sb, geo, ep)
    type(grid_t),      intent(in)    :: gr
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(in)    :: geo
    type(epot_t),      intent(inout) :: ep
  
    integer                              :: iatom, jatom
    FLOAT, allocatable                   :: rho1(:), v2(:), rho2(:)
    FLOAT                                :: temp
    FLOAT                                :: time1
    FLOAT :: rr, dd, zi, zj
    FLOAT :: sicn, chargeion, sqrtalphapi,sicn2

    PUSH_SUB(ion_interaction_sete)

    write(68,*) "ep%eii values top", ep%eii
    ep%eii=es_energy*M_HALF
    write(68,*) "ep%eii values top", ep%eii

!    SAFE _ALLOCATE(rho2(1:gr%mesh%np))
!    SAFE _ALLOCATE(rho3(1:gr%mesh%np))
!    v2(1:gr%mesh%np)= M_ZERO
!    v3(1:gr%mesh%np)= M_ZERO
!    rho2(1:gr%mesh%np)= M_ZERO
!    rho3(1:gr%mesh%np)= M_ZERO
    sicn=M_ZERO
    sicn2=M_ZERO
    !This value represents sqrt(alpha/pi)
    !On species/ps.F90, the value of sigma_erf = 0.625
    !alpha=(1/2)*(1/erf)**2
    !It is the same as alpha = CNST(1.1313708)/sqrt(M_PI)
    sqrtalphapi=sqrt(M_HALF*(CNST(1.6)**2)/M_PI)
!    do iatom=1,geo%natoms
!      call epot_local_potential(ep, gr%der, gr%dgrid, ep%poisson_solver, geo, iatom, v2, time1)
!      call species_get_density(geo%atom(iatom)%spec, geo%atom(iatom)%x, gr%mesh, rho2)
!      v3(1:gr%mesh%np)=v3(1:gr%mesh%np)+v2(1:gr%mesh%np)
!      rho3(1:gr%mesh%np)=rho3(1:gr%mesh%np)+rho2(1:gr%mesh%np)
!      !Substract this from total energy
!      write(68,*) "Getting self-interaction correction of", iatom
!      chargeion=species_zval(geo%atom(iatom)%spec)
!      sicn=sicn+sqrtalphapi*chargeion*chargeion
!      write(68,*) "chargeion, value", chargeion,sqrtalphapi*chargeion*chargeion, sicn
!    end do
!    temp=dmf_dotp(gr%mesh, rho3, v3) 
!    temp=M_HALF*temp
!    write(68,*) "temp, sicn", temp, sicn
!    ep%eii=ep%eii+temp-sicn
!    SAFE_DEALLOCATE _A(rho2)
!    SAFE_DEALLOCATE _A(v2)
!    SAFE_DEALLOCATE _A(rho3)
!    SAFE_DEALLOCATE _A(v3)

    do jatom = geo%natoms,2, -1
      SAFE_ALLOCATE(v2(1:gr%mesh%np))
      SAFE_ALLOCATE(rho2(1:gr%mesh%np))
      v2(1:gr%mesh%np)= M_ZERO
      count_atoms=jatom
      calc_gate_energy=1
      call epot_local_potential(ep, gr%der, gr%dgrid, geo, jatom, v2)
      write(68,*) "time1,", time1

      do iatom = jatom-1,1,-1
        SAFE_ALLOCATE(rho1(1:gr%mesh%np))
        call species_get_density(geo%atom(iatom)%spec, geo%atom(iatom)%x, gr%mesh, rho1)
        temp=dmf_dotp(gr%mesh, rho1, v2) 
        ep%eii = ep%eii + temp 
        SAFE_DEALLOCATE_A(rho1)
        write(68,*) "ep%eii values", iatom, jatom, ep%eii, temp
      enddo
      SAFE_DEALLOCATE_A(rho2)
      SAFE_DEALLOCATE_A(v2)
    enddo

    !Testing self-interaction approximation
    temp=0.0
    do iatom=1,geo%natoms
      SAFE_ALLOCATE(v2(1:gr%mesh%np))
      SAFE_ALLOCATE(rho2(1:gr%mesh%np))
      v2(1:gr%mesh%np)= M_ZERO
      rho2(1:gr%mesh%np)= M_ZERO
      call epot_local_potential(ep, gr%der, gr%dgrid, geo, iatom, v2)
      call species_get_density(geo%atom(iatom)%spec, geo%atom(iatom)%x, gr%mesh, rho2)
      temp=dmf_dotp(gr%mesh, rho2, v2) 
      sicn2=sicn2+M_HALF*temp
      chargeion=species_zval(geo%atom(iatom)%spec)
      sicn=sicn+sqrtalphapi*chargeion*chargeion
      write(68,*) "Getting sic of:", iatom,M_HALF*temp, sicn2, sqrtalphapi*chargeion*chargeion, sicn
      SAFE_DEALLOCATE_A(v2)
      SAFE_DEALLOCATE_A(rho2)
     enddo
     !Add off-diagonal, self-interaction and its correction
     ep%eii=ep%eii+sicn2-sicn

    write(68,*) "ep%eii values very bottom", ep%eii * CNST(2.0)*CNST(13.60569193)

    do iatom = 1, geo%natoms
      zi = species_zval(geo%atom(iatom)%spec)

      do jatom = 1, geo%natoms
        if(iatom == jatom) cycle
        zj = species_zval(geo%atom(jatom)%spec)
        rr = sqrt(sum((geo%atom(iatom)%x - geo%atom(jatom)%x)**2))
        !the force
        dd = zi * zj / rr**3
        ep%fii(1:sb%dim, iatom) = ep%fii(1:sb%dim, iatom) + &
          dd * (geo%atom(iatom)%x(1:sb%dim) - geo%atom(jatom)%x(1:sb%dim))
        ep%fii(1:sb%dim, iatom) = M_ZERO         
      enddo
      write(68,*) "SETE Force:", iatom,ep%fii(1:sb%dim,iatom)
      !calc_gate_energy=0
    enddo 

    POP_SUB(ion_interaction_sete)
   end subroutine ion_interaction_sete
      
end module epot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
