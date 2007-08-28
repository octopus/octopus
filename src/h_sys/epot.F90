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
  use double_grid_m
  use functions_m
  use global_m
  use grid_m
  use io_m
  use lib_oct_parser_m
  use lib_oct_gsl_spline_m
  use magnetic_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use multicomm_m
  use simul_box_m
  use units_m
#ifdef HAVE_FFT
  use fft_m
  use cube_function_m
#endif
  use logrid_m
  use poisson_cutoffs_m
  use ps_m
  use specie_m
  use specie_pot_m
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
    epot_generate_gauge_field, &
    epot_forces,               &
    dconmut_vnl_r,             &
    zconmut_vnl_r,             &
    build_local_part_in_real_space


  type epot_t
    ! Classic charges:
    integer :: classic_pot        ! How to include the classic charges
    FLOAT, pointer :: vclassic(:) ! We use it to store the potential of the classic charges

    ! Ions
    FLOAT,       pointer :: vpsl(:)       ! the local part of the pseudopotentials
#ifdef HAVE_FFT
    type(dcf_t), pointer :: local_cf(:)   ! for the local pseudopotential in Fourier space
    type(dcf_t), pointer :: rhocore_cf(:) ! and for the core density
#endif
    integer :: nvnl                       ! number of nonlocal operators
    type(projector_t), pointer :: p(:)    ! non-local projectors
    integer, pointer :: atomproj(:,:)     ! the range of projectors
                                          ! corresponding to an atom
    ! External e-m fields
    integer :: no_lasers                   ! number of laser pulses used
    type(laser_t), pointer :: lasers(:)    ! lasers stuff
    FLOAT, pointer :: E_field(:)           ! static electric field
    FLOAT, pointer :: v_static(:)          ! static scalar potential
    FLOAT, pointer :: B_field(:)           ! static magnetic field
    FLOAT, pointer :: A_static(:,:)        ! static vector potential
    FLOAT, pointer :: A_gauge(:)           ! gauge vector potential
    FLOAT, pointer :: A_gauge_dot(:)       ! dA_gauge/dt
    FLOAT, pointer :: A_gauge_ddot(:)      ! d^2A_gauge/dt^2
    logical :: with_gauge_field            ! true if A_gauge(:) is used
    
    ! additional arbitrary td potential defined by formula string that will be added to the local
    ! part of the potential
    character(len=1024) :: extra_td_pot

    ! The gyromagnetic ratio (-2.0 for the electron, but different if we treat
    ! *effective* electrons in a quantum dot. It affects the spin Zeeman term.
    FLOAT :: gyromagnetic_ratio
    integer :: forces
    logical :: gen_grads
#ifdef HAVE_MPI
    logical :: parallel_generate
#endif
  end type epot_t

  integer, parameter :: &
       DERIVATE_POTENTIAL = 1,      &
       DERIVATE_WAVEFUNCTION  = 2

contains

  ! ---------------------------------------------------------
  subroutine epot_init(ep, gr, geo)
    type(epot_t), intent(out) :: ep
    type(grid_t), intent(in)  :: gr
    type(geometry_t), intent(inout) :: geo

    integer :: i, nvl
    C_POINTER :: blk
    FLOAT, allocatable :: x(:)
    logical :: filter

    call push_sub('epot.epot_init')

    call loct_parse_logical(check_inp('FilterPotentials'), .false., filter)

    if(filter) then
      message(1) = 'Info: filtering the potentials.'
      call write_info(1)
    end if
    
    do i = 1, geo%nspecies
      call specie_pot_init(geo%specie(i), gr, filter)
    end do

    ! Local part of the pseudopotentials
    ALLOCATE(ep%vpsl(NP), NP)
    !$omp parallel workshare
    ep%vpsl(1:NP) = M_ZERO
    !$omp end parallel workshare

#if defined(HAVE_FFT)
    ! should we calculate the local pseudopotentials in Fourier space?
    ! This depends on wether we have periodic dimensions or not
    if(simul_box_is_periodic(gr%sb).and.(.not.geo%only_user_def)) then
      call epot_local_fourier_init(ep, gr%m, gr%sb, geo)
    end if
#endif

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
    if(ep%no_lasers>0.and.mpi_grp_is_root(mpi_world)) then
      message(1) = 'Info: Lasers'
      call write_info(1)

      call laser_write_info(ep%no_lasers, ep%lasers, stdout)
    end if

    !%Variable UserDefinedTDPotential
    !%Type string
    !%Section Hamiltonian
    !%Description
    !% The formula string defined by this variable will be used as additional time dependent
    !% potential (Note: this is available for all specie types, not only for user defined species)
    !%End
    call loct_parse_string(check_inp('UserDefinedTDPotential'), '0', ep%extra_td_pot)

    !%Variable StaticElectricField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% A static constant electrical field may be added to the usual Hamiltonian,
    !% by setting the block StaticElectricField. Atomic units will be assumed
    !% always for its magnitude, regardless of the unit system specified.
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
    
    !%Variable GaugeVectorField
    !%Type block
    !%Section Hamiltonian
    !%Description
    !% The gauge vector field is used to include a uniform (but time dependent)
    !% external electric field in a time dependent run for a periodic system
    !% By default this field is kept null.
    !%End
    ! Read the initial gauge vector field
    ep%with_gauge_field = .false.
    nullify(ep%A_gauge, ep%A_gauge_dot, ep%A_gauge_ddot)
    if(simul_box_is_periodic(gr%sb)) then
      if(loct_parse_block(check_inp('GaugeVectorField'), blk) == 0) then
        ep%with_gauge_field = .true.
        ALLOCATE(ep%A_gauge(NDIM), NDIM)
        ALLOCATE(ep%A_gauge_dot(NDIM), NDIM)
        ALLOCATE(ep%A_gauge_ddot(NDIM), NDIM)
        ep%A_gauge = M_ZERO
	ep%A_gauge_dot = M_ZERO
	ep%A_gauge_ddot = M_ZERO
	do i = 1, NDIM
          call loct_parse_block_float(blk, 0, i-1, ep%A_gauge(i))
	end do
	call loct_parse_block_end(blk)
      end if
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

    !%Variable Forces
    !%Type integer
    !%Default derivate_potential
    !%Section Hamiltonian
    !%Description
    !% The forces can be calculated either by derivating the ionic
    !% potential or the wavefunctions. This option selects how to
    !% calculate them. By default the wavefunctions are derived.
    !% Note: The option derivate_potential is deprecated and will be
    !% eliminated soon.
    !%Option derivate_potential 1
    !% Derivate the potential to calculate the forces.
    !%Option derivate_wavefunctions 2
    !% Derivate the wavefunctions to calculate the forces.
    !%End
    call loct_parse_int(check_inp('Forces'), DERIVATE_WAVEFUNCTION, ep%forces)

    ep%gen_grads = (ep%forces == DERIVATE_POTENTIAL)

#ifdef HAVE_MPI
    !%Variable ParallelPotentialGeneration
    !%Type logical
    !%Default false
    !%Section Generalities::Parallel
    !%Description
    !% If <tt>true</tt> and parallelization in states is used, the
    !% generation of the potential it is done in parallel. This is
    !% still expertimental so it is disabled dy default.
    !%End

    call loct_parse_logical(check_inp('ParallelPotentialGeneration'), .false., ep%parallel_generate)
#endif

    ! The projectors
    ep%nvnl = geometry_nvnl(geo, nvl)

    !if not periodic add also the local potential
    if(.not. simul_box_is_periodic(gr%sb)) ep%nvnl = ep%nvnl + nvl

    nullify(ep%p)
    
    if(ep%nvnl > 0) then
      ALLOCATE(ep%p(ep%nvnl), ep%nvnl)
      ALLOCATE(ep%atomproj(1:2, geo%natoms), 2*geo%natoms)
      do i = 1, ep%nvnl
        call projector_null(ep%p(i))
      end do
    end if

    call pop_sub()
  end subroutine epot_init


  ! ---------------------------------------------------------
  subroutine epot_end(ep, gr, geo)
    type(epot_t),      intent(inout) :: ep
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(inout) :: geo

    integer :: i, iproj

    call push_sub('epot.epot_end')

    do i = 1, geo%nspecies
      call specie_pot_end(geo%specie(i), gr)
    end do

#ifdef HAVE_FFT
    if(simul_box_is_periodic(gr%sb).and.(.not.geo%only_user_def)) then
      do i = 1, geo%nspecies
        call dcf_free(ep%local_cf(i))
        if(geo%specie(i)%nlcc) call dcf_free(ep%rhocore_cf(i))
      end do
      deallocate(ep%local_cf)
      if(geo%nlcc) deallocate(ep%rhocore_cf)
    end if
#endif

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
    if(associated(ep%A_gauge)) deallocate(ep%A_gauge)
    if(associated(ep%A_gauge_dot)) deallocate(ep%A_gauge_dot)
    if(associated(ep%A_gauge_ddot)) deallocate(ep%A_gauge_ddot)

    if(ep%nvnl>0) then
      do iproj = 1, ep%nvnl
        call projector_end(ep%p(iproj))
      end do

      ASSERT(associated(ep%p))
      deallocate(ep%p)
      deallocate(ep%atomproj)
    end if

    call pop_sub()

  end subroutine epot_end


#ifdef HAVE_FFT
  ! ---------------------------------------------------------
  subroutine epot_local_fourier_init(ep, m, sb, geo)
    type(epot_t),      intent(inout) :: ep
    type(mesh_t),      intent(in)    :: m
    type(simul_box_t), intent(in)    :: sb
    type(geometry_t),  intent(in)    :: geo

    integer :: vlocal_cutoff
    integer :: i, ix, iy, iz, ixx(MAX_DIM), db(MAX_DIM), c(MAX_DIM)
    FLOAT :: x(MAX_DIM)
    FLOAT :: gpar, gperp, gx, gz, modg
    FLOAT :: r_0, temp(MAX_DIM), tmp, norm

    type(specie_t), pointer :: s ! shortcuts
    type(dcf_t), pointer :: cf

    call push_sub('epot.epot_local_fourier_init')

    !%Variable VlocalCutoff
    !%Type integer
    !%Default value of PeriodicDimensions
    !%Section Hamiltonian
    !%Description
    !% Define which cutoff type is to be applied to the long range part of the local potential
    !% A cutoff is used when one wants to avoid the long range interactions
    !% among the system enclosed in the simulation box and (some of) its periodic images
    !%Option cutoff_sphere 0
    !% Cut off the interaction out of a sphere
    !%Option cutoff_cylinder 1
    !% Cut off the interaction out of a cylinder with axis parallel to the x direction
    !%Option cutoff_slab 2
    !% Cut off the interaction out of a slab in the xy plane
    !%Option cutoff_none 3
    !% Do not apply any cutoff: all the periodic images interact
    !%End
    call loct_parse_int(check_inp('VlocalCutoff'), sb%periodic_dim , vlocal_cutoff)
    if(.not.varinfo_valid_option('VlocalCutoff', vlocal_cutoff)) call input_error('VlocalCutoff')
    call messages_print_var_option(stdout, "VlocalCutoff", vlocal_cutoff)

    if (vlocal_cutoff /= sb%periodic_dim) then
      write(message(1), '(a,i1,a)')'The System is periodic in ', sb%periodic_dim, ' dimension(s),'
      write(message(2), '(a,i1,a)')'but VlocalCutoff is set for ', vlocal_cutoff, ' dimensions.'
      call write_warning(2)
    end if

    ALLOCATE(ep%local_cf(geo%nspecies), geo%nspecies)
    if(geo%nlcc) ALLOCATE(ep%rhocore_cf(geo%nspecies), geo%nspecies)

    specie: do i = 1, geo%nspecies
      s  => geo%specie(i)
      cf => ep%local_cf(i)

      !%Variable VlocalCutoffRadius
      !%Type float
      !%Default value of the largest nonperiodic box length
      !%Section Hamiltonian
      !%Description
      !% The maximum length out of which the long range part of the interaction
      !% is cut off. 
      !% It refers to the radius of the cylinder if VlocalCutoff = 1,
      !% to the thickness od the slab if VlocalCutoff = 2.
      !%End

      if(i == 1) then
        call mesh_double_box(sb, m, db)
        call dcf_new(db, cf)    ! initialize the cube
        call dcf_fft_init(cf, sb)   ! and initialize the ffts
        db = cf%n               ! dimensions of the box may have been optimized, so get them
        c(:) = db(:)/2 + 1      ! get center of double box
        if (vlocal_cutoff == 3) then
          r_0 = M_ZERO
        else
          call loct_parse_float(check_inp('VlocalCutoffRadius'),&
            maxval(db(:)*m%h(:)/M_TWO)/units_inp%length%factor , r_0)
          r_0 = r_0*units_inp%length%factor
          write(message(1),'(3a,f12.6)')'Info: Vlocal Cutoff Radius [',  &
            trim(units_out%length%abbrev), '] = ',       &
            r_0/units_out%length%factor
          call write_info(1)
        end if
      else
        call dcf_new_from(cf, ep%local_cf(1))   ! we can just copy from the first one
      end if

      if(geo%nlcc) call dcf_new_from(ep%rhocore_cf(i), ep%local_cf(1))

      call dcf_alloc_FS(cf)      ! allocate the tube in Fourier space

      norm    = M_FOUR*M_PI/m%vol_pp(1)
      temp(:) = M_TWO*M_PI/(db(:)*m%h(:))
      cf%FS   = M_Z0
      do ix = 1, cf%nx
        ixx(1) = pad_feq(ix, db(1), .true.)
        do iy = 1, db(2)
          ixx(2) = pad_feq(iy, db(2), .true.)
          do iz = 1, db(3)
            ixx(3) = pad_feq(iz, db(3), .true.)

            x = temp(:)*ixx(:)
            modg = sqrt(sum((temp(:)*ixx(:))**2))

            tmp = specie_get_local_fourier(sb%dim, s, x)
            if(modg /= M_ZERO) then
              tmp = tmp - s%z_val*exp(-(modg/(2*s%ps%a_erf))**2)/modg**2
              select case(vlocal_cutoff)
              case(0)
                cf%FS(ix, iy, iz) = tmp*poisson_cutoff_sphere(modg, r_0)

              case(1)
                gx = abs(temp(1)*ixx(1))
                gperp = sqrt((temp(2)*ixx(2))**2 + (temp(3)*ixx(3))**2)
                cf%FS(ix, iy, iz) = tmp*poisson_cutoff_infinite_cylinder(gx, gperp, r_0)

              case(2)
                gz = abs(temp(3)*ixx(3))
                gpar = sqrt((temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
                cf%FS(ix, iy, iz) = tmp*poisson_cutoff_slab(gpar, gz, r_0)

              case(3)
                cf%FS(ix, iy, iz) = tmp
              end select
            else
              select case(vlocal_cutoff)
              case(0)  ; cf%FS(ix, iy, iz) = -r_0**2/M_TWO
              case(1,2); cf%FS(ix, iy, iz) = M_ZERO
              case(3)  ; cf%FS(ix, iy, iz) = tmp
              end select
            end if

            ! multiply by normalization factor and a phase shift to get the center of the box
            ! the phase is exp(-i kR), where R denotes the center of the box
            ! note that c is a fortran index that starts at 1
            cf%FS(ix, iy, iz) = norm*        &
               exp(-M_ZI*sum(x(:)*(c(:)-1)*m%h(:)))*   &
               cf%FS(ix, iy, iz)
          end do
        end do
      end do

      ! now we built the non-local core corrections in momentum space
      nlcc: if(s%nlcc) then
        call dcf_alloc_RS(ep%rhocore_cf(i))

        do ix = 1, db(1)
          ixx(1) = ix - c(1)
          do iy = 1, db(2)
            ixx(2) = iy - c(2)
            do iz = 1, db(3)
              ixx(3) = iz - c(3)

              x(:) = m%h(:)*ixx(:)
              ep%rhocore_cf(i)%RS(ix, iy, iz) = specie_get_nlcc(s, x)
            end do
          end do
        end do
        call dcf_alloc_FS(ep%rhocore_cf(i))      ! allocate the tube in Fourier space
        call dcf_RS2FS(ep%rhocore_cf(i))         ! Fourier transform
        call dcf_free_RS(ep%rhocore_cf(i))       ! we do not need the real space any longer
      end if nlcc

    end do specie

    call pop_sub()
  end subroutine epot_local_fourier_init
#endif

  ! ---------------------------------------------------------
  subroutine epot_generate_gauge_field(ep, gr, st)
    type(epot_t),      intent(inout) :: ep
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: st
    
    integer :: i, ispin
    FLOAT :: n_el, omega2, vol

    call push_sub('epot.epot_generate_gauge_field')
    
    ASSERT(st%d%wfs_type == M_CMPLX)
    
    ! Integrate the charge density
    n_el = M_ZERO
    do ispin = 1, st%d%spin_channels
      n_el = n_el + dmf_integrate(gr%m, st%rho(1:NP,ispin))
    end do
    
    vol = M_ONE
    do i = 1, NDIM
      ! This only holds for a parallelepipedal box
      vol = vol*gr%sb%rlat(i,i)
    end do
    
    omega2 = M_FOUR*M_PI*P_c*n_el/vol
    !call calc_paramagnetic_current(gr, st, jp)
    
    ! DEBUG
    ! Harmonic oscillator
    ep%A_gauge_ddot(1:NDIM) = -omega2*ep%A_gauge(1:NDIM)
    !write(*,'(a,3f12.6)')'OUT: A = ',ep%A_gauge(1:NDIM)
    !write(*,'(a,3f12.6)')'OUT: dA/dt = ',ep%A_gauge_dot(1:NDIM)
    !write(*,'(a,3f12.6)')'OUT: d^2A/dt^2 = ',ep%A_gauge_ddot(1:NDIM)
    ! END DEBUG
       
    call pop_sub()
  
  end subroutine epot_generate_gauge_field
  
  ! ---------------------------------------------------------
  subroutine epot_generate(ep, gr, geo, mc, st, reltype, time)
    type(epot_t),      intent(inout) :: ep
    type(grid_t), target,  intent(inout) :: gr
    type(geometry_t),  intent(inout) :: geo
    type(multicomm_t), intent(in)    :: mc
    type(states_t),    intent(inout) :: st
    integer,           intent(in)    :: reltype
    FLOAT,   optional, intent(in)    :: time

    FLOAT   :: time_
    integer :: ia, l, lm, k, iproj
    type(atom_t),   pointer :: atm
#ifdef HAVE_FFT
    type(dcf_t) :: cf_loc, cf_nlcc
#endif
    type(mesh_t),      pointer :: m
    type(simul_box_t), pointer :: sb
    type(submesh_t)  :: nl_sphere

    call profiling_in(C_PROFILING_EPOT_GENERATE)
    call push_sub('epot.epot_generate')

    sb  => gr%sb
    m   => gr%m

    time_ = M_ZERO
    if (present(time)) time_ = time

    ! first we assume that we need to recalculate the ion_ion energy
    geo%eii = ion_ion_energy(geo)

#ifdef HAVE_FFT
    if(simul_box_is_periodic(sb).and.(.not.geo%only_user_def)) then
      call dcf_new_from(cf_loc, ep%local_cf(1)) ! at least one specie must exist
      call dcf_alloc_FS(cf_loc)
      cf_loc%FS = M_z0

      if(geo%nlcc) then
        call dcf_new_from(cf_nlcc, ep%local_cf(1)) ! at least one specie must exist
        call dcf_alloc_FS(cf_nlcc)
        cf_nlcc%FS = M_z0
      end if
    end if
#endif

    ! Local.
    ep%vpsl = M_ZERO
    do ia = 1, geo%natoms
      atm => geo%atom(ia) ! shortcuts

      if((.not.simul_box_is_periodic(sb)).or.geo%only_user_def) then

        call build_local_part_in_real_space(ep, gr, geo, atm, ep%vpsl, time_, st%rho_core)

#ifdef HAVE_FFT
      else ! momentum space
        call cf_phase_factor(sb, m, atm%x, ep%local_cf(atm%spec%index), cf_loc)
        if(atm%spec%nlcc) then
          call cf_phase_factor(sb, m, atm%x, ep%rhocore_cf(atm%spec%index), cf_nlcc)
        end if
#endif
      end if

    end do

    ! the pseudo potential part.
    iproj = 1
    do ia = 1, geo%natoms
      atm => geo%atom(ia)

      if(.not. specie_is_ps(atm%spec)) cycle

      ep%atomproj(1, ia) = iproj
      
      call submesh_init_sphere(nl_sphere, sb, m, atm%x, atm%spec%ps%rc_max)

      do l = 0, atm%spec%ps%l_max
        if(atm%spec%ps%l_loc == l) cycle
        do lm = -l, l

          call projector_end(ep%p(iproj))
          call submesh_copy(nl_sphere, ep%p(iproj)%sphere)
          if(simul_box_is_periodic(sb)) call projector_init_phases(ep%p(iproj), sb, m, atm, st)
          call projector_init(ep%p(iproj), atm, reltype, l, lm)

          ep%p(iproj)%iatom = ia
          iproj = iproj + 1
        end do
      end do

      call submesh_end(nl_sphere)

      if(.not. simul_box_is_periodic(gr%sb)) then
        !the local part
        ep%p(iproj)%iatom = ia

        call projector_end(ep%p(iproj))
        call submesh_init_sphere(ep%p(iproj)%sphere, &
             sb, m, atm%x, double_grid_get_rmax(gr%dgrid, atm%spec, m))
        call projector_init(ep%p(iproj), atm, force_type = M_LOCAL)

        iproj = iproj + 1
      end if

      ep%atomproj(2, ia) = iproj - 1

    end do

    call projector_build_all

#ifdef HAVE_FFT
    if(simul_box_is_periodic(sb).and.(.not.geo%only_user_def)) then
      ! first the potential
      call dcf_alloc_RS(cf_loc)
      call dcf_FS2RS(cf_loc)
      call dcf2mf(m, cf_loc, ep%vpsl)
      call dcf_free(cf_loc)

      ! and the non-local core corrections
      if(geo%nlcc) then
        call dcf_alloc_RS(cf_nlcc)
        call dcf_FS2RS(cf_nlcc)
        call dcf2mf(m, cf_nlcc, st%rho_core)
        call dcf_free(cf_nlcc)
      end if
    end if
#endif

    if (ep%classic_pot > 0) then
      ep%vpsl(1:m%np) = ep%vpsl(1:m%np) + ep%vclassic(1:m%np)
    end if

    call pop_sub()
    call profiling_out(C_PROFILING_EPOT_GENERATE)

  contains

    subroutine projector_build_all
#ifdef HAVE_MPI
      integer :: rank, size, ini, fin, mpi_err
      integer, allocatable :: rep(:)
      
      if (.not. (ep%parallel_generate .and. multicomm_strategy_is_parallel(mc, P_STRATEGY_STATES))) then
#endif

        do iproj = 1, ep%nvnl
          call projector_build(ep%p(iproj), gr, geo%atom(ep%p(iproj)%iatom), ep%gen_grads)
        end do
        
#ifdef HAVE_MPI
      else
        
        size = mc%group_sizes(P_STRATEGY_STATES)
        
        ALLOCATE(rep(1:ep%nvnl), ep%nvnl)
        
        do rank = 0, size-1
          ini = rank * ep%nvnl / size + 1
          fin = min((rank + 1 )* ep%nvnl / size, ep%nvnl)
          rep(ini:fin) = rank
        end do
        
        rank = mc%who_am_i(P_STRATEGY_STATES) - 1
        
        do iproj = 1, ep%nvnl
          if ( rep(iproj) == rank ) then 
            call projector_build(ep%p(iproj), gr, geo%atom(ep%p(iproj)%iatom), ep%gen_grads)
          end if
        end do
        
        do iproj = 1, ep%nvnl
          call projector_broadcast(ep%p(iproj), gr, mc, geo%atom(ep%p(iproj)%iatom), ep%gen_grads, rep(iproj))
        end do
        
      deallocate(rep)

    end if

#endif
    end subroutine projector_build_all

  end subroutine epot_generate

  subroutine build_local_part_in_real_space(ep, gr, geo, a, vpsl, time, rho_core)
    type(epot_t),             intent(in)    :: ep
    type(grid_t),             intent(inout) :: gr
    type(geometry_t),         intent(in)    :: geo
    type(atom_t),             intent(inout) :: a
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,                    intent(in)    :: time
    FLOAT, optional,          intent(inout) :: rho_core(:)
 
    integer :: i
    FLOAT :: x(MAX_DIM), xx(MAX_DIM), r, pot_re, pot_im
    FLOAT, allocatable  :: rho(:), vl(:)

    call push_sub('epot.build_local_part_in_real_space')

    ALLOCATE(vl(1:NP_PART), NP_PART)

    !Local potential
    call specie_get_local(a%spec, gr, a%x(1:NDIM), vl, time)
    !$omp parallel workshare
    vpsl(1:NP) = vpsl(1:NP) + vl(1:NP)
    !$omp end parallel workshare

    !Non-local core corrections
    if(present(rho_core) .and. a%spec%nlcc .and. specie_is_ps(a%spec)) then
      do i = 1, NP
        x(1:NDIM) = gr%m%x(i, 1:NDIM) - a%x(1:NDIM)
        rho_core(i) = rho_core(i) + specie_get_nlcc(a%spec, x)
      end do
    end if

    !Time dependent potential
    if(time > M_ZERO .and. (ep%extra_td_pot .ne. '0') ) then
      do i = 1, NP
        x(1:NDIM) = gr%m%x(i, 1:NDIM) - a%x(1:NDIM)
        xx(1:NDIM) = x(1:NDIM) / units_inp%length%factor   ! convert from a.u. to input units
        r = sqrt(sum(x(1:NDIM)**2)) / units_inp%length%factor
        call loct_parse_expression(pot_re, pot_im, xx(1), xx(2), xx(3), &
             r, time, ep%extra_td_pot)
        vpsl(i) = vpsl(i) + pot_re * units_inp%energy%factor  ! convert from input units to a.u.
      end do
    end if

    !Local potential from density
    if(a%spec%has_density) then

      ALLOCATE(rho(1:NP), NP)

      call specie_get_density(a%spec, a%x, gr, geo, rho)
      call dpoisson_solve(gr, vl, rho)
      vpsl(1:NP) = vpsl(1:NP) + vl(1:NP)

      deallocate(rho)

    end if

    deallocate(vl)

    call pop_sub()
  end subroutine build_local_part_in_real_space


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
    type(geometry_t), intent(inout)  :: geo
    type(epot_t),     intent(in)     :: ep
    type(states_t),   intent(inout)     :: st
    FLOAT,     optional, intent(in)    :: t

    integer :: i, j, l
    FLOAT :: d, r, zi, zj, x(MAX_DIM), time
    type(atom_t), pointer :: atm

#if defined(HAVE_MPI)
    FLOAT :: f(MAX_DIM)
#endif

    time = M_ZERO
    if(present(t)) time = t

    call profiling_in(C_PROFILING_FORCES)
    call push_sub('epot.epot_forces')

    ! init to 0
    do i = 1, geo%natoms
      geo%atom(i)%f = M_ZERO
    end do
    
    if (wfs_are_real(st) ) then 
      call dcalc_forces_from_potential(gr, geo, ep, st, time)
    else
      call zcalc_forces_from_potential(gr, geo, ep, st, time)
    end if

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      do i = 1, geo%natoms
        atm => geo%atom(i)
        call MPI_Allreduce(atm%f(1), f(1), NDIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        atm%f = f
      end do
    end if
#endif

    if(.not.geo%only_user_def) then ! exclude user defined species for the moment
      ! Now the ion, ion force term
      do i = 1, geo%natoms
        zi = geo%atom(i)%spec%Z_val
        do j = 1, geo%natoms
          if(i .ne. j) then
            zj = geo%atom(j)%spec%Z_val
            r = sqrt(sum((geo%atom(i)%x(1:NDIM) - geo%atom(j)%x(1:NDIM))**2))
            d = zi * zj/r**3

            geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
                 d*(geo%atom(i)%x(1:NDIM) - geo%atom(j)%x(1:NDIM))
          end if
        end do
      end do
    end if

#if defined(HAVE_FFT)
    if( simul_box_is_periodic(gr%sb) .and. (.not. geo%only_user_def) ) then ! fourier space
      call local_FS()
    end if
#endif

    !TODO: forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, ep%no_lasers
        select case(ep%lasers(j)%field)
        case(E_FIELD_ELECTRIC)
          call laser_field(gr%sb, ep%lasers(j), x, t)
          do i = 1, geo%natoms
            geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
              geo%atom(i)%spec%Z_val * x(1:NDIM)
          end do
        case(E_FIELD_MAGNETIC, E_FIELD_VECTOR_POTENTIAL)
          write(message(1),'(a)') 'The forces are currently not properly calculated if time-dependent'
          write(message(2),'(a)') 'magnetic fields are present.'
          call write_fatal(2)
        end select
      end do
    end if

    if(associated(ep%E_field)) then
      do i = 1, geo%natoms
        geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
             geo%atom(i)%spec%Z_val * ep%E_field(1:NDIM)
      end do
    end if


    call pop_sub()
    call profiling_out(C_PROFILING_FORCES)

#ifdef HAVE_FFT

  contains

    ! ---------------------------------------------------------
    subroutine local_FS()
      type(dcf_t) :: cf_for
      FLOAT, allocatable :: force(:)
      
      ALLOCATE(force(NP), NP)
      call dcf_new_from(cf_for, ep%local_cf(1)) ! at least one specie must exist
      call dcf_alloc_FS(cf_for)
      call dcf_alloc_RS(cf_for)

      do i = 1, geo%natoms
        atm => geo%atom(i)
        do j = 1, NDIM
          cf_for%FS = M_z0
          call cf_phase_factor(gr%sb, gr%m, atm%x, ep%local_cf(atm%spec%index), cf_for)

          call dcf_FS_grad(gr%sb, gr%m, cf_for, j)
          call dcf_FS2RS(cf_for)
          call dcf2mf(gr%m, cf_for, force)
          do l = 1, st%d%nspin
            ! FIXME: When running with partitions, vol_pp is local
            ! to the node. It is likely, that this code need changes.
            atm%f(j) = atm%f(j) + sum(force(1:NP)*st%rho(1:NP, l)*gr%m%vol_pp(1:NP))
          end do
        end do
      end do

      call dcf_free(cf_for)
      deallocate(force)
    end subroutine local_FS
#endif

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
