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
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use simul_box_m
  use units_m
#ifdef HAVE_FFT
  use fft_m
  use cube_function_m
#endif
  use logrid_m
  use ps_m
  use specie_m
  use specie_pot_m
  use geometry_m
  use states_m
  use lasers_m
  use profiling_m
  use mpi_m
  use mpi_debug_m
  use varinfo_m
  use poisson_m
  use hgh_projector_m
  use kb_projector_m
  use rkb_projector_m

  implicit none

  private
  public ::                    &
    epot_t,                    &
    epot_init,                 &
    epot_end,                  &
    epot_generate,             &
    epot_generate_gauge_field, &
    epot_laser_scalar_pot,     &
    epot_laser_field,          &
    epot_laser_vector_pot,     &
    epot_forces,               &
    dproject, zproject,        &
    projector_t

  ! The projector data type is intended to hold the non-local part of the
  ! pseudopotentials. The definition of the action of a projector (which is
  ! done through the X(project) subroutine) depends on the type of the 
  ! projector. There are three different types:
  !  -> HGH projector
  !  -> "normal" Kleinman-Bylander projector (no spin-orbit)
  !  -> "relativistc" Kleinman-Bylander projector (includes spin-orbit)
  type projector_t
    private
    integer :: type
    integer :: iatom

    integer          :: n_s     ! number of points inside the sphere
    integer, pointer :: jxyz(:) ! index of the points inside the sphere
    integer          :: nik
    CMPLX,   pointer :: phases(:,:)

    ! Only one of the following structures should be used at once
    ! The one to be use depends on the value of type variable
    type(hgh_projector_t), pointer :: hgh_p
    type(kb_projector_t),  pointer :: kb_p
    type(rkb_projector_t), pointer :: rkb_p
  end type projector_t

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

  end type epot_t

  integer, parameter :: M_HGH = 1, &
                        M_KB  = 2, &
                        M_RKB = 3

contains

  ! ---------------------------------------------------------
  subroutine epot_init(ep, gr, geo)
    type(epot_t), intent(out) :: ep
    type(grid_t), intent(in)  :: gr
    type(geometry_t), intent(inout) :: geo

    integer :: i
    C_POINTER :: blk
    FLOAT, allocatable :: x(:)

    call push_sub('epot.epot_init')

    ! Local part of the pseudopotentials
    ALLOCATE(ep%vpsl(NP), NP)
    ep%vpsl = M_ZERO

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
    !% by setting the block StaticMagneticField. Atomic units will be assumed
    !% always for its magnitude, regardless of the unit system specified.
    !% The three possible components of the block (which should only have one
    !% line) are the three components of the magnetic field vector. Note that
    !% if you are running the code in 1D mode this will not work, and if you
    !% are running the code in 2D mode the magnetic field will have to be in
    !% the z-direction, so that the first two columns should be zero.
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
    !% gyromagnetic factor that dependes on the material.
    !%End
    call loct_parse_float(check_inp('GyromagneticRatio'), P_g, ep%gyromagnetic_ratio)

    ! The projectors
    ep%nvnl = geometry_nvnl(geo)
    nullify(ep%p)
    
    if(ep%nvnl > 0) then
      ALLOCATE(ep%p(ep%nvnl), ep%nvnl)
      do i = 1, ep%nvnl
        call projector_null(ep%p(i))
      end do
    end if

    do i = 1, geo%nspecies
      call specie_pot_init(geo%specie(i), gr)
    end do

    call pop_sub()
  end subroutine epot_init


  ! ---------------------------------------------------------
  subroutine epot_end(ep, gr, geo)
    type(epot_t),      intent(inout) :: ep
    type(grid_t),      intent(in)    :: gr
    type(geometry_t),  intent(inout) :: geo

#ifdef HAVE_FFT
    integer :: i
#endif

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
      ASSERT(associated(ep%p))
      deallocate(ep%p)
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
                cf%FS(ix, iy, iz) = tmp*cutoff0(modg, r_0)
              case(1)
                gx = abs(temp(1)*ixx(1))
                gperp = sqrt((temp(2)*ixx(2))**2 + (temp(3)*ixx(3))**2)
                cf%FS(ix, iy, iz) = tmp*cutoff1(gx, gperp, r_0)
              case(2)
                gz = abs(temp(3)*ixx(3))
                gpar = sqrt((temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
                cf%FS(ix, iy, iz) = tmp*cutoff2(gpar, gz, r_0)
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
    FLOAT :: jp(NP, NDIM, st%d%nspin), n_el, omega2, vol

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
  subroutine epot_generate(ep, gr, geo, st, reltype, time, fast_generation)
    type(epot_t),      intent(inout) :: ep
    type(grid_t), target,  intent(inout) :: gr
    type(geometry_t),  intent(inout) :: geo
    type(states_t),    intent(inout) :: st
    integer,           intent(in)    :: reltype
    FLOAT,   optional, intent(in)    :: time
    logical, optional, intent(in)    :: fast_generation

    logical :: fast_generation_
    FLOAT   :: time_
    integer :: ia, i, l, lm, k, p
    type(specie_t), pointer :: s
    type(atom_t),   pointer :: a
#ifdef HAVE_FFT
    type(dcf_t) :: cf_loc, cf_nlcc
#endif
    type(mesh_t),      pointer :: m
    type(simul_box_t), pointer :: sb


    call profiling_in(C_PROFILING_EPOT_GENERATE)
    call push_sub('epot.epot_generate')

    sb  => gr%sb
    m   => gr%m

    fast_generation_ = .false.
    if (present(fast_generation)) fast_generation_ = fast_generation
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
      a => geo%atom(ia) ! shortcuts
      s => a%spec
      if (specie_is_ps(s) .and. conf%debug_level > 0) call debug_pseudo()
      call build_local_part()
    end do

    ! Nonlocal part.
    i = 1
    do ia = 1, geo%natoms
      a => geo%atom(ia)
      s => a%spec
      if(s%local) cycle
      p = 1

      do l = 0, s%ps%l_max
        if(s%ps%l_loc == l) cycle
        do lm = -l, l

          if(.not.fast_generation_) then
            ! This if is a performance hack, necessary for when the ions move.
            ! For each atom, the sphere is the same, so we just calculate it once
            if(p == 1) then
              k = i
              p = 2
              call projector_end(ep%p(i))
              call projector_build_kb_sphere(ep%p(i), sb, m, a, st)
            else
              call projector_copy_kb_sphere(ep%p(k), ep%p(i))
            end if

          end if

          call projector_init_nl_part(ep%p(i), gr, a, reltype, l, lm)

          ep%p(i)%iatom = ia
          i = i + 1
        end do
      end do
    end do

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
    ! ---------------------------------------------------------
    subroutine debug_pseudo()
      integer :: iunit, ii, nn
      FLOAT :: r, dr

      type(loct_spline_t) :: pot_corr
      
      call loct_spline_init(pot_corr)
      call dg_get_potential_correction(gr%dgrid, pot_corr)

      call io_mkdir('debug/')
      iunit = io_open('debug/pseudo-'//trim(s%label), action='write')

      dr = CNST(0.01)
      nn = CNST(10.0)/dr
      r = M_ZERO
      do ii = 1, nn
        write(iunit, '(4f12.6)') r, loct_splint(s%ps%vl, r),  &
             -s%z_val*loct_splint(pot_corr, r), &
             loct_splint(s%ps%vll, r)
        r = r + dr
      end do

      call io_close(iunit)
      
      call loct_spline_end(pot_corr)

    end subroutine debug_pseudo

    ! ---------------------------------------------------------
    subroutine build_local_part()
      integer :: i
      FLOAT :: x(MAX_DIM), xx(MAX_DIM), r, pot_re, pot_im
      FLOAT, allocatable  :: rho(:), vl(:)

      type(loct_spline_t) :: pot_corr

      call push_sub('epot.build_local_part')

      if((.not.simul_box_is_periodic(sb)).or.geo%only_user_def) then
        !Real space

        ALLOCATE(vl(1:m%np_part), m%np_part)
        
        !Local potential
        call specie_get_local(s, gr, a%x(:), vl, time_)
        ep%vpsl(1:m%np) = ep%vpsl(1:m%np) + vl(1:m%np)

        !Non-local core corrections
        if(s%nlcc .and. specie_is_ps(s)) then
          do i = 1, m%np
            x(:) = m%x(i, :) - a%x(:)
            st%rho_core(i) = st%rho_core(i) + specie_get_nlcc(s, x)
          end do
        end if
        
        !Time dependent potential
        if(time_ > M_ZERO .and. (ep%extra_td_pot .ne. '0') ) then
          do i = 1, m%np
            x(:) = m%x(i, :) - a%x(:)
            xx(:) = x(:)/units_inp%length%factor   ! convert from a.u. to input units
            r = sqrt(sum(x(:)**2))/units_inp%length%factor
            call loct_parse_expression(pot_re, pot_im, xx(1), xx(2), xx(3), &
              r, time_, ep%extra_td_pot)
            ep%vpsl(i) = ep%vpsl(i) + pot_re * units_inp%energy%factor  ! convert from input units to a.u.
          end do
        end if

        !Local potential from density
        if(s%has_density .or. &
             (specie_is_ps(s) .and. dg_add_localization_density(gr%dgrid) )) then 

          ALLOCATE(rho(1:m%np), m%np)

          call specie_get_density(s, a%x, gr, geo, rho)
          call dpoisson_solve(gr, vl, rho)
          ep%vpsl(1:m%np) = ep%vpsl(1:m%np) + vl(1:m%np)

          if (specie_is_ps(s)) then 
            
            call loct_spline_init(pot_corr)
            call dg_get_potential_correction(gr%dgrid, pot_corr)
            
            !calculate the deviation from the analitcal potential
            do i = 1, m%np
              x(:) = m%x(i, :) - a%x(:)
              r = sqrt(sum(x(:)**2))
              rho(i) = vl(i) - (-s%z_val)*loct_splint(pot_corr, r)
            end do
            call loct_spline_end(pot_corr)
            
            write(message(1),'(a, e12.6)')  'Info: Deviation from analitical potential is ', abs(dmf_integrate(m, rho))
            call write_info(1)

          end if

          deallocate(rho)

        end if

        deallocate(vl)

#ifdef HAVE_FFT
      else ! momentum space
        call cf_phase_factor(sb, m, a%x, ep%local_cf(s%index), cf_loc)
        if(s%nlcc) then
          call cf_phase_factor(sb, m, a%x, ep%rhocore_cf(s%index), cf_nlcc)
        end if
#endif
      end if

      call pop_sub()
    end subroutine build_local_part

  end subroutine epot_generate


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
  function epot_laser_scalar_pot(np, gr, ep, t) result(v)
    integer,      intent(in) :: np
    type(grid_t), intent(in) :: gr
    type(epot_t), intent(in) :: ep
    FLOAT,        intent(in) :: t

    FLOAT :: v(np)

    call laser_potential(gr%sb, ep%no_lasers, ep%lasers, t, gr%m, v)

  end function epot_laser_scalar_pot


  ! ---------------------------------------------------------
  subroutine epot_laser_vector_pot(sb, ep, t, a)
    type(simul_box_t), intent(in)  :: sb
    type(epot_t),      intent(in)  :: ep
    FLOAT,             intent(in)  :: t
    FLOAT,             intent(out) :: a(sb%dim)

    call laser_vector_field(sb, ep%no_lasers, ep%lasers, t, a)

  end subroutine epot_laser_vector_pot


  ! ---------------------------------------------------------
  subroutine epot_laser_field(sb, ep, t, e)
    type(simul_box_t), intent(in)  :: sb
    type(epot_t),      intent(in)  :: ep
    FLOAT,             intent(in)  :: t
    FLOAT,             intent(out) :: e(sb%dim)

    call laser_field(sb, ep%no_lasers, ep%lasers, t, e)

  end subroutine epot_laser_field


  ! ---------------------------------------------------------
  subroutine epot_forces(gr, geo, ep, st, t)
    type(grid_t), target, intent(inout) :: gr
    type(geometry_t), intent(inout)  :: geo
    type(epot_t),     intent(in)     :: ep
    type(states_t),   intent(in)     :: st
    FLOAT,     optional, intent(in)    :: t

    integer :: i, j, l, ist, ik, ivnl, ivnl_start, ivnl_end
    FLOAT :: d, r, zi, zj, x(MAX_DIM)
    type(atom_t), pointer :: atm

    FLOAT, allocatable :: dppsi(:, :)
    FLOAT :: dz(MAX_DIM)
    CMPLX, allocatable :: zppsi(:, :)
    CMPLX :: zz(MAX_DIM)

#if defined(HAVE_MPI)
    FLOAT :: f(MAX_DIM)
#endif

    call profiling_in(C_PROFILING_FORCES)
    call push_sub('epot.epot_forces')

    ! init to 0
    do i = 1, geo%natoms
      geo%atom(i)%f = M_ZERO
    end do

    if(.not.geo%only_user_def) then
      ! non-local component of the potential.
      if (st%d%wfs_type == M_REAL) then
        ALLOCATE(dppsi(gr%m%np, st%d%dim), gr%m%np*st%d%dim)
      else
        ALLOCATE(zppsi(gr%m%np, st%d%dim), gr%m%np*st%d%dim)
      end if

      atm_loop: do i = 1, geo%natoms
        atm => geo%atom(i)
        if(atm%spec%local) cycle

        ASSERT(NDIM == 3)

        ! Here we learn which are the projector that correspond to atom i.
        ! It assumes that the projectors of each atom are consecutive.
        ivnl_start  = - 1
        do ivnl = 1, ep%nvnl
          if(ep%p(ivnl)%iatom .eq. i) then
            ivnl_start = ivnl
            exit
          end if
        end do
        if(ivnl_start .eq. -1) cycle
        ivnl_end = ep%nvnl
        do ivnl = ivnl_start, ep%nvnl
          if(ep%p(ivnl)%iatom .ne. i) then
            ivnl_end = ivnl - 1
            exit
          end if
        end do

        ik_loop: do ik = 1, st%d%nik
          st_loop: do ist = st%st_start, st%st_end

            if (st%d%wfs_type == M_REAL) then
              dz = dpsidprojectpsi(gr%m, ep%p(ivnl_start:ivnl_end), &
                   ivnl_end - ivnl_start + 1, st%d%dim, st%dpsi(:, :, ist, ik), &
                   periodic = .false., ik = ik)
              atm%f = atm%f + M_TWO * st%occ(ist, ik) * dz
            else
              zz = zpsidprojectpsi(gr%m, ep%p(ivnl_start:ivnl_end), &
                   ivnl_end - ivnl_start + 1, st%d%dim, st%zpsi(:, :, ist, ik), &
                   periodic = .false., ik = ik)
              atm%f = atm%f + M_TWO * st%occ(ist, ik) * zz
            end if
            
          end do st_loop
        end do ik_loop

      end do atm_loop
      if (st%d%wfs_type == M_REAL) then
        deallocate(dppsi)
      else
        deallocate(zppsi)
      end if
    end if

#if defined(HAVE_MPI)
    do i = 1, geo%natoms
      atm => geo%atom(i)
      if(atm%spec%local) cycle

      if(st%parallel_in_states) then
        call MPI_Allreduce(atm%f(1), f(1), NDIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
        atm%f = f
      end if
    end do
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

    ! now comes the local part of the PP
    if(.not.simul_box_is_periodic(gr%sb).or.geo%only_user_def) then ! Real space
      call local_RS()
#if defined(HAVE_FFT)
    else ! Fourier space
      call local_FS()
#endif
    end if

    if(present(t).and.ep%no_lasers>0) then
      call laser_field(gr%sb, ep%no_lasers, ep%lasers, t, x)
      do i = 1, geo%natoms
        geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
             geo%atom(i)%spec%Z_val * x(1:NDIM)
      end do
    end if

    if(associated(ep%E_field)) then
      do i = 1, geo%natoms
        geo%atom(i)%f(1:NDIM) = geo%atom(i)%f(1:NDIM) + &
             geo%atom(i)%spec%Z_val * ep%E_field(1:NDIM)
      end do
    end if

    !TODO: forces due to the magnetic fields (static and time-dependent)

    call pop_sub()
    call profiling_out(C_PROFILING_FORCES)

  contains

    ! ---------------------------------------------------------
    subroutine local_RS()
      FLOAT :: r
      FLOAT, allocatable :: force(:,:), grho(:,:), gpot(:)
      integer  :: i, j, k, ns
      FLOAT :: x(1:MAX_DIM)
      
      ns = min(2, st%d%nspin)

      ALLOCATE(force(NP, MAX_DIM), NP*MAX_DIM)

      do i = 1, geo%natoms
        atm => geo%atom(i)

        call specie_get_glocal(atm%spec, gr, atm%x, force)

        do j = 1, NP
          force(j, 1:NDIM) = sum(st%rho(j, 1:ns))*force(j, 1:NDIM)
        end do

        do k = 1, NDIM
          atm%f(k) = atm%f(k) - dmf_integrate(gr%m, force(:, k))
        end do
      end do

      deallocate(force)

    end subroutine local_RS


#ifdef HAVE_FFT
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


  ! ---------------------------------------------------------
  subroutine projector_null(p)
    type(projector_t), intent(out) :: p

    p%type = 0
    nullify(p%phases)
    nullify(p%jxyz)
    nullify(p%hgh_p)
    nullify(p%kb_p)
    nullify(p%rkb_p)

  end subroutine projector_null

  !---------------------------------------------------------
  subroutine projector_build_kb_sphere(p, sb, m, a, st)
    type(projector_t), intent(inout) :: p
    type(simul_box_t), intent(in)    :: sb
    type(mesh_t),      intent(in)    :: m
    type(atom_t),      intent(in)    :: a
    type(states_t),    intent(in)    :: st

    integer :: i, j, k
    FLOAT :: r, X(MAX_DIM)

    call push_sub('epot.projector_build_kb_sphere')

    if (any(a%spec%ps%rc_max + m%h(1) >= sb%lsize(1:sb%periodic_dim))) then
      message(1)='KB sphere is larger than the box size'
      write(message(2),'(a,f12.6,a)')  '  rc_max+h = ', a%spec%ps%rc_max + m%h(1), ' [b]'
      write(message(3),'(a,3f12.4,a)') '  lsize    = ', sb%lsize, ' [b]'
      message(4)='Please change pseudopotential'
      call write_fatal(4)
    end if

    ! Get the total number of points inside the sphere
    j = 0
    do k = 1, m%np
      do i = 1, 3**sb%periodic_dim
        call mesh_r(m, k, r, a=a%x + sb%shift(i,:))
        if(r > a%spec%ps%rc_max + m%h(1)) cycle
        j = j + 1
        exit
      end do
    end do

    p%n_s = j
    ALLOCATE(p%jxyz(j), j)

    ! Get the index of the points inside the sphere
    j = 0
    do k = 1, m%np
      do i = 1, 3**sb%periodic_dim
        call mesh_r(m, k, r, a=a%x + sb%shift(i,:))
        ! we enlarge slightly the mesh (good for the interpolation scheme)
        if(r > a%spec%ps%rc_max + m%h(1)) cycle
        j = j + 1
        p%jxyz(j) = k
        exit
      end do
    end do

    ! and here the phases for the periodic systems
    if (sb%periodic_dim /= 0) then
      p%nik = st%d%nik

      ALLOCATE(p%phases(p%n_s, p%nik), p%n_s*p%nik)

      do j = 1, p%n_s
        x(:) = m%x(p%jxyz(j), :)
        do k = 1, p%nik
          p%phases(j, k) = exp(M_zI*sum(st%d%kpoints(:, k)*x(:)))
        end do
      end do

    end if

    call pop_sub()
  end subroutine projector_build_kb_sphere

  !---------------------------------------------------------
  subroutine projector_copy_kb_sphere(pi, po)
    type(projector_t), intent(in)    :: pi
    type(projector_t), intent(inout) :: po

    call push_sub('epot.projector_copy_kb_sphere')

    ASSERT(.not. associated(po%jxyz))

    po%n_s = pi%n_s
    ALLOCATE(po%jxyz(pi%n_s), pi%n_s)
    po%jxyz = pi%jxyz

    if (associated(pi%phases)) then
      po%nik = pi%nik
      ALLOCATE(po%phases(pi%n_s, pi%nik), pi%n_s*pi%nik)
      po%phases = pi%phases
    end if

    call pop_sub()
  end subroutine projector_copy_kb_sphere

  !---------------------------------------------------------
  subroutine projector_init_nl_part(p, gr, a, reltype, l, lm)
    type(projector_t), intent(inout) :: p
    type(grid_t),      intent(in)    :: gr
    type(atom_t),      intent(in)    :: a
    integer,           intent(in)    :: reltype
    integer,           intent(in)    :: l, lm

    call push_sub('epot.projector_init_nl_part')

    select case (a%spec%ps%kbc)
    case (1)
      p%type = M_KB
      if (reltype == 1) then
        write(message(1),'(a,a,a)') "Spin-orbit coupling for specie ", trim(a%spec%label), " is not available."
        call write_warning(1)
      end if
    case (2)
      if (l == 0 .or. reltype == 0) then
        p%type = M_KB
      else
        p%type = M_RKB
      end if
    case (3)
      p%type = M_HGH
    end select

    select case (p%type)
    case (M_HGH)
      ALLOCATE(p%hgh_p, 1)
      call hgh_projector_null(p%hgh_p)
      call hgh_projector_init(p%hgh_p, p%n_s, p%jxyz, gr, a, l, lm)
    case (M_KB)
      ALLOCATE(p%kb_p, 1)
      call kb_projector_null(p%kb_p)
      call kb_projector_init(p%kb_p, p%n_s, p%jxyz, gr, a, l, lm)
    case (M_RKB)
      ALLOCATE(p%rkb_p, 1)
      call rkb_projector_null(p%rkb_p)
      call rkb_projector_init(p%rkb_p, p%n_s, p%jxyz, gr, a, l, lm)
    end select

    call pop_sub()
  end subroutine projector_init_nl_part

  !---------------------------------------------------------
  subroutine projector_end(p)
    type(projector_t), intent(inout) :: p

    call push_sub('epot.projector_end')

    if (associated(p%jxyz)) deallocate(p%jxyz)
    if (associated(p%phases)) deallocate(p%phases)
    if (associated(p%hgh_p)) then
      call hgh_projector_end(p%hgh_p)
      deallocate(p%hgh_p)
    end if
    if (associated(p%kb_p)) then
      call kb_projector_end(p%kb_p)
      deallocate(p%kb_p)
    end if
    if (associated(p%rkb_p)) then
      call rkb_projector_end(p%rkb_p)
      deallocate(p%rkb_p)
    end if

    call pop_sub()
  end subroutine projector_end

#include "undef.F90"
#include "real.F90"
#include "epot_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "epot_inc.F90"

end module external_pot_m


