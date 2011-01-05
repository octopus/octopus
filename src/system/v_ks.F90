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

module v_ks_m
  use berry_m
  use datasets_m
  use density_m
  use derivatives_m
  use energy_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_function_m
  use lalg_basic_m
  use magnetic_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use multigrid_m
  use parser_m
  use poisson_m
  use poisson_sete_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_dim_m
  use unit_system_m
  use varinfo_m
  use xc_m
  use XC_F90(lib_m)
  use xc_functl_m
  use xc_ks_inversion_m
  use xc_OEP_m

  implicit none

  private
  public ::             &
    v_ks_t,             &
    v_ks_init,          &
    v_ks_end,           &
    v_ks_write_info,    &
    v_ks_calc,          &
    v_ks_calc_start,    &
    v_ks_calc_finish,   &
    v_ks_hartree,       &
    v_ks_freeze_hxc

  integer, parameter, public :: &
    sic_none   = 1,     &  ! no self-interaction correction
    sic_pz     = 2,     &  ! Perdew-Zunger SIC (OEP way)
    sic_amaldi = 3         ! Amaldi correction term

  type v_ks_calc_t
    logical                       :: calculating
    type(hamiltonian_t),  pointer :: hm
    logical                       :: time_present
    FLOAT                         :: time
    logical                       :: calc_berry
    FLOAT,                pointer :: rho(:, :)
    FLOAT,                pointer :: total_rho(:)
  end type v_ks_calc_t

  type v_ks_t
    integer :: theory_level

    logical :: frozen_hxc ! For RPA and SAE calculations.

    integer                  :: xc_family  ! the XC stuff
    integer                  :: sic_type   ! what kind of self-interaction correction to apply
    type(xc_t)               :: xc
    type(xc_OEP_t)           :: oep
    type(xc_ks_inversion_t)  :: ks_inversion
    type(poisson_t), pointer :: hartree_solver
    logical                  :: new_hartree
    logical                  :: tail_correction
    FLOAT                    :: tail_correction_tol
    type(grid_t), pointer    :: gr
    type(v_ks_calc_t)        :: calc
  end type v_ks_t

contains

  
  ! ---------------------------------------------------------
  subroutine v_ks_init(ks, gr, dd, geo, mc, nel)
    type(v_ks_t),         intent(out)   :: ks
    type(grid_t), target, intent(inout) :: gr
    type(states_dim_t),   intent(in)    :: dd
    type(geometry_t),     intent(inout) :: geo
    type(multicomm_t),    intent(in)    :: mc  
    FLOAT,                intent(in)    :: nel ! the total number of electrons

    PUSH_SUB(v_ks_init)

    !%Variable TheoryLevel
    !%Type integer
    !%Default dft
    !%Section Hamiltonian
    !%Description
    !% The calculations can be run with four different "theory levels":
    !%Option independent_particles 2
    !% Particles will be considered as independent, <i>i.e.</i> as non-interacting.
    !% This mode is mainly used for testing purposes, as the code is usually 
    !% much faster with <tt>independent_particles</tt>.
    !%Option hartree 1
    !% Calculation within the Hartree method. Note that, contrary to popular
    !% belief, the Hartree potential is self-interaction-free. Therefore, this run 
    !% mode will not yield the same result as <tt>dft</tt> without exchange-correlation.
    !% Development version only.
    !%Option hartree_fock 3
    !% This is the traditional Hartree-Fock scheme. Like the Hartree scheme, it is fully
    !% self-interaction-free. This mode is extremely slow. It is often more convenient
    !% to use <tt>dft</tt> within the OEP scheme to get similar (but not the same) results.
    !% Note that within this scheme you can use a correlation functional, or a hybrid
    !% functional (see <tt>XCFunctional</tt>). In the latter case, you will be following the
    !% quantum-chemistry recipe to use hybrids.
    !%Option dft 4
    !% This is the default density-functional theory scheme. Note that you can also use 
    !% hybrids in this scheme, but they will be handled the "DFT" way, <i>i.e.</i>, solving the
    !% OEP equation.
    !%End
    call parse_integer(datasets_check('TheoryLevel'), KOHN_SHAM_DFT, ks%theory_level)
    if(.not.varinfo_valid_option('TheoryLevel', ks%theory_level)) call input_error('TheoryLevel')

    call messages_obsolete_variable('NonInteractingElectrons', 'TheoryLevel')
    call messages_obsolete_variable('HartreeFock', 'TheoryLevel')

    select case(ks%theory_level)
    case(INDEPENDENT_PARTICLES)
      ks%sic_type = sic_none
    case(HARTREE)
      call messages_devel_version("Hartree theory level")
    case(HARTREE_FOCK)
      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, dd%spin_channels, dd%cdft, hartree_fock=.true.)
      ks%xc_family = ks%xc%family
      ks%sic_type = sic_none

    case(KOHN_SHAM_DFT)
      ! initialize XC modules
      call xc_init(ks%xc, gr%mesh%sb%dim, nel, dd%spin_channels, dd%cdft, hartree_fock=.false.)
      ks%xc_family = ks%xc%family

      ! check for SIC
      if(iand(ks%xc_family, XC_FAMILY_LDA + XC_FAMILY_GGA) .ne. 0) then

        !%Variable SICCorrection
        !%Type integer
        !%Default sic_none
        !%Section Hamiltonian::XC
        !%Description
        !% This variable controls which form of self-interaction correction to use. Note that
        !% this correction will be applied to the functional chosen by <tt>XCFunctional</tt>.
        !% All except <tt>sic_none</tt> are development-version only.
        !%Option sic_none 1
        !% No self-interaction correction.
        !%Option sic_pz 2
        !% Perdew-Zunger SIC, handled by the OEP technique.
        !%Option sic_amaldi 3
        !% Amaldi correction term (NOT WORKING).
        !%End
        call parse_integer(datasets_check('SICCorrection'), sic_none, ks%sic_type)
        if(.not. varinfo_valid_option('SICCorrection', ks%sic_type)) call input_error('SICCorrection')
        if(ks%sic_type .ne. sic_none) call messages_devel_version("self-interaction correction")

        ! Perdew-Zunger corrections
        if(ks%sic_type == sic_pz) ks%xc_family = ior(ks%xc_family, XC_FAMILY_OEP)
      
      else
        ks%sic_type = sic_none
      end if

      !%Variable XCTailCorrection
      !%Type logical
      !%Default no
      !%Section Hamiltonian::XC
      !%Description
      !% (Experimental) This variable enables a mechanism to correct
      !% the value of the XC functional in near-zero-density
      !% regions. This zone might have numerical noise or it might
      !% even be set to zero by <tt>libxc</tt>.
      !% The correction is performed by calculating the "XC density", 
      !% a charge density that would generate an
      !% equivalent XC potential. This XC density is a localized
      !% quantity, so we can remove the noise by cutting in the
      !% regions where the density is small.
      !%End
      call parse_logical(datasets_check('XCTailCorrection'), .false., ks%tail_correction)
      
      if(ks%tail_correction) then 
        call messages_devel_version("XC tail correction")
        write(message(1),'(a)') 'This correction shouldn''t be used with systems having nodal points of the electron density.'
        call write_info(1)
            
        !%Variable XCTailCorrectionTol
        !%Type float
        !%Default 5-e12
        !%Section Hamiltonian::XC
        !%Description
        !% This variable sets the threshold at which to cut the XC density when
        !% <tt>XCTailCorrection</tt> is enabled. The value is always assumed to
        !% be in atomic units.
        !%End
        call parse_float(datasets_check('XCTailCorrectionTol'), CNST(5e-12), ks%tail_correction_tol)
      end if

      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_init(ks%oep, ks%xc_family, gr, dd)
      endif
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_init(ks%ks_inversion, ks%xc_family, gr, geo, dd, mc)
      endif
    end select

    ks%frozen_hxc = .false.

    call v_ks_write_info(ks, stdout)

    ks%new_hartree = .false.
    nullify(ks%hartree_solver)
    if(ks%theory_level /= INDEPENDENT_PARTICLES) then
      if(gr%have_fine_mesh) then
        ks%new_hartree = .true.
        SAFE_ALLOCATE(ks%hartree_solver)
        call poisson_init(ks%hartree_solver, gr%fine%der, geo)
      else
        ks%hartree_solver => psolver
      end if
     end if
     
     ks%gr => gr
     ks%calc%calculating = .false.

    POP_SUB(v_ks_init)
  end subroutine v_ks_init
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_end(ks, gr, geo)
    type(v_ks_t),     intent(inout) :: ks
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(inout) :: geo

    PUSH_SUB(v_ks_end)

    select case(ks%theory_level)
    case(KOHN_SHAM_DFT)
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_end(ks%ks_inversion, gr, geo)
      endif
      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_end(ks%oep)
      endif
      call xc_end(ks%xc)
    end select

    if(ks%new_hartree) then
      SAFE_DEALLOCATE_P(ks%hartree_solver)
    end if

    POP_SUB(v_ks_end)
  end subroutine v_ks_end
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_write_info(ks, iunit)
    type(v_ks_t), intent(in) :: ks
    integer,      intent(in) :: iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(v_ks_write_info)

    call messages_print_stress(iunit, "Theory Level")
    call messages_print_var_option(iunit, "TheoryLevel", ks%theory_level)

    select case(ks%theory_level)
    case(HARTREE_FOCK)
      write(iunit, '(1x)')
      call xc_write_info(ks%xc, iunit)

    case(KOHN_SHAM_DFT)
      write(iunit, '(1x)')
      call xc_write_info(ks%xc, iunit)

      write(iunit, '(1x)')
      call messages_print_var_option(iunit, 'SICCorrection', ks%sic_type)

      if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
        call xc_oep_write_info(ks%oep, iunit)
      end if
      if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
        call xc_ks_inversion_write_info(ks%ks_inversion, iunit)
      end if
    end select

    call messages_print_stress(iunit)

    POP_SUB(v_ks_write_info)
  end subroutine v_ks_write_info
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_calc(ks, hm, st, calc_eigenval, time, calc_berry)
    type(v_ks_t),           intent(inout) :: ks
    type(hamiltonian_t),    intent(inout) :: hm
    type(states_t),         intent(inout) :: st
    logical,      optional, intent(in)    :: calc_eigenval
    FLOAT,        optional, intent(in)    :: time
    logical,      optional, intent(in)    :: calc_berry ! use this before wfns initialized
    
    call v_ks_calc_start(ks, hm, st, time, calc_berry)
    call v_ks_calc_finish(ks)

    if(present(calc_eigenval)) then
      if(calc_eigenval) call energy_calculate_eigenvalues(hm, ks%gr%der, st, &
        open_boundaries = ks%gr%ob_grid%open_boundaries)
    end if

  end subroutine v_ks_calc

  ! --------------------------------------------------------- 

  !> This routine starts the calculation of the Kohn-Sham
  !! potential. v_ks_calc_finish must be called to finish the
  !! calculation. The argument hm CANNOT be used until
  !! v_ks_calc_finish has been called. The st argument, on the other
  !! hand, is released inmediately and CAN be modified before
  !! finishing the calculation.

  subroutine v_ks_calc_start(ks, hm, st, time, calc_berry) 
    type(v_ks_t),                      intent(inout) :: ks 
    type(hamiltonian_t),     target,   intent(inout) :: hm 
    type(states_t),                    intent(inout) :: st
    FLOAT,                   optional, intent(in)    :: time 
    logical,                 optional, intent(in)    :: calc_berry !use this before wfns initialized
    
    FLOAT :: amaldi_factor, distance
    integer :: ip, ispin
    type(profile_t), save :: prof
    
    ! The next line is a hack to be able to perform an IP/RPA calculation
    !logical, save :: RPA_first = .true.

    PUSH_SUB(v_ks_calc_start)
    call profiling_in(prof, "KOHN_SHAM_CALC")

    ASSERT(.not. ks%calc%calculating)
    ks%calc%calculating = .true.

    if(in_debug_mode) then
      write(message(1), '(a)') 'Debug: Calculating Kohn-Sham potential.'
      call write_info(1)
    end if

    ks%calc%calc_berry = .true.
    if(present(calc_berry)) ks%calc%calc_berry = calc_berry

    ks%calc%time_present = present(time)
    if(present(time)) ks%calc%time = time

    ks%calc%hm => hm

    ! If the Hxc term is frozen, there is nothing to do (WARNING: MISSING hm%epot)
    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_start)
      return
    end if

    hm%epot = M_ZERO

    ! check whether we should introduce the Amaldi SIC correction
    amaldi_factor = M_ONE
    if(ks%sic_type == sic_amaldi) amaldi_factor = (st%qtot - M_ONE)/st%qtot

    if(ks%theory_level == INDEPENDENT_PARTICLES .or. amaldi_factor == M_ZERO) then
      hm%vhxc     = M_ZERO
      hm%epot     = M_ZERO
      hm%ehartree = M_ZERO
      hm%ex       = M_ZERO
      hm%ec       = M_ZERO
    else
      ! The next 2 lines are a hack to be able to perform an IP/RPA calculation
      !if(RPA_first) then
      !  RPA_first = .false.

        hm%ehartree = M_ZERO
        call v_ks_hartree(ks, st, hm, amaldi_factor)

        hm%vxc      = M_ZERO
        if(iand(hm%xc_family, XC_FAMILY_MGGA) .ne. 0) hm%vtau = M_ZERO
        if(hm%d%cdft) hm%axc = M_ZERO
        if(ks%theory_level .ne. HARTREE) call v_a_xc()

      !end if

      ! Build Hartree + XC potential

       forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 1) = hm%vxc(ip, 1) + hm%vhartree(ip)

      if(hm%d%ispin > UNPOLARIZED) then
        forall(ip = 1:ks%gr%mesh%np) hm%vhxc(ip, 2) = hm%vxc(ip, 2) + hm%vhartree(ip)
      end if

      if(hm%d%ispin == SPINORS) then
        forall(ispin = 3:4, ip = 1:ks%gr%mesh%np) hm%vhxc(ip, ispin) = hm%vxc(ip, ispin)
      end if

    end if

    if(associated(hm%ep%E_field) .and. simul_box_is_periodic(ks%gr%mesh%sb)) then
      if(ks%calc%calc_berry) then
        call berry_potential(st, ks%gr%mesh, hm%ep%E_field, hm%vberry)
      else
        ! before wfns are initialized, cannot calculate this term
        hm%vberry(1:ks%gr%mesh%np, 1:hm%d%nspin) = M_ZERO
      endif
    endif

    if(ks%theory_level == HARTREE .or. ks%theory_level == HARTREE_FOCK) then
      call states_end(hm%st)
      call states_copy(hm%st, st)
    end if
    if(ks%theory_level == HARTREE_FOCK) then
      hm%exx_coef = ks%xc%exx_coef
    else if (ks%theory_level==HARTREE) then
      hm%exx_coef = M_ONE
    end if

    ! Calculate the vector potential induced by the electronic current.
    ! WARNING: calculating the self-induced magnetic field here only makes
    ! sense if it is going to be used in the Hamiltonian, which does not happen
    ! now. Otherwise one could just calculate it at the end of the calculation.
    if(hm%self_induced_magnetic) call magnetic_induced(ks%gr%der, st, hm%a_ind, hm%b_ind)

    call profiling_out(prof)
    POP_SUB(v_ks_calc_start)

  contains

    ! ---------------------------------------------------------
    subroutine v_a_xc()

      FLOAT, allocatable :: rho(:, :)
      FLOAT, allocatable :: totalrho(:)
      type(profile_t), save :: prof
      FLOAT, pointer :: vxc(:, :)
      FLOAT, allocatable :: vxcc(:), nxc(:)
      integer :: ispin, ierr, ip, itmp, idim
      character(len=10) :: vxc_name


      PUSH_SUB(v_ks_calc_start.v_a_xc)
      call profiling_in(prof, "XC")

      hm%ex = M_ZERO
      hm%ec = M_ZERO
      hm%exc_j = M_ZERO

      ! get density taking into account non-linear core corrections, and the Amaldi SIC correction
      SAFE_ALLOCATE(rho(1:ks%gr%fine%mesh%np, 1:st%d%nspin))
      
      call states_total_density(st, ks%gr%fine%mesh, rho)

      ! Amaldi correction
      if(ks%sic_type == sic_amaldi) then
        rho(1:ks%gr%fine%mesh%np, 1:st%d%nspin) = amaldi_factor * rho(1:ks%gr%fine%mesh%np, 1:st%d%nspin)
      end if

      if(ks%gr%have_fine_mesh) then
        SAFE_ALLOCATE(vxc(1:ks%gr%fine%mesh%np_part, 1:st%d%nspin))
        vxc = M_ZERO
      else
        vxc => hm%vxc
      end if

      ! Get the *local* XC term
      if(hm%d%cdft) then
        call xc_get_vxc_and_axc(ks%gr%fine%der, ks%xc, st, rho, st%current, st%d%ispin, hm%vxc, hm%axc, &
             hm%ex, hm%ec, hm%exc_j, -minval(st%eigenval(st%nst, :)), st%qtot)
      else
        call xc_get_vxc(ks%gr%fine%der, ks%xc, st, rho, st%d%ispin, hm%ex, hm%ec, &
             -minval(st%eigenval(st%nst, :)), st%qtot, vxc, vtau=hm%vtau)
      end if

      if(ks%theory_level == KOHN_SHAM_DFT) then
        ! The OEP family has to be handled specially
        if(iand(ks%xc_family, XC_FAMILY_OEP) .ne. 0) then
          if (states_are_real(st)) then
            call dxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
              ks%gr, hm, st, hm%ex, hm%ec, vxc=hm%vxc)
          else
            call zxc_oep_calc(ks%oep, ks%xc, (ks%sic_type==sic_pz),  &
              ks%gr, hm, st, hm%ex, hm%ec, vxc=hm%vxc)
          end if
        endif

        if(iand(ks%xc_family, XC_FAMILY_KS_INVERSION) .ne. 0) then
          ! Also treat KS inversion separately (not part of libxc)
          call xc_ks_inversion_calc(ks%ks_inversion, ks%gr, hm, st, hm%ex, hm%ec, vxc=hm%vxc)
        endif
      end if

      if(ks%tail_correction) then

        SAFE_ALLOCATE(vxcc(1:ks%gr%fine%mesh%np_part))
        SAFE_ALLOCATE(totalrho(1:ks%gr%fine%mesh%np))
        
        totalrho(:) = M_ZERO
        do ip = 1, ks%gr%fine%mesh%np
          do ispin= 1, st%d%nspin
            totalrho(ip) = totalrho(ip) + rho(ip,ispin)
          end do
        end do
        
        do ispin = 1, st%d%nspin
                    
          vxcc(1:ks%gr%fine%mesh%np_part) = hm%vxc(1:ks%gr%fine%mesh%np_part, ispin)
          
         
          
          ! These output calls and the ones below are for debugging, XA and FB.
          
          !first we set the names of the output files

          write (vxc_name,'(i10)') ispin
          itmp = verify(vxc_name," ")
          vxc_name =  "vxc"//trim(vxc_name(itmp:))
          
                              
          !print the XC potential before the correction
          
          call doutput_function(output_axis_x, "./static", vxc_name, ks%gr%fine%mesh, vxcc, unit_one, ierr)
          call doutput_function(output_axis_y, "./static", vxc_name, ks%gr%fine%mesh, vxcc, unit_one, ierr)
          call doutput_function(output_axis_z, "./static", vxc_name, ks%gr%fine%mesh, vxcc, unit_one, ierr)

          !Performing the correction to the "XC density" 
          do ip = 1, ks%gr%fine%mesh%np
            if( (rho(ip,ispin) .ne. M_ZERO) .and. (totalrho(ip) < ks%tail_correction_tol) ) then
              distance = M_ZERO
              do idim = 1,ks%gr%mesh%sb%dim 
                distance = distance + ks%gr%fine%mesh%x(ip,idim)**2
              end do
              distance = sqrt(distance)
              vxcc(ip) =  -1/distance
            end if
          end do
          
                    
          !print the XC potential after the correction
          !call doutput_function(output_axis_x, "./static", trim(nxc_name)//trim("cut") , ks%gr%fine%mesh, nxc, unit_one, ierr)
          !call doutput_function(output_axis_x, "./static", trim(vxc_name)//trim("cut") , ks%gr%fine%mesh, vxcc, unit_one, ierr)
          
          hm%vxc(1:ks%gr%fine%mesh%np_part, ispin) = vxcc(1:ks%gr%fine%mesh%np_part)

        end do

        SAFE_DEALLOCATE_A(nxc)
        SAFE_DEALLOCATE_A(vxcc)

      end if

      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(totalrho)
      
      ! Now we calculate Int[n vxc] = hm%epot
      select case(hm%d%ispin)
      case(UNPOLARIZED)
        hm%epot = hm%epot + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), vxc(:, 1))
      case(SPIN_POLARIZED)
        hm%epot = hm%epot + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), vxc(:, 1)) &
             + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 2), vxc(:, 2))
      case(SPINORS)
        hm%epot = hm%epot + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 1), vxc(:, 1)) &
             + dmf_dotp(ks%gr%fine%mesh, st%rho(:, 2), vxc(:, 2)) &
             + M_TWO*dmf_dotp(ks%gr%fine%mesh, st%rho(:, 3), vxc(:, 3)) &
             + M_TWO*dmf_dotp(ks%gr%fine%mesh, st%rho(:, 4), vxc(:, 4))

      end select

      if(ks%gr%have_fine_mesh) then
        do ispin = 1, st%d%nspin
          call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, vxc(:, ispin), hm%vxc(:, ispin), INJECTION)
          ! some debugging output that I will keep here for the moment, XA
          !          call doutput_function(1, "./", "vxc_fine", ks%gr%fine%mesh, vxc(:, ispin), unit_one, ierr)
          !          call doutput_function(1, "./", "vxc_coarse", ks%gr%mesh, hm%vxc(:, ispin), unit_one, ierr)
        end do
        SAFE_DEALLOCATE_P(vxc)
      end if

      call profiling_out(prof)
      POP_SUB(v_ks_calc_start.v_a_xc)
    end subroutine v_a_xc
  end subroutine v_ks_calc_start
  ! ---------------------------------------------------------

  subroutine v_ks_calc_finish(ks)
    type(v_ks_t),        intent(inout) :: ks

    PUSH_SUB(v_ks_calc_finish)

    if(ks%frozen_hxc) then
      POP_SUB(v_ks_calc_finish)
      return
    end if

    ASSERT(ks%calc%calculating)
    ks%calc%calculating = .false.

    if(ks%calc%time_present) then
      call hamiltonian_update(ks%calc%hm, ks%gr%mesh, ks%calc%time)
    else
      call hamiltonian_update(ks%calc%hm, ks%gr%mesh)
    end if

    POP_SUB(v_ks_calc_finish)
  end subroutine v_ks_calc_finish

  ! ---------------------------------------------------------
  ! Hartree contribution to the XC potential
  subroutine v_ks_hartree(ks, st, hm, amaldi_factor)
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(in)    :: st
    FLOAT, optional,     intent(in)    :: amaldi_factor

    FLOAT, allocatable :: rho(:)
    FLOAT, pointer :: pot(:)
    integer :: is, ip

    PUSH_SUB(v_ks_hartree)

    ASSERT(associated(ks%hartree_solver))

    SAFE_ALLOCATE(rho(1:ks%gr%fine%mesh%np))

    ! calculate the total density
    call lalg_copy(ks%gr%fine%mesh%np, st%rho(:, 1), rho)

    do is = 2, hm%d%spin_channels
      forall(ip = 1:ks%gr%fine%mesh%np) rho(ip) = rho(ip) + st%rho(ip, is)
    end do

    ! Add, if it exists, the frozen density from the inner orbitals.
    if(associated(st%frozen_rho)) then
      do is = 1, hm%d%spin_channels
        forall(ip = 1:ks%gr%fine%mesh%np) rho(ip) = rho(ip) + st%frozen_rho(ip, is)
      end do
    end if

    ! Amaldi correction
    if(present(amaldi_factor)) rho = amaldi_factor*rho

    if(.not. ks%gr%have_fine_mesh) then
      pot => hm%vhartree
    else
      SAFE_ALLOCATE(pot(1:ks%gr%fine%mesh%np_part))
      pot = M_ZERO
    end if

    ! solve the Poisson equation
    call dpoisson_solve(ks%hartree_solver, pot, rho)
    
    ! Get the Hartree energy
    hm%ehartree = M_HALF * dmf_dotp(ks%gr%fine%mesh, rho, pot)

    if(ks%gr%have_fine_mesh) then
      ! we use injection to transfer to the fine grid, we cannot use
      ! restriction since the boundary conditions are not zero for the
      ! Hartree potential (and for some XC functionals).
      call dmultigrid_fine2coarse(ks%gr%fine%tt, ks%gr%fine%der, ks%gr%mesh, pot, hm%vhartree, INJECTION)
      ! some debugging output that I will keep here for the moment, XA
      !      call doutput_function(1, "./", "vh_fine", ks%gr%fine%mesh, pot, unit_one, is)
      !      call doutput_function(1, "./", "vh_coarse", ks%gr%mesh, hm%vhartree, unit_one, is)
      SAFE_DEALLOCATE_P(pot)
    end if

    if (poisson_get_solver(ks%hartree_solver) == POISSON_SETE) then !SEC
      hm%ehartree=hm%ehartree + poisson_energy(ks%hartree_solver)
      write(89,*) hm%ehartree*CNST(2.0*13.60569193), poisson_energy(ks%hartree_solver)*CNST(2.0*13.60569193), &
        hm%ep%eii*CNST(2.0*13.60569193)
    endif

    SAFE_DEALLOCATE_A(rho)
    POP_SUB(v_ks_hartree)
  end subroutine v_ks_hartree
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine v_ks_freeze_hxc(ks)
    type(v_ks_t), intent(inout) :: ks

    PUSH_SUB(v_ks_freeze_hxc)

    ks%frozen_hxc = .true.
    
    POP_SUB(v_ks_freeze_hxc)
  end subroutine v_ks_freeze_hxc
  ! ---------------------------------------------------------

end module v_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
