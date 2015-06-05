!! Copyright (C) 2010 H. Appel, N. Helbig
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
!! $Id$

#include "global.h"

module xc_ks_inversion_m
  use density_m
  use derivatives_m
  use eigensolver_m 
  use geometry_m 
  use global_m 
  use grid_m 
  use hamiltonian_m 
  use io_m 
  use io_function_m
  use lalg_adv_m 
  use mesh_function_m 
  use mesh_m 
  use messages_m 
  use multicomm_m
  use parser_m 
  use poisson_m 
  use profiling_m 
  use states_m 
  use states_dim_m 
  use unit_m 
  use unit_system_m 
  use varinfo_m 
  use XC_F90(lib_m) 
  use xc_m 
  use xc_functl_m

  implicit none

  private
  public ::                        &
    xc_ks_inversion_t,             &
    xc_ks_inversion_init,          &
    xc_ks_inversion_end,           &
    xc_ks_inversion_write_info,    &
    xc_ks_inversion_calc,          &
    invertks_2part,                &
    invertks_iter

  !> KS inversion methods/algorithms
  integer, public, parameter ::      &
    XC_INV_METHOD_VS_ITER      = 1,  &
    XC_INV_METHOD_TWO_PARTICLE = 2

  !> the KS inversion levels
  integer, public, parameter ::      &
    XC_KS_INVERSION_NONE      = 1,   &
    XC_KS_INVERSION_ADIABATIC = 2,   &   
    XC_KS_INVERSION_TD_EXACT  = 3

  !> asymptotic correction for v_xc
  integer, public, parameter ::      &
    XC_ASYMPTOTICS_NONE    = 1,      &
    XC_ASYMPTOTICS_SC      = 2     

  type xc_ks_inversion_t
     integer             :: method
     integer             :: level
     integer             :: asymp
     FLOAT, pointer      :: vxc_previous_step(:,:)
     type(states_t)      :: aux_st
     type(hamiltonian_t) :: aux_hm
     type(eigensolver_t) :: eigensolver
  end type xc_ks_inversion_t


contains

  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_init(ks_inv, gr, geo, st, mc)
    type(xc_ks_inversion_t), intent(out)   :: ks_inv
    type(grid_t),            intent(inout) :: gr
    type(geometry_t),        intent(inout) :: geo
    type(states_t),          intent(in)    :: st
    type(multicomm_t),       intent(in)    :: mc

    PUSH_SUB(xc_ks_inversion_init)

    if(mc%n_node > 1) &
      call messages_not_implemented("Kohn-Sham inversion in parallel")

    call messages_experimental("Kohn-Sham inversion")
    
    !%Variable InvertKSmethod
    !%Type integer
    !%Default iterative
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects whether the exact two-particle method or the iterative scheme
    !% is used to invert the density to get the KS potential.
    !%Option iterative 1
    !% Iterative scheme for <math>v_s</math>.
    !%Option two_particle 2
    !% Exact two-particle scheme.
    !%End
    call parse_variable('InvertKSmethod', XC_INV_METHOD_VS_ITER, ks_inv%method)

    if(ks_inv%method < XC_INV_METHOD_VS_ITER &
      .or. ks_inv%method > XC_INV_METHOD_TWO_PARTICLE) then
      call messages_input_error('InvertKSmethod')
      call messages_fatal(1)
    end if

    !%Variable KSInversionLevel
    !%Type integer
    !%Default ks_inversion_adiabatic
    !%Section Calculation Modes::Invert KS
    !%Description
    !% At what level <tt>Octopus</tt> shall handle the KS inversion.
    !%Option ks_inversion_none 1
    !% Do not compute KS inversion.
    !%Option ks_inversion_adiabatic 2
    !% Compute exact adiabatic <math>v_{xc}</math>.
    !%End
    call messages_obsolete_variable('KS_Inversion_Level', 'KSInversionLevel')
    call parse_variable('KSInversionLevel', XC_KS_INVERSION_ADIABATIC, ks_inv%level)
    if(.not.varinfo_valid_option('KSInversionLevel', ks_inv%level)) call messages_input_error('KSInversionLevel')

    !%Variable KSInversionAsymptotics
    !%Type integer
    !%Default xc_asymptotics_none
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Asymptotic correction applied to <math>v_{xc}</math>.
    !%Option xc_asymptotics_none 1
    !% Do not apply any correction in the asymptotic region.
    !%Option xc_asymptotics_sc 2
    !% Applies the soft-Coulomb decay of <math>-1/\sqrt{r^2+1}</math> to <math>v_{xc}</math> in the asymptotic region.
    !%End
    call parse_variable('KSInversionAsymptotics', XC_ASYMPTOTICS_NONE, ks_inv%asymp)

    if(ks_inv%level /= XC_KS_INVERSION_NONE) then
      call states_copy(ks_inv%aux_st, st, exclude_wfns = .true.)
      
      ! initialize auxiliary random wavefunctions
      call states_allocate_wfns(ks_inv%aux_st, gr%mesh)
      call states_generate_random(ks_inv%aux_st, gr%mesh)      

      ! initialize densities, hamiltonian and eigensolver
      call states_densities_init(ks_inv%aux_st, gr, geo)
      call hamiltonian_init(ks_inv%aux_hm, gr, geo, ks_inv%aux_st, INDEPENDENT_PARTICLES, XC_FAMILY_NONE)
      call eigensolver_init(ks_inv%eigensolver, gr, ks_inv%aux_st)
    end if

    POP_SUB(xc_ks_inversion_init)
  end subroutine xc_ks_inversion_init


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_end(ks_inv)
    type(xc_ks_inversion_t), intent(inout) :: ks_inv

    PUSH_SUB(xc_ks_inversion_end)

    if(ks_inv%level /= XC_KS_INVERSION_NONE) then
      ! cleanup
      call eigensolver_end(ks_inv%eigensolver)
      call hamiltonian_end(ks_inv%aux_hm)
      call states_end(ks_inv%aux_st)
    end if

    POP_SUB(xc_ks_inversion_end)
  end subroutine xc_ks_inversion_end


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_write_info(ks_inversion, iunit)
    type(xc_ks_inversion_t), intent(in) :: ks_inversion
    integer,                 intent(in) :: iunit

    if(ks_inversion%level == XC_KS_INVERSION_NONE) return

    PUSH_SUB(xc_ks_inversion_write_info)
    call messages_print_var_option(iunit, 'KSInversionLevel', ks_inversion%level)

    POP_SUB(xc_ks_inversion_write_info)
  end subroutine xc_ks_inversion_write_info


  ! ---------------------------------------------------------
  subroutine invertks_2part(target_rho, nspin, aux_hm, gr, st, eigensolver, asymptotics)
    FLOAT,               intent(in)    :: target_rho(:,:) !< (1:gr%mesh%np, 1:nspin)
    integer,             intent(in)    :: nspin
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(eigensolver_t), intent(inout) :: eigensolver
    integer,             intent(in)    :: asymptotics
           
    integer :: ii, jj, asym1, asym2 
    integer :: np
    FLOAT   :: rr, shift, smalldensity
    FLOAT, allocatable :: sqrtrho(:,:), laplace(:,:), vks(:,:)

    PUSH_SUB(invertks_2part)
    
    np = gr%mesh%np
    
    SAFE_ALLOCATE(sqrtrho(1:gr%der%mesh%np_part, 1:nspin))
    SAFE_ALLOCATE(vks(1:np, 1:nspin))
    SAFE_ALLOCATE(laplace(1:gr%der%mesh%np, 1:nspin))
    
    sqrtrho = M_ZERO
    smalldensity = 5d-6

    if(any(target_rho(:,:) < -M_EPSILON)) then
      write(message(1),*) "Target density has negative points. min value = ", minval(target_rho(:,:))
      call messages_warning(1)
    end if
    
    do jj = 1, nspin
      do ii = 1, gr%der%mesh%np
        sqrtrho(ii, jj) = sqrt(target_rho(ii, jj))
        !if (sqrtrho(ii, jj) < CNST(2.5e-6)) sqrtrho(ii, jj) = CNST(2.5e-6)
      end do
    end do   
    
    do jj = 1, nspin
      call dderivatives_lapl(gr%der, sqrtrho(:,jj), laplace(:,jj))
    end do
    
    do ii = 1, nspin
      !avoid division by zero and set parameters for asymptotics
      !only for 1D potentials at the moment
      !need to find a way to find all points from where asymptotics should start in 2 and 3D
      do jj = 1, int(np/2)
        if(target_rho(jj,ii) < smalldensity) then
          vks(jj, ii) = aux_hm%ep%vpsl(jj) + aux_hm%vhartree(jj)
          asym1 = jj
        end if
        if(target_rho(np-jj+1, ii) < smalldensity) then
          vks(np-jj+1, ii) = aux_hm%ep%vpsl(np-jj+1) + aux_hm%vhartree(np-jj+1)
          asym2 = np - jj + 1
        end if
      end do
      do jj = asym1+1, asym2-1
        vks(jj, ii) = laplace(jj, ii)/(M_TWO*sqrtrho(jj, ii))
      end do
      aux_hm%vxc(:,ii) = vks(:,ii) - aux_hm%ep%vpsl(:) - aux_hm%vhartree(:)
    end do

    !ensure correct asymptotic behavior, only for 1D potentials at the moment
    !need to find a way to find all points from where asymptotics should start in 2 and 3D
    if(asymptotics == XC_ASYMPTOTICS_SC) then
      do ii = 1, nspin
        do jj = 1, asym1
          call mesh_r(gr%mesh, jj, rr)
          aux_hm%vxc(jj, ii) = -1.0/sqrt(rr**2 + 1.0)
        end do
     
        ! calculate constant shift for correct asymptotics and shift accordingly
        call mesh_r(gr%mesh, asym1+1, rr)
        shift  = aux_hm%vxc(asym1+1, ii) + 1.0/sqrt(rr**2 + 1.0)
        do jj = asym1+1, asym2-1
          aux_hm%vxc(jj,ii) = aux_hm%vxc(jj, ii) - shift
        end do
  
        call mesh_r(gr%mesh, asym2-1, rr)
        shift  = aux_hm%vxc(asym2-1, ii) + 1.0/sqrt(rr**2 + 1.0)
        do jj = 1, asym2-1
          aux_hm%vxc(jj,ii) = aux_hm%vxc(jj, ii) - shift
        end do
        do jj = asym2, np
          call mesh_r(gr%mesh, jj, rr)
          aux_hm%vxc(jj, ii) = -1.0/sqrt(rr**2 + 1.0)
        end do
      end do 
    end if !apply asymptotic correction
    
    if(asymptotics == XC_ASYMPTOTICS_NONE) then
      do ii = 1, nspin
        ! calculate constant shift to make potential continuous
        shift  = aux_hm%vxc(asym1+1, ii)! + aux_hm%ep%vpsl(asym1+1) + aux_hm%vhartree(asym1+1)
        do jj = asym1+1, asym2-1
          aux_hm%vxc(jj,ii) = aux_hm%vxc(jj, ii) - shift
        end do
  
        shift  = aux_hm%vxc(asym2-1, ii)!+ aux_hm%ep%vpsl(asym2-1) + aux_hm%vhartree(asym2-1)
 
        do jj = 1, asym2-1
          aux_hm%vxc(jj,ii) = aux_hm%vxc(jj, ii) - shift
        end do
      end do 
    end if  

    do ii = 1, nspin
      aux_hm%vhxc(:,ii) = aux_hm%vxc(:,ii) + aux_hm%vhartree(:)
    end do
    
    call hamiltonian_update(aux_hm, gr%mesh)
    call eigensolver_run(eigensolver, gr, st, aux_hm, 1)
    call density_calc(st, gr, st%rho)
    
    SAFE_DEALLOCATE_A(sqrtrho)
    SAFE_DEALLOCATE_A(laplace)
    SAFE_DEALLOCATE_A(vks)

    POP_SUB(invertks_2part)
  end subroutine invertks_2part


  ! ---------------------------------------------------------
  subroutine invertks_iter(target_rho, nspin, aux_hm, gr, st, eigensolver, asymptotics)
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(eigensolver_t), intent(inout) :: eigensolver
    integer,             intent(in)    :: nspin
    FLOAT,               intent(in)    :: target_rho(1:gr%mesh%np, 1:nspin)
    integer,             intent(in)    :: asymptotics
        
    integer :: ii, jj, ierr, asym1, asym2
    integer :: iunit, verbosity, counter, np
    integer :: max_iter
    FLOAT :: rr, shift
    FLOAT :: alpha, beta
    FLOAT :: convdensity, diffdensity
    FLOAT, allocatable :: vhxc(:,:)

    character(len=256) :: fname

    PUSH_SUB(invertks_iter)

    np = gr%mesh%np

    !%Variable InvertKSConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Absolute difference between the calculated and the target density in the KS
    !% inversion. Has to be larger than the convergence of the density in the SCF run.
    !%End    
    call parse_variable('InvertKSConvAbsDens', CNST(1e-5), convdensity)

    !%Variable InvertKSVerbosity
    !%Type integer
    !%Default 0
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects what is output during the calculation of the KS potential.
    !%Option 0
    !% Only outputs the converged density and KS potential.
    !%Option 1
    !% Same as 0 but outputs the maximum difference to the target density in each
    !% iteration in addition.
    !%Option 2
    !% Same as 1 but outputs the density and the KS potential in each iteration in 
    !% addition.
    !%End
    call parse_variable('InvertKSVerbosity', 0, verbosity)  
    if(verbosity < 0 .or. verbosity > 2) then
      call messages_input_error('InvertKSVerbosity')
      call messages_fatal(1)
    end if

    !%Variable InvertKSMaxIter
    !%Type integer
    !%Default 200
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects how many iterations of inversion will be done in the iterative scheme
    !%End
    call parse_variable('InvertKSMaxIter', 200, max_iter)  
           
    SAFE_ALLOCATE(vhxc(1:np, 1:nspin))

    vhxc(1:np,1:nspin) = aux_hm%vhxc(1:np,1:nspin)
         
    if(verbosity == 1 .or. verbosity == 2) then
      iunit = io_open('InvertKSconvergence', action = 'write')
    end if
    
    diffdensity = M_ONE
    counter = 0
    alpha = CNST(0.05)

    do while(diffdensity > convdensity .and. counter < max_iter)
      
      counter = counter + 1 
      beta = diffdensity*0.001 !parameter to avoid numerical problems due to small denominator
  
      if(verbosity == 2) then
        write(fname,'(i6.6)') counter
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "vhxc"//fname, gr%mesh, aux_hm%vhxc(:,1), units_out%energy, ierr)
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, gr%mesh, st%rho(:,1), units_out%length**(-gr%sb%dim), ierr)
      end if

      call hamiltonian_update(aux_hm, gr%mesh)
      call eigensolver_run(eigensolver, gr, st, aux_hm, 1)
      call density_calc(st, gr, st%rho)      

      ! Inversion according to Stella/Verstraete
      do ii = 1, nspin
        do jj = 1, np
          vhxc(jj, ii) = vhxc(jj, ii) &
                       + ((st%rho(jj, ii) - target_rho(jj, ii))/(target_rho(jj, ii) + beta))*alpha
        end do
      end do

      diffdensity = M_ZERO
      do jj = 1, nspin
        do ii = 1, np
          if (abs(st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
            diffdensity = abs(st%rho(ii,jj)-target_rho(ii,jj))
          end if
        end do
      end do
            
      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity
#ifdef HAVE_FLUSH
        call flush(iunit)
#endif
      end if

      aux_hm%vhxc(1:np,1:nspin) = vhxc(1:np, 1:nspin)
      
      do jj = 1, nspin
        aux_hm%vxc(:,jj) = vhxc(:,jj) - aux_hm%vhartree(:)
      end do
    end do

    !ensure correct asymptotic behavior, only for 1D potentials at the moment
    !need to find a way to find all points from where asymptotics should start in 2 and 3D
    if(asymptotics == XC_ASYMPTOTICS_SC) then
      do ii = 1, nspin
        do jj = 1, int(np/2)
          if(target_rho(jj,ii) < convdensity*CNST(10.0)) then
            call mesh_r(gr%mesh, jj, rr)
            vhxc(jj, ii) = (st%qtot-1.0)/sqrt(rr**2 + 1.0)
            asym1 = jj
          end if
          if(target_rho(np-jj+1, ii) < convdensity*CNST(10.0)) then
            asym2 = np - jj + 1
          end if
        end do
     
        ! calculate constant shift for correct asymptotics and shift accordingly
        call mesh_r(gr%mesh, asym1+1, rr)
        shift  = vhxc(asym1+1, ii) - (st%qtot-1.0)/sqrt(rr**2 + 1.0)
        do jj = asym1+1, asym2-1
          vhxc(jj,ii) = vhxc(jj, ii) - shift
        end do
  
        call mesh_r(gr%mesh, asym2-1, rr)
        shift  = vhxc(asym2-1, ii) - (st%qtot-1.0)/sqrt(rr**2 + 1.0)
        do jj = 1, asym2-1
          vhxc(jj,ii) = vhxc(jj, ii) - shift
        end do
        do jj = asym2, np
          call mesh_r(gr%mesh, jj, rr)
          vhxc(jj, ii) = (st%qtot-1.0)/sqrt(rr**2 + 1.0)
        end do
      end do
    end if

    aux_hm%vhxc(1:np,1:nspin) = vhxc(1:np, 1:nspin)
      
    do jj = 1, nspin
      aux_hm%vxc(:,jj) = vhxc(:,jj) - aux_hm%vhartree(:)
    end do

    !calculate final density

    call hamiltonian_update(aux_hm, gr%mesh)
    call eigensolver_run(eigensolver, gr, st, aux_hm, 1)
    call density_calc(st, gr, st%rho)
    
    write(message(1),'(a,I8)') "Iterative KS inversion, iterations needed:", counter
    call messages_info(1)
    
    call io_close(iunit)      

    SAFE_DEALLOCATE_A(vhxc)

    POP_SUB(invertks_iter)

  end subroutine invertks_iter

  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_calc(ks_inversion, gr, hm, st, vxc, time)
    type(xc_ks_inversion_t),  intent(inout) :: ks_inversion
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_t),      intent(in)    :: hm
    type(states_t),           intent(inout) :: st
    FLOAT,                    intent(inout) :: vxc(:,:) !< vxc(gr%mesh%np, st%d%nspin)
    FLOAT, optional,          intent(in)    :: time

    integer :: ii

    if(ks_inversion%level == XC_KS_INVERSION_NONE) return

    PUSH_SUB(X(xc_ks_inversion_calc))

    call density_calc(st, gr, st%rho)
    
    if(present(time)) then
      write(message(1),'(A,F18.12)') 'xc_ks_inversion_calc - time:', time
      call messages_info(1)
    end if

    ks_inversion%aux_hm%energy%intnvxc     = M_ZERO
    ks_inversion%aux_hm%energy%hartree     = M_ZERO
    ks_inversion%aux_hm%energy%exchange    = M_ZERO
    ks_inversion%aux_hm%energy%correlation = M_ZERO
    
    ks_inversion%aux_hm%vhartree = hm%vhartree
 
    if (present(time) .and. time > M_ZERO) then
!      write(*,*) 'debug 1'
      do ii = 1, st%d%nspin
        ks_inversion%aux_hm%vhxc(:,ii) = ks_inversion%vxc_previous_step(:,ii)
      end do
    else 
      do ii = 1, st%d%nspin
        ks_inversion%aux_hm%vxc(:,ii)  = M_ZERO !hm%ep%vpsl(:)
        ks_inversion%aux_hm%vhxc(:,ii) = ks_inversion%aux_hm%vhartree(:) + ks_inversion%aux_hm%vxc(:,ii)
      end do
    end if
    !ks_inversion%aux_hm%ep%vpsl(:)  = M_ZERO ! hm%ep%vpsl(:)
    ks_inversion%aux_hm%ep%vpsl(:)  = hm%ep%vpsl(:)

    ! compute ks inversion, vhxc contains total KS potential
    
    select case (ks_inversion%method)
    ! adiabatic ks inversion
    case(XC_INV_METHOD_TWO_PARTICLE)
      call invertks_2part(ks_inversion%aux_st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver, ks_inversion%asymp)
    case(XC_INV_METHOD_VS_ITER)
      call invertks_iter(st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver, ks_inversion%asymp)
    end select

    !subtract external and Hartree potentials, ATTENTION: subtracts true external potential not adiabatic one 
    
    do ii = 1, st%d%nspin
      ks_inversion%aux_hm%vxc(:,ii)  = ks_inversion%aux_hm%vhxc(:,ii) - hm%vhartree(:)
    end do

    vxc = ks_inversion%aux_hm%vxc
    
    if (present(time)) then
!      write(*,*) 'debug 2'
      ks_inversion%vxc_previous_step = ks_inversion%aux_hm%vhxc
!      write(*,*) 'debug 3'
    end if

    POP_SUB(X(xc_ks_inversion_calc))

  end subroutine xc_ks_inversion_calc


end module xc_ks_inversion_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
