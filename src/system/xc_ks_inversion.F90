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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: $

#include "global.h"

module xc_ks_inversion_m
  use datasets_m
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
  use mix_m 
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
    xc_ks_inversion_messages_info,    &
    xc_ks_inversion_calc,          &
    invertks_2part,                &
    invertks_iter

  ! KS inversion methods/algorithms
  integer, public, parameter ::            &
    XC_INV_METHOD_VS_ITER      = 1,  &
    XC_INV_METHOD_TWO_PARTICLE = 2

  ! the KS inversion levels
  integer, public, parameter ::      &
    XC_KS_INVERSION_NONE      = 1,   &
    XC_KS_INVERSION_ADIABATIC = 2,   &   
    XC_KS_INVERSION_TD_EXACT  = 3

  type xc_ks_inversion_t
     integer             :: method
     integer             :: level
     FLOAT, pointer      :: vxc_previous_step(:,:)
     type(states_t)      :: aux_st
     type(hamiltonian_t) :: aux_hm
     type(eigensolver_t) :: eigensolver
  end type xc_ks_inversion_t


contains

  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_init(ks_inv, family, gr, geo, d, mc)
    type(xc_ks_inversion_t), intent(out)   :: ks_inv
    integer,                 intent(in)    :: family
    type(grid_t),            intent(inout) :: gr
    type(states_dim_t),      intent(in)    :: d
    type(geometry_t),        intent(inout) :: geo
    type(multicomm_t),       intent(in)    :: mc  

    PUSH_SUB(xc_ks_inversion_init)

    if(iand(family, XC_FAMILY_KS_INVERSION) .eq. 0) then
      ks_inv%level = XC_KS_INVERSION_NONE
      POP_SUB(xc_ks_inversion_init)
      return
    end if

#if defined(HAVE_MPI)
    message(1) = "KS Inversion currently not available in parallel. Stopping octopus."
    call messages_fatal(1)
#endif

    call messages_experimental("Kohn-Sham inversion")
    
    !%Variable InvertKSmethod
    !%Type integer
    !%Default iterative
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects whether the exact two-particle method or the iterative scheme
    !% is used to invert the density to get the KS potential.
    !%Option iterative 1
    !% Iterative scheme for v_s.
    !%Option two_particle 2
    !% Exact two-particle scheme.
    !%Option iterativevxc 3
    !% Iterative scheme for v_xc.
    !%End
    call parse_integer(datasets_check('InvertKSmethod'), &
            XC_INV_METHOD_VS_ITER, ks_inv%method)

    if(ks_inv%method < XC_INV_METHOD_VS_ITER &
      .or. ks_inv%method > XC_INV_METHOD_TWO_PARTICLE) then
      call input_error('InvertKSmethod')
      call messages_fatal(1)
    endif

    !%Variable KS_Inversion_Level
    !%Type integer
    !%Default ks_inversion_adiabatic
    !%Section Hamiltonian::XC
    !%Description
    !% At what level shall <tt>Octopus</tt> handle the KS inversion
    !%Option ks_inversion_none 1
    !% Do not compute KS inversion
    !%Option ks_inversion_adiabatic 2
    !% Compute exact adiabatic vxc
    !%End
    
    call parse_integer(datasets_check('KS_Inversion_Level'), XC_KS_INVERSION_ADIABATIC, ks_inv%level)
    if(.not.varinfo_valid_option('KS_Inversion_Level', ks_inv%level)) call input_error('KS_Inversion_Level')

    if(ks_inv%level.ne.XC_KS_INVERSION_NONE) then
      ! initialize auxilary random wavefunctions
      call states_null(ks_inv%aux_st)
      call states_init(ks_inv%aux_st, gr, geo)      
      call states_allocate_wfns(ks_inv%aux_st, gr%mesh)
      call states_generate_random(ks_inv%aux_st, gr%mesh)      
      ! initialize densities, hamiltonian and eigensolver
      call states_densities_init(ks_inv%aux_st, gr, geo, mc)
      call hamiltonian_init(ks_inv%aux_hm, gr, geo, ks_inv%aux_st, INDEPENDENT_PARTICLES, XC_FAMILY_NONE)
      call eigensolver_init(ks_inv%eigensolver, gr, ks_inv%aux_st)
    end if

    SAFE_ALLOCATE(ks_inv%vxc_previous_step(1:gr%mesh%np, 1:d%nspin))
    POP_SUB(xc_ks_inversion_init)
  end subroutine xc_ks_inversion_init


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_end(ks_inv, gr, geo)
    type(xc_ks_inversion_t), intent(inout) :: ks_inv
    type(grid_t),            intent(inout) :: gr
    type(geometry_t),        intent(inout) :: geo

    PUSH_SUB(xc_ks_inversion_end)

    if(ks_inv%level .ne. XC_KS_INVERSION_NONE) then
      ! cleanup
      call eigensolver_end(ks_inv%eigensolver)
      !call hamiltonian_end(ks_inv%aux_hm, gr, geo)
      call states_end(ks_inv%aux_st)
    end if

    SAFE_DEALLOCATE_P(ks_inv%vxc_previous_step)
    POP_SUB(xc_ks_inversion_end)
  end subroutine xc_ks_inversion_end


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_messages_info(ks_inversion, iunit)
    type(xc_ks_inversion_t), intent(in) :: ks_inversion
    integer,                 intent(in) :: iunit

    if(ks_inversion%level.eq.XC_KS_INVERSION_NONE) return

    PUSH_SUB(xc_ks_inversion_messages_info)
    call messages_print_var_option(iunit, 'KS_Inversion_Level', ks_inversion%level)

    POP_SUB(xc_ks_inversion_messages_info)
  end subroutine xc_ks_inversion_messages_info


  ! ---------------------------------------------------------
  subroutine invertks_2part(target_rho, nspin, aux_hm, gr, st, eigensolver)
    
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(eigensolver_t), intent(inout) :: eigensolver
    integer, intent(in)      :: nspin
    FLOAT,   intent(in)      :: target_rho(1:gr%mesh%np, 1:nspin)
           
    integer :: ii, jj 
    integer :: np_part, ndim, np
    FLOAT   :: spacing(1:MAX_DIM), stabilizer
    FLOAT, allocatable :: sqrtrho(:,:), laplace(:,:), vks(:,:)

    PUSH_SUB(invertks_2part)
    
    call parse_float(datasets_check('InvertKSStabilizer'), M_HALF, stabilizer)
    
    ndim = gr%sb%dim
    spacing = gr%mesh%spacing
    np = gr%mesh%np
    
    SAFE_ALLOCATE(sqrtrho(1:gr%der%mesh%np_part, 1:nspin))
    SAFE_ALLOCATE(vks(1:np, 1:nspin))
    SAFE_ALLOCATE(laplace(1:gr%der%mesh%np, 1:nspin))
    
    sqrtrho = M_ZERO
    
    do jj = 1, nspin
      do ii = 1, gr%der%mesh%np
        sqrtrho(ii, jj) = sqrt(target_rho(ii, jj))
      enddo
    enddo   
    
    do jj = 1, nspin
      call dderivatives_lapl(gr%der, sqrtrho(:,jj), laplace(:,jj))
    enddo
    
    do jj = 1, nspin
      do ii = 1, np
          vks(ii, jj) = laplace(ii, jj)/(M_TWO*sqrtrho(ii, jj)) + st%eigenval(1,jj)
      enddo
    enddo
    
    do jj = 1, nspin 
      aux_hm%vxc(:,jj) = vks(:,jj) - aux_hm%ep%vpsl(:) - aux_hm%vhartree(:)
      aux_hm%vhxc(:,jj) = aux_hm%vxc(:,jj) + aux_hm%vhartree(:)
    enddo
    
    call hamiltonian_update(aux_hm, gr%mesh)
    
    call eigensolver_run(eigensolver, gr, st, aux_hm, 1, verbose = .false.)

    call density_calc(st, gr, st%rho)
    
    SAFE_DEALLOCATE_A(sqrtrho)
    SAFE_DEALLOCATE_A(laplace)
    SAFE_DEALLOCATE_A(vks)

    POP_SUB(invertks_2part)
  end subroutine invertks_2part


  ! ---------------------------------------------------------
  subroutine invertks_iter(target_rho, nspin, aux_hm, gr, st, eigensolver)
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(inout) :: st
    type(hamiltonian_t), intent(inout) :: aux_hm
    type(eigensolver_t), intent(inout) :: eigensolver
    integer,             intent(in)    :: nspin
    FLOAT,               intent(in)    :: target_rho(1:gr%mesh%np, 1:nspin)
        
    type(mix_t) :: smix
    integer :: ii, jj, ierr, idiffmax
    integer :: iunit, iunit2, verbosity, counter, np
    FLOAT :: rr
    FLOAT :: alpha, beta, convergence, alphascale, betascale
    FLOAT :: stabilizer, convdensity, diffdensity, aa
    FLOAT, allocatable :: vhxc(:,:)

    character(len=20) :: alpha_str
    character(len=20) :: beta_str
    character(len=256) :: fname
    character(len=256) :: filename2

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
    call parse_float(datasets_check('InvertKSConvAbsDens'), CNST(1e-5), convdensity)

    !%Variable InvertKSStabilizer
    !%Type float
    !%Default 0.5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Additive constant <i>c</i> in the iterative calculation of the KS potential
    !%   (v(alpha+1)=rho(alpha)+c)/(rho_target+c)*v(alpha)
    !% ensures that very small densities do not cause numerical problems.
    !%End
    call parse_float(datasets_check('InvertKSStabilizer'), M_HALF, stabilizer)

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
    call parse_integer(datasets_check('InvertKSVerbosity'), 0, verbosity)  
    if(verbosity < 0 .or. verbosity > 2) then
      call input_error('InvertKSVerbosity')
      call messages_fatal(1)
    endif
           
    SAFE_ALLOCATE(vhxc(1:np, 1:nspin))

    vhxc(1:np,1:nspin) = aux_hm%vhxc(1:np,1:nspin)
         
    if(verbosity == 1 .or. verbosity == 2) then
      iunit = io_open('InvertKSconvergence', action = 'write')
    endif
    
    diffdensity = M_ONE
    counter = 0
    alpha = CNST(0.1)
    beta  = CNST(0.1)
    convergence = CNST(0.1)
    alphascale = 10
    betascale = 100
 
    if(verbosity == 2) then
      write(alpha_str, '(f8.4)') alpha
      write(beta_str, '(f8.4)') beta
      filename2 = "diffdens_" // trim(adjustl(alpha_str)) // "_" // trim(adjustl(beta_str)) // ".dat"
      iunit2 = io_open(trim(filename2), action='write')
    end if

    do while(diffdensity > convdensity)
      
      counter = counter + 1 
        
      if(verbosity == 2) then
        write(fname,'(i6.6)') counter
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "vhxc"//fname, gr%mesh, aux_hm%vhxc(:,1), units_out%energy, ierr)
        call dio_function_output(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, gr%mesh, st%rho(:,1), units_out%length**(-gr%sb%dim), ierr)
      endif

      call hamiltonian_update(aux_hm, gr%mesh)
      call eigensolver_run(eigensolver, gr, st, aux_hm, 1, verbose = .false.)
      call density_calc(st, gr, st%rho)      
       
      ! Inversion according to Phys. Rev. Lett. 100, 153004 (2008), Eq. (6)
      do jj = 1, np
        call mesh_r(gr%mesh, jj, rr)
        if (abs(rr).lt.CNST(1e-1)) rr = CNST(1e-1) 
        vhxc(jj, 1:nspin) = vhxc(jj, 1:nspin) + (st%rho(jj,1:nspin) - target_rho(jj,1:nspin))*alpha*rr**beta
      end do

      aux_hm%vhxc(1:np,1:nspin) = vhxc(1:np, 1:nspin)
      do jj = 1, nspin
        aux_hm%vxc(:,jj) = vhxc(:,jj) - aux_hm%vhartree(:)
      enddo
      
      diffdensity = M_ZERO
      do jj = 1, nspin
        do ii = 1, np
          if (abs(st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
            diffdensity = abs(st%rho(ii,jj)-target_rho(ii,jj))
            idiffmax=ii
          endif
        enddo
      enddo

      if (diffdensity .lt. convergence) then
        if(verbosity == 2) then
          write(iunit2,*) counter, diffdensity
        end if
        convergence = convergence / CNST(10)
      endif
      beta = min(M_ONE, diffdensity*betascale)
      alpha = M_ONE - diffdensity
            
      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity
#ifdef HAVE_FLUSH
        call flush(iunit)
#endif
      endif
    end do

    ! calculate final density
    call hamiltonian_update(aux_hm, gr%mesh)
    call eigensolver_run(eigensolver, gr, st, aux_hm, 1, verbose = .false.)
    call density_calc(st, gr, st%rho)
    
    write(message(1),'(a,I8)') "Iterative KS inversion, iterations needed:", counter
    call messages_info(1)
    
    call io_close(iunit)      
    if(verbosity == 2) then
      call io_close(iunit2)      
    end if

    SAFE_DEALLOCATE_A(vhxc)

    POP_SUB(invertks_iter)

  end subroutine invertks_iter


  ! ---------------------------------------------------------
  subroutine precond_kiks(mesh, np, nspin, st, target_rho, vhxc_out)
    type(mesh_t),   intent(in)  :: mesh
    integer,        intent(in)  :: np, nspin
    type(states_t), intent(in)  :: st
    FLOAT,          intent(in)  :: target_rho(1:np, 1:nspin)
    FLOAT,          intent(out) :: vhxc_out(1:np, 1:nspin,1:1)
    
    integer :: ip, iprime, ii, jj, ivec, jdim
    FLOAT :: numerator, diffrho, epsij, occij, inverse
    FLOAT :: vol_element
    FLOAT :: ki(1:np, 1:np)
    FLOAT :: eigenvals(1:np), inverseki(1:np,1:np)
    FLOAT, allocatable :: matrixmul(:,:), kired(:,:)
    
    PUSH_SUB(precond_kiks)

    numerator = M_ZERO
    vhxc_out = M_ZERO
    
    !do ip = 1, np
    !  diffrho = st%rho(ip, 1) - target_rho(ip, 1)
    !  numerator = numerator + diffrho**2
    !end do
    
    ki = M_ZERO
    
    do jj = 1, st%nst
      do ii = jj + 1, st%nst
        epsij = M_ONE / (st%eigenval(jj, 1) - st%eigenval(ii, 1))
	occij = st%occ(jj, 1) - st%occ(ii, 1)
        do iprime = 1, np
          do ip = 1, np
            ki(ip, iprime) = ki(ip, iprime) + occij*epsij & 
	                    * (st%dpsi(ip, 1, ii, 1)*st%dpsi(ip, 1, jj, 1)) & 
                            * (st%dpsi(iprime, 1, ii, 1)*st%dpsi(iprime, 1, jj, 1))
	  enddo
	enddo
      end do
    end do
    
    call lalg_eigensolve(np, ki, eigenvals)
    
    !do ip = 1, np
    !  if(abs(eigenvals(ip))>1d-10) then
    !    neigenval = ip
    !  endif
    !enddo
    
    SAFE_ALLOCATE(matrixmul(1:np, 1:st%nst))
    SAFE_ALLOCATE(kired(1:np, 1:st%nst))
    
    do ivec = 1, st%nst
      inverse = M_ONE/eigenvals(ivec)
      do ip = 1, np
        matrixmul(ip, ivec) = ki(ip,ivec)*inverse
	kired(ip, ivec) = ki(ip, ivec)
      enddo
    enddo
    
    inverseki = matmul(matrixmul, transpose(kired))
    
    vhxc_out = M_ZERO
    
    vol_element = M_ONE
    do jdim = 1, MAX_DIM
      if (mesh%spacing(jdim) > 1.e-10) vol_element = vol_element*mesh%spacing(jdim)
    end do
    
    do iprime = 1, np
      diffrho = target_rho(iprime, 1) - st%rho(iprime, 1)
      do ip = 1, np
	vhxc_out(ip, 1, 1) = vhxc_out(ip, 1, 1) + inverseki(ip, iprime)*diffrho
	write(200,*) ip, iprime, inverseki(ip, iprime) 
      enddo
    enddo
   
    do ip = 1, np
      write(100,*) ip, vhxc_out(ip, 1, 1),  target_rho(ip, 1) - st%rho(ip, 1)
    enddo   
    
    
    SAFE_DEALLOCATE_A(matrixmul)
    SAFE_DEALLOCATE_A(kired)
        
#ifdef HAVE_FLUSH
    !call flush(200)
#endif

    POP_SUB(precond_kiks)
  
  end subroutine precond_kiks


  ! ---------------------------------------------------------
  subroutine xc_ks_inversion_calc(ks_inversion, gr, hm, st, ex, ec, vxc, time)
    use xc_functl_m

    type(xc_ks_inversion_t),  intent(inout) :: ks_inversion
    type(grid_t),             intent(inout) :: gr
    type(hamiltonian_t),      intent(in)    :: hm
    type(states_t),           intent(inout) :: st
    FLOAT,                    intent(inout) :: ex, ec
    FLOAT,                    intent(inout) :: vxc(:,:) ! vxc(gr%mesh%np, st%d%nspin)
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
      enddo
    else 
      do ii = 1, st%d%nspin
        ks_inversion%aux_hm%vxc(:,ii)  = M_ZERO !hm%ep%vpsl(:)
        ks_inversion%aux_hm%vhxc(:,ii) = ks_inversion%aux_hm%vhartree(:) + ks_inversion%aux_hm%vxc(:,ii)
      enddo
    end if
    !ks_inversion%aux_hm%ep%vpsl(:)  = M_ZERO ! hm%ep%vpsl(:)
    ks_inversion%aux_hm%ep%vpsl(:)  = hm%ep%vpsl(:)

    ! compute ks inversion, vhxc contains total KS potential
    
    select case (ks_inversion%method)
    ! adiabatic ks inversion
    case(XC_INV_METHOD_TWO_PARTICLE)
      call invertks_2part(ks_inversion%aux_st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver)
    case(XC_INV_METHOD_VS_ITER)
      call invertks_iter(st%rho, st%d%nspin, ks_inversion%aux_hm, gr, &
                         ks_inversion%aux_st, ks_inversion%eigensolver)
    end select

    !subtract external and Hartree potentials, ATTENTION: subtracts true external potential not adiabatic one 
    
    do ii = 1, st%d%nspin
      ks_inversion%aux_hm%vxc(:,ii)  = ks_inversion%aux_hm%vhxc(:,ii) - hm%vhartree(:)
    enddo

    hm%vxc = ks_inversion%aux_hm%vxc
    vxc = ks_inversion%aux_hm%vxc
    
    !do ii = 1, st%d%nspin
    !  hm%vhxc(:,ii) = hm%vhartree(:) + hm%vxc(:,ii)
    !enddo

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
