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
!! $Id: run.F90 5113 2009-03-26 02:06:16Z dstrubbe $

#include "global.h"

module invert_ks_m
  use datasets_m
  use global_m
  use hamiltonian_m
  use h_sys_output_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use parser_m
  use poisson_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mix_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use unit_m
  use unit_system_m
  use v_ks_m
  use xc_m
  
  implicit none

  private
  public :: invert_ks_run

contains

  ! ---------------------------------------------------------
  subroutine invert_ks_run(sys, hm)
    type(system_t),              intent(inout) :: sys
    type(hamiltonian_t),         intent(inout) :: hm

    type(scf_t) :: scfv

    integer :: ii, jj, ierr, np, ndim, nspin, idiffmax
    integer :: verbosity
    FLOAT :: diffdensity
    integer :: invksmethod
    FLOAT, allocatable :: target_rho(:,:)
      
    call push_sub('invert_ks.invert_ks_run')

    call messages_devel_version("Kohn-Sham inversion")
    
    !%Variable InvertKSmethod
    !%Type integer
    !%Default iterative
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Selects whether the exact two-particle method or the iterative scheme
    !% is used to invert the density to get the KS potential.
    !%Option iterative 0
    !% Iterative scheme.
    !%Option two_particle 1
    !% Exact two-particle scheme.
    !%Option iterativevxc 2
    !% Iterative scheme for vxc
    !%End
      
    call parse_integer(datasets_check('InvertKSmethod'), 0, invksmethod)

    if(invksmethod < 0 .or. invksmethod > 2) then
      call input_error('InvertKSmethod')
      call write_fatal(1)
    endif


    ! initialize random wavefunctions
    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call states_generate_random(sys%st, sys%gr%mesh)

    !abbreviations
    np = sys%gr%mesh%np
    ndim    = sys%gr%mesh%sb%dim
    nspin   = sys%st%d%nspin

    ! read target density
    SAFE_ALLOCATE(target_rho(1:np, 1:nspin))    
    call read_target_rho()

    hm%epot     = M_ZERO
    hm%ehartree = M_ZERO
    hm%ex       = M_ZERO
    hm%ec       = M_ZERO
    hm%vhartree = M_ZERO
    hm%vxc      = M_ZERO
    !hm%ep%vpsl  = M_ZERO
    !hm%ep%local_potential = M_ZERO

    sys%ks%frozen_hxc = .true. !do not change hxc potential outside this routine
    
    call scf_init(scfv, sys%gr, sys%geo, sys%st, hm)
    
    if (invksmethod == 1) then ! 2-particle exact inversion

      call invertks_2part(target_rho, np, nspin, hm%vhxc, sys%gr%mesh%h)

    else ! iterative case
      if (invksmethod == 0) then ! iterative procedure for v_s (might be discontinued)
      
        call invertks_iter(target_rho, np, nspin, hm, sys, scfv)
      
      else
      
        call invertvxc_iter(target_rho, np, nspin, hm, sys, scfv)
      
      endif
    end if ! invksmethod

    ! this is debugging: should output quality of KS inversion
    call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, gs_run = .false., &
                 verbosity = VERB_COMPACT)
    call states_calc_dens(sys%st, sys%gr)

    diffdensity = 0d0
    do jj = 1, nspin
      do ii = 1, np
        if (abs(sys%st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
          diffdensity = abs(sys%st%rho(ii,jj)-target_rho(ii,jj))
          idiffmax=ii
        endif
      enddo
    enddo
    write (message(1),'(a,F16.6)') 'Achieved difference in densities wrt target:', &
        diffdensity
    call write_info(1)

    call scf_end(scfv)
    
    ! output for all cases    
    call h_sys_output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, STATIC_DIR)

    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vxc", sys%gr%mesh, hm%vxc(:,1), units_out%energy, ierr)
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "rho", sys%gr%mesh, sys%st%rho(:,1), units_out%length**(-sys%gr%sb%dim), ierr)

    SAFE_DEALLOCATE_A(target_rho)
    call pop_sub()
    
  contains

  ! ---------------------------------------------------------
    subroutine read_target_rho()
      character(len=256) :: filename
      integer :: pass, iunit, ierr, ii, npoints
      FLOAT   :: l_xx(MAX_DIM), l_ff(4), rr
      FLOAT, allocatable :: xx(:,:), ff(:,:)

      call push_sub('invert_ks.invert_ks_run.read_target_rho')

      !%Variable InvertKSTargetDensity
      !%Type string
      !%Default target_density.dat
      !%Section Calculation Modes::Invert KS
      !%Description
      !% Name of the file that contains the density used as the target in the 
      !% inversion of the KS equations.
      !%End
      call parse_string(datasets_check('InvertKSTargetDensity'), "target_density.dat", filename)

      iunit = io_open(filename, action='read', status='old')

      npoints = 0
      do pass = 1, 2
        ii = 0
        rewind(iunit)
        do
          read(iunit, fmt=*, iostat=ierr) l_xx(1:ndim), l_ff(1:nspin)
          if(ierr.ne.0) exit
          ii = ii + 1
          if(pass == 1) npoints = npoints + 1
          if(pass == 2) then
            xx(ii, 1:ndim)  = l_xx(1:ndim)
            ff(ii, 1:nspin) = l_ff(1:nspin)
          end if
        end do
        if(pass == 1) then
          SAFE_ALLOCATE(xx(1:npoints, 1:ndim))
          SAFE_ALLOCATE(ff(1:npoints, 1:nspin))
	  
        end if
      end do

      do ii = 1, nspin
        call dmf_interpolate_points(ndim, npoints, xx, ff(:,ii), &
             np, sys%gr%mesh%x, target_rho(:, ii))
      end do
 
      ! we now renormalize the density (necessary if we have a charged system)
      rr = M_ZERO
      do ii = 1, sys%st%d%spin_channels
        rr = rr + dmf_integrate(sys%gr%mesh, target_rho(:, ii))
      end do
      rr = sys%st%qtot/rr
      target_rho(:,:) = rr*target_rho(:,:)

      call doutput_function(io_function_fill_how("AxisX"), &
           ".", "func", sys%gr%mesh, target_rho(:,1), units_out%length**(-sys%gr%sb%dim), ierr)

      SAFE_DEALLOCATE_A(xx)
      SAFE_DEALLOCATE_A(ff)

      call pop_sub()
    end subroutine read_target_rho

  end subroutine invert_ks_run

  subroutine invertks_2part(target_rho, np, nspin, vhxc, spacing)
    integer, intent(in)  :: np, nspin
    FLOAT,   intent(in)  :: spacing(1:MAX_DIM)
    FLOAT,   intent(in)  :: target_rho(1:np, 1:nspin)
    FLOAT,   intent(out) :: vhxc(1:np, 1:nspin)
       
    integer :: ii, jj
    FLOAT, allocatable :: sqrtrho(:,:), gradrho(:,:)

    call push_sub('invert_ks.invertks_2part')
    
    !initialize the KS potential
    SAFE_ALLOCATE(sqrtrho(1:np, 1:nspin))
    SAFE_ALLOCATE(gradrho(1:np, 1:nspin))

    gradrho = 0d0
    do jj = 1, nspin
      do ii = 1, np
        sqrtrho(ii, jj) = sqrt(target_rho(ii, jj))
      enddo
      do ii = 3, np-2
        gradrho(ii, jj) = - sqrtrho(ii+2, jj) - sqrtrho(ii-2, jj) &
                        & + 16d0*(sqrtrho(ii+1, jj) + sqrtrho(ii-1, jj)) &
                        & - 30d0*sqrtrho(ii, jj)
      enddo
      
      gradrho(:, jj) = gradrho(:, jj)/(12d0*spacing(1)**2)
      do ii = 1, np
        vhxc(ii, jj) = gradrho(ii, jj)/(2d0*sqrtrho(ii, jj))
      enddo
    enddo
    SAFE_DEALLOCATE_A(sqrtrho)
    SAFE_DEALLOCATE_A(gradrho)

    call pop_sub()
  end subroutine invertks_2part


  subroutine invertks_iter(target_rho, np, nspin, hm, sys, scfv)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    type(scf_t),         intent(inout) :: scfv
    integer,             intent(in)    :: np, nspin
    FLOAT,               intent(in)    :: target_rho(1:np, 1:nspin)
        
    integer :: ii, jj, ierr, idiffmax
    integer :: iunit, verbosity, counter
    FLOAT :: stabilizer, convdensity, diffdensity
    FLOAT, allocatable :: vhxc_in(:,:,:), vhxc_out(:,:,:), vhxc_mix(:,:,:)
    type(mix_t) :: smix
    character(len=256) :: fname

    call push_sub('invert_ks.invertks_iter')
    
    !%Variable InvertKSConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Absolute difference between the calculated and the target density in the KS
    !% inversion. Has to be larger than the convergence of the density in the SCF run.
    !%End
    
    call parse_float(datasets_check('InvertKSConvAbsDens'), 1d-5, convdensity)
    
    !%Variable InvertKSStabilizer
    !%Type float
    !%Default 0.5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Additive constant c in the iterative calculation of the KS potential
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
      call write_fatal(1)
    endif
  
    
    !initialize the KS potential
    hm%vhxc = 1d0
       
    call mix_init(smix, np, nspin, 1, prefix_="InvertKS")

    SAFE_ALLOCATE(vhxc_in(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vhxc_out(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vhxc_mix(1:np, 1:nspin, 1:1))
    
    vhxc_in(1:np,1:nspin,1) = hm%vhxc(1:np,1:nspin)
         
    
    if(verbosity == 1 .or. verbosity == 2) then
      iunit = io_open('InvertKSconvergence', action = 'write')
    endif
    
    diffdensity = 1d0
    counter = 0
        
    do while(diffdensity > convdensity)
      
      counter = counter + 1 
        
      if(verbosity == 2) then
        write(fname,'(i6.6)') counter
        call doutput_function(io_function_fill_how("AxisX"), &
             ".", "vhxc"//fname, sys%gr%mesh, hm%vhxc(:,1), units_out%energy, ierr)
        call doutput_function(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, sys%gr%mesh, sys%st%rho(:,1), units_out%length**(-sys%gr%sb%dim), ierr)
      endif
    
      call hamiltonian_update_potential(hm, sys%gr%mesh)

      call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, gs_run = .false., &
                   verbosity = VERB_COMPACT)
      call states_calc_dens(sys%st, sys%gr)
      
      vhxc_out(1:np, 1:nspin, 1) = &
        (sys%st%rho(1:np,1:nspin) + stabilizer)/&
	(target_rho(1:np,1:nspin) + stabilizer) &
         * hm%vhxc(1:np, 1:nspin)

      call dmixing(smix, counter, vhxc_in, vhxc_out, vhxc_mix, dmf_dotp_aux)

      hm%vhxc(1:np,1:nspin) = vhxc_mix(1:np, 1:nspin, 1)
      vhxc_in(1:np, 1:nspin, 1) = hm%vhxc(1:np, 1:nspin)
      
      diffdensity = 0d0
      do jj = 1, nspin
        do ii = 1, np
          if (abs(sys%st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
            diffdensity = abs(sys%st%rho(ii,jj)-target_rho(ii,jj))
            idiffmax=ii
          endif
        enddo
      enddo
            
      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity
#ifdef HAVE_FLUSH
        call flush(iunit)
#endif
      endif
     
    end do

    write(message(1),'(a,I8)') "Invert KS: iterations needed:", counter
    call write_info(1)
    
    call io_close(iunit)      
    
    call mix_end(smix)

    SAFE_DEALLOCATE_A(vhxc_in)
    SAFE_DEALLOCATE_A(vhxc_out)
    SAFE_DEALLOCATE_A(vhxc_mix)

    call pop_sub()

  end subroutine invertks_iter

  subroutine invertvxc_iter(target_rho, np, nspin, hm, sys, scfv)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: hm
    type(scf_t),         intent(inout) :: scfv
    integer,             intent(in)    :: np, nspin
    FLOAT,               intent(in)    :: target_rho(1:np, 1:nspin)
        
    integer :: ii, jj, ierr, idiffmax
    integer :: iunit, verbosity, counter
    FLOAT :: stabilizer, convdensity, diffdensity
    FLOAT :: E_x, E_c
    FLOAT, allocatable :: vxc_in(:,:,:), vxc_out(:,:,:), vxc_mix(:,:,:)
    FLOAT, allocatable :: rho(:)
    type(mix_t) :: smix
    character(len=256) :: fname

    call push_sub('invert_ks.invertvxc_iter')
    
    !%Variable InvertKSConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Absolute difference between the calculated and the target density in the KS
    !% inversion. Has to be larger than the convergence of the density in the SCF run.
    !%End
    
    call parse_float(datasets_check('InvertKSConvAbsDens'), 1d-5, convdensity)
    
    !%Variable InvertKSStabilizer
    !%Type float
    !%Default 0.5
    !%Section Calculation Modes::Invert KS
    !%Description
    !% Additive constant c in the iterative calculation of the KS potential
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
      call write_fatal(1)
    endif
  
    SAFE_ALLOCATE(rho(1:np))
    
    !calculate total density
    
    rho = M_ZERO
    
    do ii = 1, nspin
      do jj = 1, np
        rho(jj) = rho(jj) + target_rho(jj, ii)
      enddo
    enddo
    
    !calculate the Hartree potential
    
    call dpoisson_solve(sys%ks%hartree_solver,hm%vhartree,rho)
    
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vhartree", sys%gr%mesh, hm%vhartree(:), units_out%energy, ierr)
    
    !initialize the KS potential
    
    call xc_get_vxc(sys%gr, sys%ks%xc, sys%st, target_rho, sys%st%d%ispin, E_x, E_c, &
      M_ZERO, sys%st%qtot, hm%vxc)
    
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vxc", sys%gr%mesh, hm%vxc(:,1), units_out%energy, ierr)
    
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vext", sys%gr%mesh, hm%ep%vpsl(:), units_out%energy, ierr)
    
    call hamiltonian_update_potential(hm, sys%gr%mesh)
    
    sys%ks%frozen_hxc = .true.
    
    call mix_init(smix, np, nspin, 1, prefix_="InvertKS")

    SAFE_ALLOCATE(vxc_in(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vxc_out(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vxc_mix(1:np, 1:nspin, 1:1))
    
    vxc_in(1:np,1:nspin,1) = hm%vxc(1:np,1:nspin)
         
    if(verbosity == 1 .or. verbosity == 2) then
      iunit = io_open('InvertKSconvergence', action = 'write')
    endif
    
    diffdensity = 1d0
    counter = 0
        
    do while(diffdensity > convdensity)
  ! do while (counter <5)   
      counter = counter + 1 
        
      if(verbosity == 2) then
        write(fname,'(i6.6)') counter
        call doutput_function(io_function_fill_how("AxisX"), &
             ".", "vxc"//fname, sys%gr%mesh, hm%vxc(:,1), units_out%energy, ierr)
        call doutput_function(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, sys%gr%mesh, sys%st%rho(:,1), units_out%length**(-sys%gr%sb%dim), ierr)
      endif
    
      call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, gs_run = .false., &
                   verbosity = VERB_COMPACT)
      call states_calc_dens(sys%st, sys%gr)
      
      vxc_out(1:np, 1:nspin, 1) = &
        (sys%st%rho(1:np,1:nspin) + stabilizer)/&
	(target_rho(1:np,1:nspin) + stabilizer) &
         * hm%vxc(1:np, 1:nspin)

      call dmixing(smix, counter, vxc_in, vxc_out, vxc_mix, dmf_dotp_aux)

      hm%vxc(1:np,1:nspin) = vxc_mix(1:np, 1:nspin, 1)
      vxc_in(1:np, 1:nspin, 1) = hm%vxc(1:np, 1:nspin)
      
      call hamiltonian_update_potential(hm, sys%gr%mesh)
      
      diffdensity = 0d0
      do jj = 1, nspin
        do ii = 1, np
          if (abs(sys%st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
            diffdensity = abs(sys%st%rho(ii,jj)-target_rho(ii,jj))
            idiffmax=ii
          endif
        enddo
      enddo
            
      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity

#ifdef HAVE_FLUSH
        call flush(iunit)
#endif
      endif
     
    end do

    write(message(1),'(a,I8)') "Invert KS: iterations needed:", counter
    call write_info(1)
    
    call io_close(iunit)      
    
    call mix_end(smix)

    SAFE_DEALLOCATE_A(vxc_in)
    SAFE_DEALLOCATE_A(vxc_out)
    SAFE_DEALLOCATE_A(vxc_mix)
    SAFE_DEALLOCATE_A(rho)
    
    call pop_sub()

  end subroutine invertvxc_iter


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
    
    call push_sub('invert_ks.precond_kiks')

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
      inverse = 1d0/eigenvals(ivec)
      do ip = 1, np
        matrixmul(ip, ivec) = ki(ip,ivec)*inverse
	kired(ip, ivec) = ki(ip, ivec)
      enddo
    enddo
    
    inverseki = matmul(matrixmul, transpose(kired))
    
    vhxc_out = M_ZERO
    
    vol_element = 1.0d0
    do jdim = 1, MAX_DIM
      if (mesh%h(jdim) > 1.e-10) vol_element = vol_element*mesh%h(jdim)
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

    call pop_sub()
  
  end subroutine precond_kiks
  
end module invert_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
