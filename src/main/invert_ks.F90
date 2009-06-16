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
  use loct_parser_m
  use mesh_function_m
  use messages_m
  use mix_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use v_ks_m
  
  implicit none

  private
  public :: invert_ks_run

contains

  ! ---------------------------------------------------------
  subroutine invert_ks_run(sys, hm)
    type(system_t),              intent(inout) :: sys
    type(hamiltonian_t),         intent(inout) :: hm

    integer :: ii, jj, kk, ierr, np, ndim, nspin, counter, idiffmax
    integer :: iunit, verbosity
    FLOAT :: diffdensity, mixing, stabilizer, convdensity
    FLOAT, allocatable :: target_rho(:,:), vhxc_in(:,:,:), vhxc_out(:,:,:), vhxc_mix(:,:,:)
    type(scf_t) :: scfv
    type(mix_t) :: smix
    character(len=256) :: fname

    ! initialize random wave-functions
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
    hm%ep%vpsl  = M_ZERO
    !hm%ep%local_potential = M_ZERO

    sys%ks%frozen_hxc = .true. !do not change hxc potential outside this routine
    
    !initialize the ks potential
    hm%vhxc = 1d0
     
    call scf_init(scfv, sys%gr, sys%geo, sys%st, hm)
    call mix_init(smix, np, nspin, 1, prefix_="InvertKS")

    SAFE_ALLOCATE(vhxc_in(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vhxc_out(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(vhxc_mix(1:np, 1:nspin, 1:1))

    vhxc_in(1:np,1:nspin,1) = hm%vhxc(1:np,1:nspin)
    
    !%Variable InvertKSConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section Main:: Invert KS
    !%Description
    !% Absolute difference between the calculated and the target density in the KS
    !% inversion, has to be larger than the convergence of the density in the scf run
    !%End
    
    call loct_parse_float(datasets_check('InvertKSConvAbsDens'), 1d-5, convdensity)
    
    !%Variable InvertKSstabilizer
    !%Type float
    !%Default 0.5
    !%Section Main:: Invert KS
    !%Description
    !% Additive constant c in the iterative calculation of the KS potential
    !%   (v(alpha+1)=rho(alpha)+c)/(rho_target+c)*v(alpha)
    !% ensures that very small densities do not cause numerical problems
    !%End

    call loct_parse_float(datasets_check('InvertKSstabilizer'), M_HALF, stabilizer)

    !%Variable InvertKSverbosity
    !%Type integer
    !%Default 0
    !%Section Main:: Invert KS
    !%Description
    !% Selects what is output during the calculation of the KS potential
    !%Option 0
    !% Only outputs the converged density and KS potential
    !%Option 1
    !% Same as 0 but outputs the maximum difference to the target density in each
    !% iteration in addition
    !%Option 2
    !% Same as 1 but outputs the density and the KS potential in each iteration in 
    !% addition  
    !%End
      
    call loct_parse_int(datasets_check('InvertKSverbosity'), 0, verbosity)  

    if(verbosity < 0 .or. verbosity > 2) then
      message(1) = 'InvertKSverbosity only has the options 0, 1, or 2'
      call write_fatal(1)
    endif
    
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
             ".", "vhxc"//fname, sys%gr%mesh, sys%gr%sb, hm%vhxc(:,1), M_ONE, ierr)
        call doutput_function(io_function_fill_how("AxisX"), &
             ".", "rho"//fname, sys%gr%mesh, sys%gr%sb, sys%st%rho(:,1), M_ONE, ierr)
      endif
    
      call hamiltonian_update_potential(hm, sys%gr%mesh)

      call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, gs_run = .false.)
      call states_calc_dens(sys%st, np)
            
      vhxc_out(1:np, 1:nspin, 1) = &
      (sys%st%rho(1:np,1:nspin) + stabilizer)/(target_rho(1:np,1:nspin) + stabilizer) &
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
      
      !write(*,*) diffdensity
      
      if(verbosity == 1 .or. verbosity == 2) then
        write(iunit,'(i6.6)', ADVANCE = 'no') counter
        write(iunit,'(es18.10)') diffdensity
#ifdef HAVE_FLUSH
        call flush(iunit)
#endif
      endif
     
    end do
    
    call io_close(iunit)      
    
    call h_sys_output_all(sys%outp, sys%gr, sys%geo, sys%st, hm, STATIC_DIR)
    
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vhxc", sys%gr%mesh, sys%gr%sb, hm%vhxc(:,1), M_ONE, ierr)
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "rho", sys%gr%mesh, sys%gr%sb, sys%st%rho(:,1), M_ONE, ierr)

    write(*,*) "iterations needed:", counter

    call mix_end(smix)
    call scf_end(scfv)

    SAFE_DEALLOCATE_A(target_rho)
    SAFE_DEALLOCATE_A(vhxc_in)
    SAFE_DEALLOCATE_A(vhxc_out)
    SAFE_DEALLOCATE_A(vhxc_mix)
    
  contains

  ! ---------------------------------------------------------
    subroutine read_target_rho()
      character(len=256) :: filename
      integer :: pass, iunit, ierr, ii, npoints
      FLOAT   :: l_xx(MAX_DIM), l_ff(4), rr
      FLOAT, allocatable :: xx(:,:), ff(:,:)

      call loct_parse_string(datasets_check('InvertKSTargetDensity'), "target_density.dat", filename)

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
           ".", "func", sys%gr%mesh, sys%gr%sb, target_rho(:,1), M_ONE, ierr)

      SAFE_DEALLOCATE_A(xx)
      SAFE_DEALLOCATE_A(ff)

    end subroutine read_target_rho

  end subroutine invert_ks_run

end module invert_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
