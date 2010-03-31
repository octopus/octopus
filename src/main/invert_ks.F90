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
  use eigensolver_m 
  use global_m
  use hamiltonian_m 
  use h_sys_output_m 
  use io_function_m 
  use io_m 
  use mesh_function_m 
  use messages_m 
  use parser_m 
  use profiling_m 
  use states_m 
  use system_m 
  use unit_m 
  use unit_system_m 
  use xc_ks_inversion_m 
  
  implicit none

  private
  public :: invert_ks_run

contains

  ! ---------------------------------------------------------
  subroutine invert_ks_run(sys, hm)
    type(system_t),              intent(inout) :: sys
    type(hamiltonian_t),         intent(inout) :: hm

    integer :: ii, jj, ierr, np, ndim, nspin, idiffmax
    integer :: verbosity
    FLOAT   :: diffdensity
    FLOAT, allocatable :: target_rho(:,:)
      
    call push_sub('invert_ks.invert_ks_run')

    ! initialize KS inversion module
    call xc_ks_inversion_init(sys%ks%ks_inversion, sys%ks%xc_family, &
         sys%gr, sys%geo, sys%st%d, sys%mc)
    call xc_ks_inversion_write_info(sys%ks%ks_inversion, stdout)

    !abbreviations
    np      = sys%gr%mesh%np
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

    ! do not change hxc potential outside this routine
    sys%ks%frozen_hxc = .true. 
    
    if (sys%ks%ks_inversion%method == XC_INVERSION_METHOD_TWO_PARTICLE) then ! 2-particle exact inversion
      call invertks_2part(target_rho, nspin, hm%vhxc, sys%gr)
    else ! iterative case
      if (sys%ks%ks_inversion%method == XC_INVERSION_METHOD_VS_ITER) then ! iterative procedure for v_s 
        call invertks_iter(target_rho, np, nspin, hm, sys%gr, &
             sys%ks%ks_inversion%aux_st, sys%ks%ks_inversion%eigensolver)
      else
        call invertvxc_iter(target_rho, np, nspin, hm, sys%ks%xc, &
             sys%ks%hartree_solver, sys%ks%frozen_hxc, sys%gr, &
             sys%ks%ks_inversion%aux_st, sys%ks%ks_inversion%eigensolver)
      endif
    end if

    ! output quality of KS inversion
    call eigensolver_run(sys%ks%ks_inversion%eigensolver, sys%gr, &
         sys%ks%ks_inversion%aux_st, hm, 1, verbose = .false.)
    call states_calc_dens(sys%ks%ks_inversion%aux_st, sys%gr)

    diffdensity = M_ZERO
    do jj = 1, nspin
      do ii = 1, np
        if (abs(sys%ks%ks_inversion%aux_st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
          diffdensity = abs(sys%ks%ks_inversion%aux_st%rho(ii,jj)-target_rho(ii,jj))
          idiffmax=ii
        endif
      enddo
    enddo
    write (message(1),'(a,F16.6)') 'Achieved difference in densities wrt target:', &
        diffdensity
    call write_info(1)

    ! output for all cases    
    call h_sys_output_all(sys%outp, sys%gr, sys%geo, sys%ks%ks_inversion%aux_st, &
         sys%ks%ks_inversion%aux_hm, STATIC_DIR)
    
    if (sys%ks%ks_inversion%method == XC_INVERSION_METHOD_VXC_ITER) then
      call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vxc", sys%gr%mesh, hm%vxc(:,1), units_out%energy, ierr)
    else
      call doutput_function(io_function_fill_how("AxisX"), &
           ".", "vks", sys%gr%mesh, hm%vhxc(:,1), units_out%energy, ierr)
    endif
    
    call doutput_function(io_function_fill_how("AxisX"), &
           ".", "rho", sys%gr%mesh, sys%ks%ks_inversion%aux_st%rho(:,1), units_out%length**(-sys%gr%sb%dim), ierr)

    SAFE_DEALLOCATE_A(target_rho)

    call xc_ks_inversion_end(sys%ks%ks_inversion, sys%gr, sys%geo)

    call pop_sub('invert_ks.invert_ks_run')
    
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
      !%Default <tt>target_density.dat</tt>
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

      call pop_sub('invert_ks.invert_ks_run.read_target_rho')
    end subroutine read_target_rho

  end subroutine invert_ks_run

  
end module invert_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
