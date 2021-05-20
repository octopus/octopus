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

module invert_ks_oct_m
  use comm_oct_m
  use density_oct_m
  use eigensolver_oct_m
  use global_oct_m
  use hamiltonian_elec_oct_m
  use output_oct_m
  use io_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multisystem_basic_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pcm_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use xc_ks_inversion_oct_m
  
  implicit none

  private
  public :: invert_ks_run

contains

  ! ---------------------------------------------------------
  subroutine invert_ks_run(system)
    class(*), intent(inout) :: system

    PUSH_SUB(invert_ks_run)

    select type (system)
    class is (multisystem_basic_t)
      message(1) = "CalculationMode = invert_ks not implemented for multi-system calculations"
      call messages_fatal(1)
    type is (electrons_t)
      call invert_ks_run_legacy(system)
    end select

    POP_SUB(invert_ks_run)
  end subroutine invert_ks_run

  ! ---------------------------------------------------------
  subroutine invert_ks_run_legacy(sys)
    type(electrons_t), intent(inout) :: sys

    integer :: ii, jj, np, ndim, nspin
    integer :: err
    FLOAT   :: diffdensity
    FLOAT, allocatable :: target_rho(:,:), rho(:)
    type(restart_t) :: restart
      
    PUSH_SUB(invert_ks_run_legacy)

    if (sys%hm%pcm%run_pcm) then
      call messages_not_implemented("PCM for CalculationMode /= gs or td")
    end if

    if (sys%kpoints%use_symmetries) then
      call messages_experimental("KPoints symmetries with CalculationMode = invert_ks")
    end if

    ! initialize KS inversion module
    call xc_ks_inversion_write_info(sys%ks%ks_inversion, stdout)

    !abbreviations
    np      = sys%gr%mesh%np
    ndim    = sys%space%dim
    nspin   = sys%st%d%nspin

    ! read target density
    SAFE_ALLOCATE(target_rho(1:np, 1:nspin)) 
    SAFE_ALLOCATE(rho(1:np)) 
       
    call read_target_rho()

    sys%hm%energy%intnvxc = M_ZERO
    sys%hm%energy%hartree = M_ZERO
    sys%hm%energy%exchange = M_ZERO
    sys%hm%energy%correlation = M_ZERO
    sys%hm%vxc = M_ZERO

    ! calculate total density
    
    rho = M_ZERO
    do ii = 1, nspin
      rho(:) = rho(:) + target_rho(:,ii)
    end do
    
    ! calculate the Hartree potential
    call dpoisson_solve(sys%hm%psolver, sys%hm%vhartree, rho)

    do ii = 1, nspin
      sys%hm%vhxc(1:np, ii) = sys%hm%vhartree(1:np)
    end do

    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)
    call eigensolver_run(sys%ks%ks_inversion%eigensolver, sys%namespace, sys%gr, sys%ks%ks_inversion%aux_st, sys%hm, 1)
    call density_calc(sys%ks%ks_inversion%aux_st, sys%gr, sys%ks%ks_inversion%aux_st%rho)
    
    write(message(1),'(a)') "Calculating KS potential"
    call messages_info(1)
       
    if (sys%ks%ks_inversion%method == XC_INV_METHOD_TWO_PARTICLE) then ! 2-particle exact inversion
     
      call invertks_2part(target_rho, nspin, sys%hm, sys%gr, &
             sys%ks%ks_inversion%aux_st, sys%ks%ks_inversion%eigensolver, sys%namespace, sys%space, &
             sys%ks%ks_inversion%asymp)
     
    else ! iterative case
      if (sys%ks%ks_inversion%method >= XC_INV_METHOD_VS_ITER .and. &
          sys%ks%ks_inversion%method <= XC_INV_METHOD_ITER_GODBY) then ! iterative procedure for v_s 
        call invertks_iter(target_rho, sys%namespace, sys%space, nspin, sys%hm, sys%gr, &
             sys%ks%ks_inversion%aux_st, sys%ks%ks_inversion%eigensolver, sys%ks%ks_inversion%asymp,&
             sys%ks%ks_inversion%method)
      end if
    end if

    ! output quality of KS inversion
    
    call hamiltonian_elec_update(sys%hm, sys%gr%mesh, sys%namespace, sys%space)
    
    call eigensolver_run(sys%ks%ks_inversion%eigensolver, sys%namespace, sys%gr, sys%ks%ks_inversion%aux_st, sys%hm, 1)
    
    call density_calc(sys%ks%ks_inversion%aux_st, sys%gr, sys%ks%ks_inversion%aux_st%rho)

    diffdensity = M_ZERO
    do jj = 1, nspin
      do ii = 1, np
        if (abs(sys%ks%ks_inversion%aux_st%rho(ii,jj)-target_rho(ii,jj)) > diffdensity) then
          diffdensity = abs(sys%ks%ks_inversion%aux_st%rho(ii,jj)-target_rho(ii,jj))
        end if
      end do
    end do
    write (message(1),'(a,F16.6)') 'Achieved difference in densities wrt target:', &
        diffdensity
    call messages_info(1)

    ! output for all cases    
    call output_all(sys%outp, sys%namespace, sys%space, STATIC_DIR, sys%gr, sys%ions, -1, &
      sys%ks%ks_inversion%aux_st, sys%hm, sys%ks)

    sys%ks%ks_inversion%aux_st%dom_st_kpt_mpi_grp = sys%st%dom_st_kpt_mpi_grp
    ! save files in restart format
    call restart_init(restart, sys%namespace, RESTART_GS, RESTART_TYPE_DUMP, sys%mc, err, mesh = sys%gr%mesh)
    call states_elec_dump(restart, sys%space, sys%ks%ks_inversion%aux_st, sys%gr%mesh, sys%kpoints, err, 0)
    if (err /= 0) then
      message(1) = "Unable to write states wavefunctions."
      call messages_warning(1)
    end if
    call restart_end(restart)

    SAFE_DEALLOCATE_A(target_rho)
    SAFE_DEALLOCATE_A(rho)
    
    POP_SUB(invert_ks_run_legacy)
    
  contains

  ! ---------------------------------------------------------
    subroutine read_target_rho()
      character(len=256) :: filename
      integer :: pass, iunit, ierr, ii, npoints
      FLOAT   :: l_xx(MAX_DIM), l_ff(4), rr
      FLOAT, allocatable :: xx(:,:), ff(:,:)
      character(len=1)   :: char

      PUSH_SUB(invert_ks_run_legacy.read_target_rho)

      ! FIXME: just use restart/gs/density*.obf file rather than needing to set this and use a different format.
      
      !%Variable InvertKSTargetDensity
      !%Type string
      !%Default <tt>target_density.dat</tt>
      !%Section Calculation Modes::Invert KS
      !%Description
      !% Name of the file that contains the density used as the target in the 
      !% inversion of the KS equations.
      !%End
      call parse_variable(sys%namespace, 'InvertKSTargetDensity', "target_density.dat", filename)

      iunit = io_open(filename, sys%namespace, action='read', status='old')

      npoints = 0
      do pass = 1, 2
        ii = -1
        rewind(iunit)
        do
          if(ii== -1) then
            read(iunit, '(a)', advance='no') char
            ii = 0
            if (char  ==  '#') read(iunit, '(a)') char
          end if
          read(iunit, fmt=*, iostat=ierr) l_xx(1:ndim), l_ff(1:nspin)
          if(ierr /= 0) exit
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
        rr = rr + dmf_integrate(sys%gr%mesh, target_rho(:, ii), reduce = .false.)
      end do
      if(sys%gr%mesh%parallel_in_domains) then
        call sys%gr%mesh%allreduce(rr)
      end if
      rr = sys%st%qtot/rr
      target_rho(:,:) = rr*target_rho(:,:)

      SAFE_DEALLOCATE_A(xx)
      SAFE_DEALLOCATE_A(ff)

      POP_SUB(invert_ks_run_legacy.read_target_rho)
    end subroutine read_target_rho

  end subroutine invert_ks_run_legacy

  
end module invert_ks_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
