!!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: kdotp.F90 4145 2008-05-02 23:29:41Z xavier $

#include "global.h"
#define RESTART_DIR "kdotp/"

module kdotp_lr_m
  use datasets_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kdotp_calc_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use h_sys_output_m
  use pert_m
  use pol_lr_m
  use restart_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m
  use v_ks_m
  
  implicit none

  private

  public :: &
       kdotp_lr_run,       &
       read_wfs

  type kdotp_t
    type(pert_t) :: perturbation

    FLOAT, pointer :: eff_mass_inv(:, :, :, :)  ! inverse effective mass tensor
                                                ! (ik, ist, idir1, idir2)

    type(lr_t), pointer :: lr(:,:) ! linear response for (NDIM,1)
                                   ! second index is dummy; should only be 1
                                   ! for compatibility with em_resp routines

    logical :: ok
  end type kdotp_t

contains

  ! ---------------------------------------------------------
  subroutine kdotp_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(kdotp_t)           :: kdotp_vars
    type(sternheimer_t)     :: sh

    integer :: idir, ierr
!    integer :: sigma, ndim, i, dir, ierr, iomega, ifactor
    character(len=100) :: dirname, str_tmp
!    logical :: complex_response
    !, have_to_calculate

    call push_sub('kdotp.static_kdotp_lr_run')

    gr => sys%gr
!    ndim = sys%gr%sb%dim

    ALLOCATE(kdotp_vars%eff_mass_inv(sys%st%d%nik, sys%st%nst, NDIM, NDIM), sys%st%d%nik * sys%st%nst * NDIM * NDIM)
    kdotp_vars%eff_mass_inv(:,:,:,:)=0

    !    call parse_input()
    call pert_init(kdotp_vars%perturbation, PERTURBATION_KDOTP, sys%gr, sys%geo)

    ALLOCATE(kdotp_vars%lr(1:NDIM, 1), NDIM)

    call read_wfs(sys%st, sys%gr, sys%geo, .true.)
    ! even if wfs are real, the response must be allowed to be complex

    kdotp_vars%lr(1:NDIM,:)%nst = sys%st%nst

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    call system_h_setup(sys, h)
    
    call sternheimer_init(sh, sys, h, "KdotP", hermitian = wfs_are_real(sys%st), ham_var_set = 0)
    ! ham_var_set = 0 results in HamiltonianVariation = V_ext_only

    do idir = 1, NDIM
      call lr_init(kdotp_vars%lr(idir, 1))
      call lr_allocate(kdotp_vars%lr(idir, 1), sys%st, sys%gr%m)

      ! load wave-functions
      if(.not.fromScratch) then
         str_tmp =  kdotp_wfs_tag(idir)
         write(dirname,'(3a, i1)') RESTART_DIR, trim(str_tmp), '_1'
         ! 1 is the sigma index which is used in em_resp
         call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
               ierr, lr=kdotp_vars%lr(idir, 1))
          
          if(ierr.ne.0) then
             message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//dirname
             call write_warning(1)
          end if
          
       end if

    end do

    call io_mkdir(trim(tmpdir)//RESTART_DIR)

    call info()

    message(1) = "Info: Calculating effective masses."
    call write_info(1)

    call io_mkdir('kdotp/')

    do idir = 1, NDIM
      kdotp_vars%ok = .true.
      write(message(1), '(a,i3)') 'Info: Calculating response for direction ', idir
      call write_info(1)
      call pert_setup_dir(kdotp_vars%perturbation, idir)
!      write(*,*) 'done with pert_setup_dir'
!      if (wfs_are_complex(sys%st)) then
!        write(*,*) 'calling zsternheimer_solve'
      call zsternheimer_solve(sh, sys, h, kdotp_vars%lr(idir,:), 1, M_Z0, &
        kdotp_vars%perturbation, RESTART_DIR, &
        kdotp_rho_tag(idir), kdotp_wfs_tag(idir), have_restart_rho=(ierr==0))
!        write(*,*) 'called zsternheimer_solve'
!      else
!        write(*,*) 'calling dsternheimer_solve'
!        call dsternheimer_solve(sh, sys, h, kdotp_vars%lr(idir,:), 1, M_ZERO, &
!          kdotp_vars%perturbation, RESTART_DIR, &
!          kdotp_rho_tag(idir), kdotp_wfs_tag(idir), have_restart_rho=(ierr==0))
!        write(*,*) 'called dsternheimer_solve'
!      end if
        
      kdotp_vars%ok = kdotp_vars%ok .and. sternheimer_has_converged(sh)
          
    end do ! idir

!    if(wfs_are_complex(sys%st)) then 
!      call zlr_calc_eff_mass_inv(sys, h, kdotp_vars)
    call zlr_calc_eff_mass_inv(sys, h, kdotp_vars%lr, &
        kdotp_vars%perturbation, kdotp_vars%eff_mass_inv)
!    else
!      call dlr_calc_eff_mass_inv(sys, h, kdotp_vars)
!      call dlr_calc_eff_mass_inv(sys, h, kdotp_vars%lr, &
!        kdotp_vars%perturbation, kdotp_vars%eff_mass_inv)
!    end if

!    call effective_masses(kdotp_vars)
    call kdotp_output(sys%st, sys%gr, kdotp_vars)

!    do iomega = 1, em_vars%nomega
!
!      em_vars%ok(1:3) = .true.
!
!      do ifactor = 1, em_vars%nfactor
!        do dir = 1, sys%gr%sb%dim
!
!          ierr = 0
!
!          str_tmp = freq2str(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor)
!          write(message(1), '(a,i1,2a)') 'Info: Calculating response for direction ', dir, &
!            ' and frequency ' , trim(str_tmp)
!          call write_info(1)
!
!          have_to_calculate = .true.
!          
!          ! if this frequency is zero and this is not the first
!          ! iteration we do not have to do anything
!          if( iomega > 1 .and. em_vars%freq_factor(ifactor) == M_ZERO) have_to_calculate = .false. 
!            
!          if(ifactor > 1) then 
!
!            ! if this frequency is the same as the previous one, just copy it
!            if( have_to_calculate .and. &
!              em_vars%freq_factor(ifactor)*em_vars%omega(iomega) == &
!              em_vars%freq_factor(ifactor-1)*em_vars%omega(iomega) ) then 
!              
!              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 1, ifactor-1), em_vars%lr(dir, 1, ifactor))
!              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 2, ifactor-1), em_vars%lr(dir, 2, ifactor))
!              
!              have_to_calculate = .false.
!              
!            end if
!
!            ! if this frequency is minus the previous one, copy it inverted
!            if( have_to_calculate .and. & 
!                 em_vars%freq_factor(ifactor) == -em_vars%freq_factor(ifactor-1) ) then 
!              
!              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 1, ifactor-1), em_vars%lr(dir, 2, ifactor))
!              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 2, ifactor-1), em_vars%lr(dir, 1, ifactor))
!              
!              have_to_calculate = .false.
!              
!            end if
!
!          end if
!
!          if(have_to_calculate) then 
!
!            if(.not. fromscratch) then 
!
!              !try to load restart density
!              if (wfs_are_complex(sys%st)) then 
!                call zrestart_read_lr_rho(kdotp_vars%lr(idir, ik), sys%gr, sys%st%d%nspin, &
!                     RESTART_DIR, &
!                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir), ierr)
!              else 
!                call drestart_read_lr_rho(kdotp_vars%lr(idir, ik), sys%gr, sys%st%d%nspin, &
!                     RESTART_DIR, &
!                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir), ierr)
!              end if
!
!              !search for the density of the closest frequency
!              if(ierr /= 0) then 
!                
!                closest_omega = em_vars%freq_factor(ifactor)*em_vars%omega(iomega)
!                call oct_search_file_lr(closest_omega, dir, ierr, trim(tmpdir)//RESTART_DIR)
!                
!                !attempt to read 
!                if(ierr == 0 ) then 
!                  if (wfs_are_complex(sys%st)) then 
!                    call zrestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
!                         RESTART_DIR, em_rho_tag(closest_omega, dir), ierr)
!                  else 
!                    call drestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
!                         RESTART_DIR, em_rho_tag(closest_omega, dir), ierr)
!                  end if
!                end if
!                
!              end if
!
!              if(ierr == 0 .and. em_vars%nsigma == 2 ) then 
!                if (wfs_are_complex(sys%st)) then 
!                  em_vars%lr(dir, 2, ifactor)%zdl_rho = conjg(em_vars%lr(dir, 1, ifactor)%zdl_rho)
!                else 
!                  em_vars%lr(dir, 2, ifactor)%ddl_rho = em_vars%lr(dir, 1, ifactor)%ddl_rho
!                end if
!              end if
!
!            end if ! .not. fromscratch
!            
!            call pert_setup_dir(em_vars%perturbation, dir)
!            if (wfs_are_complex(sys%st)) then 
!              call zsternheimer_solve(sh, sys, h, em_vars%lr(dir, :, ifactor), em_vars%nsigma , &
!                 em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, &
!                 em_vars%perturbation, RESTART_DIR,&
!                 em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir),&
!                 em_wfs_tag(dir, ifactor), have_restart_rho=(ierr==0))
!            else
!              call dsternheimer_solve(sh, sys, h, em_vars%lr(dir, :, ifactor), em_vars%nsigma , &
!                 em_vars%freq_factor(ifactor)*em_vars%omega(iomega), &
!                 em_vars%perturbation, RESTART_DIR,&
!                 em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir),&
!                 em_wfs_tag(dir, ifactor), have_restart_rho=(ierr==0))
!            end if
!            
!            em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
!            
!          end if
!
!        end do ! dir
!      end do ! ifactor

    do idir = 1, NDIM
      call lr_dealloc(kdotp_vars%lr(idir, 1))
    end do

    call sternheimer_end(sh)
    call pert_end(kdotp_vars%perturbation)

    deallocate(kdotp_vars%lr)
    call states_deallocate_wfns(sys%st)
    deallocate(kdotp_vars%eff_mass_inv)

    call pop_sub()

  contains

!    ! ---------------------------------------------------------
!    subroutine parse_input()
!      type(block_t) :: blk
!      integer   :: nrow
!      integer   :: number, j, k
!      FLOAT     :: omega_ini, omega_fin, domega
!
!      call push_sub('kdotp.parse_input')
!
!      call obsolete_variable('PolFreqs               ', 'EMFreqs             ')
!      call obsolete_variable('PolHyper               ', 'EMHyperpol          ')
!      call obsolete_variable('PolEta                 ', 'EMEta               ')
!      call obsolete_variable('PolConvAbsDens         ', 'LRConvAbsDens       ')
!      call obsolete_variable('PolHamiltonianVariation', 'HamiltonianVariation')
!
!      !%Variable EMFreqs
!      !%Type block
!      !%Section Linear Response::Polarizabilities
!      !%Description
!      !% This block defines for which frequencies the polarizabilities
!      !% will be calculated. If is not present the static (omega = 0) response
!      !% is calculated.
!      !%
!      !% Each row of the block indicates a sequence of frequency values, the
!      !% first column is an integer that indicates the number of steps, the
!      !% second number is the initial frequency, and the third number the final
!      !% frequency. If the first number is one, then only the initial value is
!      !% considered. The block can have any number of rows. Consider the next example:
!      !%
!      !% <tt>%EMFreqs
!      !% <br>31 | 0.0 | 1.0
!      !% <br> 1 | 0.32
!      !% <br>%</tt>
!      !%
!      !%End
!
!      if (loct_parse_block(check_inp('EMFreqs'), blk) == 0) then 
!
!        nrow = loct_parse_block_n(blk)
!        em_vars%nomega = 0
!
!        !count the number of frequencies
!        do i = 0, nrow-1
!          call loct_parse_block_int(blk, i, 0, number)
!          em_vars%nomega = em_vars%nomega + number
!        end do
!
!/*        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega) */
!
!        !read frequencies
!        j = 1
!        do i = 0, nrow-1
!          call loct_parse_block_int(blk, i, 0, number)
!          call loct_parse_block_float(blk, i, 1, omega_ini)
!          if(number > 1) then 
!            call loct_parse_block_float(blk, i, 2, omega_fin)
!            domega = (omega_fin - omega_ini)/(number - M_ONE)
!            do k = 0, number-1
!              em_vars%omega(j + k) = (omega_ini + domega*k)*units_inp%energy%factor
!            end do
!            j = j + number
!          else
!            em_vars%omega(j) = omega_ini*units_inp%energy%factor
!            j = j + 1
!          end if
!        end do
!
!        call loct_parse_block_end(blk)
!
!        call sort(em_vars%omega)
!
!      else
!        !there is no frequency block, we calculate response for w = 0.0
!        em_vars%nomega = 1
!/*        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega) */
!        em_vars%omega(1) = M_ZERO
!      end if
!
!      !%Variable EMEta
!      !%Type float
!      !%Default 0.0
!      !%Section Linear Response::Polarizabilities
!      !%Description
!      !% Imaginary part of the frequency.
!      !%End
!
!      call loct_parse_float(check_inp('EMEta'), M_ZERO, em_vars%eta)
!      em_vars%eta = em_vars%eta*units_inp%energy%factor
!
!      ! reset the values of these variables
!      em_vars%calc_hyperpol = .false.
!      em_vars%freq_factor(1:MAX_DIM) = M_ONE
!
!      call pert_init(em_vars%perturbation, sys%gr, sys%geo)
!
!      if(pert_type(em_vars%perturbation) == PERTURBATION_ELECTRIC) then
!        !%Variable EMHyperpol
!        !%Type block
!        !%Section Linear Response::Polarizabilities
!        !%Description
!        !% This blocks describes the multiples of the frequency used for
!        !% the dynamic hyperpolarizability.
!        !%End
!
!        if (loct_parse_block(check_inp('EMHyperpol'), blk) == 0) then 
!          call loct_parse_block_float(blk, 0, 0, em_vars%freq_factor(1))
!          call loct_parse_block_float(blk, 0, 1, em_vars%freq_factor(2))
!          call loct_parse_block_float(blk, 0, 2, em_vars%freq_factor(3))
!          
!          call loct_parse_block_end(blk)
!          
!          em_vars%calc_hyperpol = .true.
!        end if
!      end if
!
!      call pop_sub()
!
!    end subroutine parse_input
!
!
    ! ---------------------------------------------------------
    subroutine info()

      call pert_info(kdotp_vars%perturbation, stdout)

        write(message(1),'(a)') 'Effective masses'
        call messages_print_stress(stdout, trim(message(1)))

!      if (wfs_are_real(sys%st)) then 
!        message(1) = 'Wavefunctions type: Real'
!      else
     message(1) = 'Wavefunctions type: Complex'
!      end if
      call write_info(1)

!      write(message(1),'(a,i3,a)') 'Calculating response for ', em_vars%nomega, ' frequencies.'
!      call write_info(1)
!
      call messages_print_stress(stdout)

    end subroutine info

  end subroutine kdotp_lr_run

  ! ---------------------------------------------------------
  subroutine kdotp_output(st, gr, kdotp_vars)
    type(states_t),       intent(inout) :: st
    type(grid_t),         intent(inout) :: gr
    type(kdotp_t),        intent(inout) :: kdotp_vars

    integer :: iunit, ik
    FLOAT :: determinant

!    write(*,'(a)') 'kdotp_output'

!    write(*,*) 'Number of states = ', st%nst
!    write(*,*) 'Number of k-points = ', st%d%nik

    do ik = 1, st%d%nik
      call out_eff_mass(ik, st%nst)
    enddo

    contains

      subroutine out_eff_mass(ik, nst)

        integer, intent(in) :: ik, nst

        character(len=80) :: filename
        integer :: ist

        write(filename, '(a, a)') 'kdotp/kpoint_', int2str(ik)
        iunit = io_open(trim(filename), action='write')
        write(iunit,'(a, a)') '# k-point index = ', int2str(ik)
        write(iunit,'(a, 3f12.8)') '# k-point coordinates = ', st%d%kpoints(1:MAX_DIM, ik)
        if (.not.kdotp_vars%ok) write(iunit, '(a)') "# WARNING: not converged"      

        write(iunit,'(a)')
        write(iunit,'(a)') '# Inverse effective mass tensors'
        do ist = 1, nst
          write(iunit,'(a)')
          write(iunit,'(a, a, a, f12.8, a, a)') 'State #', int2str(ist), ', Energy = ', &
            st%eigenval(ist, ik)*units_out%energy%factor, ' ', units_out%energy%abbrev
          call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), NDIM, M_ONE)
        enddo

        write(iunit,'(a)')
        write(iunit,'(a)') '# Effective mass tensors'
        do ist = 1, nst
          write(iunit,'(a)')
          write(iunit,'(a, a, a, f12.8, a, a)') 'State #', int2str(ist), ', Energy = ', &
            st%eigenval(ist, ik)*units_out%energy%factor, ' ', units_out%energy%abbrev
          determinant = lalg_inverter(gr%sb%dim, kdotp_vars%eff_mass_inv(ik, ist, :, :), .true.)
          call io_output_tensor(iunit, kdotp_vars%eff_mass_inv(ik, ist, :, :), NDIM, M_ONE)
        enddo

      end subroutine out_eff_mass

      character(len=12) function int2str(i) result(str)
        integer, intent(in) :: i
      
        write(str, '(i11)') i
        str = trim(adjustl(str))

      end function int2str
            
!    ! ---------------------------------------------------------
!    subroutine out_polarizability()
!      FLOAT :: cross(MAX_DIM, MAX_DIM), crossp(MAX_DIM, MAX_DIM)
!      FLOAT :: average, anisotropy
!      
!
!      iunit = io_open(trim(dirname)//'/alpha', action='write')
!
!      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
!      write(iunit, '(2a)', advance='no') '# Polarizability tensor [', &
!        trim(units_out%length%abbrev)
!      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM
!      write(iunit, '(a)') ']'
!
!      call io_output_tensor(iunit, em_vars%alpha(:,:, ifactor), units_out%length%factor**NDIM)
!
!      call io_close(iunit)
!
!      ! CROSS SECTION (THE COMPLEX PART OF POLARIZABILITY)
!      if( wfs_are_complex(st)) then 
!        cross(1:MAX_DIM, 1:MAX_DIM) = aimag(em_vars%alpha(1:MAX_DIM, 1:MAX_DIM, ifactor)) * &
!          em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor * M_FOUR * M_PI / P_c 
!        
!        iunit = io_open(trim(dirname)//'/cross_section', action='write')
!        if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
!
!        average = M_THIRD*(cross(1, 1) + cross(2, 2) + cross(3, 3))
!        crossp(:, :) = matmul(cross(:, :), cross(:, :))
!        anisotropy = &
!             M_THIRD*(M_THREE*(crossp(1, 1) + crossp(2, 2) + crossp(3, 3)) - (cross(1, 1) + cross(2, 2) + cross(3, 3))**2)
!          
!        call cross_section_header(iunit)
!        write(iunit,'(3e20.8)', advance = 'no') &
!          em_vars%freq_factor(ifactor)*em_vars%omega(iomega) / units_out%energy%factor, average , sqrt(max(anisotropy, M_ZERO))
!        write(iunit,'(9e20.8)', advance = 'no') cross(1:3, 1:3)
!        write(iunit,'(a)', advance = 'yes')
!
!        call io_close(iunit)
!      end if
!
!    end subroutine out_polarizability

  end subroutine kdotp_output

end module kdotp_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
