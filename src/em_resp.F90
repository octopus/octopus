!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
#define RESTART_DIR "restart_pol_lr/"

module pol_lr_m
  use datasets_m
  use em_resp_calc_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use output_m
  use resp_pert_m
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
       pol_lr_run       

  public :: &
    read_wfs

  integer, parameter ::         &
     PERTURBATION_ELECTRIC = 1, &
     PERTURBATION_MAGNETIC = 2

  type em_resp_t
    type(resp_pert_t) :: perturbation

    integer :: nsigma ! 1: consider only positive values of the frequency
                      ! 2: consider both positive and negative
    integer :: nfactor! 1: only one frequency needed
                      ! 3: three frequencies (for the hyperpolarizabilities)
    integer :: nomega ! number of frequencies to consider

    FLOAT :: eta                     ! small imaginary part to add to the frequency
    FLOAT :: freq_factor(MAX_DIM)    !
    FLOAT,      pointer :: omega(:)  ! the frequencies to consider
    type(lr_t), pointer :: lr(:,:,:) ! linear response for (NDIM, nsigma, nomega)

    logical :: calc_hyperpol
    CMPLX   :: alpha(MAX_DIM, MAX_DIM, 3)        ! the linear polarizability
    CMPLX   :: beta (MAX_DIM, MAX_DIM, MAX_DIM)  ! first hyperpolarizability

    CMPLX   :: chi_para(MAX_DIM, MAX_DIM, 3)     ! The paramagnetic part of the susceptibility
    CMPLX   :: chi_dia (MAX_DIM, MAX_DIM, 3)     ! The diamagnetic  part of the susceptibility

    logical :: ok(1:3)
  end type em_resp_t

contains

  ! ---------------------------------------------------------
  subroutine pol_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(grid_t),   pointer :: gr
    type(em_resp_t)         :: em_vars
    type(sternheimer_t)     :: sh

    integer :: sigma, ndim, i, dir, ierr, iomega, ifactor
    character(len=80) :: dirname
    logical :: complex_response, have_to_calculate

    FLOAT :: closest_omega


    call push_sub('em_resp.static_pol_lr_run')

    gr => sys%gr
    ndim = sys%gr%sb%dim

    call parse_input()

    em_vars%nfactor = 1
    if(em_vars%calc_hyperpol) em_vars%nfactor=3

    em_vars%nsigma = 1  ! positive and negative values of the frequency must be considered
    if(em_vars%calc_hyperpol.or.em_vars%nomega > 1.or.em_vars%omega(1).ne.M_ZERO) em_vars%nsigma = 2

    complex_response = (em_vars%eta /= M_ZERO ) .or. wfs_are_complex(sys%st)

    ALLOCATE(em_vars%lr(1:NDIM, 1:em_vars%nsigma, 1:em_vars%nfactor), NDIM*em_vars%nsigma*em_vars%nfactor)

    call read_wfs(sys%st, sys%gr, sys%geo, complex_response)
    em_vars%lr(1:NDIM, 1:em_vars%nsigma, 1:em_vars%nfactor)%nst = sys%st%nst

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    call system_h_setup(sys, h)
    
    ! The magnetic perturbation does not change the density
    ! WARNING: This will probably not work if spin-orbit coupling is on
    if(em_vars%perturbation%resp_type == RESP_PERTURBATION_MAGNETIC) then
      h%ip_app = .true.
    end if

    call sternheimer_init(sh, sys, h, "Pol", hermitian = wfs_are_real(sys%st))

    do dir = 1, ndim
      do sigma = 1, em_vars%nsigma
        do ifactor = 1, em_vars%nfactor 

          call lr_init(em_vars%lr(dir, sigma, ifactor))
          call lr_allocate(em_vars%lr(dir, sigma, ifactor), sys%st, sys%gr%m)

          ! load wave-functions
          if(.not.fromScratch) then
            write(dirname,'(a, i1)') RESTART_DIR//trim(em_wfs_tag(dir,ifactor))//'_', sigma
            call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, sys%geo, &
              ierr, lr=em_vars%lr(dir, sigma, ifactor))

            if(ierr.ne.0) then
              message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//dirname
              call write_warning(1)
            end if

          end if

        end do
      end do
    end do

    call io_mkdir(trim(tmpdir)//RESTART_DIR)
    call info()

    message(1) = "Info: Calculating polarizabilities."
    call write_info(1)

    call io_mkdir('linear/')

    do iomega = 1, em_vars%nomega

      em_vars%ok(1:3) = .true.

      do ifactor = 1, em_vars%nfactor
        do dir = 1, sys%gr%sb%dim

          ierr = 0

          write(message(1), '(a,i1,a)') 'Info: Calculating response for direction ', dir, &
               ' and frequency '//&
               trim(freq2str(em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor))
          call write_info(1)

          have_to_calculate = .true.
          
          ! if this frequency is zero and this is not the first
          ! iteration we do not have to do anything
          if( iomega > 1 .and. em_vars%freq_factor(ifactor) == M_ZERO) then 
            have_to_calculate = .false. 
          end if
            
          if(ifactor > 1) then 

            ! if this frequency is the same as the previous one, just copy it
            if( have_to_calculate .and. &
              em_vars%freq_factor(ifactor)*em_vars%omega(iomega)==&
              em_vars%freq_factor(ifactor-1)*em_vars%omega(iomega) ) then 
              
              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 1, ifactor-1), em_vars%lr(dir, 1, ifactor))
              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 2, ifactor-1), em_vars%lr(dir, 2, ifactor))
              
              have_to_calculate = .false.
              
            end if

            ! if this frequency is minus the previous one, copy it inverted
            if( have_to_calculate .and. & 
                 em_vars%freq_factor(ifactor) == -em_vars%freq_factor(ifactor-1) ) then 
              
              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 1, ifactor-1), em_vars%lr(dir, 2, ifactor))
              call lr_copy(sys%st, sys%gr%m, em_vars%lr(dir, 2, ifactor-1), em_vars%lr(dir, 1, ifactor))
              
              have_to_calculate = .false.
              
            end if

          end if

          if(have_to_calculate) then 

            if( .not. fromscratch) then 

              !try to load restart density
              if (wfs_are_complex(sys%st)) then 
                call zrestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                     RESTART_DIR, &
                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir), ierr)
              else 
                call drestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                     RESTART_DIR, &
                     em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir), ierr)
              end if

              !search for the density of the closest frequency
              if( ierr /= 0) then 
                
                closest_omega = em_vars%freq_factor(ifactor)*em_vars%omega(iomega)
                call oct_search_file_lr(closest_omega, dir, ierr, trim(tmpdir)//RESTART_DIR)
                
                !atempt to read 
                if(ierr == 0 ) then 
                  if (wfs_are_complex(sys%st)) then 
                    call zrestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                         RESTART_DIR, em_rho_tag(closest_omega, dir), ierr)
                  else 
                    call drestart_read_lr_rho(em_vars%lr(dir, 1, ifactor), sys%gr, sys%st%d%nspin, &
                         RESTART_DIR, em_rho_tag(closest_omega, dir), ierr)
                  end if
                end if
                
              end if

              if(ierr == 0) then 
                if (wfs_are_complex(sys%st)) then 
                  em_vars%lr(dir, 2, ifactor)%zdl_rho = conjg(em_vars%lr(dir, 1, ifactor)%zdl_rho)
                else 
                  em_vars%lr(dir, 2, ifactor)%zdl_rho = em_vars%lr(dir, 1, ifactor)%zdl_rho
                end if
              end if

            end if ! .not. fromscratch
            

            call resp_pert_setup_dir(em_vars%perturbation, dir)
            if (wfs_are_complex(sys%st)) then 
              call zsternheimer_solve(sh, sys, h, em_vars%lr(dir, :, ifactor), em_vars%nsigma , &
                 em_vars%freq_factor(ifactor)*em_vars%omega(iomega) + M_zI * em_vars%eta, &
                 em_vars%perturbation, RESTART_DIR,&
                 em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir),&
                 em_wfs_tag(dir, ifactor), have_restart_rho=(ierr==0))
            else
              call dsternheimer_solve(sh, sys, h, em_vars%lr(dir, :, ifactor), em_vars%nsigma , &
                 em_vars%freq_factor(ifactor)*em_vars%omega(iomega), &
                 em_vars%perturbation, RESTART_DIR,&
                 em_rho_tag(em_vars%freq_factor(ifactor)*em_vars%omega(iomega), dir),&
                 em_wfs_tag(dir, ifactor), have_restart_rho=(ierr==0))
            end if
            
            em_vars%ok(ifactor) = em_vars%ok(ifactor) .and. sternheimer_has_converged(sh)
            
          end if

        end do ! dir
      end do ! ifactor
      
      if(em_vars%perturbation%resp_type == RESP_PERTURBATION_ELECTRIC) then
        ! calculate polarizability
        do ifactor = 1, em_vars%nfactor
          if(wfs_are_complex(sys%st)) then 
            call zlr_calc_polarizability(sys, em_vars%lr(:,:, ifactor), em_vars%alpha(:,:, ifactor))
          else
            call dlr_calc_polarizability(sys, em_vars%lr(:,:, ifactor), em_vars%alpha(:,:, ifactor))
          end if
        end do
        
        ! calculate hyperpolarizability
        if(em_vars%calc_hyperpol) then
          if(wfs_are_complex(sys%st)) then
            call zlr_calc_beta(sh, sys, h, em_vars%lr(:,:,:), em_vars%perturbation, em_vars%beta)
          else
            call dlr_calc_beta(sh, sys, h, em_vars%lr(:,:,:), em_vars%perturbation, em_vars%beta)
          end if
        end if
      
      else if(em_vars%perturbation%resp_type == RESP_PERTURBATION_MAGNETIC) then
        do ifactor = 1, em_vars%nfactor
          if(wfs_are_complex(sys%st)) then 
            call zlr_calc_susceptibility(sys, h, em_vars%lr(:,:, ifactor), em_vars%perturbation, &
               em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          else
            call dlr_calc_susceptibility(sys, h, em_vars%lr(:,:, ifactor), em_vars%perturbation, &
               em_vars%chi_para(:,:, ifactor), em_vars%chi_dia(:,:, ifactor))
          end if
        end do
      end if

      call em_resp_output(sys%st, sys%gr, sys%outp, em_vars, iomega)

    end do

    do dir = 1, ndim
      do sigma = 1, em_vars%nsigma
        do ifactor = 1, em_vars%nfactor
          call lr_dealloc(em_vars%lr(dir, sigma, ifactor))
        end do
      end do
    end do

    call sternheimer_end(sh)
    call resp_pert_end(em_vars%perturbation)

    deallocate(em_vars%omega, em_vars%lr)
    call states_deallocate_wfns(sys%st)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine parse_input()
      C_POINTER :: blk
      integer   :: nrow
      integer   :: number, j, k
      FLOAT     :: omega_ini, omega_fin, domega

      call push_sub('em_resp.parse_input')

      !%Variable PolFreqs
      !%Type block
      !%Section Linear Response::Polarizabilities
      !%Description
      !% This block defines for which frequencies the polarizabilities
      !% will be calculated. If is not present the static (omega = 0) response
      !% is calculated.
      !%
      !% Each row of the block indicates a sequence of frequency values, the
      !% first column is an integer that indicates the number of steps, the
      !% second number is the initial frequency, and the third number the final
      !% frequency. If the first number is one, then only the initial value is
      !% considered. The block can have any number of rows. Consider the next example:
      !%
      !% <tt>%PolFreqs
      !% <br>31 | 0.0 | 1.0
      !% <br> 1 | 0.32
      !% <br>%</tt>
      !%
      !%End

      if (loct_parse_block(check_inp('PolFreqs'), blk) == 0) then 

        nrow = loct_parse_block_n(blk)
        em_vars%nomega = 0

        !count the number of frequencies
        do i = 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          em_vars%nomega = em_vars%nomega + number
        end do

        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega)

        !read frequencies
        j = 1
        do i = 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          call loct_parse_block_float(blk, i, 1, omega_ini)
          if(number > 1) then 
            call loct_parse_block_float(blk, i, 2, omega_fin)
            domega = (omega_fin - omega_ini)/(number - M_ONE)
            do k = 0, number-1
              em_vars%omega(j + k) = (omega_ini + domega*k) * units_inp%energy%factor
            end do
            j = j + number
          else
            em_vars%omega(j) = omega_ini * units_inp%energy%factor
            j = j + 1
          end if
        end do

        call loct_parse_block_end(blk)

        call sort(em_vars%omega)

      else
        !there is no frequency block, we calculate response for w = 0.0
        em_vars%nomega = 1
        ALLOCATE(em_vars%omega(1:em_vars%nomega), em_vars%nomega)
        em_vars%omega(1) = M_ZERO
      end if

      !%Variable PolEta
      !%Type float
      !%Default 0.0
      !%Section Linear Response::Polarizabilities
      !%Description
      !% Imaginary part of the frequency.
      !%End

      call loct_parse_float(check_inp('PolEta'), M_ZERO, em_vars%eta)
      em_vars%eta = em_vars%eta * units_inp%energy%factor

      ! reset the values of these variables
      em_vars%calc_hyperpol = .false.
      em_vars%freq_factor(1:MAX_DIM) = M_ONE

      call resp_pert_init(em_vars%perturbation, sys%gr, sys%geo)

      if(em_vars%perturbation%resp_type == PERTURBATION_ELECTRIC) then
        !%Variable PolHyper
        !%Type block
        !%Section Linear Response::Polarizabilities
        !%Description
        !% This blocks describes the multiples of the frequency used for
        !% the dynamic hyperpolarizability.
        !%End

        if (loct_parse_block(check_inp('PolHyper'), blk) == 0) then 
          call loct_parse_block_float(blk, 0, 0, em_vars%freq_factor(1))
          call loct_parse_block_float(blk, 0, 1, em_vars%freq_factor(2))
          call loct_parse_block_float(blk, 0, 2, em_vars%freq_factor(3))
          
          call loct_parse_block_end(blk)

          em_vars%calc_hyperpol = .true.
        end if
      end if

      call pop_sub()

    end subroutine parse_input


    ! ---------------------------------------------------------
    subroutine info()

      call resp_pert_info(em_vars%perturbation, stdout)
      if(em_vars%perturbation%resp_type == RESP_PERTURBATION_ELECTRIC) then
        if(em_vars%calc_hyperpol) then 
          write(message(1),'(a)') 'Linear Reponse First Order Hyperpolarizabilities'
          call messages_print_stress(stdout, trim(message(1)))
        else 
          write(message(1),'(a)') 'Linear Reponse Polarizabilities'
          call messages_print_stress(stdout, trim(message(1)))
        end if
      else
        write(message(1),'(a)') 'Magnetic susceptibilities'
        call messages_print_stress(stdout, trim(message(1)))
      end if

      if (wfs_are_real(sys%st)) then 
        message(1) = 'Wavefunctions type: Real'
      else
        message(1) = 'Wavefunctions type: Complex'
      end if
      call write_info(1)

      write(message(1),'(a,i3,a)') 'Calculating response for ', em_vars%nomega, ' frequencies.'
      call write_info(1)

      call messages_print_stress(stdout)

    end subroutine info

  end subroutine pol_lr_run


  ! ---------------------------------------------------------
  subroutine read_wfs(st, gr, geo, complex_wfs)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    logical,          intent(in)    :: complex_wfs

    integer :: kpoints, nst, ierr, dim
      
    !check how many wfs we have

    call push_sub('em_resp.read_wfs')

    call states_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)

    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    st%nst    = nst
    st%st_end = nst
    deallocate(st%eigenval, st%occ)

    if ( complex_wfs ) then 
      call states_allocate_wfns(st, gr%m, M_CMPLX)
    else 
      call states_allocate_wfns(st, gr%m, M_REAL)
    end if

    ALLOCATE(st%eigenval(st%nst, st%d%nik), st%nst*st%d%nik)
    ALLOCATE(st%occ(st%nst, st%d%nik), st%nst*st%d%nik)

    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%spin(3, st%nst, st%d%nik), st%nst*st%d%nik)
      st%spin = M_ZERO
    end if
    st%eigenval = huge(REAL_PRECISION)
    st%occ      = M_ZERO

    ! load wave-functions
    call restart_read(trim(tmpdir)//'restart_gs', st, gr, geo, ierr)  
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    call pop_sub()

  end subroutine read_wfs

  ! ---------------------------------------------------------
  subroutine em_resp_output(st, gr, outp, em_vars, iomega)
    type(states_t),  intent(inout) :: st
    type(grid_t),    intent(inout) :: gr
    type(output_t),  intent(in)    :: outp
    type(em_resp_t), intent(inout) :: em_vars
    integer,         intent(in)    :: iomega
    
    integer :: iunit, ifactor
    character(len=80) :: dirname

    do ifactor = 1, em_vars%nfactor
      write(dirname, '(a, a)') 'linear/freq_', trim(freq2str(&
           em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor))
      call io_mkdir(trim(dirname))

      if(em_vars%perturbation%resp_type == RESP_PERTURBATION_ELECTRIC) then
        call out_polarizability()
        if(em_vars%calc_hyperpol) call out_hyperpolarizability()
      else if(em_vars%perturbation%resp_type == RESP_PERTURBATION_MAGNETIC) then
        call out_susceptibility()
      end if

      call out_projections()

      write(dirname, '(a, a)') 'linear/freq_',trim(freq2str(&
           em_vars%omega(iomega)/units_out%energy%factor))
      call io_mkdir(dirname)
      call out_wavefunctions()
    end do

  contains

    ! ---------------------------------------------------------
    ! Note: this should be in spectrum.F90
    subroutine cross_section_header(out_file)
      integer, intent(in) :: out_file

      character(len=80) :: header_string
      integer :: i, k

      !this header is the same from sprectrum.F90
      write(out_file, '(a1, a20)', advance = 'no') '#', str_center("Energy", 20)
      write(out_file, '(a20)', advance = 'no') str_center("(1/3)*Tr[sigma]", 20)
      write(out_file, '(a20)', advance = 'no') str_center("Anisotropy[sigma]", 20)

      do i = 1, 3
        do k = 1, 3
          write(header_string,'(a6,i1,a1,i1,a1)') 'sigma(',i,',',k,')'
          write(out_file, '(a20)', advance = 'no') str_center(trim(header_string), 20)
        end do
      end do

      write(out_file, *)
      write(out_file, '(a1,a20)', advance = 'no') '#', str_center('['//trim(units_out%energy%abbrev) // ']', 20)
      do i = 1, 11
        write(out_file, '(a20)', advance = 'no')  str_center('['//trim(units_out%length%abbrev) //'^2]', 20)
      end do
      write(out_file,*)
    end subroutine cross_section_header


    ! ---------------------------------------------------------
    subroutine out_polarizability()
      FLOAT :: cross(MAX_DIM, MAX_DIM), crossp(MAX_DIM, MAX_DIM)
      FLOAT :: average, anisotropy
      

      iunit = io_open(trim(dirname)//'/alpha', action='write')

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"
      write(iunit, '(2a)', advance='no') '# Polarizability tensor [', &
        trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM
      write(iunit, '(a)') ']'

      call out_tensor(iunit, em_vars%alpha(:,:, ifactor), units_out%length%factor**NDIM)

      call io_close(iunit)

      ! CROSS SECTION (THE COMPLEX PART OF POLARIZABILITY)
      if( wfs_are_complex(st)) then 
        cross(1:MAX_DIM, 1:MAX_DIM) = aimag(em_vars%alpha(1:MAX_DIM, 1:MAX_DIM, ifactor)) * &
          em_vars%freq_factor(ifactor)*em_vars%omega(iomega)/units_out%energy%factor * M_FOUR * M_PI / P_c 
        
        iunit = io_open(trim(dirname)//'/cross_section', action='write')
        if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

        average = M_THIRD* ( cross(1, 1) + cross(2, 2) + cross(3, 3) )
        crossp(:, :) = matmul(cross(:, :),cross(:, :))
        anisotropy =  M_THIRD * ( M_THREE * (crossp(1, 1) + crossp(2, 2) + crossp(3, 3)) - &
          (cross(1, 1) + cross(2, 2) + cross(3, 3))**2 )
          
        call cross_section_header(iunit)
        write(iunit,'(3e20.8)', advance = 'no') &
          em_vars%freq_factor(ifactor)*em_vars%omega(iomega) / units_out%energy%factor, &
          average , sqrt(max(anisotropy, M_ZERO)) 
        write(iunit,'(9e20.8)', advance = 'no') cross(1:3, 1:3)
        write(iunit,'(a)', advance = 'yes')

        call io_close(iunit)
      end if

    end subroutine out_polarizability


    ! ---------------------------------------------------------
    subroutine out_susceptibility()
      iunit = io_open(trim(dirname)//'/susceptibility', action='write')

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call out_tensor(iunit, em_vars%chi_para(:,:, ifactor), CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm a.u.]'
      call out_tensor(iunit, em_vars%chi_dia(:,:, ifactor), CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm a.u.]'
      call out_tensor(iunit, em_vars%chi_para(:,:, ifactor) + em_vars%chi_dia(:,:, ifactor), CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(a)') hyphens

      write(iunit, '(2a)') '# Paramagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call out_tensor(iunit, em_vars%chi_para(:,:, ifactor), M_ONE/CNST(8.9238878e-2)*CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Diamagnetic contribution to the susceptibility tensor [ppm cgs / mol]'
      call out_tensor(iunit, em_vars%chi_dia(:,:, ifactor), M_ONE/CNST(8.9238878e-2)*CNST(1e-6))
      write(iunit, '(1x)')

      write(iunit, '(2a)') '# Total susceptibility tensor [ppm cgs / mol]'
      call out_tensor(iunit, em_vars%chi_para(:,:, ifactor) + em_vars%chi_dia(:,:, ifactor), M_ONE/CNST(8.9238878e-2)*CNST(1e-6))
      write(iunit, '(1x)')

      call io_close(iunit)      
    end subroutine out_susceptibility


    ! ---------------------------------------------------------
    subroutine out_tensor(iunit, tensor, factor)
      integer, intent(in) :: iunit
      CMPLX,   intent(in) :: tensor(:,:)
      FLOAT,   intent(in) :: factor

      FLOAT :: trace
      integer :: j

      trace = M_z0
      do j = 1, NDIM
        write(iunit, '(3f12.6)') real(tensor(j, 1:NDIM)) / factor
        trace = trace + real(tensor(j, j))
      end do
      trace = trace / M_THREE

      write(iunit, '(a, f12.6)')  'Isotropic average', trace / factor
      
    end subroutine out_tensor


    ! ---------------------------------------------------------
    subroutine out_hyperpolarizability()
      character, parameter :: axis(1:3) = (/ 'x', 'y', 'z' /)

      CMPLX :: bpar(1:MAX_DIM), bper(1:MAX_DIM), bk(1:MAX_DIM)
      integer :: i, j, k

      ! Output first hyperpolarizabilty (beta)
      iunit = io_open(trim(dirname)//'/beta', action='write')

      if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

      write(iunit, '(2a)', advance='no') 'First hyperpolarizability tensor: beta [', &
        trim(units_out%length%abbrev)
      write(iunit, '(a,i1)', advance='no') '^', 5
      write(iunit, '(a)') ']'

      if (st%d%nspin /= UNPOLARIZED ) then 
        write(iunit, '(a)') 'WARNING: Hyperpolarizability has not been tested for spin polarized systems'
      end if

      write(iunit, '()')

      do i = 1, NDIM
        do j = 1, NDIM
          do k = 1, NDIM
            write(iunit,'(a,e20.8,e20.8)') 'beta '//axis(i)//axis(j)//axis(k)//' ', &
              real( em_vars%beta(i, j, k))/units_out%length%factor**(5), &
              aimag(em_vars%beta(i, j, k))/units_out%length%factor**(5)
          end do
        end do
      end do

      if (NDIM == 3) then 
        bpar = M_ZERO
        bper = M_ZERO

        do i = 1, NDIM
          do j = 1, NDIM
            bpar(i) = bpar(i) + em_vars%beta(i, j, j) + &
              em_vars%beta(j, i, j) + em_vars%beta(j, j, i)

            bper(i) = bper(i) + M_TWO*em_vars%beta(i, j, j) - &
              M_THREE*em_vars%beta(j, i, j) + M_TWO*em_vars%beta(j, j, i)
          end do
        end do

        write(iunit, '()')

        bpar = bpar / (M_FIVE * units_out%length%factor**(5))
        bper = bper / (M_FIVE * units_out%length%factor**(5))
        bk(1:NDIM) = M_THREE*M_HALF*(bpar(1:NDIM) - bper(1:NDIM))
            
        do i = 1, NDIM
          write(iunit, '(a, 2e20.8)') 'beta // '//axis(i), real(bpar(i)), aimag(bpar(i))
        end do
        
        write(iunit, '()')
            
        do i = 1, NDIM
          write(iunit, '(a, 2e20.8)') 'beta _L '//axis(i), real(bper(i)), aimag(bper(i))
        end do

        write(iunit, '()')

        do i = 1, NDIM
          write(iunit, '(a, 2e20.8)') 'beta  k '//axis(i), real(bk(i)), aimag(bk(i))
        end do
      endif

      call io_close(iunit)

    end subroutine out_hyperpolarizability


    ! ---------------------------------------------------------
    subroutine out_projections()
      CMPLX   :: proj
      integer :: ist, ivar, ik, dir, sigma
      character(len=80) :: fname

      do ik = 1, st%d%nspin
        do dir = 1, NDIM

          write(fname, '(2a,i1,a,i1)') trim(dirname), '/projection-', ik, '-', dir
          iunit = io_open(trim(fname), action='write')

          if (.not.em_vars%ok(ifactor)) write(iunit, '(a)') "# WARNING: not converged"

          write(iunit, '(a)', advance='no') '# state '
          do ivar = 1, em_vars%lr(dir, 1, 1)%nst
            do sigma = 1, em_vars%nsigma

              if( sigma == em_vars%nsigma .and. ivar == em_vars%lr(dir, 1, 1)%nst) then 
                write(iunit, '(i3)', advance='yes') (3 - 2*sigma)*ivar
              else 
                write(iunit, '(i3)', advance='no') (3 - 2*sigma)*ivar
              end if

            end do
          end do

          do ist = 1, st%nst
            write(iunit, '(i3)', advance='no') ist

            do ivar = 1, em_vars%lr(dir, 1, 1)%nst
              do sigma = 1, em_vars%nsigma

                if(wfs_are_complex(st)) then
                  proj = zstates_dotp(gr%m, st%d%dim, &
                    st%zpsi(:,:, ist, ik),                &
                    em_vars%lr(dir, sigma, ifactor)%zdl_psi(:,:, ivar, ik))
                else
                  proj = dstates_dotp(gr%m, st%d%dim, &
                    st%dpsi(:,:, ist, ik),                &
                    em_vars%lr(dir, sigma, ifactor)%ddl_psi(:,:, ivar, ik))
                end if
                  
                if( sigma == em_vars%nsigma .and. ivar == em_vars%lr(dir, 1, 1)%nst) then 
                  write(iunit, '(f12.6)', advance='yes') abs(proj)
                else 
                  write(iunit, '(f12.6,a)', advance='no') abs(proj), ' '
                end if

              end do
            end do

          end do
          call io_close(iunit)

        end do ! dir
      end do !ik

    end subroutine out_projections


    ! ---------------------------------------------------------
    subroutine out_wavefunctions()
      integer :: dir, isigma

      do dir = 1, NDIM
        if( wfs_are_complex(st) ) then 

          if(NDIM==3) then
            if(iand(outp%what, output_elf).ne.0) &
              call zlr_calc_elf(st, gr, em_vars%lr(dir, 1, ifactor), em_vars%lr(dir, 2, ifactor))
          end if
          do isigma = 1, em_vars%nsigma
            call zlr_output(st, gr, em_vars%lr(dir, isigma, ifactor), dirname, dir, isigma, outp)
          end do
        else

          if(NDIM==3) then
            if(iand(outp%what, output_elf).ne.0) &
              call dlr_calc_elf(st, gr, em_vars%lr(dir, 1, ifactor), em_vars%lr(dir, 2, ifactor))
          end if
          do isigma = 1, em_vars%nsigma
            call dlr_output(st, gr, em_vars%lr(dir, isigma, ifactor), dirname, dir, isigma, outp)
          end do

        end if
      end do

    end subroutine out_wavefunctions

  end subroutine em_resp_output

  
#include "undef.F90"
#include "complex.F90"
#include "em_resp_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "em_resp_inc.F90"

end module pol_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
