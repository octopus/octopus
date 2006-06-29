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

module static_pol_lr_m
  use global_m
  use messages_m
  use units_m
  use mesh_m
  use mesh_function_m
  use system_m
  use restart_m
  use hamiltonian_m
  use mix_m
  use poisson_m
  use linear_response_m
  use io_m
  use lib_oct_parser_m
  use math_m

  implicit none

  private
  public :: &
    static_pol_lr_run
  
  type pol_props_t
    logical :: complex_response
    logical :: add_fxc
    logical :: add_hartree
    logical :: dynamic
  end type pol_props_t

  type status_t
    logical :: ok
  end type status_t
  
contains

  ! ---------------------------------------------------------
  subroutine static_pol_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(lr_t), allocatable :: lr(:,:,:) ! lr(NDIM,NS,NFREQ)
    type(grid_t),   pointer :: gr
    type(status_t)          :: status

    FLOAT ::  pol(1:MAX_DIM, 1:MAX_DIM)
    CMPLX :: zpol(1:MAX_DIM, 1:MAX_DIM)
    FLOAT :: hpol(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM)

    FLOAT, allocatable :: hpol_density(:,:,:,:)

    character(len=30) :: dirname

    FLOAT :: eta
    integer :: nfreq, nsigma, ndim, i, j, sigma, ierr, kpoints, dim, nst, k
    character(len=80) :: fname

    type(pol_props_t) :: props
    
    integer :: nomega
    FLOAT, allocatable :: omega(:)

    integer ::iunit

    call push_sub('em_resp.static_pol_lr_run')
    
    gr => sys%gr
    ndim = sys%gr%sb%dim

    call parse_input()

    if(.not.props%dynamic) then
      props%complex_response = .false.
      nfreq = 1
      nsigma = 1
    else 
      props%complex_response = (eta /= M_ZERO )
      nfreq = 1   ! for now only polarizability, so only one frequency
      nsigma = 2  ! but positive and negative values of the frequency must be considered
    end if

    !if wfs are complex then we will calculate complex_response
    props%complex_response = props%complex_response .or. (sys%st%d%wfs_type == M_CMPLX)

    call read_wfs()

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)

    call system_h_setup(sys, h)

    ALLOCATE(lr(1:ndim, 1:nsigma, 1:nfreq), ndim*nfreq*nsigma)
    
    do i = 1, ndim
      do sigma = 1, nsigma
        do j = 1, nfreq
          
          if (props%complex_response) then 
            call lr_init(lr(i, sigma, j), "Pol", def_solver=LR_BICGSTAB)
          else
            call lr_init(lr(i, sigma, j), "Pol", def_solver=LR_CG)            
          end if

          call lr_alloc_fHxc (sys%st, gr%m, lr(i, sigma, j))

          if(.not.props%complex_response ) then
            ierr = dlr_alloc_psi(sys%st, gr%m, lr(i, sigma, j))
            lr(i, sigma, j)%ddl_rho = M_ZERO
          else
            ierr = zlr_alloc_psi(sys%st, gr%m, lr(i, sigma, j))
            lr(i, sigma, j)%zdl_rho = M_ZERO
          end if

          if(props%add_fxc) then 
            call lr_build_fxc(gr%m, sys%st, sys%ks%xc, lr(i, sigma, j)%dl_Vxc)
          else 
            lr(i, sigma, j)%dl_Vxc=M_ZERO
          end if
          
        end do
      end do
    end do

    if(.not.fromScratch) then
      ! load wave-functions
      do i = 1, ndim
        do sigma = 1, nsigma
          do j = 1, nfreq
            write(dirname,'(a,i1,a,i1)') RESTART_DIR, i, "_", sigma
            call restart_read(trim(tmpdir)//dirname, sys%st, sys%gr, ierr, lr=lr(i, sigma, j))
            if(ierr.ne.0) then
              message(1) = "Could not load response wave-functions from '"//trim(tmpdir)//dirname
              call write_warning(1)
            end if
          end do
        end do
      end do
    end if
    call io_mkdir(trim(tmpdir)//RESTART_DIR)

    call info()

    if(props%dynamic) then 

      !!DYNAMIC

      if(.not. fromScratch) then
        iunit = io_open('linear/dynpols', action='read', status='old', die=.false.)
        if(iunit > 0) then
          call io_close(iunit)
        else
          fromScratch = .true.
        end if
      end if

      if(fromScratch) then 
        call io_mkdir('linear')
        iunit = io_open('linear/dynpols', action='write')

        write(iunit, '(a)')       '# Frequency dependent polarizability, real and imaginary parts are included'
        write(iunit, '(3a)')      '# Frequency units: [',  trim(units_out%energy%abbrev),'/hbar]'
        write(iunit, '(3a,i1,a)') '# Polarizability units: [',  trim(units_out%length%abbrev),'^',NDIM,']'
        write(iunit, '(8a)')      '# ', 'Frequency ', &
             'alpha_11 ' , 'Alpha_22 ' , 'alpha_33 ', 'alpha_12 ', 'alpha_13 ', 'alpha_23 '
        call io_close(iunit)
      end if

      message(1) = "Info: Calculating dynamic polarizabilities."
      call write_info(1)
           
      do i = 1, nomega

        write(message(1), '(a,f12.6,3a)') 'Info: Calculating polarizability for frequency: ', & 
             omega(i)/units_out%energy%factor, ' [',trim(units_out%energy%abbrev),']'

        call write_info(1)

        if( props%complex_response ) then 
          call zdynamic_response(sys, h, lr, props, zpol(1:ndim, 1:ndim), cmplx(omega(i), eta, PRECISION),status)
        else
          call ddynamic_response(sys, h, lr, props, zpol(1:ndim, 1:ndim), real(omega(i), PRECISION),status)
        end if

        iunit = io_open('linear/dynpols', action='write', position='append' )
        if(status%ok) then           
          !convert units 
          zpol = zpol/units_out%length%factor**NDIM
          write(iunit, '(13f12.6)') omega(i), zpol(1,1), zpol(2,2), zpol(3,3), zpol(1,2), zpol(1,3), zpol(2,3)
        else
          write(iunit, '(a,f12.6)') '#calculation didnt converge for frequency ', omega(i)
        end if
        
        call io_close(iunit)

        if(1 == nomega ) then 
          
          if( props%complex_response ) then 
            do j = 1, NDIM
              call zlr_calc_elf(sys%st,sys%gr, lr(j, 1, 1), lr(j,2,1))
              call zlr_output(sys%st, sys%gr, lr(j, 1, 1) ,"linear", j, sys%outp)
            end do
          else
            do j = 1, NDIM
              call dlr_calc_elf(sys%st,sys%gr, lr(j, 1, 1), lr(j,2,1))
              call dlr_output(sys%st, sys%gr, lr(j, 1, 1) ,"linear", j, sys%outp)
            end do
          end if
        end if
        
      end do
      
    else 

      !!! STATIC
      ALLOCATE(hpol_density(sys%NP,MAX_DIM,MAX_DIM,MAX_DIM),sys%NP*MAX_DIM**3)

      if ( .not. props%complex_response) then
        call dstatic_response(sys, h, lr(:, :, 1), props, &
           pol(1:ndim, 1:ndim), hpol(1:ndim, 1:ndim, 1:ndim),hpol_density)
        call output()
        do i = 1, ndim
          call dlr_calc_elf(sys%st,sys%gr, lr(i, 1, 1))          
          call dlr_output(sys%st, sys%gr, lr(i, 1, 1) ,"linear", i, sys%outp)
        end do

      else
      
        call zstatic_response(sys, h, lr(:, :, 1), props, &
           pol(1:ndim, 1:ndim), hpol(1:ndim, 1:ndim, 1:ndim),hpol_density)
        call output()
        do i = 1, ndim
          call zlr_calc_elf(sys%st,sys%gr, lr(i, 1, 1))
          call zlr_output(sys%st, sys%gr, lr(i, 1, 1) ,"linear", i, sys%outp)
        end do

      end if
      
      if(iand(sys%outp%what, output_pol_density).ne.0) then
        do i=1,NDIM
          do j=1,NDIM
            do k=1,NDIM
              write(fname, '(a,i1,a,i1,a,i1)') 'beta-', i, '-', j, '-', k
              call doutput_function (sys%outp%how, "linear", fname, gr%m, gr%sb, hpol_density(:,i,j,k), M_ONE, ierr)
            end do
          end do
        end do
      end if

      deallocate(hpol_density)
    end if
    
    do i = 1, ndim
      do sigma = 1, nsigma
        do j = 1, nfreq
          call lr_dealloc(lr(i, sigma, j))
        end do
      end do
    end do

    if(props%dynamic) deallocate(omega)
    deallocate(lr)
    call states_deallocate_wfns(sys%st)
    call pop_sub()
      
  contains

    ! ---------------------------------------------------------
    subroutine parse_input()
      integer(POINTER_SIZE) :: blk
      integer :: nrow
      integer :: number, j, k, ham_var
      FLOAT   :: omega_ini, omega_fin, domega

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

        props%dynamic = .true.

        nrow = loct_parse_block_n(blk)
        nomega = 0

        !count the number of frequencies
        do i= 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          nomega = nomega + number
        end do

        ALLOCATE(omega(1:nomega), nomega)
        
        !read frequencies
        j = 1
        do i = 0, nrow-1
          call loct_parse_block_int(blk, i, 0, number)
          call loct_parse_block_float(blk, i, 1, omega_ini)
          if(number > 1) then 
            call loct_parse_block_float(blk, i, 2, omega_fin)
            domega = (omega_fin - omega_ini)/(number - M_ONE)
            do k = 0, number-1
              omega(j + k) = (omega_ini + domega*k) * units_inp%energy%factor
            end do
            j = j + number
          else
            omega(j) = omega_ini * units_inp%energy%factor
            j = j + 1
          end if
        end do

        call loct_parse_block_end(blk)

        call sort(omega)
      else 
        props%dynamic = .false. 
      end if

    !%Variable PolEta
    !%Type float
    !%Default 0.0
    !%Section Linear Response::Polarizabilities
    !%Description
    !% Imaginary part of the frequency.
    !%End

    call loct_parse_float(check_inp('PolEta'), M_ZERO, eta)
    eta = eta * units_inp%energy%factor

    !%Variable PolHamiltonianVariation
    !%Type integer
    !%Default hartree+fxc
    !%Section Linear Response::Polarizabilities
    !%Description
    !% The terms are considered in the variation of the
    !% hamiltonian. V_ext is always considered. The default is to include
    !% the fxc and hartree terms. If you want to do RPA only include
    !% hartree, it will not be faster though.
    !%Option hartree 1 
    !% The variation of the hartree potential.
    !%Option fxc 2
    !% The exchange and correlation kernel, the variation of the
    !% exchange and correlation potential.
    !%End

    if(.not. h%ip_app) then 
      call loct_parse_int(check_inp('PolHamiltonianVariation'), 3, ham_var)    
      props%add_fxc = ((ham_var/2) == 1)
      props%add_hartree = (mod(ham_var, 2) == 1)
    else
      props%add_fxc = .false. 
      props%add_hartree = .false.
    end if

    call pop_sub()

  end subroutine parse_input


  ! ---------------------------------------------------------
  subroutine read_wfs()    
    !check how many wfs we have

    call push_sub('em_resp.parse_input')

    call restart_look(trim(tmpdir)//'restart_gs', sys%gr%m, kpoints, dim, nst, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    ! load wave-functions
    if ( props%complex_response ) then 
      call states_allocate_wfns(sys%st, gr%m, M_CMPLX)
    else 
      call states_allocate_wfns(sys%st, gr%m, M_REAL)
    end if
    
    call restart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)   
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    call pop_sub()

  end subroutine read_wfs


  ! ---------------------------------------------------------
  subroutine output()
    integer :: j, iunit, i, k
    FLOAT :: msp, bpar(MAX_DIM)

    call io_mkdir('linear')
    
    ! Output polarizabilty
    iunit = io_open('linear/polarizability_lr', action='write' )
    write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
           trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM
      write(iunit, '(a)') ']'

      msp = M_ZERO
      do j = 1, NDIM
        write(iunit, '(3f12.6)') pol(j, 1:NDIM) &
             / units_out%length%factor**NDIM
        msp = msp + pol(j, j)
      end do
      msp = msp / M_THREE

      write(iunit, '(a, f12.6)')  'Mean static polarizability', msp &
           / units_out%length%factor**NDIM

      call io_close(iunit)

      ! Output first hyperpolarizabilty (beta)
      iunit = io_open('linear/beta_lr', action='write')
      write(iunit, '(2a)', advance='no') '# Static hyperpolarizability tensor [', &
           trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM+2
      write(iunit, '(a)') ']'

      if (sys%st%d%nspin /= UNPOLARIZED ) then 
        write(iunit, '(a)') 'WARNING: Hyperpolarizability has not been tested for spin polarized systems'
      end if

      do i = 1, NDIM
        do j = 1, NDIM
          do k = 1, NDIM
            write(iunit,'(3i2,e20.8)') i, j, k, hpol(i, j, k)/units_out%length%factor**(5)
          end do
        end do
      end do

      if (NDIM == 3) then 

        bpar = M_ZERO
        do i = 1, NDIM
          do j = 1, NDIM
            if( i < j ) then 
              bpar(i) = bpar(i) + M_THREE*hpol(i, j, j)
            else 
              bpar(i) = bpar(i) + M_THREE*hpol(j, j, i)
            endif
          end do
        end do

        bpar = bpar / (M_FIVE * units_out%length%factor**(NDIM+2))
        write(iunit, '(a, 3f12.6,a)') 'B||', bpar(1:NDIM),&
               '  ( B||_i = 1/5 \sum_j(B_ijj+B_jij+B_jji) ) '

      endif

      call io_close(iunit)

    end subroutine output


    subroutine info()

      write(message(1),'(a)') 'Linear Reponse Polarizabilities'
      call messages_print_stress(stdout, trim(message(1)))

      if (props%complex_response) then 
        message(1) = 'Wavefunctions type: Complex'
      else
        message(1) = 'Wavefunctions type: Real'
      end if
      call write_info(1)

      if (props%add_hartree .and. props%add_fxc) then 
        message(1)='Hamiltonian variation: V_ext + hartree + fxc'
      else
        message(1)='Hamiltonian variation: V_ext'
        if (props%add_fxc) message(1)='Hamiltonian variation: V_ext + fxc'
        if (props%add_hartree) message(1)='Hamiltonian variation: V_ext + hartree'
      end if
      call write_info(1)
      
      if (props%dynamic) then 
        write(message(1),'(a,i3,a)') 'Calculating dynamic polarizability tensor for ', nomega, ' frequencies.'
      else
        message(1)='Calculating static polarizability and first static hyperpolarizability tensors.'
      end if
      call write_info(1)

      call messages_print_stress(stdout)

    end subroutine info

  end subroutine static_pol_lr_run

  subroutine lr_calc_current(st, gr, lr, lr_m)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    type(lr_t),       intent(inout) :: lr
    type(lr_t), optional, intent(inout) :: lr_m

    integer :: k, ist, ispin, idim, ndim, np

    CMPLX, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)

    ALLOCATE(lr%dl_j(gr%m%np, MAX_DIM, st%d%nspin), gr%m%np*MAX_DIM*st%d%nspin)

    np = NP
    ndim = NDIM

    ALLOCATE(   gpsi(1:np, 1:ndim), np*ndim)
    ALLOCATE(gdl_psi(1:np, 1:ndim), np*ndim)
    if(present(lr_m)) ALLOCATE(gdl_psi_m(1:np, 1:ndim), np*ndim)
   
    lr%dl_j = M_ZERO

    do ispin = 1, st%d%nspin
      do ist = 1, st%nst
        do idim = 1, st%d%dim

          call zf_gradient(gr%sb, gr%f_der, lr%zdl_psi(:,idim,ist,ispin), gdl_psi)
          call zf_gradient(gr%sb, gr%f_der, st%zpsi(:,idim,ist,ispin), gpsi)
          
          if(present(lr_m)) then               
            

            call zf_gradient(gr%sb, gr%f_der, lr_m%zdl_psi(:,idim,ist,ispin), gdl_psi_m)

            do k = 1, NDIM 
              
              lr%dl_j(1:np,k,ispin) = lr%dl_j(1:np, k, ispin) + (           &
                   + conjg(st%zpsi(1:np, idim, ist, ispin)) *       gdl_psi(1:np,k)   &
                   -       st%zpsi(1:np, idim, ist, ispin) * conjg(gdl_psi_m(1:np,k))  &
                   + conjg(lr_m%zdl_psi(1:np, idim, ist, ispin)) *     gpsi(1:np,k)   & 
                   -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np,k))  &
                   )/(M_TWO*M_zI)
            end do
            
          else 
            
            do k = 1, NDIM 
              
              lr%dl_j(1:np,k,ispin) = lr%dl_j(1:np, k, ispin) + (           &
                   + conjg(st%zpsi(1:np, idim, ist, ispin)) *       gdl_psi(1:np,k)   &
                   -       st%zpsi(1:np, idim, ist, ispin)  * conjg(gdl_psi(1:np,k))  &
                   + conjg(lr%zdl_psi(1:np, idim, ist, ispin)) *       gpsi(1:np,k)   & 
                   -       lr%zdl_psi(1:np, idim, ist, ispin)  * conjg(gpsi(1:np,k))  &
                   )/(M_TWO*M_zI)
              
            end do
            
          end if
          
        end do
      end do
    end do
    
    deallocate(gpsi)
    deallocate(gdl_psi)
    
  end subroutine lr_calc_current

#include "undef.F90"
#include "complex.F90"

#include "em_resp_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "em_resp_inc.F90"

end module static_pol_lr_m
