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
  
  type pol_props 
     logical :: complex_response
     logical :: add_fxc
     logical :: use_unoccupied
     logical :: dynamic
     logical :: ort_each_step
  end type pol_props

  
contains

  ! ---------------------------------------------------------
  subroutine static_pol_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(lr_t), allocatable :: lr(:,:,:) ! lr(NDIM,NS,NFREQ)
    type(grid_t),   pointer :: gr

    FLOAT ::  pol(1:MAX_DIM, 1:MAX_DIM)
    CMPLX :: zpol(1:MAX_DIM, 1:MAX_DIM)
    FLOAT :: hpol(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM)

    FLOAT :: delta
    integer :: nfreq, nsigma, ndim, i, j, sigma, ierr, kpoints, dim, nst
    type(pol_props) :: props
    
    integer :: nomega
    FLOAT, allocatable :: omega(:)

    integer ::iunit

    call push_sub('em_resp.static_pol_lr_run')
    
    gr => sys%gr
    ndim = sys%gr%sb%dim
    call parse_freq_blk()

    call loct_parse_logical(check_inp('PolAddFXC'), .true., props%add_fxc)
    call loct_parse_logical(check_inp('PolOrtEachStep'), .false., props%ort_each_step)

    if(.not.props%dynamic) then
      props%complex_response = .false.
      nfreq = 1
      nsigma = 1
    else 
      call loct_parse_float(check_inp('Delta'), M_ZERO, delta)
      props%complex_response = (delta /= M_ZERO ) 
      nfreq = 1   ! for now only polarizability, so only one frequency
      nsigma = 2  ! but positive and negative values of the frequency must be considered
    end if

    call read_wfs()

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    
    if(props%complex_response) then 
      call zsystem_h_setup(sys, h)
    else
      call X(system_h_setup) (sys, h)
    endif
    
    fromScratch = .true.
    
    ALLOCATE(lr(1:ndim,1:nsigma,1:nfreq), ndim*nfreq*nsigma)
    
    do i = 1, ndim
      do sigma = 1, nsigma
        do j = 1, nfreq

          call lr_init(lr(i, sigma, j), "SP")

          if(.not.props%complex_response) then 
            call X(lr_alloc_fHxc)  (sys%st, gr%m, lr(i, sigma, j))
            ierr = X(lr_alloc_psi) (sys%st, gr%m, lr(i, sigma, j))
            lr(i, sigma, j)%X(dl_psi) = M_ZERO
          else
            call zlr_alloc_fHxc(sys%st, gr%m, lr(i, sigma, j))
            ierr = zlr_alloc_psi(sys%st, gr%m, lr(i, sigma, j))
            lr(i, sigma, j)%zdl_psi = M_ZERO
            lr(i, sigma, j)%zdl_rho = M_ZERO
          end if

          if(props%add_fxc) then 
            call lr_build_fxc(gr%m, sys%st, sys%ks%xc, lr(i, sigma, j)%dl_Vxc)
          else 
            lr(i, sigma, j)%dl_Vxc=M_ZERO
          end if

        end do
      end do
    enddo

    if(props%dynamic) then 
      
      call io_mkdir('linear')
      iunit = io_open('linear/dynpols', action='write')

      ! todo: write header
      write(iunit, '(8a)')  '# ', 'freq ', '11 ' , '22 ' , '33 ', '12 ', '13 ', '23 '
      call io_close(iunit)

      ! the dynamic case
      message(1) = "Info: Calculating dynamic polarizabilities."
      call write_info(1)
           
      do i= 1, nomega

        write(message(1), '(a,f12.6)') 'Info: Calculating polarizability for frequency: ', omega(i)
        call write_info(1)

        if( props%complex_response ) then 
          call zdynamic_response(sys, h, lr, props, zpol(1:ndim, 1:ndim), cmplx(omega(i), delta, PRECISION))
        else 
          call X(dynamic_response)(sys, h, lr, props, zpol(1:ndim, 1:ndim), R_TOTYPE(omega(i)))
        end if
        iunit = io_open('linear/dynpols', action='write', position='append' )
        write(iunit, '(13f12.6)') omega(i), zpol(1,1), zpol(2,2), zpol(3,3), zpol(1,2), zpol(1,3), zpol(2,3)
        call io_close(iunit)
      end do

    else 

      ! the static case
      call static_response(sys, h, lr(:, :, 1), props, &
           pol(1:ndim, 1:ndim), hpol(1:ndim, 1:ndim, 1:ndim))
      call output()
      do i = 1,ndim
        call X(lr_output) (sys%st, sys%gr, lr(i, 1, 1) ,"linear", i, sys%outp)
      end do
    end if
    
    do i = 1, ndim
      do sigma = 1, nsigma
        do j= 1, nfreq
          call lr_dealloc(lr(i, sigma, j))
        end do
      end do
    end do

    if(props%dynamic) deallocate(omega)
    deallocate(lr)
    if(props%complex_response) then 
      deallocate(sys%st%zpsi)
    else
      deallocate(sys%st%X(psi))
    endif
    call pop_sub()
      
  contains

    ! ---------------------------------------------------------
    subroutine parse_freq_blk()
      
      integer(POINTER_SIZE) :: blk
      integer :: nrow
      integer :: number, j, k
      FLOAT   :: omega_ini, omega_fin, domega
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
            domega = (omega_fin-omega_ini)/(number-M_ONE)
            do k = 0, number-1
              omega(j+k) = omega_ini + domega*k
            end do
            j = j + number
          else
            omega(j) = omega_ini
            j = j + 1
          end if
        end do

        call loct_parse_block_end(blk)

        call sort(omega)
    else 
      props%dynamic = .false. 
    end if

  end subroutine parse_freq_blk


  ! ---------------------------------------------------------
  subroutine read_wfs()    
    !check how many wfs we have
    
    call restart_look(trim(tmpdir)//'restart_gs', sys%gr%m, kpoints, dim, nst, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    !if there are unoccupied wfs, read them
    if (nst > sys%st%nst) then 
      
      call loct_parse_logical(check_inp('LRUseUnoccupied'), .false. , props%use_unoccupied)

      if (props%use_unoccupied) then 
      
        write(message(1), '(a,i2,a)') 'Info: Found ', (nst - sys%st%nst), ' unoccupied wavefunctions.'
        call write_info(1)
        
        sys%st%nst    = nst
        sys%st%st_end = nst
        deallocate(sys%st%eigenval, sys%st%occ)
        
        ALLOCATE(sys%st%eigenval(sys%st%nst, sys%st%d%nik), sys%st%nst*sys%st%d%nik)
        ALLOCATE(     sys%st%occ(sys%st%nst, sys%st%d%nik), sys%st%nst*sys%st%d%nik)
      else
        write(message(1), '(a)') 'Info: Not using unoccupied wavefunctions.'
        call write_info(1)
      end if

    end if

    ! load wave-functions
    if(props%complex_response) then 
      call zstates_allocate_wfns(sys%st, gr%m)
      call zrestart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    else
      call X(states_allocate_wfns)(sys%st, gr%m)
      call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    endif
    
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if
    
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
        msp = msp + pol(j,j)
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
            write(iunit,'(3i2,f12.6)') i, j, k, hpol(i,j,k)/units_out%length%factor**(5)
          end do
        end do
      end do

      if (NDIM == 3) then 

        bpar = M_ZERO
        do i = 1, NDIM
          do j = 1, NDIM
            if( i < j ) then 
              bpar(i) = bpar(i) + M_THREE*hpol(i,j,j)
            else 
              bpar(i) = bpar(i) + M_THREE*hpol(j,j,i)
            endif
!              bpar(i)=bpar(i)+(hpol(i,j,j)+hpol(j,i,j)+hpol(j,j,i))
          end do
        end do

        bpar = bpar / (M_FIVE * units_out%length%factor**(NDIM+2))
        write(iunit, '(a, 3f12.6,a)') 'B||', bpar(1:NDIM),&
               '  ( B||_i = 1/5 \sum_j(B_ijj+B_jij+B_jji) ) '

      endif

      call io_close(iunit)

    end subroutine output

  end subroutine static_pol_lr_run


  ! ---------------------------------------------------------
  subroutine static_response(sys, h, lr, props, pol, hpol)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(lr_t),          intent(inout) :: lr(:,:) ! lr(NDIM,1,1)
    type(pol_props),     intent(in)    :: props
    FLOAT,               intent(out)   :: pol(:,:)
    FLOAT,               intent(out)   :: hpol(:,:,:)

    integer :: i, j, k
    integer :: ispin, ist, ispin2, ist2, n, np, dim

    R_TYPE :: prod

    R_TYPE, allocatable :: tmp(:,:), dVde(:,:,:), drhode(:,:)
    FLOAT,  allocatable :: kxc(:,:,:,:)
    FLOAT               :: spinfactor

    np  = sys%gr%m%np
    dim = sys%gr%sb%dim

    call push_sub('em_resp.static_response')

    message(1) = "Info:  Calculating static properties"
    call write_info(1)

    ALLOCATE(tmp(1:np, 1), np)
    ALLOCATE(dVde(1:np, 1:sys%st%d%nspin, 1:dim), np*sys%st%d%nspin*dim)
    ALLOCATE(drhode(1:np, 1:dim), np*dim)
    ALLOCATE(kxc(1:np, 1:sys%st%d%nspin, 1:sys%st%d%nspin, 1:sys%st%d%nspin), np*sys%st%d%nspin**3)

    call xc_get_kxc(sys%ks%xc, sys%gr%m, sys%st%rho, sys%st%d%ispin, kxc)

    ! first the derivatives in all directions are calculated and stored
    write(message(1), '(a)') 'Info: Calculating derivatives of the orbitals:'
    call write_info(1)

    do i = 1, dim
      write(message(1), '(a,i1)') 'Info: Derivative direction: ', i
      call write_info(1)

      call mix_init(lr(i,1)%mixer, sys%gr%m, 1, sys%st%d%nspin)

      call X(get_response_e)(sys, h, lr(:,:), i, 1, R_TOTYPE(M_ZERO), props)
      call mix_end(lr(i, 1)%mixer)

      do ispin = 1, sys%st%d%nspin
        ! the potential derivatives

        ! Hartree and the potential and the derivative of the
        !  potential associated at the electric field

        dVde(1:np,ispin, i) = lr(i,1)%X(dl_Vhar)(1:np) + sys%gr%m%x(1:np, i)

        ! xc
        do ispin2 = 1, sys%st%d%nspin
          dVde(1:np,ispin,i) = dVde(1:np,ispin,i) + &
             lr(i,1)%dl_Vxc(1:np, ispin, ispin2)*lr(i,1)%X(dl_rho)(1:np, ispin2)
        end do
      end do

      ! the density
      do n = 1, np
        drhode(n, i) = sum(lr(i,1)%X(dl_rho)(n, 1:sys%st%d%nspin))
      end do
      
    end do

    write(message(1), '(a)') 'Info: Calculating polarizability tensor'
    call write_info(1)

    pol = M_ZERO
    do i = 1, dim
      do j = 1, np
        pol(i, 1:dim) = pol(i, 1:dim) - &
           sys%gr%m%x(j, 1:dim) * drhode(j, i) * sys%gr%m%vol_pp(j)
      end do
    end do

    write(message(1), '(a)') 'Info: Calculating hyperpolarizability tensor'
    call write_info(1)
    
    if (sys%st%d%nspin /= UNPOLARIZED ) then 
      write(message(1), '(a)') 'WARNING: Hyperpolarizability has not been tested for spin polarized systems'
      call write_warning(1)
      spinfactor = M_ONE
    else 
      spinfactor = M_TWO
    end if

    hpol = M_ZERO

    do i = 1, dim
      do j = 1, dim
        do k = 1, dim
          
          do ispin = 1, sys%st%d%nspin
            do ist = 1, sys%st%nst

              ! <D\psi_n | P_c DV_scf P_c | D\psi_n >

              tmp(1:np, 1)  = dVde(1:np, ispin, j) * lr(k,1)%X(dl_psi)(1:np,1,ist,ispin)
              hpol(i, j, k) = hpol(i, j, k) + &
                 spinfactor * sum(R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * &
                 tmp(1:np, 1) * sys%gr%m%vol_pp(1:np))

              do ispin2 = 1, sys%st%d%nspin
                do ist2 = 1, sys%st%nst

                  prod = sum(R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * &
                     lr(k, 1)%X(dl_psi)(1:np, 1, ist2, ispin2) * sys%gr%m%vol_pp(1:np))

                  prod = prod * sum( &
                     R_CONJ(sys%st%X(psi)(1:np, 1, ist2, ispin2)) * &
                     dVde(1:np, ispin, j) * sys%st%X(psi)(1:np, 1, ist, ispin) * &
                     sys%gr%m%vol_pp(1:np))

                  hpol(i, j, k) = hpol(i, j, k) - spinfactor*prod
                end do ! ist2
              end do ! ispin2

            end do ! ist
          end do ! ispin

          hpol(i, j, k) = hpol(i, j, k) + &
             sum(kxc(1:np, 1, 1, 1) * drhode(1:np, i) * drhode(1:np, j)*drhode(1:np, k) * &
             sys%gr%m%vol_pp(1:np))/CNST(6.0)

        end do ! k
      end do ! j
    end do ! i

    do k = 1, dim
      do j = 1, k
        do i = 1, j
          hpol(i,j,k) = -(hpol(i,j,k) + hpol(j,k,i) + hpol(k,i,j) + hpol(k,j,i) + hpol(j,i,k) + hpol(i,k,j))
        end do ! k
      end do ! j
    end do ! i

    deallocate(tmp)
    deallocate(dVde)
    deallocate(drhode)
    deallocate(kxc)

    call pop_sub()
  end subroutine static_response


#include "undef.F90"
#include "complex.F90"

#include "em_resp_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "em_resp_inc.F90"

end module static_pol_lr_m
