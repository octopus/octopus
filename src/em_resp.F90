!! Copyright (C) 2004 Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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

  implicit none

  private
  public :: &
       static_pol_lr_run

contains

  ! ---------------------------------------------------------
  subroutine static_pol_lr_run(sys, h, fromScratch)
    type(system_t), target, intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                   intent(inout) :: fromScratch

    type(lr_t), allocatable :: lr(:,:,:) ! lr(NDIM,NS,NFREQ)
    type(grid_t),   pointer :: gr

    FLOAT :: pol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim)
    FLOAT :: hpol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim, 1:sys%gr%sb%dim)

    FLOAT :: w, delta
    integer :: nfreq, nsigma, ndim, i, j, sigma, ierr
    logical :: dynamic_pol, complex_response


    call push_sub('em_resp.static_pol_lr_run')
    
    gr => sys%gr

    !FIRST WE CHECK WHAT WE WANT TO DO

    ndim=sys%gr%sb%dim
    if ( conf%devel_version ) then 
      !! Read frequency
      call loct_parse_float(check_inp('PolBaseFrequency'), M_ZERO , w)
    else
      w=M_ZERO
    endif

    if( w == M_ZERO ) then
      dynamic_pol=.false.
      complex_response=.false.
      nfreq=1
      nsigma=1
    else 
      dynamic_pol=.true.
      complex_response=.true.
      call loct_parse_float(check_inp('Delta'), M_ZERO , delta)
!      complex_response=.false.
      nfreq=1  !for now only polarizability, so only one frequency
      nsigma=2 !but positive and negative values of the frequency must be considered
    end if

    !!

    if(complex_response) then 
      ! allocate wfs
      allocate(sys%st%zpsi(gr%m%np_part, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
      ! load wave-functions
      call zrestart_read(trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    else
      allocate(sys%st%X(psi)(gr%m%np_part, sys%st%d%dim, sys%st%nst, sys%st%d%nik))
      call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%gr, ierr)
    endif

    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)

    if(complex_response) then 
      call zsystem_h_setup(sys, h)
    else
      call X(system_h_setup) (sys, h)
    endif

    !if(.not.fromScratch) then ! try to load delta_psi
    !  if(X(restart_read) (trim(tmpdir)//'restart_lr_static_pol', sys%st, m).ne.0) then
    fromScratch = .true.



    ALLOCATE(lr(1:ndim,1:nsigma,1:nfreq),ndim*nfreq*nsigma)

    do i=1,ndim
      do sigma=1,nsigma
        do j=1,nfreq

          call lr_init(lr(i,sigma,j), "SP")

          if( .not. complex_response ) then 
            call X(lr_alloc_fHxc)  (sys%st, gr%m, lr(i,sigma,j))
            ierr = X(lr_alloc_psi) (sys%st, gr%m, lr(i,sigma,j))
            lr(i,sigma,j)%X(dl_psi)=M_ZERO
          else
            call zlr_alloc_fHxc(sys%st, gr%m, lr(i,sigma,j))
            ierr = zlr_alloc_psi(sys%st, gr%m, lr(i,sigma,j))
            lr(i,sigma,j)%zdl_psi=M_ZERO
            lr(i,sigma,j)%zdl_rho=M_ZERO
          end if

          call lr_build_fxc(gr%m, sys%st, sys%ks%xc, lr(i,sigma,j)%dl_Vxc)


        end do
      end do
    enddo

    if( dynamic_pol ) then 
      print*, "CALCULATING DYNAMIC POLARIZABILITIES :-D"
      print*, "WARNING, NOT WORKING!!!!!!!!!!!!!!!"
      if( complex_response ) then 
        call zdynamic(sys, h, lr, pol, cmplx(w,delta,PRECISION))
      else 
        call X(dynamic)(sys, h, lr, pol, w)
      end if
    else 
      call static(sys, h, lr(:,:,1), pol, hpol)
      call output()
    end if


    do i=1,ndim
      do sigma=1,nsigma
        do j=1,nfreq
          call lr_dealloc(lr(i,sigma,j))
        end do
      end do
    end do

    deallocate(lr)
    if(complex_response) then 
      deallocate(sys%st%zpsi)
    else
      deallocate(sys%st%X(psi))
    endif
    call pop_sub()
      
      
  contains

    ! ---------------------------------------------------------
    subroutine output()
      integer :: j, iunit,i,k
      FLOAT :: msp, bpar(3)
      call io_mkdir('linear')

      !! Output polarizabilty
      iunit = io_open('linear/polarizability_lr', action='write')
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

      !! Output first hyperpolarizabilty (beta)
      iunit = io_open('linear/beta_lr', action='write')
      write(iunit, '(2a)', advance='no') '# Static hyperpolarizability tensor [', &
           trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM+2
      write(iunit, '(a)') ']'

      if (sys%st%d%nspin /= UNPOLARIZED ) then 
        write(iunit, '(a)') 'WARNING: Hyperpolarizability has not been tested for spin polarized systems'
      end if

      do i=1,NDIM
        do j=1,NDIM
          do k=1,NDIM
            write(iunit,'(3i2,f12.6)') i, j, k, hpol(i,j,k)/units_out%length%factor**(5)
          end do
        end do
      end do

      if (NDIM == 3) then 

        bpar=M_ZERO
        do i=1,NDIM
          do j=1,NDIM
            if( i < j ) then 
              bpar(i)=bpar(i)+M_THREE*hpol(i,j,j)
            else 
              bpar(i)=bpar(i)+M_THREE*hpol(j,j,i)
            endif
!              bpar(i)=bpar(i)+(hpol(i,j,j)+hpol(j,i,j)+hpol(j,j,i))
          end do
        end do

        bpar=bpar/(M_FIVE * units_out%length%factor**(NDIM+2))
        write(iunit, '(a, 3f12.6,a)') 'B||', bpar(1:NDIM),&
               '  ( B||_i = 1/5 \sum_j(B_ijj+B_jij+B_jji) ) '

      endif


      call io_close(iunit)

    end subroutine output

  end subroutine static_pol_lr_run

  ! ---------------------------------------------------------
  subroutine static(sys, h, lr, pol, hpol)
    type(system_t),      intent(inout) :: sys
    type(hamiltonian_t), intent(inout) :: h
    type(lr_t),          intent(inout) :: lr(:,:) ! lr(NDIM,1,1)
    FLOAT,                  intent(out)   :: pol(:,:)
    FLOAT,                  intent(out)   :: hpol(:,:,:)

    integer :: i, j, k
    integer :: ispin, ist, ispin2, ist2,n,np, dim

    R_TYPE :: prod

    R_TYPE, allocatable :: tmp(:,:), dVde(:,:,:), drhode(:,:)
    FLOAT,  allocatable :: kxc(:,:,:,:)
    FLOAT               :: spinfactor

    np  = sys%gr%m%np
    dim = sys%gr%sb%dim

    call push_sub('em_resp.static')

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

      call X(get_response_e)(sys, h, lr(:,:), i, nsigma=1, omega=R_TOTYPE(M_ZERO))
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
      do j = 1,np
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
          hpol(i,j,k)=-(hpol(i,j,k)+hpol(j,k,i)+hpol(k,i,j)+hpol(k,j,i)+hpol(j,i,k)+hpol(i,k,j))
        end do ! k
      end do ! j
    end do ! i

    deallocate(tmp)
    deallocate(dVde)
    deallocate(drhode)
    deallocate(kxc)

    call pop_sub()
  end subroutine static



#include "undef.F90"
#include "complex.F90"

#include "em_resp_inc.F90"

#include "undef.F90"
#include "real.F90"
#include "em_resp_inc.F90"



end module static_pol_lr_m
