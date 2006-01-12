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

module static_pol_lr
  use global
  use messages
  use units
  use mesh
  use mesh_function
  use system
  use restart
  use hamiltonian
  use mix
  use poisson
  use linear_response
  use io
  use lib_oct_parser

  implicit none

  private
  public :: &
       static_pol_lr_run

contains

  ! ---------------------------------------------------------
  subroutine static_pol_lr_run(sys, h, fromScratch)
    type(system_type), target, intent(inout) :: sys
    type(hamiltonian_type),    intent(inout) :: h
    logical,                   intent(inout) :: fromScratch

    type(lr_type), allocatable :: lr(:,:,:) ! lr(NDIM,NS,NFREQ)
    type(grid_type),   pointer :: gr

    FLOAT :: pol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim)
    FLOAT :: hpol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim, 1:sys%gr%sb%dim)

    FLOAT :: w 
    integer :: nfreq, nsigma, ndim, i, j, sigma, ierr
    logical :: dynamic_pol

    call init_()

    ! load wave-functions
    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, gr%m, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)

    call X(system_h_setup) (sys, h)

    !if(.not.fromScratch) then ! try to load delta_psi
    !  if(X(restart_read) (trim(tmpdir)//'restart_lr_static_pol', sys%st, m).ne.0) then
    fromScratch = .true.

    ndim=sys%gr%sb%dim

    if ( conf%devel_version ) then 
      !! Read frequency
      call loct_parse_float(check_inp('PolBaseFrequency'), M_ZERO , w)
    else
      w=M_ZERO
    endif

    if( w == M_ZERO ) then
      dynamic_pol=.false.
      nfreq=1
      nsigma=1
    else 
      dynamic_pol=.true.
      nfreq=1  !for now only polarizability, so only one frequency
      nsigma=2 !but positive and negative values of the frequency must be considered
    end if

    ALLOCATE(lr(1:ndim,1:nsigma,1:nfreq),ndim*nfreq*nsigma)

    do i=1,ndim
      do sigma=1,nsigma
        do j=1,nfreq

          call lr_init(lr(i,sigma,j), "SP")

          call X(lr_alloc_fHxc)  (sys%st, gr%m, lr(i,sigma,j))
          ierr = X(lr_alloc_psi) (sys%st, gr%m, lr(i,sigma,j))
      
          call lr_build_fxc(gr%m, sys%st, sys%ks%xc, lr(i,sigma,j)%dl_Vxc)

        end do
      end do
    enddo

    if( dynamic_pol ) then 
      print*, "CALCULATING DYNAMIC POLARIZABILITIES :-D"
      print*, "WARNING, NOT TESTED"
      call dynamic(sys, h, lr, pol, w)
    else 
      call static(sys, h, lr(:,:,1), pol, hpol)
    end if

    call output()

    do i=1,ndim
      do sigma=1,nsigma
        do j=1,nfreq
          call lr_dealloc(lr(i,sigma,j))
        end do
      end do
    end do

    call end_()
    
    deallocate(lr)

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol_lr.static_pol_lr_run')

      gr => sys%gr

      ! allocate wfs
      allocate(sys%st%X(psi)(gr%m%np_part, sys%st%d%dim, sys%st%nst, sys%st%d%nik))

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      deallocate(sys%st%X(psi))

      call pop_sub()
    end subroutine end_


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

      do k=1,NDIM
        do j=1,k
          do i=1,j
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
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr(:,:) ! lr(NDIM,1,1)
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

    call push_sub('static_pol_lr.static')

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

      call get_response_e(sys, h, lr(:,:), i, nsigma=1, omega=R_TOTYPE(M_ZERO))
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

    hpol = -M_SIX*hpol

    do k = 1, dim
      do j = 1, k
        do i = 1, j
!          print*,i,j,k
          hpol(i,j,k)=(hpol(i,j,k)+hpol(j,k,i)+hpol(k,i,j)+hpol(k,j,i)+hpol(j,i,k)+hpol(i,k,j))/M_SIX
        end do ! k
      end do ! j
    end do ! i

    deallocate(tmp)
    deallocate(dVde)
    deallocate(drhode)
    deallocate(kxc)

    call pop_sub()
  end subroutine static


  ! ---------------------------------------------------------
  subroutine dynamic(sys, h, lr, pol, inw)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr(:,:,:) ! dim, nsigma(=2), nfreq(=1)
    FLOAT,                  intent(out)   :: pol(:,:)
    FLOAT,                  intent(in)    :: inw 

    integer :: dir, j, freq, sigma

    FLOAT :: rhov
    R_TYPE :: w

    call push_sub('static_pol_lr.dynamic')

    freq = 1
    w    = R_TOTYPE(inw)

    pol = M_ZERO

    do dir = 1, sys%gr%sb%dim
      write(message(1), '(a,i1,a,f12.6)') 'Info: Calculating polarizability for direction ', &
         dir, ' and frequency ', w

      call write_info(1)
      
      do sigma = 1, 2
        call mix_init(lr(dir,sigma,freq)%mixer, sys%gr%m, 1, sys%st%d%nspin)
      end do
      
      call get_response_e(sys, h, lr(:,:,freq), dir, nsigma=2 , omega=w)
      
      do sigma = 1, 2
        call mix_end(lr(dir,sigma,freq)%mixer)
      end do
      
      do sigma = 1, 2
        do j = 1, sys%gr%m%np
          rhov = sum(lr(dir,sigma,freq)%X(dl_rho)(j,1:sys%st%d%nspin))*sys%gr%m%vol_pp(j)
          pol(dir, :) = pol(dir, :) - M_HALF*sys%gr%m%x(j,:)*rhov
        end do
      end do
    end do
    
    print*, pol

    call pop_sub()

  end subroutine dynamic


  ! ---------------------------------------------------------
  subroutine get_response_e(sys, h, lr, dir, nsigma, omega)
    type(system_type), target, intent(inout) :: sys
    type(hamiltonian_type),    intent(inout) :: h
    type(lr_type),             intent(inout) :: lr(:,:) !ndim,nsigma
    integer,                   intent(in)    :: dir 
    integer,                   intent(in)    :: nsigma 
    R_TYPE,                    intent(in)    :: omega
    FLOAT, allocatable :: diff(:,:,:)
    FLOAT :: dpsimod,freq_sign
    integer :: iter, sigma, ik, ik2, ist, i
    FLOAT, allocatable :: tmp(:), dl_rhoin(:,:,:,:), dl_rhonew(:,:,:,:), dl_rhotmp(:,:,:,:)
    R_TYPE, allocatable :: Y(:, :),dV(:,:)
    R_TYPE, allocatable :: a(:,:)
    logical :: finish(2)

    type(mesh_type), pointer :: m
    type(states_type), pointer :: st

    call push_sub('static_pol_lr.get_response_e')

    ASSERT( nsigma==1 .or. nsigma ==2 )

    m => sys%gr%m
    st => sys%st

    ALLOCATE(tmp(m%np),m%np)
    ALLOCATE(Y(m%np,1),m%np*1)
    ALLOCATE(dV(m%np,st%d%nspin),m%np*st%d%nspin)
    ALLOCATE(dl_rhoin(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
    ALLOCATE(dl_rhonew(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
    ALLOCATE(dl_rhotmp(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
    ALLOCATE(diff(st%nst,st%d%nspin,nsigma),st%nst*st%d%nspin*nsigma)
    ALLOCATE(a(st%nst,st%d%nspin),st%nst*st%d%nspin)

    diff=M_ZERO

    call init_response_e()
    finish = .false.
    iter_loop: do iter=1, lr(dir,1)%max_iter
      
      do sigma=1,nsigma
        dl_rhoin(1,:,:,sigma) = lr(dir,sigma)%X(dl_rho)(:,:)
      end do

      if(.not.h%ip_app) then
        do sigma=1,nsigma
          do i = 1, m%np
            tmp(i) = sum(lr(dir,sigma)%X(dl_rho)(i,:))
          end do
          call dpoisson_solve(sys%gr, lr(dir,sigma)%ddl_Vhar, tmp)
        end do
      end if


      do ik = 1, st%d%nspin
        do sigma = 1, nsigma
          
          if(finish(sigma)) cycle !if we are ready with this sign, do nothing

          if(sigma==1) then 
            freq_sign =  M_ONE
          else
            freq_sign = -M_ONE
          end if

          dV(1:m%np,ik) = (lr(dir,sigma)%ddl_Vhar(1:m%np) + m%x(1:m%np,dir))
          
          do ik2 = 1, st%d%nspin
            dV(1:m%np,ik) = dV(1:m%np,ik) +&
                 lr(dir,sigma)%dl_Vxc(1:m%np, ik, ik2)*lr(dir,sigma)%X(dl_rho)(1:m%np,ik2)
          end do
          
          do ist = 1, st%nst
            if(st%occ(ist, ik) <= M_ZERO) cycle
            
            Y(1:m%np,1) = -dV(1:m%np,ik)*st%X(psi)(1:m%np, 1, ist, ik)

            
            call X(lr_orth_vector)(m, st, Y, ik)

            call X(lr_solve_HXeY) (lr(dir,1), h, sys%gr, sys%st%d, ik, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), Y, &
                 -sys%st%eigenval(ist, ik) + freq_sign*omega)

            !altough dl_psi should be orthogonal to psi
            !a re-orthogonalization is sometimes necessary 
            call X(lr_orth_vector)(m, st, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), ik)

            dpsimod = sum(lr(dir,sigma)%X(dl_psi)(1:m%np, 1, ist, ik)**2 * sys%gr%m%vol_pp(1:m%np))

            print*, iter, (3-2*sigma)*ist, dpsimod, (dpsimod-diff(ist,ik,sigma))

            diff(ist,ik,sigma)=dpsimod

          end do
        end do
      end do

      ! calculate dl_rho

      if(nsigma == 2 ) then
        call build_rho_dynamic()
      else ! static case
        lr(dir, 1)%X(dl_rho) = M_ZERO
        call X(lr_build_dl_rho)(m, st, lr(dir,1), type=3)
      end if
      
      ! mix to get new density
      dl_rhonew(1,:,:,:) = M_ZERO

      finish=.true.

      do sigma = 1, nsigma
        dl_rhotmp(1,:,:,sigma) = lr(dir,sigma)%X(dl_rho)(:,:)

        call mixing(lr(dir,sigma)%mixer, m, iter, 1, st%d%nspin, &
             dl_rhoin(:,:,:,sigma), dl_rhotmp(:,:,:,sigma), dl_rhonew(:,:,:,sigma))

        ! check for convergence
        lr(dir, sigma)%abs_dens = M_ZERO
      
        do ik = 1, st%d%nspin
          tmp(:) = (dl_rhoin(1,:,ik,sigma) - lr(dir,sigma)%X(dl_rho)(:,ik))**2
          lr(dir,sigma)%abs_dens = lr(dir,sigma)%abs_dens + dmf_integrate(m, tmp)
        end do
        lr(dir,sigma)%abs_dens = sqrt(lr(dir,sigma)%abs_dens)
        print*, lr(dir,sigma)%abs_dens
        ! are we finished?
        finish(sigma) = (lr(dir,sigma)%abs_dens <= lr(dir,sigma)%conv_abs_dens)
        lr(dir,sigma)%abs_dens=M_ZERO
      end do

      if( finish(1) .and. finish(2) ) then
        write(message(1), '(a, i4, a)')        &
             'Info: SCF for response converged in ', &
             iter, ' iterations'
        call write_info(1)
        exit
      else
        do sigma=1,nsigma
          lr(dir,sigma)%X(dl_rho)(:,:) = dl_rhonew(1,:,:,sigma)
        end do
      end if
        
    end do iter_loop

    deallocate(dl_rhoin, dl_rhonew, dl_rhotmp)
    deallocate(tmp, Y,dV)
    deallocate(diff,a)
    call pop_sub()

  contains

    !------------------------------------------------------------
    subroutine init_response_e()
      integer :: ik, ist, i
      FLOAT :: rd

      call push_sub('static_pol_lr.init_response_e')

      do ik = 1, st%d%nspin
        do ist = 1, st%nst
          if (st%occ(ist, ik) > M_ZERO) then
            do i = 1, m%np
              call mesh_r(m, i, rd)
              do sigma = 1, nsigma
                lr(dir, sigma)%X(dl_psi)(i, 1, ist, ik) = st%X(psi)(i, 1, ist, ik)*rd*exp(-rd)
              end do
            end do
          end if
        end do
      end do

      do sigma=1,nsigma
        lr(dir,sigma)%ddl_Vhar(:) = M_ZERO
        call lr_orth_response(m, st, lr(dir,sigma))
      end do

      if(nsigma == 2 ) then 
        call build_rho_dynamic() 
      else
        call X(lr_build_dl_rho)(m, st, lr(dir,1), 3)
      end if

      call pop_sub()
    end subroutine init_response_e

    !------------------------------------------------------------
    subroutine build_rho_dynamic()
      integer :: msigma

      !!the dinamic case
      do sigma=1,nsigma 
        lr(dir,sigma)%X(dl_rho)=M_ZERO
        
        !!the oposite sign
        if (sigma==1) msigma=2
        if (sigma==2) msigma=1

        do ik = 1, st%d%nspin
          do ist = 1, st%nst
            lr(dir,sigma)%X(dl_rho)(1:m%np,ik)=lr(dir,sigma)%X(dl_rho)(1:m%np,ik)+&
                 st%d%kweights(ik)*st%occ(ist, ik) * (&
                 R_CONJ(st%X(psi)(1:m%np, 1, ist, ik))*lr(dir,sigma)%X(dl_psi)(1:m%np,1,ist,ik) +&
                 R_CONJ(lr(dir,msigma)%X(dl_psi)(1:m%np,1,ist,ik))*st%X(psi)(1:m%np, 1, ist, ik)&
                 )
          end do
        end do

      end do
        
    end  subroutine build_rho_dynamic
    
  end subroutine get_response_e
  
end module static_pol_lr
