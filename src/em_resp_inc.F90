!! Copyright (C) 2005 M. Marques, X. Andrade
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


! ---------------------------------------------------------
subroutine X(dynamic_response)(sys, h, lr, props, pol, w, status)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:,:) ! dim, nsigma(=2), nfreq(=1)
  type(pol_props_t),     intent(in)    :: props
  CMPLX,               intent(inout) :: pol(:,:)
  R_TYPE,              intent(in)    :: w 
  type(status_t),      intent(out)   :: status
  
  integer :: dir, j, freq, sigma

  R_TYPE :: rhov
  R_TYPE :: X(pol)(MAX_DIM, MAX_DIM)

  call push_sub('static_pol_lr.dynamic')

  freq = 1
  X(pol) = M_ZERO

  do dir = 1, sys%gr%sb%dim
    write(message(1), '(a,i1)') 'Info: Calculating polarizability for direction ', dir

    call write_info(1)

    do sigma = 1, 2
      call mix_init(lr(dir,sigma,freq)%mixer, sys%gr%m, 1, sys%st%d%nspin)
    end do

    call X(get_response_e)(sys, h, lr(:,:,freq), dir, 2 , w, props, status)

    do sigma = 1, 2
      call mix_end(lr(dir,sigma,freq)%mixer)
    end do

    do j = 1, sys%gr%m%np
      rhov = sum(lr(dir,2,freq)%X(dl_rho)(j,1:sys%st%d%nspin))*sys%gr%m%vol_pp(j)
      X(pol)(dir, :) = X(pol)(dir, :) - sys%gr%m%x(j,:)*rhov
    end do

  end do

  pol=X(pol)

  call pop_sub()

end subroutine X(dynamic_response)


! ---------------------------------------------------------
subroutine X(get_response_e)(sys, h, lr, dir, nsigma, omega, props, status)
  type(system_t), target, intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: h
  type(lr_t),             intent(inout) :: lr(:,:) !ndim,nsigma
  integer,                intent(in)    :: dir 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(in)    :: omega
  type(pol_props_t),      intent(in)    :: props
  type(status_t),optional,intent(out)   :: status

  FLOAT, allocatable :: diff(:,:,:)
  FLOAT :: dpsimod,freq_sign
  integer :: iter, sigma, ik, ik2, ist, i, ist2
  FLOAT, allocatable :: dl_rhoin(:,:,:,:), dl_rhonew(:,:,:,:), dl_rhotmp(:,:,:,:), dtmp(:)
  R_TYPE, allocatable :: Y(:, :),dV(:,:), tmp(:), Hp(:,:)
  R_TYPE, allocatable :: a(:,:)
  logical :: finish(2)

  type(mesh_t), pointer :: m
  type(states_t), pointer :: st
  
  call push_sub('static_pol_lr.get_response_e')

  ASSERT( nsigma==1 .or. nsigma ==2 )
  
  m => sys%gr%m
  st => sys%st
  
  ALLOCATE(tmp(m%np),m%np)
  ALLOCATE(dtmp(m%np),m%np)
  ALLOCATE(Y(m%np,1),m%np*1)
  ALLOCATE(Hp(m%np,1),m%np*1)
  ALLOCATE(dV(m%np,st%d%nspin),m%np*st%d%nspin)
  ALLOCATE(dl_rhoin(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(dl_rhonew(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(dl_rhotmp(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(diff(st%nst,st%d%nspin,nsigma),st%nst*st%d%nspin*nsigma)
  ALLOCATE(a(st%nst,st%d%nspin),st%nst*st%d%nspin)

  diff=M_ZERO
!  print*, "OMEGA",  omega
  call init_response_e()
  finish = .false.
  iter_loop: do iter=1, lr(dir,1)%max_iter

!    do sigma=1,nsigma
!      dl_rhoin(1,:,:,sigma) = lr(dir,sigma)%X(dl_rho)(:,:)
!    end do
    !we will use the two mixers, one for the real part and one for the complex

    dl_rhoin(1,:,:,1) = R_REAL(lr(dir,1)%X(dl_rho)(:,:))

    if(nsigma==2) dl_rhoin(1,:,:,2) = R_AIMAG(lr(dir,1)%X(dl_rho)(:,:))

    if (props%add_hartree .and. (.not. h%ip_app) ) then

      do sigma=1,nsigma
        do i = 1, m%np
          tmp(i) = sum(lr(dir,sigma)%X(dl_rho)(i,:))
        end do
        call X(poisson_solve)(sys%gr, lr(dir,sigma)%X(dl_Vhar), tmp)
      end do
    end if
    
    do ik = 1, st%d%nspin
      do sigma = 1, nsigma

        if(sigma==1) then 
          freq_sign =  M_ONE
        else
          freq_sign = -M_ONE
        end if
        
        dV(1:m%np,ik) = m%x(1:m%np,dir)

        if (props%add_hartree) dV(1:m%np,ik) = dV(1:m%np,ik) + lr(dir,sigma)%X(dl_Vhar)(1:m%np) 

        if(props%add_fxc) then 
          do ik2 = 1, st%d%nspin
            dV(1:m%np,ik) = dV(1:m%np,ik) +&
                 lr(dir,sigma)%dl_Vxc(1:m%np, ik, ik2)*lr(dir,sigma)%X(dl_rho)(1:m%np,ik2)
          end do
          
        end if

        do ist = 1, st%nst
          if(st%occ(ist, ik) <= lr_min_occ) cycle
          Y(1:m%np,1) = -dV(1:m%np,ik)*st%X(psi)(1:m%np, 1, ist, ik)
          
          !          if (nsigma == 2 ) then
!          do ist2 = 1, st%nst
!            if(st%occ(ist, ik) > min_occ) then 
!              a(ist2,ik)=sum(R_CONJ(st%X(psi)(1:m%np, 1, ist2, ik))*Y(1:m%np ,1)*&
!                   sys%gr%m%vol_pp(1:m%np))
!              Y(1:m%np,1)=Y(1:m%np,1)-a(ist2,ik)*st%X(psi)(1:m%np, 1, ist2, ik)
!            end if
            !              print*,sum(Y(1:m%np,1)*st%X(psi)(1:m%np, 1, ist2, ik))
!          end do
            
          !          else
          
          call X(lr_orth_vector)(m, st, Y, ik)
          
          !          do ist2 = 1, st%nst
          !            print*,sum(R_CONJ(Y(1:m%np,1))*st%X(psi)(1:m%np, 1, ist2, ik))
          !          end do

          if(props%ort_each_step) then 
            call X(lr_solve_HXeY) (lr(dir,1), h, sys%gr, sys%st%d, ik, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), Y, &
                 -sys%st%eigenval(ist, ik) + freq_sign*omega, sys%st)
          else
            call X(lr_solve_HXeY) (lr(dir,1), h, sys%gr, sys%st%d, ik, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), Y, &
                 -sys%st%eigenval(ist, ik) + freq_sign*omega)
          end if

!          call X(Hpsi)(h, sys%gr, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik) , Hp, ik)
!          Hp(1:m%np,1)=Hp(1:m%np,1)+&
!               (-sys%st%eigenval(ist, ik)+freq_sign*omega)*lr(dir,sigma)%X(dl_psi)(1:m%np,1, ist, ik)&
!               -Y(1:m%np,1)
         
!          print*,"RES", sum(R_REAL(R_CONJ(Hp(1:m%np,1))*Hp(1:m%np,1)))

!          do ist2 = 1, st%nst
!            print*,sum(R_CONJ(lr(dir,sigma)%X(dl_psi)(1:m%np,1, ist, ik))*st%X(psi)(1:m%np, 1, ist2, ik))
!          end do

          !altough dl_psi should be orthogonal to psi
          !a re-orthogonalization is sometimes necessary 
          call X(lr_orth_vector)(m, st, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), ik)

!          if (nsigma == 2 ) then
            
!            do ist2 = 1, st%nst
!              if(sys%st%occ(ist2,ik) == M_ZERO ) then 
!                lr(dir,sigma)%X(dl_psi)(1:m%np,1, ist, ik)=lr(dir,sigma)%X(dl_psi)(1:m%np, 1, ist, ik)&
!                     +a(ist2,ik)/(sys%st%eigenval(ist2, ik)-sys%st%eigenval(ist, ik)+freq_sign*omega)*&
!                     st%X(psi)(1:m%np, 1, ist2, ik)
!              end if
!            end do
            
!          endif
         
          dpsimod = sum(R_ABS(lr(dir,sigma)%X(dl_psi)(1:m%np, 1, ist, ik))**2 * sys%gr%m%vol_pp(1:m%np))
          
          print*, iter, (3-2*sigma)*ist, dpsimod, (dpsimod-diff(ist,ik,sigma))

          diff(ist,ik,sigma)=dpsimod

        end do
      end do
    end do

    if( lr(dir,1)%max_iter == iter  ) then 
      message(1) = "Self-consistent iteration for response did not converge"
      if(present(status)) then 
        status%ok = .false. 
        call write_warning(1);
      else
        call write_fatal(1);
      end if 
    end if
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

    !we will use the two mixers, one for the real part and one for the complex
    dl_rhotmp(1,:,:,1) = R_REAL(lr(dir,1)%X(dl_rho)(:,:))
    if(nsigma==2) dl_rhotmp(1,:,:,2) = R_AIMAG(lr(dir,1)%X(dl_rho)(:,:))


    do sigma = 1, nsigma

      call mixing(lr(dir,sigma)%mixer, m, iter, 1, st%d%nspin, &
           dl_rhoin(:,:,:,sigma), dl_rhotmp(:,:,:,sigma), dl_rhonew(:,:,:,sigma))


      ! check for convergence
      lr(dir, sigma)%abs_dens = M_ZERO

      do ik = 1, st%d%nspin
        dtmp(:) = (dl_rhoin(1,:,ik,sigma) - dl_rhotmp(1,:,ik,sigma))**2
        lr(dir,sigma)%abs_dens = lr(dir,sigma)%abs_dens + dmf_integrate(m, dtmp)
      end do

      lr(dir,sigma)%abs_dens = sqrt(lr(dir,sigma)%abs_dens)
      print*, "SCF RES", sigma, lr(dir,sigma)%abs_dens
      ! are we finished?
      finish(sigma) = (lr(dir,sigma)%abs_dens <= lr(dir,sigma)%conv_abs_dens)
      lr(dir,sigma)%abs_dens=M_ZERO
    end do

    if( finish(1) .and. finish(2) ) then 
      if(present(status)) status%ok = .true.
      write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, ' iterations'
      call write_info(1)
      exit
    else

#ifdef R_TCOMPLEX
        lr(dir,1)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),dl_rhonew(1,:,:,2),PRECISION)
        lr(dir,2)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),-dl_rhonew(1,:,:,2),PRECISION)
#else
      do sigma=1,nsigma
        lr(dir,sigma)%X(dl_rho)(:,:) = dl_rhonew(1,:,:,1)
      end do
#endif
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
              !                lr(dir, sigma)%X(dl_psi)(i, 1, ist, ik) = M_ZERO
            end do
          end do
        end if
      end do
    end do

    do sigma=1,nsigma
      lr(dir,sigma)%X(dl_Vhar)(:) = M_ZERO
      call X(lr_orth_response)(m, st, lr(dir,sigma))
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
    R_TYPE  :: a


    do sigma=1,nsigma 
      lr(dir,sigma)%X(dl_rho)=M_ZERO
    end do

    do ik = 1, st%d%nspin
      do ist = 1, st%nst
        do i=1,m%np
          
          a=st%d%kweights(ik)*st%occ(ist, ik)*(&
               R_CONJ(st%X(psi)(i, 1, ist, ik))*lr(dir,1)%X(dl_psi)(i,1,ist,ik) +&
               R_CONJ(lr(dir,2)%X(dl_psi)(i,1,ist,ik))*st%X(psi)(i, 1, ist, ik))

          lr(dir,1)%X(dl_rho)(i,ik)=lr(dir,1)%X(dl_rho)(i,ik)+ a

          lr(dir,2)%X(dl_rho)(i,ik)=lr(dir,2)%X(dl_rho)(i,ik)+R_CONJ(a)

        end do
      end do
    end do
    
  end  subroutine build_rho_dynamic

end subroutine X(get_response_e)


! ---------------------------------------------------------
subroutine X(static_response) (sys, h, lr, props, pol, hpol)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:) ! lr(NDIM,1,1)
  type(pol_props_t),     intent(in)    :: props
  FLOAT,               intent(out)   :: pol(:,:)
  FLOAT,               intent(out)   :: hpol(:,:,:)

  integer :: i, j, k
  integer :: ispin, ist, ispin2, ist2, n, np, dim, ik

  R_TYPE :: prod

  R_TYPE, allocatable :: tmp(:,:), dVde(:,:,:), drhode(:,:)
  FLOAT,  allocatable :: kxc(:,:,:,:)

  np  = sys%gr%m%np
  dim = sys%gr%sb%dim

  call push_sub('em_resp.static_response')

  message(1) = "Info: Calculating static properties"
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

    call mix_init(lr(i, 1)%mixer, sys%gr%m, 1, sys%st%d%nspin)

    call X(get_response_e)(sys, h, lr(:,:), i, 1, R_TOTYPE(M_ZERO), props)
    call mix_end(lr(i, 1)%mixer)

    do ispin = 1, sys%st%d%nspin
      ! the potential derivatives

      ! Hartree and the potential and the derivative of the
      !  potential associated at the electric field

      dVde(1:np,ispin, i) = sys%gr%m%x(1:np, i)

      if(props%add_hartree) dVde(1:np, ispin, i) = dVde(1:np, ispin, i) + lr(i,1)%X(dl_Vhar)(1:np) 

      if(props%add_fxc) then 
        ! xc
        do ispin2 = 1, sys%st%d%nspin
          dVde(1:np, ispin,i) = dVde(1:np, ispin,i) + &
               lr(i, 1)%dl_Vxc(1:np, ispin, ispin2)*lr(i,1)%X(dl_rho)(1:np, ispin2)
        end do
      end if

    end do

    ! the density
    do n = 1, np
      drhode(n, i) = sum(lr(i, 1)%X(dl_rho)(n, 1:sys%st%d%nspin))
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
  end if

  hpol = M_ZERO

  do i = 1, dim
    do j = 1, dim
      do k = 1, dim

        do ik=1, sys%st%d%nik
          do ispin = 1, sys%st%d%nspin
            do ist = 1, sys%st%nst

              if( sys%st%occ(ist, ik) > lr_min_occ ) then 
                ! <D\psi_n | P_c DV_scf P_c | D\psi_n >
                
                tmp(1:np, 1)  = dVde(1:np, ispin, j) * lr(k,1)%X(dl_psi)(1:np, 1, ist, ispin)
                hpol(i, j, k) = hpol(i, j, k) + &
                     sys%st%d%kweights(ik)*sys%st%occ(ist, ik)* sum(R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * &
                     tmp(1:np, 1) * sys%gr%m%vol_pp(1:np))
                
                do ispin2 = 1, sys%st%d%nspin
                  do ist2 = 1, sys%st%nst
                    if( sys%st%occ(ist2, ik) > lr_min_occ ) then 
                      prod = sum(R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * &
                           lr(k, 1)%X(dl_psi)(1:np, 1, ist2, ispin2) * sys%gr%m%vol_pp(1:np))
                      
                      prod = prod * sum( &
                           R_CONJ(sys%st%X(psi)(1:np, 1, ist2, ispin2)) * &
                           dVde(1:np, ispin, j) * sys%st%X(psi)(1:np, 1, ist, ispin) * &
                         sys%gr%m%vol_pp(1:np))
                      
                      hpol(i, j, k) = hpol(i, j, k) - sys%st%d%kweights(ik)*sys%st%occ(ist, ik)*prod
                    end if
                  end do ! ist2
                end do ! ispin2

              end if
            
            end do ! ist
          end do ! ispin
        end do !ik

        if(props%add_fxc) then 
          hpol(i, j, k) = hpol(i, j, k) + &
               sum(kxc(1:np, 1, 1, 1) * drhode(1:np, i) * drhode(1:np, j)*drhode(1:np, k) * &
               sys%gr%m%vol_pp(1:np))/CNST(6.0)
        end if
        
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
end subroutine X(static_response)
