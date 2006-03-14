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
subroutine X(dynamic)(sys, h, lr, pol, w)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:,:) ! dim, nsigma(=2), nfreq(=1)
  FLOAT,                  intent(out)   :: pol(:,:)
  R_TYPE,                  intent(in)   :: w 

  
  integer :: dir, j, freq, sigma

  R_TYPE :: rhov
  R_TYPE :: X(pol)(3,3)

  call push_sub('static_pol_lr.dynamic')

  freq = 1
  X(pol) = M_ZERO

  do dir = 1, sys%gr%sb%dim
    write(message(1), '(a,i1,a,f12.6)') 'Info: Calculating polarizability for direction ', &
         dir, ' and frequency ', dble(w)

    call write_info(1)

    do sigma = 1, 2
      call mix_init(lr(dir,sigma,freq)%mixer, sys%gr%m, 1, sys%st%d%nspin)
    end do
    if( dir==3 ) then
      call X(get_response_e)(sys, h, lr(:,:,freq), dir, nsigma=2 , omega=w)
    else
      print*, "WARNING: SETTING TO ZERO DIR ", dir
      lr(dir,2,freq)%X(dl_rho)=M_ZERO
    endif
    do sigma = 1, 2
      call mix_end(lr(dir,sigma,freq)%mixer)
    end do

    do j = 1, sys%gr%m%np
      rhov = sum(lr(dir,2,freq)%X(dl_rho)(j,1:sys%st%d%nspin))*sys%gr%m%vol_pp(j)
      X(pol)(dir, :) = X(pol)(dir, :) - sys%gr%m%x(j,:)*rhov
    end do

  end do

  print*, X(pol)

  call pop_sub()

end subroutine X(dynamic)



! ---------------------------------------------------------
subroutine X(get_response_e)(sys, h, lr, dir, nsigma, omega)
  type(system_t), target, intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: h
  type(lr_t),             intent(inout) :: lr(:,:) !ndim,nsigma
  integer,                   intent(in)    :: dir 
  integer,                   intent(in)    :: nsigma 
  R_TYPE,                   intent(in)    :: omega
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
  print*, "OMEGA",  omega
  call init_response_e()
  finish = .false.
  iter_loop: do iter=1, lr(dir,1)%max_iter

!    do sigma=1,nsigma
!      dl_rhoin(1,:,:,sigma) = lr(dir,sigma)%X(dl_rho)(:,:)
!    end do
    !we will use the two mixers, one for the real part and one for the complex

    dl_rhoin(1,:,:,1) = R_REAL(lr(dir,1)%X(dl_rho)(:,:))

    if(nsigma==2) dl_rhoin(1,:,:,2) = R_AIMAG(lr(dir,1)%X(dl_rho)(:,:))

    if(.not.h%ip_app) then
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
        
        dV(1:m%np,ik) = lr(dir,sigma)%X(dl_Vhar)(1:m%np) + m%x(1:m%np,dir)

        do ik2 = 1, st%d%nspin
          dV(1:m%np,ik) = dV(1:m%np,ik) +&
               lr(dir,sigma)%dl_Vxc(1:m%np, ik, ik2)*lr(dir,sigma)%X(dl_rho)(1:m%np,ik2)
        end do
        
        do ist = 1, st%nst
          if(st%occ(ist, ik) <= M_ZERO) cycle
          
          Y(1:m%np,1) = -dV(1:m%np,ik)*st%X(psi)(1:m%np, 1, ist, ik)
          
          
          !          if (nsigma == 2 ) then
          
          do ist2 = 1, st%nst
            a(ist2,ik)=sum(R_CONJ(st%X(psi)(1:m%np, 1, ist2, ik))*Y(1:m%np ,1)*&
                 sys%gr%m%vol_pp(1:m%np))
            Y(1:m%np,1)=Y(1:m%np,1)-a(ist2,ik)*st%X(psi)(1:m%np, 1, ist2, ik)
            !              print*,sum(Y(1:m%np,1)*st%X(psi)(1:m%np, 1, ist2, ik))
          end do
            
          !          else
          
          !          call X(lr_orth_vector)(m, st, Y, ik)
          
!          do ist2 = 1, st%nst
            !            print*,sum(R_CONJ(Y(1:m%np,1))*st%X(psi)(1:m%np, 1, ist2, ik))
!          end do

          
          call X(lr_solve_HXeY) (lr(dir,1), h, sys%gr, sys%st%d, ik, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), Y, &
               -sys%st%eigenval(ist, ik) + freq_sign*omega)
          
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
!              lr(dir,sigma)%X(dl_psi)(1:m%np,1, ist, ik)=lr(dir,sigma)%X(dl_psi)(1:m%np, 1, ist, ik)&
!                   +a(ist2,ik)/(sys%st%eigenval(ist2, ik)-sys%st%eigenval(ist, ik)+freq_sign*omega)*&
!                   st%X(psi)(1:m%np, 1, ist2, ik)
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
      call write_fatal(1);
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
      print*, sigma, lr(dir,sigma)%abs_dens
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

#ifdef R_TCOMPLEX
        lr(dir,1)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),dl_rhonew(1,:,:,2),PRECISION)
        lr(dir,2)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),-dl_rhonew(1,:,:,2),PRECISION)
#else
      do sigma=1,nsigma
        lr(dir,sigma)%X(dl_rho)(:,:) = dl_rhonew(1,:,:,1)
      end do
#endif
    end if

!    print*, "SUM", sum(dl_rhonew(1,:,:,2))
    
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
    integer :: msigma
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


