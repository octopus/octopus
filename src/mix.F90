#include "config.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
! Module mixing 0.2 (25/06/2001)                                              !
!                                                                             !
! Performs linear, Broyden or (not yet) Anderson-Pulay mixing  of the         !
!       densities.                                                            !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mix
use global
use fdf
use mesh
use states

implicit none

private
public :: mix_type, mix_init, mix_end, mix_dens, calcdens

integer, parameter :: LINEAR   = 0, &
                      ANDERSON = 1, &
                      BROYDEN  = 2

! Anderson Mixing
integer, parameter :: A_IN = 1, A_OUT = 2

type mix_type
  integer  :: type_of_mixing                
  integer  :: norbitals               !  Number of orbitals
  real(r8) :: weight                  !  Weight of each component to perform 
                                      !  the integrals.

  real(r8) :: alpha                   !  rho_new = (1-a)*rho_in + a*rho_out

  ! Anderson Mixing
  !real(r8) :: anderson_mixing_alpha
  !real(r8), pointer :: densities(:,:,:,:)
  !integer :: anderson_number

  ! Broyden Mixing
  integer :: broyden_number
  real(r8), pointer :: df(:,:,:), dv(:,:,:)
end type mix_type

contains

! Initialization...
subroutine mix_init(smix, m, st)
  type(mix_type), intent(out) :: smix
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st

  ! check input parameters
  smix%type_of_mixing = fdf_integer("TypeOfMixing", 0)  
  if(smix%type_of_mixing == 1) then
    message(1) = 'Anderson mixing (TypeOfMixing = 1) currently not implemented)'
    call write_fatal(1)
  end if
  if(smix%type_of_mixing < 0 .or. smix%type_of_mixing > 2) then
    message(1) = 'Type of mixing passed to mix_init not allowed'
    call write_fatal(1)
  end if

  smix%alpha = fdf_double("Mixing", 0.3_r8);
  if(smix%alpha <= 0.0_r8 .or. smix%alpha > 1.0_r8) then
    write(message(1), '(a, f14.6,a)') "Input: '", smix%alpha, &
         "' is not a valid Mixing"
    message(2) = '(0 < Mixing <= 1)'
    call write_fatal(2)
  end if

  if(smix%type_of_mixing == 2) then
    smix%broyden_number = fdf_integer("BroydenNumber",3)
    if(smix%broyden_number <= 1 .or. smix%broyden_number > 5) then
      write(message(1), '(a, i4,a)') "Input: '", smix%broyden_number, &
           "' is not a valid BroydenNumber"
      message(2) = '(1 < BroydenNumber <= 5)'
      call write_fatal(2)
    end if
  end if

  if(smix%type_of_mixing == BROYDEN) then
    allocate(smix%df(m%np, smix%broyden_number, st%ispin))
    allocate(smix%dv(m%np, smix%broyden_number, st%ispin))
  end if

  return
end subroutine mix_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mix_end(smix)
  type(mix_type), intent(inout) :: smix

  if(smix%type_of_mixing == BROYDEN) then
    if(associated(smix%df)) then
      deallocate(smix%df); nullify(smix%df);
      deallocate(smix%dv); nullify(smix%dv);
    end if
  end if
end subroutine mix_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mix_dens(smix, iter, st, m, dist)
  type(mix_type), intent(inout) :: smix
  integer, intent(in)      :: iter
  type(states_type), intent(inout) :: st
  type(mesh_type), intent(IN) :: m
  real(r8), intent(out)    :: dist
  
  real(r8), allocatable :: rhoout(:,:), dummy(:)
  integer :: is, errorflag
    
  allocate(rhoout(m%np, st%ispin))
  call calcdens(st, m%np, rhoout)

  dist = 0._r8
  allocate(dummy(m%np))
  do is = 1, st%ispin
    dummy = (st%rho(:,is) - rhoout(:,is))**2
    dist = dist + dmesh_integrate(m, dummy)
  end do
  dist = sqrt(dist)
  deallocate(dummy)
    
  select case(smix%type_of_mixing)
  case(LINEAR)
    call mix_linear(smix, m%np, st%ispin, st%rho, rhoout)
!!$    case(ANDERSON)
!!$      call anderson_mix(rho,rhoout, iter, errorflag)
  case(BROYDEN)
    do is = 1, st%ispin
!      call mix_broyden(is, smix, rho(:, is),rhoout(:, is), iter, errorflag)
    end do
  end select
    
  deallocate(rhoout)
    
end subroutine mix_dens

! Calculates the new density out the wavefunctions and occupations...
subroutine calcdens(st, np, rho)
  type(states_type), intent(IN) :: st
  integer, intent(in) :: np
  real(r8), intent(out) :: rho(np, st%ispin)

  integer :: i, ik, p, sp

  if(st%ispin == 2) then
    sp = 2
  else
    sp = 1
  end if

  ! TODO: for polymers, do not forget to introduce integration factor
  ! in momentum space
  rho = 0._r8
  do ik = 1, st%nik, sp
    do p  = 1, st%nst
      do i = 1, np
        ! spin-up density
        rho(i, 1) = rho(i, 1) + st%occ(p, ik)*R_ABS(st%R_FUNC(psi)(i, 1, p, ik))**2

        ! spin-down density
        if(st%ispin == 2) then
          rho(i, 2) = rho(i, 2) + st%occ(p, ik+1)*R_ABS(st%R_FUNC(psi)(i, 1, p, ik+1))**2
        end if

        ! off-diagonal densities
        if(st%ispin == 4) then
          rho(i, 2) = rho(i, 2) + st%occ(p, ik)*R_ABS(st%R_FUNC(psi)(i, 2, p, ik))**2
!          rho(i, 3) = st%occ(p, ik)*&
!               st%R_FUNC(psi)(i, 1, p, ik)*R_CONJ(st%R_FUNC(psi)(i, 2, p, ik))
!          rho(i, 4) = R_CONJ(rho(i, 3)) ! this is in principle not necessary!
        end if

      end do
    end do
  end do

  return
end subroutine calcdens

! Performs the linear mixing...
subroutine mix_linear(smix, np, ispin, rho, rhoout)
  type(mix_type), intent(IN) :: smix
  integer, intent(in) :: np, ispin
  real(r8), intent(inout) :: rho(np, ispin)
  real(r8), intent(IN) :: rhoout(np, ispin)

  call dscal(np*ispin, 1.0_r8-smix%alpha, rho, 1)
  call daxpy(np*ispin, smix%alpha, rhoout, 1, rho, 1)

  return
end subroutine mix_linear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! Performs the Anderson-Pulay mixing...                                       !
!!$!                                                                             !
!!$!       (Ths subroutine is just an interface which takes care of the saving   !
!!$!       
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$  subroutine anderson_mix(rho,rhoout,iter,errorflag)
!!$
!!$  implicit none
!!$
!!$  integer, intent(in) :: iter
!!$  real(r8), intent(inout) :: rho(dim,nspin)
!!$  real(r8), intent(in) :: rhoout(dim,nspin)
!!$  integer, intent(out) :: errorflag
!!$
!!$  integer :: niter,j,is
!!$
!!$  if(iter<1) then
!!$    errorflag=-1
!!$    return
!!$  endif
!!$
!!$  if(iter==1) then
!!$    densities(1:dim,1:nspin,0,A_IN)=rho(1:dim,1:nspin)
!!$    densities(1:dim,1:nspin,0,A_OUT)=rhoout(1:dim,1:nspin)
!!$    errorflag=-2
!!$    return
!!$  endif
!!$
!!$  niter=min(iter,anderson_number)
!!$
!!$  do j=-(anderson_number-1),-(anderson_number-1)+niter-1
!!$     densities(1:dim,1:nspin,j,:)=densities(1:dim,1:nspin,j+1,:)
!!$  enddo
!!$  densities(1:dim,1:nspin,0,A_IN)=rho(1:dim,1:nspin)
!!$  densities(1:dim,1:nspin,0,A_OUT)=rhoout(1:dim,1:nspin)
!!$
!!$  do is=1,nspin
!!$     call do_anderson_mix(densities(1:dim,is,-(niter-1):0,A_IN),    &
!!$                          densities(1:dim,is,-(niter-1):0,A_OUT),   &
!!$                          niter,rho(1:dim,is),errorflag)
!!$  enddo
!!$
!!$  return
!!$  end subroutine anderson_mix
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$! It receives as input the niter input densities nin, the niter output        !
!!$!       densities nout. It outputs the new density nnew. The dot product is   !
!!$!       defined through the weight w.                                         !
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$  subroutine do_anderson_mix(nin,nout,niter,nnew,errorflag)
!!$
!!$  implicit none
!!$
!!$  integer, intent(in) :: niter
!!$  real(r8), intent(in) :: nin(dim,-(niter-1):0), nout(dim,-(niter-1):0)
!!$  real(r8), intent(out) :: nnew(dim)
!!$  integer, intent(out) :: errorflag
!!$
!!$  real(r8) :: kappa(-(niter-1):0,-(niter-1):0), betha(-(niter-1):0 )
!!$  integer :: i,j,l
!!$
!!$
!!$! Calculates the kappa matrix
!!$
!!$  kappa = 0.0_r8
!!$  do i=-(niter-1),0
!!$  do j=i,0
!!$     do l=1,dim
!!$        kappa(i,j)=kappa(i,j) + (nout(l,i)-nin(l,i))*(nout(l,j)-nin(l,j))
!!$     enddo
!!$     kappa(i,j) = kappa(i,j)*weight
!!$     !kappa(i,j)=weight*sum( (nout(1:dim,i)-nin(1:dim,i)) * (nout(1:dim,j)-nin(1:dim,j)) )
!!$     kappa(j,i)=kappa(i,j)
!!$  enddo
!!$  enddo
!!$  
!!$! Solves the system for the bethas through a call to solve_system...  
!!$
!!$  errorflag = 0
!!$  call solve_system(dim,niter,kappa,betha,errorflag)
!!$
!!$
!!$! Finds out nnew
!!$
!!$
!!$  nnew = 0.0_r8
!!$  do j = - (niter-1), 0
!!$     nnew(:) = (1.0_r8-alpha)*betha(j)*nin(:,j) + &
!!$               alpha*betha(j)*nout(:,j)
!!$  enddo
!!$
!!$! A little test to check  out normalization...
!!$  write(*,*) 'CHECK OF NORMALIZATION IN anderson...',sum(nnew(:))*weight
!!$
!!$  return
!!$  end subroutine do_anderson_mix
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$
!!$  subroutine solve_system(n,niter,kappa,betha,errorflag)
!!$
!!$  implicit none
!!$
!!$  integer, intent(in) :: n,niter
!!$  real(r8), intent(in) :: kappa(n,1:niter)
!!$  real(r8), intent(out) :: betha(1:niter)
!!$  integer, intent(out) :: errorflag
!!$!! This is for the check:
!!$  integer :: i
!!$  real(r8) :: eps=1.0e-8_r8,val
!!$!!
!!$
!!$  integer :: ipiv,info
!!$  real(r8) :: m(niter+1,niter+1),vec(niter+1),res(niter+1)
!!$
!!$  m(1:niter,1:niter) = kappa(1:niter,1:niter)
!!$  m(1:niter,niter+1) = -0.5_r8
!!$  m(niter+1,1:niter) = 1.0_r8
!!$  m(niter+1,niter+1) = 0.0_r8
!!$  vec(1:niter) = 0.0_r8
!!$  vec(niter+1) = 1.0_r8
!!$
!!$  errorflag = 0
!!$  call dgesv(niter,1,m,niter,ipiv,vec,niter,info)
!!$  if(info.ne.0) then
!!$     errorflag = -99
!!$     return
!!$  endif
!!$
!!$  betha(1:niter)=res(1:niter)
!!$  ! A little check to see if the system has been successfully solved...
!!$  do i=-(niter-1),0
!!$     val = sum(kappa(i,:)*betha(:))-0.5_r8*res(niter+1)
!!$     if(abs(val)>eps) stop 'System not solved'
!!$  enddo
!!$  val = sum(betha)
!!$  if(abs(val-1.0_r8)>eps) stop 'System not solved'
!!$
!!$  end subroutine solve_system
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Broyden mixing...                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine mix_broyden(is, smix, vin, vout, iter, errorflag)
!!$    integer, intent(in) :: is
!!$    type(mix_type), intent(inout) :: smix
!!$    integer, intent(in)     :: iter              ! counter of the iterations
!!$    integer, intent(out)    :: errorflag         ! not really used yet.
!!$    real(r8), intent(inout) :: vout(smix%dim), & ! the new density; 
!!$                               vin(smix%dim)     ! the old H+xc density
!!$
!!$    integer  ::        &
!!$         ndim,         & ! input:  leading dimens. of arrays vout, vin, df, dv
!!$         n_iter          ! input:  numb. of iterations used in potential mixing
!!$
!!$    real(r8) ::        & 
!!$         alphamix,     & ! input:  mixing factor
!!$         dr2             ! output: parameter for the control convergence 
!!$
!!$    integer, parameter :: maxter = 5 
!!$
!!$    integer  ::        &
!!$         iunit,        & ! counter on I/O unit numbers  
!!$         n,            & ! do-loop counter on array-dimension ndim
!!$         i,            & ! first do-loop counter on number of (part.) iterat. 
!!$         j,            & ! second do-loop counter on number of (part.) iterat. 
!!$         iwork(maxter),& ! dummy array used as output by libr. routines 
!!$         info,         & ! flag saying if the exec. of libr. routines was ok
!!$         iter_used,    & ! actual number of iteractions used
!!$         ipos,         & ! index of the present iteraction
!!$         inext           ! index of the next iteraction
!!$
!!$    integer :: ndimtot
!!$    real(r8) ::                    &
!!$         vinsave(smix%dim),        & ! work space
!!$         beta(maxter, maxter),     & !
!!$         gamma,                    & !
!!$         work(maxter),             & !
!!$         w(maxter),                & !
!!$         w0,                       & ! a constant value 
!!$         norm                        ! Euclidean norm of ndim vector df(i,j) 
!!$                                     ! (i=1,ndim ; j fixed)
!!$
!!$    real(r8), external ::  & 
!!$         ddot,             & ! function which computes the dot product of two vectors
!!$         dnrm2               ! function which computes the Euclidean norm of a vector
!!$
!!$
!!$    ndim = smix%dim
!!$    alphamix = smix%alpha
!!$    n_iter = smix%broyden_number
!!$
!!$    w0 = 0.01_r8
!!$    w  = maxter*1.0_r8
!!$
!!$    if (iter.lt.1) then
!!$      errorflag = -1
!!$      return
!!$    endif
!!$    if (n_iter.gt.maxter) then
!!$      errorflag = -2
!!$      return
!!$    endif
!!$    if (ndim.le.0) then    
!!$      errorflag = -3
!!$      return
!!$    endif
!!$
!!$    call daxpy( ndim, -1.0_r8, vin, 1, vout, 1)
!!$
!!$    ndimtot = ndim
!!$
!!$    iter_used = min(iter - 1, n_iter)
!!$
!!$    ! ipos is the position in which results from the present iteraction 
!!$    ! are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
!!$
!!$    ipos = iter - 1 - ((iter - 2)/n_iter)*n_iter
!!$
!!$    if(iter.gt.1) then
!!$      do n = 1, ndim
!!$        smix%df(n, ipos, is) = vout(n) - smix%df(n, ipos, is)
!!$        smix%dv(n, ipos, is) = vin(n)  - smix%dv(n, ipos, is)
!!$      end do
!!$      norm = DNRM2(ndim, smix%df(1, ipos, is), 1)
!!$      call dscal (ndim, 1.d0/norm, smix%df(1, ipos, is), 1)
!!$      call dscal (ndim, 1.d0/norm, smix%dv(1, ipos, is), 1)
!!$    endif
!!$
!!$    call dcopy(ndim, vin ,1,vinsave,1)
!!$    
!!$    do i = 1,iter_used
!!$      do j = i + 1, iter_used
!!$        beta(i, j) = w(i)*w(j)*DDOT(ndim, smix%df(1, j, is), 1, smix%df(1, i, is), 1)
!!$      end do
!!$      beta(i, i) = w0**2 + w(i)**2
!!$    end do
!!$    
!!$    call dsytrf('u', iter_used, beta, maxter, iwork, work, maxter, info)
!!$    if(info .ne. 0) then
!!$      errorflag = 1
!!$      return
!!$    endif
!!$    call dsytri('u', iter_used, beta, maxter, iwork, work, info)
!!$    if(info .ne. 0) then
!!$      errorflag = 2
!!$      return
!!$    endif
!!$    
!!$    do i = 1, iter_used
!!$      do j = i + 1, iter_used
!!$        beta(j, i) = beta(i, j)
!!$      end do
!!$    end do
!!$    
!!$    do i = 1, iter_used
!!$      work(i) = ddot(ndim, smix%df(1, i, is), 1, vout, 1)
!!$    end do
!!$    
!!$    call daxpy(ndim, alphamix, vout, 1, vin, 1)
!!$    
!!$    do i = 1, iter_used
!!$      gamma = 0.d0
!!$      do j = 1, iter_used
!!$        gamma = gamma + beta(j, i)*w(j)*work(j)
!!$      end do
!!$      do n = 1, ndim
!!$        vin(n) = vin(n) - w(i)*gamma*(alphamix*smix%df(n, i, is) + smix%dv(n, i, is))
!!$      end do
!!$    end do
!!$    
!!$    inext = iter - ((iter - 1)/n_iter)*n_iter
!!$    call dcopy(ndim, vout, 1, smix%df(1, inext, is), 1)
!!$    call dcopy(ndim, vinsave, 1, smix%dv(1, inext, is), 1)
!!$    
!!$    return
!!$  end subroutine mix_broyden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  end module mix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "undef.F90"
