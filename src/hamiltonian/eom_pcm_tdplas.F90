!! Copyright (C) 2017 Gabriel Gil, Stefano Corni, Silvio Pipolo, Carlo Andrea Rozzi
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

module eom_pcm_tdplas_m
  use global_oct_m
  use messages_oct_m
  private
  public :: pcm_charges_propagation, pcm_tessera_t, debye_param_t, drude_param_t 
  save

  type :: pcm_tessera_t
    FLOAT :: point(1:3)  !< representative point of the tessera  
    FLOAT :: area        !< area of the tessera
    FLOAT :: normal(1:3) !< unitary outgoing vector normal to the tessera surface  
    FLOAT :: r_sphere    !< radius of the sphere to which the tessera belongscd
  end type pcm_tessera_t

  type :: debye_param_t
    FLOAT :: eps_d
    FLOAT :: tau
  end type debye_param_t

  type :: drude_param_t
    FLOAT :: aa
    FLOAT :: gm
    FLOAT :: w0
  end type drude_param_t

  type(pcm_tessera_t), allocatable :: cts_act(:)

  type(debye_param_t) :: deb
  type(drude_param_t) :: drl

  integer :: nts_act

  FLOAT :: eps_0
  character(3) :: which_eps

  FLOAT, parameter :: zero=0.d0
  FLOAT, parameter :: one=1.d0
  FLOAT, parameter :: two=2.d0
  FLOAT, parameter :: pi=3.141592653589793D+00
  FLOAT, parameter :: twp=two*pi

  FLOAT :: dt

  FLOAT, allocatable :: q_tp(:),dq_tp(:)
  FLOAT, allocatable :: pot_tp(:)
  FLOAT :: f1,f2,f3,f4,f5
  FLOAT, allocatable :: force(:),force_p(:)

  FLOAT, allocatable :: cals(:,:),cald(:,:)
  FLOAT, allocatable :: sm1(:,:)
  FLOAT, allocatable :: eigv(:),eigt(:,:)
  FLOAT, allocatable :: sm12(:,:),sp12(:,:)
  FLOAT, allocatable :: matq0(:,:),matqd(:,:)
  FLOAT, allocatable :: matqv(:,:),matqq(:,:)

  contains

  subroutine pcm_charges_propagation(q_t, pot_t, this_dt, this_cts_act, this_eps_0, this_eps, this_deb, this_drl) 
   save
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   FLOAT, intent(in) :: this_dt

   type(pcm_tessera_t), intent(in) :: this_cts_act(:)

   FLOAT, intent(in) :: this_eps_0
   character(*), intent(in) :: this_eps
   type(debye_param_t), optional, intent(in) :: this_deb ! By now, these arguments are both optional,
   type(drude_param_t), optional, intent(in) :: this_drl ! and at least one of them should be required.
							 ! A tentative solution: using interface driving
							 ! directly to each case and adding them to the module. E.g.
							 !interface pcm_charges_propagation
    							 !module procedure pcm_charges_propagation_deb, pcm_charges_propagation_drl
  							 !end interface

   logical :: firsttime=.true.


   PUSH_SUB(pcm_charges_propagation)

   if(firsttime) then
    dt=this_dt
    nts_act=size(this_cts_act)
    if (size(q_t) /= nts_act) stop
    allocate(cts_act(nts_act))
    cts_act=this_cts_act
    eps_0=this_eps_0
    which_eps=this_eps
    if (which_eps=='deb' .and. .not.(present(this_deb))) then
     stop 'Debye dielectric function requires three parameters!'
    else if (which_eps=='deb' .and. present(this_deb)) then
     deb=this_deb
    endif
    if (which_eps=='drl' .and. .not.(present(this_drl))) then
     stop 'Drude-Lorentz dielectric function requires three parameters!'
    else if (which_eps=='drl' .and. present(this_drl)) then
     drl=this_drl
    endif
    call init_BEM
    call init_charges(pot_t)
    q_t=q_tp
    firsttime=.false.
    POP_SUB(pcm_charges_propagation)
    return
   endif

   if( which_eps .eq. "deb") then
    if( deb%tau /= zero ) then
     call pcm_ief_prop_deb(pot_t,q_t)
    else
     POP_SUB(pcm_charges_propagation)
     return
    endif
   elseif( which_eps .eq. "drl" ) then
    call pcm_ief_prop_vv_ief_drl(pot_t,q_t)
   else
    stop 'Only Debye or Drude-Lorent is allowed!'
   endif
   q_tp=q_t
   pot_tp=pot_t

   POP_SUB(pcm_charges_propagation)
  end subroutine pcm_charges_propagation

  subroutine init_charges(pot_t)
   FLOAT, intent(in)  :: pot_t(:)

   PUSH_SUB(init_charges)

   allocate(q_tp(nts_act))
   allocate(pot_tp(nts_act))
   q_tp=matmul(matq0,pot_t)
   pot_tp=pot_t
   if( which_eps .eq. "drl" ) then
    allocate(dq_tp(nts_act))
    allocate(force_p(nts_act))
    allocate(force(nts_act))
    dq_tp=zero
    force_p=zero
    call init_vv_propagator
   endif

   POP_SUB(init_charges)
  end subroutine

  subroutine init_vv_propagator
   PUSH_SUB(init_vv_propagator)

   f1=dt*(1.d0-dt*0.5d0*drl%gm)
   f2=dt*dt*0.5d0
   f3=1.d0-dt*drl%gm*(1.d0-dt*0.5*drl%gm)
   f4=0.5d0*dt
   f5=drl%gm*f2

   POP_SUB(init_vv_propagator)
  end subroutine init_vv_propagator

  subroutine pcm_ief_prop_deb(pot_t,q_t)
   FLOAT, intent(in)  :: pot_t(:)
   FLOAT, intent(out) :: q_t(:)

   PUSH_SUB(pcm_ief_prop_deb)

   q_t(:) = q_tp(:) - dt*matmul(matqq,q_tp)   + &
		      dt*matmul(matqv,pot_tp) + &
			 matmul(matqd,pot_t-pot_tp)

   POP_SUB(pcm_ief_prop_deb)
  end subroutine pcm_ief_prop_deb

  subroutine pcm_ief_prop_vv_ief_drl(pot_t,q_t)
   FLOAT, intent(in)  :: pot_t(:)
   FLOAT, intent(out) :: q_t(:)

   PUSH_SUB(pcm_ief_prop_vv_ief_drl)

   q_t = q_tp + f1*dq_tp + f2*force_p
   force = -matmul(matqq,q_t) + matmul(matqv,pot_t)
   dq_tp = f3*dq_tp + f4*(force+force_p) -f5*force_p
   force_p = force
   q_tp = q_t

   POP_SUB(pcm_ief_prop_vv_ief_drl)
  end subroutine pcm_ief_prop_vv_ief_drl

  subroutine init_BEM
   PUSH_SUB(init_BEM)

   allocate(matqv(nts_act,nts_act))
   allocate(matqq(nts_act,nts_act))
   allocate(matq0(nts_act,nts_act))
   allocate(matqd(nts_act,nts_act))
   call do_PCM_propMat

   POP_SUB(init_BEM)
  end subroutine init_BEM

  subroutine do_PCM_propMat
   integer :: i,j
   FLOAT, allocatable :: scr4(:,:),scr1(:,:)
   FLOAT, allocatable :: scr2(:,:),scr3(:,:)
   FLOAT, allocatable :: fact1(:),fact2(:)
   FLOAT, allocatable :: Kdiag0(:),Kdiagd(:)
   FLOAT :: sgn,fac_eps0,fac_epsd
   FLOAT :: temp

   PUSH_SUB(do_PCM_propMat)

   sgn=one

   allocate(cals(nts_act,nts_act),cald(nts_act,nts_act))
   allocate(Kdiag0(nts_act),Kdiagd(nts_act))
   allocate(fact1(nts_act),fact2(nts_act))
   do i=1,nts_act
    do j=1,nts_act
     call green_s(i,j,temp)
     cals(i,j)=temp
     call green_d(i,j,temp)
     cald(i,j)=temp
    enddo
   enddo
   call allocate_TS_matrix
   call do_TS_matrix
   !deallocate(cals,cald)
   allocate(scr4(nts_act,nts_act),scr1(nts_act,nts_act))
   allocate(scr2(nts_act,nts_act),scr3(nts_act,nts_act))
   if (which_eps .eq. "deb") then
    fac_eps0=(eps_0+1)/(eps_0-1)
    Kdiag0(:)=(twp-sgn*eigv(:))/(twp*fac_eps0-sgn*eigv(:))
    fac_epsd=(deb%eps_d+1)/(deb%eps_d-1)
    Kdiagd(:)=(twp-sgn*eigv(:))/(twp*fac_epsd-sgn*eigv(:))
    fact1(:)=((twp-sgn*eigv(:))*eps_0+twp+eigv(:))/ &
    ((twp-sgn*eigv(:))*deb%eps_d+twp+eigv(:))/deb%tau
    fact2(:)=Kdiag0(:)*fact1(:)
   elseif (which_eps .eq. "drl") then
    Kdiagd(:)=0.d0
    fact2(:)=(twp-sgn*eigv(:))*drl%aa/(two*twp)
    do i=1,nts_act
     if(fact2(i).lt.0.d0) fact2(i)=0d0
    enddo
    if (drl%w0.eq.0.d0) drl%w0=1.d-8
    fact1(:)=fact2(:)+drl%w0*drl%w0
    Kdiag0(:)=fact2(:)/fact1(:)
   endif
   scr3=matmul(sm12,eigt)
   scr2=matmul(transpose(eigt),sp12)
   scr4=transpose(scr3)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact1(i)
   enddo
   matqq=matmul(scr1,scr2)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact2(i)
   enddo
   matqv=-matmul(scr1,scr4)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiag0(i)
   enddo
   matq0=-matmul(scr1,scr4)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiagd(i)
   enddo
   matqd=-matmul(scr1,scr4)
   deallocate(scr4,scr1)
   deallocate(scr2,scr3)
   deallocate(fact1,fact2)
   deallocate(Kdiag0,Kdiagd)
   call deallocate_TS_matrix

   POP_SUB(do_PCM_propMat)
  end subroutine do_PCM_propMat

  subroutine green_d(i,j,value)
   integer, intent(in)  :: i,j
   FLOAT   , intent(out) :: value
   FLOAT :: dist,diff(3)

   PUSH_SUB(green_d)

   if (i.ne.j) then
    diff=cts_act(i)%point-cts_act(j)%point
    dist=sqrt(dot_product(diff,diff))
    value=dot_product(cts_act(j)%normal,diff)/dist**3
   else
    value=-1.0694*sqrt(two*twp*cts_act(i)%area)
    value=value/(2.d0*cts_act(i)%r_sphere)/cts_act(i)%area
   endif

   POP_SUB(green_d)
  end subroutine green_d

  subroutine green_s(i,j,value)
   integer, intent(in):: i,j
   FLOAT, intent(out) :: value
   FLOAT:: dist,diff(3)

   PUSH_SUB(green_s)

   if (i.ne.j) then
    diff=cts_act(i)%point-cts_act(j)%point
    dist=sqrt(dot_product(diff,diff))
    value=one/dist
   else
    value=1.0694*sqrt(two*twp/cts_act(i)%area)
   endif

   POP_SUB(green_s)
  end subroutine green_s

  subroutine allocate_TS_matrix
   PUSH_SUB(allocate_TS_matrix)

   allocate(eigv(nts_act))
   allocate(eigt(nts_act,nts_act))
   allocate(sm12(nts_act,nts_act))
   allocate(sp12(nts_act,nts_act))

   POP_SUB(allocate_TS_matrix)
  end subroutine allocate_TS_matrix

  subroutine deallocate_TS_matrix
   PUSH_SUB(deallocate_TS_matrix)

   deallocate(eigv)
   deallocate(eigt)
   deallocate(sm12)
   deallocate(sp12)

   POP_SUB(deallocate_TS_matrix)
  end subroutine deallocate_TS_matrix

  subroutine do_TS_matrix
   integer :: i,j
   integer :: info,lwork,liwork
   FLOAT, allocatable :: scr1(:,:),scr2(:,:),eigt_t(:,:)
   FLOAT :: sgn,fac_eps0,fac_epsd
   FLOAT:: temp,fact1,fact2
   character jobz,uplo
   integer, allocatable :: iwork(:)
   FLOAT,allocatable :: work(:)

   PUSH_SUB(do_TS_matrix)

   allocate(scr1(nts_act,nts_act),scr2(nts_act,nts_act))
   allocate(eigt_t(nts_act,nts_act))
   allocate(work(1+6*nts_act+2*nts_act*nts_act))
   allocate(iwork(3+5*nts_act))

   sgn=one

   jobz = 'V'
   uplo = 'U'
   lwork = 1+6*nts_act+2*nts_act*nts_act
   liwork = 3+5*nts_act
   eigt = cals
   call dsyevd(jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
               iwork,liwork,info)
   do i=1,nts_act
    if(eigv(i).lt.0.d0) then
     write(6,*) "WARNING:",i," eig of S is negative!"
     write(6,*) "   I put it to 1e-8"
     eigv(i)=1.d-8
    endif
    scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
   enddo
   eigt_t=transpose(eigt)
   sp12=matmul(scr1,eigt_t)
   do i=1,nts_act
    scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
   enddo
   sm12=matmul(scr1,eigt_t)
   do i=1,nts_act
    scr1(:,i)=cald(:,i)*cts_act(i)%area
   enddo
   scr2=matmul(matmul(sm12,scr1),sp12)
   do j=1,nts_act
    do i=1,nts_act
     eigt(i,j)=0.5*(scr2(i,j)+scr2(j,i))
    enddo
   enddo
   call dsyevd(jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
               iwork,liwork,info)

   if (which_eps .eq. 'drl') then
    scr1=matmul(sm12,eigt)
    scr2=transpose(scr1)
    do i=1,nts_act
     fact2=(twp-sgn*eigv(i))*drl%aa/(two*twp)
     fact1=fact2+drl%w0*drl%w0
     if (fact1.lt.0.d0) eigv(i)=-twp
     scr1(:,i)=scr1(:,i)*fact2/fact1
    enddo
    eigt_t=-matmul(scr1,scr2)
   endif
   if(allocated(cals).and.allocated(cald)) deallocate(cals,cald)
   deallocate(work)
   deallocate(scr1,scr2,eigt_t)

   POP_SUB(do_TS_matrix)
  end subroutine do_TS_matrix

end module eom_pcm_tdplas_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

