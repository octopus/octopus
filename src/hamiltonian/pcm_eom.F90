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

module pcm_eom_oct_m
  use global_oct_m
  use messages_oct_m
  private
  public :: pcm_charges_propagation, pcm_tessera_t, debye_param_t, drude_param_t 
  save

  !> tesselation derived type
  type :: pcm_tessera_t
    FLOAT :: point(1:3)  !< representative point of the tessera  
    FLOAT :: area        !< area of the tessera
    FLOAT :: normal(1:3) !< unitary outgoing vector normal to the tessera surface  
    FLOAT :: r_sphere    !< radius of the sphere to which the tessera belongscd
  end type pcm_tessera_t

  type(pcm_tessera_t), allocatable :: cts_act(:) !< tesselation arrays (nts_act)
  integer :: nts_act				 !< number of tesserae

  !> set of parameters for Debye dielectric model
  type :: debye_param_t
    FLOAT :: eps_0 !< static  dielectric constant (at ZERO     frequency of the field)
    FLOAT :: eps_d !< dynamic dielectric constant (at INFINITE frequency of the field)
    FLOAT :: tau   !< Debye relaxation time
  end type debye_param_t

  !> set of parameters for Drude-Lorentz dielectric model
  !> no parameter allows the match with the dynamic dielectric constant
  type :: drude_param_t
    FLOAT :: aa !< parameter matching the static dielectric constant
    FLOAT :: gm !< damping of the plasma oscillation
    FLOAT :: w0 !< plasma frequency
  end type drude_param_t

  type(debye_param_t) :: deb
  type(drude_param_t) :: drl

  character(3) :: which_eps			 !< character flag for Debye ('deb') and Drude-Lorentz models ('drl')

  FLOAT :: dt 					 !< time-step of the propagation

  FLOAT, allocatable :: q_tp(:), dq_tp(:)	 !< polarization charges and variation of in it in the previous iteration
  FLOAT, allocatable :: pot_tp(:)		 !< potential in previous iteration
                                                 !< See Chem.Phys.Lett. 429 (2006) 310-316 for Velocity-Verlet (VV) algorithm... 
  FLOAT :: f1, f2, f3, f4, f5			 !< auxiliar constants for VV
  FLOAT, allocatable :: force(:),force_p(:)	 !< analogous to force in the equation of motion for the pol.charges

  						 !< In J.Phys.Chem.A 2015, 119, 5405-5416...
  FLOAT, allocatable :: cals(:,:),cald(:,:)	 !< Calderon matrices S and D from Eq.(5), ibid.
  FLOAT, allocatable :: eigv(:),eigt(:,:)        !< \Lambda and T matrices from Eq.(10), ibid.
  FLOAT, allocatable :: sm12(:,:),sp12(:,:)      !< S^{-1/2} and S^{1/2}
  FLOAT, allocatable :: matq0(:,:),matqd(:,:)    !< Q^{IEF(d)}_0 (not used in ref.) and Q^{IEF(d)}_d from Eq.(18) with eps_0/eps_d
  FLOAT, allocatable :: matqv(:,:),matqq(:,:)    !< \tilde{Q} and R matrices from Eq.(38)-(39), respectively
                                                 !< Q^{IEF(d)}_d, \tilde{Q} and R matrices are those that enter the EOM eq.(37)
  !> mathematical constants
  FLOAT, parameter :: twopi=M_TWO*M_Pi
  FLOAT, parameter :: fourpi=M_FOUR*M_Pi

  contains

  !------------------------------------------------------------------------------------------------------------------------------
  !> Driving subroutine for the Equation of Motion (EOM) propagation of the polarization charges
  !> within the Integral Equation Formalism (IEF) formulation of the Polarization Continuum Model (PCM).
  subroutine pcm_charges_propagation(q_t, pot_t, this_dt, this_cts_act, this_eps, this_deb, this_drl) 
   save
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   FLOAT, intent(in) :: this_dt

   type(pcm_tessera_t), intent(in) :: this_cts_act(:)

   character(*), intent(in) :: this_eps !< type of dielectric model to be used, either Debye ('deb') or Drude-Lorentz ('drl')
   type(debye_param_t), optional, intent(in) :: this_deb
   type(drude_param_t), optional, intent(in) :: this_drl 

   logical :: firsttime=.true.


   PUSH_SUB(pcm_charges_propagation)

   if(firsttime) then
    dt=this_dt
    nts_act=size(this_cts_act)
    if (size(q_t) /= nts_act) then
      message(1) = "pcm_charges_propagation: Number of tesserae do not coincide with size of PCM charges array."
      call messages_fatal(1)     
    endif 
    allocate(cts_act(nts_act))
    cts_act=this_cts_act
    which_eps=this_eps
    if (which_eps=='deb' .and. .not.(present(this_deb))) then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Debye dielectric function requires three parameters."
     call messages_fatal(1)
    else if (which_eps=='deb' .and. present(this_deb)) then
     deb=this_deb
    endif
    if (which_eps=='drl' .and. .not.(present(this_drl))) then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Drude-Lorentz dielectric function requires three parameters."
     call messages_fatal(1)
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
    if( deb%tau /= M_ZERO ) then
     call pcm_ief_prop_deb(pot_t,q_t)
    else
     POP_SUB(pcm_charges_propagation)
     return
    endif
   elseif( which_eps .eq. "drl" ) then
    call pcm_ief_prop_vv_ief_drl(pot_t,q_t)
   else
    message(1) = "pcm_charges_propagation: EOM-PCM error. Only Debye or Drude-Lorent dielectric models are allowed."
    call messages_fatal(1)
   endif
   q_tp   = q_t
   pot_tp = pot_t

   POP_SUB(pcm_charges_propagation)
  end subroutine pcm_charges_propagation

  !------------------------------------------------------------------------------------------------------------------------------
  !> Polarization charges initialization in equilibrium with the initial potential.
  subroutine init_charges(pot_t)
   FLOAT, intent(in)  :: pot_t(:)

   PUSH_SUB(init_charges)

   allocate(q_tp(nts_act))
   allocate(pot_tp(nts_act))
   q_tp=matmul(matq0,pot_t) !< applying the static IEF-PCM response matrix (corresponging to epsilon_0)
   pot_tp=pot_t             !< to the initial potential, which is assumed to be the ground state potential (!)
   if( which_eps .eq. "drl" ) then
    allocate(dq_tp(nts_act))
    allocate(force_p(nts_act))
    allocate(force(nts_act))
    dq_tp=M_ZERO
    force_p=M_ZERO
    call init_vv_propagator !< initializing Velocity-Verlet algorithm for the integration of EOM-PCM for Drude-Lorentz
   endif

   POP_SUB(init_charges)
  end subroutine

  !------------------------------------------------------------------------------------------------------------------------------
  !> Euler method for integrating first order EOM for the polarization charges within IEF-PCM
  !> in the case of Debye dielectric functions.
  subroutine pcm_ief_prop_deb(pot_t,q_t)
   FLOAT, intent(in)  :: pot_t(:)
   FLOAT, intent(out) :: q_t(:)

   PUSH_SUB(pcm_ief_prop_deb)

   !> Eq.(47) in S. Corni, S. Pipolo and R. Cammi, J.Phys.Chem.A 2015, 119, 5405-5416.
   q_t(:) = q_tp(:) - dt*matmul(matqq,q_tp)   + &
		      dt*matmul(matqv,pot_tp) + &
			 matmul(matqd,pot_t-pot_tp)

   POP_SUB(pcm_ief_prop_deb)
  end subroutine pcm_ief_prop_deb

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to initialize numerical constants required by the Velocity-Verlet (VV) algorithm.
  subroutine init_vv_propagator
   PUSH_SUB(init_vv_propagator)

   !> See E. Vanden-Eijnden, G. Ciccotti, Chem.Phys.Lett. 429 (2006) 310-316.
   !> Using the scheme in Eq.(21) and (17), ibid.
   !> See subroutine pcm_ief_prop_vv_ief_drl
   f1=dt*(1.d0-dt*0.5d0*drl%gm)
   f2=dt*dt*0.5d0
   f3=1.d0-dt*drl%gm*(1.d0-dt*0.5*drl%gm)
   f4=0.5d0*dt
   f5=drl%gm*f2

   POP_SUB(init_vv_propagator)
  end subroutine init_vv_propagator

  !------------------------------------------------------------------------------------------------------------------------------
  !> VV algorithm for integrating second order EOM for the polarization charges within IEF-PCM
  !> in the case of Drude-Lorentz dielectric functions.
  subroutine pcm_ief_prop_vv_ief_drl(pot_t,q_t)
   FLOAT, intent(in)  :: pot_t(:)
   FLOAT, intent(out) :: q_t(:)

   PUSH_SUB(pcm_ief_prop_vv_ief_drl)

   !> From Eq.(15) S. Pipolo and S. Corni, J.Phys.Chem.C 2016, 120, 28774-28781.
   !> Using the scheme in Eq.(21) and (17) of E. Vanden-Eijnden, G. Ciccotti, Chem.Phys.Lett. 429 (2006) 310-316.
   q_t = q_tp + f1*dq_tp + f2*force_p
   force = -matmul(matqq,q_t) + matmul(matqv,pot_t)
   dq_tp = f3*dq_tp + f4*(force+force_p) -f5*force_p
   force_p = force
   q_tp = q_t

   POP_SUB(pcm_ief_prop_vv_ief_drl)
  end subroutine pcm_ief_prop_vv_ief_drl

  !------------------------------------------------------------------------------------------------------------------------------
  !> Boundary Element Method (BEM) EOM-IEF-PCM matrices initialization.
  subroutine init_BEM
   PUSH_SUB(init_BEM)

   allocate(matqv(nts_act,nts_act))
   allocate(matqq(nts_act,nts_act))
   allocate(matq0(nts_act,nts_act))
   allocate(matqd(nts_act,nts_act))
   call do_PCM_propMat

   POP_SUB(init_BEM)
  end subroutine init_BEM

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build the required BEM matrices for the EOM-IEF-PCM for Debye and Drude-Lorentz cases.
  !> Following from
  !> Ref.1 - J.Phys.Chem.A 2015, 119, 5405-5416.    - Debye case
  !> Ref.2 - J.Phys.Chem.C 2016, 120, 28-774-28781. - Drude-Lorentz case
  !> The matrices are required for the EOMs, eq.(37) in Ref.1 and eq.(15) in Ref.2 
  subroutine do_PCM_propMat
   integer :: i,j
   FLOAT, allocatable :: scr4(:,:),scr1(:,:)
   FLOAT, allocatable :: scr2(:,:),scr3(:,:)
   FLOAT, allocatable :: fact1(:),fact2(:)
   FLOAT, allocatable :: Kdiag0(:),Kdiagd(:)
   FLOAT :: sgn,fac_eps0,fac_epsd
   FLOAT :: temp

   PUSH_SUB(do_PCM_propMat)

   sgn=M_ONE ! In the case of NP this is -one because the normal to the cavity is always pointing outward 
             ! and the field created by a positive unit charge outside the cavity is directed inward.

   allocate(cals(nts_act,nts_act),cald(nts_act,nts_act))
   allocate(Kdiag0(nts_act),Kdiagd(nts_act))
   allocate(fact1(nts_act),fact2(nts_act))
   !> generate Calderon S and D matrices
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
    fac_eps0=(deb%eps_0+M_ONE)/(deb%eps_0-M_ONE)			 
    Kdiag0(:)=(twopi-sgn*eigv(:))/(twopi*fac_eps0-sgn*eigv(:)) !< Eq.(14) with eps_0 in Ref.1
    fac_epsd=(deb%eps_d+M_ONE)/(deb%eps_d-M_ONE)
    Kdiagd(:)=(twopi-sgn*eigv(:))/(twopi*fac_epsd-sgn*eigv(:)) !< Eq.(14) with eps_d, ibid.
    fact1(:)=((twopi-sgn*eigv(:))*deb%eps_0+twopi+eigv(:))/ &  !< inverse of Eq.(32), ibid.
    ((twopi-sgn*eigv(:))*deb%eps_d+twopi+eigv(:))/deb%tau
    fact2(:)=Kdiag0(:)*fact1(:)				       !< tau^{-1}K_0 in Eq.(38), ibid.
   elseif (which_eps .eq. "drl") then
    Kdiagd(:)=M_ZERO					       !< from Eq.(10) up in Ref.2
    fact2(:)=(twopi-sgn*eigv(:))*drl%aa/fourpi		       !< Eq.(10) down
    do i=1,nts_act
     if(fact2(i).lt.M_ZERO) fact2(i)=M_ZERO
    enddo
    if (drl%w0.eq.M_ZERO) drl%w0=1.d-8
    fact1(:)=fact2(:)+drl%w0*drl%w0			       !< Eq.(19), ibid.
    Kdiag0(:)=fact2(:)/fact1(:)				       !< from Eq.(10) up, ibid.
   endif
   scr3=matmul(sm12,eigt)
   scr2=matmul(transpose(eigt),sp12)
   scr4=transpose(scr3)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact1(i)
   enddo
   matqq=matmul(scr1,scr2)				       !< Eq.(39) in Ref.1 and Eq.(16) in Ref.2
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact2(i)
   enddo
   matqv=-matmul(scr1,scr4)				       !< Eq.(38) in Ref.1 and Eq.(17) in Ref.2
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiag0(i)
   enddo
   matq0=-matmul(scr1,scr4)				       !< from Eq.(14) and (18) for eps_0 in Ref.1
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiagd(i)
   enddo
   matqd=-matmul(scr1,scr4)				       !< from Eq.(14) and (18) for eps_d, ibid.
   deallocate(scr4,scr1)
   deallocate(scr2,scr3)
   deallocate(fact1,fact2)
   deallocate(Kdiag0,Kdiagd)
   call deallocate_TS_matrix

   POP_SUB(do_PCM_propMat)
  end subroutine do_PCM_propMat

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build BEM matrix corresponding to the Calderon operator D, using the Green function of an isotropic medium (!)
  !> Only solvents can be treated with this. To be changed for surfaces, spherical nanoparticles, etc.
  subroutine green_d(i,j,value)
   integer, intent(in)  :: i,j
   FLOAT   , intent(out) :: value
   FLOAT :: dist,diff(3)

   PUSH_SUB(green_d)

   if (i.ne.j) then
    diff=cts_act(i)%point-cts_act(j)%point
    dist=sqrt(dot_product(diff,diff))
    value=dot_product(cts_act(j)%normal,diff)/dist**3	    !< Eq.(5) in Refs.1-2
   else
    value=-1.0694*sqrt(fourpi*cts_act(i)%area)
    value=value/(M_TWO*cts_act(i)%r_sphere)/cts_act(i)%area !< diagonal part is a bit cumbersome
   endif

   POP_SUB(green_d)
  end subroutine green_d

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build BEM matrix corresponding to the Calderon operator S, using the Green function of an isotropic medium (!)
  !> Only solvents can be treated with this. To be changed for surfaces, spherical nanoparticles, etc.    
  subroutine green_s(i,j,value)
   integer, intent(in):: i,j
   FLOAT, intent(out) :: value
   FLOAT:: dist,diff(3)

   PUSH_SUB(green_s)

   if (i.ne.j) then
    diff=cts_act(i)%point-cts_act(j)%point
    dist=sqrt(dot_product(diff,diff))
    value=M_ONE/dist					    !< Eq.(5) in Refs.1-2
   else
    value=1.0694*sqrt(fourpi/cts_act(i)%area)		    !< diagonal part is a bit cumbersome
   endif

   POP_SUB(green_s)
  end subroutine green_s

  !------------------------------------------------------------------------------------------------------------------------------
  subroutine allocate_TS_matrix
   PUSH_SUB(allocate_TS_matrix)

   allocate(eigv(nts_act))
   allocate(eigt(nts_act,nts_act))
   allocate(sm12(nts_act,nts_act))
   allocate(sp12(nts_act,nts_act))

   POP_SUB(allocate_TS_matrix)
  end subroutine allocate_TS_matrix

  !------------------------------------------------------------------------------------------------------------------------------
  subroutine deallocate_TS_matrix
   PUSH_SUB(deallocate_TS_matrix)

   deallocate(eigv)
   deallocate(eigt)
   deallocate(sm12)
   deallocate(sp12)

   POP_SUB(deallocate_TS_matrix)
  end subroutine deallocate_TS_matrix

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build matrices S^{1/2}, S^{-1/2}, T and \Lambda (notation of Refs.1-2)
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

   sgn=M_ONE

   jobz = 'V'
   uplo = 'U'
   lwork = 1+6*nts_act+2*nts_act*nts_act
   liwork = 3+5*nts_act
   eigt = cals
   call dsyevd(jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
               iwork,liwork,info)
   do i=1,nts_act
    if(eigv(i).lt.M_ZERO) then
     write(6,*) "WARNING:",i," eig of S is negative!"
     write(6,*) "   I put it to 1e-8"
     eigv(i)=1.d-8
    endif
    scr1(:,i)=eigt(:,i)*sqrt(eigv(i))
   enddo
   eigt_t=transpose(eigt)
   sp12=matmul(scr1,eigt_t) 		 !< building S^{1/2} to be used in R (Ref.1) and Q_{\omega} (Ref.2)
   do i=1,nts_act
    scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
   enddo
   sm12=matmul(scr1,eigt_t) 		 !< building S^{-1/2} to be used in R and \tilde{Q} (Ref.1) and Q_{\omega} and Q_f (Ref.2)
   do i=1,nts_act
    scr1(:,i)=cald(:,i)*cts_act(i)%area
   enddo
   scr2=matmul(matmul(sm12,scr1),sp12)   !< Eq.(10) in Ref.1 (paragraph after Eq.(9), Ref.2)
   do j=1,nts_act
    do i=1,nts_act
     eigt(i,j)=0.5*(scr2(i,j)+scr2(j,i)) !< re-symmetrizing for numerical reasons
    enddo
   enddo
   call dsyevd(jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
               iwork,liwork,info)	 !< obtaining \Lambda (eigv) and T (eigt), Eq.(10), ibid.

   if(allocated(cals).and.allocated(cald)) deallocate(cals,cald)
   deallocate(work)
   deallocate(scr1,scr2,eigt_t)

   POP_SUB(do_TS_matrix)
  end subroutine do_TS_matrix

end module pcm_eom_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

