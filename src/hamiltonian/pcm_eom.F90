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
  use io_oct_m
  use messages_oct_m
  use profiling_oct_m

  private
  public :: pcm_charges_propagation, pcm_eom_enough_initial, pcm_eom_end, pcm_tessera_t, debye_param_t, drude_param_t 
  save

  !> tesselation derived type
  type :: pcm_tessera_t
    FLOAT :: point(1:3)                           !< representative point of the tessera  
    FLOAT :: area                                 !< area of the tessera
    FLOAT :: normal(1:3)                          !< unitary outgoing vector normal to the tessera surface  
    FLOAT :: r_sphere                             !< radius of the sphere to which the tessera belongscd
  end type pcm_tessera_t

  type(pcm_tessera_t), allocatable :: cts_act(:)  !< tesselation arrays (nts_act)
  integer :: nts_act				                      !< number of tesserae

  !> set of parameters for Debye dielectric model
  type :: debye_param_t
    FLOAT :: eps_0                                !< static  dielectric constant (at ZERO     frequency of the field)
    FLOAT :: eps_d                                !< dynamic dielectric constant (at INFINITE frequency of the field)
    FLOAT :: tau                                  !< Debye relaxation time
  end type debye_param_t

  !> set of parameters for Drude-Lorentz dielectric model
  !> no parameter allows the match with the dynamic dielectric constant
  type :: drude_param_t
    FLOAT :: aa                                   !< parameter matching the static dielectric constant
    FLOAT :: gm                                   !< damping of the plasma oscillation
    FLOAT :: w0                                   !< resonance frequency
  end type drude_param_t

  type(debye_param_t) :: deb
  type(drude_param_t) :: drl

  character(8) :: which_eom			                  !< character flag for PCM charges due to:
                                                  !< electrons ('electron'), external potential ('external') or kick ('justkick')

  character(3) :: which_eps			                  !< character flag for Debye ('deb') and Drude-Lorentz models ('drl')

  FLOAT :: dt 					                          !< time-step of the propagation

						                                      !< polarization charges and variation of it in the previous iteration:
  FLOAT, allocatable :: q_tp(:), dq_tp(:)	        !< 			 due to solute electrons
  FLOAT, allocatable :: qext_tp(:), dqext_tp(:)	  !< 			 due to external potential
  FLOAT, allocatable :: qkick_tp(:), dqkick_tp(:)	!< 			 due to kick

  FLOAT, allocatable :: pot_tp(:)		              !< Hartree (electrons) potential in previous iteration
  FLOAT, allocatable :: potext_tp(:)		          !< external potential in previous iteration

                                                  !< See Chem.Phys.Lett. 429 (2006) 310-316 for Velocity-Verlet (VV) algorithm... 
  FLOAT :: f1, f2, f3, f4, f5			                !< auxiliar constants for VV
                                                  !< analogous to force in the equation of motion for the pol.charges, prev. iter.
  FLOAT, allocatable :: force_tp(:)	              !< 			 due to solute electrons
  FLOAT, allocatable :: force_qext_tp(:)          !< 			 due to external potential
  FLOAT, allocatable :: force_qkick_tp(:)         !< 			 due to kick

  						                                    !< In Ref.1 - J.Phys.Chem.A 2015, 119, 5405-5416... - Debye dielectric func. (eps)
                                                  !< In Ref.2 - In J. Phys. Chem. C 2016, 120, 28774−28781... - Drude-Lorentz eps
  FLOAT, allocatable :: cals(:,:),cald(:,:)	      !< Calderon matrices S and D from Eq.(5), Ref.1
  FLOAT, allocatable :: eigv(:),eigt(:,:)         !< \Lambda and T matrices from Eq.(10), Ref.1
  FLOAT, allocatable :: sm12(:,:),sp12(:,:)       !< S^{-1/2} and S^{1/2}
                                                  !< Q^{IEF(d)}_0 (not used in ref.) and Q^{IEF(d)}_d from Eq.(18) with eps_0/eps_d
  FLOAT, allocatable :: matq0(:,:),matqd(:,:)       !< 			 for solute
  FLOAT, allocatable :: matq0_lf(:,:),matqd_lf(:,:) !< 			 for external potential
                                                  !< \tilde{Q} and R matrices from Eq.(38)-(39) (Ref.1), respectively
                                                  !<  Q_f and Q_{\omega} from Eq.(17)-(16) (Ref.2), respectively
  FLOAT, allocatable :: matqv(:,:),matqq(:,:)     !< 			 for solute
  FLOAT, allocatable :: matqv_lf(:,:)							!< 			 for external potential
                                                  !< Q^{IEF(d)}_d, \tilde{Q} and R matrices enter the EOM in eq.(37), Ref.1
                                                  !< Q_f and Q_{\omega} matrices enter the EOM in eq.(15), Ref.2
                                                  !< N.B.: matrices R (in case of Debye) and Q_{\omega} (in case of Drude-Lorentz), 
                                                  !< i.e., matqq in this implementation, 
                                                  !< are the same in the EOMs for polarization charges due to the solute or ext.pot.

  !> mathematical constants
  FLOAT, parameter :: twopi=M_TWO*M_Pi
  FLOAT, parameter :: fourpi=M_FOUR*M_Pi

  logical :: enough_initial = .false.

  contains

  !------------------------------------------------------------------------------------------------------------------------------
  !> Driving subroutine for the Equation of Motion (EOM) propagation of the polarization charges
  !> within the Integral Equation Formalism (IEF) formulation of the Polarization Continuum Model (PCM).
  subroutine pcm_charges_propagation(q_t, pot_t, this_dt, this_cts_act, input_asc, this_eom, this_eps, this_deb, this_drl) 
   save
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   FLOAT, intent(in) :: this_dt

   type(pcm_tessera_t), intent(in) :: this_cts_act(:)

   logical, intent(in) :: input_asc

   character(*), intent(in) :: this_eom !< EOM case, either due to electrons ('electron') or due to external potential ('external')
   character(*), intent(in) :: this_eps !< type of dielectric model to be used, either Debye ('deb') or Drude-Lorentz ('drl')
   type(debye_param_t), optional, intent(in) :: this_deb
   type(drude_param_t), optional, intent(in) :: this_drl

   logical :: firsttime=.true.
   logical :: initial_electron=.true.
   logical :: initial_external=.true.
   logical :: initial_kick=.true.

   PUSH_SUB(pcm_charges_propagation)

   which_eom=this_eom
   if (which_eom /= 'electron' .and. which_eom /= 'external' .and. which_eom /= 'justkick' ) then
    message(1) = "pcm_charges_propagation: EOM evolution of PCM charges can only be due to solute electrons"
    message(2) = "or external potential (including the kick)."
    call messages_fatal(2)
   endif

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
    if (which_eps=='deb' .and. (.not.(present(this_deb)))) then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Debye dielectric function requires three parameters."
     call messages_fatal(1)
    else if (which_eps=='deb' .and. present(this_deb)) then
     deb=this_deb
    endif
    if (which_eps=='drl' .and. (.not.(present(this_drl)))) then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Drude-Lorentz dielectric function requires three parameters."
     call messages_fatal(1)
    else if (which_eps=='drl' .and. present(this_drl)) then
     drl=this_drl
    endif
    if (which_eps/='deb' .and. which_eps/='drl') then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Only Debye or Drude-Lorent dielectric models are allowed."
     call messages_fatal(1)
    endif
    if( abs(deb%tau) <= M_EPSILON ) then
     message(1) = "pcm_charges_propagation: EOM-PCM error. Debye EOM-PCM require a non-null Debye relaxation time."
     call messages_fatal(1)
    endif
    firsttime=.false.
   endif


   if( input_asc .and. which_eps == 'deb' ) then

    !> initialize pcm charges due to electrons, external potential or kick

    call pcm_charges_from_input_file(q_t,pot_t)

    POP_SUB(pcm_charges_propagation)
    return

   else if( input_asc .and. which_eps /= 'deb' ) then

    message(1) = "pcm_charges_propagation: EOM-PCM error. Only Debye EOM-PCM can startup from input charges."
    call messages_fatal(1)

   endif

   if( ( initial_electron .and. which_eom == 'electron' ) .or. & 
       ( initial_external .and. which_eom == 'external' ) .or. & 
       ( initial_kick     .and. which_eom == 'justkick' )      ) then

     !> initialize pcm matrices

     call pcm_bem_init

     !> initialize pcm charges due to electrons, external potential or kick

     call init_charges(q_t, pot_t)

     if(initial_electron .and. which_eom == 'electron' ) initial_electron=.false.
     if(initial_external .and. which_eom == 'external' ) initial_external=.false.
     if(initial_kick     .and. which_eom == 'justkick' ) initial_kick=.false.

   else

     !> propagate pcm charges due to electrons or external potential (including possible kick)

     if( which_eps .eq. "deb") then
      call pcm_ief_prop_deb(q_t,pot_t)
     else if( which_eps .eq. "drl" ) then
      call pcm_ief_prop_vv_ief_drl(q_t,pot_t)
     endif

   endif
  
   POP_SUB(pcm_charges_propagation)
  end subroutine pcm_charges_propagation

  !------------------------------------------------------------------------------------------------------------------------------
  !> Polarization charges initialization from input file
  subroutine pcm_charges_from_input_file(q_t, pot_t)
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   FLOAT :: aux1(3)
   integer :: aux2

   integer :: asc_unit

   integer :: ia

   PUSH_SUB(pcm_charges_from_input_file)

   if( which_eom == 'electron' ) then
    SAFE_ALLOCATE(q_tp(nts_act))
    asc_unit = io_open(PCM_DIR//'ASC_e.dat', action='read')
    do ia = 1, nts_act
      read(asc_unit,*) aux1, q_t(ia), aux2
    end do
    q_tp=q_t
   else if( which_eom == 'external' ) then
    SAFE_ALLOCATE(qext_tp(nts_act))
    asc_unit = io_open(PCM_DIR//'ASC_ext.dat', action='read')
    do ia = 1, nts_act
      read(asc_unit,*) aux1, q_t(ia), aux2
    end do
    qext_tp = q_t
   else if( which_eom == 'justkick' ) then
    SAFE_ALLOCATE(qkick_tp(nts_act))
    asc_unit = io_open(PCM_DIR//'ASC_kick.dat', action='read')
    do ia = 1, nts_act
      read(asc_unit,*) aux1, q_t(ia), aux2
    end do
    qkick_tp = q_t  
   endif
   call io_close(asc_unit)   

   SAFE_ALLOCATE(pot_tp(nts_act))
   pot_tp = pot_t

   POP_SUB(pcm_charges_from_input_file)
  end subroutine pcm_charges_from_input_file

  !------------------------------------------------------------------------------------------------------------------------------
  !> Polarization charges initialization (in equilibrium with the initial potential for electrons)
  subroutine init_charges(q_t,pot_t)
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   PUSH_SUB(init_charges)

   if( which_eom == 'electron' ) then
    !< Here we consider the potential at any earlier time equal to the initial potential.
    !< Therefore, we can suppose that the solvent is already in equilibrium with the initial solute potential.
    !< The latter is valid when starting the propagation from the ground state but does not hold in general.
    message(1) = 'EOM-PCM for solvent polarization due to solute electrons considers that you start from a ground state run.'
    call messages_warning(1)
  
    SAFE_ALLOCATE(pot_tp(nts_act))
    pot_tp = pot_t

    !< applying the static IEF-PCM response matrix (corresponging to epsilon_0) to the initial potential
    q_t=matmul(matq0,pot_t) 

    SAFE_ALLOCATE(q_tp(nts_act))
    q_tp = q_t

    if( which_eps .eq. "drl" ) then
     SAFE_ALLOCATE(dq_tp(nts_act))
     SAFE_ALLOCATE(force_tp(nts_act))
     dq_tp=M_ZERO
     force_tp=M_ZERO
    endif
			    
   else if( which_eom == 'external' ) then
    !< Here (instead) we consider zero the potential at any earlier time.
    !< Therefore, the solvent is not initially in equilibrium with the external potential unless its initial value is zero.

    SAFE_ALLOCATE(potext_tp(nts_act))
    potext_tp = pot_t

    !< applying the dynamic IEF-PCM response matrix (for epsilon_d) to the initial potential
    q_t=matmul(matqd_lf,pot_t)

    SAFE_ALLOCATE(qext_tp(nts_act))
    qext_tp = q_t

    if( which_eps .eq. "drl" ) then
     SAFE_ALLOCATE(dqext_tp(nts_act))
     SAFE_ALLOCATE(force_qext_tp(nts_act))
     dqext_tp=M_ZERO
     force_qext_tp=M_ZERO
    endif

   else if( which_eom == 'justkick' ) then

    if( which_eps .eq. "drl" ) then
     SAFE_ALLOCATE(dqkick_tp(nts_act))
     SAFE_ALLOCATE(force_qkick_tp(nts_act))
     dqkick_tp=M_ZERO
     force_qkick_tp=M_ZERO
    endif

    if( which_eps .eq. "drl" ) then
      dqkick_tp=matmul(matqv_lf,pot_t)*dt
    else
      q_t = matmul(matqv_lf-matmul(matqq,matqd_lf),pot_t)
    endif

    SAFE_ALLOCATE(qkick_tp(nts_act))
    qkick_tp = q_t

   endif

   !< initializing Velocity-Verlet algorithm for the integration of EOM-PCM for Drude-Lorentz
   if( which_eps == 'drl' ) call init_vv_propagator

   POP_SUB(init_charges)
  end subroutine init_charges

  !------------------------------------------------------------------------------------------------------------------------------
  !> Euler method for integrating first order EOM for the polarization charges within IEF-PCM
  !> in the case of Debye dielectric functions.
  subroutine pcm_ief_prop_deb(q_t,pot_t)
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   PUSH_SUB(pcm_ief_prop_deb)

   if( which_eom == 'electron' ) then
    !> Eq.(47) in S. Corni, S. Pipolo and R. Cammi, J.Phys.Chem.A 2015, 119, 5405-5416.
    q_t(:) = q_tp(:) - dt*matmul(matqq,q_tp)   + &
	  	       dt*matmul(matqv,pot_tp) + &
			  matmul(matqd,pot_t-pot_tp)

    q_tp = q_t

    pot_tp = pot_t

   else if( which_eom == 'external' ) then 
    !> analogous to Eq.(47) ibid. with different matrices
    q_t(:) = qext_tp(:) - dt*matmul(matqq,qext_tp)   + &
	  	          dt*matmul(matqv_lf,potext_tp) + &
			     matmul(matqd_lf,pot_t-potext_tp)             !< N.B. matqq

    qext_tp = q_t

    potext_tp = pot_t

   else if( which_eom == 'justkick' ) then
    !> simplifying for kick
    q_t(:) = qkick_tp(:) - dt*matmul(matqq,qkick_tp)    !< N.B. matqq

    qkick_tp = q_t

   endif

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
   f3=1.d0-dt*drl%gm*(1.d0-dt*M_HALF*drl%gm)
   f4=0.5d0*dt
   f5=drl%gm*f2

   POP_SUB(init_vv_propagator)
  end subroutine init_vv_propagator

  !------------------------------------------------------------------------------------------------------------------------------
  !> VV algorithm for integrating second order EOM for the polarization charges within IEF-PCM
  !> in the case of Drude-Lorentz dielectric functions.
  subroutine pcm_ief_prop_vv_ief_drl(q_t,pot_t)
   FLOAT, intent(out) :: q_t(:)
   FLOAT, intent(in)  :: pot_t(:)

   FLOAT :: force(nts_act)

   PUSH_SUB(pcm_ief_prop_vv_ief_drl)

   if( which_eom == 'electron' ) then 
    !> From Eq.(15) S. Pipolo and S. Corni, J.Phys.Chem.C 2016, 120, 28774-28781.
    !> Using the scheme in Eq.(21) and (17) of E. Vanden-Eijnden, G. Ciccotti, Chem.Phys.Lett. 429 (2006) 310-316.
    q_t = q_tp + f1*dq_tp + f2*force_tp
    force = -matmul(matqq,q_t) + matmul(matqv,pot_t)
    dq_tp = f3*dq_tp + f4*(force+force_tp) -f5*force_tp
    force_tp = force
    q_tp = q_t

   else if( which_eom == 'external' ) then
    !> analagous to Eq.(15) ibid. with different matrices
    q_t = qext_tp + f1*dqext_tp + f2*force_qext_tp
    force = -matmul(matqq,q_t) + matmul(matqv_lf,pot_t)	!< N.B. matqq
    dqext_tp = f3*dqext_tp + f4*(force+force_qext_tp) -f5*force_qext_tp
    force_qext_tp = force
    q_tp = q_t

   else if( which_eom == 'justkick' ) then
    !> simplifying for kick
    q_t = qkick_tp + f1*dqkick_tp + f2*force_qkick_tp
    force = -matmul(matqq,q_t)												  !< N.B. matqq
    dqkick_tp = f3*dqkick_tp + f4*(force+force_qkick_tp) -f5*force_qkick_tp
    force_qkick_tp = force
    q_tp = q_t

   endif

   pot_tp = pot_t

   POP_SUB(pcm_ief_prop_vv_ief_drl)
  end subroutine pcm_ief_prop_vv_ief_drl

  !------------------------------------------------------------------------------------------------------------------------------
  !> Boundary Element Method (BEM) EOM-IEF-PCM matrices initialization.
  subroutine pcm_bem_init

   integer :: itess, jtess
   integer :: pcmmat0_unit, pcmmatd_unit

   PUSH_SUB(pcm_bem_init)

   if( which_eom == 'electron' ) then
    SAFE_ALLOCATE(matq0(nts_act,nts_act))
    SAFE_ALLOCATE(matqd(nts_act,nts_act))
    SAFE_ALLOCATE(matqv(nts_act,nts_act))
    if( .not.allocated(matqq) ) then
     SAFE_ALLOCATE(matqq(nts_act,nts_act))
    endif
   else if( which_eom == 'external' .or. which_eom == 'justkick' ) then
    SAFE_ALLOCATE(matq0_lf(nts_act,nts_act)) !< not used yet
    SAFE_ALLOCATE(matqd_lf(nts_act,nts_act))
    SAFE_ALLOCATE(matqv_lf(nts_act,nts_act))
    if( .not.allocated(matqq) ) then
     SAFE_ALLOCATE(matqq(nts_act,nts_act))
    endif
   endif  
   call do_PCM_propMat

   if( which_eom == 'electron' ) then
    pcmmat0_unit = io_open(PCM_DIR//'pcm_matrix_static_from_eom.out', action='write')
    pcmmatd_unit = io_open(PCM_DIR//'pcm_matrix_dynamic_from_eom.out', action='write')
     do jtess = 1, nts_act
      do itess = 1, nts_act
       write(pcmmat0_unit,*) matq0(itess,jtess)
       write(pcmmatd_unit,*) matqd(itess,jtess)
      end do
     end do
    call io_close(pcmmat0_unit)
    call io_close(pcmmatd_unit)
   else if( which_eom == 'external' ) then
    pcmmat0_unit = io_open(PCM_DIR//'pcm_matrix_static_lf_from_eom.out', action='write')
    pcmmatd_unit = io_open(PCM_DIR//'pcm_matrix_dynamic_lf_from_eom.out', action='write')
     do jtess = 1, nts_act
      do itess = 1, nts_act
       write(pcmmat0_unit,*) matq0_lf(itess,jtess)
       write(pcmmatd_unit,*) matqd_lf(itess,jtess)
      end do
     end do
    call io_close(pcmmat0_unit)
    call io_close(pcmmatd_unit)
   endif

   POP_SUB(pcm_bem_init)
  end subroutine pcm_bem_init

  !------------------------------------------------------------------------------------------------------------------------------

  subroutine pcm_eom_end
   PUSH_SUB(pcm_eom_end)

   ! pcm charges
   if( allocated(q_tp) ) then
    SAFE_DEALLOCATE_A(q_tp)
   endif
   if( allocated(qext_tp) ) then
    SAFE_DEALLOCATE_A(qext_tp)
   endif
   if( allocated(qkick_tp) ) then
    SAFE_DEALLOCATE_A(qkick_tp)
   endif

   ! increment pcm charges
   if( allocated(dq_tp) ) then
    SAFE_DEALLOCATE_A(dq_tp)
   endif
   if( allocated(dqext_tp) ) then
    SAFE_DEALLOCATE_A(dqext_tp)
   endif
   if( allocated(dqkick_tp) ) then
    SAFE_DEALLOCATE_A(dqkick_tp)
   endif

   ! force-like term in pcm-eom equation
   if( allocated(force_tp) ) then
    SAFE_DEALLOCATE_A(force_tp)
   endif
   if( allocated(force_qext_tp) ) then
    SAFE_DEALLOCATE_A(force_qext_tp)
   endif
   if( allocated(force_qkick_tp) ) then
    SAFE_DEALLOCATE_A(force_qkick_tp)
   endif

   ! pcm-eom bem matrices
   if( allocated(matq0) ) then
    SAFE_DEALLOCATE_A(matq0)
   endif
   if( allocated(matqd) ) then
    SAFE_DEALLOCATE_A(matqd)
   endif
   if( allocated(matq0_lf) ) then
    SAFE_DEALLOCATE_A(matq0_lf)
   endif
   if( allocated(matqd_lf) ) then
    SAFE_DEALLOCATE_A(matqd_lf)
   endif
   if( allocated(matqv) ) then
    SAFE_DEALLOCATE_A(matqv)
   endif
   if( allocated(matqv_lf) ) then
    SAFE_DEALLOCATE_A(matqv_lf)
   endif
   if( allocated(matqq) ) then
    SAFE_DEALLOCATE_A(matqv)
   endif

   POP_SUB(pcm_eom_end)
  end subroutine pcm_eom_end

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build the required BEM matrices for the EOM-IEF-PCM for Debye and Drude-Lorentz cases.
  !> Following from
  !> Ref.1 - J.Phys.Chem.A 2015, 119, 5405-5416.    - Debye case
  !> Ref.2 - J.Phys.Chem.C 2016, 120, 28-774-28781. - Drude-Lorentz case
  !> The matrices are required for the EOMs, eq.(37) in Ref.1 and eq.(15) in Ref.2 
  subroutine do_PCM_propMat
   save
   integer :: i,j
   FLOAT, allocatable :: scr4(:,:),scr1(:,:)
   FLOAT, allocatable :: scr2(:,:),scr3(:,:)
   FLOAT, allocatable :: fact1(:),fact2(:)   !< tau^{-1} and tau^{-1}K_0 from (38)-(39), respectively, with Eq.(32), Ref.1
                                             !< K_{\omega} and K_f from (19)-(10), respectively, Ref.2
   FLOAT, allocatable :: Kdiag0(:),Kdiagd(:) !< correspond to K_0 (static dielec. const.) and K_d (dynamic dielec. const.) in Ref.1
                                             !< e.g. Eq.(38) -^                               ^- e.g. implicitly in Q_d^{IEF(d)}, 
                                             !<                                                       see Eq.(18) and (37)
                                             !< not used in Ref.2
   FLOAT :: sgn,sgn_lf,fac_eps0,fac_epsd
   FLOAT :: temp

   logical :: firsttime=.true.

   PUSH_SUB(do_PCM_propMat)

   sgn=M_ONE ! In the case of NP this is -1 because the normal to the cavity is always pointing outward 
             ! and the field created by a positive unit charge outside the cavity is directed inward.

   sgn_lf=M_ONE
   !> 'local field' differ from 'standard' PCM response matrix in some sign changes
   if( which_eom == 'external' .or. which_eom == 'justkick' ) sgn_lf=-M_ONE 

   SAFE_ALLOCATE(cals(nts_act,nts_act))
   SAFE_ALLOCATE(cald(nts_act,nts_act))
   SAFE_ALLOCATE(Kdiag0(nts_act))
   SAFE_ALLOCATE(Kdiagd(nts_act))
   SAFE_ALLOCATE(fact1(nts_act))
   SAFE_ALLOCATE(fact2(nts_act))
   !> generate Calderon S and D matrices
   do i=1,nts_act
    do j=1,nts_act
     call green_s(i,j,temp)
     cals(i,j)=temp
     call green_d(i,j,temp)
     cald(i,j)=temp
    enddo
   enddo
   if ( firsttime ) then 
    call allocate_TS_matrix
    call do_TS_matrix
    firsttime=.false.
   endif
   SAFE_DEALLOCATE_A(cals)
   SAFE_DEALLOCATE_A(cald)
   SAFE_ALLOCATE(scr4(nts_act,nts_act))
   SAFE_ALLOCATE(scr1(nts_act,nts_act))
   SAFE_ALLOCATE(scr2(nts_act,nts_act))
   SAFE_ALLOCATE(scr3(nts_act,nts_act))
   if (which_eps .eq. "deb") then
    if( deb%eps_0 /= M_ONE ) then
     fac_eps0=(deb%eps_0+M_ONE)/(deb%eps_0-M_ONE)			 
     Kdiag0(:)=sgn_lf*(twopi-sgn*sgn_lf*eigv(:))/(twopi*fac_eps0-sgn*eigv(:)) !< Eq.(14) with eps_0 in Ref.1
    else
     Kdiag0(:)=M_ZERO
    endif
    if( deb%eps_d /= M_ONE ) then
     fac_epsd=(deb%eps_d+M_ONE)/(deb%eps_d-M_ONE)
     Kdiagd(:)=sgn_lf*(twopi-sgn*sgn_lf*eigv(:))/(twopi*fac_epsd-sgn*eigv(:)) !< Eq.(14) with eps_d, ibid.
    else
     Kdiagd(:)=M_ZERO
    endif
    fact1(:)=((twopi-sgn*eigv(:))*deb%eps_0+twopi+eigv(:))/ &                 !< inverse of Eq.(32), ibid.
    ((twopi-sgn*eigv(:))*deb%eps_d+twopi+eigv(:))/deb%tau
    fact2(:)=Kdiag0(:)*fact1(:)				                                        !< tau^{-1}K_0 in Eq.(38), ibid.
   elseif (which_eps .eq. "drl") then
    Kdiagd(:)=M_ZERO					                                                !< from Eq.(10) up in Ref.2
    fact2(:)=(twopi-sgn*eigv(:))*drl%aa/fourpi                                !< Eq.(10) down
    do i=1,nts_act
     if(fact2(i).lt.M_ZERO) fact2(i)=M_ZERO !< check out
    enddo
    if (abs(drl%w0)<=M_EPSILON) drl%w0=1.d-8      !< check out
    fact1(:)=fact2(:)+drl%w0*drl%w0			                                      !< Eq.(19), ibid.
    fact2(:)=sgn_lf*(twopi-sgn*sgn_lf*eigv(:))*drl%aa/fourpi                  !< Eq.(10) down, local field analogous
    Kdiag0(:)=fact2(:)/fact1(:)				                                        !< from Eq.(10) up, ibid.
   endif
   scr3=matmul(sm12,eigt)
   scr2=matmul(transpose(eigt),sp12)
   scr4=transpose(scr3)
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact1(i)
   enddo
   matqq=matmul(scr1,scr2)				                                            !< Eq.(39) in Ref.1 and Eq.(16) in Ref.2
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*fact2(i)
   enddo
   if( which_eom == 'electron' ) then
    matqv=-matmul(scr1,scr4)				                                          !< Eq.(38) in Ref.1 and Eq.(17) in Ref.2
   else if( which_eom == 'external' .or. which_eom == 'justkick' ) then
    matqv_lf=-matmul(scr1,scr4)				                                        !< local field analogous
   endif
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiag0(i)
   enddo
   if( which_eom == 'electron' ) then
    matq0=-matmul(scr1,scr4)				                                          !< from Eq.(14) and (18) for eps_0 in Ref.1
   else if( which_eom == 'external' .or. which_eom == 'justkick' ) then
    matq0_lf=-matmul(scr1,scr4)			                                          !< local field analogous !< not used yet
   endif
   do i=1,nts_act
    scr1(:,i)=scr3(:,i)*Kdiagd(i)
   enddo
   if( which_eom == 'electron' ) then
    matqd=-matmul(scr1,scr4)				                                          !< from Eq.(14) and (18) for eps_d, ibid.
   else if( which_eom == 'external' .or. which_eom == 'justkick' ) then
    matqd_lf=-matmul(scr1,scr4)				                                        !< local field analogous
   endif
   SAFE_DEALLOCATE_A(scr4)
   SAFE_DEALLOCATE_A(scr1)
   SAFE_DEALLOCATE_A(scr2)
   SAFE_DEALLOCATE_A(scr3)
   SAFE_DEALLOCATE_A(fact1)
   SAFE_DEALLOCATE_A(fact2)
   SAFE_DEALLOCATE_A(Kdiag0)
   SAFE_DEALLOCATE_A(Kdiagd)

   if( enough_initial ) then
    call deallocate_TS_matrix
   endif

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
    value=dot_product(cts_act(j)%normal,diff)/dist**3       !< Eq.(5) in Refs.1-2
   else
    value=CNST(-1.0694)*sqrt(fourpi*cts_act(i)%area)
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
    value=M_ONE/dist					                !< Eq.(5) in Refs.1-2
   else
    value=CNST(1.0694)*sqrt(fourpi/cts_act(i)%area) !< diagonal part is a bit cumbersome
   endif

   POP_SUB(green_s)
  end subroutine green_s

  !------------------------------------------------------------------------------------------------------------------------------
  subroutine allocate_TS_matrix
   PUSH_SUB(allocate_TS_matrix)

   SAFE_ALLOCATE(eigv(nts_act))
   SAFE_ALLOCATE(eigt(nts_act,nts_act))
   SAFE_ALLOCATE(sm12(nts_act,nts_act))
   SAFE_ALLOCATE(sp12(nts_act,nts_act))

   POP_SUB(allocate_TS_matrix)
  end subroutine allocate_TS_matrix

  !------------------------------------------------------------------------------------------------------------------------------
  subroutine deallocate_TS_matrix
   PUSH_SUB(deallocate_TS_matrix)

   SAFE_DEALLOCATE_A(eigv)
   SAFE_DEALLOCATE_A(eigt)
   SAFE_DEALLOCATE_A(sm12)
   SAFE_DEALLOCATE_A(sp12)

   POP_SUB(deallocate_TS_matrix)
  end subroutine deallocate_TS_matrix

  !------------------------------------------------------------------------------------------------------------------------------
  !> Subroutine to build matrices S^{1/2}, S^{-1/2}, T and \Lambda (notation of Refs.1-2)
  subroutine do_TS_matrix
   integer :: i,j
   integer :: info,lwork,liwork
   FLOAT, allocatable :: scr1(:,:),scr2(:,:),eigt_t(:,:)
   FLOAT :: sgn
   character jobz,uplo
   integer, allocatable :: iwork(:)
   FLOAT,allocatable :: work(:)

   PUSH_SUB(do_TS_matrix)

   SAFE_ALLOCATE(scr1(nts_act,nts_act))
   SAFE_ALLOCATE(scr2(nts_act,nts_act))
   SAFE_ALLOCATE(eigt_t(nts_act,nts_act))
   SAFE_ALLOCATE(work(1+6*nts_act+2*nts_act*nts_act))
   SAFE_ALLOCATE(iwork(3+5*nts_act))

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
   sp12=matmul(scr1,eigt_t) 		         !< building S^{1/2} to be used in R (Ref.1) and Q_{\omega} (Ref.2)
   do i=1,nts_act
    scr1(:,i)=eigt(:,i)/sqrt(eigv(i))
   enddo
   sm12=matmul(scr1,eigt_t) 		         !< building S^{-1/2} to be used in R and \tilde{Q} (Ref.1) and Q_{\omega} and Q_f (Ref.2)
   do i=1,nts_act
    scr1(:,i)=cald(:,i)*cts_act(i)%area
   enddo
   scr2=matmul(matmul(sm12,scr1),sp12)   !< Eq.(10) in Ref.1 (paragraph after Eq.(9), Ref.2)
   do j=1,nts_act
    do i=1,nts_act
     eigt(i,j)=M_HALF*(scr2(i,j)+scr2(j,i)) !< re-symmetrizing for numerical reasons
    enddo
   enddo
   call dsyevd(jobz,uplo,nts_act,eigt,nts_act,eigv,work,lwork, &
               iwork,liwork,info)        !< obtaining \Lambda (eigv) and T (eigt), Eq.(10), ibid.

   SAFE_DEALLOCATE_A(work)
   SAFE_DEALLOCATE_A(iwork)
   SAFE_DEALLOCATE_A(scr1)
   SAFE_DEALLOCATE_A(scr2)
   SAFE_DEALLOCATE_A(eigt_t)

   POP_SUB(do_TS_matrix)
  end subroutine do_TS_matrix

  ! -----------------------------------------------------------------------------
  subroutine pcm_eom_enough_initial(not_yet_called)
    logical, intent(out) :: not_yet_called
    PUSH_SUB(pcm_eom_enough_initial)
    enough_initial = .true.
    not_yet_called = .false.
    POP_SUB(pcm_eom_enough_initial)
  end subroutine pcm_eom_enough_initial

end module pcm_eom_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

