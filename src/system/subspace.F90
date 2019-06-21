!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module subspace_oct_m
  use accel_oct_m
  use accel_blas_oct_m
  use batch_oct_m
  use batch_ops_oct_m
  use blas_oct_m
  use blacs_oct_m
  use blacs_proc_grid_oct_m
  use comm_oct_m
  use derivatives_oct_m
#ifdef HAVE_ELPA
  use elpa
#endif
  use global_oct_m
  use hamiltonian_oct_m
  use hardware_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_batch_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use pblas_oct_m
  use profiling_oct_m
  use scalapack_oct_m
  use sort_oct_m
  use states_oct_m
  use states_calc_oct_m
  use states_parallel_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none
  
  private

  public ::               &
    subspace_t,           &
    subspace_init,        &
    dsubspace_diag,       &
    zsubspace_diag, &
    zsubspace_diag_hamiltonian, &
    zfloquet_FBZ_subspace_diag

  type subspace_t
    integer :: method
  end type subspace_t

  type(profile_t),     save    :: diagon_prof, hamiltonian_prof
  
contains

  subroutine subspace_init(this, st, no_sd)
    type(subspace_t),  intent(out) :: this
    type(states_t),    intent(in)  :: st
    logical,           intent(in)  :: no_sd

    integer :: default

    PUSH_SUB(subspace_init)

    if(no_sd) then

      this%method = OPTION__SUBSPACEDIAGONALIZATION__NONE

    else

      !%Variable SubspaceDiagonalization
      !%Type integer
      !%Default standard
      !%Section SCF::Eigensolver
      !%Description
      !% Selects the method to perform subspace diagonalization. The
      !% default is <tt>standard</tt>, unless states parallelization is used,
      !% when the default is <tt>scalapack</tt>.
      !%Option none 0
      !% No subspace diagonalization. WARNING: this will generally give incorrect results.
      !%Option standard 1
      !% The standard routine. Can be used with domain parallelization but not
      !% state parallelization.
      !%Option scalapack 3
      !% State-parallelized version using ScaLAPACK. (Requires that
      !% Octopus was compiled with ScaLAPACK support.)
      !%Option floquet_ss 4
      !% (default and only option for floquet calculation mode, do not use!)
      !%End

      default = OPTION__SUBSPACEDIAGONALIZATION__STANDARD

#ifdef HAVE_SCALAPACK
      if(st%parallel_in_states) default = OPTION__SUBSPACEDIAGONALIZATION__SCALAPACK
#endif

      call parse_variable('SubspaceDiagonalization', default, this%method)

      if(.not.varinfo_valid_option('SubspaceDiagonalization', this%method)) call messages_input_error('SubspaceDiagonalization')
    end if

    call messages_print_var_option(stdout, 'SubspaceDiagonalization', this%method)

    ! some checks for ingenious users
    if(this%method == OPTION__SUBSPACEDIAGONALIZATION__SCALAPACK) then
#ifndef HAVE_MPI
      message(1) = 'The scalapack subspace diagonalization can only be used in parallel.'
      call messages_fatal(1, only_root_writes = .true.)
#else
#ifndef HAVE_SCALAPACK
      message(1) = 'The scalapack subspace diagonalization requires scalapack.'
      call messages_fatal(1, only_root_writes = .true.)
#endif
      if(st%dom_st_mpi_grp%size == 1) then
        message(1) = 'The scalapack subspace diagonalization is designed to be used with domain or state parallelization.'
        call messages_warning(1)
      end if

      if(st%d%kpt%parallel) then
        message(1) = 'Currently the scalapack subspace diagonalization does not use k-point parallelization.'
        call messages_warning(1)
      end if
#endif
    end if

    POP_SUB(subspace_init)
  end subroutine subspace_init




  subroutine zfloquet_FBZ_subspace_diag(der, st,  hm, ik, start)
    type(derivatives_t),    intent(in)    :: der
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(inout) :: hm
    integer,                intent(in)    :: ik
    logical, optional,      intent(in)    :: start

    CMPLX, allocatable :: evecs(:, :), rdiff(:), H0(:,:), one(:,:), HF(:,:), P(:,:), Pd(:,:), Heff(:,:)
    CMPLX, allocatable :: rot_state(:,:,:), state(:,:),  evecs_temp(:,:)
    FLOAT, allocatable :: evalues_full(:), mix(:),norms(:), evalues0(:), evals_diff(:,:), evecs_reshape(:,:)
    integer            :: block0, nst0, maxiter, nmix, ist,jst, im,in, jj,ii, Fdim, it, pos, idim, imm, inn
    integer, allocatable :: idx(:)
    logical            :: start_
    FLOAT, allocatable :: overlap(:)    
FLOAT :: temp    
    PUSH_SUB(floquet_FBZ_subspace_diag)

    ASSERT(hm%F%order(1) == -hm%F%order(2))

    start_ = optional_default(start,.false.)

    Fdim=hm%F%floquet_dim

    SAFE_ALLOCATE(H0(1:st%nst, 1:st%nst)) ! the zero-block, this could be gs Hamiltonian(?)
    SAFE_ALLOCATE( P(1:st%nst, 1:st%nst)) ! the off-diagonal intereaction block
    SAFE_ALLOCATE(Pd(1:st%nst, 1:st%nst)) ! P^\dagger  
    SAFE_ALLOCATE(HF(1:st%nst*Fdim, 1:st%nst*Fdim)) ! the full Floquet matrix
    SAFE_ALLOCATE(one(1:st%nst, 1:st%nst)) 
    SAFE_ALLOCATE(evecs(1:st%nst*Fdim, 1:st%nst))
SAFE_ALLOCATE(evecs_temp(1:st%nst*Fdim, 1:st%nst))
    SAFE_ALLOCATE(evalues_full(1:st%nst*Fdim))
    SAFE_ALLOCATE(evalues0(1:st%nst))
    SAFE_ALLOCATE(Heff(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(evals_diff(1:st%nst*Fdim, 1:st%nst))
    SAFE_ALLOCATE(evecs_reshape(1:st%nst,1:Fdim))
    SAFE_ALLOCATE(norms(1:st%nst*Fdim))
    SAFE_ALLOCATE(idx(1:st%nst*Fdim))
    SAFE_ALLOCATE(overlap(1:st%nst*Fdim))
    evecs = M_ZERO

    one = M_ZERO
    do ist=1,st%nst
      one(ist,ist) = M_ONE
    end do

    call sort(st%eigenval(1:st%nst,ik), idx)

    hm%F%floquet_apply = .false.
    hm%d%dim = hm%F%spindim !st%d%dim
    ! get the action of H0 on the states
    call zsubspace_diag_hamiltonian(der,  st, hm, ik, H0)
    P = M_ZERO
    do it=1,hm%F%nT
      Pd = M_ZERO
      call zsubspace_diag_hamiltonian(der, st, hm%td_hm(it), ik, Pd)
      P = P + Pd*exp(-M_zI*hm%F%omega*it*hm%F%dt)/hm%F%nT
    end do
    Pd(1:st%nst, 1:st%nst) = transpose(conjg(P))
    ! construct the full matrix by shifting the H0 action
    HF = M_ZERO
    do im=1,Fdim
      HF((im-1)*st%nst+1:im*st%nst , (im-1)*st%nst+1:im*st%nst) = H0+(im-hm%F%order(2)-1)*hm%F%omega*one
      if(im>1) HF((im-2)*st%nst+1:(im-1)*st%nst , (im-1)*st%nst+1:im*st%nst) = Pd(1:st%nst,1:st%nst)
      if(im<Fdim) HF(im*st%nst+1:(im+1)*st%nst , (im-1)*st%nst+1:im*st%nst) = P(1:st%nst,1:st%nst)
    end do
    hm%F%floquet_apply = .true.
    hm%d%dim = st%d%dim

    call lalg_eigensolve(st%nst*Fdim, HF, evalues_full)

    ! if this is the initial step filter the zero harmonics by downfolding residual
    if (start_) then
      do ist=1,st%nst
        ! SCF cycle that should find a state approximating the zero sector of the zero harmonic (full state)
        do ii=1,hm%F%cf_nsteps
          call continued_fraction(st%nst, H0, P, Pd,st%eigenval(ist,ik), hm%F%omega, Fdim/2,Heff)
          call lalg_eigensolve(st%nst, Heff, evalues0)
          do jst=1,st%nst
            evals_diff(ist,jst) = abs(st%eigenval(ist,ik) - evalues0(jst))
          end do
          pos = minloc(evals_diff(ist,:),dim=1)
          st%eigenval(ist,ik) = evalues0(pos)   
        end do ! end cycle 
        ! find state that has a zero sector that most resembles the downfolding result
        do jst=1,st%nst*Fdim
          overlap(jst)= abs(dot_product(HF((hm%F%order(2))*st%nst+1:(hm%F%order(2)+1)*st%nst,jst),Heff(1:st%nst,ist))) 
        end do
        pos = maxloc(overlap,dim=1)
        evecs_temp(1:st%nst*Fdim,ist) = HF(1:st%nst*Fdim,pos)
        st%eigenval(ist,ik) = evalues_full(pos)

      end do

      call sort(st%eigenval(1:st%nst,ik), idx)
      
      do jst=1,st%nst
        evecs(1:st%nst*Fdim,jst) =  evecs_temp(1:st%nst*Fdim,idx(jst))
      end do
    else
!      ! for all regular iterations we filter the FBZ by comparison with input
!      do ist=1,st%nst*Fdim
!        do jst=1,st%nst
!          evals_diff(ist,jst) = abs(evalues_full(ist) - st%eigenval(jst,ik))
!        end do
!      end do
!
!      ! keep states with those eigenvalues that have the smallest difference to input
!      do jst=1,st%nst
!        pos = minloc(evals_diff(:,jst),dim=1)
!        st%eigenval(jst,ik) = evalues_full(pos)
!        evecs(1:st%nst*Fdim,jst) = HF(1:st%nst*Fdim,pos) 
!      end do
     norms = M_ONE
     do ist=1,st%nst*Fdim 
       evecs_reshape = reshape(HF(:,ist),(/st%nst,Fdim/))
       norms(ist) = sqrt(dot_product(evecs_reshape(1:st%nst,Fdim/2+1),evecs_reshape(1:st%nst,Fdim/2+1)))
     end do

     call sort(norms, idx)

     do jst=1,st%nst 
        pos = idx(st%nst*Fdim+1-jst)
        st%eigenval(jst,ik) = evalues_full(pos) 
        evecs_temp(1:st%nst*Fdim,jst) = HF(1:st%nst*Fdim,pos)
     end do

     call sort(st%eigenval(1:st%nst,ik), idx)

      do jst=1,st%nst
        evecs(1:st%nst*Fdim,jst) =  evecs_temp(1:st%nst*Fdim,idx(jst))
      end do

   end if


    ! rotate states by hand
    SAFE_ALLOCATE(    state(1:der%mesh%np,1:st%d%dim))
    SAFE_ALLOCATE(rot_state(1:st%nst,1:der%mesh%np,1:st%d%dim))
    rot_state = M_z0
    do ist=1,st%nst
       evecs_reshape = reshape(evecs(:,ist),(/st%nst,Fdim/))
       do jst=1,st%nst
          call states_get_state(st,der%mesh, jst, ik, state)
          do im=hm%F%order(1),hm%F%order(2)
             imm = im - hm%F%order(1) + 1
            do in=hm%F%order(1),hm%F%order(2)
               inn = in - hm%F%order(1) + 1
!               if(im-in > Fdim .or. im-in < 1) cycle
              if(im-in-hm%F%order(1) > Fdim .or. im-in-hm%F%order(1) < 0) cycle
              do idim=1,hm%F%spindim
                rot_state(ist,1:der%mesh%np,(imm-1)*hm%F%spindim+idim) =  rot_state(ist,1:der%mesh%np,(imm-1)*hm%F%spindim+idim)+ & 
                                                                   evecs_reshape(jst,imm)*state(1:der%mesh%np,(im-in-hm%F%order(1)-1)*hm%F%spindim+idim)
              end do
            end do
          end do

       end do
    end do
  
    do ist=1,st%nst
       state(1:der%mesh%np,1:st%d%dim) = rot_state(ist,1:der%mesh%np,1:st%d%dim)
       state(1:der%mesh%np,1:st%d%dim) = state(1:der%mesh%np,1:st%d%dim)/zmf_nrm2(der%mesh,st%d%dim,state(1:der%mesh%np,1:st%d%dim))
!       print *,ist, zmf_nrm2(der%mesh,st%d%dim,state(1:der%mesh%np,1:st%d%dim))
      call states_set_state(st,der%mesh, ist, ik, state(1:der%mesh%np,1:st%d%dim))
    end do

    SAFE_DEALLOCATE_A(state)
    SAFE_DEALLOCATE_A(rot_state)


    
    SAFE_DEALLOCATE_A(H0)
    SAFE_DEALLOCATE_A( P)
    SAFE_DEALLOCATE_A(Pd)
    SAFE_DEALLOCATE_A(Heff)
    SAFE_DEALLOCATE_A(evecs)
    SAFE_DEALLOCATE_A(evalues_full)
    SAFE_DEALLOCATE_A(evalues0)
    SAFE_DEALLOCATE_A(evals_diff)
    SAFE_DEALLOCATE_A(HF)
    SAFE_DEALLOCATE_A(one)
    SAFE_DEALLOCATE_A(overlap)
    SAFE_DEALLOCATE_A(norms)
    SAFE_DEALLOCATE_A(idx)
    SAFE_DEALLOCATE_A(evecs_temp)

    POP_SUB(floquet_FBZ_subspace_diag)

    contains

      subroutine continued_fraction(nn,H0, P, Pd, Q, omega, order,Heff)
        integer :: nn
        CMPLX  :: H0(nn,nn), P(nn,nn), Pd(nn,nn), Heff(nn,nn)
        FLOAT  ::  Q, omega
        integer :: order, ii

        FLOAT   :: one(nn,nn)
        CMPLX  :: Heff_plus(nn,nn)
        CMPLX  :: Heff_minus(nn,nn)

        one = M_ZERO
        do ii=1,nn
           one(ii,ii) = M_ONE
        end do

        Heff_minus(:,:) = M_ZERO
        Heff_plus(:,:)  = M_ZERO

        do ii=order,1,-1
           Heff_minus(:,:) = (Q-ii*omega)*one(:,:)-H0(:,:)-Heff_minus(:,:)
           Heff_plus(:,:)  = (Q+ii*omega)*one(:,:)-H0(:,:)-Heff_plus(:,:)

           call lalg_sym_inverter('u',nn, Heff_minus)
           call lalg_sym_inverter('u',nn, Heff_plus)
 
           Heff_minus = matmul(matmul(Pd,Heff_minus),P)
           Heff_plus = matmul(matmul(P,Heff_plus),Pd)
        end do

        Heff = H0+Heff_minus+Heff_plus

      end subroutine continued_fraction
  
  end subroutine zfloquet_FBZ_subspace_diag





#include "undef.F90"
#include "real.F90"
#include "subspace_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "subspace_inc.F90"

end module subspace_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
