!! Copyright (C) 2016 N. Tancogne-Dejean 
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
!! $Id$


subroutine X(lda_u_apply)(this, d, ik, psib, hpsib, has_phase)
  type(lda_u_t),      intent(in) :: this
  integer,            intent(in) :: ik
  type(batch_t),      intent(in) :: psib
  type(batch_t),   intent(inout) :: hpsib
  type(states_dim_t), intent(in) :: d
  logical,            intent(in) :: has_phase !True if the wavefunction has an associated phase

  integer :: ibatch, ios, imp, im, ispin, bind1, bind2
  R_TYPE, allocatable :: dot(:,:,:), reduced(:,:)
  type(orbitalset_t), pointer  :: os
  type(profile_t), save :: prof

  call profiling_in(prof, "DFTU_APPLY")

  PUSH_SUB(lda_u_apply)

  SAFE_ALLOCATE(reduced(1:this%max_np,1:psib%nst_linear))
  SAFE_ALLOCATE(dot(1:d%dim,1:this%maxnorbs, 1:psib%nst))

  ispin = states_dim_get_spin_index(d, ik)

  ! We have to compute 
  ! hpsi> += sum_m |phi m> sum_mp Vmmp <phi mp | psi >
  !
  ! We first compute <phi m | psi> for all orbitals of the atom
  !
  do ios = 1, this%norbsets
    os => this%orbsets(ios)
    call X(orbitalset_get_coeff_batch)(os, d%dim, psib, ik, has_phase, dot(1:d%dim,1:os%norbs,1:psib%nst))

    !
    reduced(:,:) = R_TOTYPE(M_ZERO) 
    !
    do ibatch = 1, psib%nst
      bind1 = batch_ist_idim_to_linear(psib, (/ibatch, 1/))
      bind2 = batch_ist_idim_to_linear(psib, (/ibatch, 2/))
      do im = 1,this%orbsets(ios)%norbs
        ! sum_mp Vmmp <phi mp | psi >
        do imp = 1, this%orbsets(ios)%norbs

          !Note here that V_{mmp} =U/2(delta_{mmp}-2n_{mpm})
          if(d%ispin /= SPINORS) then
            reduced(im,ibatch) = reduced(im,ibatch) + this%X(V)(im,imp,ispin,ios)*dot(1,imp,ibatch)
          else
            reduced(im,bind1) = reduced(im,bind1) + this%X(V)(im,imp,1,ios)*dot(1,imp,ibatch)
            reduced(im,bind1) = reduced(im,bind1) + this%X(V)(im,imp,3,ios)*dot(2,imp,ibatch)
            reduced(im,bind2) = reduced(im,bind2) + this%X(V)(im,imp,4,ios)*dot(1,imp,ibatch)
            reduced(im,bind2) = reduced(im,bind2) + this%X(V)(im,imp,2,ios)*dot(2,imp,ibatch)
          end if
        end do
      end do      
    end do !ibatch
 
    !We add the orbitals properly weighted to hpsi
    call X(orbitalset_add_to_batch)(this%orbsets(ios), d%dim, hpsib, ik, has_phase, reduced)
  end do
 
  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(reduced)

  POP_SUB(lda_u_apply)
  call profiling_out(prof)

end subroutine X(lda_u_apply)

! ---------------------------------------------------------
!> This routine computes the values of the occupation matrices
! ---------------------------------------------------------
subroutine X(update_occ_matrices)(this, mesh, st, lda_u_energy, phase)
  type(lda_u_t), intent(inout)         :: this
  type(mesh_t),     intent(in)         :: mesh
  type(states_t),  intent(in)          :: st
  FLOAT, intent(inout)                 :: lda_u_energy
  CMPLX, pointer, optional             :: phase(:,:) 

  integer :: ios, im, ik, ist, ispin, norbs, ip, idim
  integer :: ios2, im2
  R_TYPE, allocatable :: psi(:,:) 
  R_TYPE, allocatable :: dot(:,:,:)
  FLOAT   :: weight
  R_TYPE  :: renorm_weight
  type(orbitalset_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, "DFTU_OCC_MATRICES")
  
  PUSH_SUB(update_occ_matrices)

  this%X(n)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(dot(1:st%d%dim,1:this%maxnorbs,1:this%norbsets))

  if(this%useACBN0) then
    this%X(n_alt)(1:this%maxnorbs,1:this%maxnorbs,1:st%d%nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
  end if

  !TODO: use symmetries of the occupation matrices
  do ik = st%d%kpt%start, st%d%kpt%end
    ispin =  states_dim_get_spin_index(st%d,ik)

    !We recompute the overlap matrices here
    if(this%useACBN0) then
  !    print *, '=========== ik ', ik,  '========='
      call X(build_overlap_matrices)(this, ik, present(phase))
  !    print *, '================================='
    end if

    do ist = st%st_start, st%st_end

      weight = st%d%kweights(ik)*st%occ(ist, ik)     
      call states_get_state(st, mesh, ist, ik, psi )

      if(present(phase)) then
#ifdef R_TCOMPLEX
        ! Apply the phase that contains both the k-point and vector-potential terms.
        do idim = 1, st%d%dim 
          !$omp parallel do
          do ip = 1, mesh%np
            psi(ip, idim) = phase(ip, ik)*psi(ip, idim)
          end do
          !$omp end parallel do
        end do
#endif
      end if

      do ios = 1, this%norbsets
        os => this%orbsets(ios)
        !We first compute the matrix elemets <\psi | orb_m>
        !taking into account phase correction if needed 
        call X(orbitalset_get_coefficients)(os, st%d%dim, psi, ik, present(phase), &
                            dot(1:st%d%dim,1:os%norbs,ios))
      end do !ios


      !We compute the on-site occupation of the site, if needed 
      if(this%useACBN0) then
        this%X(renorm_occ)(:,:,:,ist,ik) = R_TOTYPE(M_ZERO)
        do ios = 1, this%norbsets
          os => this%orbsets(ios)
          norbs = this%orbsets(ios)%norbs
          do im = 1, norbs 
            !This sum runs over all the orbitals of the same type as im
            do ios2 = 1, this%norbsets
              do im2 = 1, this%orbsets(ios2)%norbs
                do idim = 1, st%d%dim
                  this%X(renorm_occ)(species_index(os%spec),os%nn,os%ll,ist,ik) = &
                   this%X(renorm_occ)(species_index(os%spec),os%nn,os%ll,ist,ik) &
                     + R_CONJ(dot(idim,im,ios))*dot(idim,im2,ios2)*os%X(S)(im,im2,ios2)
                end do
              end do
            end do
          end do
        end do
      end if
     

      if(st%d%ispin /= SPINORS) then !Collinear case

        !We can compute the (renormalized) occupation matrices
        do ios = 1, this%norbsets
          os => this%orbsets(ios)
          norbs = this%orbsets(ios)%norbs
          do im = 1, norbs
              this%X(n)(1:norbs,im,ispin,ios) = this%X(n)(1:norbs,im,ispin,ios) &
                                           + weight*dot(1,1:norbs,ios)*R_CONJ(dot(1,im,ios))
              !We compute the renomalized occupation matrices
              if(this%useACBN0) then
                renorm_weight = this%X(renorm_occ)(species_index(os%spec),os%nn,os%ll,ist,ik)*weight
                this%X(n_alt)(1:norbs,im,ispin,ios) = this%X(n_alt)(1:norbs,im,ispin,ios) &
                                           + renorm_weight*dot(1,1:norbs,ios)*R_CONJ(dot(1,im,ios))
              end if 
          end do !im
        end do !ios

      else !Noncollinear case

        !We can compute the (renormalized) occupation matrices
        do ios = 1, this%norbsets
          os => this%orbsets(ios)
          norbs = this%orbsets(ios)%norbs
          do im = 1, norbs
            this%X(n)(1:norbs,im,1,ios) = this%X(n)(1:norbs,im,1,ios) &
                                      + weight*dot(1,1:norbs,ios)*R_CONJ(dot(1,im,ios))
            this%X(n)(1:norbs,im,2,ios) = this%X(n)(1:norbs,im,2,ios) &
                                      + weight*dot(2,1:norbs,ios)*R_CONJ(dot(2,im,ios))
            this%X(n)(1:norbs,im,3,ios) = this%X(n)(1:norbs,im,3,ios) &
                                      + weight*dot(1,1:norbs,ios)*R_CONJ(dot(2,im,ios))
            this%X(n)(1:norbs,im,4,ios) = this%X(n)(1:norbs,im,4,ios) &
                                      + weight*dot(2,1:norbs,ios)*R_CONJ(dot(1,im,ios))
            !We compute the renomalized occupation matrices
            if(this%useACBN0) then
              renorm_weight = this%X(renorm_occ)(species_index(os%spec),os%nn,os%ll,ist,ik)*weight
              this%X(n_alt)(1:norbs,im,1,ios) = this%X(n_alt)(1:norbs,im,1,ios) &
                                         + renorm_weight*dot(1,1:norbs,ios)*R_CONJ(dot(1,im,ios))
              this%X(n_alt)(1:norbs,im,2,ios) = this%X(n_alt)(1:norbs,im,2,ios) &
                                         + renorm_weight*dot(2,1:norbs,ios)*R_CONJ(dot(2,im,ios))
              this%X(n_alt)(1:norbs,im,3,ios) = this%X(n_alt)(1:norbs,im,3,ios) &
                                         + renorm_weight*dot(1,1:norbs,ios)*R_CONJ(dot(2,im,ios))
              this%X(n_alt)(1:norbs,im,4,ios) = this%X(n_alt)(1:norbs,im,4,ios) &
                                         + renorm_weight*dot(2,1:norbs,ios)*R_CONJ(dot(1,im,ios))
            end if
          end do !im
        end do !ios 

      end if !Spinors
    end do
  end do

  
  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(psi)

#if defined(HAVE_MPI)        
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp%comm, this%X(n))
    if(this%useACBN0) &
      call comm_allreduce(st%st_kpt_mpi_grp%comm, this%X(n_alt))
  end if
#endif      

  if(this%useACBN0 .and. .not.this%freeze_u) then
    if(this%nspins > 1 ) then
      do ios = 1, this%norbsets
        if(this%orbsets(ios)%ndim  == 1) then
          call X(compute_ACBNO_U)(this, ios)
        else
          call compute_ACBNO_U_noncollinear(this, ios)
        end if
      end do
    else
      call X(compute_ACBNO_U_restricted)(this)
    end if
  end if


  call X(compute_dftu_energy)(this, lda_u_energy, st)
  call X(lda_u_update_potential)(this,st)

  POP_SUB(update_occ_matrices)
  call profiling_out(prof)
end subroutine X(update_occ_matrices)

! ---------------------------------------------------------
!> This routine computes the value of the double counting term in the LDA+U energy
! ---------------------------------------------------------
subroutine X(compute_dftu_energy)(this, energy, st)
  type(lda_u_t), intent(inout)    :: this
  FLOAT, intent(inout)            :: energy
  type(states_t), intent(in)      :: st 

  integer :: ios, imp, im, ispin

  PUSH_SUB(compute_dftu_energy)

  energy = M_ZERO

  do ios = 1, this%norbsets
    do ispin = 1, this%nspins
      !TODO: These are matrix operations, that could be optimized
      do im = 1, this%orbsets(ios)%norbs
        do imp = 1, this%orbsets(ios)%norbs
          energy = energy - CNST(0.5)*this%orbsets(ios)%Ueff*abs(this%X(n)(im,imp,ispin,ios))**2/st%smear%el_per_state
        end do
        if(ispin <= this%spin_channels) &
          energy = energy + CNST(0.5)*this%orbsets(ios)%Ueff*real(this%X(n)(im,im,ispin,ios))
      end do
    end do
  end do

  POP_SUB(compute_dftu_energy)
end subroutine X(compute_dftu_energy)


! ---------------------------------------------------------
!> This routine computes the potential that, once multiplied
!> by the projector Pmm' and summed over m and m' for all the atoms
!> gives the full Hubbard potential
! ---------------------------------------------------------
subroutine X(lda_u_update_potential)(this, st)
  type(lda_u_t), intent(inout)    :: this
  type(states_t), intent(in)      :: st

  integer :: ios, im, ispin, norbs
  type(profile_t), save :: prof

  call profiling_in(prof, "DFTU_POTENTIAL")

  PUSH_SUB(lda_u_update_potential)

  this%X(V) = M_ZERO

  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    do ispin = 1, this%nspins
      do im = 1, norbs
        this%X(V)(1:norbs,im,ispin,ios) = - this%orbsets(ios)%Ueff*this%X(n)(1:norbs,im,ispin,ios)/st%smear%el_per_state
        ! Only the diagonal part in spin space (for spinors)
        if(ispin <= this%spin_channels) &
          this%X(V)(im,im,ispin,ios) = this%X(V)(im,im,ispin,ios) + CNST(0.5)*this%orbsets(ios)%Ueff
      end do
    end do
  end do

  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp%comm, this%X(V))
  end if

  POP_SUB(lda_u_update_potential)
  call profiling_out(prof)
end subroutine X(lda_u_update_potential)

! ---------------------------------------------------------
!> This routine computes the effective U following the expression 
!> given in Agapito et al., Phys. Rev. X 5, 011006 (2015)
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U)(this, ios)
  type(lda_u_t), intent(inout)    :: this
  integer,       intent(in)       :: ios
  
  integer :: im, imp, impp, imppp, ispin1, ispin2, norbs
  FLOAT   :: numU, numJ, denomU, denomJ, tmpU, tmpJ

  PUSH_SUB(compute_ACBNO_U)

  ASSERT(this%orbsets(ios)%ndim == 1)

  norbs = this%orbsets(ios)%norbs
  numU = M_ZERO
  numJ = M_ZERO
  denomU = M_ZERO
  denomJ = M_ZERO

  if(norbs > 1) then

    do im = 1, norbs
    do imp = 1,norbs
      do impp = 1, norbs
      do imppp = 1, norbs
        ! We first compute the terms
        ! sum_{alpha,beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
        ! sum_{alpha} P^alpha_{mmp}P^alpha_{mpp,mppp}
        tmpU = M_ZERO
        tmpJ = M_ZERO
    
        do ispin1 = 1, this%spin_channels
          do ispin2 = 1, this%spin_channels
            tmpU = tmpU + real(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin2,ios))
          end do
          tmpJ = tmpJ + real(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin1,ios))
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          tmpJ = tmpJ + real(this%X(n_alt)(im,imp,3,ios)*this%X(n_alt)(impp,imppp,4,ios) &
                            +this%X(n_alt)(im,imp,4,ios)*this%X(n_alt)(impp,imppp,3,ios))
        end if
        ! These are the numerator of the ACBN0 U and J
        numU = numU + tmpU*this%coulomb(im,imp,impp,imppp,ios)
        numJ = numJ + tmpJ*this%coulomb(im,imppp,impp,imp,ios)
      end do
      end do

      ! We compute the term
      ! sum_{alpha} sum_{m,mp/=m} N^alpha_{m}N^alpha_{mp}
      tmpJ = M_ZERO
      tmpU = M_ZERO
      if(imp/=im) then
        do ispin1 = 1, this%spin_channels
          tmpJ = tmpJ + real(this%X(n)(im,im,ispin1,ios)*this%X(n)(imp,imp,ispin1,ios))
          tmpU = tmpU + real(this%X(n)(im,im,ispin1,ios)*this%X(n)(imp,imp,ispin1,ios))
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          tmpJ = tmpJ + real(this%X(n)(im,im,3,ios)*this%X(n)(imp,imp,4,ios) &
                            +this%X(n)(im,im,4,ios)*this%X(n)(imp,imp,3,ios))
        end if
      end if
      denomJ = denomJ + tmpJ

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels
          if(ispin1 /= ispin2) then
            tmpU = tmpU + real(this%X(n)(im,im,ispin1,ios)*this%X(n)(imp,imp,ispin2,ios))
          end if
        end do
      end do

      if(this%nspins>this%spin_channels) then !Spinors
        if(im == imp) then
          tmpU = tmpU - real(this%X(n)(im,im,3,ios)*this%X(n)(im,im,4,ios) &
                            +this%X(n)(im,im,4,ios)*this%X(n)(im,im,3,ios))
        end if
      end if 

      denomU = denomU + tmpU

    end do
    end do
    this%orbsets(ios)%Ueff = numU/denomU - numJ/denomJ
    this%orbsets(ios)%Ubar = numU/denomU
    this%orbsets(ios)%Jbar = numJ/denomJ
 
  else !In the case of s orbitals, the expression is different
    ! sum_{alpha/=beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
    ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
    numU = M_ZERO
    denomU = M_ZERO
    do ispin1 = 1, this%spin_channels
      do ispin2 = 1, this%spin_channels
        if(ispin1 /= ispin2) then
          numU = numU + real(this%X(n_alt)(1,1,ispin1,ios)*this%X(n_alt)(1,1,ispin2,ios)) 
          denomU = denomU + real(this%X(n)(1,1,ispin1,ios)*this%X(n)(1,1,ispin2,ios))
        end if
      end do
    end do

    if(this%nspins>this%spin_channels) then !Spinors
      denomU = denomU + real(this%X(n)(1,1,3,ios)*this%X(n)(1,1,4,ios) & 
                            +this%X(n)(1,1,4,ios)*this%X(n)(1,1,3,ios))
    end if

    ! We have to be careful in the case of hydrogen atom for instance 
    if(abs(denomU)> CNST(1.0e-3)) then
      this%orbsets(ios)%Ubar = (numU/denomU)*this%coulomb(1,1,1,1,ios)
    else
      if( abs(numU-denomU) < CNST(1.0e-3)) then
        this%orbsets(ios)%Ubar = this%coulomb(1,1,1,1,ios)
      else
        this%orbsets(ios)%Ubar = (numU/denomU)
        write(message(1),'(a,a)')' Small denominator value for the s orbital ', this%orbsets(ios)%Ubar
        write(message(2),'(a,a)')' to be multiplied by ',  this%coulomb(1,1,1,1,ios)
        call messages_warning(2) 
        this%orbsets(ios)%Ubar = this%orbsets(ios)%Ubar*this%coulomb(1,1,1,1,ios)
      end if
    end if
    
    this%orbsets(ios)%Jbar = 0
    this%orbsets(ios)%Ueff = this%orbsets(ios)%Ubar
  end if

  POP_SUB(compute_ACBNO_U)  
end subroutine X(compute_ACBNO_U)


! ---------------------------------------------------------
!> This routine computes the effective Uin the spin-unpolarised case
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U_restricted)(this)
  type(lda_u_t), intent(inout)    :: this
  
  integer :: ios, im, imp, impp, imppp, norbs
  FLOAT   :: numU, numJ, denomU, denomJ

  PUSH_SUB(compute_ACBNO_U_restricted)

  do ios = 1, this%norbsets
    norbs = this%orbsets(ios)%norbs
    numU = M_ZERO
    numJ = M_ZERO
    denomU = M_ZERO
    denomJ = M_ZERO

    if(norbs > 1) then

      do im = 1, norbs
      do imp = 1,norbs
        do impp = 1, norbs
        do imppp = 1, norbs
          numU = numU + real(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios)) &
                             *this%coulomb(im,imp,impp,imppp,ios)
          numJ = numJ + real(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios)) &
                             *this%coulomb(im,imppp,impp,imp,ios)
        end do
        end do
        ! We compute the term
        ! sum_{m,mp/=m} N_{m}N_{mp}
        if(imp/=im) then
          denomJ = denomJ + real(this%X(n)(im,im,1,ios)*this%X(n)(imp,imp,1,ios))
          denomU = denomU + real(this%X(n)(im,im,1,ios)*this%X(n)(imp,imp,1,ios))
        end if
        ! We compute the term
        ! sum_{m,mp} N_{m}N_{mp}
        denomU = denomU + real(this%X(n)(im,im,1,ios)*this%X(n)(imp,imp,1,ios))
      end do
      end do
      this%orbsets(ios)%Ueff = M_TWO*numU/denomU - numJ/denomJ
      this%orbsets(ios)%Ubar = M_TWO*numU/denomU
      this%orbsets(ios)%Jbar = numJ/denomJ
 
    else !In the case of s orbitals, the expression is different
      ! P_{mmp}P_{mpp,mppp}(m,mp|mpp,mppp)  
      numU = real(this%X(n_alt)(1,1,1,ios)*this%X(n_alt)(1,1,1,ios))*this%coulomb(1,1,1,1,ios)

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      denomU = real(this%X(n)(1,1,1,ios)*this%X(n)(1,1,1,ios))

      this%orbsets(ios)%Ueff = 2*numU/denomU
      this%orbsets(ios)%Ubar = 2*numU/denomU
      this%orbsets(ios)%Jbar = 0
    end if
  end do

  POP_SUB(compute_ACBNO_U_restricted)  
end subroutine X(compute_ACBNO_U_restricted)

! ---------------------------------------------------------
! TODO: Merge this with the two_body routine in system/output_me_inc.F90
subroutine X(compute_coulomb_integrals) (this, mesh, der)
  type(lda_u_t),   intent(inout)  :: this
  type(mesh_t),       intent(in)  :: mesh
  type(derivatives_t), intent(in) :: der

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: norbs, np_sphere, ios, ip
  integer :: idone, ntodo
  FLOAT, allocatable :: tmp(:), vv(:), nn(:)
  type(orbitalset_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, "DFTU_COULOMB_INTEGRALS")

  PUSH_SUB(X(compute_coulomb_integrals))

  ASSERT(.not. mesh%parallel_in_domains)
  
  SAFE_ALLOCATE(nn(1:this%max_np))
  SAFE_ALLOCATE(vv(1:this%max_np))
  SAFE_ALLOCATE(tmp(1:this%max_np))

  SAFE_ALLOCATE(this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs, 1:this%norbsets))
  this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%norbsets) = M_ZERO

  !Lets counts the number of orbital to treat, to display a progress bar
  ntodo = 0
  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    ntodo= ntodo + ((norbs+1)*norbs/2)*((norbs+1)*norbs/2+1)/2
  end do 
  idone = 0
  if(mpi_world%rank == 0) call loct_progress_bar(-1, ntodo)

  do ios = this%orbs_dist%start, this%orbs_dist%end
    os => this%orbsets(ios)
    norbs = os%norbs
    np_sphere = os%sphere%np  

    call poisson_init_sm(os%poisson, psolver, der, os%sphere) 
 
    ijst=0
    do ist = 1, norbs
      
      do jst = 1, norbs
        if(jst > ist) cycle
        ijst=ijst+1

        !$omp parallel do
        do ip=1,np_sphere
          nn(ip)  = real(os%X(orb)(ip,1,ist))*real(os%X(orb)(ip,1,jst))
        end do
        !$omp end parallel do    

        !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
        call dpoisson_solve_sm(os%poisson, os%sphere, vv(1:np_sphere), nn(1:np_sphere))

        klst=0
        do kst = 1, norbs
          do lst = 1, norbs
            if(lst > kst) cycle
            klst=klst+1
            if(klst > ijst) cycle

            !$omp parallel do
            do ip=1,np_sphere
              tmp(ip) = vv(ip)*real(os%X(orb)(ip,1,lst))*real(os%X(orb)(ip,1,kst))
            end do
            !$omp end parallel do

            this%coulomb(ist,jst,kst,lst,ios) = dsm_integrate(mesh, os%sphere, tmp(1:np_sphere))

            if(abs(this%coulomb(ist,jst,kst,lst,ios))<CNST(1.0e-12)) then
              this%coulomb(ist,jst,kst,lst,ios) = M_ZERO
            else
              this%coulomb(kst,lst,ist,jst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(jst,ist,lst,kst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(lst,kst,jst,ist,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(jst,ist,kst,lst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(lst,kst,ist,jst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(ist,jst,lst,kst,ios) = this%coulomb(ist,jst,kst,lst,ios)
              this%coulomb(kst,lst,jst,ist,ios) = this%coulomb(ist,jst,kst,lst,ios)
            end if
              
          !Update the progress bar
          idone = idone + 1
          if(mpi_world%rank == 0) call loct_progress_bar(idone, ntodo)
          end do !lst
        end do !kst
      end do !jst
    end do !ist
    call poisson_end(os%poisson)
  end do !iorb

  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp%comm, this%coulomb)
  end if
 
  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(compute_coulomb_integrals))
  call profiling_out(prof)
end subroutine X(compute_coulomb_integrals)

 ! ---------------------------------------------------------
 !> This routine computes [r,V_lda+u].
 ! ---------------------------------------------------------
 subroutine X(lda_u_commute_r)(this, mesh, d, ik, ist, psi, gpsi, has_phase)
   type(lda_u_t),      intent(in) :: this
   type(mesh_t),    intent(in)    :: mesh
   type(states_dim_t), intent(in) :: d
   R_TYPE,          intent(in)    :: psi(:,:)
   integer,         intent(in)    :: ik, ist
   R_TYPE,          intent(inout) :: gpsi(:, :, :)
   logical,            intent(in) :: has_phase !True if the wavefunction has an associated phase 

   integer :: ios, idim, idir, im, imp, is, ispin
   integer :: idim_orb
   R_TYPE, allocatable :: dot(:,:)
   R_TYPE, allocatable :: epsi(:,:)
   type(orbitalset_t), pointer  :: os 
   R_TYPE, allocatable  :: reduced(:,:)
   type(profile_t), save :: prof

   call profiling_in(prof, "DFTU_COMMUTE_R")

   PUSH_SUB(lda_u_commute_r)

   if(simul_box_is_periodic(mesh%sb) .and. .not. this%submeshforperiodic) then
     SAFE_ALLOCATE(epsi(1:mesh%np,1:d%dim))
   else
     SAFE_ALLOCATE(epsi(1:this%max_np,1:d%dim))
   end if
   SAFE_ALLOCATE(dot(1:d%dim,1:this%maxnorbs))
   SAFE_ALLOCATE(reduced(1:d%dim,1:this%maxnorbs))

   ispin = states_dim_get_spin_index(d, ik)

   do ios = 1, this%norbsets
      ! We have to compute 
      ! hpsi> += r sum_m |phi m> sum_mp Vmmp <phi mp | psi >
      !
      ! We first compute <phi m | psi> for all orbitals of the atom
      !
      os => this%orbsets(ios)
      ! 
      call X(orbitalset_get_coefficients)(os, d%dim, psi, ik, has_phase, dot)
      !
      reduced(:,:) = M_ZERO
      do im = 1, os%norbs
        ! sum_mp Vmmp <phi mp | psi >
        do imp = 1, os%norbs
          if(d%ispin /= SPINORS) then
            reduced(1,im) = reduced(1,im) + this%X(V)(im,imp,ispin,ios)*dot(1,imp)
          else
            reduced(1,im) = reduced(1,im) + this%X(V)(im,imp,1,ios)*dot(1,imp)
            reduced(1,im) = reduced(1,im) + this%X(V)(im,imp,3,ios)*dot(2,imp)
            reduced(2,im) = reduced(2,im) + this%X(V)(im,imp,4,ios)*dot(1,imp)
            reduced(2,im) = reduced(2,im) + this%X(V)(im,imp,2,ios)*dot(2,imp)
          end if
        end do
      end do

      do idir = 1, mesh%sb%dim
        do idim = 1, d%dim
          idim_orb = min(idim,os%ndim)
          do im = 1, os%norbs
          !In case of phase, we have to apply the conjugate of the phase here
          if(has_phase) then
#ifdef R_TCOMPLEX
            if(simul_box_is_periodic(mesh%sb) .and. .not. this%submeshforperiodic) then
              epsi(:,idim) = R_TOTYPE(M_ZERO)
              !$omp parallel do
              do is = 1, os%sphere%np
                  epsi(os%sphere%map(is),idim) = epsi(os%sphere%map(is),idim) &
                     + os%sphere%x(is,idir)*os%zorb(is,idim_orb,im)*os%phase(is,ik)
              end do
              !$omp end parallel do
              call lalg_axpy(mesh%np, reduced(idim,im), epsi(1:mesh%np,idim), &
                                  gpsi(1:mesh%np,idir,idim))
            else
              !$omp parallel do
              do is = 1, os%sphere%np
                epsi(is,idim) = os%sphere%x(is,idir)*os%eorb_submesh(is,idim_orb,im,ik)
              end do
              !$omp end parallel do
              call submesh_add_to_mesh(os%sphere, epsi(1:os%sphere%np,idim), &
                                  gpsi(1:mesh%np,idir,idim), reduced(idim,im))
            end if
#endif
            else
              !$omp parallel do
              do is = 1, os%sphere%np
                epsi(is,idim) = os%sphere%x(is,idir)*os%X(orb)(is,idim_orb,im)
              end do
              !$omp end parallel do
              call submesh_add_to_mesh(os%sphere, epsi(1:os%sphere%np,idim), &
                                    gpsi(1:mesh%np,idir,idim), reduced(idim,im))
            end if
          end do !im
        end do !idim
      end do !idir

     do idir = 1, mesh%sb%dim
       ! We have to compute 
       ! hpsi> -= sum_m |phi m> sum_mp Vmmp <phi mp| r | psi >
       !
       ! We first compute <phi m| r | psi> for all orbitals of the atom
       !
       !
       if(simul_box_is_periodic(mesh%sb) .and. .not. this%submeshforperiodic) then
         forall(is = 1:mesh%np)
           epsi(is,1) = mesh%x(is,idir)*psi(is,1)
         end forall
       else
         !$omp parallel do 
         do is = 1, os%sphere%np
           epsi(is,1) = os%sphere%x(is,idir)*psi(os%sphere%map(is),1)
         end do
         !$omp end parallel do
       end if     

       if(has_phase) then
#ifdef R_TCOMPLEX
         if(simul_box_is_periodic(mesh%sb) .and. .not. this%submeshforperiodic) then
           do im = 1, os%norbs
             do idim = 1, d%dim
               dot(idim,im) = X(mf_dotp)(mesh, os%eorb_mesh(1:mesh%np,idim,im,ik),&
                               epsi(1:mesh%np,idim), reduce = .false.)
             end do
           end do
         else
           do im = 1, os%norbs
             do idim = 1, d%dim
               dot(idim,im) = X(mf_dotp)(mesh, os%eorb_submesh(1:os%sphere%np,idim,im,ik),&
                               epsi(1:os%sphere%np,idim), reduce = .false., np = os%sphere%np)
             end do
           end do
         end if
#endif
       else
         do im = 1, os%norbs
           do idim = 1, d%dim
             dot(idim,im) = X(mf_dotp)(mesh, os%X(orb)(1:os%sphere%np,idim,im),&
                               epsi(1:os%sphere%np,idim), reduce = .false., np = os%sphere%np)
           end do
         end do
       end if
 
       do im = 1, os%norbs
         ! sum_mp Vmmp <phi mp|r| psi >
         do imp = 1, os%norbs
           if(d%ispin /= SPINORS) then
             reduced(1,im) = reduced(1,im) - this%X(V)(im,imp,ispin,ios)*dot(1,imp)
            else
              reduced(1,im) = reduced(1,im) - this%X(V)(im,imp,1,ios)*dot(1,imp)
              reduced(1,im) = reduced(1,im) - this%X(V)(im,imp,3,ios)*dot(2,imp)
              reduced(2,im) = reduced(2,im) - this%X(V)(im,imp,4,ios)*dot(1,imp)
              reduced(2,im) = reduced(2,im) - this%X(V)(im,imp,2,ios)*dot(2,imp)
            end if
         end do

       end do

       call X(orbitalset_add_to_psi)(os, d%dim, gpsi(1:mesh%np,idir,1:d%dim), ik, has_phase, &
                                         reduced(1:d%dim,1:os%norbs)) 
     end do !idir

   end do !ios

   SAFE_DEALLOCATE_A(epsi)
   SAFE_DEALLOCATE_A(dot)

   POP_SUB(lda_u_commute_r)
   call profiling_out(prof)
 end subroutine X(lda_u_commute_r)

 subroutine X(lda_u_force)(this, mesh, st, iq, ndim, psib, grad_psib, force, phase)
   type(lda_u_t),             intent(in)    :: this
   type(mesh_t),              intent(in)    :: mesh 
   type(states_t),            intent(in)    :: st
   integer,                   intent(in)    :: iq, ndim
   type(batch_t),             intent(in)    :: psib
   type(batch_t),             intent(in)    :: grad_psib(:)
   FLOAT,                     intent(inout) :: force(:, :)
   logical,                   intent(in)    :: phase

   integer :: ios, iatom, ibatch, ist, im, imp, ispin, idir
   type(orbitalset_t), pointer  :: os
   R_TYPE :: ff(1:ndim)
   R_TYPE, allocatable :: psi(:,:), gpsi(:,:)
   R_TYPE, allocatable :: dot(:,:), gdot(:,:,:), gradn(:,:,:,:)
   FLOAT :: weight

   if(.not. this%apply) return
   if(st%d%ispin == SPINORS) call messages_not_implemented("Hubbard forces with spinors")

   PUSH_SUB(X(lda_u_force))

   SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
   SAFE_ALLOCATE(gpsi(1:mesh%np, 1:st%d%dim))
   SAFE_ALLOCATE(dot(1:st%d%dim, 1:this%maxnorbs))
   SAFE_ALLOCATE(gdot(1:st%d%dim, 1:this%maxnorbs,1:ndim))
   SAFE_ALLOCATE(gradn(1:this%maxnorbs,1:this%maxnorbs,1:this%nspins,1:ndim))

   ispin = states_dim_get_spin_index(st%d, iq)

   do ios = 1, this%norbsets 
     os => this%orbsets(ios)
     iatom = os%iatom

     gradn(1:os%norbs,1:os%norbs,1:this%nspins,1:ndim) = M_ZERO

     do ibatch = 1, psib%nst_linear
       ist = batch_linear_to_ist(psib, ibatch) 
       weight = st%d%kweights(iq)*st%occ(ist, iq)

       call batch_get_state(psib, ibatch, mesh%np, psi)

       !We first compute the matrix elemets <\psi | orb_m>
       !taking into account phase correction if needed   
       ! 
       call X(orbitalset_get_coefficients)(os, st%d%dim, psi, iq, phase, dot)

       do idir = 1, ndim
         call batch_get_state(grad_psib(idir), ibatch, mesh%np, gpsi)     
         !We first compute the matrix elemets <\psi | orb_m>
         !taking into account phase correction if needed 
         ! 
         call X(orbitalset_get_coefficients)(os, st%d%dim, gpsi, iq, phase, gdot(1:st%d%dim,1:os%norbs,idir))

         do im = 1, os%norbs
           gradn(1:os%norbs,im,ispin,idir) = gradn(1:os%norbs,im,ispin,idir) &
                                        + weight*(R_CONJ(gdot(1,1:os%norbs,idir))*dot(1,im) &
                                                 +gdot(1,im,idir)*R_CONJ(dot(1,1:os%norbs)))
         end do
       end do !idir
       
     end do !ibatch

     ff(1:ndim) = M_ZERO
     do im = 1, os%norbs
       do imp = 1, os%norbs
        ff(1:ndim) = ff(1:ndim) - this%X(n)(im,imp,ispin,ios)/st%smear%el_per_state*gradn(im,imp,ispin,1:ndim)
       end do !imp
     ff(1:ndim) = ff(1:ndim) + CNST(0.5)*gradn(im, im, ispin,1:ndim)
     end do !im

     force(1:ndim, iatom) = force(1:ndim, iatom) - os%Ueff*real(ff(1:ndim))
   end do !ios

   SAFE_DEALLOCATE_A(psi)
   SAFE_DEALLOCATE_A(gpsi)
   SAFE_DEALLOCATE_A(dot)
   SAFE_DEALLOCATE_A(gdot)
   SAFE_DEALLOCATE_A(gradn)

   POP_SUB(X(lda_u_force))
 end subroutine X(lda_u_force)

! ---------------------------------------------------------
!> This routine is an interface for constructing the orbital basis.
! ---------------------------------------------------------
subroutine X(construct_orbital_basis)(this, geo, mesh, st)
  type(lda_u_t),             intent(inout)    :: this
  type(geometry_t), target,  intent(in)       :: geo
  type(mesh_t),              intent(in)       :: mesh
  type(states_t),            intent(in)       :: st 

  integer :: ia, iorb, norb, ntotorb, offset, ios, idim
  integer ::  hubbardl, ii, nn, ll, mm, work, work2, iorbset
  FLOAT   :: norm, hubbardj, radius, jj
  FLOAT, allocatable :: minradii(:)
  integer :: nSorbitals
  type(orbitalset_t), pointer :: os
  logical :: hasjdependence

  PUSH_SUB(X(construct_orbital_basis))

  write(message(1),'(a)')    'Building the LDA+U localized orbital basis.'
  call messages_info(1)

  !We first count the number of orbital sets we have to treat
  norb = 0
  if( .not. this%useAllOrbitals ) then
    do ia = 1, geo%natoms
      hubbardl = species_hubbard_l(geo%atom(ia)%species)
      hubbardj = species_hubbard_j(geo%atom(ia)%species)
      if( hubbardl .eq. -1 ) cycle
  
      !This is a dirty way to detect if the pseudopotential has j-dependent atomic wavefunctions
      hasjdependence = .false.
      call species_iwf_j(geo%atom(ia)%species, 1, jj)
      if(jj /= M_ZERO) hasjdependence = .true.

      if(hasjdependence .and. hubbardj == M_ZERO) then
        norb = norb+2
      else
        norb = norb+1
      end if

      if( hubbardj /= 0 .and. .not. hasjdependence) then
        write(message(1),'(a,i1,a)') 'Atom ', ia, ' has no j-dependent atomic wavefunctions.'
        write(message(2),'(a)') 'This is not compatible with the hubbard_j option.'
        call messages_fatal(2)
      end if
    end do
  else
    do ia = 1, geo%natoms
      work = 0
      nSorbitals = 0
      hubbardj = species_hubbard_j(geo%atom(ia)%species)
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm) 
        call species_iwf_j(geo%atom(ia)%species, iorb, jj)
        if(ll == 0) nSorbitals = nSorbitals + 1
        work = max(work, ii)

        if( hubbardj /= 0 .and. jj == 0 ) then
          write(message(1),'(a,i1,a)') 'Atom ', ia, ' has no j-dependent atomic wavefunction.'
          write(message(2),'(a)') 'This is not compatible with the hubbard_j option.'
          call messages_fatal(2)  
        end if
      end do
      if(this%skipSOrbitals) work = work-nSorbitals
      norb = norb + work
    end do
  end if


  write(message(1),'(a, i3, a)')    'Found ', norb, ' orbital sets.'
  call messages_info(1)

  this%norbsets = norb
  SAFE_ALLOCATE(this%orbsets(1:norb))
  do iorbset = 1, this%norbsets
    call orbitalset_nullify(this%orbsets(iorbset))
  end do

  if( this%useAllOrbitals .and. this%minimalAtomicSphere ) then
   SAFE_ALLOCATE(minradii(1:geo%natoms))
   call find_minimal_atomic_spheres(geo, mesh, minradii, this%truncation, this%orbitals_threshold)
  end if

  iorbset = 0
  do ia = 1, geo%natoms
    hubbardj = species_hubbard_j(geo%atom(ia)%species)
    !This is a dirty way to detect if the pseudopotential has j-dependent atomic wavefunctions
    hasjdependence = .false.
    call species_iwf_j(geo%atom(ia)%species, 1, jj)
    if(jj /= M_ZERO) hasjdependence = .true.
    if (debug%info .and. hasjdependence) then
      write(message(1),'(a,i3,a)')  'Debug: Atom ', ia, ' has j-dependent pseudo-wavefunctions.'
      call messages_info(1)
    end if 

    if(.not. this%useAllOrbitals) then
      hubbardl = species_hubbard_l(geo%atom(ia)%species)
      if( hubbardl .eq. -1 ) cycle
      !In this case we only have one orbital or we only want one
      if(.not. hasjdependence .or. hubbardj /= M_ZERO &
            .or. hubbardl == 0 ) then
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(geo%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. hubbardj == jj ) then
            norb = norb + 1
            os%ll = hubbardl
            os%nn = nn
            os%jj = jj 
            os%ii = ii
            os%radius = atomic_orbital_get_radius(geo, mesh, ia, iorb, 1, &
                          this%truncation, this%orbitals_threshold)
          end if
        end do
        if( hasjdependence ) then
          os%ndim = st%d%dim
          os%norbs = norb + int((os%jj - os%ll)*2)
        else
          os%ndim = 1
          os%norbs = norb
        end if
        os%Ueff = species_hubbard_u(geo%atom(ia)%species)
        os%submeshforperiodic = this%submeshforperiodic
        os%spec => geo%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, geo, mesh)
      else
        !j = l-1/2
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(geo%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. jj < ll ) then
            norb = norb + 1
            os%ll = hubbardl
            os%nn = nn
            os%jj = jj
            os%ii = ii
            os%radius = atomic_orbital_get_radius(geo, mesh, ia, iorb, 1, &
                        this%truncation, this%orbitals_threshold)
          end if
        end do
        os%ndim = st%d%dim
        os%norbs = norb-1
        os%Ueff = species_hubbard_u(geo%atom(ia)%species)
        os%submeshforperiodic = this%submeshforperiodic
        os%spec => geo%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, geo, mesh)

        !j = l+1/2
        iorbset = iorbset + 1
        os => this%orbsets(iorbset)
        norb = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(geo%atom(ia)%species, iorb, jj)
          if(ll .eq. hubbardl .and. jj > ll ) then
            norb = norb + 1
            os%ll = hubbardl
            os%nn = nn
            os%jj = jj
            os%ii = ii
            os%radius = atomic_orbital_get_radius(geo, mesh, ia, iorb, 1, &
                        this%truncation, this%orbitals_threshold)
          end if
        end do
        os%ndim = st%d%dim
        os%norbs = norb+1
        os%Ueff = species_hubbard_u(geo%atom(ia)%species)
        os%submeshforperiodic = this%submeshforperiodic
        os%spec => geo%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, geo, mesh)
      end if
    else !useAllOrbitals
      ASSERT(.not.hasjdependence)
      work = 0
      nSorbitals = 0
      do iorb = 1, species_niwfs(geo%atom(ia)%species)
        call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
        if(ll == 0) nSorbitals = nSorbitals + 1
        work = max(work, ii)          
      end do
      if(this%skipSOrbitals) work = work-nSorbitals
      offset = 0
      if(this%skipSOrbitals) offset = 1

      !We loop over the orbital sets of the atom ia
      do norb = 1, work
        os => this%orbsets(iorbset+norb)
        !We count the number of orbital for this orbital set
        work2 = 0
        do iorb = 1, species_niwfs(geo%atom(ia)%species)
          call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
          call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
          call species_iwf_j(geo%atom(ia)%species, iorb, jj )
          if(ii == norb+offset .and. hubbardj == jj) then
            work2 = work2 + 1
            os%ll = ll
            os%nn = nn
            os%jj = jj
            os%ii = ii
            if(this%minimalAtomicSphere) then
              radius = minradii(ia)
            else
              os%radius = atomic_orbital_get_radius(geo, mesh, ia, iorb, 1, &
                               this%truncation, this%orbitals_threshold)
            end if
          end if
        end do
        os%norbs = work2
        os%ndim = 1
        os%Ueff = species_hubbard_u(geo%atom(ia)%species)
        os%submeshforperiodic = this%submeshforperiodic
        os%spec => geo%atom(ia)%species
        os%iatom = ia
        call X(orbitalset_utils_getorbitals)(os, geo, mesh)
      end do !norb
      iorbset = iorbset + work
    end if
  end do

  ! We have to normalize the orbitals, 
  ! in case the orbitals that comes out of the pseudo are not properly normalised
  do iorbset = 1, this%norbsets
    os => this%orbsets(iorbset)
    do iorb = 1, os%norbs
      norm = M_ZERO
      do idim = 1, os%ndim
        norm = norm + X(sm_nrm2)(os%sphere, os%X(orb)(1:os%sphere%np,idim,iorb))**2
      end do
      norm = sqrt(norm)
      if(this%normalizeOrbitals) then
        do idim = 1, os%ndim
          os%X(orb)(1:os%sphere%np,idim,iorb) =  &
            os%X(orb)(1:os%sphere%np,idim,iorb)/norm
        end do
      else
        write(message(1),'(a,i3,a,i3,a,f8.5)') 'Info: Orbset ', iorbset, ' Orbital ', iorb, &
                           ' norm= ',  norm
        call messages_info(1)
      end if
    end do
  end do

  this%maxnorbs = 0
  this%max_np = 0
  do iorbset = 1, this%norbsets
    os => this%orbsets(iorbset)
    if( os%norbs > this%maxnorbs ) this%maxnorbs = os%norbs

    nullify(os%phase)
    ! In case of complex wavefunction, we allocate the array for the phase correction
  #ifdef R_TCOMPLEX
    SAFE_ALLOCATE(os%phase(1:os%sphere%np, st%d%kpt%start:st%d%kpt%end))
    os%phase(:,:) = M_ZERO
    if(simul_box_is_periodic(mesh%sb) .and. .not. this%submeshforperiodic) then 
      SAFE_ALLOCATE(os%eorb_mesh(1:mesh%np, 1:os%ndim, 1:os%norbs, st%d%kpt%start:st%d%kpt%end))
      os%eorb_mesh(:,:,:,:) = M_ZERO
    else
      SAFE_ALLOCATE(os%eorb_submesh(1:os%sphere%np, 1:os%ndim, 1:os%norbs, st%d%kpt%start:st%d%kpt%end))
      os%eorb_submesh(:,:,:,:) = M_ZERO
    end if
  #endif

    ! We need to know the maximum number of points in order to allocate a temporary array
    ! to apply the phase in lda_u_apply
    if(os%sphere%np > this%max_np) this%max_np = os%sphere%np

  end do  

  do ios = 1, this%norbsets
    if(this%orbsets(ios)%sphere%np == -1) then
       write(message(1),'(a,a4,i1,a1,a)')    'The orbital ',trim(species_label(this%orbsets(ios)%spec)), &
                      this%orbsets(ios)%nn, l_notation(this%orbsets(ios)%ll), ' has no grid point.'
       write(message(2),'(a)') 'Change the input file or use a pseudopotential that contains these orbitals.'
       call messages_fatal(1)
    end if

    write(message(1),'(a,i2,a,f8.5,a)')    'Orbital set ', ios, ' has a value of U of ',&
                         this%orbsets(ios)%Ueff   , ' Ha.'
    write(message(2),'(a,i2,a)')    'It contains ', this%orbsets(ios)%norbs, ' orbitals.'
    write(message(3),'(a,f8.5,a,i6,a)') 'The radius is ', this%orbsets(ios)%sphere%radius, &
                        ' Bohr,  with ', this%orbsets(ios)%sphere%np, ' grid points.'
     call messages_info(3)
  end do 

  if(this%UseACBN0) then
    do iorbset = 1, this%norbsets
      os => this%orbsets(iorbset)
      if(states_are_real(st)) then
        SAFE_ALLOCATE(os%dS(os%norbs,this%maxnorbs,this%norbsets))
      else
        SAFE_ALLOCATE(os%zS(os%norbs,this%maxnorbs,this%norbsets))
      end if
    end do
  end if

  if( this%useAllOrbitals .and. this%minimalAtomicSphere ) then
   SAFE_DEALLOCATE_A(minradii)
  end if

 
  POP_SUB(X(construct_orbital_basis))

end subroutine X(construct_orbital_basis)

subroutine X(build_overlap_matrices)(this, ik, has_phase)
  type(lda_u_t),    intent(inout)    :: this
  integer,          intent(in)       :: ik
  logical,          intent(in)       :: has_phase

  integer :: ios, ios2, im, im2, norbs, np
  type(orbitalset_t), pointer :: os, os2
!  R_TYPE, allocatable :: orb1(:,:), orb2(:)

  if(this%nspins > this%spin_channels .and. this%IncludeOverlap) &
      call messages_not_implemented("Overlap matrices with spinors")

  PUSH_SUB(X(build_overlap_matrices))

  np = this%orbsets(1)%sphere%mesh%np

  do ios = 1, this%norbsets
    os => this%orbsets(ios)
    norbs = this%orbsets(ios)%norbs
    os%X(S)(1:norbs,1:this%maxnorbs,1:this%norbsets) = R_TOTYPE(M_ZERO)

    !TODO: Use symmetry of the overlap matrices
    norbs = this%orbsets(ios)%norbs

    do im = 1, norbs

      do ios2 = 1, this%norbsets
        os2 => this%orbsets(ios2)

        if(ios2 == ios) then
          os%X(S)(im,im,ios2) = M_ONE
        else
          if(this%IncludeOverlap) then  
          do im2 = 1, os2%norbs
            if(has_phase) then
 #ifdef R_TCOMPLEX
              if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. this%submeshforperiodic) then
                os%X(S)(im,im2,ios2) = zmf_dotp(os%sphere%mesh, &
                                          os2%eorb_mesh(1:os2%sphere%mesh%np,1,im2,ik), &
                                          os%eorb_mesh(1:os%sphere%mesh%np,1,im,ik))
              else
                os%X(S)(im,im2,ios2) = zsubmesh_to_submesh_dotp(os2%sphere, &
                                         os2%eorb_submesh(1:os2%sphere%np,1,im2,ik), &
                                         os%sphere, os%eorb_submesh(1:os%sphere%np,1,im,ik))
              end if
 #endif
            else
              os%X(S)(im,im2,ios2) = X(submesh_to_submesh_dotp)(os2%sphere, os2%X(orb)(1:os2%sphere%np,1,im2), &
                   os%sphere, os%X(orb)(1:os%sphere%np,1,im))
            end if
         !   if(ios == ios2 ) then
         !     print *, ios, ios2, im, im2, os%X(S)(im,im2,ios2)
         !    end if
            end do ! im2
          end if
        end if
      end do !im
    end do !ios2
  end do !ios 

  POP_SUB(X(build_overlap_matrices))
end subroutine X(build_overlap_matrices)

 ! ---------------------------------------------------------
  subroutine X(lda_u_set_occupations)(this, occ)
    type(lda_u_t),  intent(inout) :: this
    R_TYPE,         intent(in)    :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs

    PUSH_SUB(X(lda_u_set_occupations))

    ind = 0
    do ios = 1, this%norbsets
      norbs = this%orbsets(ios)%norbs
      do ispin = 1, this%nspins
        do im = 1, norbs
          do imp = 1,norbs
            ind = ind + 1
            this%X(n)(im,imp,ispin,ios) = occ(ind)
          end do
        end do
      end do
    end do 

    if(this%useACBN0) then
      do ios = 1, this%norbsets
        norbs = this%orbsets(ios)%norbs
        do ispin = 1, this%nspins
          do im = 1, norbs
            do imp = 1,norbs
              ind = ind + 1
              this%X(n_alt)(im,imp,ispin,ios) = occ(ind)
            end do
          end do
        end do
      end do
    end if

    POP_SUB(X(lda_u_set_occupations))
  end subroutine X(lda_u_set_occupations)

   ! ---------------------------------------------------------
  subroutine X(lda_u_get_occupations)(this, occ)
    type(lda_u_t),  intent(in) :: this
    R_TYPE,      intent(inout) :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs

    PUSH_SUB(X(lda_u_get_occupations))

    ind = 0
    do ios = 1, this%norbsets
      norbs = this%orbsets(ios)%norbs
      do ispin = 1, this%nspins
        do im = 1, norbs
          do imp = 1,norbs
            ind = ind + 1
            occ(ind) = this%X(n)(im,imp,ispin,ios)
          end do
        end do
      end do
    end do

    if(this%useACBN0) then
      do ios = 1, this%norbsets
        norbs = this%orbsets(ios)%norbs
        do ispin = 1, this%nspins
          do im = 1, norbs
            do imp = 1,norbs
              ind = ind + 1
              occ(ind) = this%X(n_alt)(im,imp,ispin,ios)
            end do
          end do
        end do
      end do
    end if

    POP_SUB(X(lda_u_get_occupations))
  end subroutine X(lda_u_get_occupations)

  ! ---------------------------------------------------------
  subroutine X(lda_u_allocate)(this, st)
    type(lda_u_t),  intent(inout) :: this
    type(states_t), intent(in)    :: st

    integer :: maxorbs, nspin

    PUSH_SUB(X(lda_u_allocate))

    maxorbs = this%maxnorbs
    nspin = this%nspins

    SAFE_ALLOCATE(this%X(n)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
    this%X(n)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
    SAFE_ALLOCATE(this%X(V)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
    this%X(V)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%useACBN0) then
      SAFE_ALLOCATE(this%X(n_alt)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
      this%X(n_alt)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
      SAFE_ALLOCATE(this%X(renorm_occ)(this%nspecies,0:5,0:(MAX_L-1),st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end))
      this%X(renorm_occ)(this%nspecies,0:5,0:(MAX_L-1),st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end) = R_TOTYPE(M_ZERO)
    end if

    POP_SUB(X(lda_u_allocate))
  end subroutine X(lda_u_allocate)
