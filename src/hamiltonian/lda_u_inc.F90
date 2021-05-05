!! Copyright (C) 2016-2019 N. Tancogne-Dejean
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


subroutine X(lda_u_apply)(this, d, mesh, psib, hpsib)
  type(lda_u_t),           intent(in)    :: this
  type(mesh_t),            intent(in)    :: mesh
  type(wfs_elec_t),        intent(in)    :: psib
  type(wfs_elec_t),        intent(inout) :: hpsib
  type(states_elec_dim_t), intent(in)    :: d

  integer :: ibatch, ios, imp, im, ispin, bind1, bind2, inn, ios2
  R_TYPE, allocatable :: dot(:,:,:,:), reduced(:,:,:), psi(:,:)
  type(orbitalset_t), pointer  :: os
  type(profile_t), save :: prof
  integer :: el_per_state
  R_TYPE :: weight

  call profiling_in(prof, TOSTRING(X(DFTU_APPLY)))

  PUSH_SUB(lda_u_apply)

  SAFE_ALLOCATE(reduced(1:this%maxnorbs, 1:psib%nst_linear, 1:this%norbsets))
  SAFE_ALLOCATE(dot(1:d%dim,1:this%maxnorbs, 1:this%norbsets, 1:psib%nst))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:d%dim))

  ispin = d%get_spin_index(psib%ik)
  if(d%ispin == UNPOLARIZED) then
    el_per_state = 2
  else
    el_per_state = 1
  end if


  ! We have to compute 
  ! hpsi> += sum_m |phi m> sum_mp Vmmp <phi mp | psi >
  !
  ! We first compute <phi m | psi> for all orbitals of the atom
  !
  do ibatch = 1, psib%nst
    call batch_get_state(psib, ibatch, mesh%np, psi)
    do ios = 1, this%norbsets
      os => this%orbsets(ios)
      call X(orbitalset_get_coefficients)(os, d%dim, psi, psib%ik, psib%has_phase, &
                                                  dot(1:d%dim,1:os%norbs,ios,ibatch))
    end do
  end do

  reduced(:,:,:) = R_TOTYPE(M_ZERO)

  do ios = 1, this%norbsets
    !
    os => this%orbsets(ios)
    do ibatch = 1, psib%nst
      bind1 = psib%ist_idim_to_linear((/ibatch, 1/))
      bind2 = psib%ist_idim_to_linear((/ibatch, 2/))
      do im = 1,os%norbs
        ! sum_mp Vmmp <phi mp | psi >
        do imp = 1, os%norbs
          !Note here that V_{mmp} =U/2(delta_{mmp}-2n_{mpm})
          if(d%ispin /= SPINORS) then
            reduced(im, ibatch, ios) = reduced(im, ibatch, ios) + this%X(V)(im, imp, ispin, ios)*dot(1,imp, ios, ibatch)
          else
            reduced(im, bind1, ios) = reduced(im, bind1, ios) + this%X(V)(im, imp, 1, ios)*dot(1, imp, ios, ibatch)
            reduced(im, bind1, ios) = reduced(im, bind1, ios) + this%X(V)(im, imp, 3, ios)*dot(2, imp, ios, ibatch)
            reduced(im, bind2, ios) = reduced(im, bind2, ios) + this%X(V)(im, imp, 4, ios)*dot(1, imp, ios, ibatch)
            reduced(im, bind2, ios) = reduced(im, bind2, ios) + this%X(V)(im, imp, 2, ios)*dot(2, imp, ios, ibatch)
          end if
        end do
      end do      
    end do !ibatch

    !We add the intersite interaction
    if(this%intersite) then
      !Loop over nearest neighbors
      do inn = 1, os%nneighbors
        ios2 = os%map_os(inn) 

        if(psib%has_phase) then
#ifdef R_TCOMPLEX
          weight = os%phase_shift(inn, psib%ik)*os%V_ij(inn, 0)/el_per_state
#else
          !Phase can only be applied to complex wavefunctions
          ASSERT(.false.)
#endif
        else
          weight = os%V_ij(inn, 0)/el_per_state
        end if

        do ibatch = 1, psib%nst
          bind1 = psib%ist_idim_to_linear((/ibatch, 1/))
          bind2 = psib%ist_idim_to_linear((/ibatch, 2/))
          do im = 1, os%norbs
            do imp = 1, this%orbsets(ios2)%norbs
              if(d%ispin /= SPINORS) then
                reduced(im, ibatch, ios) = reduced(im, ibatch, ios) - dot(1, imp, ios2, ibatch) &
                         *this%X(n_ij)(im, imp, ispin, ios, inn)*weight
              else ! Spinors
                reduced(im, bind1, ios) = reduced(im, bind1, ios) - dot(1, imp, ios2, ibatch) &
                         *this%X(n_ij)(im, imp, 1, ios, inn)*weight
                reduced(im, bind1, ios) = reduced(im, bind1, ios) - dot(2, imp, ios2, ibatch) &
                         *this%X(n_ij)(im, imp, 3, ios, inn)*weight
                reduced(im, bind2, ios) = reduced(im, bind2, ios) - dot(1, imp, ios2, ibatch) &
                         *this%X(n_ij)(im, imp, 4, ios, inn)*weight
                reduced(im, bind2, ios) = reduced(im, bind2, ios) - dot(2, imp, ios2, ibatch) &
                         *this%X(n_ij)(im, imp, 2, ios, inn)*weight
              end if
            end do !imp
          end do !im
        end do !ibatch
      end do !inn

    end if
  end do

  !We add the orbitals properly weighted to hpsi
  do ios = 1, this%norbsets 
    os => this%orbsets(ios)
    call X(orbitalset_add_to_batch)(os, d%dim, hpsib, reduced(:,:,ios))
  end do
 
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(reduced)

  POP_SUB(lda_u_apply)
  call profiling_out(prof)

end subroutine X(lda_u_apply)

! ---------------------------------------------------------
!> This routine computes the values of the occupation matrices
! ---------------------------------------------------------
subroutine X(update_occ_matrices)(this, namespace, mesh, st, lda_u_energy, phase)
  type(lda_u_t),       intent(inout) :: this
  type(namespace_t),   intent(in)    :: namespace
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in)    :: st
  FLOAT,               intent(inout) :: lda_u_energy
  CMPLX,     optional, intent(in)    :: phase(1:mesh%np_part, st%d%kpt%start:st%d%kpt%end)

  integer :: ios, im, ik, ist, ispin, norbs, idim, inn, im2, ios2
  R_TYPE, allocatable :: psi(:,:) 
  R_TYPE, allocatable :: dot(:,:,:)
  FLOAT   :: weight, renorm_weight
  R_TYPE  :: weight_phase, renorm_weight_phase
  type(orbitalset_t), pointer :: os, os2
  type(profile_t), save :: prof
  integer :: spec_ind, spec_ind2
  FLOAT, allocatable :: muliken_charge(:)

  call profiling_in(prof, TOSTRING(X(DFTU_OCC_MATRICES)))
  
  PUSH_SUB(update_occ_matrices)

  this%X(n)(1:this%maxnorbs, 1:this%maxnorbs, 1:st%d%nspin, 1:this%norbsets) = R_TOTYPE(M_ZERO)
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(dot(1:st%d%dim, 1:this%maxnorbs, 1:this%norbsets))
  if(this%level == DFT_U_ACBN0) then
    this%X(n_alt)(1:this%maxnorbs, 1:this%maxnorbs, 1:st%d%nspin, 1:this%norbsets) = R_TOTYPE(M_ZERO)
    if(this%intersite) then
      this%X(n_ij)(1:this%maxnorbs, 1:this%maxnorbs, 1:st%d%nspin, &
                     1:this%norbsets, 1:this%maxneighbors) = R_TOTYPE(M_ZERO)
      this%X(n_alt_ij)(1:this%maxnorbs, 1:this%maxnorbs, 1:st%d%nspin, &
                     1:this%norbsets, 1:this%maxneighbors) = R_TOTYPE(M_ZERO)
      this%X(n_alt_ii)(1:2, 1:this%maxnorbs, 1:st%d%nspin, &
                     1:this%norbsets, 1:this%maxneighbors) = R_TOTYPE(M_ZERO)
    end if
    SAFE_ALLOCATE(muliken_charge(this%norbsets))
  end if

  !TODO: use symmetries of the occupation matrices
  do ik = st%d%kpt%start, st%d%kpt%end
    ispin =  st%d%get_spin_index(ik)

    do ist = st%st_start, st%st_end

      weight = st%d%kweights(ik)*st%occ(ist, ik)     
      if(weight < CNST(1.0e-10)) cycle
 
      call states_elec_get_state(st, mesh, ist, ik, psi)

      if(present(phase)) then
#ifdef R_TCOMPLEX
        call states_elec_set_phase(st%d, psi, phase(:,ik), mesh%np, .false.)
#endif
      end if

      do ios = 1, this%norbsets
        os => this%orbsets(ios)
        !We first compute the matrix elemets <orb_m |\psi>
        !taking into account phase correction if needed 
        call X(orbitalset_get_coefficients)(os, st%d%dim, psi, ik, present(phase), &
                          dot(1:st%d%dim, 1:os%norbs, ios))
      end do !ios

      !We compute the on-site occupation of the site, if needed 
      if(this%level == DFT_U_ACBN0) then
        this%renorm_occ(:,:,:,ist,ik) = M_ONE-this%acbn0_screening
        muliken_charge = M_ZERO
        do ios = 1, this%norbsets
          os => this%orbsets(ios)
          if(this%basisfromstates) then
            spec_ind = 1
          else
            spec_ind = species_index(os%spec)
          end if

          do im = 1, os%norbs 
            do idim = 1, st%d%dim
              muliken_charge(ios) = muliken_charge(ios) + abs(dot(idim,im,ios))**2 * this%acbn0_screening
            end do
          end do

          this%renorm_occ(spec_ind,os%nn,os%ll,ist,ik) = &
            this%renorm_occ(spec_ind,os%nn,os%ll,ist,ik) + muliken_charge(ios)
        end do

      end if


      !We can compute the (renormalized) occupation matrices
      do ios = 1, this%norbsets
        os => this%orbsets(ios)
        if(this%basisfromstates) then
          spec_ind = 1
        else
          spec_ind = species_index(os%spec)
        end if
        norbs = os%norbs

        do im = 1, norbs
        
          if(st%d%ispin /= SPINORS) then !Collinear case
            this%X(n)(1:norbs, im, ispin, ios) = this%X(n)(1:norbs, im, ispin, ios) &
              + weight * dot(1, im, ios) * R_CONJ(dot(1, 1:norbs, ios))
            !We compute the renomalized occupation matrices
            if(this%level == DFT_U_ACBN0) then
              renorm_weight = this%renorm_occ(spec_ind, os%nn, os%ll, ist, ik) * weight
              this%X(n_alt)(1:norbs, im, ispin, ios) = this%X(n_alt)(1:norbs, im, ispin, ios) &
                      + renorm_weight * dot(1, im, ios) * R_CONJ(dot(1, 1:norbs, ios))

              !Generalized occupation matrices
              if(this%intersite) then

                do inn = 1, os%nneighbors
                  ios2 = os%map_os(inn)
                  os2 => this%orbsets(ios2)

                  if(this%basisfromstates) then
                    spec_ind2 = 1
                  else
                    spec_ind2 = species_index(os2%spec)
                  end if
                  
                  renorm_weight = (this%renorm_occ(spec_ind, os%nn, os%ll, ist, ik) + & 
                                   this%renorm_occ(spec_ind2, os2%nn, os2%ll, ist, ik)) * weight
                  if(spec_ind == spec_ind2) then 
                    renorm_weight = M_HALF* renorm_weight
                  end if
                       
                  this%X(n_alt_ii)(1, im, ispin, ios, inn) = this%X(n_alt_ii)(1, im, ispin, ios, inn) &
                    + renorm_weight * dot(1, im, ios) * R_CONJ(dot(1, im, ios))

                  weight_phase = weight
                  renorm_weight_phase = renorm_weight
                  if(present(phase)) then
#ifdef R_TCOMPLEX
                    weight_phase = weight_phase * os%phase_shift(inn, ik)
                    renorm_weight_phase = renorm_weight_phase * os%phase_shift(inn, ik)
#else
                    ! Phase can only by applied when having complex wavefunctions
                    ASSERT(.false.)
#endif
                  end if


                  do im2 = 1, os2%norbs
              
                    if(im == 1) then
                       this%X(n_alt_ii)(2, im2, ispin, ios, inn) = this%X(n_alt_ii)(2, im2, ispin, ios, inn) &
                       + renorm_weight * dot(1, im2, ios2) * R_CONJ(dot(1, im2, ios2))
                    end if

                    this%X(n_ij)(im, im2, ispin, ios, inn) = this%X(n_ij)(im, im2, ispin, ios, inn) &
                           + weight_phase * dot(1, im2, ios2) * R_CONJ(dot(1, im, ios))
                    this%X(n_alt_ij)(im, im2, ispin, ios, inn) = this%X(n_alt_ij)(im, im2, ispin, ios, inn) &
                           + renorm_weight_phase * dot(1, im2, ios2) * R_CONJ(dot(1, im, ios))
                  end do !im2
                end do !inn
              end if
            end if

          else !Noncollinear case

            this%X(n)(1:norbs, im, 1, ios) = this%X(n)(1:norbs, im, 1, ios) &
              + weight*dot(1, im, ios) * R_CONJ(dot(1, 1:norbs, ios))
            this%X(n)(1:norbs, im, 2, ios) = this%X(n)(1:norbs, im, 2, ios) &
              + weight*dot(2, im, ios) * R_CONJ(dot(2, 1:norbs, ios))
            this%X(n)(1:norbs, im, 3, ios) = this%X(n)(1:norbs, im, 3, ios) &
              + weight*dot(1, im, ios) * R_CONJ(dot(2, 1:norbs, ios))
            this%X(n)(1:norbs, im, 4, ios) = this%X(n)(1:norbs, im, 4, ios) &
              + weight*dot(2, im, ios) * R_CONJ(dot(1, 1:norbs, ios))
            !We compute the renomalized occupation matrices
            if(this%level == DFT_U_ACBN0) then
              renorm_weight = this%renorm_occ(spec_ind,os%nn,os%ll,ist,ik)*weight
              this%X(n_alt)(1:norbs,im,1,ios) = this%X(n_alt)(1:norbs, im, 1, ios) &
                            + renorm_weight * dot(1, im, ios) * R_CONJ(dot(1, 1:norbs, ios))
              this%X(n_alt)(1:norbs,im,2,ios) = this%X(n_alt)(1:norbs, im, 2, ios) &
                            + renorm_weight * dot(2, im, ios) * R_CONJ(dot(2, 1:norbs, ios))
              this%X(n_alt)(1:norbs,im,3,ios) = this%X(n_alt)(1:norbs, im, 3, ios) &
                            + renorm_weight * dot(1, im, ios) * R_CONJ(dot(2, 1:norbs, ios))
              this%X(n_alt)(1:norbs,im,4,ios) = this%X(n_alt)(1:norbs, im, 4, ios) &
                            + renorm_weight * dot(2, im, ios) * R_CONJ(dot(1, 1:norbs, ios))

              if(this%intersite) then
                do inn = 1, os%nneighbors
                  ios2 = os%map_os(inn)
                  os2 => this%orbsets(ios2)

                  if(this%basisfromstates) then
                    spec_ind2 = 1
                  else
                    spec_ind2 = species_index(os2%spec)
                  end if

                  renorm_weight = (this%renorm_occ(spec_ind, os%nn, os%ll, ist, ik) + &
                                   this%renorm_occ(spec_ind2, os2%nn, os2%ll, ist, ik) ) * weight
                  if(spec_ind == spec_ind2) then
                    renorm_weight = M_HALF * renorm_weight
                  end if

                  weight_phase = weight
                  renorm_weight_phase = renorm_weight
                  if(present(phase)) then
#ifdef R_TCOMPLEX
                    weight_phase = weight_phase * os%phase_shift(inn, ik)
                    renorm_weight_phase = renorm_weight_phase * os%phase_shift(inn, ik)
#else
                    ! Phase can only by applied when having complex wavefunctions
                    ASSERT(.false.)
#endif
                  end if

                  this%X(n_alt_ii)(1, im, 1, ios, inn) = this%X(n_alt_ii)(1, im, 1, ios, inn) &
                    + renorm_weight * dot(1, im, ios) * R_CONJ(dot(1, im, ios))
                  this%X(n_alt_ii)(1, im, 2, ios, inn) = this%X(n_alt_ii)(1, im, 2, ios, inn) &
                    + renorm_weight * dot(2, im, ios) * R_CONJ(dot(2, im, ios))
                  this%X(n_alt_ii)(1, im, 3, ios, inn) = this%X(n_alt_ii)(1, im, 3, ios, inn) &
                    + renorm_weight * dot(1, im, ios) * R_CONJ(dot(2, im, ios))
                  this%X(n_alt_ii)(1, im, 4, ios, inn) = this%X(n_alt_ii)(1, im, 4, ios, inn) &
                    + renorm_weight * dot(2, im, ios) * R_CONJ(dot(1, im, ios))

                  do im2 = 1, os2%norbs

                    if(im == 1) then
                      this%X(n_alt_ii)(2, im2, 1, ios, inn) = this%X(n_alt_ii)(2, im2, 1, ios, inn) &  
                        + renorm_weight * dot(1, im2, ios2) * R_CONJ(dot(1, im2, ios2))
                      this%X(n_alt_ii)(2, im2, 2, ios, inn) = this%X(n_alt_ii)(2, im2, 2, ios, inn) &
                        + renorm_weight * dot(2, im2, ios2) * R_CONJ(dot(2, im2, ios2))
                      this%X(n_alt_ii)(2, im2, 3, ios, inn) = this%X(n_alt_ii)(2, im2, 3, ios, inn) & 
                        + renorm_weight * dot(1, im2, ios2) * R_CONJ(dot(2, im2, ios2))
                      this%X(n_alt_ii)(2, im2, 4, ios, inn) = this%X(n_alt_ii)(2, im2, 4, ios, inn) &  
                        + renorm_weight * dot(2, im2, ios2) * R_CONJ(dot(1, im2, ios2))
                    end if


                    this%X(n_ij)(im, im2, 1, ios, inn) = this%X(n_ij)(im, im2, 1, ios, inn) &
                           + weight_phase*dot(1, im2, ios2)*R_CONJ(dot(1, im, ios))
                    this%X(n_ij)(im, im2, 2, ios, inn) = this%X(n_ij)(im, im2, 2, ios, inn) &
                           + weight_phase*dot(2, im2, ios2)*R_CONJ(dot(2, im, ios))
                    this%X(n_ij)(im, im2, 3, ios, inn) = this%X(n_ij)(im, im2, 3, ios, inn) &
                           + weight_phase*dot(1, im2, ios2)*R_CONJ(dot(2, im, ios))
                    this%X(n_ij)(im, im2, 4, ios, inn) = this%X(n_ij)(im, im2, 4, ios, inn) &
                           + weight_phase*dot(2, im2, ios2)*R_CONJ(dot(1, im, ios))


                    this%X(n_alt_ij)(im, im2, 1, ios, inn) = this%X(n_alt_ij)(im, im2, 1, ios, inn) &
                           + renorm_weight_phase * dot(1, im2, ios2) * R_CONJ(dot(1, im, ios))
                    this%X(n_alt_ij)(im, im2, 2, ios, inn) = this%X(n_alt_ij)(im, im2, 2, ios, inn) &
                           + renorm_weight_phase * dot(2, im2, ios2) * R_CONJ(dot(2, im, ios))
                    this%X(n_alt_ij)(im, im2, 3, ios, inn) = this%X(n_alt_ij)(im, im2, 3, ios, inn) &
                           + renorm_weight_phase * dot(1, im2, ios2) * R_CONJ(dot(2, im, ios))
                    this%X(n_alt_ij)(im, im2, 4, ios, inn) = this%X(n_alt_ij)(im, im2, 4, ios, inn) &
                           + renorm_weight_phase * dot(2, im2, ios2) * R_CONJ(dot(1, im, ios))
                   end do !im2
                 end do !inn
               end if
            end if
          end if !Spinors
        end do !im
      end do !ios
    end do
  end do

  SAFE_DEALLOCATE_A(dot)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(muliken_charge)

  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp, this%X(n))
    if(this%level == DFT_U_ACBN0) then
      call comm_allreduce(st%st_kpt_mpi_grp, this%X(n_alt))
      if(this%intersite) then
        call comm_allreduce(st%st_kpt_mpi_grp, this%X(n_ij))
        call comm_allreduce(st%st_kpt_mpi_grp, this%X(n_alt_ij))
        call comm_allreduce(st%st_kpt_mpi_grp, this%X(n_alt_ii))
      end if
    end if
  end if

  if(this%level == DFT_U_ACBN0 .and. .not.this%freeze_u) then
    if(this%nspins > 1 ) then
      do ios = 1, this%norbsets
        if(this%orbsets(ios)%ndim  == 1) then
          call X(compute_ACBNO_U)(this, ios, namespace)
          call X(compute_ACBNO_V)(this, ios)
        else
          call compute_ACBNO_U_noncollinear(this, ios, namespace)
          if(this%intersite) then
            call messages_not_implemented("Intersite interaction with spinors orbitals.", namespace=namespace)
          end if
        end if
      end do
    else
      call X(compute_ACBNO_U_restricted)(this)
      call X(compute_ACBNO_V_restricted)(this)
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
  type(states_elec_t), intent(in) :: st 

  integer :: ios, imp, im, ispin, inn
  type(orbitalset_t), pointer :: os
  FLOAT :: nsigma

  PUSH_SUB(compute_dftu_energy)

  energy = M_ZERO

  do ios = 1, this%norbsets
    do ispin = 1, this%nspins
      !TODO: These are matrix operations, that could be optimized
      do im = 1, this%orbsets(ios)%norbs
        do imp = 1, this%orbsets(ios)%norbs
          energy = energy - M_HALF*this%orbsets(ios)%Ueff &
                                        *abs(this%X(n)(im, imp, ispin, ios))**2/st%smear%el_per_state
        end do
        if(ispin <= this%spin_channels) &
          energy = energy + M_HALF*this%orbsets(ios)%Ueff*TOFLOAT(this%X(n)(im, im, ispin, ios))
      end do
    end do
  end do

  if(this%intersite) then
    do ios = 1, this%norbsets
      os => this%orbsets(ios)
      do inn = 1, os%nneighbors
        do ispin = 1, this%nspins
          do im = 1, os%norbs
            do imp = 1, this%orbsets(os%map_os(inn))%norbs
              energy = energy - M_HALF*os%V_ij(inn,0)/st%smear%el_per_state &
                               *abs(this%X(n_ij)(im, imp, ispin, ios, inn))**2
            end do
          end do
        end do
      end do
    end do
  end if

  !See for instance PRB 67, 153106 (2003), Eq.(4)
  if(this%double_couting == DFT_U_AMF) then
    ASSERT(st%d%ispin /= SPINORS)
    do ios = 1, this%norbsets
      do ispin = 1, this%nspins
        nsigma = M_ZERO
        do im = 1, this%orbsets(ios)%norbs
          nsigma = nsigma + R_REAL(this%X(n)(im,im,ispin,ios))/st%smear%el_per_state
        end do

        do im = 1, this%orbsets(ios)%norbs
          energy = energy + M_HALF*this%orbsets(ios)%Ueff*nsigma*(M_ONE-nsigma/this%orbsets(ios)%norbs)
          energy = energy - M_HALF*this%orbsets(ios)%Ueff*R_REAL(this%X(n)(im, im, ispin, ios))
        end do
      end do
    end do
  end if

  POP_SUB(compute_dftu_energy)
end subroutine X(compute_dftu_energy)


! ---------------------------------------------------------
!> This routine computes the potential that, once multiplied
!> by the projector Pmm' and summed over m and m' for all the atoms
!> gives the full Hubbard potential
! ---------------------------------------------------------
subroutine X(lda_u_update_potential)(this, st)
  type(lda_u_t), intent(inout)    :: this
  type(states_elec_t), intent(in) :: st

  integer :: ios, im, ispin, norbs
  type(profile_t), save :: prof
  FLOAT :: nsigma

  call profiling_in(prof, TOSTRING(X(DFTU_POTENTIAL)))

  PUSH_SUB(lda_u_update_potential)

  this%X(V) = M_ZERO

  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    do ispin = 1, this%nspins
      do im = 1, norbs
        this%X(V)(1:norbs,im,ispin,ios) = - this%orbsets(ios)%Ueff*this%X(n)(1:norbs,im,ispin,ios)/st%smear%el_per_state
        ! Only the diagonal part in spin space (for spinors)
        if(ispin <= this%spin_channels) &
          this%X(V)(im,im,ispin,ios) = this%X(V)(im,im,ispin,ios) + M_HALF*this%orbsets(ios)%Ueff

        if(this%orbsets(ios)%alpha > CNST(1.0e-6)) then
          this%X(V)(im,im,ispin,ios) = this%X(V)(im,im,ispin,ios) + this%orbsets(ios)%alpha
        end if
      end do
    end do

    !See for instance PRB 67, 153106 (2003), Eq.(4)
    if(this%double_couting == DFT_U_AMF) then
      ASSERT(st%d%ispin /= SPINORS)
      do ispin = 1, this%nspins
        nsigma = M_ZERO
        do im = 1, norbs
          nsigma = nsigma + R_REAL(this%X(n)(im,im,ispin,ios))/st%smear%el_per_state
        end do
        do im = 1, norbs
          this%X(V)(im,im,ispin,ios) = this%X(V)(im,im,ispin,ios) + this%orbsets(ios)%Ueff &
                            *(nsigma/norbs - M_HALF)
        end do
      end do
    end if
  end do


  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp, this%X(V))
  end if

  POP_SUB(lda_u_update_potential)
  call profiling_out(prof)
end subroutine X(lda_u_update_potential)

! ---------------------------------------------------------
!> This routine computes the effective U following the expression 
!> given in Agapito et al., Phys. Rev. X 5, 011006 (2015)
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U)(this, ios, namespace)
  type(lda_u_t),     intent(inout) :: this
  integer,           intent(in)    :: ios
  type(namespace_t), intent(in)    :: namespace
  
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
            tmpU = tmpU + R_REAL(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin2,ios))
          end do
          tmpJ = tmpJ + R_REAL(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin1,ios))
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          tmpJ = tmpJ + R_REAL(this%X(n_alt)(im,imp,3,ios)*this%X(n_alt)(impp,imppp,4,ios)) &
                         +R_REAL(this%X(n_alt)(im,imp,4,ios)*this%X(n_alt)(impp,imppp,3,ios))
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
          tmpJ = tmpJ + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin1,ios))
          tmpU = tmpU + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin1,ios))
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          tmpJ = tmpJ + R_REAL(this%X(n)(im,im,3,ios)*this%X(n)(imp,imp,4,ios)) &
                            +R_REAL(this%X(n)(im,im,4,ios)*this%X(n)(imp,imp,3,ios))
        end if
      end if

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels
          if(ispin1 /= ispin2) then
            tmpU = tmpU + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin2,ios))
          end if
        end do
      end do

      if(this%nspins>this%spin_channels) then !Spinors
        if(im == imp) then
          tmpU = tmpU -(R_REAL(this%X(n)(im,im,3,ios)*this%X(n)(im,im,4,ios)) &
                       +R_REAL(this%X(n)(im,im,4,ios)*this%X(n)(im,im,3,ios)))
        end if
      end if 

      if(this%rot_inv) then
        !Rotationally invariance term
        !sum_{alpha} sum_{m,mp/=m} n^alpha_{mmp}n^alpha_{mpm}
        if(imp/=im) then
          do ispin1 = 1, this%spin_channels
            tmpJ = tmpJ + R_REAL(this%X(n)(im,imp,ispin1,ios)*this%X(n)(imp,im,ispin1,ios))
            tmpU = tmpU + R_REAL(this%X(n)(im,imp,ispin1,ios)*this%X(n)(imp,im,ispin1,ios))
          end do
          ASSERT(this%nspins==this%spin_channels)!Spinors not yet implemented
        end if
      end if

      denomJ = denomJ + tmpJ
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
          numU = numU + R_REAL(this%X(n_alt)(1,1,ispin1,ios))*R_REAL(this%X(n_alt)(1,1,ispin2,ios)) 
          denomU = denomU + R_REAL(this%X(n)(1,1,ispin1,ios))*R_REAL(this%X(n)(1,1,ispin2,ios))
        end if
      end do
    end do

    if(this%nspins>this%spin_channels) then !Spinors
      denomU = denomU - (R_REAL(this%X(n)(1,1,3,ios)*this%X(n)(1,1,4,ios)) & 
                        +R_REAL(this%X(n)(1,1,4,ios)*this%X(n)(1,1,3,ios)))
    end if

    ! We have to be careful in the case of hydrogen atom for instance 
    if(abs(denomU)> CNST(1.0e-3)) then
      this%orbsets(ios)%Ubar = (numU/denomU)*R_REAL(this%coulomb(1,1,1,1,ios))
    else
      if( abs(numU-denomU) < CNST(1.0e-3)) then
        this%orbsets(ios)%Ubar = R_REAL(this%coulomb(1,1,1,1,ios))
      else
        this%orbsets(ios)%Ubar = (numU/denomU)
        write(message(1),'(a,a)')' Small denominator value for the s orbital ', this%orbsets(ios)%Ubar
        write(message(2),'(a,a)')' to be multiplied by ',  this%coulomb(1,1,1,1,ios)
        call messages_warning(2, namespace=namespace)
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
          numU = numU + R_REAL(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios)) &
                             *this%coulomb(im,imp,impp,imppp,ios)
          numJ = numJ + R_REAL(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios)) &
                             *this%coulomb(im,imppp,impp,imp,ios)
        end do
        end do
        ! We compute the term
        ! sum_{m,mp/=m} N_{m}N_{mp}
        if(imp/=im) then
          denomJ = denomJ + R_REAL(this%X(n)(im,im,1,ios))*R_REAL(this%X(n)(imp,imp,1,ios))
          denomU = denomU + R_REAL(this%X(n)(im,im,1,ios))*R_REAL(this%X(n)(imp,imp,1,ios))
        end if
        ! We compute the term
        ! sum_{m,mp} N_{m}N_{mp}
        denomU = denomU + R_REAL(this%X(n)(im,im,1,ios))*R_REAL(this%X(n)(imp,imp,1,ios))

        if(this%rot_inv) then
          !Rotationally invariance term
          !sum_{m,mp/=m} n^alpha_{mmp}n^alpha_{mpm}
          if(imp/=im) then
            denomJ = denomJ + R_REAL(this%X(n)(im,imp,1,ios)*this%X(n)(imp,im,1,ios))
            denomU = denomU + R_REAL(this%X(n)(im,imp,1,ios)*this%X(n)(imp,im,1,ios))
          end if
        end if
      end do
      end do

      this%orbsets(ios)%Ueff = M_TWO*numU/denomU - numJ/denomJ
      this%orbsets(ios)%Ubar = M_TWO*numU/denomU
      this%orbsets(ios)%Jbar = numJ/denomJ

    else !In the case of s orbitals, the expression is different
      ! P_{mmp}P_{mpp,mppp}(m,mp|mpp,mppp)  
      numU = R_REAL(this%X(n_alt)(1,1,1,ios))**2*this%coulomb(1,1,1,1,ios)

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      denomU = R_REAL(this%X(n)(1,1,1,ios))**2

      this%orbsets(ios)%Ueff = numU/denomU
      this%orbsets(ios)%Ubar = numU/denomU
      this%orbsets(ios)%Jbar = 0
    end if
  end do

  POP_SUB(compute_ACBNO_U_restricted)  
end subroutine X(compute_ACBNO_U_restricted)

! ---------------------------------------------------------
!> This routine computes the effective V in the spin-polarized case
! ---------------------------------------------------------
subroutine X(compute_ACBNO_V)(this, ios)
  type(lda_u_t), intent(inout)    :: this
  integer,       intent(in)       :: ios

  integer :: ios2, im, imp
  integer :: inn, norbs, norbs2
  integer :: ispin1, ispin2
  FLOAT   :: numV, denomV

  if(.not.this%intersite) return

  PUSH_SUB(compute_ACBNO_V)

  norbs = this%orbsets(ios)%norbs
  do inn = 1, this%orbsets(ios)%nneighbors
    numV = M_ZERO
    denomV = M_ZERO

    ios2 = this%orbsets(ios)%map_os(inn)
    norbs2 = this%orbsets(ios2)%norbs

    do im = 1, norbs
      do imp= 1, norbs2
        do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels
          numV = numV + R_REAL(this%X(n_alt_ii)(1,im,ispin1,ios,inn))*R_REAL(this%X(n_alt_ii)(2,imp,ispin2,ios,inn))   &
                       *this%orbsets(ios)%coulomb_IIJJ(im,im,imp,imp,inn)
          if(ispin1 == ispin2) then
            numV = numV - R_REAL(this%X(n_alt_ij)(im,imp,ispin1,ios,inn)*R_CONJ(this%X(n_alt_ij)(im,imp,ispin1,ios,inn)))&
                       *this%orbsets(ios)%coulomb_IIJJ(im,im,imp,imp,inn)
          end if
        end do
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          numV = numV -(R_REAL(this%X(n_alt_ij)(im,imp,3,ios,inn)*R_CONJ(this%X(n_alt_ij)(im,imp,3,ios,inn))) &
                       +R_REAL(this%X(n_alt_ij)(im,imp,4,ios,inn)*R_CONJ(this%X(n_alt_ij)(im,imp,4,ios,inn)))) &
                        *this%orbsets(ios)%coulomb_IIJJ(im,im,imp,imp,inn)
        end if
      end do !imp

      do imp = 1,norbs2
        ! We compute the term
        ! sum_{m,mp} ( 2*N_{m}N_{mp} - n_{mmp}n_{mpm})
        do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels 
          denomV = denomV + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin2,ios2))
          if(ispin1 == ispin2) denomV = denomV - abs(this%X(n_ij)(im,imp,ispin1,ios,inn))**2
        end do
        end do
        if(this%nspins>this%spin_channels) then !Spinors
          denomV = denomV - abs(this%X(n_ij)(im,imp,3,ios,inn))**2 - abs(this%X(n_ij)(im,imp,4,ios,inn))**2
        end if
      end do !imp
    end do !im


    this%orbsets(ios)%V_ij(inn,0) = numV/denomV*M_HALF
  end do !inn

  POP_SUB(compute_ACBNO_V)
end subroutine X(compute_ACBNO_V)



! ---------------------------------------------------------
!> This routine computes the effective V in the spin-unpolarised case
! ---------------------------------------------------------
subroutine X(compute_ACBNO_V_restricted)(this)
  type(lda_u_t), intent(inout)    :: this
  
  integer :: ios, ios2, im, imp
  integer :: inn, norbs, norbs2
  FLOAT   :: numV, denomV

  if(.not.this%intersite) return

  PUSH_SUB(compute_ACBNO_V_restricted)

  do ios = 1, this%norbsets
    norbs = this%orbsets(ios)%norbs
    do inn = 1, this%orbsets(ios)%nneighbors
      numV = M_ZERO
      denomV = M_ZERO

      ios2 = this%orbsets(ios)%map_os(inn)
      norbs2 = this%orbsets(ios2)%norbs

      do im = 1, norbs
        do imp = 1, norbs2
          numV = numV + (M_TWO*R_REAL(this%X(n_alt_ii)(1,im,1,ios,inn)*this%X(n_alt_ii)(2,imp,1,ios,inn))   &
               - R_REAL(this%X(n_alt_ij)(im,imp,1,ios,inn)*R_CONJ(this%X(n_alt_ij)(im,imp,1,ios,inn))))&
                         *this%orbsets(ios)%coulomb_IIJJ(im,im,imp,imp,inn)
        end do

        do imp = 1, norbs2
          ! We compute the term
          ! sum_{m,mp} ( 2*N_{m}N_{mp} - n_{mmp}n_{mpm})
          denomV = denomV + M_TWO*R_REAL(this%X(n)(im,im,1,ios))*R_REAL(this%X(n)(imp,imp,1,ios2)) &
                        - abs(this%X(n_ij)(im,imp,1,ios,inn))**2
        end do
      end do

      this%orbsets(ios)%V_ij(inn,0) = numV/denomV * M_HALF
    end do !inn
  end do !ios

  POP_SUB(compute_ACBNO_V_restricted)  
end subroutine X(compute_ACBNO_V_restricted)

! ---------------------------------------------------------
!> This routine computes the Kanamori U, Up, and J
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U_kanamori)(this, kanamori)
  type(lda_u_t), intent(in)       :: this
  FLOAT,         intent(out)      :: kanamori(:,:)
  
  integer :: ios, im, imp, impp, imppp, norbs
  integer :: ispin1, ispin2
  FLOAT   :: numU, numUp, numJ, denomU, denomUp, denomJ
  FLOAT   :: tmpU, tmpUp, tmpJ

  PUSH_SUB(compute_ACBNO_U_kanamori)

  ASSERT(this%nspins == this%spin_channels)

  do ios = 1, this%norbsets
    norbs = this%orbsets(ios)%norbs
    numU = M_ZERO
    denomU = M_ZERO
    numUp = M_ZERO
    denomUp = M_ZERO
    numJ = M_ZERO
    denomJ = M_ZERO

    if(norbs > 1) then

    do im = 1, norbs
    do imp = 1,norbs
      do impp = 1, norbs
      do imppp = 1, norbs
        tmpU = M_ZERO
        tmpUp = M_ZERO
        tmpJ = M_ZERO

        do ispin1 = 1, this%spin_channels
          do ispin2 = 1, this%spin_channels
            tmpUp = tmpUp + R_REAL(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin2,ios))
            if(ispin1 /= ispin2) then
               tmpU = tmpU + R_REAL(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin2,ios))
            end if
          end do
          tmpJ = tmpJ + R_REAL(this%X(n_alt)(im,imp,ispin1,ios)*this%X(n_alt)(impp,imppp,ispin1,ios))
        end do

        if(im == imp .and. impp == imppp .and. im == impp) then
          numU = numU + tmpU*this%coulomb(im,imp,impp,imppp,ios)
        else
          numUp = numUp + tmpUp*this%coulomb(im,imp,impp,imppp,ios)
          numJ = numJ + tmpJ*this%coulomb(im,imppp,impp,imp,ios)
        end if
      end do
      end do

      tmpU = M_ZERO
      tmpUp = M_ZERO
      tmpJ = M_ZERO
      if(imp/=im) then
        do ispin1 = 1, this%spin_channels
          tmpUp = tmpUp + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin1,ios))
          tmpJ = tmpJ   + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin1,ios))
        end do
      end if

      do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels
          if(ispin1 /= ispin2) then
            if(im /= imp) then
              tmpUp = tmpUp + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin2,ios))
            else
              tmpU = tmpU + R_REAL(this%X(n)(im,im,ispin1,ios))*R_REAL(this%X(n)(imp,imp,ispin2,ios))
            end if
          end if
        end do
      end do

      denomU = denomU + tmpU
      denomUp = denomUp + tmpUp
      denomJ = denomJ + tmpJ

    end do
    end do

    kanamori(1,ios) = numU/denomU 
    kanamori(2,ios) = numUp/denomUp
    kanamori(3,ios) = numJ/denomJ

  else !In the case of s orbitals, the expression is different
    kanamori(1,ios) = this%orbsets(ios)%Ubar
    kanamori(2,ios) = M_ZERO
    kanamori(3,ios) = M_ZERO
  end if


  end do

  POP_SUB(compute_ACBNO_U_kanamori)  
end subroutine X(compute_ACBNO_U_kanamori)

! ---------------------------------------------------------
!> This routine computes the Kanamori U, Up, and J
! ---------------------------------------------------------
subroutine X(compute_ACBNO_U_kanamori_restricted)(this, kanamori)
  type(lda_u_t), intent(in)       :: this
  FLOAT,         intent(out)      :: kanamori(3)
  
  integer :: ios, im, imp, impp, imppp, norbs
  FLOAT   :: numU, numUp, numJ, denomU, denomUp, denomJ
  FLOAT   :: tmpU, tmpUp, tmpJ

  PUSH_SUB(compute_ACBNO_U_kanamori_restricted)

  ASSERT(this%nspins == 1)

  do ios = 1, this%norbsets
    norbs = this%orbsets(ios)%norbs
    numU = M_ZERO
    denomU = M_ZERO
    numUp = M_ZERO
    denomUp = M_ZERO
    numJ = M_ZERO
    denomJ = M_ZERO

    if(norbs > 1) then

    do im = 1, norbs
    do imp = 1,norbs
      do impp = 1, norbs
      do imppp = 1, norbs
        tmpU = M_ZERO
        tmpUp = M_ZERO
        tmpJ = M_ZERO

        tmpUp = tmpUp + M_TWO*R_REAL(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios))
        tmpU = tmpU + R_REAL(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios))
        tmpJ = tmpJ + R_REAL(this%X(n_alt)(im,imp,1,ios)*this%X(n_alt)(impp,imppp,1,ios))

        ! These are the numerator of the ACBN0 U and J
        if(im == imp .and. impp == imppp .and. im == impp) then
          numU = numU + tmpU*this%coulomb(im,imp,impp,imppp,ios)
        else
          numUp = numUp + tmpUp*this%coulomb(im,imp,impp,imppp,ios)
          numJ = numJ + tmpU*this%coulomb(im,imppp,impp,imp,ios)
        end if
      end do
      end do

      if(im /= imp) then
        denomUp = denomUp + M_TWO*R_REAL(this%X(n)(im,im,1,ios)*this%X(n)(imp,imp,1,ios))
        denomJ = denomJ + R_REAL(this%X(n)(im,im,1,ios))*R_REAL(this%X(n)(imp,imp,1,ios))
      else
        denomU = denomU + R_REAL(this%X(n)(im,im,1,ios)*this%X(n)(imp,imp,1,ios))
      end if

    end do
    end do

    kanamori(1) = numU/denomU 
    kanamori(2) = numUp/denomUp
    kanamori(3) = numJ/denomJ


  else !In the case of s orbitals, the expression is different
    kanamori(1) = this%orbsets(ios)%Ubar
    kanamori(2) = M_ZERO
    kanamori(3) = M_ZERO
  end if


  end do

  POP_SUB(compute_ACBNO_U_kanamori_restricted)  
end subroutine X(compute_ACBNO_U_kanamori_restricted)


! ---------------------------------------------------------
! TODO: Merge this with the two_body routine in system/output_me_inc.F90
subroutine X(compute_coulomb_integrals) (this, namespace, space, mesh, der, psolver)
  type(lda_u_t),       intent(inout)  :: this
  type(namespace_t),   intent(in)     :: namespace
  type(space_t),       intent(in)     :: space
  type(mesh_t),        intent(in)     :: mesh
  type(derivatives_t), intent(in)     :: der
  type(poisson_t),     intent(in)     :: psolver

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: norbs, np_sphere, ios, ip
  integer :: idone, ntodo
  FLOAT, allocatable :: tmp(:), vv(:), nn(:)
  type(orbitalset_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(DFTU_COULOMB_INTEGRALS)))

  PUSH_SUB(X(compute_coulomb_integrals))

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

    call submesh_build_global(os%sphere, space)

    select case (this%sm_poisson)
    case(SM_POISSON_DIRECT)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_DIRECT_SUM) 
    case(SM_POISSON_ISF)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_ISF)
    case(SM_POISSON_PSOLVER)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_PSOLVER)
    case(SM_POISSON_FFT)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_FFT)
    end select
 
    ijst=0
    do ist = 1, norbs
      
      do jst = 1, norbs
        if(jst > ist) cycle
        ijst=ijst+1

        !$omp parallel do
        do ip = 1,np_sphere
          nn(ip)  = TOFLOAT(os%X(orb)(ip,1,ist))*TOFLOAT(os%X(orb)(ip,1,jst))
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
            do ip = 1,np_sphere
              tmp(ip) = vv(ip)*TOFLOAT(os%X(orb)(ip,1,kst))*TOFLOAT(os%X(orb)(ip,1,lst))
            end do
            !$omp end parallel do

            this%coulomb(ist,jst,kst,lst,ios) = dsm_integrate(mesh, os%sphere, tmp(1:np_sphere), reduce = .false.)

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
    call submesh_end_cube_map(os%sphere)
    call submesh_end_global(os%sphere)
  end do !iorb

  if(mesh%parallel_in_domains) then 
    call mesh%allreduce(this%coulomb)
  end if 

  if(this%orbs_dist%parallel) then
    call comm_allreduce(this%orbs_dist%mpi_grp, this%coulomb)
  end if

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(compute_coulomb_integrals))
  call profiling_out(prof)
end subroutine X(compute_coulomb_integrals)

subroutine X(compute_periodic_coulomb_integrals)(this, namespace, space, der, mc)
  type(lda_u_t),       intent(inout)  :: this
  type(namespace_t),   intent(in)     :: namespace
  type(space_t),       intent(in)     :: space
  type(derivatives_t), intent(in)     :: der
  type(multicomm_t),   intent(in)     :: mc

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: norbs, np, ip
  integer :: idone, ntodo
  FLOAT, allocatable :: tmp(:), vv(:), nn(:)
  type(orbitalset_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(DFTU_PER_COULOMB)))

  !At the moment the basis is not spin polarized
  ASSERT(this%nspins == 1)

  PUSH_SUB(X(compute_periodic_coulomb_integrals))

  SAFE_ALLOCATE(nn(1:der%mesh%np))
  SAFE_ALLOCATE(vv(1:der%mesh%np))
  SAFE_ALLOCATE(tmp(1:der%mesh%np))

  SAFE_ALLOCATE(this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs, 1:this%norbsets))
  this%coulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%norbsets) = M_ZERO

  !Lets counts the number of orbital to treat, to display a progress bar
  ntodo = 0
  norbs = this%orbsets(1)%norbs
  ntodo= ntodo + ((norbs+1)*norbs/2)*((norbs+1)*norbs/2+1)/2
  idone = 0
  if(mpi_world%rank == 0) call loct_progress_bar(-1, ntodo)

  os => this%orbsets(1)
  norbs = os%norbs
  np = der%mesh%np  

  call poisson_init(os%poisson, namespace, space, der, mc, M_ZERO, solver=POISSON_DIRECT_SUM) !POISSON_ISF)

  ijst=0
  do ist = 1, norbs
    do jst = 1, norbs
      if(jst > ist) cycle
      ijst=ijst+1

      !$omp parallel do
      do ip = 1,np
        nn(ip)  = TOFLOAT(os%X(orb)(ip,1,ist))*TOFLOAT(os%X(orb)(ip,1,jst))
      end do
      !$omp end parallel do    

      !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
      call dpoisson_solve(os%poisson, vv, nn, all_nodes=.true.)

      klst=0
      do kst = 1, norbs
        do lst = 1, norbs
          if(lst > kst) cycle
          klst=klst+1
          if(klst > ijst) cycle

          !$omp parallel do
          do ip = 1,np
            tmp(ip) = vv(ip)*TOFLOAT(os%X(orb)(ip,1,lst))*TOFLOAT(os%X(orb)(ip,1,kst))
          end do
          !$omp end parallel do

          this%coulomb(ist,jst,kst,lst,1) = dmf_integrate(der%mesh, tmp)

          if(abs(this%coulomb(ist,jst,kst,lst,1))<CNST(1.0e-12)) then
            this%coulomb(ist,jst,kst,lst,1) = M_ZERO
          else
            this%coulomb(kst,lst,ist,jst,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(jst,ist,lst,kst,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(lst,kst,jst,ist,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(jst,ist,kst,lst,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(lst,kst,ist,jst,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(ist,jst,lst,kst,1) = this%coulomb(ist,jst,kst,lst,1)
            this%coulomb(kst,lst,jst,ist,1) = this%coulomb(ist,jst,kst,lst,1)
          end if
              
          !Update the progress bar
          idone = idone + 1
          if(mpi_world%rank == 0) call loct_progress_bar(idone, ntodo)
        end do !lst
      end do !kst
    end do !jst
  end do !ist

  call poisson_end(os%poisson)

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(compute_periodic_coulomb_integrals))
  call profiling_out(prof)
end subroutine X(compute_periodic_coulomb_integrals)

 ! ---------------------------------------------------------
 !> This routine computes [r,V_lda+u].
 ! ---------------------------------------------------------
 subroutine X(lda_u_commute_r)(this, mesh, d, namespace, ik, psi, gpsi, has_phase)
   type(lda_u_t),           intent(in) :: this
   type(mesh_t),            intent(in) :: mesh
   type(states_elec_dim_t), intent(in) :: d
   type(namespace_t),       intent(in) :: namespace
   R_TYPE,                  intent(in) :: psi(:,:)
   integer,                 intent(in) :: ik
   R_TYPE,               intent(inout) :: gpsi(:, :, :)
   logical,                 intent(in) :: has_phase !True if the wavefunction has an associated phase 

   integer :: ios, idim, idir, im, imp, is, ispin
   integer :: idim_orb, inn, ios2, el_per_state
   R_TYPE, allocatable :: dot(:,:,:), epsi(:,:), reduced(:,:,:)
   type(orbitalset_t), pointer  :: os 
   type(profile_t), save :: prof

   call profiling_in(prof, TOSTRING(X(DFTU_COMMUTE_R)))

   PUSH_SUB(lda_u_commute_r)

   if(this%double_couting /= DFT_U_FLL) then
    call messages_not_implemented("AMF double couting and commutator [r,V_u]", namespace=namespace)
   end if

   if(this%intersite .and. d%ispin == SPINORS) then
     call messages_not_implemented("Intersite interaction, spinors, and commutator [r,V_u]", namespace=namespace)
   end if

   if(.not.this%basis%submesh) then
     SAFE_ALLOCATE(epsi(1:mesh%np,1:d%dim))
   else
     SAFE_ALLOCATE(epsi(1:this%max_np,1:d%dim))
   end if
   SAFE_ALLOCATE(dot(1:d%dim,1:this%maxnorbs, 1:this%norbsets))
   SAFE_ALLOCATE(reduced(1:d%dim,1:this%maxnorbs, 1:this%norbsets))

   ispin = d%get_spin_index(ik)
   if(d%ispin == UNPOLARIZED) then
     el_per_state = 2
   else
     el_per_state = 1
   end if

   do ios = 1, this%norbsets
      ! We have to compute 
      ! hpsi> += r sum_m |phi m> sum_mp Vmmp <phi mp | psi >
      !
      ! We first compute <phi m | psi> for all orbitals of the atom
      !
      os => this%orbsets(ios)
      ! 
      call X(orbitalset_get_coefficients)(os, d%dim, psi, ik, has_phase, dot(:,:,ios))
   end do
   !
   reduced(:,:,:) = M_ZERO
   !
   do ios = 1, this%norbsets
     os => this%orbsets(ios)
     do im = 1, os%norbs
       ! sum_mp Vmmp <phi mp | psi >
       do imp = 1, os%norbs
         if(d%ispin /= SPINORS) then
           reduced(1,im, ios) = reduced(1,im, ios) + this%X(V)(im,imp,ispin,ios)*dot(1,imp, ios)
         else
           reduced(1,im, ios) = reduced(1,im, ios) + this%X(V)(im,imp,1,ios)*dot(1,imp, ios)
           reduced(1,im, ios) = reduced(1,im, ios) + this%X(V)(im,imp,3,ios)*dot(2,imp, ios)
           reduced(2,im, ios) = reduced(2,im, ios) + this%X(V)(im,imp,4,ios)*dot(1,imp, ios)
           reduced(2,im, ios) = reduced(2,im, ios) + this%X(V)(im,imp,2,ios)*dot(2,imp, ios)
         end if
       end do
     end do

     !We add the intersite interaction
     if(this%intersite) then
       ASSERT(d%ispin /= SPINORS)

       !Loop over nearest neighbors
       do inn = 1, os%nneighbors
         ios2 = os%map_os(inn)

         do im = 1,os%norbs
           do imp = 1, this%orbsets(ios2)%norbs
             if(has_phase) then
#ifdef R_TCOMPLEX
               reduced(1,im,ios) = reduced(1,im,ios) - dot(1,imp,ios2)*os%phase_shift(inn, ik) &
                 *this%zn_ij(im,imp,ispin,ios,inn)*M_HALF*os%V_ij(inn,0)/el_per_state
               reduced(1,imp,ios2) = reduced(1,imp,ios2) - dot(1, im, ios)*R_CONJ(os%phase_shift(inn, ik)) &
                 *R_CONJ(this%zn_ij(im,imp,ispin,ios,inn))*M_HALF*os%V_ij(inn,0)/el_per_state
#else
               ! Phase can only be applied to complex wavefunctions
               ASSERT(.false.)
#endif
              else
                reduced(1,im,ios) = reduced(1,im,ios) - dot(1,imp,ios2) &
                           *this%X(n_ij)(im,imp,ispin,ios,inn)*M_HALF*os%V_ij(inn,0)/el_per_state
                reduced(1,imp,ios2) = reduced(1,imp,ios2) - dot(1,im,ios) &
                           *R_CONJ(this%X(n_ij)(im,imp,ispin,ios,inn))*M_HALF*os%V_ij(inn,0)/el_per_state
              end if
            end do !imp
          end do !im
        end do !inn
      end if
    end do !ios


    do ios = 1, this%norbsets
      os => this%orbsets(ios)
      do idir = 1, mesh%sb%dim
        do idim = 1, d%dim
          idim_orb = min(idim,os%ndim)
          do im = 1, os%norbs
          !In case of phase, we have to apply the conjugate of the phase here
          if(has_phase) then
#ifdef R_TCOMPLEX
            if(.not. this%basis%submesh) then
              !If we orthogonalize, the orbital is not anymore zorb*phase
              if(.not.this%basis%orthogonalization) then
                epsi(:,idim) = R_TOTYPE(M_ZERO)
                do is = 1, os%sphere%np
                    epsi(os%sphere%map(is),idim) = epsi(os%sphere%map(is),idim) &
                       + os%sphere%x(is,idir)*os%zorb(is,idim_orb,im)*os%phase(is,ik)
                end do
                call lalg_axpy(mesh%np, reduced(idim,im,ios), epsi(1:mesh%np,idim), &
                                    gpsi(1:mesh%np,idir,idim))
              else
                !$omp parallel do
                do is = 1, mesh%np
                  epsi(is,idim) = os%sphere%x(is,idir)*os%eorb_mesh(is,im,idim_orb, ik)
                end do
                call lalg_axpy(mesh%np, reduced(idim,im,ios), epsi(1:mesh%np,idim), &
                                  gpsi(1:mesh%np,idir,idim))
              end if
            else
              !$omp parallel do
              do is = 1, os%sphere%np
                epsi(is,idim) = os%sphere%x(is,idir)*os%eorb_submesh(is,idim_orb,im,ik)
              end do
              !$omp end parallel do
              call submesh_add_to_mesh(os%sphere, epsi(1:os%sphere%np,idim), &
                                  gpsi(1:mesh%np,idir,idim), reduced(idim,im,ios))
            end if
#endif
            else
              if(.not. this%basis%submesh) then
                 !$omp parallel do
                 do is = 1, mesh%np
                    epsi(is,idim) = os%sphere%x(is,idir)*os%X(orb)(is,idim_orb,im)
                 end do
                 call lalg_axpy(mesh%np, reduced(idim,im,ios), epsi(1:mesh%np,idim), &
                                  gpsi(1:mesh%np,idir,idim))
              else
                !$omp parallel do
                do is = 1, os%sphere%np
                  epsi(is,idim) = os%sphere%x(is,idir)*os%X(orb)(is,idim_orb,im)
                end do
                call submesh_add_to_mesh(os%sphere, epsi(1:os%sphere%np,idim), &
                                    gpsi(1:mesh%np,idir,idim), reduced(idim,im,ios))
              end if
            end if
          end do !im
        end do !idim
      end do !idir
    end do !ios

   do idir = 1, mesh%sb%dim
     reduced = M_ZERO

     ! We have to compute 
     ! hpsi> -= sum_m |phi m> sum_mp Vmmp <phi mp| r | psi >
     !
     ! We first compute <phi m| r | psi> for all orbitals of the atom
     !
     !
     if(.not. this%basis%submesh) then
       do is = 1, mesh%np
         epsi(is,1) = mesh%x(is,idir)*psi(is,1)
       end do
     else
       !$omp parallel do 
       do is = 1, os%sphere%np
         epsi(is,1) = os%sphere%x(is,idir)*psi(os%sphere%map(is),1)
       end do
       !$omp end parallel do
     end if

     do ios = 1, this%norbsets
       os => this%orbsets(ios) 

       if(has_phase) then
#ifdef R_TCOMPLEX
         if(.not. this%basis%submesh) then
           do im = 1, os%norbs
             do idim = 1, d%dim
               idim_orb = min(idim,os%ndim)
               dot(idim,im,ios) = X(mf_dotp)(mesh, os%eorb_mesh(1:mesh%np,im, idim_orb,ik),&
                               epsi(1:mesh%np,idim), reduce = .false.)
             end do
           end do
         else
           do im = 1, os%norbs
             do idim = 1, d%dim
               idim_orb = min(idim,os%ndim)
               dot(idim,im,ios) = X(mf_dotp)(mesh, os%eorb_submesh(1:os%sphere%np,idim_orb,im,ik),&
                               epsi(1:os%sphere%np,idim), reduce = .false., np = os%sphere%np)
             end do
           end do
         end if
#else
         ! Phase can only be applied to complex wavefunctions
         ASSERT(.false.)
#endif
       else
         if(.not. this%basis%submesh) then
           do idim = 1, d%dim
             idim_orb = min(idim,os%ndim)
             dot(idim,im,ios) = X(mf_dotp)(mesh, os%X(orb)(1:mesh%np,idim_orb,im),&
                              epsi(1:mesh%np,idim), reduce=.false.)
           end do
         else
           do im = 1, os%norbs
             do idim = 1, d%dim
               idim_orb = min(idim,os%ndim)
               dot(idim,im,ios) = X(mf_dotp)(mesh, os%X(orb)(1:os%sphere%np,idim_orb,im),&
                               epsi(1:os%sphere%np,idim), reduce = .false., np = os%sphere%np)
             end do
           end do
         end if
      end if
    end do !ios
 
    if(mesh%parallel_in_domains) call mesh%allreduce(dot)

    do ios = 1, this%norbsets
       os => this%orbsets(ios)
       do im = 1, os%norbs
         ! sum_mp Vmmp <phi mp|r| psi >
         do imp = 1, os%norbs
           if(d%ispin /= SPINORS) then
             reduced(1,im,ios) = reduced(1,im,ios) - this%X(V)(im,imp,ispin,ios)*dot(1,imp,ios)
            else
              reduced(1,im,ios) = reduced(1,im,ios) - this%X(V)(im,imp,1,ios)*dot(1,imp,ios)
              reduced(1,im,ios) = reduced(1,im,ios) - this%X(V)(im,imp,3,ios)*dot(2,imp,ios)
              reduced(2,im,ios) = reduced(2,im,ios) - this%X(V)(im,imp,4,ios)*dot(1,imp,ios)
              reduced(2,im,ios) = reduced(2,im,ios) - this%X(V)(im,imp,2,ios)*dot(2,imp,ios)
            end if
         end do
       end do
 
       !We add the intersite interaction
       if(this%intersite) then
         ASSERT(d%ispin /= SPINORS)

         !Loop over nearest neighbors
         do inn = 1, os%nneighbors

           ios2 = os%map_os(inn)
           do im = 1,os%norbs
             do imp = 1, this%orbsets(ios2)%norbs
               if(has_phase) then
#ifdef R_TCOMPLEX
                 reduced(1,im, ios) = reduced(1,im,ios) - dot(1,imp,ios2)*os%phase_shift(inn, ik) &
                       *this%zn_ij(im,imp,ispin,ios,inn) * os%V_ij(inn,0) / el_per_state
#else
                 ! Phase can only be applied to complex wavefunctions
               ASSERT(.false.)
#endif
               else
                 reduced(1,im,ios) = reduced(1,im,ios) - dot(1,imp,ios2) &
                         *this%X(n_ij)(im,imp,ispin,ios,inn) * os%V_ij(inn,0) / el_per_state
               end if
             end do !imp
           end do !im
         end do !inn
       end if
    end do
    
     do ios = 1, this%norbsets
       os => this%orbsets(ios)
       call X(orbitalset_add_to_psi)(os, d%dim, gpsi(1:mesh%np, idir, 1:d%dim), ik, has_phase, &
                                     reduced(1:d%dim, 1:os%norbs, ios)) 
     end do
   end do !idir

   SAFE_DEALLOCATE_A(epsi)
   SAFE_DEALLOCATE_A(dot)

   POP_SUB(lda_u_commute_r)
   call profiling_out(prof)
 end subroutine X(lda_u_commute_r)

 subroutine X(lda_u_force)(this, namespace, space, mesh, st, iq, psib, grad_psib, force, phase)
   type(lda_u_t),             intent(in)    :: this
   type(namespace_t),         intent(in)    :: namespace
   type(space_t),             intent(in)    :: space
   type(mesh_t),              intent(in)    :: mesh 
   type(states_elec_t),       intent(in)    :: st
   integer,                   intent(in)    :: iq
   type(wfs_elec_t),          intent(in)    :: psib
   type(wfs_elec_t),          intent(in)    :: grad_psib(:)
   FLOAT,                     intent(inout) :: force(:, :)
   logical,                   intent(in)    :: phase

   integer :: ios, iatom, ibatch, ist, im, imp, ispin, idir
   type(orbitalset_t), pointer  :: os
   R_TYPE :: ff(1:space%dim)
   R_TYPE, allocatable :: psi(:,:), gpsi(:,:)
   R_TYPE, allocatable :: dot(:,:), gdot(:,:,:), gradn(:,:,:,:)
   FLOAT :: weight
   type(profile_t), save :: prof

   if(this%level == DFT_U_NONE) return
   !In this casem there is no contribution to the force
   if(this%basisfromstates) return

   if(this%double_couting /= DFT_U_FLL) then
    call messages_not_implemented("AMF double couting with forces", namespace=namespace)
   end if

   PUSH_SUB(X(lda_u_force))

   call profiling_in(prof, TOSTRING(X(FORCES_DFTU)))

   !TODO: Implement
   if(this%intersite) then
     message(1) = "Intersite V forces are not implemented."
     call messages_warning(1, namespace=namespace)
   end if

   SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
   SAFE_ALLOCATE(gpsi(1:mesh%np, 1:st%d%dim))
   SAFE_ALLOCATE(dot(1:st%d%dim, 1:this%maxnorbs))
   SAFE_ALLOCATE(gdot(1:st%d%dim, 1:this%maxnorbs,1:space%dim))
   SAFE_ALLOCATE(gradn(1:this%maxnorbs,1:this%maxnorbs,1:this%nspins,1:space%dim))

   ispin = st%d%get_spin_index(iq)

   do ios = 1, this%norbsets 
     os => this%orbsets(ios)
     iatom = os%iatom

     gradn(1:os%norbs,1:os%norbs,1:this%nspins,1:space%dim) = R_TOTYPE(M_ZERO)

     do ibatch = 1, psib%nst
       ist = psib%ist(ibatch)
       weight = st%d%kweights(iq)*st%occ(ist, iq)
       if(weight < CNST(1.0e-10)) cycle

       call batch_get_state(psib, ibatch, mesh%np, psi)
       !No phase here, this is already added

       !We first compute the matrix elemets <\psi | orb_m>
       !taking into account phase correction if needed   
       ! 
       call X(orbitalset_get_coefficients)(os, st%d%dim, psi, iq, phase, dot)

       do idir = 1, space%dim
         call batch_get_state(grad_psib(idir), ibatch, mesh%np, gpsi)     
         !We first compute the matrix elemets <\psi | orb_m>
         !taking into account phase correction if needed 
         ! 
         !No phase here, this is already added

         call X(orbitalset_get_coefficients)(os, st%d%dim, gpsi, iq, phase, &
                     gdot(1:st%d%dim,1:os%norbs,idir))

         if(st%d%ispin /= SPINORS) then
           do im = 1, os%norbs
             gradn(1:os%norbs,im,ispin,idir) = gradn(1:os%norbs,im,ispin,idir) &
                                          + weight*(gdot(1,1:os%norbs,idir)*R_CONJ(dot(1,im)) &
                                                   +R_CONJ(gdot(1,im,idir))*dot(1,1:os%norbs))
           end do
         else
           do im = 1, os%norbs
             do ispin = 1, this%spin_channels
               gradn(1:os%norbs,im,ispin,idir) = gradn(1:os%norbs,im,ispin,idir) &
                                            + weight*(gdot(ispin,1:os%norbs,idir)*R_CONJ(dot(ispin,im)) &
                                                     +R_CONJ(gdot(ispin,im,idir))*dot(ispin,1:os%norbs))
             end do
             gradn(1:os%norbs,im,3,idir) = gradn(1:os%norbs,im,3,idir) &
                                            + weight*(gdot(1,1:os%norbs,idir)*R_CONJ(dot(2,im)) &
                                                     +R_CONJ(gdot(2,im,idir))*dot(1,1:os%norbs))
             gradn(1:os%norbs,im,4,idir) = gradn(1:os%norbs,im,4,idir) &
                                            + weight*(gdot(2,1:os%norbs,idir)*R_CONJ(dot(1,im)) &
                                                     +R_CONJ(gdot(1,im,idir))*dot(2,1:os%norbs))
           end do
           
         end if
       end do !idir
       
     end do !ibatch

     if(st%d%ispin /= SPINORS) then
       ff(1:space%dim) = M_ZERO
       do im = 1, os%norbs
         do imp = 1, os%norbs
          ff(1:space%dim) = ff(1:space%dim) - this%X(n)(imp,im,ispin,ios)/st%smear%el_per_state*gradn(im,imp,ispin,1:space%dim)
         end do !imp
       ff(1:space%dim) = ff(1:space%dim) + CNST(0.5)*gradn(im, im, ispin,1:space%dim)
       end do !im
     else
       ff(1:space%dim) = M_ZERO
       do ispin = 1, st%d%nspin
         do im = 1, os%norbs
           do imp = 1, os%norbs
             !We use R_CONJ to get n(imp,im, sigmap, sigma) from n(im,imp, sigma,sigmap)
             ff(1:space%dim) = ff(1:space%dim) - R_CONJ(this%X(n)(im,imp,ispin,ios))/st%smear%el_per_state &
                                           *gradn(im,imp,ispin,1:space%dim)
           end do !imp 
           if(ispin <= this%spin_channels) &
             ff(1:space%dim) = ff(1:space%dim) + CNST(0.5)*gradn(im, im, ispin,1:space%dim)
         end do !im
       end do !ispin
     end if

     ! We convert the force to Cartesian coordinates
     ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
     if (space%is_periodic() .and. mesh%sb%latt%nonorthogonal) then
        ff(1:space%dim) = matmul(mesh%sb%latt%klattice_primitive(1:space%dim, 1:space%dim), ff(1:space%dim))
      end if

     force(1:space%dim, iatom) = force(1:space%dim, iatom) - os%Ueff*TOFLOAT(ff(1:space%dim))
   end do !ios

   SAFE_DEALLOCATE_A(psi)
   SAFE_DEALLOCATE_A(gpsi)
   SAFE_DEALLOCATE_A(dot)
   SAFE_DEALLOCATE_A(gdot)
   SAFE_DEALLOCATE_A(gradn)

   call profiling_out(prof)

   POP_SUB(X(lda_u_force))
 end subroutine X(lda_u_force)

 ! ---------------------------------------------------------
  subroutine X(lda_u_set_occupations)(this, occ)
    type(lda_u_t),  intent(inout) :: this
    R_TYPE,         intent(in)    :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs, inn, ios2

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

    if(this%level == DFT_U_ACBN0) then
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

      if(this%intersite) then
        do ios = 1, this%norbsets
          do inn = 1, this%orbsets(ios)%nneighbors
            ios2 = this%orbsets(ios)%map_os(inn)
            do ispin = 1, this%nspins
              do im = 1, this%orbsets(ios)%norbs
                do imp = 1,this%orbsets(ios2)%norbs
                  ind = ind + 1
                  this%X(n_ij)(im,imp,ispin,ios,inn) = occ(ind)
                  ind = ind + 1
                  this%X(n_alt_ij)(im,imp,ispin,ios,inn) = occ(ind)
                end do
              end do
            end do
          end do
        end do
      end if
    end if

    POP_SUB(X(lda_u_set_occupations))
  end subroutine X(lda_u_set_occupations)

   ! ---------------------------------------------------------
  subroutine X(lda_u_get_occupations)(this, occ)
    type(lda_u_t),  intent(in) :: this
    R_TYPE,      intent(inout) :: occ(:)

    integer :: ios, ispin, im, imp, ind, norbs, inn, ios2

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

    if(this%level == DFT_U_ACBN0) then
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

      if(this%intersite) then
        do ios = 1, this%norbsets
          do inn = 1, this%orbsets(ios)%nneighbors
            ios2 = this%orbsets(ios)%map_os(inn)
            do ispin = 1, this%nspins
              do im = 1, this%orbsets(ios)%norbs
                do imp = 1,this%orbsets(ios2)%norbs
                  ind = ind + 1
                  occ(ind) = this%X(n_ij)(im,imp,ispin,ios,inn)
                  ind = ind + 1
                  occ(ind) = this%X(n_alt_ij)(im,imp,ispin,ios,inn)
                end do
              end do
            end do
          end do
        end do 
      end if
    end if

    POP_SUB(X(lda_u_get_occupations))
  end subroutine X(lda_u_get_occupations)

  ! ---------------------------------------------------------
  subroutine X(lda_u_allocate)(this, st)
    type(lda_u_t),  intent(inout)  :: this
    type(states_elec_t), intent(in):: st

    integer :: maxorbs, nspin

    PUSH_SUB(X(lda_u_allocate))

    maxorbs = this%maxnorbs
    nspin = this%nspins

    SAFE_ALLOCATE(this%X(n)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
    this%X(n)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
    SAFE_ALLOCATE(this%X(V)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
    this%X(V)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)

    !In case we use the ab-initio scheme, we need to allocate extra resources
    if(this%level == DFT_U_ACBN0) then
      SAFE_ALLOCATE(this%X(n_alt)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets))
      this%X(n_alt)(1:maxorbs,1:maxorbs,1:nspin,1:this%norbsets) = R_TOTYPE(M_ZERO)
      SAFE_ALLOCATE(this%renorm_occ(this%nspecies,0:5,0:(MAX_L-1),st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end))
      this%renorm_occ(this%nspecies,0:5,0:(MAX_L-1),st%st_start:st%st_end,st%d%kpt%start:st%d%kpt%end) = M_ZERO
    end if

    POP_SUB(X(lda_u_allocate))
  end subroutine X(lda_u_allocate)
