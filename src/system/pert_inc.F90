!! Copyright (C) 2007 the octopus team
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
!! $Id: pert_inc.F90 2548 2006-11-06 21:42:27Z xavier $

! --------------------------------------------------------------------------
subroutine X(pert_apply)(this, gr, geo, hm, ik, f_in, f_out)
! Returns f_out = H' f_in, where H' is perturbation Hamiltonian
! Note that e^ikr phase is applied to f_in, then is removed afterward
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(inout) :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:, :)
  R_TYPE,               intent(out)   :: f_out(:, :)

  R_TYPE, allocatable :: f_in_copy(:, :)
  logical :: apply_kpoint
  integer :: ip, idim

  PUSH_SUB(X(pert_apply))

  call profiling_in(prof, "PERT_APPLY")

  ASSERT(this%dir .ne. -1)

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
    SAFE_ALLOCATE(f_in_copy(1:gr%mesh%np_part, 1:hm%d%dim))
    do idim = 1, hm%d%dim
      call lalg_copy(gr%mesh%np, f_in(:, idim), f_in_copy(:, idim))
      call X(derivatives_set_bc(gr%der, f_in_copy(:, idim)))
    end do
  endif
  ! no derivatives in electric, so ghost points not needed

  apply_kpoint = associated(hm%phase) .and. this%pert_type /= PERTURBATION_ELECTRIC
  ! electric does not need it since (e^-ikr)r(e^ikr) = r

  if (apply_kpoint) then
    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np_part)
      f_in_copy(ip, idim) = hm%phase(ip, ik) * f_in_copy(ip, idim)
    end forall
  endif

  select case(this%pert_type)
    case(PERTURBATION_ELECTRIC)
      call electric()

    case(PERTURBATION_MAGNETIC)
      call magnetic()

    case(PERTURBATION_IONIC)
      call ionic()

    case(PERTURBATION_KDOTP)
      call kdotp()

    case(PERTURBATION_NONE)
      call none()

  end select
  
  if (apply_kpoint) then
    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np)
      f_out(ip, idim) = conjg(hm%phase(ip, ik)) * f_out(ip, idim)
    end forall
  endif

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
    SAFE_DEALLOCATE_A(f_in_copy)
  endif

  call profiling_out(prof)
  POP_SUB(X(pert_apply))

contains

  ! --------------------------------------------------------------------------

  subroutine none()

    PUSH_SUB(X(pert_apply).none)
    f_out(1:gr%mesh%np, 1:hm%d%dim) = M_ZERO
    POP_SUB(X(pert_apply).none)

  end subroutine none

  ! --------------------------------------------------------------------------
  subroutine electric()

    PUSH_SUB(X(pert_apply).electric)

    if(this%dir <= gr%sb%periodic_dim) then
      message(1) = "Internal error: cannot apply electric perturbation in periodic direction."
      call messages_fatal(1)
    endif

    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np)
      f_out(ip, idim) = f_in(ip, idim) * gr%mesh%x(ip, this%dir)
    end forall

    POP_SUB(X(pert_apply).electric)

  end subroutine electric

  ! --------------------------------------------------------------------------
  subroutine kdotp()
  ! perturbation is grad + [V,r]
    R_TYPE, allocatable :: grad(:, :, :)
    integer :: iatom

    PUSH_SUB(X(pert_apply).kdotp)

    SAFE_ALLOCATE(grad(1:gr%mesh%np, 1:gr%sb%dim, 1:hm%d%dim))

    do idim = 1, hm%d%dim
      call X(derivatives_grad) (gr%der, f_in_copy(:, idim), grad(:, :, idim), set_bc = .false.)
      ! set_bc done already separately
    end do

    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) f_out(ip, idim) = grad(ip, this%dir, idim)

    ! i delta_H_k = i (-i*grad + k) . delta_k
    ! representation on psi is just grad . delta_k
    ! note that second-order term is left out
    if(this%use_nonlocalpps) then
      do iatom = 1, geo%natoms
        if(species_is_ps(geo%atom(iatom)%spec)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, this%dir, ik, f_in_copy, f_out)
        end if
      end do
    endif
    
    SAFE_DEALLOCATE_A(grad)
    POP_SUB(X(pert_apply).kdotp)
    
  end subroutine kdotp

  ! --------------------------------------------------------------------------
  subroutine magnetic()
    R_TYPE, allocatable :: lf(:,:), vrnl(:,:,:)
    R_TYPE :: cross(1:MAX_DIM), vv(1:MAX_DIM)
    FLOAT :: xx(1:MAX_DIM)
    integer :: iatom, idir
    
    PUSH_SUB(X(pert_apply).magnetic)

    SAFE_ALLOCATE(lf(1:gr%mesh%np, 1:gr%sb%dim))

    do idim = 1, hm%d%dim
      ! Note that we leave out the term 1/P_c
      call X(physics_op_L)(gr%der, f_in_copy(:, idim), lf, set_bc = .false.)
      f_out(1:gr%mesh%np, idim) = M_HALF*lf(1:gr%mesh%np, this%dir)
    end do

    SAFE_DEALLOCATE_A(lf)

    if(this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
      SAFE_ALLOCATE(vrnl(1:gr%mesh%np, 1:hm%d%dim, 1:gr%sb%dim))
      
      do iatom = 1, geo%natoms

        vrnl = M_ZERO
        do idir = 1, gr%mesh%sb%dim
          if(this%dir == idir) cycle ! this direction is not used in the cross product
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir, ik, f_in_copy, vrnl(:, :, idir))
        end do
        
        xx(1:gr%mesh%sb%dim) = geo%atom(iatom)%x(1:gr%mesh%sb%dim)
        
        do idim = 1, hm%d%dim
          do ip = 1, gr%mesh%np
            
            if(this%gauge == GAUGE_ICL) xx(1:gr%mesh%sb%dim) = gr%mesh%x(ip, 1:gr%mesh%sb%dim)
            
            vv(1:gr%mesh%sb%dim) = vrnl(ip, idim, 1:gr%mesh%sb%dim)
            
            cross(1) = xx(2) * vv(3) - xx(3) * vv(2)
            cross(2) = xx(3) * vv(1) - xx(1) * vv(3)
            cross(3) = xx(1) * vv(2) - xx(2) * vv(1)
            
#if !defined(R_TCOMPLEX)
            f_out(ip, idim) = f_out(ip, idim) + M_HALF*cross(this%dir)
#else
            f_out(ip, idim) = f_out(ip, idim) - M_zI*M_HALF*cross(this%dir)
#endif
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(vrnl)
    end if
    
    POP_SUB(X(pert_apply).magnetic)
    
  end subroutine magnetic
  
  ! --------------------------------------------------------------------------
  subroutine ionic
    integer :: iatom, idir
    R_TYPE, allocatable  :: tmp(:)

    PUSH_SUB(X(pert_apply).ionic)
    SAFE_ALLOCATE(tmp(1:gr%mesh%np))

    ASSERT(hm%d%dim == 1)

    f_out(1:gr%mesh%np, 1) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim

        if (this%ionic%pure_dir .and. iatom /= this%atom1 .and. idir /= this%dir) cycle

        call X(ionic_perturbation)(this, gr, geo, hm, ik, f_in_copy(:, 1), tmp, iatom, idir)
        
        call lalg_axpy(gr%mesh%np, this%ionic%mix1(iatom, idir), tmp, f_out(:, 1))

      end do
    end do
      
    SAFE_DEALLOCATE_A(tmp)
    POP_SUB(X(pert_apply).ionic)

  end subroutine ionic

end subroutine X(pert_apply)

  ! --------------------------------------------------------------------------
subroutine X(ionic_perturbation)(this, gr, geo, hm, ik, f_in, f_out, iatom, idir)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(inout) :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)
  integer,              intent(in)    :: iatom, idir    

  ! FIX ME: may need to tell derivatives_perform not to apply boundary conditions 
  ! more things about ghost points may need to be done

  R_TYPE, allocatable :: grad(:,:), fin(:, :), fout(:, :)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  PUSH_SUB(X(ionic_perturbation))

  SAFE_ALLOCATE(vloc(1:gr%mesh%np))
  vloc(1:gr%mesh%np) = M_ZERO
  call epot_local_potential(hm%ep, gr%der, gr%dgrid, geo, iatom, vloc, CNST(0.0))

  SAFE_ALLOCATE(fin(1:gr%mesh%np_part, 1:1))
  call lalg_copy(gr%mesh%np_part, f_in, fin(:, 1))

  !d^T v |f>
  SAFE_ALLOCATE(fout(1:gr%mesh%np_part, 1:1))
  forall(ip = 1:gr%mesh%np) fout(ip, 1) = vloc(ip)*fin(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, fin, fout, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fout(:,1), f_out)

  !v d |f>
  SAFE_ALLOCATE(grad(1:gr%mesh%np, 1:1))
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), grad(:,1))
  forall(ip = 1:gr%mesh%np) fout(ip, 1) = vloc(ip)*grad(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, grad, fout, ik)
  forall(ip = 1:gr%mesh%np) f_out(ip) = -f_out(ip) + fout(ip, 1)

  SAFE_DEALLOCATE_A(grad)
  SAFE_DEALLOCATE_A(fin)
  SAFE_DEALLOCATE_A(fout)
  SAFE_DEALLOCATE_A(vloc)
  POP_SUB(X(ionic_perturbation))

end subroutine X(ionic_perturbation)

! --------------------------------------------------------------------------
subroutine X(pert_apply_order_2) (this, gr, geo, hm, ik, f_in, f_out)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(inout) :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:, :)
  R_TYPE,               intent(out)   :: f_out(:, :)

  integer :: ip, idim
  R_TYPE, allocatable :: f_in_copy(:,:)
  logical :: apply_kpoint

  PUSH_SUB(X(pert_apply_order_2))

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
    SAFE_ALLOCATE(f_in_copy(1:gr%mesh%np_part, 1:hm%d%dim))
    do idim = 1, hm%d%dim
      call lalg_copy(gr%mesh%np, f_in(:, idim), f_in_copy(:, idim))
      call X(derivatives_set_bc(gr%der, f_in_copy(:, idim)))
    end do
  endif
  ! no derivatives in electric, so ghost points not needed

  apply_kpoint = associated(hm%phase) .and. this%pert_type /= PERTURBATION_ELECTRIC &
    .and. this%pert_type /= PERTURBATION_KDOTP
  ! electric does not need it since (e^-ikr)r(e^ikr) = r
  ! kdotp has the perturbation written in terms of the periodic part with the phase

  if (apply_kpoint) then
    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np_part)
      f_in_copy(ip, idim) = hm%phase(ip, ik) * f_in_copy(ip, idim)
    end forall
  endif

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:gr%mesh%np, 1:hm%d%dim) = R_TOTYPE(M_ZERO)
  case(PERTURBATION_IONIC)
    call ionic()
  case(PERTURBATION_MAGNETIC)
    call magnetic()
  case(PERTURBATION_KDOTP)
    call kdotp()
  case(PERTURBATION_NONE)
    f_out(1:gr%mesh%np, 1:hm%d%dim) = R_TOTYPE(M_ZERO)
  end select

  if (apply_kpoint) then
    forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np)
      f_out(ip, idim) = conjg(hm%phase(ip, ik)) * f_out(ip, idim)
    end forall
  endif

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
    SAFE_DEALLOCATE_A(f_in_copy)
  endif

  POP_SUB(X(pert_apply_order_2))

contains

  ! --------------------------------------------------------------------------
  subroutine magnetic()
    R_TYPE, allocatable :: f_in2(:, :, :), dnl(:, :, :), vrnl(:,:), xf(:, :)
    R_TYPE :: cross1(1:3), bdir(1:MAX_DIM, 2)
    FLOAT  :: rdelta
    R_TYPE :: contr

    integer :: iatom, idir, idir2, ip, idim

    PUSH_SUB(X(pert_apply_order_2).magnetic)

    do idim = 1, hm%d%dim
      do ip = 1, gr%mesh%np
        rdelta = sum(gr%mesh%x(ip, 1:MAX_DIM)**2)*ddelta(this%dir, this%dir2)
        f_out(ip, idim) = M_FOURTH*(rdelta - gr%mesh%x(ip, this%dir)*gr%mesh%x(ip, this%dir2))*f_in_copy(ip, idim)
      end do
    end do

    ! gauge correction
    apply_gauge: if(this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
      bdir(1:MAX_DIM, 1:2) = M_ZERO
      bdir(this%dir,  1)   = M_ONE
      bdir(this%dir2, 2)   = M_ONE

      SAFE_ALLOCATE(f_in2(1:gr%mesh%np_part, 1:hm%d%dim, 1:gr%mesh%sb%dim))
      SAFE_ALLOCATE( vrnl(1:gr%mesh%np_part, 1:hm%d%dim))
      SAFE_ALLOCATE(  dnl(1:gr%mesh%np, 1:hm%d%dim, 1:gr%mesh%sb%dim))
      SAFE_ALLOCATE(   xf(1:gr%mesh%np, 1:hm%d%dim))

      f_in2  = R_TOTYPE(M_ZERO)
      atoms: do iatom = 1, geo%natoms

        ! This calculates f_in2 = (B x r) f_in_copy
        do ip = 1, gr%mesh%np
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(gr%mesh%x(ip, :)))
          end select

          do idim = 1,hm%d%dim
            f_in2(ip, idim, 1:gr%sb%dim) = cross1(1:gr%sb%dim)*f_in_copy(ip, idim)
          end do
        end do

        ! let us now get sum_beta Dnl f_in2
        dnl = R_TOTYPE(M_ZERO)
        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            vrnl = M_ZERO
            !calculate dnl |f_in2> = -[x,vnl] |f_in2>
            call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir2, ik, f_in2(:, :, idir2), vrnl)

            ! -x vnl |f>
            forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) 
              dnl(ip, idim, idir) = dnl(ip, idim, idir) - gr%mesh%x(ip, idir)*vrnl(ip, idim)
            end forall

            ! vnl x |f>
            forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) xf(ip, idim) = gr%mesh%x(ip, idir) * f_in2(ip, idim, idir2)

            vrnl = M_ZERO
            call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir2, ik, xf, vrnl)

            forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) dnl(ip, idim, idir) = dnl(ip, idim, idir) + vrnl(ip, idim)
          end do
        end do

        do ip = 1, gr%mesh%np
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(gr%mesh%x(ip, :)))
          end select

          do idim = 1, hm%d%dim
            contr = M_ZERO
            do idir = 1, gr%sb%dim
              contr = contr + cross1(idir)*dnl(ip, idim, idir)
            end do
            f_out(ip, idim) = f_out(ip, idim) + M_FOURTH * contr
          end do
        end do
      end do atoms

      SAFE_DEALLOCATE_A(f_in2)
      SAFE_DEALLOCATE_A(vrnl)
      SAFE_DEALLOCATE_A(dnl)
      SAFE_DEALLOCATE_A(xf)
    end if apply_gauge

    POP_SUB(X(pert_apply_order_2).magnetic)

  end subroutine magnetic

  ! --------------------------------------------------------------------------
  subroutine ionic
    integer :: iatom, idir, jdir
    R_TYPE, allocatable  :: tmp(:)
    
    PUSH_SUB(X(pert_apply_order_2).ionic)

    ASSERT(hm%d%dim == 1)

    SAFE_ALLOCATE(tmp(1:gr%mesh%np))
    
    f_out(1:gr%mesh%np, 1) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        do jdir = 1, gr%mesh%sb%dim
          
          if (this%ionic%pure_dir &
               .and. iatom /= this%atom1 .and. idir /= this%dir &
               .and. iatom /= this%atom2 .and. jdir /= this%dir2) cycle

          call X(ionic_perturbation_order_2)(this, gr, geo, hm, ik, f_in_copy(:, 1), tmp, iatom, idir, jdir)
          
          call lalg_axpy(gr%mesh%np, this%ionic%mix1(iatom, idir)*this%ionic%mix2(iatom, jdir), tmp, f_out(:, 1))
          
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
    POP_SUB(X(pert_apply_order_2).ionic)

  end subroutine ionic


  ! --------------------------------------------------------------------------
  ! d^2/dki dkj (-(1/2) ki kj [ri,[rj,H]])
  ! for i  = j : 1 - [ri,[rj,Vnl]]
  ! for i != j : -(1/2) [ri,[rj,Vnl]]
  ! Ref: Eq. 3 from M Cardona and FH Pollak, Phys. Rev. 142, 530-543 (1966)
  subroutine kdotp
    integer :: iatom
    R_TYPE, allocatable :: cpsi(:,:)

    PUSH_SUB(X(pert_apply_order_2).kdotp)

    f_out(1:gr%mesh%np, 1:hm%d%dim) = M_ZERO
    SAFE_ALLOCATE(cpsi(1:gr%mesh%np_part, 1:hm%d%dim))
    cpsi(1:gr%mesh%np_part, 1:hm%d%dim) = M_ZERO
    
    if(this%use_nonlocalpps) then
      do iatom = 1, geo%natoms
        if(species_is_ps(geo%atom(iatom)%spec)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, this%dir, ik, f_in_copy, cpsi(:, :))
        endif
      enddo

      forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np)
        f_out(ip, idim) = f_out(ip, idim) + gr%mesh%x(ip, this%dir2) * cpsi(ip, idim) - cpsi(ip, idim) * gr%mesh%x(ip, this%dir2)
      end forall
    endif

    if(this%dir == this%dir2) then
      ! add delta_ij
      forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) f_out(ip, idim) = - f_out(ip, idim) + f_in_copy(ip, idim)
    else
      forall(idim = 1:hm%d%dim, ip = 1:gr%mesh%np) f_out(ip, idim) = - M_HALF * f_out(ip, idim)
    endif

    SAFE_DEALLOCATE_A(cpsi)
    POP_SUB(X(pert_apply_order_2).kdotp)
  end subroutine kdotp

end subroutine X(pert_apply_order_2)


! --------------------------------------------------------------------------
subroutine X(ionic_perturbation_order_2) (this, gr, geo, hm, ik, f_in, f_out, iatom, idir, jdir)
  type(pert_t),        intent(in)    :: this
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: f_in(:)
  R_TYPE,              intent(out)   :: f_out(:)
  integer,             intent(in)    :: iatom, idir, jdir

  ! FIX ME: may need to tell derivatives_oper not to apply boundary conditions

  R_TYPE, allocatable :: fin(:, :)
  R_TYPE, allocatable :: tmp1(:, :), tmp2(:,:)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  PUSH_SUB(X(ionic_perturbation_order_2))

  SAFE_ALLOCATE( fin(1:gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmp1(1:gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(tmp2(1:gr%mesh%np_part, 1:1))
  SAFE_ALLOCATE(vloc(1:gr%mesh%np))

  forall(ip = 1:gr%mesh%np) vloc(ip) = M_ZERO
  call epot_local_potential(hm%ep, gr%der, gr%dgrid, geo, iatom, vloc, CNST(0.0))

  call lalg_copy(gr%mesh%np_part, f_in, fin(:, 1))
   
  !di^T dj^T v |f>
  forall(ip = 1:gr%mesh%np) tmp1(ip, 1) = vloc(ip)*fin(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, fin, tmp1, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, tmp1(:,1), tmp2(:,1))
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp2(:,1), f_out)

  !di^T v dj |f>
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, fin(:,1), tmp1(:,1))
  forall(ip = 1:gr%mesh%np) tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, tmp2(:,1), tmp1(:,1))
  forall(ip = 1:gr%mesh%np) f_out(ip) = f_out(ip) - tmp1(ip, 1)

  !dj^T v di |f>
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  forall(ip = 1:gr%mesh%np) tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp2(:,1), tmp1(:,1))
  forall(ip = 1:gr%mesh%np) f_out(ip) = f_out(ip) - tmp1(ip, 1)

  !v di dj |f>
  call X(derivatives_perform)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  call X(derivatives_perform)(gr%der%grad(jdir), gr%der, tmp1(:,1), tmp2(:,1))
  forall(ip = 1:gr%mesh%np) tmp1(ip, 1) = vloc(ip)*tmp2(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp2, tmp1, ik)
  forall(ip = 1:gr%mesh%np) f_out(ip) = f_out(ip) + tmp1(ip, 1)

  POP_SUB(X(ionic_perturbation_order_2))

end subroutine X(ionic_perturbation_order_2)

! --------------------------------------------------------------------------
subroutine X(ionic_pert_matrix_elements_2)(this, gr, geo, hm, ik, st, psi, vib, factor, matrix)
  type(pert_t),        intent(in)    :: this
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(inout) :: hm
  integer,             intent(in)    :: ik
  type(states_t),      intent(in)    :: st
  R_TYPE,              intent(inout) :: psi(:, :, :)
  type(vibrations_t),  intent(in)    :: vib
  R_TYPE,              intent(in)    :: factor
  R_TYPE,              intent(inout) :: matrix(:, :)

  integer :: ist, idim, ip
  integer :: imat, jmat, iatom, idir, jdir
  FLOAT, allocatable :: vloc(:)
  R_TYPE, allocatable :: gpsi(:, :, :), g2psi(:, :, :, :), tmp1(:, :)
  R_TYPE :: dot

  PUSH_SUB(X(ionic_pert_matrix_elements_2))

  SAFE_ALLOCATE( vloc(1:gr%mesh%np))
  SAFE_ALLOCATE( gpsi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim))
  SAFE_ALLOCATE(g2psi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim, 1:gr%sb%dim))
  SAFE_ALLOCATE( tmp1(1:gr%mesh%np, 1:st%d%dim))

  do ist = 1, st%nst
    do idim = 1, st%d%dim
      call X(derivatives_grad)(gr%der, psi(:, idim, ist), gpsi(:, idim, :))
      do idir = 1, gr%sb%dim
        call X(derivatives_grad)(gr%der, gpsi(:, idim, idir), g2psi(:, idim, idir, :))
      end do
    end do
    
    do imat = 1, vib%num_modes
      iatom = vibrations_get_atom(vib, imat)
      idir  = vibrations_get_dir (vib, imat)

      forall(ip = 1:gr%mesh%np) vloc(ip) = M_ZERO
      call epot_local_potential(hm%ep, gr%der, gr%dgrid, geo, iatom, vloc, CNST(0.0))

      do jdir = 1, gr%sb%dim
        jmat = vibrations_get_index(vib, iatom, jdir)

        dot = M_ZERO

        !2<f|dj^T v di |f>
        forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np) tmp1(ip, idim) = vloc(ip)*gpsi(ip, idim, idir)
        call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, st%d%dim, gpsi(:, :, idir), tmp1, ik)
        dot = dot + CNST(2.0)*X(mf_dotp)(gr%mesh, st%d%dim, gpsi(:, :, jdir), tmp1)

        !2<f|di^T dj^T v |f> 
        forall (idim = 1:st%d%dim, ip = 1:gr%mesh%np) tmp1(ip, idim) = vloc(ip)*psi(ip, idim, ist)
        call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, st%d%dim, psi(:, :, ist), tmp1, ik)
        dot = dot + CNST(2.0)*X(mf_dotp)(gr%mesh, st%d%dim, g2psi(:, :, idir, jdir), tmp1)

        dot = dot*st%occ(ist, ik)*factor
        
        matrix(imat, jmat) = matrix(imat, jmat) + dot

      end do

    end do
  end do
  
  POP_SUB(X(ionic_pert_matrix_elements_2))
end subroutine X(ionic_pert_matrix_elements_2)

! --------------------------------------------------------------------------
! This routine includes occupations for psib if pert_order == 2, correct if used
! as <psi(0)|H(2)|psi(0)>. It does not include occupations if pert_order == 1,
! correct if used as <psi(0)|H(1)|psi(1)> since the LR wavefunctions include the
! occupation. This routine must be modified if used differently than these two
! ways.
subroutine X(pert_expectation_density) (this, gr, geo, hm, st, psia, psib, density, pert_order)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(inout) :: hm
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  R_TYPE,               intent(out)   :: density(:)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: pertpsib(:, :)
  integer :: ik, ist, idim, order
  FLOAT   :: ikweight

  PUSH_SUB(X(pert_expectation_density))

  SAFE_ALLOCATE(pertpsib(1:gr%mesh%np, 1:st%d%dim))

  order = 1
  if(present(pert_order)) order = pert_order
  ASSERT(order == 1 .or. order == 2)

  density(1:gr%mesh%np) = R_TOTYPE(M_ZERO)

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      if(order == 1) then 
        call X(pert_apply)(this, gr, geo, hm, ik, psib(:, :, ist, ik), pertpsib)
        ikweight = st%d%kweights(ik) * st%smear%el_per_state
      else
        call X(pert_apply_order_2)(this, gr, geo, hm, ik, psib(:, :, ist, ik), pertpsib)
        ikweight = st%d%kweights(ik) * st%occ(ist, ik)
      end if

      do idim = 1, st%d%dim
        density(1:gr%mesh%np) = density(1:gr%mesh%np) + ikweight * &
          R_CONJ(psia(1:gr%mesh%np, idim, ist, ik)) * pertpsib(1:gr%mesh%np, idim)
      end do

    end do
  end do

  SAFE_DEALLOCATE_A(pertpsib)
  POP_SUB(X(pert_expectation_density))

end subroutine X(pert_expectation_density)

! --------------------------------------------------------------------------
R_TYPE function X(pert_expectation_value) (this, gr, geo, hm, st, psia, psib, pert_order) result(expval)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(inout) :: hm
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: density(:)
#ifdef HAVE_MPI
  R_TYPE :: expval_tmp
#endif
  integer :: order

  PUSH_SUB(X(pert_expectation_value))

  order = 1
  if(present(pert_order)) order = pert_order

  ASSERT(order == 1 .or. order == 2)

  SAFE_ALLOCATE(density(1:gr%mesh%np))

  call X(pert_expectation_density)(this, gr, geo, hm, st, psia, psib, density, pert_order = order)

  expval = X(mf_integrate)(gr%mesh, density)

#ifdef HAVE_MPI
    ! reduce density
    if(st%parallel_in_states) then
        call MPI_Allreduce(expval, expval_tmp, 1, R_MPITYPE, MPI_SUM, st%mpi_grp%comm, mpi_err)
        expval = expval_tmp
    end if
#endif

  SAFE_DEALLOCATE_A(density)
  POP_SUB(X(pert_expectation_value))

end function X(pert_expectation_value)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
