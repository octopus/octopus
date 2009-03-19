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
subroutine X(pert_apply) (this, gr, geo, hm, ik, f_in, f_out)
! Returns f_out = H' f_in, where H' is perturbation Hamiltonian
! Note that e^ikr phase is applied to f_in, then is removed afterward
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)

  R_TYPE, allocatable :: f_in_copy(:)
   logical :: apply_kpoint

  call push_sub('pert_inc.Xpert_apply')

  call profiling_in(prof, "PERT_APPLY")

  ASSERT(this%dir .ne. -1)

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
     ALLOCATE(f_in_copy(1:NP_PART), NP_PART)
     call lalg_copy(NP_PART, f_in, f_in_copy)
     call X(set_bc(gr%der, f_in_copy(:)))
  endif
  ! no derivatives in electric, so ghost points not needed

  apply_kpoint = simul_box_is_periodic(gr%sb) .and. .not. kpoint_is_gamma(hm%d, ik) &
    .and. this%pert_type /= PERTURBATION_ELECTRIC
  ! electric does not need it since (e^-ikr)r(e^ikr) = r

  if (apply_kpoint) then
    f_in_copy(1:NP_PART) = hm%phase(1:NP_PART, ik) * f_in_copy(1:NP_PART)
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
    f_out(1:NP) = conjg(hm%phase(1:NP, ik)) * f_out(1:NP)
  endif

  if (this%pert_type /= PERTURBATION_ELECTRIC) then
     deallocate(f_in_copy)
  endif

  call profiling_out(prof)
  call pop_sub()

contains

  ! --------------------------------------------------------------------------

  subroutine none()

    f_out(1:NP) = M_ZERO

  end subroutine none

  ! --------------------------------------------------------------------------
  subroutine electric()

    call push_sub('pert_inc.X(pert_apply).electric')
    f_out(1:NP) = f_in(1:NP) * gr%mesh%x(1:NP, this%dir)
    call pop_sub()

  end subroutine electric

  ! --------------------------------------------------------------------------
  subroutine kdotp()
    R_TYPE, allocatable :: grad(:,:), cpsi(:,:)
    integer iatom

    call push_sub('pert_inc.X(pert_apply).kdotp')

    ALLOCATE(cpsi(1:NP_PART, hm%d%dim), NP_PART * hm%d%dim)
    ALLOCATE(grad(gr%mesh%np, gr%sb%dim), gr%mesh%np * gr%sb%dim)

    call X(derivatives_grad) (gr%der, f_in_copy, grad, set_bc = .false.)
    ! set_bc done already separately
    f_out(1:NP) = - M_zI * (grad(1:NP, this%dir))
    ! delta_H_k = (-i*grad + k) . delta_k
    ! representation on psi is just -i*grad . delta_k
    ! note that second-order term is left out

    if (this%use_nonlocalpps) then
      do iatom = 1, geo%natoms
        if(species_is_ps(geo%atom(iatom)%spec)) then
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, this%dir, ik, f_in_copy, cpsi(:, :))
          f_out(1:NP) = f_out(1:NP) - M_zI * cpsi(1:NP, 1)
          ! using only the first spinor component
        end if
      end do
    endif
 
    deallocate(grad, cpsi)
    call pop_sub()
    
  end subroutine kdotp

  ! --------------------------------------------------------------------------
  subroutine magnetic()
    R_TYPE, allocatable :: lf(:,:), vrnl(:,:,:)
    R_TYPE :: cross(1:MAX_DIM), vv(1:MAX_DIM)
    FLOAT :: xx(1:MAX_DIM)
    integer :: iatom, idir, ip

    call push_sub('pert_inc.X(pert_apply).magnetic')

    ALLOCATE(lf(gr%mesh%np, gr%sb%dim), gr%mesh%np * gr%sb%dim)

    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%mesh, gr%der, f_in_copy, lf, set_bc = .false.)
    f_out(1:NP) = M_HALF * lf(1:NP, this%dir)

    deallocate(lf)

    if(this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
      ALLOCATE(vrnl(NP, hm%d%dim, gr%sb%dim), NP * hm%d%dim * gr%sb%dim)
      vrnl(1:NP, 1:hm%d%dim, this%dir) = M_ZERO

      do iatom = 1, geo%natoms

        do idir = 1, gr%mesh%sb%dim
          if(this%dir == idir) cycle ! this direction is not used in the cross product
          call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir, ik, f_in_copy, vrnl(:, :, idir))
        end do

        xx(1:gr%mesh%sb%dim) = geo%atom(iatom)%x(1:gr%mesh%sb%dim)

        do ip = 1, NP

          if(this%gauge == GAUGE_ICL) xx(1:gr%mesh%sb%dim) = gr%mesh%x(ip, 1:gr%mesh%sb%dim)
         
          vv(1:gr%mesh%sb%dim) = vrnl(ip, 1, 1:gr%mesh%sb%dim)

          cross(1) = xx(2) * vv(3) - xx(3) * vv(2)
          cross(2) = xx(3) * vv(1) - xx(1) * vv(3)
          cross(3) = xx(1) * vv(2) - xx(2) * vv(1)

#if !defined(R_TCOMPLEX)
          f_out(ip) = f_out(ip) + M_HALF * cross(this%dir)
#else
          f_out(ip) = f_out(ip) - M_zI * M_HALF * cross(this%dir)
#endif
        end do
      end do

      deallocate(vrnl)
    end if

    call pop_sub()

  end subroutine magnetic

  ! --------------------------------------------------------------------------
  subroutine ionic
    integer :: iatom, idir
    R_TYPE, allocatable  :: tmp(:)

    call push_sub('pert_inc.Xpert_apply.ionic')
    ALLOCATE(tmp(1:NP), NP*1)
    
    f_out(1:NP) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim

        if (this%ionic%pure_dir .and. iatom /= this%atom1 .and. idir /= this%dir) cycle

        call X(ionic_perturbation)(this, gr, geo, hm, ik, f_in_copy, tmp, iatom, idir)
        
        call lalg_axpy(NP, this%ionic%mix1(iatom, idir), tmp, f_out)

      end do
    end do
      
    deallocate(tmp)
    call pop_sub()

  end subroutine ionic

end subroutine X(pert_apply)

  ! --------------------------------------------------------------------------
subroutine X(ionic_perturbation)(this, gr, geo, hm, ik, f_in, f_out, iatom, idir)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)
  integer,              intent(in)    :: iatom, idir    

  ! FIX ME: may need to tell derivatives_oper not to apply boundary conditions 
  ! more things about ghost points may need to be done

  R_TYPE, allocatable :: grad(:,:), fin(:, :), fout(:, :)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  call push_sub('pert_inc.Xionic_perturbation')

  ALLOCATE(vloc(1:NP), NP)
  vloc(1:NP) = M_ZERO
  call epot_local_potential(hm%ep, gr, gr%mesh, geo, iatom, vloc, CNST(0.0))

  ALLOCATE(fin(1:NP_PART, 1), NP_PART)
  call lalg_copy(NP_PART, f_in, fin(:, 1))

  !d^T v |f>
  ALLOCATE(fout(1:NP_PART, 1), NP_PART)
  forall(ip = 1:NP) fout(ip, 1) = vloc(ip)*fin(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, fin, fout, ik)
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, fout(:,1), f_out)

  !v d |f>
  ALLOCATE(grad(1:NP, 1), NP)
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, fin(:,1), grad(:,1))
  forall(ip = 1:NP) fout(ip, 1) = vloc(ip)*grad(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, grad, fout, ik)
  forall(ip = 1:NP) f_out(ip) = -f_out(ip) + fout(ip, 1)

  deallocate(grad, fin, fout, vloc)
  call pop_sub()

end subroutine X(ionic_perturbation)

! --------------------------------------------------------------------------
subroutine X(pert_apply_order_2) (this, gr, geo, hm, ik, f_in, f_out)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: hm
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)

! FIX ME: need to apply phases here

  call push_sub('pert_inc.Xpert_apply_order2')

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:NP) = R_TOTYPE(M_ZERO)
  case(PERTURBATION_IONIC)
    call ionic()
  case(PERTURBATION_MAGNETIC)
    call magnetic()
  end select

  call pop_sub()

contains

  ! --------------------------------------------------------------------------
  subroutine magnetic()
    R_TYPE, allocatable :: f_in2(:,:), dnl(:,:), vrnl(:,:), xf(:)
    R_TYPE :: cross1(1:3), bdir(1:MAX_DIM, 2)
    FLOAT  :: rdelta
    R_TYPE :: contr

    integer :: iatom, idir, idir2, ip

    call push_sub('pert_inc.Xpert_apply_order2.magnetic')

    do ip = 1, NP
      rdelta = sum(gr%mesh%x(ip, 1:MAX_DIM)**2) * ddelta(this%dir, this%dir2)
      f_out(ip) = M_FOURTH * (rdelta - gr%mesh%x(ip, this%dir)*gr%mesh%x(ip, this%dir2)) * f_in(ip)
    end do

    ! gauge correction
    apply_gauge: if(this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
      bdir(1:MAX_DIM, 1:2) = M_ZERO
      bdir(this%dir,  1)   = M_ONE
      bdir(this%dir2, 2)   = M_ONE

      ALLOCATE(f_in2(NP_PART, gr%mesh%sb%dim), NP_PART * gr%mesh%sb%dim)
      ALLOCATE(vrnl(NP_PART, hm%d%dim), NP_PART * hm%d%dim)
      ALLOCATE(dnl(NP, gr%mesh%sb%dim), NP * gr%mesh%sb%dim)
      ALLOCATE(xf(NP), NP)

      f_in2(NP:NP_PART,:) = R_TOTYPE(M_ZERO)
      atoms: do iatom = 1, geo%natoms

        ! This calculates f_in2 = (B x r) f_in
        do ip = 1, NP
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(gr%mesh%x(ip, :)))
          end select

          f_in2(ip, 1:gr%sb%dim) = cross1(1:gr%sb%dim) * f_in(ip)
        end do

        ! let us now get sum_beta Dnl f_in2
        dnl(1:NP, 1:gr%mesh%sb%dim) = R_TOTYPE(M_ZERO)
        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            !calculate dnl |f_in2> = -[x,vnl] |f_in2>
            call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir2, ik, f_in2(:, idir2), vrnl(:, :))

            ! -x vnl |f>
            dnl(1:NP, idir) = dnl(1:NP, idir) - gr%mesh%x(1:NP, idir) * vrnl(1:NP, 1)

            ! vnl x |f>
            xf(1:NP) = gr%mesh%x(1:NP, idir) * f_in2(1:NP, idir2)
            call X(projector_commute_r)(hm%ep%proj(iatom), gr, hm%d%dim, idir2, ik, xf, vrnl(:, :))

            dnl(1:NP, idir) = dnl(1:NP, idir) + vrnl(1:NP, 1)
          end do
        end do

        do ip = 1, NP
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(gr%mesh%x(ip, :)))
          end select

          contr = M_ZERO
          do idir = 1, gr%sb%dim
            contr = contr + cross1(idir) * dnl(ip, idir)
          end do
          f_out(ip) = f_out(ip) + M_FOURTH * contr
        end do

      end do atoms

      deallocate(f_in2, vrnl, dnl, xf)
    end if apply_gauge

    call pop_sub()

  end subroutine magnetic

  ! --------------------------------------------------------------------------
  subroutine ionic
    integer :: iatom, idir, jdir
    R_TYPE, allocatable  :: tmp(:)
    
    call push_sub('pert_inc.Xpert_apply_order2.ionic')

    ALLOCATE(tmp(1:NP), NP*1)
    
    f_out(1:NP) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        do jdir = 1, gr%mesh%sb%dim
          
          if (this%ionic%pure_dir &
               .and. iatom /= this%atom1 .and. idir /= this%dir &
               .and. iatom /= this%atom2 .and. jdir /= this%dir2) cycle

          call X(ionic_perturbation_order_2)(this, gr, geo, hm, ik, f_in, tmp, iatom, idir, jdir)
          
          call lalg_axpy(NP, this%ionic%mix1(iatom, idir)*this%ionic%mix2(iatom, jdir), tmp, f_out)
          
        end do
      end do
    end do

    deallocate(tmp)
    call pop_sub()

  end subroutine ionic

end subroutine X(pert_apply_order_2)

! --------------------------------------------------------------------------
subroutine X(ionic_perturbation_order_2) (this, gr, geo, hm, ik, f_in, f_out, iatom, idir, jdir)
  type(pert_t),        intent(in)    :: this
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: f_in(:)
  R_TYPE,              intent(out)   :: f_out(:)
  integer,             intent(in)    :: iatom, idir, jdir

  ! FIX ME: may need to tell derivatives_oper not to apply boundary conditions

  R_TYPE, allocatable :: fin(:, :)
  R_TYPE, allocatable :: tmp1(:, :), tmp2(:,:)
  FLOAT,  allocatable :: vloc(:)
  integer :: ip

  call push_sub('pert_inc.Xionic_perturbation_order2')

  ALLOCATE(fin(1:NP_PART, 1), NP_PART)
  ALLOCATE(tmp1(1:NP_PART, 1), NP_PART)
  ALLOCATE(tmp2(1:NP_PART, 1), NP_PART)
  ALLOCATE(vloc(1:NP), NP)

  forall(ip = 1:NP) vloc(ip) = M_ZERO
  call epot_local_potential(hm%ep, gr, gr%mesh, geo, iatom, vloc, CNST(0.0))

  call lalg_copy(NP_PART, f_in, fin(:, 1))
   
  !di^T dj^T v |f>
  forall(ip = 1:NP) tmp1(ip, 1) = vloc(ip)*fin(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, fin, tmp1, ik)
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, tmp1(:,1), tmp2(:,1))
  call X(derivatives_oper)(gr%der%grad(jdir), gr%der, tmp2(:,1), f_out)

  !di^T v dj |f>
  call X(derivatives_oper)(gr%der%grad(jdir), gr%der, fin(:,1), tmp1(:,1))
  forall(ip = 1:NP) tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, tmp2(:,1), tmp1(:,1))
  forall(ip = 1:NP) f_out(ip) = f_out(ip) - tmp1(ip, 1)

  !dj^T v di |f>
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  forall(ip = 1:NP) tmp2(ip, 1) = vloc(ip)*tmp1(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_oper)(gr%der%grad(jdir), gr%der, tmp2(:,1), tmp1(:,1))
  forall(ip = 1:NP) f_out(ip) = f_out(ip) - tmp1(ip, 1)

  !v di dj |f>
  call X(derivatives_oper)(gr%der%grad(idir), gr%der, fin(:,1), tmp1(:,1))
  call X(derivatives_oper)(gr%der%grad(jdir), gr%der, tmp1(:,1), tmp2(:,1))
  forall(ip = 1:NP) tmp1(ip, 1) = vloc(ip)*tmp2(ip, 1)
  call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, 1, tmp2, tmp1, ik)
  forall(ip = 1:NP) f_out(ip) = f_out(ip) + tmp1(ip, 1)

  call pop_sub()

end subroutine X(ionic_perturbation_order_2)

! --------------------------------------------------------------------------
subroutine X(ionic_pert_matrix_elements_2)(this, gr, geo, hm, ik, st, psi, vib, factor, matrix)
  type(pert_t),        intent(in)    :: this
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: hm
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

  ALLOCATE(vloc(1:NP), NP)
  ALLOCATE(gpsi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim), gr%mesh%np_part*st%d%dim*gr%sb%dim)
  ALLOCATE(g2psi(1:gr%mesh%np_part, 1:st%d%dim, 1:gr%sb%dim, 1:gr%sb%dim), gr%mesh%np_part*st%d%dim*gr%sb%dim**2)
  ALLOCATE(tmp1(1:gr%mesh%np, 1:st%d%dim), gr%mesh%np*st%d%dim)

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

      forall(ip = 1:NP) vloc(ip) = M_ZERO
      call epot_local_potential(hm%ep, gr, gr%mesh, geo, iatom, vloc, CNST(0.0))

      do jdir = 1, gr%sb%dim
        jmat = vibrations_get_index(vib, iatom, jdir)

        dot = M_ZERO

        !2<f|dj^T v di |f>
        forall (idim = 1:st%d%dim, ip = 1:NP) tmp1(ip, idim) = vloc(ip)*gpsi(ip, idim, idir)
        call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, st%d%dim, gpsi(:, :, idir), tmp1, ik)
        dot = dot + CNST(2.0)*X(mf_dotp)(gr%mesh, st%d%dim, gpsi(:, :, jdir), tmp1)

        !2<f|di^T dj^T v |f> 
        forall (idim = 1:st%d%dim, ip = 1:NP) tmp1(ip, idim) = vloc(ip)*psi(ip, idim, ist)
        call X(project_psi)(gr%mesh, hm%ep%proj(iatom:iatom), 1, st%d%dim, psi(:, :, ist), tmp1, ik)
        dot = dot + CNST(2.0)*X(mf_dotp)(gr%mesh, st%d%dim, g2psi(:, :, idir, jdir), tmp1)

        dot = dot*st%occ(ist, ik)*factor
        
        matrix(imat, jmat) = matrix(imat, jmat) + dot

      end do

    end do
  end do
  
end subroutine X(ionic_pert_matrix_elements_2)

! --------------------------------------------------------------------------
subroutine X(pert_expectation_density) (this, gr, geo, hm, st, psia, psib, density, pert_order)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: hm
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  R_TYPE,               intent(out)   :: density(:)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: pertpsib(:)
  integer :: ik, ist, idim, order
  FLOAT   :: ikweight

  call push_sub('pert_inc.Xpert_expectation_density')

  ALLOCATE(pertpsib(1:NP), NP)

  order = 1
  if(present(pert_order)) order = pert_order
  ASSERT(order == 1 .or. order == 2)

  density(1:NP) = R_TOTYPE(M_ZERO)

  do ik = 1, st%d%nik
    do ist  = st%st_start, st%st_end
      do idim = 1, st%d%dim

        if(order == 1) then 
          call X(pert_apply) (this, gr, geo, hm, ik, psib(:, idim, ist, ik), pertpsib)
          ikweight = st%d%kweights(ik)*st%smear%el_per_state
        else
          call X(pert_apply_order_2) (this, gr, geo, hm, ik, psib(:, idim, ist, ik), pertpsib)
          ikweight = st%d%kweights(ik)*st%occ(ist, ik)
        end if

        density(1:NP) = density(1:NP) + ikweight * &
             R_CONJ(psia(1:NP, idim, ist, ik))*pertpsib(1:NP)

      end do
    end do
  end do

  deallocate(pertpsib)
  call pop_sub()

end subroutine X(pert_expectation_density)

! --------------------------------------------------------------------------
R_TYPE function X(pert_expectation_value) (this, gr, geo, hm, st, psia, psib, pert_order) result(expval)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: hm
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: density(:)
#ifdef HAVE_MPI
  R_TYPE :: expval_tmp
#endif
  integer :: order

  call push_sub('pert_inc.Xpert_expectation_value')

  order = 1
  if(present(pert_order)) order = pert_order

  ASSERT(order == 1 .or. order == 2)

  ALLOCATE(density(1:NP), NP)

  call X(pert_expectation_density)(this, gr, geo, hm, st, psia, psib, density, pert_order = order)

  expval = X(mf_integrate)(gr%mesh, density)

#ifdef HAVE_MPI
    ! reduce density
    if(st%parallel_in_states) then
        call MPI_Allreduce(expval, expval_tmp, 1, R_MPITYPE, MPI_SUM, st%mpi_grp%comm, mpi_err)
        expval = expval_tmp
    end if
#endif

  deallocate(density)
  call pop_sub()

end function X(pert_expectation_value)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
