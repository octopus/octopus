!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: em_resp.F90 2548 2006-11-06 21:42:27Z xavier $


! ---------------------------------------------------------
! \warning: This subroutine is clearly broken after the changes
! to include temperature in linear response
subroutine X(lr_calc_elf)(st, gr, lr, lr_m)
  type(states_t),       intent(inout) :: st
  type(grid_t),         intent(inout) :: gr
  type(lr_t),           intent(inout) :: lr
  type(lr_t), optional, intent(inout) :: lr_m !when this argument is present, we are doing dynamical response

  integer :: ip, idir, is, ist, idim, ik

  R_TYPE, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)
  FLOAT,  allocatable :: rho(:), grho(:,:)
  R_TYPE, allocatable :: dl_rho(:), gdl_rho(:,:)
  FLOAT,  allocatable :: elf(:,:), de(:,:), current(:, :, :)
  FLOAT :: dl_d0, d0
  FLOAT :: factor, spin_factor

  FLOAT, parameter :: dmin = CNST(1e-10)
  FLOAT :: ik_weight

  PUSH_SUB(X(lr_calc_elf))

  ASSERT(.false.)

  SAFE_ALLOCATE(   gpsi(1:gr%mesh%np, 1:gr%mesh%sb%dim))
  SAFE_ALLOCATE(gdl_psi(1:gr%mesh%np, 1:gr%mesh%sb%dim))

  if(present(lr_m)) then
    SAFE_ALLOCATE(gdl_psi_m(1:gr%mesh%np, 1:gr%mesh%sb%dim))
  end if

  SAFE_ALLOCATE(   grho(1:gr%mesh%np, 1:gr%mesh%sb%dim))
  SAFE_ALLOCATE(gdl_rho(1:gr%mesh%np, 1:gr%mesh%sb%dim))

  SAFE_ALLOCATE(   rho(1:gr%mesh%np_part))
  SAFE_ALLOCATE(dl_rho(1:gr%mesh%np_part))

  SAFE_ALLOCATE(   elf(1:gr%mesh%np, 1:st%d%nspin))
  SAFE_ALLOCATE(    de(1:gr%mesh%np, 1:st%d%nspin))
  SAFE_ALLOCATE(current(1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:st%d%nspin))

  if( .not. associated(lr%X(dl_de)) ) then
    SAFE_ALLOCATE(lr%X(dl_de)(1:gr%mesh%np, 1:st%d%nspin)) 
  end if
  if( .not. associated(lr%X(dl_elf))) then
    SAFE_ALLOCATE(lr%X(dl_elf)(1:gr%mesh%np, 1:st%d%nspin))
  end if

  !calculate the gs elf
  call elf_calc(st, gr, elf, de)

  !calculate current and its variation
  if(states_are_complex(st)) then 
    call calc_physical_current(gr%der, st, current)
    if(present(lr_m)) then 
      call lr_calc_current(st, gr, lr, lr_m)
    else
      call lr_calc_current(st, gr, lr)
    end if
  end if

  ! single or double occupancy
  if(st%d%nspin == 1) then
    spin_factor = M_TWO
  else
    spin_factor = M_ONE
  end if

  lr%X(dl_de) = M_ZERO

  do is = 1, st%d%nspin
    rho  = M_ZERO
    grho = M_ZERO

    dl_rho  = M_ZERO
    gdl_rho = M_ZERO

    !first we calculate the densities and their gradients; this could
    !be done directly, but it is less precise numerically
    do ik = is, st%d%nik, st%d%nspin
      do ist = 1, st%nst
        ik_weight = st%d%kweights(ik) * st%occ(ist, ik) / spin_factor

        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          ! sum over states to obtain the spin-density
          rho(1:gr%mesh%np) = rho(1:gr%mesh%np) + ik_weight * abs(st%X(psi)(1:gr%mesh%np, idim, ist, is))**2

          !the gradient of the density
          do idir = 1, gr%mesh%sb%dim
            grho(1:gr%mesh%np, idir) = grho(1:gr%mesh%np, idir) + &
                 ik_weight * M_TWO * R_REAL(R_CONJ(st%X(psi)(1:gr%mesh%np, idim, ist, is)) * gpsi(1:gr%mesh%np, idir))
          end do

          !the variation of the density and its gradient

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            dl_rho(1:gr%mesh%np) = dl_rho(1:gr%mesh%np) + ik_weight * ( &
                 R_CONJ(st%X(psi)(1:gr%mesh%np, idim, ist, is)) * lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is)+ & 
                 st%X(psi)(1:gr%mesh%np, idim, ist, is) * R_CONJ(lr_m%X(dl_psi)(1:gr%mesh%np, idim, ist, is)) )

            do idir = 1, gr%mesh%sb%dim

              gdl_rho(1:gr%mesh%np, idir) = gdl_rho(1:gr%mesh%np, idir) + ik_weight * ( &
                   R_CONJ(st%X(psi)(1:gr%mesh%np, idim, ist, is)) * gdl_psi(1:gr%mesh%np, idir) +      &
                   R_CONJ(gpsi(1:gr%mesh%np, idir))* lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is)  +      &
                   st%X(psi)(1:gr%mesh%np, idim, ist, is) * R_CONJ(gdl_psi_m(1:gr%mesh%np, idir)) +      &
                   gpsi(1:gr%mesh%np, idir) * R_CONJ(lr_m%X(dl_psi)(1:gr%mesh%np, idim, ist, is))  )

            end do

          else
            
            dl_rho(1:gr%mesh%np) = dl_rho(1:gr%mesh%np) + ik_weight*M_TWO* &
                 R_REAL(R_CONJ(st%X(psi)(1:gr%mesh%np, idim, ist, is))*lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is))

            do idir = 1, gr%mesh%sb%dim
              gdl_rho(1:gr%mesh%np, idir) = gdl_rho(1:gr%mesh%np, idir) + ik_weight * M_TWO * ( &
                   R_CONJ(st%X(psi)(1:gr%mesh%np, idim, ist, is)) * gdl_psi(1:gr%mesh%np, idir) + &
                   gpsi(1:gr%mesh%np, idir)*R_CONJ(lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is))  )
            end do

          end if

        end do !idim
      end do !ist
    end do !ik

    !now we start to calculate the ELF

    !first the term that depends on the orbitals
    !this is the only term that is different for the dynamical case
    do ik = is, st%d%nik, st%d%nspin
      do ist = 1, st%nst
        ik_weight = st%d%kweights(ik) * st%occ(ist, ik) / spin_factor
        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            do ip = 1, gr%mesh%np
              lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) + &
                dl_rho(ip) * ik_weight * sum(R_ABS(gpsi(ip, 1:gr%mesh%sb%dim))**2) + &
                rho(ip) * ik_weight * &
                sum(R_CONJ(gpsi(ip, 1:gr%mesh%sb%dim)) * gdl_psi(ip, 1:gr%mesh%sb%dim) + &
                gpsi(ip, 1:gr%mesh%sb%dim) * R_CONJ(gdl_psi_m(ip, 1:gr%mesh%sb%dim)))
            end do

          else 

            do ip = 1, gr%mesh%np
              lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) + &
                dl_rho(ip) * ik_weight * sum(R_ABS(gpsi(ip, 1:gr%mesh%sb%dim))**2) + &
                rho(ip) * ik_weight * M_TWO * &
                (sum(R_CONJ(gpsi(ip, 1:gr%mesh%sb%dim)) * gdl_psi(ip, 1:gr%mesh%sb%dim)))
            end do

          end if

        end do
      end do
    end do

    !the density term
    do ip = 1, gr%mesh%np
      if(abs(st%rho(ip, is)) >= dmin) then
        lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) - &
          M_HALF * sum(grho(ip, 1:gr%mesh%sb%dim) * gdl_rho(ip, 1:gr%mesh%sb%dim))
      end if
    end do

    !the current term
    if(states_are_complex(st)) then       
      do ip = 1, gr%mesh%np
        if(abs(st%rho(ip, is)) >= dmin) then
          lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) + &
            M_TWO * sum(current(ip, 1:gr%mesh%sb%dim, is) * lr%dl_j(ip, 1:gr%mesh%sb%dim, is))
        end if
      end do
    end if

    !now the normalization 
    factor = M_THREE/M_FIVE * (M_SIX * M_PI**2)**M_TWOTHIRD
    do ip = 1, gr%mesh%np

      if(abs(st%rho(ip, is)) >= dmin) then
        d0    = factor * rho(ip)**(M_EIGHT / M_THREE)
        dl_d0 = M_EIGHT/M_THREE * factor * dl_rho(ip) * rho(ip)**(M_FIVE/M_THREE)

        lr%X(dl_elf)(ip, is) = &
             M_TWO * d0 * dl_d0/(d0**2 + de(ip, is)**2) * (M_ONE - elf(ip,is)) &
             - M_TWO * de(ip, is) * lr%X(dl_de)(ip, is) / (d0**2 + de(ip, is)**2) * elf(ip, is)
      else
        lr%X(dl_elf)(ip, is) = M_ZERO
      end if

    end do

  end do

  SAFE_DEALLOCATE_A(current)
  SAFE_DEALLOCATE_A(gpsi)
  SAFE_DEALLOCATE_A(gdl_psi)
  if(present(lr_m)) then
    SAFE_DEALLOCATE_A(gdl_psi_m)
  end if

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(dl_rho)

  SAFE_DEALLOCATE_A(grho)
  SAFE_DEALLOCATE_A(gdl_rho)

  SAFE_DEALLOCATE_A(elf)
  SAFE_DEALLOCATE_A(de)

  POP_SUB(X(lr_calc_elf))

end subroutine X(lr_calc_elf)


! ---------------------------------------------------------
subroutine X(calc_polarizability_periodic)(sys, em_lr, kdotp_lr, nsigma, zpol, ndir)
  type(system_t),         intent(inout) :: sys
  type(lr_t),             intent(inout) :: em_lr(:,:)
  type(lr_t),             intent(inout) :: kdotp_lr(:)
  integer,                intent(in)    :: nsigma
  CMPLX,                  intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)
  integer, optional,      intent(in)    :: ndir

  integer :: dir1, dir2, ndir_, ist, ik
  CMPLX :: term, subterm
  type(mesh_t), pointer :: mesh
#ifdef HAVE_MPI
  CMPLX :: zpol_temp(1:MAX_DIM, 1:MAX_DIM)
#endif

  PUSH_SUB(X(calc_polarizability_periodic))

  mesh => sys%gr%mesh

  ndir_ = mesh%sb%periodic_dim
  if(present(ndir)) ndir_ = ndir

  ! alpha_ij(w) = -e sum(m occ, k) [(<u_mk(0)|-id/dk_i)|u_mkj(1)(w)> + <u_mkj(1)(-w)|(-id/dk_i|u_mk(0)>)]

  do dir1 = 1, ndir_
    do dir2 = 1, sys%gr%sb%dim

      zpol(dir1, dir2) = M_ZERO

      do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
        do ist = sys%st%st_start, sys%st%st_end
          term = M_ZERO
        
          subterm = - X(mf_dotp)(mesh, sys%st%d%dim, &
            kdotp_lr(dir1)%X(dl_psi)(1:mesh%np, 1:sys%st%d%dim, ist, ik), &
            em_lr(dir2, 1)%X(dl_psi)(1:mesh%np, 1:sys%st%d%dim, ist, ik))
          term = term + subterm

          if(nsigma == 1) then
            term = term + conjg(subterm)
          else
            term = term - X(mf_dotp)(mesh, sys%st%d%dim, &
              em_lr(dir2, 2)%X(dl_psi)(1:mesh%np, 1:sys%st%d%dim, ist, ik), & 
              kdotp_lr(dir1)%X(dl_psi)(1:mesh%np, 1:sys%st%d%dim, ist, ik))
          end if

          zpol(dir1, dir2) = zpol(dir1, dir2) + &
            term * sys%st%d%kweights(ik) * sys%st%occ(ist, ik)
        enddo

      enddo
    enddo
  enddo

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    zpol(1:mesh%sb%periodic_dim, 1:mesh%sb%dim) = zpol_temp(1:mesh%sb%periodic_dim, 1:mesh%sb%dim)
  endif
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    zpol(1:mesh%sb%periodic_dim, 1:mesh%sb%dim) = zpol_temp(1:mesh%sb%periodic_dim, 1:mesh%sb%dim)
  endif
#endif

  POP_SUB(X(calc_polarizability_periodic))

end subroutine X(calc_polarizability_periodic)


! ---------------------------------------------------------
subroutine X(calc_polarizability_finite)(sys, hm, lr, nsigma, perturbation, zpol, doalldirs, ndir)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr(:,:)
  integer,                intent(in)    :: nsigma
  type(pert_t),           intent(inout) :: perturbation
  CMPLX,                  intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)
  logical, optional,      intent(in)    :: doalldirs
  integer, optional,      intent(in)    :: ndir

  integer :: dir1, dir2, ndir_, startdir

  PUSH_SUB(X(calc_polarizability_finite))

  ndir_ = sys%gr%sb%dim
  if(present(ndir)) ndir_ = ndir

  startdir = sys%gr%sb%periodic_dim + 1
  if(present(doalldirs)) then
    if(doalldirs) startdir = 1
  endif

  ! alpha_ij(w) = - sum(m occ) [<psi_m(0)|r_i|psi_mj(1)(w)> + <psi_mj(1)(-w)|r_i|psi_m(0)>]
  ! minus sign is from electronic charge -e

  do dir1 = startdir, ndir_
    do dir2 = 1, sys%gr%sb%dim
      call pert_setup_dir(perturbation, dir1)
      zpol(dir1, dir2) = -X(pert_expectation_value)(perturbation, sys%gr, sys%geo, hm, &
        sys%st, sys%st%X(psi), lr(dir2, 1)%X(dl_psi))

      if(nsigma == 1) then
        zpol(dir1, dir2) = zpol(dir1, dir2) + R_CONJ(zpol(dir1, dir2))
      else
        zpol(dir1, dir2) = zpol(dir1, dir2) &
          - X(pert_expectation_value)(perturbation, sys%gr, sys%geo, hm, &
          sys%st, lr(dir2, 2)%X(dl_psi), sys%st%X(psi))
      end if
    end do
  end do

  POP_SUB(X(calc_polarizability_finite))

end subroutine X(calc_polarizability_finite)


! ---------------------------------------------------------
subroutine X(lr_calc_susceptibility)(sys, hm, lr, nsigma, perturbation, chi_para, chi_dia)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr(:,:)
  integer,                intent(in)    :: nsigma
  type(pert_t),           intent(inout) :: perturbation
  CMPLX,                  intent(out)   :: chi_para(:,:), chi_dia(:,:)

  integer :: dir1, dir2
  R_TYPE  :: trace

  PUSH_SUB(X(lr_calc_susceptibility))

  chi_para = M_ZERO
  chi_dia  = M_ZERO

  do dir1 = 1, sys%gr%sb%dim
    do dir2 = 1, sys%gr%sb%dim

      call pert_setup_dir(perturbation, dir1, dir2)

      trace = X(pert_expectation_value)(perturbation, sys%gr, sys%geo, hm, sys%st, sys%st%X(psi), lr(dir2, 1)%X(dl_psi))
      
      if (nsigma == 1) then 
        trace = trace + R_CONJ(trace)
      else
        trace = trace + &
          X(pert_expectation_value)(perturbation, sys%gr, sys%geo, hm, sys%st, lr(dir2, 2)%X(dl_psi), sys%st%X(psi))
      end if
     
      ! first the paramagnetic term 
      chi_para(dir1, dir2) = chi_para(dir1, dir2) + trace

      chi_dia(dir1, dir2) = chi_dia(dir1, dir2) + &
        X(pert_expectation_value)(perturbation, sys%gr, sys%geo, hm, sys%st, sys%st%X(psi), sys%st%X(psi), pert_order=2)

    end do
  end do


  ! We now add the minus sign from the definition of the susceptibility (chi = -d^2 E /d B^2)
  chi_para(:,:) = -chi_para(:,:)/P_C**2
  chi_dia (:,:) = -chi_dia (:,:)/P_C**2

#if defined(R_TREAL)
  ! When using real wavefunctions there is an extra factor of (-i)*(-i) = -1
  chi_para(:,:) = -chi_para(:,:)
#endif

  POP_SUB(X(lr_calc_susceptibility))

end subroutine X(lr_calc_susceptibility)

! ---------------------------------------------------------
subroutine X(lr_calc_beta) (sh, sys, hm, em_lr, dipole, beta, kdotp_lr, kdotp_em_lr, occ_response)
! See (16) in X Andrade et al., J. Chem. Phys. 126, 184106 (2006) for finite systems
! and (10) in A Dal Corso et al., Phys. Rev. B 15, 15638 (1996) for periodic systems
! Supply only em_lr for finite systems, and both kdotp_lr and kdotp_em_lr for periodic
! em_lr(dir, sigma, omega) = electric perturbation of ground-state wavefunctions
! kdotp_lr(dir) = kdotp perturbation of ground-state wavefunctions
! kdotp_em_lr(dir1, dir2, sigma, omega) = kdotp perturbation of electric-perturbed wfns
  type(sternheimer_t),     intent(inout) :: sh
  type(system_t), target,  intent(inout) :: sys
  type(hamiltonian_t),     intent(inout) :: hm
  type(lr_t),              intent(inout) :: em_lr(:,:,:)
  type(pert_t),            intent(inout) :: dipole
  CMPLX,                   intent(out)   :: beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM)
  type(lr_t),    optional, intent(in)    :: kdotp_lr(:)
  type(lr_t),    optional, intent(in)    :: kdotp_em_lr(:,:,:,:) ! kdotp dir, em dir, sigma, factor
  logical,       optional, intent(in)    :: occ_response ! do the wfns include the occ subspace?

  ! occ_response = yes is based on Baroni et al., RMP 73, 515 (2001), eqn 122
  ! occ_response = no  is based on Baroni et al., RMP 73, 515 (2001), eqn 123
  !   The occ_response = no version can be used even if the wfns do include the
  !   occupied subspace, it is just more efficient to use the other formula.

  type(states_t), pointer :: st
  type(mesh_t),   pointer :: mesh

  integer :: ifreq, jfreq, isigma, idim, ispin, np, ndir, idir, ist, ik
  integer :: ii, jj, kk, iperm, op_sigma, ist2, ip, is1, is2, is3
  integer :: perm(1:3), u(1:3), w(1:3), ijk(1:3)

  R_TYPE, allocatable :: hvar(:, :, :, :, :, :)
  R_TYPE, allocatable :: tmp(:, :), ppsi(:, :), me010(:, :, :, :, :)
  type(matrix_t), allocatable :: me11(:, :, :, :, :, :)
  FLOAT,  allocatable :: rho(:,:), kxc(:, :, :, :)
  R_TYPE, allocatable :: hpol_density(:)
  logical :: occ_response_

  PUSH_SUB(X(lr_calc_beta))

  call profiling_in(beta_prof, "CALC_BETA")

  st   => sys%st
  mesh => sys%gr%mesh
  np   =  mesh%np
  ndir =  sys%gr%sb%dim

  occ_response_ = .false.
  if(present(occ_response)) occ_response_ = occ_response

  ASSERT(present(kdotp_lr) .eqv. present(kdotp_em_lr))
  ! either both are absent for finite, or both present for periodic

  if(sternheimer_add_fxc(sh)) then
    !calculate kxc, the derivative of fxc
    SAFE_ALLOCATE(kxc(1:np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
    kxc = M_ZERO
    
    SAFE_ALLOCATE(rho(1:np, 1:st%d%nspin))
    call states_total_density(st, mesh, rho)
    
    call xc_get_kxc(sys%ks%xc, mesh, rho, st%d%ispin, kxc)
    SAFE_DEALLOCATE_A(rho)
    SAFE_ALLOCATE(hpol_density(1:np))
  endif

  SAFE_ALLOCATE(tmp(1:np, 1:st%d%dim))
  SAFE_ALLOCATE(ppsi(1:np, 1:st%d%dim))
  SAFE_ALLOCATE(hvar(1:np, 1:st%d%nspin, 1:2, 1:st%d%dim, 1:ndir, 1:3))
  SAFE_ALLOCATE(me010(1:st%nst, 1:st%nst, 1:mesh%sb%dim, 1:3, st%d%kpt%start:st%d%kpt%end))
  SAFE_ALLOCATE(me11(1:mesh%sb%dim, 1:mesh%sb%dim, 1:3, 1:3, 1:2, st%d%kpt%start:st%d%kpt%end))

  do ifreq = 1, 3
    do idir = 1, ndir
      do idim = 1, st%d%dim
        call X(sternheimer_calc_hvar)(sh, sys, hm, em_lr(idir, :, ifreq), 2, hvar(:, :, :, idim, idir, ifreq))
      end do !idim
    end do !idir
  end do !ifreq

  call get_matrix_elements()

  beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM) = M_ZERO

  do ii = 1, ndir
    do jj = 1, ndir
      do kk = 1, ndir

        ijk = (/ ii, jj, kk /)

        do isigma = 1, 2

          op_sigma = 2 
          if(isigma == 2) op_sigma = 1

          do iperm = 1, 6

            call get_permutation(iperm, perm)

            u(1:3) = ijk(perm(1:3))
            w(1:3) = perm(1:3)

            do ik = st%d%kpt%start, st%d%kpt%end
              do ist = 1, st%nst
                do idim = 1, st%d%dim

                  ispin = states_dim_get_spin_index(sys%st%d, ik)

                  if (present(kdotp_em_lr) .and. u(2) <= sys%gr%sb%periodic_dim) then
                    tmp(1:np, 1) = - kdotp_em_lr(u(2), u(3), isigma, w(3))%X(dl_psi)(1:np, idim, ist, ik)
                  else
                    call pert_setup_dir(dipole, u(2))
                    call X(pert_apply) &
                      (dipole, sys%gr, sys%geo, hm, ik, em_lr(u(3), isigma, w(3))%X(dl_psi)(:, :, ist, ik), tmp)
                  endif

                  do ip = 1, np
                    tmp(ip, 1) = tmp(ip, 1) + hvar(ip, ispin, isigma, idim, u(2), w(2)) &
                      * em_lr(u(3), isigma, w(3))%X(dl_psi)(ip, idim, ist, ik)
                  enddo

                  beta(ii, jj, kk) = beta(ii, jj, kk) &
                    - M_HALF * st%d%kweights(ik) * st%smear%el_per_state &
                    * X(mf_dotp)(mesh, em_lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, idim, ist, ik), tmp(1:np, 1))
                  
                  do ist2 = 1, st%nst
                    if(occ_response_ .and. ist2 .ne. ist) cycle
                    beta(ii, jj, kk) = beta(ii, jj, kk) + & 
                      M_HALF * st%d%kweights(ik)*st%smear%el_per_state*me010(ist2, ist, u(2), w(2), ik)*&
                      me11(u(1), u(3), w(1), w(3), isigma, ik)%X(matrix)(ist, ist2)
                  end do ! ist2

                end do !idim
              end do !ist
            end do !ik

            if(sternheimer_add_fxc(sh)) then 
              
              hpol_density(:) = M_ZERO
              do ip = 1, np
                do is1 = 1, st%d%nspin
                  do is2 = 1, st%d%nspin
                    do is3 = 1, st%d%nspin
                      hpol_density(ip) = hpol_density(ip) + kxc(ip, is1, is2, is3) &
                        * em_lr(u(1), isigma, w(1))%X(dl_rho)(ip, is1) & 
                        * em_lr(u(2), isigma, w(2))%X(dl_rho)(ip, is2) & 
                        * em_lr(u(3), isigma, w(3))%X(dl_rho)(ip, is3)
                    enddo
                  enddo
                enddo
              end do

              beta(ii, jj, kk) = beta(ii, jj, kk) - M_HALF * X(mf_integrate)(mesh, hpol_density / CNST(6.0))

            end if

          end do ! iperm

        end do !isigma

      end do !kk
    end do !jj
  end do !ii

  do ik = st%d%kpt%start, st%d%kpt%end
    do ii = 1, ndir
      do jj = 1, ndir
        do ifreq = 1, 3
          do jfreq = 1, 3
            do isigma = 1,2
              SAFE_DEALLOCATE_P(me11(ii, jj, ifreq, jfreq, isigma, ik)%X(matrix))
            end do
          end do
        end do
      end do
    end do
  end do

  if(sternheimer_add_fxc(sh)) then
    SAFE_DEALLOCATE_A(hpol_density)
    SAFE_DEALLOCATE_A(kxc)
  endif
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(hvar)
  SAFE_DEALLOCATE_A(me010)
  SAFE_DEALLOCATE_A(me11)
  SAFE_DEALLOCATE_A(hvar)

  POP_SUB(X(lr_calc_beta))

  call profiling_out(beta_prof)

contains

  subroutine get_permutation(ii, pp)
    integer, intent(in)  :: ii
    integer, intent(out) :: pp(1:3)

    PUSH_SUB(X(lr_calc_beta).get_permutation)

    ASSERT( ii >= 1 .and. ii <= 6)

    select case(ii)

    case(1) ; pp(1)=1 ; pp(2)=2 ; pp(3)=3
    case(2) ; pp(1)=2 ; pp(2)=3 ; pp(3)=1
    case(3) ; pp(1)=3 ; pp(2)=1 ; pp(3)=2
    case(4) ; pp(1)=3 ; pp(2)=2 ; pp(3)=1
    case(5) ; pp(1)=1 ; pp(2)=3 ; pp(3)=2
    case(6) ; pp(1)=2 ; pp(2)=1 ; pp(3)=3

    end select

    POP_SUB(X(lr_calc_beta).get_permutation)
  end subroutine get_permutation

  subroutine get_matrix_elements()
    PUSH_SUB(X(lr_calc_beta).get_matrix_elements)

    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_dim_get_spin_index(sys%st%d, ik)
      do ist = 1, st%nst
        do ist2 = 1, st%nst

          if(occ_response_ .and. ist2 .ne. ist) cycle
          do ii = 1, ndir
            call pert_setup_dir(dipole, ii)

            do ifreq = 1, 3

              if (present(kdotp_lr) .and. ii <= sys%gr%sb%periodic_dim) then
                forall (idim = 1:st%d%dim, ip = 1:np) ppsi(ip, idim) = -kdotp_lr(ii)%X(dl_psi)(ip, idim, ist, ik)
              else
                call X(pert_apply)(dipole, sys%gr, sys%geo, hm, ik, st%X(psi)(:, :, ist, ik), ppsi)
              endif

              isigma = 1
              forall (idim = 1:st%d%dim, ip = 1:np)
                ppsi(ip, idim) = ppsi(ip, idim) + hvar(ip, ispin, isigma, idim, ii, ifreq)*st%X(psi)(ip, idim, ist, ik)
              end forall

              me010(ist2, ist, ii, ifreq, ik) = X(mf_dotp)(mesh, st%d%dim, st%X(psi)(:, :, ist2, ik), ppsi)

            end do
          end do
        end do
      end do
    end do

    do ik = st%d%kpt%start, st%d%kpt%end
      do ii = 1, ndir
        do jj = 1, ndir
          do ifreq = 1, 3
            do jfreq = 1, 3
              do isigma = 1, 2
                op_sigma = 2 
                if(isigma == 2) op_sigma = 1
                SAFE_ALLOCATE(me11(ii, jj, ifreq, jfreq, isigma, ik)%X(matrix)(1:st%nst, 1:st%nst))

                if(occ_response_) then
                  do ist = 1, st%nst
                    me11(ii, jj, ifreq, jfreq, isigma, ik)%X(matrix)(ist, ist) = &
                      X(mf_dotp)(mesh, st%d%dim, em_lr(ii, op_sigma, ifreq)%X(dl_psi)(:, :, ist, ik), &
                                                 em_lr(jj, isigma,   jfreq)%X(dl_psi)(:, :, ist, ik))
                  enddo
                else
                  call states_blockt_mul(mesh, st, st%st_start, st%st_end, st%st_start, st%st_end, &
                    em_lr(ii, op_sigma, ifreq)%X(dl_psi)(:, :, :, ik), &
                    em_lr(jj, isigma, jfreq)%X(dl_psi)(:, :, :, ik), &
                    me11(ii, jj, ifreq, jfreq, isigma, ik)%X(matrix))
                endif

              end do
            end do
          end do
        end do
      end do
    end do

    POP_SUB(X(lr_calc_beta).get_matrix_elements)
  end subroutine get_matrix_elements

end subroutine X(lr_calc_beta)

! ---------------------------------------------------------
subroutine X(lr_calc_2np1) (sh, hm, st, geo, gr, lr1, lr2, lr3, pert1, pert2, pert3, val)
  type(sternheimer_t),     intent(inout) :: sh
  type(hamiltonian_t),     intent(inout) :: hm
  type(states_t),          intent(in)    :: st
  type(geometry_t),        intent(in)    :: geo
  type(grid_t),            intent(inout) :: gr
  type(lr_t),              intent(in)    :: lr1
  type(lr_t),              intent(in)    :: lr2
  type(lr_t),              intent(in)    :: lr3
  type(pert_t),            intent(in)    :: pert1
  type(pert_t),            intent(in)    :: pert2
  type(pert_t),            intent(in)    :: pert3
  R_TYPE,                  intent(out)   :: val

  integer :: ik, ist, jst
  R_TYPE :: term
  R_TYPE, allocatable :: tmp(:, :), me23(:, :)

  PUSH_SUB(X(lr_calc_2np1))

  SAFE_ALLOCATE(tmp(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(me23(1:st%nst, 1:st%nst))

  val = R_TOTYPE(M_ZERO)

  do ik = st%d%kpt%start, st%d%kpt%end

    !precalculate these matrix elements (it could be done in a better way)
    do ist = 1, st%nst
      do jst = 1, st%nst
        me23(ist, jst) = X(mf_dotp)(gr%mesh, st%d%dim, lr2%X(dl_psi)(:, :, ist, ik), lr3%X(dl_psi)(:, :, jst, ik))
      end do
    end do

    do ist = 1, st%nst
      
      call X(pert_apply)(pert2, gr, geo, hm, ik, lr3%X(dl_psi)(:, :, ist, ik), tmp)
      term = X(mf_dotp)(gr%mesh, st%d%dim, lr1%X(dl_psi)(:, :, ist, ik), tmp)
      
      do jst = 1, st%nst
        call X(pert_apply)(pert1, gr, geo, hm, ik, st%X(psi)(:, :, jst, ik), tmp)
        term = term - X(mf_dotp)(gr%mesh, st%d%dim, st%X(psi)(:, :, jst, ik), tmp)*me23(jst, ist)
      end do

      val = val + st%d%kweights(ik)*st%smear%el_per_state*term
      
    end do
  end do
  
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(me23)

  POP_SUB(X(lr_calc_2np1))
end subroutine X(lr_calc_2np1)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
