!! Copyright (C) 2004-2012 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca), David Strubbe
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


! ---------------------------------------------------------
! \warning: This subroutine is clearly broken after the changes
! to include temperature in linear response
subroutine X(lr_calc_elf)(st, gr, lr, lr_m)
  type(states_elec_t),  intent(inout) :: st
  type(grid_t),         intent(inout) :: gr
  type(lr_t),           intent(inout) :: lr
  type(lr_t), optional, intent(inout) :: lr_m !< when this argument is present, we are doing dynamical response

  integer :: ip, idir, is, ist, idim, ik

  R_TYPE, allocatable :: psi(:, :), gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)
  FLOAT,  allocatable :: rho(:), grho(:,:)
  R_TYPE, allocatable :: dl_rho(:), gdl_rho(:,:)
  FLOAT,  allocatable :: elf(:,:), de(:,:), current(:, :, :)
  R_TYPE :: dl_d0
  FLOAT :: d0, factor, spin_factor

  FLOAT, parameter :: dmin = CNST(1e-10)
  FLOAT :: ik_weight

  PUSH_SUB(X(lr_calc_elf))

  ASSERT(.false.)

  SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1:gr%mesh%sb%dim))
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

        call states_elec_get_state(st, gr%mesh, ist, is, psi)
        
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/spin_factor

        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, psi(:, idim), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          ! sum over states to obtain the spin-density
          rho(1:gr%mesh%np) = rho(1:gr%mesh%np) + ik_weight * abs(psi(1:gr%mesh%np, idim))**2

          !the gradient of the density
          do idir = 1, gr%mesh%sb%dim
            grho(1:gr%mesh%np, idir) = grho(1:gr%mesh%np, idir) + &
                 ik_weight*M_TWO*R_REAL(R_CONJ(psi(1:gr%mesh%np, idim))*gpsi(1:gr%mesh%np, idir))
          end do

          !the variation of the density and its gradient

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            dl_rho(1:gr%mesh%np) = dl_rho(1:gr%mesh%np) + ik_weight*( &
              R_CONJ(psi(1:gr%mesh%np, idim))*lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is) + &
              psi(1:gr%mesh%np, idim)*R_CONJ(lr_m%X(dl_psi)(1:gr%mesh%np, idim, ist, is)))

            do idir = 1, gr%mesh%sb%dim

              gdl_rho(1:gr%mesh%np, idir) = gdl_rho(1:gr%mesh%np, idir) + ik_weight*( &
                   R_CONJ(psi(1:gr%mesh%np, idim))*gdl_psi(1:gr%mesh%np, idir) + &
                   R_CONJ(gpsi(1:gr%mesh%np, idir))* lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is) + &
                   psi(1:gr%mesh%np, idim)*R_CONJ(gdl_psi_m(1:gr%mesh%np, idir)) + &
                   gpsi(1:gr%mesh%np, idir)*R_CONJ(lr_m%X(dl_psi)(1:gr%mesh%np, idim, ist, is)))

            end do

          else
            
            dl_rho(1:gr%mesh%np) = dl_rho(1:gr%mesh%np) + ik_weight*M_TWO* &
                 R_REAL(R_CONJ(psi(1:gr%mesh%np, idim))*lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is))

            do idir = 1, gr%mesh%sb%dim
              gdl_rho(1:gr%mesh%np, idir) = gdl_rho(1:gr%mesh%np, idir) + ik_weight*M_TWO*(&
                R_CONJ(psi(1:gr%mesh%np, idim))*gdl_psi(1:gr%mesh%np, idir) + &
                gpsi(1:gr%mesh%np, idir)*R_CONJ(lr%X(dl_psi)(1:gr%mesh%np, idim, ist, is)))
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

        call states_elec_get_state(st, gr%mesh, ist, is, psi)
        
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/spin_factor

        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, psi(:, idim), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            do ip = 1, gr%mesh%np
              lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) + &
                dl_rho(ip) * ik_weight * sum(abs(gpsi(ip, 1:gr%mesh%sb%dim))**2) + &
                rho(ip) * ik_weight * &
                sum(R_CONJ(gpsi(ip, 1:gr%mesh%sb%dim))*gdl_psi(ip, 1:gr%mesh%sb%dim) + &
                gpsi(ip, 1:gr%mesh%sb%dim) * R_CONJ(gdl_psi_m(ip, 1:gr%mesh%sb%dim)))
            end do

          else 

            do ip = 1, gr%mesh%np
              lr%X(dl_de)(ip, is) = lr%X(dl_de)(ip, is) + &
                dl_rho(ip) * ik_weight * sum(abs(gpsi(ip, 1:gr%mesh%sb%dim))**2) + &
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
          lr%zdl_de(ip, is) = lr%zdl_de(ip, is) + &
            M_TWO * sum(current(ip, 1:gr%mesh%sb%dim, is) * lr%dl_j(ip, 1:gr%mesh%sb%dim, is))
        end if
      end do
    end if

    !now the normalization 
    factor = M_THREE/M_FIVE * (CNST(6.0) * M_PI**2)**M_TWOTHIRD
    do ip = 1, gr%mesh%np

      if(abs(st%rho(ip, is)) >= dmin) then
        d0    = factor * rho(ip)**(CNST(8.0) / M_THREE)
        dl_d0 = CNST(8.0)/M_THREE * factor * dl_rho(ip) * rho(ip)**(M_FIVE/M_THREE)

        lr%X(dl_elf)(ip, is) = &
             M_TWO * d0 * dl_d0/(d0**2 + de(ip, is)**2) * (M_ONE - elf(ip,is)) &
             - M_TWO * de(ip, is) * lr%X(dl_de)(ip, is) / (d0**2 + de(ip, is)**2) * elf(ip, is)
      else
        lr%X(dl_elf)(ip, is) = M_ZERO
      end if

    end do

  end do

  SAFE_DEALLOCATE_A(psi)
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
!> alpha_ij(w) = -e sum(m occ, k) [(<u_mk(0)|-id/dk_i)|u_mkj(1)(w)> + <u_mkj(1)(-w)|(-id/dk_i|u_mk(0)>)]
subroutine X(calc_polarizability_periodic)(sys, em_lr, kdotp_lr, nsigma, zpol, ndir, zpol_k)
  type(electrons_t),   target, intent(inout) :: sys
  type(lr_t),                  intent(inout) :: em_lr(:,:)
  type(lr_t),                  intent(inout) :: kdotp_lr(:)
  integer,                     intent(in)    :: nsigma
  CMPLX,                       intent(out)   :: zpol(:, :) !< (sb%dim, sb%dim)
  integer,           optional, intent(in)    :: ndir
  CMPLX,             optional, intent(out)   :: zpol_k(:, :, :)

  integer :: dir1, dir2, ndir_, ist, ik
  CMPLX :: term, subterm
  type(mesh_t), pointer :: mesh
  logical :: kpt_output
#ifdef HAVE_MPI
  CMPLX :: zpol_temp(1:MAX_DIM, 1:MAX_DIM)
#endif

  PUSH_SUB(X(calc_polarizability_periodic))

  mesh => sys%gr%mesh

  ndir_ = mesh%sb%periodic_dim
  if(present(ndir)) ndir_ = ndir
 
  kpt_output = present(zpol_k)

  if(kpt_output) zpol_k(:, :, :) = M_ZERO

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
          if(kpt_output) zpol_k(dir1, dir2, ik) = zpol_k(dir1, dir2, ik) + term * sys%st%occ(ist, ik)
        end do

      end do
    end do
  end do

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    zpol(1:mesh%sb%periodic_dim, 1:mesh%sb%dim) = zpol_temp(1:mesh%sb%periodic_dim, 1:mesh%sb%dim)
  end if
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    zpol(1:mesh%sb%periodic_dim, 1:mesh%sb%dim) = zpol_temp(1:mesh%sb%periodic_dim, 1:mesh%sb%dim)
  end if
  if(kpt_output) then
    if(sys%st%parallel_in_states) then
      do ik = 1, sys%st%d%nik
        zpol_temp(:, :) = M_ZERO
        call MPI_Allreduce(zpol_k(1:MAX_DIM, 1:MAX_DIM, ik), zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, &
          sys%st%mpi_grp%comm, mpi_err)
        do dir1 = 1, ndir_
          do dir2 = 1, sys%gr%sb%dim
            zpol_k(dir1, dir2, ik) = zpol_temp(dir1, dir2)
          end do
        end do
      end do
    end if
    if(sys%st%d%kpt%parallel) then
      do ik = 1, sys%st%d%nik
        zpol_temp(:, :) = M_ZERO
        call MPI_Allreduce(zpol_k(1:MAX_DIM, 1:MAX_DIM, ik), zpol_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, &
          sys%st%d%kpt%mpi_grp%comm, mpi_err)
        do dir1 = 1, ndir_
          do dir2 = 1, sys%gr%sb%dim
            zpol_k(dir1, dir2, ik) = zpol_temp(dir1, dir2)
          end do
        end do
      end do
    end if
  end if
#endif

  call zsymmetrize_tensor_cart(mesh%sb%symm, zpol)
  
  POP_SUB(X(calc_polarizability_periodic))

end subroutine X(calc_polarizability_periodic)


! ---------------------------------------------------------
!> alpha_ij(w) = - sum(m occ) [<psi_m(0)|r_i|psi_mj(1)(w)> + <psi_mj(1)(-w)|r_i|psi_m(0)>]
!! minus sign is from electronic charge -e
subroutine X(calc_polarizability_finite)(sys, lr, nsigma, perturbation, zpol, doalldirs, ndir)
  type(electrons_t),   target, intent(inout) :: sys
  type(lr_t),                  intent(inout) :: lr(:,:)
  integer,                     intent(in)    :: nsigma
  type(pert_t),                intent(inout) :: perturbation
  CMPLX,                       intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)
  logical,           optional, intent(in)    :: doalldirs
  integer,           optional, intent(in)    :: ndir

  type(states_elec_t), pointer :: st
  integer :: dir1, dir2, ndir_, startdir, np
  R_TYPE, pointer :: psi(:, :, :, :)
  
  PUSH_SUB(X(calc_polarizability_finite))

  st => sys%st
  np = sys%gr%mesh%np
  
  ndir_ = sys%gr%sb%dim
  if(present(ndir)) ndir_ = ndir

  startdir = sys%gr%sb%periodic_dim + 1
  if(present(doalldirs)) then
    if(doalldirs) startdir = 1
  end if

  SAFE_ALLOCATE(psi(1:np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call states_elec_get_state(sys%st, sys%gr%mesh, psi)
  
  do dir1 = startdir, ndir_
    do dir2 = 1, sys%gr%sb%dim
      call pert_setup_dir(perturbation, dir1)
      zpol(dir1, dir2) = -X(pert_expectation_value)(perturbation, sys%namespace, sys%gr, sys%geo, &
        sys%hm, sys%st, psi, lr(dir2, 1)%X(dl_psi))
      if(nsigma == 1) then
        zpol(dir1, dir2) = zpol(dir1, dir2) + R_CONJ(zpol(dir1, dir2))
      else
        zpol(dir1, dir2) = zpol(dir1, dir2) - X(pert_expectation_value)(perturbation, sys%namespace, sys%gr, sys%geo, &
          sys%hm, sys%st, lr(dir2, 2)%X(dl_psi), psi)
      end if
    end do
  end do

  SAFE_DEALLOCATE_P(psi)
  
  POP_SUB(X(calc_polarizability_finite))

end subroutine X(calc_polarizability_finite)


! ---------------------------------------------------------
subroutine X(lr_calc_susceptibility)(sys, lr, nsigma, perturbation, chi_para, chi_dia)
  type(electrons_t),   target, intent(inout) :: sys
  type(lr_t),                  intent(inout) :: lr(:,:)
  integer,                     intent(in)    :: nsigma
  type(pert_t),                intent(inout) :: perturbation
  CMPLX,                       intent(out)   :: chi_para(:,:), chi_dia(:,:)

  type(states_elec_t), pointer :: st
  integer :: dir1, dir2, np
  R_TYPE  :: trace
  R_TYPE, pointer :: psi(:, :, :, :)
  
  PUSH_SUB(X(lr_calc_susceptibility))

  st => sys%st
  np = sys%gr%mesh%np
  
  chi_para = M_ZERO
  chi_dia  = M_ZERO

  SAFE_ALLOCATE(psi(1:np, 1:st%d%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call states_elec_get_state(sys%st, sys%gr%mesh, psi)

  do dir1 = 1, sys%gr%sb%dim
    do dir2 = 1, sys%gr%sb%dim

      call pert_setup_dir(perturbation, dir1, dir2)

      trace = X(pert_expectation_value)(perturbation, sys%namespace, sys%gr, sys%geo, sys%hm, sys%st, psi, &
        lr(dir2, 1)%X(dl_psi))
      
      if (nsigma == 1) then 
        trace = trace + R_CONJ(trace)
      else
        trace = trace + X(pert_expectation_value)(perturbation, sys%namespace, sys%gr, sys%geo, sys%hm, sys%st, &
          lr(dir2, 2)%X(dl_psi), psi)
      end if
     
      ! first the paramagnetic term 
      chi_para(dir1, dir2) = chi_para(dir1, dir2) + trace

      chi_dia(dir1, dir2) = chi_dia(dir1, dir2) + X(pert_expectation_value)(perturbation, sys%namespace, sys%gr, sys%geo, &
        sys%hm, sys%st, psi, psi, pert_order=2)

    end do
  end do
  
  SAFE_DEALLOCATE_P(psi)
  
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
!> See (16) in X Andrade et al., J. Chem. Phys. 126, 184106 (2006) for finite systems
!! and (10) in A Dal Corso et al., Phys. Rev. B 15, 15638 (1996) for periodic systems
!! Supply only em_lr for finite systems, and both kdotp_lr and kdotp_em_lr for periodic
!! em_lr(dir, sigma, omega) = electric perturbation of ground-state wavefunctions
!! kdotp_lr(dir) = kdotp perturbation of ground-state wavefunctions
!! kdotp_em_lr(dir1, dir2, sigma, omega) = kdotp perturbation of electric-perturbed wfns
subroutine X(lr_calc_beta) (sh, sys, em_lr, dipole, beta, kdotp_lr, kdotp_em_lr, occ_response, dl_eig)
  type(sternheimer_t),         intent(inout) :: sh
  type(electrons_t),   target, intent(inout) :: sys
  type(lr_t),                  intent(inout) :: em_lr(:,:,:)
  type(pert_t),                intent(inout) :: dipole
  CMPLX,                       intent(out)   :: beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM)
  type(lr_t),        optional, intent(in)    :: kdotp_lr(:)
  type(lr_t),        optional, intent(in)    :: kdotp_em_lr(:,:,:,:) !< kdotp dir, em dir, sigma, factor
  logical,           optional, intent(in)    :: occ_response !< do the wfns include the occ subspace?
  !< occ_response = yes is based on Baroni et al., RMP 73, 515 (2001), eqn 122
  !! occ_response = no  is based on Baroni et al., RMP 73, 515 (2001), eqn 123
  !!   The occ_response = no version can be used even if the wfns do include the
  !!   occupied subspace, it is just more efficient to use the other formula.
  FLOAT,         optional, intent(in)    :: dl_eig(:,:,:) !< state, kpt, dir

  type(states_elec_t), pointer :: st
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

  call profiling_in(beta_prof, "X(CALC_BETA)")

  st   => sys%st
  mesh => sys%gr%mesh
  np   =  mesh%np
  ndir =  sys%gr%sb%dim

  occ_response_ = .false.
  if(present(occ_response)) occ_response_ = occ_response

  ASSERT(present(kdotp_lr) .eqv. present(kdotp_em_lr))
  ASSERT(present(dl_eig) .eqv. present(kdotp_em_lr))
  ! either both are absent for finite, or both present for periodic

  if(sternheimer_add_fxc(sh)) then
    !calculate kxc, the derivative of fxc
    SAFE_ALLOCATE(kxc(1:np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
    kxc = M_ZERO
    
    SAFE_ALLOCATE(rho(1:np, 1:st%d%nspin))
    call states_elec_total_density(st, mesh, rho)
    
    call xc_get_kxc(sys%ks%xc, mesh, sys%namespace, rho, st%d%ispin, kxc)
    SAFE_DEALLOCATE_A(rho)
    SAFE_ALLOCATE(hpol_density(1:np))
  end if

  SAFE_ALLOCATE(tmp(1:np, 1:st%d%dim))
  SAFE_ALLOCATE(ppsi(1:np, 1:st%d%dim))
  SAFE_ALLOCATE(hvar(1:np, 1:st%d%nspin, 1:2, 1:st%d%dim, 1:ndir, 1:3))
  SAFE_ALLOCATE(me010(1:st%nst, 1:st%nst, 1:mesh%sb%dim, 1:3, st%d%kpt%start:st%d%kpt%end))
  SAFE_ALLOCATE(me11(1:mesh%sb%dim, 1:mesh%sb%dim, 1:3, 1:3, 1:2, st%d%kpt%start:st%d%kpt%end))

  do ifreq = 1, 3
    do idir = 1, ndir
      do idim = 1, st%d%dim
        call X(sternheimer_calc_hvar)(sh, sys, em_lr(idir, :, ifreq), 2, hvar(:, :, :, idim, idir, ifreq))
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

                  ispin = states_elec_dim_get_spin_index(sys%st%d, ik)

                  if (present(kdotp_em_lr) .and. u(2) <= sys%gr%sb%periodic_dim) then
                    tmp(1:np, idim) = kdotp_em_lr(u(2), u(3), isigma, w(3))%X(dl_psi)(1:np, idim, ist, ik)
                  else
                    call pert_setup_dir(dipole, u(2))
                    call X(pert_apply)(dipole, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
                      em_lr(u(3), isigma, w(3))%X(dl_psi)(:, :, ist, ik), tmp)
                  end if

                  do ip = 1, np
                    tmp(ip, idim) = tmp(ip, idim) + hvar(ip, ispin, isigma, idim, u(2), w(2)) &
                      * em_lr(u(3), isigma, w(3))%X(dl_psi)(ip, idim, ist, ik)
                  end do

                  beta(ii, jj, kk) = beta(ii, jj, kk) &
                    - M_HALF * st%d%kweights(ik) * st%smear%el_per_state &
                    * X(mf_dotp)(mesh, em_lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, idim, ist, ik), tmp(1:np, idim))
                  
                  do ist2 = 1, st%nst
                    if(occ_response_ .and. ist2 /= ist) cycle
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
                    end do
                  end do
                end do
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
  end if
  SAFE_DEALLOCATE_A(tmp)
  SAFE_DEALLOCATE_A(ppsi)
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

    R_TYPE, allocatable :: psi1(:, :), psi2(:, :)

    PUSH_SUB(X(lr_calc_beta).get_matrix_elements)
    
    SAFE_ALLOCATE(psi1(1:sys%gr%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(psi2(1:sys%gr%mesh%np_part, 1:st%d%dim))
    
    do ik = st%d%kpt%start, st%d%kpt%end
      ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
      do ist = 1, st%nst

        call states_elec_get_state(st, sys%gr%mesh, ist, ik, psi1)
        
        do ist2 = 1, st%nst

          if(occ_response_ .and. ist2 /= ist) cycle

          call states_elec_get_state(st, sys%gr%mesh, ist2, ik, psi2)

          do ii = 1, ndir
            call pert_setup_dir(dipole, ii)

            do ifreq = 1, 3

              ! ist = ist2 term cannot be captured by k.p perturbation
              if (present(kdotp_lr) .and. ii <= sys%gr%sb%periodic_dim) then
                if(ist == ist2) then
                  ppsi = M_ZERO ! will put in eigenvalue directly later
                else
                  do idim = 1, st%d%dim
                    do ip = 1, np
                      ppsi(ip, idim) = kdotp_lr(ii)%X(dl_psi)(ip, idim, ist, ik)
                    end do
                  end do
                end if
              else
                call X(pert_apply)(dipole, sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi1, ppsi)
              end if

              isigma = 1
              do idim = 1, st%d%dim
                do ip = 1, np
                  ppsi(ip, idim) = ppsi(ip, idim) + hvar(ip, ispin, isigma, idim, ii, ifreq)*psi1(ip, idim)
                end do
              end do

              me010(ist2, ist, ii, ifreq, ik) = X(mf_dotp)(mesh, st%d%dim, psi2, ppsi)

              if (present(kdotp_lr) .and. ii <= sys%gr%sb%periodic_dim .and. ist == ist2) then
                me010(ist, ist, ii, ifreq, ik) = me010(ist, ist, ii, ifreq, ik) + dl_eig(ist, ik, ii)
              end if

            end do
!            if(mpi_grp_is_root(mpi_world)) write(10,'(4i3,f10.6)') ist2, ist, ii, ik, me010(ist2, ist, ii, 1, ik)
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi1)
    SAFE_DEALLOCATE_A(psi2)    

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
                  end do
                else
                  call states_elec_blockt_mul(mesh, st, st%st_start, st%st_start, &
                    em_lr(ii, op_sigma, ifreq)%X(dl_psi)(:, :, :, ik), &
                    em_lr(jj, isigma, jfreq)%X(dl_psi)(:, :, :, ik), &
                    me11(ii, jj, ifreq, jfreq, isigma, ik)%X(matrix))
                end if

!                if(ifreq == 1 .and. jfreq == 1 .and. mpi_grp_is_root(mpi_world)) then
!                  do ist = 1, st%nst
!                    do ist2 = 1, st%nst
!                      write(11,'(6i3,f10.6)') ii, jj, isigma, ik, ist, ist2, me11(ii, jj, 1, 1, isigma, ik)%X(matrix)(ist, ist2)
!                    end do
!                  end do
!                end if
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
subroutine X(post_orthogonalize)(sys, nfactor, nsigma, freq_factor, omega, eta, em_lr, kdotp_em_lr)
  type(electrons_t),   intent(in)    :: sys
  integer,             intent(in)    :: nfactor
  integer,             intent(in)    :: nsigma
  FLOAT,               intent(in)    :: freq_factor(:)
  FLOAT,               intent(in)    :: omega
  FLOAT,               intent(in)    :: eta                  !< should be zero when wfns are real
  type(lr_t),          intent(inout) :: em_lr(:,:,:)         !< em dir, sigma, factor
  type(lr_t),          intent(inout) :: kdotp_em_lr(:,:,:,:) !< kdotp dir, em dir, sigma, factor

  integer :: kdotp_dir, em_dir, isigma, ifactor
  R_TYPE :: frequency

  PUSH_SUB(X(post_orthogonalize))

#ifdef R_TREAL
  if(abs(eta) > M_EPSILON) then
    message(1) = "Internal error: dpost_orthogonalize cannot be called with argument eta != 0"
    call messages_fatal(1)
  end if
#endif

  do ifactor = 1, nfactor
    do isigma = 1, nsigma
      frequency = R_TOPREC(omega * freq_factor(ifactor) + M_zI * eta)
      if(isigma == 2) frequency = -R_CONJ(frequency)
      
      do em_dir = 1, sys%gr%sb%dim
        call X(lr_orth_response)(sys%gr%mesh, sys%st, em_lr(em_dir, isigma, ifactor), frequency)
        
        do kdotp_dir = 1, sys%gr%sb%periodic_dim
          call X(lr_orth_response)(sys%gr%mesh, sys%st, kdotp_em_lr(kdotp_dir, em_dir, isigma, ifactor), frequency)
        end do
      end do
    end do
  end do

  POP_SUB(X(post_orthogonalize))
end subroutine X(post_orthogonalize)


! --------------------------------------------------------- 
! <\psi|-i d/dk|\psi> cannot be calculated from kdotp perturbation which gives only diagonal matrix elements
! but it can be approximated for large cells along the lines of the "single-point Berry phase":
! <\psi|-i d/dk|\psi> =~ (L/2 \pi) Im <\psi|exp(2 \pi i x / L)|\psi> ( =~ <\psi|x|\psi> )
subroutine X(em_resp_calc_eigenvalues)(sys, dl_eig)
  type(electrons_t),   intent(in)  :: sys
  FLOAT,               intent(out) :: dl_eig(:,:,:) !< (ist, ik, idir)

  integer :: ik, ist, ip, idim, idir
  R_TYPE, allocatable :: psi(:,:)
  CMPLX, allocatable :: integrand(:)
#ifdef HAVE_MPI
  FLOAT, allocatable :: dl_eig_temp(:,:,:)
#endif

  PUSH_SUB(X(em_resp_calc_eigenvalues))

  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  SAFE_ALLOCATE(integrand(1:sys%gr%mesh%np))

  dl_eig(:, :, :) = M_ZERO

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    do ist = sys%st%st_start, sys%st%st_end

      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)

      do idir = 1, sys%gr%sb%periodic_dim
        do idim = 1, sys%st%d%dim
          do ip = 1, sys%gr%mesh%np
            integrand(ip) = exp(M_zI*(M_PI/sys%gr%mesh%sb%lsize(idir))*sys%gr%mesh%x(ip, idir)) * abs(psi(ip, idim))**2
          end do
          dl_eig(ist, ik, idir) = dl_eig(ist, ik, idir) + &
            (sys%gr%mesh%sb%lsize(idir)/M_PI) * aimag(zmf_integrate(sys%gr%mesh, integrand))
        end do
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(integrand)

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states .or. sys%st%d%kpt%parallel) then
    SAFE_ALLOCATE(dl_eig_temp(1:sys%st%nst, 1:sys%st%d%nik, 1:sys%gr%sb%periodic_dim))

    call MPI_Allreduce(dl_eig, dl_eig_temp, sys%st%nst * sys%st%d%nik * sys%gr%sb%periodic_dim, &
      MPI_FLOAT, MPI_SUM, sys%st%st_kpt_mpi_grp%comm, mpi_err)

    dl_eig(:,:,:) = dl_eig_temp(:,:,:)
    SAFE_DEALLOCATE_A(dl_eig_temp)
  end if
#endif

  POP_SUB(X(em_resp_calc_eigenvalues))
end subroutine X(em_resp_calc_eigenvalues)

! ---------------------------------------------------------
subroutine X(lr_calc_magneto_optics_finite)(sh, sh_mo, sys, nsigma, nfactor, lr_e, &
  lr_b, chi)
  type(sternheimer_t),    intent(inout) :: sh, sh_mo
  type(electrons_t),      intent(inout) :: sys
  integer,                intent(in)    :: nsigma
  integer,                intent(in)    :: nfactor
  type(lr_t),             intent(inout) :: lr_e(:,:,:)
  type(lr_t),             intent(inout) :: lr_b(:,:)
  CMPLX,                  intent(out)   :: chi(:,:,:)

  integer :: dir1, dir2, dir3, ist, ist_occ, dir(nfactor)
  integer :: ik, idim, ip, isigma, ispin, ii
  type(pert_t) :: pert_m, pert_e(nfactor)
  R_TYPE, allocatable :: pertpsi_e(:,:,:), pertpsi_b(:,:)
  FLOAT :: weight
  R_TYPE, allocatable :: hpol_density(:), hvar_e(:,:,:,:,:), hvar_b(:,:,:,:)
  type(matrix_t) :: prod_eb(nsigma), prod_ee 
  R_TYPE, allocatable :: psi(:,:), psi1(:,:)
  CMPLX :: factor
  FLOAT :: factor1
  logical :: calc_var, calc_var_mo
  R_TYPE :: term(nsigma)
    
  PUSH_SUB(X(lr_calc_magneto_optics_finite))

  SAFE_ALLOCATE(pertpsi_e(1:sys%gr%mesh%np, 1:sys%hm%d%dim, 1:nsigma))
  SAFE_ALLOCATE(pertpsi_b(1:sys%gr%mesh%np, 1:sys%hm%d%dim))
  
  calc_var = sternheimer_add_fxc(sh) .or. sternheimer_add_hartree(sh)
  calc_var_mo = sternheimer_add_fxc(sh_mo) .or. sternheimer_add_hartree(sh_mo)
  if(sternheimer_add_fxc(sh)) then 
    SAFE_ALLOCATE(hpol_density(1:sys%gr%mesh%np))
  end if
  
  do isigma = 1, nsigma
    SAFE_ALLOCATE(prod_eb(isigma)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  end do
  SAFE_ALLOCATE(prod_ee%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  if(calc_var) then
    SAFE_ALLOCATE(hvar_e(1:sys%gr%sb%dim, 1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:nsigma, 1:nfactor))
  end if
  if(calc_var_mo) then
    SAFE_ALLOCATE(hvar_b(1:sys%gr%sb%dim, 1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))
  end if
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  SAFE_ALLOCATE(psi1(1:sys%gr%mesh%np, 1:sys%st%d%dim))

#if defined(R_TCOMPLEX)
  factor = M_ONE        
  factor1 = M_ONE
#else
  factor = -M_zI
  factor1 = -M_ONE
#endif

  chi = M_ZERO

  pertpsi_e(:,:,:) = M_ZERO
  pertpsi_b(:,:) = M_ZERO
  
  call pert_init(pert_m, sys%namespace, PERTURBATION_MAGNETIC, sys%gr, sys%geo)
  do ii = 1, nfactor
    call pert_init(pert_e(ii), sys%namespace, PERTURBATION_ELECTRIC, sys%gr, sys%geo)
  end do
  
  psi(:,:) = M_ZERO
  psi1(:,:) = M_ZERO
  
  if(calc_var) then
    hvar_e(:,:,:,:,:) = M_ZERO
    do dir1 = 1, sys%gr%sb%dim
      do ii = 1, nfactor
        call X(sternheimer_calc_hvar)(sh, sys, lr_e(dir1, :, ii), nsigma, hvar_e(dir1, :, :, :, ii))
      end do
    end do
  end if
  if(calc_var_mo) then
    hvar_b(:,:,:,:) = M_ZERO
    do dir1 = 1, sys%gr%sb%dim
      call X(sternheimer_calc_hvar)(sh_mo, sys, lr_b(dir1, :), 1, hvar_b(dir1, :, :, :))
    end do
  end if
  
  do dir1 = 1, sys%gr%sb%dim
    do dir2 = 1, sys%gr%sb%dim
      do dir3 = 1, sys%gr%sb%dim
        dir(1) = dir2
        dir(nfactor) = dir1
        call pert_setup_dir(pert_m, dir3)
        do ii = 1, nfactor
          call pert_setup_dir(pert_e(ii), dir(ii))
        end do
        do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
          ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
          weight = sys%st%d%kweights(ik) * sys%st%smear%el_per_state
          do ist = 1, sys%st%nst
            if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
              do ii = 1, nfactor
                do isigma = 1, nsigma
                  call X(pert_apply)(pert_e(swap_sigma(ii)), sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
                    lr_e(dir(ii), isigma, ii)%X(dl_psi)(:, :, ist, ik), pertpsi_e(:, :, isigma))
                  if(calc_var) then
                    do idim = 1, sys%hm%d%dim
                      do ip = 1, sys%gr%mesh%np
                        pertpsi_e(ip, idim, isigma) = pertpsi_e(ip, idim, isigma) + &
                          hvar_e(dir(swap_sigma(ii)), ip, ispin, isigma, swap_sigma(ii)) * &
                          lr_e(dir(ii), isigma, ii)%X(dl_psi)(ip, idim, ist, ik)
                      end do
                    end do
                  end if
                  term(isigma) = weight * factor1 * X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, &
                    lr_b(dir3, 1)%X(dl_psi)(:, :, ist, ik), pertpsi_e(:, :, isigma))
                end do
                chi(dir(nfactor), dir(1), dir3) = chi(dir(nfactor), dir(1), dir3) + term(1) + R_CONJ(term(nsigma))

                call X(pert_apply)(pert_m, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
                  lr_e(dir(ii), 1, ii)%X(dl_psi)(:, :, ist, ik), pertpsi_b(:, :))
                if(calc_var_mo) then
                  do idim = 1, sys%hm%d%dim
                    do ip = 1, sys%gr%mesh%np
                      pertpsi_b(ip, idim) = pertpsi_b(ip, idim) + hvar_b(dir3, ip, ispin, 1) * &
                        lr_e(dir(ii), 1, ii)%X(dl_psi)(ip, idim, ist, ik)
                    end do
                  end do
                end if
                term(1) = weight * factor1 * X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, &
                  lr_e(dir(swap_sigma(ii)), nsigma, swap_sigma(ii))%X(dl_psi)(:, :, ist, ik), &
                  pertpsi_b(:, :))
                chi(dir(nfactor), dir(1), dir3) = chi(dir(nfactor), dir(1), dir3) + term(1) 
              end do
            end if
          end do
          do ii = 1, nfactor 
            do isigma = 1, nsigma
              call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
                lr_e(dir(swap_sigma(ii)), swap_sigma(isigma), swap_sigma(ii))%X(dl_psi)(:, :, :, ik), &
                lr_b(dir3, 1)%X(dl_psi)(:, :, :, ik), prod_eb(isigma)%X(matrix))
            end do
            call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
              lr_e(dir(ii), nsigma, ii)%X(dl_psi)(:, :, :, ik), lr_e(dir(swap_sigma(ii)), 1, &
              swap_sigma(ii))%X(dl_psi)(:, :, :, ik), prod_ee%X(matrix))

            do ist = 1, sys%st%nst
              if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
                call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
                call X(pert_apply)(pert_e(ii), sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi, pertpsi_e(:, :, 1))
                if(nfactor > 1) pertpsi_e(:, :, nfactor) = pertpsi_e(:, :, 1)
                if(calc_var) then
                  do idim = 1, sys%hm%d%dim
                    do ip = 1, sys%gr%mesh%np
                      do isigma = 1, nsigma
                        pertpsi_e(ip, idim, isigma) = pertpsi_e(ip, idim, isigma) + &
                          hvar_e(dir(ii), ip, ispin, isigma, ii) * psi(ip, idim)
                      end do
                    end do
                  end do
                end if
                call X(pert_apply)(pert_m, sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi, pertpsi_b(:, :))
                if(calc_var_mo) then
                  do idim = 1, sys%hm%d%dim
                    do ip = 1, sys%gr%mesh%np
                      pertpsi_b(ip, idim) = pertpsi_b(ip, idim) + hvar_b(dir3, ip, ispin, 1) * psi(ip, idim)
                    end do
                  end do
                end if
                do ist_occ = 1, sys%st%nst
                  if(abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then
                    call states_elec_get_state(sys%st, sys%gr%mesh, ist_occ, ik, psi1)
                    do isigma = 1, nsigma
                      term(isigma) = - weight * X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi1, pertpsi_e(:, :, isigma)) &
                        * prod_eb(isigma)%X(matrix)(ist, ist_occ) 
                    end do
                    chi(dir(nfactor), dir(1), dir3) = chi(dir(nfactor), dir(1), dir3) + term(1) + R_CONJ(term(nsigma))

                    term(1) = - weight * X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi1, pertpsi_b(:, :))* &
                      prod_ee%X(matrix)(ist, ist_occ) 
                    chi(dir(nfactor), dir(1), dir3) = chi(dir(nfactor), dir(1), dir3) + term(1) 
                  end if
                end do
              end if
            end do
          end do
        end do
      end do
    end do
  end do
  
  do dir1 = 1, sys%gr%sb%dim
    do dir2 = 1, sys%gr%sb%dim
      dir(1) = dir2
      dir(nfactor) = dir1
      do dir3 = 1, sys%gr%sb%dim
        if(sternheimer_add_fxc(sh_mo) .and. sternheimer_add_fxc(sh)) then 
          hpol_density(:) = M_ZERO
          do ii = 1, nfactor
            call X(calc_kvar_energy)(sh_mo, sys, lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), &
              lr_e(dir(ii), 1, ii), lr_b(dir3, 1), hpol_density)
            call X(calc_kvar_energy)(sh_mo, sys, lr_e(dir(ii), 1, ii), lr_b(dir3, 1), &
              lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), hpol_density)
            call X(calc_kvar_energy)(sh_mo, sys, lr_b(dir3, 1), &
              lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), lr_e(dir(ii), 1, ii), hpol_density)
          end do
          chi(dir(nfactor), dir(1), dir3) = chi(dir(nfactor), dir(1), dir3) + &
            M_ONE/CNST(6.0) * X(mf_integrate)(sys%gr%mesh, hpol_density) ! the coefficient should be checked
        end if
      end do
    end do
  end do
  
  chi(:,:,:) = chi(:,:,:) / P_C * factor
  
  do isigma = 1, nsigma
    SAFE_DEALLOCATE_P(prod_eb(isigma)%X(matrix))
  end do
  SAFE_DEALLOCATE_P(prod_ee%X(matrix))
  
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi1)
  SAFE_DEALLOCATE_A(pertpsi_e)
  SAFE_DEALLOCATE_A(pertpsi_b)
  if(calc_var) then
    SAFE_DEALLOCATE_A(hvar_e)
  end if
  if(calc_var_mo) then
    SAFE_DEALLOCATE_A(hvar_b)
  end if
  
  if(sternheimer_add_fxc(sh)) then 
    SAFE_DEALLOCATE_A(hpol_density)
  end if
  POP_SUB(X(lr_calc_magneto_optics_finite))
end subroutine X(lr_calc_magneto_optics_finite)

! ---------------------------------------------------------
subroutine X(calc_kvar_energy)(sh_mo, sys, lr1, lr2, lr3, hpol_density)
  type(sternheimer_t),    intent(inout) :: sh_mo
  type(electrons_t),      intent(inout) :: sys
  type(lr_t),             intent(inout) :: lr1, lr2, lr3
  R_TYPE,                 intent(inout) :: hpol_density(:)
    
  R_TYPE, allocatable :: kvar(:,:,:) 
  integer :: is, ip
    
  PUSH_SUB(X(calc_kvar_energy))
  
  SAFE_ALLOCATE(kvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))
  
  call X(calc_kvar)(sh_mo, sys, lr1%X(dl_rho), lr2%X(dl_rho), 1, kvar)
  do ip = 1, sys%gr%mesh%np
    do is = 1, sys%st%d%nspin
      hpol_density(ip) = hpol_density(ip) + kvar(ip,is,1) * lr3%X(dl_rho)(ip,is)
    end do
  end do
  
  SAFE_DEALLOCATE_A(kvar)
  
  POP_SUB(X(calc_kvar_energy))
end subroutine X(calc_kvar_energy)

! ---------------------------------------------------------
subroutine X(lr_calc_magneto_optics_periodic)(sh, sh2, sys, nsigma, nfactor, &
  nfactor_ke, freq_factor, lr_e, lr_b, lr_k, lr_k2, lr_ke, lr_kb, frequency, zpol,&
  zpol_kout)
  type(sternheimer_t),  intent(inout) :: sh, sh2
  type(electrons_t),    intent(inout) :: sys 
  integer,              intent(in)    :: nsigma
  integer,              intent(in)    :: nfactor
  integer,              intent(in)    :: nfactor_ke
  FLOAT,                intent(in)    :: freq_factor(:)
  type(lr_t),           intent(inout) :: lr_e(:,:,:)
  type(lr_t),           intent(inout) :: lr_b(:,:)
  type(lr_t),           intent(inout) :: lr_k(:,:)
  type(lr_t),           intent(inout) :: lr_k2(:,:,:)
  type(lr_t),           intent(inout) :: lr_ke(:,:,:,:)
  type(lr_t),           intent(inout) :: lr_kb(:,:,:)
  CMPLX,                intent(in)    :: frequency
  CMPLX,                intent(inout) :: zpol(:,:,:)
  CMPLX, optional,      intent(inout) :: zpol_kout(:,:,:,:)
  
  integer :: idir1, idir2, idir3, ist, ii, dir(nfactor), &
    ispin, idim, ndim, np, ik, ndir, ist_occ, isigma, isigma_alt, ip
  FLOAT :: weight
  R_TYPE, allocatable :: gpsi(:,:,:,:), gdl_e(:,:,:), gdl_k(:,:,:,:), gdl_b(:,:,:), &
    gdl_ke(:,:,:,:), hvar(:,:,:,:,:), hvar2(:,:,:,:)
  R_TYPE:: factor, factor0, factor1, factor_e, factor_rho, factor_sum
  type(matrix_t) :: mat_be, mat_eb
  type(matrix_t), allocatable:: mat_kek(:,:), mat_kke(:,:), mat_g(:,:)
  R_TYPE, allocatable :: psi_be(:,:,:,:), density(:)
  logical :: add_hartree1, add_fxc1, add_hartree2, add_fxc2, kpt_output
  type(lr_t) :: lr0(1)
  CMPLX :: zpol0(1:MAX_DIM,1:MAX_DIM,1:MAX_DIM), zpol_k(1:MAX_DIM,1:MAX_DIM,1:MAX_DIM)  

#ifdef HAVE_MPI
  CMPLX :: zpol_temp(1:MAX_DIM,1:MAX_DIM,1:MAX_DIM)
#endif

#if defined(R_TCOMPLEX)
  factor = -M_zI
  factor0 = -M_zI
  factor1 = M_ONE
  factor_rho = M_zI * M_HALF
#else
  factor = -M_ONE
  factor0 = M_ONE
  factor1 = -M_ONE
  factor_rho = M_HALF
#endif
  factor_e = -M_ONE 
  factor_sum = -M_ONE 

  PUSH_SUB(X(lr_calc_magneto_optics_periodic))
  
  ASSERT(abs(frequency) .gt. M_EPSILON)
  ASSERT(nfactor == 2)
  
  np = sys%gr%mesh%np
  ndir = sys%gr%mesh%sb%dim
  ndim = sys%st%d%dim

  SAFE_ALLOCATE(gpsi(1:np, 1:ndim, 1:sys%st%nst, 1:nsigma))
  SAFE_ALLOCATE(gdl_e(1:np, 1:ndim, 1:nsigma))
  SAFE_ALLOCATE(gdl_k(1:ndir, 1:np, 1:ndim, 1:nsigma))
  SAFE_ALLOCATE(gdl_b(1:np, 1:ndim, 1:nsigma))
  SAFE_ALLOCATE(gdl_ke(1:ndir, 1:np, 1:ndim, 1:nsigma))
  SAFE_ALLOCATE(mat_kke(1:ndir, 1:ndir))
  SAFE_ALLOCATE(mat_kek(1:ndir, 1:ndir))
  SAFE_ALLOCATE(mat_g(1:ndir, 1:nsigma))
  SAFE_ALLOCATE(mat_eb%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(mat_be%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi_be(1:np, 1:ndim, 1:sys%st%nst, 1:nsigma))
  SAFE_ALLOCATE(hvar(1:np, 1:sys%st%d%nspin, 1:nsigma, 1:ndir, 1:nfactor))
  SAFE_ALLOCATE(hvar2(1:np, 1:sys%st%d%nspin, 1:1, 1:ndir))

  kpt_output = present(zpol_kout)
  do idir1 = 1, ndir
    do isigma = 1, nsigma
      SAFE_ALLOCATE(mat_g(idir1, isigma)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
    end do
    do idir2 = 1, ndir
      SAFE_ALLOCATE(mat_kek(idir1,idir2)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
      SAFE_ALLOCATE(mat_kke(idir1,idir2)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
    end do
  end do

  hvar(:,:,:,:,:) = M_ZERO
  hvar2(:,:,:,:) = M_ZERO
  zpol(:,:,:) = M_ZERO 
  zpol0(:,:,:) = M_ZERO 
  if(kpt_output) zpol_kout(:,:,:,:) = M_ZERO  

  add_hartree1 = sternheimer_add_hartree(sh) 
  add_fxc1 = sternheimer_add_fxc(sh)
  add_hartree2 = sternheimer_add_hartree(sh2) 
  add_fxc2 = sternheimer_add_fxc(sh2)
  
  call lr_init(lr0(1))
  if(add_hartree2 .or. add_fxc2) then
    call lr_allocate(lr0(1), sys%st, sys%gr%mesh, allocate_rho = .true.)
    if(add_fxc2) then
      SAFE_ALLOCATE(density(1:np))
    end if
    do idir1 = 1, ndir
      lr0(1)%X(dl_rho) = M_ZERO 
      call X(calc_rho)(sys, factor_rho, factor_sum, factor, &
        factor, lr_k(magn_dir(idir1, 1), 1), lr_k(magn_dir(idir1, 2), 1), lr0(1))
      call X(calc_rho)(sys, factor_sum * factor_rho, factor_sum, &
        factor, factor, lr_k(magn_dir(idir1, 2), 1), lr_k(magn_dir(idir1, 1), 1), lr0(1))
      
      if(sys%st%parallel_in_states .or. sys%st%d%kpt%parallel) &
        call comm_allreduce(sys%st%st_kpt_mpi_grp%comm, lr0(1)%X(dl_rho))

      do ip = 1, sys%gr%mesh%np 
        do ispin = 1, sys%st%d%nspin 
          lr0(1)%X(dl_rho)(ip, ispin) = lr0(1)%X(dl_rho)(ip, ispin) + lr_b(idir1, 1)%X(dl_rho)(ip, ispin)
        end do  
      end do
      call X(sternheimer_calc_hvar)(sh2, sys, lr0, 1, hvar2(:,:,:, idir1))
      if(add_fxc2 .and. add_fxc1) then  
        do idir2 = 1, ndir
          do idir3 = 1, ndir 
            dir(1) = idir2
            dir(nfactor) = idir3   
            density(:) = M_ZERO
            do ii = 1, nfactor
              call X(calc_kvar_energy)(sh2, sys, lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), &
                lr_e(dir(ii), 1, ii), lr0(1), density)
              call X(calc_kvar_energy)(sh2, sys, lr_e(dir(ii), 1, ii), lr0(1), &
                lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), density)
              call X(calc_kvar_energy)(sh2, sys, lr0(1), lr_e(dir(swap_sigma(ii)), 1, swap_sigma(ii)), &
                lr_e(dir(ii), 1, ii), density)
            end do
            if(nfactor_ke > 1) then  
              do ii = 1, nfactor
                call X(calc_kvar_energy)(sh2, sys, lr_e(dir(ii), 1, swap_sigma(ii)), &
                  lr_e(dir(swap_sigma(ii)), 1, ii), lr0(1), density)
                call X(calc_kvar_energy)(sh2, sys, lr_e(dir(swap_sigma(ii)), 1, ii), lr0(1), &
                  lr_e(ii, 1, swap_sigma(ii)), density)
                call X(calc_kvar_energy)(sh2, sys, lr0(1), lr_e(dir(ii), 1, swap_sigma(ii)), &
                  lr_e(dir(swap_sigma(ii)), 1, ii), density)
              end do
            end if
            zpol0(idir3, idir2, idir1) = zpol0(idir3, idir2, idir1) - & 
              frequency/P_C * factor0/CNST(6.0) * X(mf_integrate)(sys%gr%mesh, density) 
          end do 
        end do
      end if
    end do
    if(add_fxc2) then 
      SAFE_DEALLOCATE_A(density) 
    end if
  else
    call lr_allocate(lr0(1), sys%st, sys%gr%mesh, allocate_rho = .false.)
  end if   
  
  if(add_hartree1 .or. add_fxc1) then
    do idir1 = 1, ndir
      do ii = 1, nfactor
        call X(sternheimer_calc_hvar)(sh, sys, lr_e(idir1, :, ii), nsigma, hvar(:,:,:, idir1, ii))
      end do
    end do
  end if
  
  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
    weight = sys%st%d%kweights(ik) * sys%st%smear%el_per_state
    zpol_k(:,:,:) = M_ZERO
    do ii = 1, nfactor_ke
      do idir1 = 1, ndir
        gpsi(:,:,:,:) = M_ZERO 
        do idir2 = 1, ndir
          psi_be(:,:,:,:) = M_ZERO
          do isigma = 1, nsigma
            if(nsigma == 2) isigma_alt = swap_sigma(isigma)
            call X(inhomog_EB)(sys, ik, add_hartree1, add_fxc1, hvar(:, ispin, isigma, idir2, ii), &
              lr_b(idir1, 1)%X(dl_psi)(:,:,:, ik), lr_kb(idir2, idir1, 1)%X(dl_psi)(:,:,:, ik), &
              factor1, psi_be(:,:,:, isigma), lr_k(magn_dir(idir1, 1), 1)%X(dl_psi)(:,:,:, ik), &
              lr_k(magn_dir(idir1, 2), 1)%X(dl_psi)(:,:,:, ik))

            call X(inhomog_BE)(sys, magn_dir(idir1, 1), magn_dir(idir1, 2), ik, & 
              add_hartree2, add_fxc2, hvar2(:, ispin, 1, idir1),  & 
              lr_e(idir2, isigma, ii)%X(dl_psi)(:,:,:, ik), lr_e(idir2, isigma_alt, ii)%X(dl_psi)(:,:,:, ik), &
              lr_ke(magn_dir(idir1, 1), idir2, isigma, ii)%X(dl_psi)(:,:,:, ik), &
              lr_ke(magn_dir(idir1, 2), idir2, isigma, ii)%X(dl_psi)(:,:,:, ik), &
              lr_k(magn_dir(idir1, 1), 1)%X(dl_psi)(:,:,:, ik), &
              lr_k(magn_dir(idir1, 2), 1)%X(dl_psi)(:,:,:, ik), &
              factor_e, factor_e, psi_be(:,:,:, isigma))
          end do
          do ist = 1, sys%st%nst
            if (abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
              do idir3 = 1, ndir 
                do idim = 1, ndim    
                  zpol_k(idir3, idir2, idir1) = zpol_k(idir3, idir2, idir1) - freq_factor(ii) * & 
                    frequency / P_C * factor_e * (X(mf_dotp)(sys%gr%mesh, psi_be(:, idim, ist, nsigma), &
                    factor0 * lr_e(idir3, 1, swap_sigma(ii))%X(dl_psi)(:, idim, ist, ik))&
                    - X(mf_dotp)(sys%gr%mesh, factor0 * lr_e(idir3, nsigma, swap_sigma(ii))%X(dl_psi)(:, idim, ist, ik), &
                    psi_be(:, idim, ist, 1)))

                  zpol_k(idir3, idir2, idir1) = zpol_k(idir3, idir2, idir1) - M_ONE / P_C *&
                    (X(mf_dotp)(sys%gr%mesh, psi_be(:, idim, ist, nsigma), &
                    factor0 * lr_k(idir3, 1)%X(dl_psi)(:, idim, ist, ik))&
                    + X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(idir3, 1)%X(dl_psi)(:, idim, ist, ik), &
                    psi_be(:, idim, ist, 1)))
                end do
              end do
            end if
          end do
        end do
        do ist = 1, sys%st%nst
          if (abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
            call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, lr0(1)%X(dl_psi)(:, :, ist, ik))
            call apply_v(add_hartree1, add_fxc1, nsigma, 1, idir1, ist, ik, freq_factor(ii) * frequency, &
              hvar(:, ispin, :, idir1, swap_sigma(ii)), lr0(:), gpsi(:, :, ist, :))
            do idir2 = 1, ndir
              call apply_v(add_hartree1, add_fxc1, nsigma, 1, idir1, ist, ik, freq_factor(ii) * frequency, &
                hvar(:, ispin, :, idir1, swap_sigma(ii)), lr_k(idir2, :), gdl_k(idir2, :, :, :))
            end do
            do idir2 = 1, ndir
              call apply_v(add_hartree1, add_fxc1, nsigma, 1, idir1, ist, ik, freq_factor(ii) * frequency, &
                hvar(:, ispin, :, idir1, swap_sigma(ii)), lr_b(idir2, :), gdl_b(:, :, :))
              call apply_v(add_hartree1, add_fxc1, nsigma, nsigma, idir1, ist, ik, freq_factor(ii) * frequency, &
                hvar(:, ispin, :, idir1, swap_sigma(ii)), lr_e(idir2, :, ii), gdl_e(:, :, :))
              do idir3 = 1, ndir
                call apply_v(add_hartree1, add_fxc1, nsigma, nsigma, idir1, ist, ik, freq_factor(ii) * frequency, &
                  hvar(:, ispin, :, idir1, swap_sigma(ii)), lr_ke(idir3, idir2, :, ii), gdl_ke(idir3, :, :, :))
              end do
              do idir3 = 1, ndir
                do idim = 1, ndim 
                  zpol_k(idir1, idir3, idir2) = zpol_k(idir1, idir3, idir2) + M_HALF / P_C * factor_e * &
                    (X(mf_dotp)(sys%gr%mesh, lr_e(idir3, nsigma, ii)%X(dl_psi)(:, idim, ist, ik), factor0 * gdl_b(:, idim, 1)) &
                    + X(mf_dotp)(sys%gr%mesh, factor0 * gdl_b(:, idim, nsigma), lr_e(idir3, 1, ii)%X(dl_psi)(:, idim, ist, ik)))
  
                  zpol_k(idir1, idir2, idir3) = zpol_k(idir1, idir2, idir3) + M_HALF / P_C * factor_e * &
                    (X(mf_dotp)(sys%gr%mesh, factor * gdl_e(:, idim, nsigma), lr_b(idir3, 1)%X(dl_psi)(:, idim, ist, ik))&
                    + X(mf_dotp)(sys%gr%mesh, lr_b(idir3, 1)%X(dl_psi)(:, idim, ist, ik), factor * gdl_e(:, idim, 1)))

                  zpol_k(idir1, idir2, idir3) = zpol_k(idir1, idir2, idir3) - M_FOURTH / P_C * &
                    (X(mf_dotp)(sys%gr%mesh, lr_ke(magn_dir(idir3, 1), idir2, nsigma, ii)%X(dl_psi)(:, idim, ist, ik),&
                      factor0 * gdl_k(magn_dir(idir3, 2), :, idim, 1)) &
                    + X(mf_dotp)(sys%gr%mesh, factor * gdl_ke(magn_dir(idir3, 1), :, idim, nsigma),&
                      lr_k(magn_dir(idir3, 2), 1)%X(dl_psi)(:, idim, ist, ik))&
                    + X(mf_dotp)(sys%gr%mesh, -lr_k(magn_dir(idir3, 1), 1)%X(dl_psi)(:, idim, ist, ik),&
                      factor0 * gdl_ke(magn_dir(idir3, 2), :, idim, 1)) &
                    + X(mf_dotp)(sys%gr%mesh, -gdl_k(magn_dir(idir3, 1), :, idim, nsigma),&
                      -factor0 * lr_ke(magn_dir(idir3, 2), idir2, 1, ii)%X(dl_psi)(:, idim, ist, ik))&
                    - X(mf_dotp)(sys%gr%mesh, lr_ke(magn_dir(idir3, 2), idir2, nsigma, ii)%X(dl_psi)(:, idim, ist, ik),&
                      factor0 * gdl_k(magn_dir(idir3, 1), :, idim, 1)) &
                    - X(mf_dotp)(sys%gr%mesh, factor * gdl_ke(magn_dir(idir3, 2), :, idim, nsigma),&
                      lr_k(magn_dir(idir3, 1), 1)%X(dl_psi)(:,idim, ist, ik))&
                    - X(mf_dotp)(sys%gr%mesh, -lr_k(magn_dir(idir3, 2), 1)%X(dl_psi)(:, idim, ist, ik),&
                      factor0 * gdl_ke(magn_dir(idir3, 1), :, idim, 1))&
                    - X(mf_dotp)(sys%gr%mesh, -gdl_k(magn_dir(idir3, 2), :, idim, nsigma),&
                      -factor0 * lr_ke(magn_dir(idir3, 1), idir2, 1, ii)%X(dl_psi)(:, idim, ist, ik)))
                end do
              end do 
            end do
          end if
        end do
        do isigma = 1, nsigma
          call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
            lr0(1)%X(dl_psi)(:, :, :, ik), factor * gpsi(:, :, :, isigma), mat_g(idir1, isigma)%X(matrix))
        end do
      end do
      do idir1 = 1, ndir
        do idir2 = 1, ndir
          do idir3 = 1, ndir  
            call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
              factor0 * lr_k(idir3, 1)%X(dl_psi)(:,:,:, ik), lr_ke(idir2, idir1, 1, ii)%X(dl_psi)(:,:,:, ik),&
              mat_kke(idir3, idir2)%X(matrix))
            call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
              lr_ke(idir3, idir1, nsigma, ii)%X(dl_psi)(:,:,:, ik), factor0 * lr_k(idir2, 1)%X(dl_psi)(:,:,:, ik),&
              mat_kek(idir3, idir2)%X(matrix))
          end do
        end do
        do idir2 = 1, ndir
          call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
            factor1 * lr_b(idir2, 1)%X(dl_psi)(:,:,:, ik), lr_e(idir1, 1, ii)%X(dl_psi)(:,:,:, ik),&
            mat_be%X(matrix))
          call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
            lr_e(idir1, nsigma, ii)%X(dl_psi)(:,:,:, ik), lr_b(idir2, 1)%X(dl_psi)(:,:,:, ik),&
            mat_eb%X(matrix))
          do idir3 = 1, ndir
            do ist = 1, sys%st%nst
              if (abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
                do ist_occ  = sys%st%st_start, sys%st%st_end
                  if (abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then

                    zpol_k(idir3, idir1, idir2) = zpol_k(idir3, idir1, idir2) &  
                      - M_HALF / P_C * factor_e * (factor1 * mat_g(idir3, 1)%X(matrix)(ist_occ, ist)&
                      + R_CONJ(mat_g(idir3, nsigma)%X(matrix)(ist, ist_occ))) * (mat_be%X(matrix)(ist, ist_occ)&
                      + mat_eb%X(matrix)(ist, ist_occ))

                    zpol_k(idir3, idir1, idir2) = zpol_k(idir3, idir1, idir2) &
                      - M_FOURTH / P_C * factor * (factor1 * mat_g(idir3, 1)%X(matrix)(ist_occ, ist)&
                      + R_CONJ(mat_g(idir3, nsigma)%X(matrix)(ist, ist_occ))) * (&
                      mat_kke(magn_dir(idir2, 2),magn_dir(idir2, 1))%X(matrix)(ist, ist_occ)&
                      - mat_kke(magn_dir(idir2, 1),magn_dir(idir2, 2))%X(matrix)(ist, ist_occ)&
                      + mat_kek(magn_dir(idir2, 2),magn_dir(idir2, 1))%X(matrix)(ist, ist_occ)&
                      - mat_kek(magn_dir(idir2, 1),magn_dir(idir2, 2))%X(matrix)(ist, ist_occ))
                  end if
                end do
              end if
            end do
          end do
        end do
      end do
    end do
    do idir1 = 1, ndir
      do idir2 = 1, ndir
        do idir3 = 1, ndir
          zpol(idir3, idir1, idir2) = zpol(idir3, idir1, idir2) + weight * zpol_k(idir3, idir1, idir2)
          if(kpt_output) zpol_kout(idir3, idir1, idir2, ik) = zpol_kout(idir3, idir1, idir2, ik) + &
            zpol_k(idir3, idir1, idir2) * sys%st%smear%el_per_state
        end do
      end do
    end do
  end do
  
  call lr_dealloc(lr0(1))

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**3, MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    zpol(1:ndir, 1:ndir, 1:ndir) = zpol_temp(1:ndir, 1:ndir, 1:ndir)
  endif
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(zpol, zpol_temp, MAX_DIM**3, MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    zpol(1:ndir, 1:ndir, 1:ndir) = zpol_temp(1:ndir, 1:ndir, 1:ndir)
  endif
  if(kpt_output) then
    if(sys%st%parallel_in_states) then
      do ik = 1, sys%st%d%nik
        zpol_temp(:, :, :) = M_ZERO
        call MPI_Allreduce(zpol_kout(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM, ik), zpol_temp, MAX_DIM**3, &
          MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
        do idir1 = 1, ndir
          do idir2 = 1, ndir
            do idir3 = 1, ndir
              zpol_kout(idir1, idir2, idir3, ik) = zpol_temp(idir1, idir2, idir3)
            end do
          end do
        end do
      end do
    endif
    if(sys%st%d%kpt%parallel) then
      do ik = 1, sys%st%d%nik
        zpol_temp(:, :, :) = M_ZERO
        call MPI_Allreduce(zpol_kout(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM, ik), zpol_temp, MAX_DIM**3, &
          MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
        do idir1 = 1, ndir
          do idir2 = 1, ndir
            do idir3 = 1, ndir
              zpol_kout(idir1, idir2, idir3, ik) = zpol_temp(idir1, idir2, idir3)
            end do
          end do
        end do
      end do
    end if
  end if
#endif

  do idir1 = 1, ndir
    do idir2 = 1, ndir
      do idir3 = 1, ndir
        zpol(idir3, idir1, idir2) = zpol(idir3, idir1, idir2) + zpol0(idir3, idir1, idir2)   
      end do
    end do
  end do

  zpol(:,:,:) = -M_zI / (frequency) * zpol(:,:,:) 
  if(nfactor_ke > 1) zpol(:,:,:) = M_HALF * zpol(:,:,:)
  call zsymmetrize_magneto_optics_cart(sys%gr%mesh%sb%symm, zpol(:,:,:))
  
  if(kpt_output) then
    zpol_kout(:,:,:,:) = -M_zI / (frequency) * zpol_kout(:,:,:,:)
    if(nfactor_ke > 1) zpol_kout(:,:,:,:) = M_HALF * zpol_kout(:,:,:,:) 
  end if
  SAFE_DEALLOCATE_P(mat_eb%X(matrix))
  SAFE_DEALLOCATE_P(mat_be%X(matrix))
  do idir2 = 1, ndir
    do isigma = 1, nsigma
      SAFE_DEALLOCATE_P(mat_g(idir2, isigma)%X(matrix))
    end do
    do idir3 = 1, ndir
      SAFE_DEALLOCATE_P(mat_kek(idir2, idir3)%X(matrix))
      SAFE_DEALLOCATE_P(mat_kke(idir2, idir3)%X(matrix))
    end do
  end do
  
  SAFE_DEALLOCATE_A(mat_kke)
  SAFE_DEALLOCATE_A(mat_kek)
  SAFE_DEALLOCATE_A(gpsi)
  SAFE_DEALLOCATE_A(gdl_e)
  SAFE_DEALLOCATE_A(gdl_k)
  SAFE_DEALLOCATE_A(gdl_b)
  SAFE_DEALLOCATE_A(gdl_ke)
  SAFE_DEALLOCATE_A(psi_be)
  SAFE_DEALLOCATE_A(hvar)
  SAFE_DEALLOCATE_A(hvar2)

  POP_SUB(X(lr_calc_magneto_optics_periodic))

contains
  subroutine apply_v(add_hartree, add_fxc, nsigma_h, nsigma_in, dir, ist0, ik0, frequency0, &
    hvar_in, lr_in, psi_out)
    
    logical,      intent(in)    :: add_hartree, add_fxc
    integer,      intent(in)    :: nsigma_h, nsigma_in
    integer,      intent(in)    :: dir, ist0, ik0
    CMPLX,        intent(in)    :: frequency0
    R_TYPE,       intent(inout) :: hvar_in(:, :) 
    type(lr_t),   intent(inout) :: lr_in(:)  
    R_TYPE,       intent(inout) :: psi_out(:, :, :)  
    
    type(pert_t)  :: pert_kdotp
    integer       :: ip, idim

    PUSH_SUB(X(lr_calc_magneto_optics_periodic).apply_v)

    ASSERT(nsigma_h .ge. nsigma_in)

    call pert_init(pert_kdotp, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%geo)
    call pert_setup_dir(pert_kdotp, dir)

    psi_out(:, :, :) = M_ZERO
    do isigma = 1, nsigma_in
      call X(pert_apply)(pert_kdotp, sys%namespace, sys%gr, sys%geo, sys%hm, ik0, &
        lr_in(isigma)%X(dl_psi)(:, :, ist0, ik0), psi_out(:, :, isigma))
    end do
    if((nsigma_h == 2) .and. (nsigma_in == 1)) psi_out(:, :, nsigma_h) = psi_out(:, :, 1)

    if(add_hartree .or. add_fxc) then
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, 1) = psi_out(ip, idim, 1) + frequency0 * &
            factor_e * hvar_in(ip, 1) * lr_in(1)%X(dl_psi)(ip, idim, ist0, ik0)
        end do
      end do
      if(nsigma_h == 2) then 
        do idim = 1, sys%hm%d%dim
          do ip = 1, sys%gr%mesh%np
            psi_out(ip, idim, nsigma_h) = psi_out(ip, idim, nsigma_h) - R_CONJ(frequency0) * &
              factor_e * hvar_in(ip, nsigma_h) * lr_in(nsigma_in)%X(dl_psi)(ip, idim, ist0, ik0)
          end do
        end do
      end if
    end if

    call pert_end(pert_kdotp)

    POP_SUB(X(lr_calc_magneto_optics_periodic).apply_v)
  end subroutine apply_v
end subroutine X(lr_calc_magneto_optics_periodic)

! ---------------------------------------------------------
! See papers Shi et. al Phys. Rev. Lett. 99, 197202 (2007)
! K.-T. Chen and P. A. Lee, Phys. Rev. B 84, 205137 (2011)
subroutine X(lr_calc_magnetization_periodic)(sys, lr_k, magn)
  type(electrons_t),    intent(inout) :: sys 
  type(lr_t),           intent(inout) :: lr_k(:) 
  CMPLX,                intent(out)   :: magn(:)

  integer :: idir1, idir2, idir, ist, idim, ndim, ik, ndir, np
  R_TYPE :: factor
  R_TYPE, allocatable :: Hdl_psi(:,:,:)
  FLOAT :: weight
  
#ifdef HAVE_MPI
  CMPLX :: magn_temp(1:MAX_DIM)
#endif

  PUSH_SUB(X(lr_calc_magnetization_periodic))
  
#if defined(R_TCOMPLEX)
  factor = -M_zI
#else
  factor = M_ONE
#endif
  
  np = sys%gr%mesh%np
  ndir = sys%gr%mesh%sb%dim
  ndim = sys%st%d%dim

  SAFE_ALLOCATE(Hdl_psi(1:np, 1:ndim, 1:ndir))
  
  magn(:) = M_ZERO
  
  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    weight = sys%st%d%kweights(ik) * sys%st%smear%el_per_state
    do ist = 1, sys%st%nst
      if (abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
        do idir = 1, ndir
          call X(hamiltonian_elec_apply_single)(sys%hm, sys%namespace, sys%gr%mesh, lr_k(idir)%X(dl_psi)(:, :, ist, ik), &
            Hdl_psi(:, :, idir), ist, ik)
        end do
   
        do idir = 1, ndir
          idir1 = magn_dir(idir, 1)
          idir2 = magn_dir(idir, 2)
          do idim = 1, ndim
            magn(idir) = magn(idir) + M_zI * M_HALF * M_HALF * weight / P_C * &
              (X(mf_dotp)(sys%gr%mesh, lr_k(idir2)%X(dl_psi)(1:np, idim, ist, ik), Hdl_psi(1:np, idim, idir1)) &
              - X(mf_dotp)(sys%gr%mesh, lr_k(idir1)%X(dl_psi)(1:np, idim, ist, ik), Hdl_psi(1:np, idim, idir2)) &
              + X(mf_dotp)(sys%gr%mesh, Hdl_psi(1:np, idim, idir2), lr_k(idir1)%X(dl_psi)(1:np, idim, ist, ik)) &
              - X(mf_dotp)(sys%gr%mesh, Hdl_psi(1:np, idim, idir1), lr_k(idir2)%X(dl_psi)(1:np, idim, ist, ik)))

            magn(idir) = magn(idir) + M_zI * M_HALF * weight * sys%st%eigenval(ist, ik) / P_C * (&
              X(mf_dotp)(sys%gr%mesh, lr_k(idir2)%X(dl_psi)(1:np, idim, ist, ik), &
                lr_k(idir1)%X(dl_psi)(1:np, idim, ist, ik))&
              - X(mf_dotp)(sys%gr%mesh, lr_k(idir1)%X(dl_psi)(1:np, idim, ist, ik), &
                lr_k(idir2)%X(dl_psi)(1:np, idim, ist, ik))) 
          end do
        end do
      end if
    end do
  end do

  SAFE_DEALLOCATE_A(Hdl_psi)

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(magn, magn_temp, MAX_DIM, MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    magn(1:ndir) = magn_temp(1:ndir)
  endif
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(magn, magn_temp, MAX_DIM, MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    magn(1:ndir) = magn_temp(1:ndir)
  endif
#endif

  POP_SUB(X(lr_calc_magnetization_periodic))
end subroutine X(lr_calc_magnetization_periodic)

!--------------------------------------------------------
! According to paper X. Gonze and J. W. Zwanziger, Phys. Rev. B 84, 064445 (2011)
subroutine X(lr_calc_susceptibility_periodic)(sys, nsigma, lr_k, lr_b, &
  lr_kk, lr_kb, magn)
  type(electrons_t),    intent(inout) :: sys 
  integer,              intent(in)    :: nsigma
  type(lr_t),           intent(inout) :: lr_k(:) 
  type(lr_t),           intent(inout) :: lr_b(:) 
  type(lr_t),           intent(inout) :: lr_kk(:,:)
  type(lr_t),           intent(inout) :: lr_kb(:,:) 
  CMPLX,                intent(out)   :: magn(:,:)

  integer :: idir1, idir2, idir, ist, ispin, idim, ndim 
  integer :: dir1, dir2, dir3, dir4, np, ik, ist_occ, ndir

  FLOAT :: weight
  R_TYPE:: factor, factor0
  R_TYPE, allocatable :: Hdl_k(:,:,:,:), Hdl_b(:,:,:) 
  R_TYPE, allocatable :: Hdl_kb(:,:,:,:), Hdl_kk(:,:,:,:)
  type(matrix_t), allocatable :: kH_mat(:,:), Hk_mat(:,:), kk_mat(:,:)

#ifdef HAVE_MPI
  CMPLX :: magn_temp(1:MAX_DIM,1:MAX_DIM)
#endif 

#if defined(R_TCOMPLEX)
  factor = -M_zI
  factor0 = -M_zI
#else
  factor = -M_ONE
  factor0 = M_ONE
#endif
  
  np = sys%gr%mesh%np
  ndir = sys%gr%mesh%sb%dim
  ndim = sys%st%d%dim
  
  PUSH_SUB(X(lr_calc_susceptibility_periodic))

  SAFE_ALLOCATE(Hdl_b(1:np, 1:ndim, 1:ndir))
  SAFE_ALLOCATE(Hdl_kb(1:np, 1:ndim, 1:ndir, 1:ndir))
  SAFE_ALLOCATE(Hdl_k(1:np, 1:ndim, 1:sys%st%nst, 1:ndir))
  SAFE_ALLOCATE(Hdl_kk(1:np, 1:ndim, 1:ndir, 1:ndir))
  SAFE_ALLOCATE(Hk_mat(1:ndir, 1:ndir))
  SAFE_ALLOCATE(kH_mat(1:ndir, 1:ndir))
  SAFE_ALLOCATE(kk_mat(1:ndir, 1:ndir))
  
  do idir1 = 1, ndir
    do idir2 = 1, ndir
      SAFE_ALLOCATE(Hk_mat(idir1, idir2)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
      SAFE_ALLOCATE(kH_mat(idir1, idir2)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
      SAFE_ALLOCATE(kk_mat(idir1, idir2)%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
    end do
  end do

  
  magn(:,:) = M_ZERO  
  
  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    weight = sys%st%d%kweights(ik)*sys%st%smear%el_per_state
    ispin = states_elec_dim_get_spin_index(sys%st%d,ik)
    do ist = 1, sys%st%nst
      if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
        do idir = 1, ndir
          call X(hamiltonian_elec_apply_single)(sys%hm, sys%namespace, sys%gr%mesh, lr_b(idir)%X(dl_psi)(:,:,ist,ik), &
            Hdl_b(:,:,idir), ist, ik)
          call X(hamiltonian_elec_apply_single)(sys%hm, sys%namespace, sys%gr%mesh, lr_k(idir)%X(dl_psi)(:,:,ist,ik), &
            Hdl_k(:,:,ist,idir), ist, ik)
        end do

        do idir1 = 1, ndir
          do idir2 = 1, ndir
            call X(hamiltonian_elec_apply_single)(sys%hm, sys%namespace, sys%gr%mesh, &
              lr_kb(idir1, idir2)%X(dl_psi)(:,:, ist, ik), Hdl_kb(:,:, idir1, idir2), ist, ik)
            call X(hamiltonian_elec_apply_single)(sys%hm, sys%namespace, sys%gr%mesh, &
              lr_kk(max(idir1, idir2), min(idir1, idir2))%X(dl_psi)(:,:, ist, ik), Hdl_kk(:,:, idir1, idir2), ist, ik)
          end do
        end do
        do idir1 = 1, ndir
          do idir2 = 1, ndir
            do idim = 1, ndim

              magn(idir1, idir2) = magn(idir1, idir2) + M_HALF * weight / (P_C**2) * (&
                X(mf_dotp)(sys%gr%mesh, lr_b(idir2)%X(dl_psi)(:, idim, ist, ik), Hdl_b(:, idim, idir1))&
                + X(mf_dotp)(sys%gr%mesh, lr_b(idir1)%X(dl_psi)(:, idim, ist, ik), Hdl_b(:, idim, idir2))&
                + X(mf_dotp)(sys%gr%mesh, Hdl_b(:,idim, idir2), lr_b(idir1)%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, Hdl_b(:,idim, idir1), lr_b(idir2)%X(dl_psi)(:, idim, ist, ik)))
        
              magn(idir1, idir2) = magn(idir1, idir2) - weight * sys%st%eigenval(ist, ik) / (P_C**2) * (&
                X(mf_dotp)(sys%gr%mesh, lr_b(idir2)%X(dl_psi)(:, idim, ist, ik), &
                  lr_b(idir1)%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, lr_b(idir1)%X(dl_psi)(:, idim, ist, ik), &
                  lr_b(idir2)%X(dl_psi)(:, idim, ist, ik)))
  
              magn(idir1, idir2) = magn(idir1, idir2) - M_HALF * factor0 / (P_C**2) * weight * sys%st%eigenval(ist, ik) * (&
                X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir2, 2))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kb(magn_dir(idir2, 1), idir1)%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir2, 1))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kb(magn_dir(idir2, 2),idir1)%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir1, 2))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kb(magn_dir(idir1, 1), idir2)%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir1, 1))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kb(magn_dir(idir1, 2), idir2)%X(dl_psi)(:, idim, ist, ik))& 
                + X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir2, 2), idir1)%X(dl_psi)(:, idim, ist, ik), &
                  factor * lr_k(magn_dir(idir2, 1))%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir2, 1), idir1)%X(dl_psi)(:, idim, ist, ik), &
                  factor * lr_k(magn_dir(idir2, 2))%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir1, 2), idir2)%X(dl_psi)(:, idim, ist, ik), &
                  factor * lr_k(magn_dir(idir1, 1))%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir1, 1), idir2)%X(dl_psi)(:, idim, ist, ik), &
                  factor * lr_k(magn_dir(idir1, 2))%X(dl_psi)(:, idim, ist, ik)))
 
              magn(idir1, idir2) = magn(idir1, idir2) + M_FOURTH * weight / (P_C**2) * (&
                X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir2, 1))%X(dl_psi)(:, idim, ist, ik), &
                  factor0 * Hdl_kb(:, idim, magn_dir(idir2, 2), idir1))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir2, 2))%X(dl_psi)(:, idim, ist, ik), &
                  factor0 * Hdl_kb(:, idim, magn_dir(idir2, 1), idir1))&
                + X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir1, 1))%X(dl_psi)(:, idim, ist, ik), &
                  factor0 * Hdl_kb(:, idim, magn_dir(idir1, 2), idir2))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * lr_k(magn_dir(idir1, 2))%X(dl_psi)(:, idim, ist, ik), &
                  factor0 * Hdl_kb(:, idim, magn_dir(idir1, 1), idir2))&   
                + X(mf_dotp)(sys%gr%mesh, factor0 * Hdl_k(:, idim, ist, magn_dir(idir2, 1)), &
                  factor0 * lr_kb(magn_dir(idir2, 2), idir1)%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * Hdl_k(:, idim, ist, magn_dir(idir2, 2)), &
                  factor0 * lr_kb(magn_dir(idir2, 1), idir1)%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, factor0 * Hdl_k(:, idim, ist, magn_dir(idir1, 1)), &
                  factor0 * lr_kb(magn_dir(idir1, 2), idir2)%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, factor0 * Hdl_k(:, idim, ist, magn_dir(idir1, 2)), &
                  factor0 * lr_kb(magn_dir(idir1, 1), idir2)%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir2, 1), idir1)%X(dl_psi)(:, idim, ist, ik),&
                  -Hdl_k(:, idim, ist, magn_dir(idir2, 2)))&
                - X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir2, 2), idir1)%X(dl_psi)(:, idim, ist, ik),&
                  -Hdl_k(:, idim, ist, magn_dir(idir2, 1)))&
                + X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir1, 1), idir2)%X(dl_psi)(:, idim, ist, ik),&
                  -Hdl_k(:, idim, ist, magn_dir(idir1, 2)))&
                - X(mf_dotp)(sys%gr%mesh, lr_kb(magn_dir(idir1, 2), idir2)%X(dl_psi)(:, idim, ist, ik),&
                  -Hdl_k(:, idim, ist, magn_dir(idir1, 1)))&  
                + X(mf_dotp)(sys%gr%mesh, Hdl_kb(:, idim, magn_dir(idir2, 1), idir1), &
                  -lr_k(magn_dir(idir2, 2))%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, Hdl_kb(:, idim, magn_dir(idir2, 2), idir1), &
                  -lr_k(magn_dir(idir2, 1))%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, Hdl_kb(:, idim, magn_dir(idir1, 1), idir2), &
                  -lr_k(magn_dir(idir1, 2))%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, Hdl_kb(:, idim, magn_dir(idir1, 2), idir2), &
                  -lr_k(magn_dir(idir1, 1))%X(dl_psi)(:, idim, ist, ik)))
        
              magn(idir1, idir2) = magn(idir1, idir2) - M_FOURTH * M_HALF * weight / (P_C**2) * (&
                X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik), &
                  Hdl_kk(:, idim, magn_dir(idir1,2), magn_dir(idir2,2))) &
                - X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik), &
                  Hdl_kk(:, idim, magn_dir(idir1, 1), magn_dir(idir2, 2))) &
                - X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik), &
                  Hdl_kk(:, idim, magn_dir(idir1, 2), magn_dir(idir2, 1))) &
                + X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik), &
                  Hdl_kk(:, idim, magn_dir(idir1, 1), magn_dir(idir2, 1))) &
                + X(mf_dotp)(sys%gr%mesh, Hdl_kk(:, idim, magn_dir(idir1, 1), magn_dir(idir2, 1)), &
                  lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik)) &
                - X(mf_dotp)(sys%gr%mesh, Hdl_kk(:, idim, magn_dir(idir1, 2), magn_dir(idir2, 1)),&
                  lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 2)), &
                   min(magn_dir(idir1, 1), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik))&
                - X(mf_dotp)(sys%gr%mesh, Hdl_kk(:, idim, magn_dir(idir1, 1), magn_dir(idir2, 2)),&
                  lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, Hdl_kk(:, idim, magn_dir(idir1, 2), magn_dir(idir2, 2)),&
                  lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik)))

              magn(idir1, idir2) = magn(idir1, idir2) + M_FOURTH / (P_C**2) * weight * sys%st%eigenval(ist, ik) * (&
                + X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 2) ,magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik)) &
                - X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik)) &
                - X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik))&
                + X(mf_dotp)(sys%gr%mesh, lr_kk(max(magn_dir(idir1, 1), magn_dir(idir2, 1)), &
                  min(magn_dir(idir1, 1), magn_dir(idir2, 1)))%X(dl_psi)(:, idim, ist, ik), &
                  lr_kk(max(magn_dir(idir1, 2), magn_dir(idir2, 2)), &
                  min(magn_dir(idir1, 2), magn_dir(idir2, 2)))%X(dl_psi)(:, idim, ist, ik)))
            end do
          end do
        end do
      end if
    end do

    do idir1 = 1, ndir
      do idir2 = 1, ndir
        call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
          lr_k(idir1)%X(dl_psi)(:,:,:, ik), Hdl_k(:,:,:, idir2), kH_mat(idir1, idir2)%X(matrix))
    
        call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
          Hdl_k(:,:,:, idir1), lr_k(idir2)%X(dl_psi)(:,:,:, ik), Hk_mat(idir1, idir2)%X(matrix))
    
        call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
          lr_k(idir1)%X(dl_psi)(:,:,:, ik), lr_k(idir2)%X(dl_psi)(:,:,:, ik), kk_mat(idir1, idir2)%X(matrix))
      end do
    end do
    do idir1 = 1, ndir
      do idir2 = 1, ndir
        dir1 = magn_dir(idir1, 1)
        dir2 = magn_dir(idir2, 1)
        dir3 = magn_dir(idir1, 2)
        dir4 = magn_dir(idir2, 2)

        do ist = 1, sys%st%nst
          if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
            do ist_occ = 1, sys%st%nst
              if (abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then
                magn(idir1, idir2) = magn(idir1, idir2) - M_FOURTH * M_HALF / (P_C**2) * weight * (&
                  kH_mat(dir4, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir3)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir4, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir1)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir2, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir3)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir2, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir1)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir3, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir3, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir1, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir1, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir4, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir4, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir2, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir2, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir3, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir4)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir1, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir4)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir3, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir2)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir1, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir2)%X(matrix)(ist, ist_occ))

                magn(idir1, idir2) = magn(idir1, idir2) - M_FOURTH * M_HALF / (P_C**2) * weight * (&
                  Hk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir3)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir1)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir2, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir3)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir2, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir1)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir3, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir3, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir1, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir1, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir4, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir4, dir2)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir2, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir2, dir4)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir3, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir4)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir1, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir2, dir4)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir3, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir2)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir1, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir4, dir2)%X(matrix)(ist, ist_occ))

                magn(idir1, idir2) = magn(idir1, idir2) + M_FOURTH * sys%st%eigenval(ist, ik) / (P_C**2) * weight * (&
                  kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir3, dir4)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir1, dir4)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir3, dir2)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir1, dir2)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir2, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir4, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir3)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir2, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir4, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir1)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir2, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir3, dir4)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir4, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir3, dir2)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir2, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir1, dir4)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir4, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir1, dir2)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir3)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir1)%X(matrix)(ist_occ, ist))
  
                magn(idir1, idir2) = magn(idir1, idir2) - M_HALF * M_FOURTH / (P_C**2) * weight * (&
                  kH_mat(dir4, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir4, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) &
                  - kH_mat(dir2, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) &
                  + kH_mat(dir2, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) &
                  - Hk_mat(dir2, dir1)%X(matrix)(ist_occ, ist) * kk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) &
                  + Hk_mat(dir2, dir3)%X(matrix)(ist_occ, ist) * kk_mat(dir1, dir4)%X(matrix)(ist, ist_occ))

                magn(idir1, idir2) = magn(idir1, idir2) - M_HALF * M_FOURTH / (P_C**2) * weight * (&
                  kH_mat(dir3, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) &
                  - kH_mat(dir3, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir1)%X(matrix)(ist_occ, ist) &
                  - kH_mat(dir1, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) &
                  + kH_mat(dir1, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir3)%X(matrix)(ist_occ, ist) &
                  + Hk_mat(dir3, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir1)%X(matrix)(ist_occ, ist) &
                  - Hk_mat(dir3, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir1)%X(matrix)(ist_occ, ist) &
                  - Hk_mat(dir1, dir2)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir3)%X(matrix)(ist_occ, ist) &
                  + Hk_mat(dir1, dir4)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir3)%X(matrix)(ist_occ, ist))

                magn(idir1, idir2) = magn(idir1, idir2) + M_FOURTH * (&
                  sys%st%eigenval(ist, ik) + sys%st%eigenval(ist_occ, ik)) / (P_C**2) * weight * (&
                  kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir4)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir2, dir4)%X(matrix)(ist_occ, ist) &
                  - kk_mat(dir1, dir3)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir2)%X(matrix)(ist_occ, ist) &
                  + kk_mat(dir3, dir1)%X(matrix)(ist, ist_occ) * kk_mat(dir4, dir2)%X(matrix)(ist_occ, ist))
              end if 
            end do
          end if
        end do
      end do
    end do
  end do
  
  magn(:,:) = -magn(:,:)
  
  do idir1 = 1, ndir
    do idir2 = 1, ndir
      SAFE_DEALLOCATE_P(Hk_mat(idir1, idir2)%X(matrix))
      SAFE_DEALLOCATE_P(kH_mat(idir1, idir2)%X(matrix))
      SAFE_DEALLOCATE_P(kk_mat(idir1, idir2)%X(matrix))
    end do
  end do

  SAFE_DEALLOCATE_A(Hdl_k)
  SAFE_DEALLOCATE_A(Hdl_b)
  SAFE_DEALLOCATE_A(Hdl_kb)
  SAFE_DEALLOCATE_A(Hdl_kk) 

  call zsymmetrize_tensor_cart(sys%gr%mesh%sb%symm, magn(:,:))

#ifdef HAVE_MPI
  if(sys%st%parallel_in_states) then
    call MPI_Allreduce(magn, magn_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%mpi_grp%comm, mpi_err)
    magn(1:ndir, 1:ndir) = magn_temp(1:ndir, 1:ndir)
  endif
  if(sys%st%d%kpt%parallel) then
    call MPI_Allreduce(magn, magn_temp, MAX_DIM**2, MPI_CMPLX, MPI_SUM, sys%st%d%kpt%mpi_grp%comm, mpi_err)
    magn(1:ndir, 1:ndir) = magn_temp(1:ndir, 1:ndir)
  endif
#endif
  
  POP_SUB(X(lr_calc_susceptibility_periodic))
end subroutine X(lr_calc_susceptibility_periodic)

! ---------------------------------------------------------
!  -Pc{V1!dl_dk2>+!dn_dk2><n!V1!l>}, n,l=occ   
! This subroutine can be used for setting the right-hand side of
! the Sternheimer equation for magnetic and second-order kdotp 
! perturbations according to the density matrix formulation.
subroutine X(inhomog_per_component)(sys, idir, ik, & 
  psi_k2, psi_out, factor_tot, factor_k, factor_second)
  type(electrons_t),    intent(inout) :: sys 
  integer,              intent(in)    :: idir
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: psi_k2(:,:,:)    
  R_TYPE,               intent(inout) :: psi_out(:,:,:)  
  R_TYPE,               intent(in)    :: factor_tot, factor_k, factor_second
  
  R_TYPE, allocatable :: f_out(:,:), vel(:,:,:)
  R_TYPE :: factor
  integer :: ip, ist, idim, ist_occ
  type(matrix_t):: vel_mat
  type(pert_t)  :: pert_kdotp
  R_TYPE, allocatable :: psi(:,:,:)

  PUSH_SUB(X(inhomog_per_component))
  
  SAFE_ALLOCATE(f_out(1:sys%gr%mesh%np,1:sys%hm%d%dim))
  SAFE_ALLOCATE(vel(1:sys%gr%mesh%np,1:sys%hm%d%dim, 1:sys%st%nst))
  SAFE_ALLOCATE(vel_mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim, 1:sys%st%nst))
  
#if defined(R_TCOMPLEX)
  factor = -M_zI
#else
  factor = M_ONE
#endif

  f_out(:,:) = M_ZERO
  
  call pert_init(pert_kdotp, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%geo)
  call pert_setup_dir(pert_kdotp, idir)
  
  vel(:,:,:) = M_ZERO
  psi(:,:,:) = M_ZERO
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi(:,:,ist))
      call X(pert_apply)(pert_kdotp, sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi(:,:,ist), f_out) 
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          vel(ip,idim,ist) = factor * f_out(ip,idim)
        end do
      end do
    end if
  end do  
    
  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    psi(:,:,:), vel(:,:,:), vel_mat%X(matrix))

  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then

      call X(pert_apply)(pert_kdotp, sys%namespace, sys%gr, sys%geo, sys%hm, ik, factor_k * psi_k2(:, :, ist), f_out)
        
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, ist) = psi_out(ip, idim, ist) - &
            factor_tot * factor * f_out(ip, idim)
        end do
      end do

      do ist_occ = 1, sys%st%nst
        if(abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then 
           do idim = 1, sys%hm%d%dim
            do ip = 1, sys%gr%mesh%np
              psi_out(ip, idim, ist) = psi_out(ip, idim, ist) - factor_second * factor_tot *&
                factor_k * psi_k2(ip, idim, ist_occ) * vel_mat%X(matrix)(ist_occ, ist)   
            end do
          end do
        end if
      end do 
    end if
  end do 
  
  call pert_end(pert_kdotp)
  SAFE_DEALLOCATE_P(vel_mat%X(matrix))
 
  SAFE_DEALLOCATE_A(f_out)
  SAFE_DEALLOCATE_A(vel)
  SAFE_DEALLOCATE_A(psi)
  POP_SUB(X(inhomog_per_component))
end subroutine X(inhomog_per_component)

!--------------------------------------------------------------------------
!  -Pc(!dn_dk2><dn_de!V1!l>-V1!m><dm_e!dl_k2>)
! This subroutine can be used for setting the right-hand side of
! the Sternheimer equation for second-order magnetic
! perturbations according to the density matrix formulation.  
subroutine X(inhomog_per_component_2nd_order)(sys, idir, ik, & 
  psi_k2, psi_e, psi_out, factor_tot, factor_k, factor_e)
  type(electrons_t),    intent(inout) :: sys 
  integer,              intent(in)    :: idir
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: psi_k2(:,:,:)   
  R_TYPE,               intent(inout) :: psi_e(:,:,:)
  R_TYPE,               intent(inout) :: psi_out(:,:,:) 
  R_TYPE,               intent(in)    :: factor_tot, factor_k, factor_e
  
  R_TYPE, allocatable :: f_out(:,:)
  R_TYPE, allocatable :: vel(:,:,:)
  R_TYPE :: factor
  integer :: ip, ist, idim, ist_occ
  type(matrix_t):: vel_mat, prod_mat
  type(pert_t)  :: pert_kdotp
  R_TYPE, allocatable :: psi(:,:)

  PUSH_SUB(X(inhomog_per_component_2nd_order))

  SAFE_ALLOCATE(f_out(1:sys%gr%mesh%np, 1:sys%hm%d%dim))
  SAFE_ALLOCATE(vel(1:sys%gr%mesh%np,1:sys%hm%d%dim, 1:sys%st%nst))
  SAFE_ALLOCATE(vel_mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(prod_mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
    
#if defined(R_TCOMPLEX)
  factor = -M_zI
#else
  factor = M_ONE
#endif

  f_out(:,:) = M_ZERO
  
  call pert_init(pert_kdotp, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%geo)
  call pert_setup_dir(pert_kdotp, idir)

  vel(:,:,:) = M_ZERO
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
      call X(pert_apply)(pert_kdotp, sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi, f_out)
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          vel(ip, idim, ist) = factor * f_out(ip, idim)
        end do
      end do
    end if
  end do  

  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    factor_e * psi_e(:,:,:), vel(:,:,:), vel_mat%X(matrix))

  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    factor_e * psi_e(:,:,:), factor_k * psi_k2(:,:,:), prod_mat%X(matrix))

  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      do ist_occ = 1, sys%st%nst
        if(abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then
          do idim = 1, sys%hm%d%dim
            do ip = 1, sys%gr%mesh%np
              psi_out(ip, idim, ist) = psi_out(ip, idim, ist) - factor_tot * &
                vel_mat%X(matrix)(ist_occ, ist) * factor_k * psi_k2(ip, idim, ist_occ)
              psi_out(ip, idim, ist) = psi_out(ip, idim, ist) + factor_tot * &
                prod_mat%X(matrix)(ist_occ, ist) * vel(ip, idim, ist_occ)
            end do
          end do
        end if
      end do
    end if
  end do

  call pert_end(pert_kdotp)

  SAFE_DEALLOCATE_P(vel_mat%X(matrix))
  SAFE_DEALLOCATE_P(prod_mat%X(matrix))
  SAFE_DEALLOCATE_A(f_out)
  SAFE_DEALLOCATE_A(vel)
  SAFE_DEALLOCATE_A(psi)
  POP_SUB(X(inhomog_per_component_2nd_order))
end subroutine X(inhomog_per_component_2nd_order)


  ! --------------------------------------------------------------------------
subroutine X(inhomog_B)(sh, sys, idir1, idir2, & 
  lr_k1, lr_k2, psi_out)
  type(sternheimer_t),  intent(inout) :: sh
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir1
  integer,              intent(in)    :: idir2
  type(lr_t),           intent(inout) :: lr_k1(:) 
  type(lr_t),           intent(inout) :: lr_k2(:)   
  R_TYPE,               intent(inout) :: psi_out(:,:,:,:,:)
    
  R_TYPE :: factor_plus, factor_minus, factor_k, &
            factor_magn, factor_rho, factor_sum
  type(lr_t) :: lr0(1)
  R_TYPE, allocatable :: hvar(:,:,:)
  integer :: ispin, ik, ik0
  
  PUSH_SUB(X(inhomog_B))
  
#if defined(R_TCOMPLEX)
  factor_rho = M_zI * M_HALF
  factor_plus = M_HALF * M_zI
  factor_k = -M_zI
#else
  factor_rho = M_HALF
  factor_plus = M_HALF
  factor_k = -M_ONE
#endif
  factor_minus = -factor_plus
  factor_magn = M_ONE

  psi_out(:,:,:,:,:) = M_ZERO
 
  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ik0 = ik-sys%st%d%kpt%start+1
    call X(inhomog_per_component)(sys, idir1, ik, lr_k2(1)%X(dl_psi)(:, :, :, ik), &
      psi_out(:, :, :, ik0, 1), factor_plus, factor_k, factor_magn)
    call X(inhomog_per_component)(sys, idir2, ik, lr_k1(1)%X(dl_psi)(:, :, :, ik), &
      psi_out(:, :, :, ik0, 1), factor_minus, factor_k, factor_magn)
  end do
  if(sternheimer_add_hartree(sh) .or. sternheimer_add_fxc(sh)) then
    SAFE_ALLOCATE(hvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))
    call lr_init(lr0(1))
    call lr_allocate(lr0(1), sys%st, sys%gr%mesh, allocate_rho = .true.)
    lr0(1)%X(dl_rho)(1:sys%gr%mesh%np, 1:sys%st%d%nspin) = M_ZERO
    factor_sum = -M_ONE
    call X(calc_rho)(sys, factor_rho, factor_sum, factor_k, factor_k,&
      lr_k1(1),lr_k2(1),lr0(1))
    call X(calc_rho)(sys, factor_sum * factor_rho, factor_sum, factor_k,&
      factor_k,lr_k2(1),lr_k1(1),lr0(1))
    if(sys%st%parallel_in_states .or. sys%st%d%kpt%parallel) &
      call comm_allreduce(sys%st%st_kpt_mpi_grp%comm, lr0(1)%X(dl_rho))
    call X(sternheimer_calc_hvar)(sh, sys, lr0(1:1), 1, hvar) 
    call lr_dealloc(lr0(1))

    do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
      ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
      ik0 = ik-sys%st%d%kpt%start+1
      call X(calc_hvar_psi)(sys, ik, hvar(:, ispin, 1), psi_out(:, :, :, ik0, 1))
    end do
    SAFE_DEALLOCATE_A(hvar)
  end if

  POP_SUB(X(inhomog_B))
end subroutine X(inhomog_B)

! --------------------------------------------------------------------------
subroutine X(inhomog_EB)(sys, ik, add_hartree, add_fxc, & 
  hvar, psi_b, psi_kb, factor_b, psi_out, psi_k1, psi_k2)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: ik
  logical,              intent(in)    :: add_hartree, add_fxc
  R_TYPE,               intent(inout) :: hvar(:)
  R_TYPE,               intent(inout) :: psi_b(:,:,:)
  R_TYPE,               intent(inout) :: psi_kb(:,:,:)  
  R_TYPE,               intent(in)    :: factor_b 
  R_TYPE,               intent(inout) :: psi_out(:,:,:) 
  R_TYPE,     optional, intent(inout) :: psi_k1(:,:,:), psi_k2(:,:,:)
  
  R_TYPE :: factor, factor_e, factor_sum
  integer :: ip, idim, ist

  PUSH_SUB(X(inhomog_EB))
  
#if defined(R_TCOMPLEX)
  factor = M_zI
#else
  factor = -M_ONE
#endif
  factor_e = -M_ONE

  do ist = 1, sys%st%nst
    do idim = 1, sys%hm%d%dim
      do ip = 1, sys%gr%mesh%np
        psi_out(ip, idim, ist) = psi_out(ip, idim, ist) + factor * psi_kb(ip, idim, ist)
      end do
    end do
  end do
  
  if(add_hartree .or. add_fxc) then
    call X(calc_hvar_lr)(sys, ik, hvar(:), psi_b, factor_e, factor_b, psi_out) 
    if((present(psi_k1)) .and. (present(psi_k2))) then
      factor_sum = M_ONE
      call calc_hvar_lr2(psi_k1, psi_k2, hvar, factor_sum)
      factor_sum = -M_ONE
      call calc_hvar_lr2(psi_k2, psi_k1, hvar, factor_sum)
    end if 
  end if 
 
  POP_SUB(X(inhomog_EB))
  
contains
  subroutine calc_hvar_lr2(tlr_k1, tlr_k2, hvar_e, factor_tot)
    R_TYPE,       intent(inout) :: tlr_k1(:,:,:), tlr_k2(:,:,:)
    R_TYPE,       intent(inout) :: hvar_e(:)
    R_TYPE,       intent(in)    :: factor_tot
    
    R_TYPE, allocatable :: psi(:,:,:), psi0(:,:)
    integer :: ist1
    R_TYPE :: factor0, factor_k
    type(matrix_t):: prod_mat2, prod_mat
    
    PUSH_SUB(X(inhomog_EB).calc_hvar_lr2)
  
#if defined(R_TCOMPLEX)
  factor0 = M_zI * M_HALF
  factor_k = -M_zI
#else
  factor0 = -M_HALF
  factor_k = M_ONE
#endif

    SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, sys%hm%d%dim, 1:sys%st%nst))
    SAFE_ALLOCATE(prod_mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
    SAFE_ALLOCATE(prod_mat2%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
    SAFE_ALLOCATE(psi0(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  
    psi(:,:,:) = M_ZERO
    do ist = 1, sys%st%nst
      if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
        call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi0)
        do idim = 1, sys%hm%d%dim
          do ip = 1, sys%gr%mesh%np
            psi(ip, idim, ist) = factor_e * hvar_e(ip) * psi0(ip,idim)
          end do
        end do
      end if
    end do
      
    call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
      factor_k * tlr_k1(:,:,:), factor_k * tlr_k2(:,:,:), prod_mat%X(matrix))

    call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
      factor_k * tlr_k1(:,:,:), psi(:,:,:), prod_mat2%X(matrix))
  
    do ist1 = 1, sys%st%nst
      if(abs(sys%st%occ(ist1, ik)) .gt. M_EPSILON) then
        call states_elec_get_state(sys%st, sys%gr%mesh, ist1, ik, psi0)
        do ist  = 1, sys%st%nst
          if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
            do idim = 1, sys%hm%d%dim
              do ip = 1, sys%gr%mesh%np 
                psi_out(ip, idim, ist) = psi_out(ip, idim, ist) + factor_tot * factor0 * (&
                  prod_mat%X(matrix)(ist1, ist) * factor_e * hvar_e(ip) * psi0(ip, idim) - &
                  prod_mat2%X(matrix)(ist1, ist) * factor_k * tlr_k2(ip, idim, ist1))                  
              end do
            end do
          end if 
        end do
      end if
    end do
  
    SAFE_DEALLOCATE_P(prod_mat%X(matrix))
    SAFE_DEALLOCATE_P(prod_mat2%X(matrix))
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(psi0)
  
    POP_SUB(X(inhomog_EB).calc_hvar_lr2)
  end subroutine calc_hvar_lr2
end subroutine X(inhomog_EB)

! --------------------------------------------------------------------------
subroutine X(inhomog_BE)(sys, idir1, idir2, ik, & 
  add_hartree, add_fxc, hvar, psi_e1, psi_e2, psi_ek1, &
  psi_ek2, psi_k1, psi_k2, factor_e1, factor_e2, psi_out)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir1, idir2
  integer,              intent(in)    :: ik
  logical,              intent(in)    :: add_hartree, add_fxc
  R_TYPE,               intent(inout) :: hvar(:)
  R_TYPE,               intent(inout) :: psi_e1(:,:,:), psi_e2(:,:,:)       
  R_TYPE,               intent(inout) :: psi_ek1(:,:,:), psi_ek2(:,:,:)  
  R_TYPE,               intent(inout) :: psi_k1(:,:,:) 
  R_TYPE,               intent(inout) :: psi_k2(:,:,:)   
  R_TYPE,               intent(in)    :: factor_e1, factor_e2
  R_TYPE,               intent(inout) :: psi_out(:,:,:) 
  
  R_TYPE :: factor_plus, factor_minus, factor_magn
  R_TYPE :: factor_k, factor_k0, factor_b
  
  PUSH_SUB(X(inhomog_BE))

#if defined(R_TCOMPLEX)
  factor_plus = M_HALF * M_zI
  factor_minus = -M_HALF * M_zI
  factor_k = -M_zI
  factor_k0 = -M_zI
  factor_b = M_ONE
#else
  factor_plus = M_HALF
  factor_minus = -M_HALF
  factor_k = M_ONE
  factor_k0 = -M_ONE
  factor_b = -M_ONE
#endif
  factor_magn = M_ONE

  call X(inhomog_per_component)(sys, idir1, ik, & 
    psi_ek2, psi_out, factor_plus, factor_magn, factor_magn)
  call X(inhomog_per_component)(sys, idir2, ik,& 
    psi_ek1, psi_out, factor_minus, factor_magn, factor_magn)
    
  call X(inhomog_per_component_2nd_order)(sys, idir1, ik, & 
    psi_k2, psi_e2, psi_out, factor_plus, factor_k, factor_e1)
  call X(inhomog_per_component_2nd_order)(sys, idir2, ik, &
    psi_k1, psi_e2, psi_out, factor_minus, factor_k, factor_e1)
  
  call X(inhomog_per_component_2nd_order)(sys, idir1, ik, &
    psi_e1, psi_k2, psi_out, factor_plus, factor_e2, factor_k0)
  call X(inhomog_per_component_2nd_order)(sys, idir2, ik, &
    psi_e1, psi_k1, psi_out, factor_minus, factor_e2, factor_k0)

  if(add_hartree .or. add_fxc) then 
    call X(calc_hvar_lr)(sys, ik, hvar, psi_e1, &
      factor_b, factor_e2, psi_out)
  end if
  
  POP_SUB(X(inhomog_BE))
end subroutine X(inhomog_BE)

! --------------------------------------------------------------------------
subroutine X(inhomog_KE)(sys, idir, ik, & 
  psi_e, psi_out)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir, ik
  R_TYPE,               intent(inout) :: psi_e(:,:,:) 
  R_TYPE,               intent(inout) :: psi_out(:,:,:) 

  R_TYPE :: factor_plus, factor_magn, factor_e

  PUSH_SUB(X(inhomog_KE))

  factor_plus = M_ONE
  factor_magn = -M_ONE
  factor_e = -M_ONE
   
  call X(inhomog_per_component)(sys, idir, ik, & 
    psi_e, psi_out, factor_plus, factor_e, factor_magn)

  POP_SUB(X(inhomog_KE))
end subroutine X(inhomog_KE)

! --------------------------------------------------------------------------

subroutine X(inhomog_K2)(sys, idir1, idir2, ik, & 
  psi_k1, psi_k2, psi_out)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir1, idir2
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: psi_k1(:,:,:) 
  R_TYPE,               intent(inout) :: psi_k2(:,:,:) 
  R_TYPE,               intent(inout) :: psi_out(:,:,:)

  R_TYPE :: factor_plus, factor_magn, factor_k
  integer :: ip, ist, idim
  R_TYPE, allocatable :: f_out(:,:), psi(:,:)
  type(pert_t) :: pert_kdotp2
  
  PUSH_SUB(X(inhomog_K2))

  SAFE_ALLOCATE(f_out(1:sys%gr%mesh%np, 1:sys%hm%d%dim))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  factor_plus = M_ONE
#if defined(R_TCOMPLEX)
  factor_k = -M_zI
#else
  factor_k = -M_ONE
#endif 
  factor_magn = -M_ONE

  call X(inhomog_per_component)(sys, idir1, ik, &
    psi_k2, psi_out, factor_plus, factor_k, factor_magn)
  
  call pert_init(pert_kdotp2, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%geo)
  call pert_setup_dir(pert_kdotp2,idir1,idir2)
  
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi) 
      call X(pert_apply_order_2)(pert_kdotp2, sys%namespace, sys%gr, sys%geo, sys%hm, ik, psi, f_out)
      if(idir1 == idir2) f_out(:,:) = M_HALF * f_out(:,:)
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, ist) = psi_out(ip, idim, ist) + &
            factor_plus * f_out(ip, idim)
        end do
      end do
    end if 
  end do 

  call pert_end(pert_kdotp2)
  SAFE_DEALLOCATE_A(f_out)
  SAFE_DEALLOCATE_A(psi)

  POP_SUB(X(inhomog_K2))
end subroutine X(inhomog_K2)

! --------------------------------------------------------------------------
subroutine X(inhomog_KB)(sys, idir, idir1, idir2, ik, & 
  psi_b, psi_k, psi_k1, psi_k2, psi_out)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir, idir1, idir2
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: psi_k1(:,:,:)
  R_TYPE,               intent(inout) :: psi_k2(:,:,:)
  R_TYPE,               intent(inout) :: psi_k(:,:,:) 
  R_TYPE,               intent(inout) :: psi_b(:,:,:)
  R_TYPE,               intent(inout) :: psi_out(:,:,:) 

  R_TYPE :: factor_plus, factor_minus, factor_magn, factor, factor_k
  integer :: ist_occ, ist, idim, ip
  R_TYPE, allocatable :: f_out1(:,:,:), f_out2(:,:,:)
  type(matrix_t):: vel_mat1, vel_mat2
  type(pert_t)  :: pert_kdotp2
  R_TYPE, allocatable :: psi(:,:,:)

  PUSH_SUB(X(inhomog_KB))
  
  ASSERT(sys%gr%sb%dim .ne. -1)
  ASSERT(idir .ne. -1)
  
  SAFE_ALLOCATE(f_out1(1:sys%gr%mesh%np,1:sys%hm%d%dim, 1:sys%st%nst))
  SAFE_ALLOCATE(f_out2(1:sys%gr%mesh%np,1:sys%hm%d%dim, 1:sys%st%nst))
  SAFE_ALLOCATE(vel_mat1%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(vel_mat2%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim, 1:sys%st%nst))
  
#if defined(R_TCOMPLEX)
  factor_plus = M_zI * M_HALF
  factor_minus = -M_zI * M_HALF
  factor_k = -M_zI
#else
  factor_plus = M_HALF
  factor_minus = -M_HALF
  factor_k = M_ONE
#endif

  factor_magn = -M_ONE
  factor = M_ONE

  call X(inhomog_per_component_2nd_order)(sys, idir, ik, &
    psi_k2, psi_k1, psi_out, factor_plus, factor_k, factor_k)
  call X(inhomog_per_component_2nd_order)(sys, idir, ik, &
    psi_k1, psi_k2, psi_out, factor_minus, factor_k, factor_k)
  call X(inhomog_per_component)(sys, idir, ik, & 
    psi_b, psi_out, factor, factor, factor_magn)

  call pert_init(pert_kdotp2, sys%namespace, PERTURBATION_KDOTP, sys%gr, sys%geo)
  call pert_setup_dir(pert_kdotp2, idir1, idir2)


  f_out1(:,:,:) = M_ZERO
  f_out2(:,:,:) = M_ZERO
  psi(:,:,:) = M_ZERO
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi(:, :, ist))
      call pert_setup_dir(pert_kdotp2, idir, idir1)
      call X(pert_apply_order_2)(pert_kdotp2, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
        psi_k2(:, :, ist), f_out1(:, :, ist))
      if(idir .ne. idir1) f_out1(:,:,ist) = M_TWO * f_out1(:, :, ist)   
      call pert_setup_dir(pert_kdotp2,idir,idir2)
      call X(pert_apply_order_2)(pert_kdotp2, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
        psi_k1(:, :, ist), f_out2(:, :, ist))
      if(idir .ne. idir2) f_out2(:,:,ist) = M_TWO * f_out2(:, :, ist)

      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, ist) = psi_out(ip, idim, ist) +  &
            factor_k * (factor_plus * f_out1(ip, idim, ist) + factor_minus * f_out2(ip, idim, ist))
        end do
      end do

      call pert_setup_dir(pert_kdotp2, idir, idir1)
      call X(pert_apply_order_2)(pert_kdotp2, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
        psi(:, :, ist), f_out1(:, :, ist))
      if(idir .ne. idir1) f_out1(:, :, ist) = M_TWO * f_out1(:, :, ist)  
      call pert_setup_dir(pert_kdotp2, idir, idir2)
      call X(pert_apply_order_2)(pert_kdotp2, sys%namespace, sys%gr, sys%geo, sys%hm, ik, &
        psi(:, :, ist), f_out2(:, :,ist))
      if(idir .ne. idir2) f_out2(:, :, ist) = M_TWO * f_out2(:, :, ist)
    end if
  end do
    
  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    psi, f_out1(:,:,:), vel_mat1%X(matrix))

  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    psi, f_out2(:,:,:), vel_mat2%X(matrix))

  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      do ist_occ = 1, sys%st%nst
        if(abs(sys%st%occ(ist_occ, ik)) .gt. M_EPSILON) then
          do idim = 1, sys%hm%d%dim
            do ip = 1, sys%gr%mesh%np
              psi_out(ip, idim, ist) = psi_out(ip, idim, ist) +  &
                factor_k * (factor_plus * psi_k2(ip, idim, ist_occ) * vel_mat1%X(matrix)(ist_occ, ist) + &
                factor_minus * psi_k1(ip, idim, ist_occ) * vel_mat2%X(matrix)(ist_occ, ist))
            end do  
          end do
        end if
      end do
    end if
  end do

  call pert_end(pert_kdotp2)
  
  SAFE_DEALLOCATE_P(vel_mat1%X(matrix))
  SAFE_DEALLOCATE_P(vel_mat2%X(matrix))
  SAFE_DEALLOCATE_A(f_out1)
  SAFE_DEALLOCATE_A(f_out2)
  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(inhomog_KB))
end subroutine X(inhomog_KB)


 ! --------------------------------------------------------------------------
subroutine X(inhomog_KB_tot)(sh, sys, idir, idir1, idir2, & 
  lr_k, lr_b, lr_k1, lr_k2, lr_kk1, lr_kk2, psi_out)
  type(sternheimer_t),  intent(inout) :: sh
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir, idir1, idir2
  type(lr_t),           intent(inout) :: lr_k(:)
  type(lr_t),           intent(inout) :: lr_b(:)
  type(lr_t),           intent(inout) :: lr_kk1(:)
  type(lr_t),           intent(inout) :: lr_kk2(:)
  type(lr_t),           intent(inout) :: lr_k1(:)
  type(lr_t),           intent(inout) :: lr_k2(:)
  R_TYPE,               intent(inout) :: psi_out(:,:,:,:,:) 

  R_TYPE :: factor1, factor2, factor_sum, factor_rho
  R_TYPE, allocatable :: hvar(:,:,:)
  type(lr_t) :: lr0(1)
  integer :: ispin, ik, ip, ik0
  logical :: add_hartree, add_fxc

  PUSH_SUB(X(inhomog_KB_tot))

  SAFE_ALLOCATE(hvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))

#if defined(R_TCOMPLEX)
  factor1 = -M_zI
  factor2 = -M_zI
  factor_rho = M_zI * M_HALF
#else 
  factor1 = M_ONE
  factor2 = -M_ONE
  factor_rho = M_HALF
#endif 
  factor_sum = -M_ONE 

  psi_out(:,:,:,:,:) = M_ZERO
  hvar(:,:,:) = M_ZERO

  add_hartree = sternheimer_add_hartree(sh) 
  add_fxc = sternheimer_add_fxc(sh)
  
  if(add_hartree .or. add_fxc) then
    call lr_init(lr0(1))
    call lr_allocate(lr0(1), sys%st, sys%gr%mesh, allocate_rho = .true.)
    lr0(1)%X(dl_rho) = M_ZERO
  
    call X(calc_rho)(sys, factor_rho, factor_sum, factor1, &
      factor1, lr_k1(1), lr_k2(1), lr0(1))
    call X(calc_rho)(sys, factor_sum * factor_rho, factor_sum, &
      factor1, factor1, lr_k2(1), lr_k1(1), lr0(1))
    if(sys%st%parallel_in_states .or. sys%st%d%kpt%parallel) &  
      call comm_allreduce(sys%st%st_kpt_mpi_grp%comm, lr0(1)%X(dl_rho))

    do ip = 1, sys%gr%mesh%np 
      do ispin = 1, sys%st%d%nspin 
       lr0(1)%X(dl_rho)(ip, ispin) = lr0(1)%X(dl_rho)(ip, ispin) + lr_b(1)%X(dl_rho)(ip, ispin)
      end do  
    end do
   
    call X(sternheimer_calc_hvar)(sh, sys, lr0, 1, hvar)
    call lr_dealloc(lr0(1))
  end if
  
  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ik0 = ik-sys%st%d%kpt%start+1
    call X(inhomog_KB)(sys, idir, idir1, idir2, ik, & 
      lr_b(1)%X(dl_psi)(:, :, :, ik), lr_k(1)%X(dl_psi)(:, :, :, ik), &
      lr_k1(1)%X(dl_psi)(:, :, :, ik), lr_k2(1)%X(dl_psi)(:, :, :, ik), &
      psi_out(:, :, :, ik0, 1))
    
    ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
    call X(inhomog_BE)(sys, idir1, idir2, ik, & 
      add_hartree, add_fxc, hvar(:, ispin, 1),  &
      lr_k(1)%X(dl_psi)(:, :, :, ik), lr_k(1)%X(dl_psi)(:, :, :, ik), &
      lr_kk1(1)%X(dl_psi)(:, :, :, ik), lr_kk2(1)%X(dl_psi)(:, :, :, ik), &
      lr_k1(1)%X(dl_psi)(:, :, :, ik), lr_k2(1)%X(dl_psi)(:, :, :, ik), &
      factor1, factor2, psi_out(:, :, :, ik0, 1))
  end do
  
  SAFE_DEALLOCATE_A(hvar)

  POP_SUB(X(inhomog_KB_tot))
end subroutine X(inhomog_KB_tot)

! --------------------------------------------------------------------------
subroutine X(inhomog_KE_tot)(sh, sys, idir, nsigma, & 
  lr_k, lr_e, lr_kk, psi_out)
  type(sternheimer_t),  intent(inout) :: sh
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir, nsigma
  type(lr_t),           intent(inout) :: lr_e(:) 
  type(lr_t),           intent(inout) :: lr_k(:)  
  type(lr_t),           intent(inout) :: lr_kk(:)  
  R_TYPE,               intent(inout) :: psi_out(:,:,:,:,:) 
   
  R_TYPE :: factor_k
  integer :: isigma, ispin, ik, ik0
  logical :: add_hartree, add_fxc
  R_TYPE, allocatable :: hvar(:,:,:)

  PUSH_SUB(X(inhomog_KE_tot))

  SAFE_ALLOCATE(hvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:nsigma))

#if defined(R_TCOMPLEX)
  factor_k = -M_zI
#else 
  factor_k = M_ONE
#endif 

  psi_out(:,:,:,:,:) = M_ZERO
  hvar(:,:,:) = M_ZERO
 
  add_hartree = sternheimer_add_hartree(sh) 
  add_fxc = sternheimer_add_fxc(sh)
  
  if(add_hartree .or. add_fxc) then
    call X(sternheimer_calc_hvar)(sh, sys, lr_e, nsigma, hvar)
  end if

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ik0 = ik-sys%st%d%kpt%start+1
    ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
    do isigma = 1, nsigma
      call X(inhomog_EB)(sys, ik, add_hartree, add_fxc, &
        hvar(:, ispin, isigma), lr_k(1)%X(dl_psi)(:, :, :, ik), &
        lr_kk(1)%X(dl_psi)(:, :, :, ik), factor_k, psi_out(:, :, :, ik0, isigma))

      call X(inhomog_KE)(sys, idir, ik, & 
        lr_e(isigma)%X(dl_psi)(:, :, :, ik), psi_out(:, :, :, ik0, isigma))
    end do
  end do
  
  SAFE_DEALLOCATE_A(hvar)

  POP_SUB(X(inhomog_KE_tot))
end subroutine X(inhomog_KE_tot)
 
! --------------------------------------------------------------------------
subroutine X(inhomog_K2_tot)(sh, sys, idir1, idir2, & 
  lr_k1, lr_k2, psi_out)
  type(sternheimer_t),  intent(inout) :: sh
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: idir1, idir2
  type(lr_t),           intent(inout) :: lr_k1(:) 
  type(lr_t),           intent(inout) :: lr_k2(:) 
  R_TYPE,               intent(inout) :: psi_out(:,:,:,:,:) 

  integer :: ik, ik0
  
  PUSH_SUB(X(inhomog_K2_tot))

  psi_out(:,:,:,:,:) = M_ZERO

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ik0 = ik-sys%st%d%kpt%start+1
    call X(inhomog_K2)(sys, idir1, idir2, ik, & 
      lr_k1(1)%X(dl_psi)(:,:,:, ik), lr_k2(1)%X(dl_psi)(:,:,:, ik), &
      psi_out(:,:,:, ik0, 1))

    if(idir1 == idir2) then
      psi_out(:,:,:, ik0, 1) = M_TWO * psi_out(:,:,:, ik0, 1)
    else
      call X(inhomog_K2)(sys, idir2, idir1, ik, & 
        lr_k2(1)%X(dl_psi)(:,:,:, ik), lr_k1(1)%X(dl_psi)(:,:,:, ik), &
        psi_out(:,:,:, ik0, 1))
    end if
  end do

  POP_SUB(X(inhomog_K2_tot))
end subroutine X(inhomog_K2_tot)

!----------------------------------------------------------
! Calculation of contribution to the density for second-order perturbations
! or magnetic perturbations with kdotp that come from elements of the density 
! matrix within occupied and unoccupied subspaces
subroutine X(calc_rho)(sys, factor, factor_sum, factor_e, &
  factor_k, lr_e, lr_k, lr0)
  type(electrons_t),    intent(inout) :: sys
  R_TYPE,               intent(in)    :: factor
  R_TYPE,               intent(in)    :: factor_sum
  R_TYPE,               intent(in)    :: factor_e
  R_TYPE,               intent(in)    :: factor_k
  type(lr_t),           intent(inout) :: lr_e 
  type(lr_t),           intent(inout) :: lr_k 
  type(lr_t),           intent(inout) :: lr0 

  integer :: ip, ik, ist, idim, ist_occ, ispin
  FLOAT :: weight
  type(matrix_t):: mat
  R_TYPE, allocatable :: psi(:,:), psi1(:,:)

  PUSH_SUB(X(calc_rho))

  SAFE_ALLOCATE(mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  SAFE_ALLOCATE(psi1(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  do ik = sys%st%d%kpt%start, sys%st%d%kpt%end
    ispin = states_elec_dim_get_spin_index(sys%st%d, ik)
    call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
      factor_k * lr_k%X(dl_psi)(:,:,:, ik), factor_e*lr_e%X(dl_psi)(:,:,:, ik), mat%X(matrix))

    do ist  = 1, sys%st%nst
      if (sys%st%occ(ist, ik) .gt. M_EPSILON) then
        call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
        weight = sys%st%d%kweights(ik) * sys%st%smear%el_per_state
        do idim = 1, sys%hm%d%dim
          do ip = 1, sys%gr%mesh%np
            lr0%X(dl_rho)(ip, ispin) = lr0%X(dl_rho)(ip, ispin) + factor * weight * factor_e * &
              lr_e%X(dl_psi)(ip, idim, ist, ik) * R_CONJ(factor_k * lr_k%X(dl_psi)(ip, idim, ist, ik))
          end do
        end do
        do ist_occ  = 1, sys%st%nst
          if (sys%st%occ(ist_occ, ik) .gt. M_EPSILON) then
            call states_elec_get_state(sys%st, sys%gr%mesh, ist_occ, ik, psi1)
            do idim = 1, sys%hm%d%dim
              do ip = 1, sys%gr%mesh%np
                lr0%X(dl_rho)(ip, ispin) = lr0%X(dl_rho)(ip, ispin) - factor * factor_sum * weight * &
                  mat%X(matrix)(ist, ist_occ) * psi(ip, idim) * &
                  R_CONJ(psi1(ip, idim))
              end do
            end do
          end if
        end do
      end if
    end do
  end do

  SAFE_DEALLOCATE_P(mat%X(matrix))
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi1)

  POP_SUB(X(calc_rho))
end subroutine X(calc_rho)

! --------------------------------------------------------------------------
! Calculation of V_{hxc}[n^{(1)}] | \psi^{(0)}>  for magnetic perturbations 
! with kdotp coming from the contribution to the density from elements of
! the density matrix within occupied and unoccupied subspaces
subroutine X(calc_hvar_psi)(sys, ik, hvar, psi_out)
  type(electrons_t),    intent(inout) :: sys 
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: hvar(:) 
  R_TYPE,               intent(inout) :: psi_out(:,:,:)   
    
  R_TYPE, allocatable :: psi(:,:)
  integer :: ip, ist, idim
    
  PUSH_SUB(X(calc_hvar_psi))

  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
  
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi)
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, ist) = psi_out(ip, idim, ist) &
            - hvar(ip) * psi(ip, idim)
        end do
      end do
    end if
  end do

  SAFE_DEALLOCATE_A(psi)
  
  POP_SUB(X(calc_hvar_psi))
end subroutine X(calc_hvar_psi)

! --------------------------------------------------------------------------
! Calculation of V_{hxc}[n^{(1)}] |  \psi^{(1)}> for second-order perturbations
subroutine X(calc_hvar_lr)(sys, ik, hvar, psi_in, &
  factor1, factor2, psi_out)
  type(electrons_t),    intent(inout) :: sys
  integer,              intent(in)    :: ik
  R_TYPE,               intent(inout) :: hvar(:)
  R_TYPE,               intent(inout) :: psi_in(:,:,:)
  R_TYPE,               intent(in)    :: factor1, factor2
  R_TYPE,               intent(inout) :: psi_out(:,:,:)
    
  R_TYPE, allocatable :: psi(:,:,:)
  integer :: ip, ist, idim, ist1
  type(matrix_t):: mat
  R_TYPE, allocatable :: psi0(:,:,:)
    
  PUSH_SUB(X(calc_hvar_lr))
  
  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np,1:sys%hm%d%dim, 1:sys%st%nst))
  SAFE_ALLOCATE(mat%X(matrix)(1:sys%st%nst, 1:sys%st%nst))
  SAFE_ALLOCATE(psi0(1:sys%gr%mesh%np, 1:sys%st%d%dim, 1:sys%st%nst))

  psi(:,:,:) = M_ZERO
  psi0(:,:,:) = M_ZERO
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi0(:, :, ist))
      do idim = 1, sys%hm%d%dim
        do ip = 1, sys%gr%mesh%np
          psi_out(ip, idim, ist) = psi_out(ip, idim, ist) &
            - factor1 * hvar(ip) * factor2 * psi_in(ip, idim, ist)
            psi(ip, idim, ist) = factor1 * hvar(ip) * psi0(ip, idim, ist)
        end do
      end do
    end if
  end do

  call states_elec_blockt_mul(sys%gr%mesh, sys%st, sys%st%st_start, sys%st%st_start, &
    psi0, psi(:,:,:), mat%X(matrix))
        
  do ist = 1, sys%st%nst
    if(abs(sys%st%occ(ist, ik)) .gt. M_EPSILON) then
      do ist1 = 1, sys%st%nst
        if(abs(sys%st%occ(ist1, ik)) .gt. M_EPSILON) then
          do idim = 1, sys%hm%d%dim
            do ip = 1, sys%gr%mesh%np
              psi_out(ip, idim, ist) = psi_out(ip, idim, ist) + &
                mat%X(matrix)(ist1, ist) * factor2 * psi_in(ip, idim, ist1)
            end do
          end do
        end if 
      end do
    end if
  end do
   
  SAFE_DEALLOCATE_P(mat%X(matrix))
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(psi0)
  
  POP_SUB(X(calc_hvar_lr))
end subroutine X(calc_hvar_lr)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
