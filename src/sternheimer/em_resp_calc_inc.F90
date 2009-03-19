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
! WARNING: This subroutine is clearly broken after the changes
! to include temperature in linear response
subroutine X(lr_calc_elf)(st, gr, lr, lr_m)
  type(states_t),   intent(inout) :: st
  type(grid_t),     intent(inout) :: gr
  type(lr_t),       intent(inout) :: lr
  type(lr_t), optional, intent(inout) :: lr_m !when this argument is present, we are doing dynamical response

  integer :: i, is, ist, idim, ik

  R_TYPE, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)
  FLOAT,  allocatable :: rho(:), grho(:,:)
  R_TYPE, allocatable :: dl_rho(:), gdl_rho(:,:)
  FLOAT,  allocatable :: elf(:,:), de(:,:)
  FLOAT :: dl_d0, d0
  FLOAT :: f, s

  FLOAT, parameter :: dmin = CNST(1e-10)
  FLOAT :: ik_weight

  call push_sub('em_resp_inc.Xcalc_lr_elf')

  ASSERT(.false.)

  ALLOCATE(   gpsi(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)
  ALLOCATE(gdl_psi(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)

  if(present(lr_m)) ALLOCATE(gdl_psi_m(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)

  ALLOCATE(   grho(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)
  ALLOCATE(gdl_rho(NP, gr%mesh%sb%dim), NP*gr%mesh%sb%dim)

  ALLOCATE(   rho(gr%mesh%np_part), NP)
  ALLOCATE(dl_rho(gr%mesh%np_part), NP)

  ALLOCATE(   elf(NP, st%d%nspin), NP*st%d%nspin)
  ALLOCATE(    de(NP, st%d%nspin), NP*st%d%nspin)

  if( .not. associated(lr%X(dl_de)) ) ALLOCATE(lr%X(dl_de)(NP, st%d%nspin), NP*st%d%nspin) 
  if( .not. associated(lr%X(dl_elf))) ALLOCATE(lr%X(dl_elf)(NP, st%d%nspin), NP*st%d%nspin)

  !calculate the gs elf
  call elf_calc(st, gr, elf, de)

  !calculate current and its variation
  if(st%wfs_type == M_CMPLX) then 
    call calc_physical_current(gr, st, st%j)
    if(present(lr_m)) then 
      call lr_calc_current(st, gr, lr, lr_m)
    else
      call lr_calc_current(st, gr, lr)
    end if
  end if

  ! single or double occupancy
  if(st%d%nspin == 1) then
    s = M_TWO
  else
    s = M_ONE
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
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/s

        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          ! sum over states to obtain the spin-density
          rho(1:NP) = rho(1:NP) + ik_weight*abs(st%X(psi)(1:NP, idim, ist, is))**2

          !the gradient of the density
          do i = 1, gr%mesh%sb%dim
            grho(1:NP, i) = grho(1:NP, i) + &
                 ik_weight*M_TWO*R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is))*gpsi(1:NP, i))
          end do

          !the variation of the density and its gradient

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            dl_rho(1:NP) = dl_rho(1:NP) + ik_weight * ( &
                 R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * lr%X(dl_psi)(1:NP, idim, ist, is)+ & 
                 st%X(psi)(1:NP, idim, ist, is) * R_CONJ(lr_m%X(dl_psi)(1:NP, idim, ist, is)) )

            do i=1, gr%mesh%sb%dim

              gdl_rho(1:NP,i) = gdl_rho(1:NP,i) + ik_weight * ( &
                   R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gdl_psi(1:NP,i) +      &
                   R_CONJ(gpsi(1:NP,i))* lr%X(dl_psi)(1:NP, idim, ist, is)  +      &
                   st%X(psi)(1:NP, idim, ist, is) * R_CONJ(gdl_psi_m(1:NP,i)) +      &
                   gpsi(1:NP,i) * R_CONJ(lr_m%X(dl_psi)(1:NP, idim, ist, is))  )

            end do

          else
            
            dl_rho(1:NP) = dl_rho(1:NP) + ik_weight*M_TWO* &
                 R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is))*lr%X(dl_psi)(1:NP, idim, ist, is))

            do i=1,gr%mesh%sb%dim
              gdl_rho(1:NP,i) = gdl_rho(1:NP,i) + ik_weight * M_TWO * ( &
                   R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gdl_psi(1:NP, i) + &
                   gpsi(1:NP, i)*R_CONJ(lr%X(dl_psi)(1:NP, idim, ist, is))  )
            end do

          end if

        end do !idim
      end do !ist
    end do !ik

    !now we start to calculate the elf

    !first the term that depends on the orbitals
    !this is the only term that is different for the dynamical case
    do ik = is, st%d%nik, st%d%nspin
      do ist = 1, st%nst
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/s
        do idim = 1, st%d%dim

          call X(derivatives_grad)(gr%der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(derivatives_grad)(gr%der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          if(present(lr_m)) then 

            call X(derivatives_grad)(gr%der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            do i = 1, NP
              lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) + &
                   dl_rho(i)*ik_weight*sum(R_ABS(gpsi(i, 1:gr%mesh%sb%dim))**2) + &
                   rho(i)*ik_weight*&
                   sum(R_CONJ(gpsi(i,1:gr%mesh%sb%dim))*gdl_psi(i, 1:gr%mesh%sb%dim) + &
                   gpsi(i,1:gr%mesh%sb%dim)*R_CONJ(gdl_psi_m(i, 1:gr%mesh%sb%dim)))
            end do

          else 

            do i = 1, NP
              lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) + &
                   dl_rho(i)*ik_weight*sum(R_ABS(gpsi(i, 1:gr%mesh%sb%dim))**2) + &
                   rho(i)*ik_weight*M_TWO*(sum(R_CONJ(gpsi(i, 1:gr%mesh%sb%dim))*gdl_psi(i, 1:gr%mesh%sb%dim)))
            end do

          end if

        end do
      end do
    end do

    !the density term
    do i = 1, NP
      if(abs(st%rho(i, is)) >= dmin) then
        lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) - M_HALF*sum(grho(i, 1:gr%mesh%sb%dim)*gdl_rho(i, 1:gr%mesh%sb%dim))
      end if
    end do

    !the current term
    if(st%wfs_type == M_CMPLX) then       
      do i = 1, NP
        if(abs(st%rho(i, is)) >= dmin) then
          lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) + M_TWO*sum(st%j(i, 1:gr%mesh%sb%dim, is)*lr%dl_j(i, 1:gr%mesh%sb%dim, is))
        end if
      end do
    end if

    !now the normalization 
    f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
    do i = 1, NP

      if(abs(st%rho(i, is)) >= dmin) then
        d0    = f*rho(i)**(M_EIGHT/M_THREE)
        dl_d0 = M_EIGHT/M_THREE*f*dl_rho(i)*rho(i)**(M_FIVE/M_THREE)

        lr%X(dl_elf)(i, is) = &
             M_TWO*d0*dl_d0/(d0**2 + de(i, is)**2)*(M_ONE - elf(i,is)) &
             - M_TWO*de(i, is)*lr%X(dl_de)(i, is)/(d0**2 + de(i, is)**2)*elf(i, is)
      else
        lr%X(dl_elf)(i, is) = M_ZERO
      end if

    end do

  end do

  deallocate(gpsi)
  deallocate(gdl_psi)
  if(present(lr_m)) deallocate(gdl_psi_m)

  deallocate(rho)
  deallocate(dl_rho)

  deallocate(grho)
  deallocate(gdl_rho)

  deallocate(elf)
  deallocate(de)

  call pop_sub()

end subroutine X(lr_calc_elf)


! ---------------------------------------------------------
subroutine X(calc_polarizability_finite)(sys, hm, lr, nsigma, perturbation, zpol, ndir)
  type(system_t),         intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: hm
  type(lr_t),             intent(inout) :: lr(:,:)
  integer,                intent(in)    :: nsigma
  type(pert_t),           intent(inout) :: perturbation
  CMPLX,                  intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)
  integer, optional,      intent(in)    :: ndir

  integer :: dir1, dir2, ndir_

  call push_sub('em_resp_calc_inc.Xcalc_polarizability_finite')

  ndir_ = sys%gr%mesh%sb%dim
  if(present(ndir)) ndir_ = ndir

  ! alpha_ij(w) = - sum(m occ) [<psi_m(0)|r_i|psi_mj(1)(w)> + <psi_mj(1)(-w)|r_i|psi_m(0)>]
  ! minus sign is from electronic charge -e

  do dir1 = 1, ndir_
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

  call pop_sub()

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

  call push_sub('em_resp_calc_inc.Xlr_calc_susceptibility')

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
  ! When using real wave-functions there is an extra factor of (-i)*(-i) = -1
  chi_para(:,:) = -chi_para(:,:)
#endif

  call pop_sub()

end subroutine X(lr_calc_susceptibility)


! ---------------------------------------------------------
subroutine X(lr_calc_beta) (sh, sys, hm, em_lr, dipole, beta, kdotp_lr, kdotp_em_lr)
! Note: correctness for spinors is unclear
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
  type(lr_t),    optional, intent(in)    :: kdotp_em_lr(:,:,:,:)

  type(states_t), pointer :: st
  type(mesh_t),   pointer :: mesh

  integer :: ifreq, isigma, idim, ispin, np, ndim, idir, ist, ik
  integer :: ii, jj, kk, iperm, op_sigma, idim2, ist2, ip
  integer :: perm(1:3), u(1:3), w(1:3), ijk(1:3)

  R_TYPE :: prod

  R_TYPE, allocatable :: hvar(:, :, :, :, :, :)
  R_TYPE, allocatable :: tmp(:)
  FLOAT,  allocatable :: rho(:,:), kxc(:, :, :, :)
  R_TYPE, allocatable :: hpol_density(:)

  call push_sub('em_resp_inc.Xlr_calc_beta')

  call profiling_in(beta_prof, "CALC_BETA")

  np   =  sys%gr%mesh%np
  ndim =  sys%gr%sb%dim
  st   => sys%st
  mesh => sys%gr%mesh

  write(message(1), '(a)') 'Info: Calculating hyperpolarizability tensor'
  call write_info(1)

  ASSERT(present(kdotp_lr) .eqv. present(kdotp_em_lr))
  ! either both are absent for finite, or both present for periodic

  !calculate kxc, the derivative of fxc
  ALLOCATE(kxc(1:np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin), np*st%d%nspin**3)
  kxc = M_ZERO

  ALLOCATE(rho(np, st%d%nspin), np*st%d%nspin)
  call states_total_density(st, mesh, rho)

  call xc_get_kxc(sys%ks%xc, mesh, rho, st%d%ispin, kxc)
  deallocate(rho)

  ALLOCATE(tmp(1:np), np)
  ALLOCATE(hvar(1:np, 1:st%d%nspin, 1:2, 1:st%d%dim, 1:ndim, 1:3), np*st%d%nspin*2*st%d%dim*ndim*3)
  ALLOCATE(hpol_density(1:np), np)

  do ifreq = 1, 3
    do idir = 1, ndim
      do idim = 1, st%d%dim
        call X(sternheimer_calc_hvar)(sh, sys, hm, em_lr(idir, :, ifreq), 2, hvar(:, :, :, idim, idir, ifreq))
      end do !idim
    end do !idir
  end do !ifreq

  beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM) = M_ZERO

  do ii = 1, ndim
    do jj = 1, ndim
      do kk = 1, ndim

        ijk = (/ ii, jj, kk /)

        do isigma = 1, 2

          op_sigma = 2 
          if(isigma == 2) op_sigma = 1

          do iperm = 1, 6

            call get_permutation(iperm, perm)

            u(1:3) = ijk(perm(1:3))
            w(1:3) = perm(1:3)

            do ik = 1, st%d%nik
              do ist = 1, st%nst
                do idim = 1, st%d%dim

                  ispin = states_dim_get_spin_index(sys%st%d, ik)

                  if (present(kdotp_em_lr)) then
                    tmp(1:np) = - M_zI * kdotp_em_lr(u(2), u(3), isigma, w(3))%X(dl_psi)(1:np, idim, ist, ik)
                  else
                    call pert_setup_dir(dipole, u(2))
                    call X(pert_apply) &
                      (dipole, sys%gr, sys%geo, hm, ik, em_lr(u(3), isigma, w(3))%X(dl_psi)(:, idim, ist, ik), tmp)
                  endif

                  do ip = 1, np
                    tmp(ip) = tmp(ip) + R_REAL(hvar(ip, ispin, isigma, idim, u(2), w(2))) &
                      * em_lr(u(3), isigma, w(3))%X(dl_psi)(ip, idim, ist, ik)
                  enddo

                  beta(ii, jj, kk) = beta(ii, jj, kk) &
                    - M_HALF * st%d%kweights(ik) * st%smear%el_per_state &
                    * X(mf_dotp)(mesh, em_lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, idim, ist, ik), tmp(1:np))
                  
                  do ist2 = 1, st%nst
                    do idim2 = 1, st%d%dim
                      ! there is no coupling between states with different k or spin

                      if (present(kdotp_lr)) then
                        tmp(1:np) = - M_zI * kdotp_lr(u(2))%X(dl_psi)(1:np, idim, ist, ik) 
                      else
                        call pert_setup_dir(dipole, u(2))
                        call X(pert_apply)(dipole, sys%gr, sys%geo, hm, ik, st%X(psi)(:, idim, ist, ik), tmp)
                      endif

                      do ip = 1, np
                        tmp(ip) = tmp(ip) + &
                          R_REAL(hvar(ip, ispin, isigma, idim2, u(2), w(2))) * st%X(psi)(ip, idim, ist, ik)
                      enddo

                      prod = X(mf_dotp)(mesh, st%X(psi)(:, idim2, ist2, ik), tmp)

                      beta(ii, jj, kk) = beta(ii, jj, kk) + & 
                        M_HALF * st%d%kweights(ik) * st%smear%el_per_state * prod * X(mf_dotp)(mesh, &
                        em_lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, idim, ist, ik), &
                        em_lr(u(3), isigma,   w(3))%X(dl_psi)(1:np, idim2, ist2, ik))

                      end do !idim2
                    end do ! ist2

                end do !idim
              end do !ist
            end do !ik

            if(sternheimer_add_fxc(sh)) then 
              do ip = 1, np 
                hpol_density(ip) = kxc(ip, 1, 1, 1) & 
                  * sum(em_lr(u(1), isigma, w(1))%X(dl_rho)(ip, 1:st%d%nspin)) & 
                  * sum(em_lr(u(2), isigma, w(2))%X(dl_rho)(ip, 1:st%d%nspin)) & 
                  * sum(em_lr(u(3), isigma, w(3))%X(dl_rho)(ip, 1:st%d%nspin)) & 
                  / CNST(6.0) 
              end do
            end if

            beta(ii, jj, kk) = beta(ii, jj, kk) - M_HALF * X(mf_integrate)(mesh, hpol_density)


          end do ! iperm

        end do !isigma

      end do !kk
    end do !jj
  end do !ii

  deallocate(hpol_density)
  deallocate(tmp)
  deallocate(hvar)

  call pop_sub()

  call profiling_out(beta_prof)

contains

  subroutine get_permutation(i, p)
    integer, intent(in)  :: i
    integer, intent(out) :: p(1:3)

    ASSERT( i>=1 .and. i <= 6)

    select case(i)

    case(1) ; p(1)=1 ; p(2)=2 ; p(3)=3
    case(2) ; p(1)=2 ; p(2)=3 ; p(3)=1
    case(3) ; p(1)=3 ; p(2)=1 ; p(3)=2
    case(4) ; p(1)=3 ; p(2)=2 ; p(3)=1
    case(5) ; p(1)=1 ; p(2)=3 ; p(3)=2
    case(6) ; p(1)=2 ; p(2)=1 ; p(3)=3

    end select

  end subroutine get_permutation

end subroutine X(lr_calc_beta)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
